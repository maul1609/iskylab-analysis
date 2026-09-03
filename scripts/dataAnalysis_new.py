"""Initialise, batch-run and analyse iSKYLAB experiments with the BMM.

This is the maintained iSKYLAB batch driver.  It reads the supplied chamber
measurements and initial aerosol PSDs, creates one BMM namelist per experiment,
optionally runs ``main.exe``, and supports both historical batch-number
comparisons and detailed single-experiment time-series/PSD validation.

Important implementation notes
------------------------------
* Namelists are edited by variable name, never by matching a copied historical
  value.  Missing variables therefore raise an error instead of silently using
  the template value.
* The whole ``&chamber_spec`` block is regenerated for each experiment, so the
  chamber data length is no longer tied to an old 1853-point template.
* Pressure, gas temperature and qtot forcing are independent model controls.
* qtot data can be written without forcing the model.  By default qtot is
  constructed as dew-point vapour + OPC liquid water and forcing is OFF.
* Measured pressure, gas temperature and wall temperature can be smoothed
  with time-based Savitzky-Golay windows before namelist export.  The raw
  observations are retained and diagnostic plots/summary statistics can be
  written for every experiment.
* The measured ``Tww_mean`` series is written automatically when the coupled
  coupled wall-BL treatment is enabled with non-zero thermal coupling, or when
  wall-temperature BL mode explicitly requires it.
* Aerosol diameters from the PNSD fitting code are in micrometres internally
  and are converted to metres exactly once when the namelist is written.
* ``--experiment Exp005`` runs/analyses one experiment against ql, Nd, Deff,
  relative dispersion and the full OPC size distribution.  WELAS and native
  BMM moving particles are conservatively rebinned onto a configurable common
  coarse log-D grid for a readable like-for-like PSD comparison.
* ``--cloud-formation-time-min`` (alias ``--saturation-time-min``) can move the
  model saturation target for a single run without moving the raw observed
  cloud timing.
* ``--shift-cloud-measurements`` optionally applies the diagnosed WELAS/BMM ql
  lag to the complete WELAS time coordinate, so liquid water, number, Deff,
  dispersion, ice number and the full observed PSD are shifted consistently.
* ``--bl-tau-s`` and ``--bl-temp-offset-k`` override the chamber BL mixing
  timescale and sensible-temperature offset for a run without editing config.
* ``--bl-wall-water-mode 1`` enables a finite prognostic chamber-wall water
  reservoir.  The wall starts dry by default, accumulates liquid/frost only by
  condensation/deposition, and can subsequently return no more water than it
  has stored.  Reservoir masses are physical kg, not kg kg-1 mixing ratios.
* ``--bl-wall-water-mode 1`` retains the finite-reservoir fractional-relaxation
  closure for reproducibility. ``--bl-wall-water-mode 2`` uses the same finite
  prognostic reservoir but computes a physical wall vapour mass flux from a
  transfer velocity and the vapour-pressure disequilibrium with measured Twall.
  In mode 2 wall vapour exchange is independent of ``chamber_bl_tau``.
* ``--sce-bins``/``--scebins`` creates a temporary copy of ``sce/namelist.in``,
  changes ``n_binsc`` only in that copy, and points the generated BMM namelist
  ``scefile`` at it.  The repository SCE namelist is never modified.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, savgol_filter, sosfiltfilt

import experiment_metadata as meta
import iskylab_config as cfg
import readMeteoCPC
import readOPC_Merged
import readPNSD_Mrg_new
import svp
from namelist_utils import (
    read_text, replace_group, set_array, set_literal_array, set_or_insert_value,
    set_value, write_text,
)

R_GAS = 8.314
M_WATER = 18.0e-3
M_DRY_AIR = 28.96e-3
R_D = R_GAS / M_DRY_AIR
R_V = R_GAS / M_WATER
EPSILON = R_D / R_V



def _smoothing_cfg(name, default):
    """Return a smoothing config value without requiring an updated config file."""
    return getattr(cfg, name, default)


def _time_sampling_info(time, *, name="series"):
    """Validate a time coordinate and return its representative sampling interval."""
    time = np.asarray(time, dtype=float)
    if time.ndim != 1 or time.size < 3:
        raise ValueError(f"{name}: time coordinate must contain at least three samples")
    if not np.all(np.isfinite(time)):
        raise ValueError(f"{name}: time coordinate contains non-finite values")
    dt_all = np.diff(time)
    if np.any(dt_all <= 0.0):
        raise ValueError(f"{name}: time coordinate must be strictly increasing")
    dt = float(np.median(dt_all))
    rel_irregularity = float(np.max(np.abs(dt_all - dt)) / dt)
    return dt, rel_irregularity


def _window_samples(time, window_seconds, *, minimum=3, name="series"):
    """Convert a physical smoothing/despiking window to a sensible odd sample count."""
    dt, _ = _time_sampling_info(time, name=name)
    if window_seconds <= 0.0:
        raise ValueError(f"{name}: smoothing window must be positive")
    nwin = max(int(round(float(window_seconds) / dt)), int(minimum))
    if nwin % 2 == 0:
        nwin += 1
    nmax = len(time) if len(time) % 2 == 1 else len(time) - 1
    nwin = min(nwin, nmax)
    if nwin < minimum:
        raise ValueError(f"{name}: record is too short for requested filtering")
    return int(nwin), dt


def _hampel_time_series(time, values, window_seconds, nsigma, *, name="series"):
    """Robustly replace isolated spikes using a local median/MAD criterion.

    This is deliberately a despiker, not a smoother.  Genuine slower chamber
    structure is retained; only points that are strong local outliers are
    replaced by the local median.
    """
    time = np.asarray(time, dtype=float)
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or values.size != time.size:
        raise ValueError(f"{name}: time and values must be one-dimensional and equal length")

    nwin, dt = _window_samples(time, window_seconds, minimum=3, name=name)
    half = nwin // 2
    padded = np.pad(values, (half, half), mode="edge")
    windows = np.lib.stride_tricks.sliding_window_view(padded, nwin)
    med = np.median(windows, axis=1)
    mad = np.median(np.abs(windows - med[:, None]), axis=1)

    # 1.4826*MAD estimates sigma for Gaussian noise.  The small numerical floor
    # still allows an isolated spike in an otherwise flat segment to be caught.
    scale_floor = 50.0 * np.finfo(float).eps * np.maximum(1.0, np.abs(med))
    sigma = np.maximum(1.4826 * mad, scale_floor)
    flagged = np.abs(values - med) > float(nsigma) * sigma

    out = values.copy()
    out[flagged] = med[flagged]
    return out, {
        "despike_window_seconds": float(window_seconds),
        "despike_window_samples": int(nwin),
        "despike_nsigma": float(nsigma),
        "n_despiked": int(np.count_nonzero(flagged)),
        "median_dt_s": float(dt),
    }


def _savgol_time_series(time, values, window_seconds, polyorder, *, name="series"):
    """Return zero-lag Savitzky-Golay-smoothed data using a physical-time window."""
    time = np.asarray(time, dtype=float)
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or len(time) != len(values):
        raise ValueError(f"{name}: time and values must be one-dimensional and equal length")
    if polyorder < 0:
        raise ValueError(f"{name}: polynomial order must be non-negative")

    dt, rel_irregularity = _time_sampling_info(time, name=name)
    nwin, _ = _window_samples(
        time, window_seconds, minimum=max(polyorder + 2, 3), name=name
    )
    if nwin <= polyorder:
        raise ValueError(
            f"{name}: effective window ({nwin} samples) must exceed polyorder={polyorder}"
        )

    if rel_irregularity <= 0.05:
        smoothed = savgol_filter(values, nwin, polyorder, mode="interp")
    else:
        uniform_time = np.arange(time[0], time[-1] + 0.5 * dt, dt)
        uniform_values = np.interp(uniform_time, time, values)
        nwin_uniform = max(int(round(window_seconds / dt)), polyorder + 2)
        if nwin_uniform % 2 == 0:
            nwin_uniform += 1
        nmax = len(uniform_time) if len(uniform_time) % 2 == 1 else len(uniform_time) - 1
        nwin_uniform = min(nwin_uniform, nmax)
        if nwin_uniform <= polyorder:
            raise ValueError(f"{name}: effective uniform-grid window is too short")
        uniform_smooth = savgol_filter(
            uniform_values, nwin_uniform, polyorder, mode="interp"
        )
        smoothed = np.interp(time, uniform_time, uniform_smooth)

    return np.asarray(smoothed, dtype=float), {
        "method": "savgol",
        "median_dt_s": dt,
        "window_samples": int(nwin),
        "timescale_seconds": float(window_seconds),
        "max_dt_irregularity_fraction": rel_irregularity,
        "polyorder": int(polyorder),
    }


def _butterworth_time_series(time, values, timescale_seconds, order=2, *, name="series"):
    """Return zero-phase Butterworth low-pass-smoothed data.

    ``timescale_seconds`` is interpreted as the cutoff period: variability with
    periods much shorter than this is increasingly attenuated.  Filtering is
    performed on a uniform median-dt grid and with ``sosfiltfilt`` so no phase
    lag is introduced into chamber cooling/warming features.
    """
    time = np.asarray(time, dtype=float)
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or len(time) != len(values):
        raise ValueError(f"{name}: time and values must be one-dimensional and equal length")
    if timescale_seconds <= 0.0:
        raise ValueError(f"{name}: Butterworth cutoff period must be positive")
    if order < 1:
        raise ValueError(f"{name}: Butterworth order must be >= 1")

    dt, rel_irregularity = _time_sampling_info(time, name=name)
    uniform_time = np.arange(time[0], time[-1] + 0.5 * dt, dt)
    uniform_values = np.interp(uniform_time, time, values)

    fs = 1.0 / dt
    nyquist = 0.5 * fs
    requested_cutoff_hz = 1.0 / float(timescale_seconds)
    cutoff_hz = min(requested_cutoff_hz, 0.90 * nyquist)

    sos = butter(int(order), cutoff_hz, btype="lowpass", fs=fs, output="sos")
    uniform_smooth = sosfiltfilt(sos, uniform_values)
    smoothed = np.interp(time, uniform_time, uniform_smooth)

    return np.asarray(smoothed, dtype=float), {
        "method": "butterworth",
        "median_dt_s": dt,
        "timescale_seconds": float(timescale_seconds),
        "cutoff_hz": float(cutoff_hz),
        "requested_cutoff_hz": float(requested_cutoff_hz),
        "order": int(order),
        "max_dt_irregularity_fraction": rel_irregularity,
    }


def _smoothing_candidates(time, values, timescale_seconds, *, name):
    """Return raw/despiked/Savitzky-Golay/Butterworth candidates for one series."""
    time = np.asarray(time, dtype=float)
    filled = np.asarray(values, dtype=float)

    do_despike = bool(_smoothing_cfg("CHAMBER_DESPIKE", True))
    despike_window = float(_smoothing_cfg("CHAMBER_DESPIKE_WINDOW", 10.0))
    despike_nsigma = float(_smoothing_cfg("CHAMBER_DESPIKE_NSIGMA", 4.0))
    polyorder = int(_smoothing_cfg("CHAMBER_SMOOTH_POLYORDER", 2))
    butter_order = int(_smoothing_cfg("CHAMBER_BUTTERWORTH_ORDER", 2))

    if do_despike:
        despiked, despike_info = _hampel_time_series(
            time, filled, despike_window, despike_nsigma, name=name
        )
    else:
        despiked = filled.copy()
        dt, _ = _time_sampling_info(time, name=name)
        despike_info = {
            "despike_window_seconds": 0.0,
            "despike_window_samples": 1,
            "despike_nsigma": despike_nsigma,
            "n_despiked": 0,
            "median_dt_s": dt,
        }

    savgol, savgol_info = _savgol_time_series(
        time, despiked, timescale_seconds, polyorder, name=name
    )
    butterworth, butter_info = _butterworth_time_series(
        time, despiked, timescale_seconds, butter_order, name=name
    )

    return {
        "raw": filled,
        "despiked": despiked,
        "savgol": savgol,
        "butterworth": butterworth,
        "info": {
            "despike": despike_info,
            "savgol": savgol_info,
            "butterworth": butter_info,
        },
    }


def _selected_smoothing_method():
    """Resolve the forcing method while retaining backward compatibility."""
    if not bool(_smoothing_cfg("SMOOTH_CHAMBER_FORCING", True)):
        return "raw"
    method = str(_smoothing_cfg("CHAMBER_SMOOTH_METHOD", "butterworth")).strip().lower()
    aliases = {
        "none": "raw",
        "off": "raw",
        "sgolay": "savgol",
        "savitzky-golay": "savgol",
        "butter": "butterworth",
    }
    method = aliases.get(method, method)
    if method not in {"raw", "savgol", "butterworth"}:
        raise ValueError(
            "CHAMBER_SMOOTH_METHOD must be 'raw', 'savgol' or 'butterworth'"
        )
    return method


def _record_selected_smoothing(block, key, candidates, method, timescale_seconds):
    """Store all candidates and make the selected one the working forcing series."""
    block[f"{key}_despiked"] = np.asarray(candidates["despiked"], dtype=float)
    block[f"{key}_savgol"] = np.asarray(candidates["savgol"], dtype=float)
    block[f"{key}_butterworth"] = np.asarray(candidates["butterworth"], dtype=float)
    block[key] = np.asarray(candidates[method], dtype=float).copy()

    filled = np.asarray(candidates["raw"], dtype=float)
    selected = np.asarray(block[key], dtype=float)
    correction = selected - filled
    info = {
        "selected_method": method,
        "timescale_seconds": float(timescale_seconds),
        "rms_correction": float(np.sqrt(np.mean(correction**2))),
        "max_abs_correction": float(np.max(np.abs(correction))),
        "mean_correction": float(np.mean(correction)),
        "n_despiked": int(candidates["info"]["despike"]["n_despiked"]),
        "despike_window_seconds": float(
            candidates["info"]["despike"]["despike_window_seconds"]
        ),
        "despike_nsigma": float(candidates["info"]["despike"]["despike_nsigma"]),
        "savgol": dict(candidates["info"]["savgol"]),
        "butterworth": dict(candidates["info"]["butterworth"]),
    }
    block["forcing_smoothing"][key] = info


def _smooth_chamber_forcing_block(block, *, experiment):
    """Prepare robust chamber P/T forcing and retain all smoothing alternatives.

    Gas temperature and pressure are smoothed directly.  By default the wall
    temperature is reconstructed as

        T_wall,filtered = T_gas,filtered + filtered(T_wall - T_gas)

    so filtering does not accidentally distort the wall-minus-gas temperature
    contrast that drives the chamber boundary-layer treatment.
    """
    time = np.asarray(block["Time"], dtype=float)
    method = _selected_smoothing_method()
    block.setdefault("forcing_smoothing", {})

    # Preserve raw instrument observations once and fill only missing values.
    for key, label in (
        ("Tgw mean", "gas temperature"),
        ("Tww mean", "wall temperature"),
        ("Pressure", "pressure"),
    ):
        raw_key = f"{key}_raw"
        if raw_key not in block:
            block[raw_key] = np.asarray(block[key], dtype=float).copy()
        filled = _fill_nan_linear(time, block[key], name=f"{experiment} {label}")
        block[f"{key}_filled"] = filled.copy()

    temp_window = float(_smoothing_cfg("CHAMBER_SMOOTH_TEMP_WINDOW", 30.0))
    wall_window = float(_smoothing_cfg("CHAMBER_SMOOTH_WALL_TEMP_WINDOW", 45.0))
    pressure_window = float(_smoothing_cfg("CHAMBER_SMOOTH_PRESSURE_WINDOW", 15.0))

    gas_candidates = _smoothing_candidates(
        time, block["Tgw mean_filled"], temp_window,
        name=f"{experiment} gas temperature",
    )
    pressure_candidates = _smoothing_candidates(
        time, block["Pressure_filled"], pressure_window,
        name=f"{experiment} pressure",
    )

    use_delta_t = bool(
        _smoothing_cfg("CHAMBER_SMOOTH_WALL_AS_DELTA_T", True)
    )
    if use_delta_t:
        delta_filled = (
            np.asarray(block["Tww mean_filled"], dtype=float)
            - np.asarray(block["Tgw mean_filled"], dtype=float)
        )
        delta_candidates = _smoothing_candidates(
            time, delta_filled, wall_window,
            name=f"{experiment} wall-minus-gas temperature",
        )
        wall_candidates = {
            "raw": np.asarray(block["Tww mean_filled"], dtype=float),
            "despiked": gas_candidates["despiked"] + delta_candidates["despiked"],
            "savgol": gas_candidates["savgol"] + delta_candidates["savgol"],
            "butterworth": (
                gas_candidates["butterworth"] + delta_candidates["butterworth"]
            ),
            "info": {
                "despike": delta_candidates["info"]["despike"],
                "savgol": delta_candidates["info"]["savgol"],
                "butterworth": delta_candidates["info"]["butterworth"],
            },
        }
        block["wall_minus_gas_raw"] = delta_filled
        block["wall_minus_gas_despiked"] = delta_candidates["despiked"]
        block["wall_minus_gas_savgol"] = delta_candidates["savgol"]
        block["wall_minus_gas_butterworth"] = delta_candidates["butterworth"]
    else:
        wall_candidates = _smoothing_candidates(
            time, block["Tww mean_filled"], wall_window,
            name=f"{experiment} wall temperature",
        )
        block["wall_minus_gas_raw"] = (
            np.asarray(block["Tww mean_filled"])
            - np.asarray(block["Tgw mean_filled"])
        )
        for candidate_method in ("despiked", "savgol", "butterworth"):
            block[f"wall_minus_gas_{candidate_method}"] = (
                np.asarray(wall_candidates[candidate_method])
                - np.asarray(gas_candidates[candidate_method])
            )

    _record_selected_smoothing(
        block, "Tgw mean", gas_candidates, method, temp_window
    )
    _record_selected_smoothing(
        block, "Tww mean", wall_candidates, method, wall_window
    )
    _record_selected_smoothing(
        block, "Pressure", pressure_candidates, method, pressure_window
    )

    block["wall_minus_gas"] = (
        np.asarray(block["Tww mean"], dtype=float)
        - np.asarray(block["Tgw mean"], dtype=float)
    )
    block["forcing_smoothing"]["wall_as_delta_t"] = use_delta_t
    block["forcing_smoothing"]["selected_method"] = method


def _series_derivative(time, values):
    """Return d(values)/dt on the native time grid."""
    time = np.asarray(time, dtype=float)
    values = np.asarray(values, dtype=float)
    return np.gradient(values, time)


def _write_forcing_smoothing_diagnostics(data):
    """Compare raw, despiked, Savitzky-Golay and Butterworth chamber forcing."""
    if not bool(_smoothing_cfg("SAVE_FORCING_SMOOTHING_DIAGNOSTICS", False)):
        return

    outdir = cfg.OUTPUT_ROOT / "forcing_smoothing"
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []

    method_colours = {
        "raw": "#7F7F7F",
        "despiked": "#CC79A7",
        "savgol": "#0072B2",
        "butterworth": "#009E73",
    }

    for key in readMeteoCPC.metStr:
        if key not in data:
            continue
        met = data[key]
        exp = key.split("-")[-1]
        time = np.asarray(met["Time"], dtype=float)
        selected = met.get("forcing_smoothing", {}).get("selected_method", "raw")

        field_specs = (
            ("Tgw mean", "Gas temperature", "degC"),
            ("Tww mean", "Wall temperature", "degC"),
            ("Pressure", "Pressure", "hPa"),
        )
        for field, label, units in field_specs:
            info = met.get("forcing_smoothing", {}).get(field, {})
            raw = np.asarray(met[f"{field}_filled"], dtype=float)
            for method in ("raw", "despiked", "savgol", "butterworth"):
                series_key = (
                    f"{field}_filled" if method == "raw" else f"{field}_{method}"
                )
                series = np.asarray(met[series_key], dtype=float)
                correction = series - raw
                deriv = _series_derivative(time, series)
                rows.append(
                    {
                        "experiment": exp,
                        "field": field,
                        "units": units,
                        "method": method,
                        "selected": method == selected,
                        "timescale_seconds": info.get("timescale_seconds", 0.0),
                        "n_despiked": info.get("n_despiked", 0),
                        "rms_correction": float(np.sqrt(np.mean(correction**2))),
                        "max_abs_correction": float(np.max(np.abs(correction))),
                        "derivative_rms_per_s": float(np.sqrt(np.mean(deriv**2))),
                        "derivative_max_abs_per_s": float(np.max(np.abs(deriv))),
                    }
                )

        fig, axes = plt.subplots(4, 2, figsize=(14, 14), sharex=True)
        ax = axes.ravel()

        # Gas temperature and its time derivative.
        for method in ("raw", "despiked", "savgol", "butterworth"):
            key_name = (
                "Tgw mean_filled" if method == "raw" else f"Tgw mean_{method}"
            )
            linewidth = 2.4 if method == selected else 1.3
            alpha = 1.0 if method == selected else 0.72
            label = f"{method}" + (" [FORCING]" if method == selected else "")
            y = np.asarray(met[key_name], dtype=float)
            ax[0].plot(
                time / 60.0, y, color=method_colours[method],
                linewidth=linewidth, alpha=alpha, label=label,
            )
            ax[1].plot(
                time / 60.0, _series_derivative(time, y),
                color=method_colours[method], linewidth=linewidth,
                alpha=alpha, label=label,
            )
        ax[0].set_title("Gas temperature smoothing comparison")
        ax[0].set_ylabel(r"$T_{gas}$ ($^\circ$C)")
        ax[1].set_title("Gas-temperature derivative")
        ax[1].set_ylabel(r"$dT_{gas}/dt$ (K s$^{-1}$)")

        # Wall temperature and the physically important wall-gas contrast.
        for method in ("raw", "despiked", "savgol", "butterworth"):
            wall_key = (
                "Tww mean_filled" if method == "raw" else f"Tww mean_{method}"
            )
            delta_key = (
                "wall_minus_gas_raw"
                if method == "raw"
                else f"wall_minus_gas_{method}"
            )
            linewidth = 2.4 if method == selected else 1.3
            alpha = 1.0 if method == selected else 0.72
            label = f"{method}" + (" [FORCING]" if method == selected else "")
            ax[2].plot(
                time / 60.0, np.asarray(met[wall_key]),
                color=method_colours[method], linewidth=linewidth,
                alpha=alpha, label=label,
            )
            ax[3].plot(
                time / 60.0, np.asarray(met[delta_key]),
                color=method_colours[method], linewidth=linewidth,
                alpha=alpha, label=label,
            )
        ax[2].set_title("Wall temperature smoothing comparison")
        ax[2].set_ylabel(r"$T_{wall}$ ($^\circ$C)")
        ax[3].set_title(r"Wall–gas contrast (important for BL forcing)")
        ax[3].set_ylabel(r"$T_{wall}-T_{gas}$ (K)")
        ax[3].axhline(0.0, color="0.5", linewidth=0.8)

        # Pressure and derivative.
        for method in ("raw", "despiked", "savgol", "butterworth"):
            pkey = "Pressure_filled" if method == "raw" else f"Pressure_{method}"
            linewidth = 2.4 if method == selected else 1.3
            alpha = 1.0 if method == selected else 0.72
            label = f"{method}" + (" [FORCING]" if method == selected else "")
            p = np.asarray(met[pkey], dtype=float)
            ax[4].plot(
                time / 60.0, p, color=method_colours[method],
                linewidth=linewidth, alpha=alpha, label=label,
            )
            ax[5].plot(
                time / 60.0, _series_derivative(time, p),
                color=method_colours[method], linewidth=linewidth,
                alpha=alpha, label=label,
            )
        ax[4].set_title("Pressure smoothing comparison")
        ax[4].set_ylabel("Pressure (hPa)")
        ax[5].set_title("Pressure derivative")
        ax[5].set_ylabel(r"$dP/dt$ (hPa s$^{-1}$)")

        # Residuals show exactly what each smoother removes.
        gas_raw = np.asarray(met["Tgw mean_filled"], dtype=float)
        p_raw = np.asarray(met["Pressure_filled"], dtype=float)
        for method in ("despiked", "savgol", "butterworth"):
            ax[6].plot(
                time / 60.0,
                np.asarray(met[f"Tgw mean_{method}"]) - gas_raw,
                color=method_colours[method], label=method,
            )
            ax[7].plot(
                time / 60.0,
                np.asarray(met[f"Pressure_{method}"]) - p_raw,
                color=method_colours[method], label=method,
            )
        ax[6].set_title("Gas-temperature correction")
        ax[6].set_ylabel(r"Filtered - raw $T$ (K)")
        ax[7].set_title("Pressure correction")
        ax[7].set_ylabel("Filtered - raw P (hPa)")

        for a in ax:
            a.grid(alpha=0.3)
            a.legend(fontsize=8, loc="best")
            a.set_xlabel("Experiment time (min)")

        gas_info = met["forcing_smoothing"]["Tgw mean"]
        wall_info = met["forcing_smoothing"]["Tww mean"]
        p_info = met["forcing_smoothing"]["Pressure"]
        delta_mode = met["forcing_smoothing"].get("wall_as_delta_t", False)
        fig.suptitle(
            f"{exp}: chamber-forcing smoothing\n"
            f"selected={selected}; T={gas_info['timescale_seconds']:.0f}s, "
            f"wall/ΔT={wall_info['timescale_seconds']:.0f}s, "
            f"P={p_info['timescale_seconds']:.0f}s; "
            f"wall-as-ΔT={delta_mode}",
            fontsize=13,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        fig.savefig(outdir / f"forcing-smoothing-comparison-{exp}.png", dpi=180)
        plt.close(fig)

    if rows:
        csv_path = outdir / "forcing-smoothing-method-comparison.csv"
        with csv_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"Forcing-smoothing diagnostics written to {outdir}")


# Historical batch group to run when executing this file directly.
THIS_RUN = 6
RUN_MODEL = True
DO_ANALYSIS = True
DO_PLOT = True
READ_DATA = True


def _svp_liq(temperature_k):
    return np.asarray(svp.svp(np.asarray(temperature_k), "buck2", "liq"), dtype=float)


def _mixing_ratio_from_vapour_pressure(e_pa, p_pa):
    return EPSILON * e_pa / (p_pa - e_pa)


def _fill_nan_linear(time, values, *, name="series"):
    """Fill isolated/missing observations by linear interpolation.

    ``np.interp`` also extends the nearest valid value across missing end
    points.  A completely missing series is rejected because silently filling
    it would create an artificial chamber forcing.
    """
    time = np.asarray(time, dtype=float)
    values = np.asarray(values, dtype=float).copy()
    good = np.isfinite(time) & np.isfinite(values)
    if not np.any(good):
        raise ValueError(f"No finite data available for {name}")
    values[~np.isfinite(values)] = np.interp(
        time[~np.isfinite(values)], time[good], values[good]
    )
    return values


def _interp_at(time, values, t, *, name="series"):
    time = np.asarray(time, dtype=float)
    values = np.asarray(values, dtype=float)
    good = np.isfinite(time) & np.isfinite(values)
    if not np.any(good):
        raise ValueError(f"No finite data available for {name}")
    if t < np.nanmin(time[good]) or t > np.nanmax(time[good]):
        raise ValueError(
            f"Requested {name} at t={t:g} s outside observed range "
            f"[{np.nanmin(time[good]):g}, {np.nanmax(time[good]):g}] s"
        )
    return float(np.interp(t, time[good], values[good]))


def load_all_data(experiments=None):
    """Read only the requested experiments, or every configured experiment.

    ``experiments`` contains canonical experiment names such as ``Exp021``.
    This is important for single-run mode: a missing file for some unrelated
    experiment (for example Exp020) must not prevent ``--experiment Exp021``
    from running.
    """
    data = {}
    wanted = None if experiments is None else {
        _normalise_experiment_name(exp) for exp in experiments
    }

    def _wanted(key):
        return wanted is None or key.split("-")[-1] in wanted

    for i, key in enumerate(readMeteoCPC.metStr):
        if not _wanted(key):
            continue
        block = readMeteoCPC.readData(readThis=i, metStr=key)[key]

        # Preserve and smooth only the variables used as fast chamber forcing.
        # Dew point is deliberately *not* smoothed here: it is a separate water
        # measurement/diagnostic whose uncertainty should not be hidden by the
        # thermal-forcing filter.
        _smooth_chamber_forcing_block(block, experiment=key)
        block["TDew"] = _fill_nan_linear(
            block["Time"], block["TDew"], name=f"{key} dew point"
        )
        data[key] = block

    for i, key in enumerate(readOPC_Merged.opcStr):
        if not _wanted(key):
            continue
        data[key] = readOPC_Merged.readData(readThis=i, opcStr=key)[key]

    for i, key in enumerate(readPNSD_Mrg_new.npsdStr):
        if not _wanted(key):
            continue
        data[key] = readPNSD_Mrg_new.readData(readThis=i, npsdStr=key)[key]

    if wanted is not None:
        for exp in wanted:
            required = (
                f"MeteoCPC-{exp}",
                f"MergedOPC-{exp}",
                f"InitialPNSD-{exp}",
            )
            missing = [key for key in required if key not in data]
            if missing:
                raise KeyError(
                    f"{exp}: requested experiment is missing configured reader "
                    f"entries/data for: {', '.join(missing)}"
                )

    _derive_water_series(data)
    return data


def _derive_water_series(data):
    """Add qv, qsat, RH and a configurable qtot diagnostic to meteo blocks.

    Dew point determines vapour mixing ratio.  For the recommended
    ``vapour_plus_opc`` mode, OPC liquid water is interpolated onto the meteo
    time grid and converted from g m-3 to kg kg-1 dry air before forming qtot.
    This matches the interpretation used elsewhere in this repository and
    avoids labelling dew-point vapour alone as total water.
    """
    for key, met in list(data.items()):
        if not key.startswith("MeteoCPC-"):
            continue
        exp = key.split("-")[-1]
        p_pa = met["Pressure"] * 100.0
        t_k = met["Tgw mean"] + 273.15
        e = _svp_liq(met["TDew"] + 273.15)
        es = _svp_liq(t_k)

        qv = _mixing_ratio_from_vapour_pressure(e, p_pa)
        qsat = _mixing_ratio_from_vapour_pressure(es, p_pa)
        met["qv_dewpoint"] = qv
        met["qsat"] = qsat
        met["rh_dewpoint"] = qv / qsat

        if cfg.QTOT_DATA_MODE == "dewpoint_only":
            qtot = qv.copy()
        elif cfg.QTOT_DATA_MODE == "vapour_plus_opc":
            opc_key = f"MergedOPC-{exp}"
            if opc_key not in data:
                raise KeyError(f"Missing OPC data required for qtot: {opc_key}")
            opc = data[opc_key]
            lwc = np.asarray(opc["lwc"], dtype=float)
            finite = np.isfinite(opc["Time"]) & np.isfinite(lwc)
            if not np.any(finite):
                raise ValueError(f"No finite OPC LWC for {exp}")
            lwc_interp = np.interp(met["Time"], opc["Time"][finite], lwc[finite])
            rho_d = p_pa / (R_D * t_k)
            ql = lwc_interp * 1.0e-3 / rho_d  # g m-3 -> kg m-3 -> kg kg-1
            qtot = qv + ql
            met["ql_opc"] = ql
        else:
            raise ValueError(f"Unknown QTOT_DATA_MODE={cfg.QTOT_DATA_MODE!r}")

        met["qtot_forcing"] = _fill_nan_linear(met["Time"], qtot, name=f"{exp} qtot")


def _normalise_experiment_name(experiment):
    """Return canonical ExpNNN form for an experiment identifier."""
    exp = str(experiment)
    if exp.startswith("Exp"):
        suffix = exp[3:]
        if suffix.isdigit():
            return f"Exp{int(suffix):03d}"
        return exp
    return f"Exp{int(exp):03d}"


def _model_saturation_time(exp, observed_cloud_time, runtime_time, *, cli_override_s=None):
    """Return the model saturation target without changing observational timing.

    The historical target is ``experiment_metadata.CLOUD_ONSET[exp]``.  A
    global config shift may move all experiments and a per-experiment absolute
    config override may replace it.  For a single run, an explicit CLI
    override has highest precedence.
    """
    target = float(observed_cloud_time) + float(cfg.MODEL_SATURATION_TIME_SHIFT_S)
    overrides = dict(getattr(cfg, "MODEL_SATURATION_TIME_OVERRIDES_S", {}))
    if exp in overrides:
        target = float(overrides[exp])
    if cli_override_s is not None:
        target = float(cli_override_s)

    tmin = float(np.nanmin(runtime_time))
    tmax = float(np.nanmax(runtime_time))
    if target < tmin or target > tmax:
        raise ValueError(
            f"{exp}: model saturation target {target:.3f} s is outside the "
            f"measured forcing interval [{tmin:.3f}, {tmax:.3f}] s"
        )
    return target


def _cloud_onset_for_experiment(exp, data):
    """Return curated cloud onset, or estimate it from OPC Nd when absent."""
    if exp in meta.CLOUD_ONSET:
        return float(meta.CLOUD_ONSET[exp])
    opc_key = f"MergedOPC-{exp}"
    if opc_key not in data:
        raise KeyError(f"No cloud-onset metadata or OPC data for {exp}")
    opc = data[opc_key]
    time = np.asarray(opc["Time"], dtype=float)
    nd = np.asarray(opc["ndrop"], dtype=float)
    good = np.isfinite(time) & np.isfinite(nd)
    if not np.any(good):
        raise ValueError(f"Cannot estimate cloud onset for {exp}: no finite OPC Nd")
    peak = float(np.nanmax(nd[good]))
    threshold = max(1.0, 0.05 * peak)
    active = good & (nd >= threshold)
    if not np.any(active):
        raise ValueError(f"Cannot estimate cloud onset for {exp}: no OPC cloud signal")
    onset = float(time[np.flatnonzero(active)[0]])
    print(f"{exp}: no curated cloud onset; using OPC-derived onset {onset/60.0:.2f} min")
    return onset


def build_initial_state(data, *, saturation_time_overrides_s=None):
    """Calculate initial thermodynamics and aerosol normalisation metadata.

    ``saturation_time_overrides_s`` is intended for one-off single-run CLI
    sensitivity tests.  It changes only the model saturation target used to
    infer the initial vapour amount; observed cloud-onset metadata and plotting
    windows remain fixed.
    """
    states = {}
    saturation_time_overrides_s = saturation_time_overrides_s or {}
    for met_key in readMeteoCPC.metStr:
        exp = met_key.split("-")[-1]
        if met_key not in data:
            continue
        observed_cloud_time = _cloud_onset_for_experiment(exp, data)
        met = data[met_key]
        saturation_time = _model_saturation_time(
            exp,
            observed_cloud_time,
            met["Time"],
            cli_override_s=saturation_time_overrides_s.get(exp),
        )
        t0 = _interp_at(met["Time"], met["Tgw mean"], 0.0, name=f"{exp} Tgas") + 273.15
        p0 = _interp_at(met["Time"], met["Pressure"], 0.0, name=f"{exp} pressure") * 100.0
        tsat = _interp_at(
            met["Time"], met["Tgw mean"], saturation_time,
            name=f"{exp} saturation-target T",
        ) + 273.15
        psat = _interp_at(
            met["Time"], met["Pressure"], saturation_time,
            name=f"{exp} saturation-target p",
        ) * 100.0

        if cfg.INITIAL_RH_METHOD == "cloud_onset":
            # Historical initialisation: choose qv so that, in the absence of
            # pre-cloud water exchange, the measured P/T trajectory reaches
            # liquid saturation at the prescribed cloud-onset time.
            es_target = float(_svp_liq([tsat])[0])
            qv0 = _mixing_ratio_from_vapour_pressure(es_target, psat)
            es0 = float(_svp_liq([t0])[0])
            rh0 = qv0 / _mixing_ratio_from_vapour_pressure(es0, p0)
        elif cfg.INITIAL_RH_METHOD == "dewpoint":
            qv0 = _interp_at(met["Time"], met["qv_dewpoint"], 0.0, name=f"{exp} qv")
            es0 = float(_svp_liq([t0])[0])
            rh0 = qv0 / _mixing_ratio_from_vapour_pressure(es0, p0)
        else:
            raise ValueError(f"Unknown INITIAL_RH_METHOD={cfg.INITIAL_RH_METHOD!r}")

        if cfg.AEROSOL_INIT_TIME == "t0":
            aerosol_time = 0.0
        elif cfg.AEROSOL_INIT_TIME == "cloud_onset":
            # Aerosol normalisation remains tied to the observed cloud onset,
            # not to a model saturation-time sensitivity setting.
            aerosol_time = observed_cloud_time
        else:
            raise ValueError(f"Unknown AEROSOL_INIT_TIME={cfg.AEROSOL_INIT_TIME!r}")

        cpc_cm3 = _interp_at(met["Time"], met["CPC_TotBot"], aerosol_time, name=f"{exp} CPC")
        ta = _interp_at(met["Time"], met["Tgw mean"], aerosol_time, name=f"{exp} aerosol T") + 273.15
        pa = _interp_at(met["Time"], met["Pressure"], aerosol_time, name=f"{exp} aerosol p") * 100.0
        rho_d = pa / (R_D * ta)
        aerosol_number_mixing_ratio = cpc_cm3 * 1.0e6 / rho_d  # # kg-1 dry air

        states[exp] = {
            "cloud_time": observed_cloud_time,
            "saturation_time": saturation_time,
            "initT": t0,
            "initP": p0,
            "cloudT": tsat,
            "cloudP": psat,
            "initRH": rh0,
            "aerConc": aerosol_number_mixing_ratio,
            "aerosol_norm_time": aerosol_time,
        }

    if cfg.AEROSOL_INIT_TIME == "cloud_onset" and (
        cfg.CHAMBER_FAN_LOSS or cfg.CHAMBER_WALL_LOSS or cfg.FALLOUT_FLAG
    ):
        print(
            "WARNING: aerosol is normalised at cloud onset while particle-loss "
            "processes are enabled.  Pre-cloud aerosol loss will therefore be "
            "double represented.  AEROSOL_INIT_TIME='t0' is recommended."
        )
    if cfg.INITIAL_RH_METHOD == "cloud_onset" and cfg.CHAMBER_BL_MIX:
        print(
            "WARNING: INITIAL_RH_METHOD='cloud_onset' assumes no pre-cloud water "
            "exchange, but chamber BL mixing is enabled.  Cloud onset is then a "
            "diagnostic prediction rather than a guaranteed constraint."
        )
    return states


def _format_fortran_array(name, values, n, values_per_line=8):
    values = np.asarray(values, dtype=float)
    if len(values) != n:
        raise ValueError(f"{name}: expected {n} values, got {len(values)}")
    pieces = [f"{v:.12g}" for v in values]
    lines = []
    for start in range(0, n, values_per_line):
        chunk = ", ".join(pieces[start : start + values_per_line]) + ","
        if start == 0:
            lines.append(f"    {name}(1:{n}) = {chunk}")
        else:
            lines.append(f"        {chunk}")
    return "\n".join(lines)


def _set_dynamic_array_slice(text, name, values):
    """Replace a 1-D namelist array assignment regardless of old upper bound.

    For example, if the template contains ``inp_temp(1:16)`` and 50 values are
    supplied, this finds the existing slice and asks ``set_array`` to replace
    that exact template object with ``inp_temp(1:50)``.
    """
    values = np.asarray(values, dtype=float).reshape(-1)
    if values.size < 1:
        raise ValueError(f"{name} must contain at least one value")

    pattern = re.compile(
        rf"(?im)^[ \t]*({re.escape(name)}\s*\(\s*1\s*:\s*\d+\s*\))\s*="
    )
    match = pattern.search(text)
    if match is None:
        raise ValueError(
            f"Could not find an existing {name}(1:N) assignment in the BMM template"
        )

    old_object = match.group(1)
    updated = set_array(text, old_object, values)

    # set_array replaces values but may retain the old object text/bounds,
    # depending on namelist_utils version. Ensure the LHS bound matches the
    # actual number of supplied classes.
    updated = re.sub(
        rf"(?im)^([ \t]*){re.escape(old_object)}\s*=",
        lambda m: f"{m.group(1)}{name}(1:{values.size}) =",
        updated,
        count=1,
    )
    return updated


def chamber_spec_body(met, *, require_wall_temperature=False):
    """Create the variable-length &chamber_spec body for one experiment.

    ``Pressure``, ``Tgw mean`` and ``Tww mean`` are the working series from
    :func:`load_all_data`; when chamber smoothing is enabled these are the
    smoothed trajectories, so raw instrument jitter is never exported by
    accident.
    """
    n = len(met["Time"])
    blocks = [
        "    ! Chamber observations on a common experiment time grid.",
        _format_fortran_array("time_chamber", met["Time"], n),
        _format_fortran_array("press_chamber", met["Pressure"] * 100.0, n),
        _format_fortran_array("temp_chamber", met["Tgw mean"] + 273.15, n),
    ]
    if cfg.WRITE_QTOT_DATA:
        blocks.extend(
            [
                "    ! Diagnostic/optional forcing total-water series.",
                _format_fortran_array("qtot_chamber", met["qtot_forcing"], n),
            ]
        )
    if require_wall_temperature:
        blocks.extend(
            [
                "    ! Measured wall temperature (Tww_mean) for wall-coupled BL processing.",
                _format_fortran_array("wall_temp_chamber", met["Tww mean"] + 273.15, n),
            ]
        )
    return "\n".join(blocks)


def _aerosol_mode_arrays(exp, state, data):
    """Return BMM aerosol mode arrays derived from the fitted measured PSD."""
    key = f"InitialPNSD-{exp}"
    index = readPNSD_Mrg_new.npsdStr.index(key)
    psd = data[key]

    comp1 = readPNSD_Mrg_new.comp[index][0]
    n1_fit = np.asarray(psd[f"Nfit_{comp1}_dNdlogDve_cc"], dtype=float)
    d1_um = np.asarray(psd[f"dfit_{comp1}_dNdlogDve_cc"], dtype=float)
    sig1 = np.asarray(psd[f"lnsigfit_{comp1}_dNdlogDve_cc"], dtype=float)

    n2_fit = np.zeros(3)
    d2_um = np.full(3, 0.1)   # fitting/reader diameter unit is micrometres
    sig2 = np.full(3, 0.5)
    if psd["num1"] > 1:
        comp2 = readPNSD_Mrg_new.comp[index][1]
        n2_fit = np.asarray(psd[f"Nfit_{comp2}_dNdlogDve_cc"], dtype=float)
        d2_um = np.asarray(psd[f"dfit_{comp2}_dNdlogDve_cc"], dtype=float)
        sig2 = np.asarray(psd[f"lnsigfit_{comp2}_dNdlogDve_cc"], dtype=float)

    total_fit = np.sum(n1_fit) + np.sum(n2_fit)
    if total_fit <= 0.0:
        raise ValueError(f"Fitted aerosol concentration is zero for {exp}")

    # Scale fitted component fractions to the observed CPC total number mixing
    # ratio.  The result is directly in BMM number-per-kg-dry-air units.
    scale = state["aerConc"] / total_fit
    n1 = np.zeros(3)
    n2 = np.zeros(3)
    n1[: len(n1_fit)] = scale * n1_fit
    n2[: len(n2_fit)] = scale * n2_fit

    # Keep tiny nonzero placeholders only where the BMM template historically
    # expected a third submode.  They are physically negligible but avoid an
    # exactly empty mode in workflows that assume three internal submodes.
    n1[len(n1_fit) :] = 1.0e3
    if psd["num1"] <= 1:
        n2[:] = [0.0, 0.0, 1.0e3]
    elif len(n2_fit) < 3:
        n2[len(n2_fit) :] = 1.0e3

    d1_out_um = np.full(3, 0.1)
    sig1_out = np.full(3, 0.5)
    d1_out_um[: len(d1_um)] = d1_um
    sig1_out[: len(sig1)] = sig1

    d2_out_um = np.full(3, 0.1)
    sig2_out = np.full(3, 0.5)
    if psd["num1"] > 1:
        d2_out_um[: len(d2_um)] = d2_um
        sig2_out[: len(sig2)] = sig2

    density = np.asarray(readPNSD_Mrg_new.density[index], dtype=float)
    kappa = np.asarray(readPNSD_Mrg_new.kappa[index], dtype=float)
    density4 = np.array([1770.0, 1770.0, 1770.0, 2160.0])
    kappa4 = np.array([0.1, 0.1, 0.1, 0.1])
    density4[: len(density)] = density
    kappa4[: len(kappa)] = kappa

    afhh4 = np.zeros(4, dtype=float)
    bfhh4 = np.zeros(4, dtype=float)
    inp_category4 = ["none"] * 4
    components = list(readPNSD_Mrg_new.comp[index])
    for j, component in enumerate(components):
        if component in cfg.DUST_FHH:
            # Dust adsorption is represented with FHH, not a negative-kappa
            # sentinel.  Its IASD identity is passed separately to BMM.
            kappa4[j] = 0.0
            afhh4[j] = float(cfg.DUST_FHH[component]["A"])
            bfhh4[j] = float(cfg.DUST_FHH[component]["B"])
            inp_category4[j] = cfg.DUST_INP_CATEGORY[component]

    return {
        "n1": n1,
        "d1_m": d1_out_um * 1.0e-6,
        "sig1": sig1_out,
        "n2": n2,
        "d2_m": d2_out_um * 1.0e-6,
        "sig2": sig2_out,
        "density": density4,
        "kappa": kappa4,
        "afhh": afhh4,
        "bfhh": bfhh4,
        "inp_category": inp_category4,
    }


def make_namelist(exp, state, data, output_file, *, winit=1.3, scefile=None,
                  bl_tau_s=None, bl_temp_offset_k=None,
                  bl_wall_water_mode=None, wall_liquid_water_init_kg=None,
                  wall_ice_water_init_kg=None, wall_water_efficiency=None,
                  wall_vapour_transfer_velocity_ms=None,
                  bl_evap_size_exp=None):
    """Return a complete BMM namelist for one iSKYLAB experiment.

    If ``scefile`` is supplied, point the main BMM namelist at that SCE
    namelist.  This lets per-run SCE sensitivities use a temporary copy
    without modifying the repository's ``sce/namelist.in``.
    """
    text = read_text(cfg.BMM_TEMPLATE)
    met = data[f"MeteoCPC-{exp}"]
    n = len(met["Time"])
    runtime = min(float(np.nanmax(met["Time"])), 30.0 * 60.0)

    # General initial/run controls.
    text = set_value(text, "outputfile", str(output_file))
    text = set_value(text, "runtime", runtime)
    text = set_value(text, "tinit", state["initT"])
    text = set_value(text, "pinit", state["initP"])
    text = set_value(text, "rhinit", state["initRH"])
    text = set_value(text, "winit", float(winit))
    if scefile is not None:
        text = set_value(text, "scefile", str(scefile))
    is_dust = exp in getattr(meta, "DUST_EXPERIMENTS", set())
    if is_dust and cfg.DUST_ENABLE_ICE:
        text = set_value(text, "ice_flag", 1)
        text = set_literal_array(
            text, "ice_nucleation_mech(1:4)", cfg.DUST_ICE_NUCLEATION_MECH
        )
    text = set_value(text, "fallout_flag", bool(cfg.FALLOUT_FLAG))
    text = set_value(text, "residence_depth", float(cfg.RESIDENCE_DEPTH))

    # Chamber observation/physics controls.
    force_p = cfg.FORCE_PRESSURE and not cfg.SYNTHETIC_UPDRAFT
    force_t = cfg.FORCE_TEMPERATURE and not cfg.SYNTHETIC_UPDRAFT
    force_q = cfg.FORCE_QTOT and not cfg.SYNTHETIC_UPDRAFT
    bl_mix = 0 if cfg.SYNTHETIC_UPDRAFT else cfg.CHAMBER_BL_MIX
    wall_water_mode = int(
        getattr(cfg, "CHAMBER_BL_WALL_WATER_MODE", 0)
        if bl_wall_water_mode is None else bl_wall_water_mode
    )
    wall_liq_init_kg = float(
        getattr(cfg, "CHAMBER_WALL_LIQUID_WATER_INIT_KG", 0.0)
        if wall_liquid_water_init_kg is None else wall_liquid_water_init_kg
    )
    wall_ice_init_kg = float(
        getattr(cfg, "CHAMBER_WALL_ICE_WATER_INIT_KG", 0.0)
        if wall_ice_water_init_kg is None else wall_ice_water_init_kg
    )
    wall_eff = float(
        getattr(cfg, "CHAMBER_WALL_WATER_EFFICIENCY", 1.0)
        if wall_water_efficiency is None else wall_water_efficiency
    )
    wall_km = float(
        getattr(cfg, "CHAMBER_WALL_VAPOUR_TRANSFER_VELOCITY", 1.0e-3)
        if wall_vapour_transfer_velocity_ms is None
        else wall_vapour_transfer_velocity_ms
    )
    evap_size_exp = float(
        getattr(
            cfg,
            "CHAMBER_BL_EVAP_SIZE_EXP",
            getattr(cfg, "CHAMBER_BL_INHOM_SIZE_EXP", 2.0),
        )
        if bl_evap_size_exp is None else bl_evap_size_exp
    )
    fan_loss = 0 if cfg.SYNTHETIC_UPDRAFT else cfg.CHAMBER_FAN_LOSS
    wall_loss = 0 if cfg.SYNTHETIC_UPDRAFT else cfg.CHAMBER_WALL_LOSS

    if force_q and not cfg.WRITE_QTOT_DATA:
        raise ValueError("FORCE_QTOT=True requires WRITE_QTOT_DATA=True")
    if bl_mix not in (0, 1):
        raise ValueError("CHAMBER_BL_MIX must be 0 (off) or 1 (on)")
    if wall_water_mode not in (0, 1, 2):
        raise ValueError("chamber BL wall-water mode must be 0, 1, or 2")
    if wall_liq_init_kg < 0.0 or wall_ice_init_kg < 0.0:
        raise ValueError("initial wall liquid/ice reservoirs must be non-negative kg")
    if not 0.0 <= wall_eff <= 1.0:
        raise ValueError("wall water efficiency must be between 0 and 1")
    if wall_km < 0.0:
        raise ValueError("wall vapour transfer velocity must be >= 0 m s-1")
    if evap_size_exp < 0.0:
        raise ValueError("chamber BL evaporation size exponent must be >= 0")
    need_wall_temperature = bool(
        bl_mix and (cfg.CHAMBER_BL_ALPHA_T > 0.0 or wall_water_mode >= 1)
    )
    if need_wall_temperature and "Tww mean" not in met:
        raise ValueError("BL thermal/wall-water processing requires measured Tww_mean")

    chamber_values = {
        "n_levels_c": n,
        "chamber_force_pressure": force_p,
        "chamber_force_temperature": force_t,
        "chamber_force_qtot": force_q,
        "chamber_bl_mix": bl_mix,
        "chamber_bl_tau": (cfg.CHAMBER_BL_TAU if bl_tau_s is None else float(bl_tau_s)),
        "chamber_bl_alpha_t": cfg.CHAMBER_BL_ALPHA_T,
        "chamber_bl_evap_mode": cfg.CHAMBER_BL_EVAP_MODE,
        "chamber_bl_evap_size_exp": evap_size_exp,
        # Keep the deprecated Fortran alias neutral.  This prevents an older
        # template value from overriding the new common exponent.
        "chamber_bl_inhom_size_exp": -1.0,
        "chamber_bl_temp_offset": (
            cfg.CHAMBER_BL_TEMP_OFFSET
            if bl_temp_offset_k is None else float(bl_temp_offset_k)
        ),
        "chamber_bl_wall_water_mode": wall_water_mode,
        "chamber_wall_water_efficiency": wall_eff,
        "chamber_wall_vapour_transfer_velocity": wall_km,
        "chamber_wall_liquid_water_init": wall_liq_init_kg,
        "chamber_wall_ice_water_init": wall_ice_init_kg,
        "chamber_fan_loss": fan_loss,
        "chamber_fan_loss_kmax": cfg.CHAMBER_FAN_LOSS_KMAX,
        "chamber_fan_loss_d50_ref": cfg.CHAMBER_FAN_LOSS_D50_REF,
        "chamber_fan_loss_exp": cfg.CHAMBER_FAN_LOSS_EXP,
        "chamber_fan_rpm": cfg.CHAMBER_FAN_RPM,
        "chamber_fan_rpm_ref": cfg.CHAMBER_FAN_RPM_REF,
        "chamber_wall_loss": wall_loss,
        "chamber_wall_ustar": cfg.CHAMBER_WALL_USTAR,
        "chamber_diameter": cfg.CHAMBER_DIAMETER,
        "chamber_height": cfg.CHAMBER_HEIGHT,
    }
    # The four wall-water controls post-date some iSKYLAB BMM templates.
    # Migrate only these known schema additions when an older template is
    # encountered; all established variables remain strict and must already
    # exist.  Anchored insertion keeps the additions inside &chamber_options.
    wall_water_schema = {
        "chamber_bl_evap_size_exp": (
            "chamber_bl_evap_mode",
            "Common evaporation size exponent p: homogeneous changes size, inhomogeneous changes number; p=2 is D2-based.",
        ),
        "chamber_bl_inhom_size_exp": (
            "chamber_bl_evap_size_exp",
            "Deprecated compatibility alias; keep at -1 so it cannot override the new common exponent.",
        ),
        "chamber_bl_wall_water_mode": (
            "chamber_bl_temp_offset",
            "Wall-vapour closure: 0=legacy; 1=reservoir relaxation; 2=physical finite-rate mass transfer.",
        ),
        "chamber_wall_water_efficiency": (
            "chamber_bl_wall_water_mode",
            "Fractional vapour equilibration during one BL wall encounter [0..1].",
        ),
        "chamber_wall_vapour_transfer_velocity": (
            "chamber_wall_water_efficiency",
            "Physical wall-vapour mass-transfer velocity [m s-1], used by wall-water mode 2.",
        ),
        "chamber_wall_liquid_water_init": (
            "chamber_wall_vapour_transfer_velocity",
            "Initial physical liquid-water mass stored on the wall [kg]; default dry wall = 0.",
        ),
        "chamber_wall_ice_water_init": (
            "chamber_wall_liquid_water_init",
            "Initial physical ice/frost mass stored on the wall [kg]; default dry wall = 0.",
        ),
    }
    for name, value in chamber_values.items():
        if name in wall_water_schema:
            after, comment = wall_water_schema[name]
            text = set_or_insert_value(
                text, name, value, after=after, group="chamber_options",
                comment=comment,
            )
        else:
            text = set_value(text, name, value)
    text = replace_group(
        text, "chamber_spec",
        chamber_spec_body(met, require_wall_temperature=need_wall_temperature),
    )

    # Measured/fitted aerosol initial conditions.
    aer = _aerosol_mode_arrays(exp, state, data)

    inp_temp_c = np.asarray(cfg.INP_TEMP_C, dtype=float).reshape(-1)
    if inp_temp_c.size < 1:
        raise ValueError("INP_TEMP_C must contain at least one temperature threshold")
    if not np.all(np.isfinite(inp_temp_c)):
        raise ValueError("INP_TEMP_C contains non-finite values")
    if np.any(np.diff(inp_temp_c) >= 0.0):
        raise ValueError(
            "INP_TEMP_C must be strictly ordered from warm to cold "
            "(e.g. -18, -18.3, ..., -33 degC)"
        )

    for name, values in [
        ("n_aer1(1:3,1:1)", aer["n1"]),
        ("d_aer1(1:3,1:1)", aer["d1_m"]),
        ("sig_aer1(1:3,1:1)", aer["sig1"]),
        ("n_aer1(1:3,2:2)", aer["n2"]),
        ("d_aer1(1:3,2:2)", aer["d2_m"]),
        ("sig_aer1(1:3,2:2)", aer["sig2"]),
        ("density_core1(1:4)", aer["density"]),
        ("kappa_core1(1:4)", aer["kappa"]),
        ("afhh_core1(1:4)", aer["afhh"]),
        ("bfhh_core1(1:4)", aer["bfhh"]),
    ]:
        text = set_array(text, name, values)

    text = _set_dynamic_array_slice(text, "inp_temp", inp_temp_c)
    text = set_literal_array(text, "inp_category(1:4)", aer["inp_category"])
    text = set_value(text, "n_inp_classes", int(inp_temp_c.size))

    return text, aer


def _make_temp_sce_namelist(sce_bins, tmpdir):
    """Create a temporary SCE namelist with an optional ``n_binsc`` override.

    The repository copy of ``sce/namelist.in`` is never modified.  When
    ``sce_bins`` is supplied, the base SCE namelist is copied into ``tmpdir``,
    ``n_binsc`` is changed in the copy, and the returned path should be written
    into the main BMM namelist via its ``scefile`` variable.

    If ``sce_bins`` is None, return None so the main BMM template keeps its
    existing ``scefile`` setting unchanged.
    """
    if sce_bins is None:
        return None

    sce_bins = int(sce_bins)
    if sce_bins < 0:
        raise ValueError("sce_bins must be >= 0")

    source = Path(
        getattr(cfg, "SCE_NAMELIST", cfg.BMM_MODEL_FOLDER / "sce" / "namelist.in")
    )
    if not source.exists():
        raise FileNotFoundError(f"BMM SCE namelist not found: {source}")

    tmpdir = Path(tmpdir)
    tmpdir.mkdir(parents=True, exist_ok=True)
    target = tmpdir / f"sce-namelist-nbinsc{sce_bins}.in"

    text = read_text(source)
    text = set_value(text, "n_binsc", sce_bins)
    write_text(target, text)

    print(
        f"Temporary SCE namelist: n_binsc={sce_bins}; "
        f"source={source}; run copy={target}"
    )
    return target


def run_batch(batch_sims, states, data, winit, *, sce_bins=None,
              bl_tau_s=None, bl_temp_offset_k=None,
              bl_wall_water_mode=None, wall_liquid_water_init_kg=None,
              wall_ice_water_init_kg=None, wall_water_efficiency=None,
              wall_vapour_transfer_velocity_ms=None, bl_evap_size_exp=None):
    """Generate namelists, run the BMM and return aerosol mode totals."""
    cfg.OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    if not cfg.BMM_EXECUTABLE.exists():
        raise FileNotFoundError(f"BMM executable not found: {cfg.BMM_EXECUTABLE}")
    if not cfg.BMM_TEMPLATE.exists():
        raise FileNotFoundError(f"BMM template not found: {cfg.BMM_TEMPLATE}")

    mode1 = np.zeros(len(batch_sims))
    mode2 = np.zeros(len(batch_sims))

    with tempfile.TemporaryDirectory(prefix="iskylab_bmm_") as tmpdir:
        tmpdir = Path(tmpdir)
        scefile = _make_temp_sce_namelist(sce_bins, tmpdir)
        for nn, exp in enumerate(batch_sims):
            output_file = cfg.OUTPUT_ROOT / f"output{nn:03d}.nc"
            namelist_text, aer = make_namelist(
                exp, states[exp], data, output_file, winit=winit[nn],
                scefile=scefile, bl_tau_s=bl_tau_s,
                bl_temp_offset_k=bl_temp_offset_k,
                bl_wall_water_mode=bl_wall_water_mode,
                wall_liquid_water_init_kg=wall_liquid_water_init_kg,
                wall_ice_water_init_kg=wall_ice_water_init_kg,
                wall_water_efficiency=wall_water_efficiency,
                wall_vapour_transfer_velocity_ms=wall_vapour_transfer_velocity_ms,
                bl_evap_size_exp=bl_evap_size_exp,
            )
            namelist_file = tmpdir / f"namelist-{exp}.in"
            write_text(namelist_file, namelist_text)

            if cfg.SAVE_GENERATED_NAMELISTS:
                shutil.copy2(namelist_file, cfg.OUTPUT_ROOT / f"namelist-{exp}.in")

            print(f"Run {nn:03d}: {exp} -> {output_file}")
            try:
                completed = subprocess.run(
                    [str(cfg.BMM_EXECUTABLE), str(namelist_file)],
                    cwd=cfg.BMM_MODEL_FOLDER,
                    check=True,
                    text=True,
                    capture_output=True,
                )
            except subprocess.CalledProcessError as exc:
                # Preserve the model diagnostics that are otherwise hidden by
                # capture_output when a batch member fails.
                if exc.stdout:
                    print(exc.stdout)
                if exc.stderr:
                    print(exc.stderr)
                raise
            if completed.stdout.strip():
                print(completed.stdout, end="" if completed.stdout.endswith("\n") else "\n")
            if completed.stderr.strip():
                print(completed.stderr, end="" if completed.stderr.endswith("\n") else "\n")

            mode1[nn] = np.sum(aer["n1"])
            mode2[nn] = np.sum(aer["n2"])

    return mode1, mode2


def peak_cdnc(cloud_window, opc):
    """Return peak observed CDNC inside the prescribed OPC analysis window."""
    mask = (opc["Time"] >= cloud_window[0]) & (opc["Time"] < cloud_window[1])
    if not np.any(mask):
        return np.nan
    return float(np.nanmax(opc["ndrop"][mask]))


def analyse_batch(batch_sims, states, data, vals1, vals2, group_type):
    """Reproduce the historical activation-number comparison plots."""
    from netCDF4 import Dataset
    ndrop_model = np.zeros(len(batch_sims))
    for i, exp in enumerate(batch_sims):
        path = cfg.OUTPUT_ROOT / f"output{i:03d}.nc"
        with Dataset(path) as nc:
            nwat = np.asarray(nc["nwat"][:])
            mwat = np.asarray(nc["mwat"][:])
            if "dwet" in nc.variables:
                dwet = np.asarray(nc["dwet"][:])
            else:
                dwet = np.where(
                    mwat > 0.0,
                    (mwat / (np.pi / 6.0 * 1000.0)) ** (1.0 / 3.0),
                    0.0,
                )
            time = np.asarray(nc["time"][:])
            p = np.asarray(nc["p"][:])
            t = np.asarray(nc["t"][:])
            batch_model = {
                "time": time,
                "nwat": nwat,
                "mwat": mwat,
                "dwet": dwet,
            }
            if "mbinedges" in nc.variables:
                batch_model["mbinedges"] = np.asarray(nc["mbinedges"][:])
            if "bin_scheme_flag" in nc.variables:
                batch_model["bin_scheme_flag"] = np.asarray(nc["bin_scheme_flag"][:])
            elif "bin_scheme_flag" in nc.ncattrs():
                batch_model["bin_scheme_flag"] = np.asarray(
                    [nc.getncattr("bin_scheme_flag")]
                )
            fraction = _warm_fraction_above_diameter(
                batch_model, cfg.SINGLE_COMPARE_DROP_MIN_UM
            )
            conc = np.sum(
                nwat * fraction, axis=tuple(range(1, nwat.ndim))
            )
            window = meta.CLOUD_WINDOWS[exp]
            mask = (time > window[0]) & (time <= window[1])
            if not np.any(mask):
                ndrop_model[i] = np.nan
            else:
                rho_d = p[mask] / (R_D * t[mask])
                ndrop_model[i] = np.nanmax(conc[mask] * rho_d) / 1.0e6

    cdnc_obs = np.zeros(len(batch_sims))
    for i, exp in enumerate(batch_sims):
        opc = data[f"MergedOPC-{exp}"]
        # Recalculate observed drop number above the configured diameter cutoff using the historical OPC
        # calibration factor; retain the supplied ndrop field separately.
        diam_mask = opc["Dp"] > float(cfg.SINGLE_COMPARE_DROP_MIN_UM)
        opc["ndrop"] = np.nansum(opc["Conc"][:, diam_mask], axis=1) * 0.009839
        cdnc_obs[i] = peak_cdnc(meta.CLOUD_WINDOWS[exp], opc)

    # Convert the initial BMM number mixing ratios back to cm-3 at the chosen
    # aerosol-normalisation state for plotting.
    vals1_cm3 = np.empty_like(vals1)
    vals2_cm3 = np.empty_like(vals2)
    for i, exp in enumerate(batch_sims):
        state = states[exp]
        met = data[f"MeteoCPC-{exp}"]
        ta = _interp_at(met["Time"], met["Tgw mean"], state["aerosol_norm_time"]) + 273.15
        pa = _interp_at(met["Time"], met["Pressure"], state["aerosol_norm_time"]) * 100.0
        rho_d = pa / (R_D * ta)
        vals1_cm3[i] = vals1[i] * rho_d / 1.0e6
        vals2_cm3[i] = vals2[i] * rho_d / 1.0e6

    x = vals1_cm3 if group_type == 1 else vals2_cm3
    xlabel = "First aerosol-component number (cm$^{-3}$)" if group_type == 1 else "Second aerosol-component number (cm$^{-3}$)"

    fig = plt.figure(figsize=(15, 4))
    ax = fig.add_subplot(131)
    ax.plot(x, ndrop_model / ndrop_model[0], ".-", ms=10, color="#0072B2", label="BMM")
    ax.plot(x, cdnc_obs / cdnc_obs[0], ".-", ms=10, color="#000000", label="iSKYLAB/WELAS")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Relative number of drops")
    ax.legend()
    ax.grid()

    ax = fig.add_subplot(132)
    denom = vals1_cm3 + vals2_cm3
    ax.plot(x, ndrop_model / denom, ".-", ms=10, color="#0072B2", label="BMM")
    ax.plot(x, cdnc_obs / denom, ".-", ms=10, color="#000000", label="iSKYLAB/WELAS")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Activated fraction")
    ax.legend()
    ax.set_ylim((0, 1))
    ax.grid()

    ax = fig.add_subplot(133)
    ax.plot(x, ndrop_model, ".-", ms=10, color="#0072B2", label="BMM")
    ax.plot(x, cdnc_obs, ".-", ms=10, color="#000000", label="iSKYLAB/WELAS")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Number of drops (cm$^{-3}$)")
    ax.legend()
    ax.grid()
    fig.tight_layout()
    return fig


def _time_edges(time):
    """Return plotting edges for a strictly increasing time-centre coordinate."""
    time = np.asarray(time, dtype=float)
    if time.ndim != 1 or len(time) < 2:
        raise ValueError("time coordinate must contain at least two samples")
    edges = np.empty(len(time) + 1, dtype=float)
    edges[1:-1] = 0.5 * (time[:-1] + time[1:])
    edges[0] = time[0] - 0.5 * (time[1] - time[0])
    edges[-1] = time[-1] + 0.5 * (time[-1] - time[-2])
    return edges


def _opc_geometry_and_moments(opc, drop_min_um):
    """Reconstruct OPC bin geometry and cloud-drop moments from the PSD itself.

    This intentionally does not require ``readOPC_Merged.py`` to provide the
    newer ``Dp_edges``, ``dlogD`` or bulk-moment keys.  Older iSKYLAB reader
    versions only return ``Dp`` and ``Conc``.  ``Conc`` is dN/dlog10(D), so the
    represented number in bin j is Conc_j * Delta log10(D)_j.
    """
    dp = np.asarray(opc["Dp"], dtype=float)
    conc = np.asarray(opc["Conc"], dtype=float)
    if dp.ndim != 1 or dp.size < 2:
        raise ValueError("OPC Dp must contain at least two diameter-bin centres")
    if conc.ndim != 2 or conc.shape[1] != dp.size:
        raise ValueError(
            f"OPC Conc/Dp shape mismatch: Conc={conc.shape}, Dp={dp.shape}"
        )
    if np.any(~np.isfinite(dp)) or np.any(dp <= 0.0) or np.any(np.diff(dp) <= 0.0):
        raise ValueError("OPC Dp centres must be finite, positive and strictly increasing")

    # Prefer reader-supplied edges when present and valid, otherwise reconstruct
    # logarithmic edges from adjacent centre midpoints.  The reconstruction is
    # exactly what the updated readOPC_Merged.py uses.
    edges = np.asarray(opc.get("Dp_edges", []), dtype=float)
    if edges.shape != (dp.size + 1,) or np.any(~np.isfinite(edges)) or np.any(edges <= 0.0):
        logc = np.log10(dp)
        loge = np.empty(dp.size + 1, dtype=float)
        loge[1:-1] = 0.5 * (logc[:-1] + logc[1:])
        loge[0] = logc[0] - 0.5 * (logc[1] - logc[0])
        loge[-1] = logc[-1] + 0.5 * (logc[-1] - logc[-2])
        edges = 10.0 ** loge

    dlog = np.diff(np.log10(edges))
    if np.any(~np.isfinite(dlog)) or np.any(dlog <= 0.0):
        raise ValueError("Reconstructed OPC logarithmic bin widths are invalid")

    nbin = conc * dlog[None, :]
    drop_mask = dp > float(drop_min_um)
    w = np.where(drop_mask[None, :], nbin, 0.0)
    m0 = np.nansum(w, axis=1)
    m1 = np.nansum(w * dp[None, :], axis=1)
    m2 = np.nansum(w * dp[None, :]**2, axis=1)
    m3 = np.nansum(w * dp[None, :]**3, axis=1)

    dmean = np.full(m0.shape, np.nan, dtype=float)
    dvol = np.full(m0.shape, np.nan, dtype=float)
    deff = np.full(m0.shape, np.nan, dtype=float)
    rel_disp = np.full(m0.shape, np.nan, dtype=float)
    good0 = m0 > 0.0
    good2 = m2 > 0.0
    dmean[good0] = m1[good0] / m0[good0]
    dvol[good0] = (m3[good0] / m0[good0]) ** (1.0 / 3.0)
    deff[good2] = m3[good2] / m2[good2]
    var = np.zeros_like(m0)
    var[good0] = np.maximum(m2[good0] / m0[good0] - dmean[good0]**2, 0.0)
    good_disp = good0 & (dmean > 0.0)
    rel_disp[good_disp] = np.sqrt(var[good_disp]) / dmean[good_disp]

    # WELAS-equivalent liquid water represented by the same D > drop_min_um
    # population.  The instrument interprets each optical/wet diameter as a
    # pure spherical water drop, so deliberately apply that same retrieval
    # rather than trying to subtract dry aerosol volume.  nbin is # cm-3 and
    # dp is micrometres; convert both to SI before summing water mass/air volume.
    rho_w = 1000.0
    lwc_equiv_kg_m3 = (
        rho_w * np.pi / 6.0
        * np.nansum(w * 1.0e6 * (dp[None, :] * 1.0e-6) ** 3, axis=1)
    )

    return {
        "Dp_um": dp,
        "Dp_edges_um": edges,
        "dlogD": dlog,
        "psd": conc,
        "ndrop_psd": m0,
        "dmean_um": dmean,
        "dvol_um": dvol,
        "deff_um": deff,
        "rel_disp": rel_disp,
        "lwc_equiv_kg_m3": lwc_equiv_kg_m3,
    }


def _observed_bulk_series(exp, data):
    """Return OPC bulk variables in units directly comparable with BMM output."""
    opc = data[f"MergedOPC-{exp}"]
    met = data[f"MeteoCPC-{exp}"]
    tobs = np.asarray(opc["Time"], dtype=float)
    p_hpa = np.interp(tobs, met["Time"], met["Pressure"])
    t_k = np.interp(tobs, met["Time"], met["Tgw mean"]) + 273.15
    rho_d = p_hpa * 100.0 / (R_D * t_k)
    ql_total = np.asarray(opc["lwc"], dtype=float) * 1.0e-3 / rho_d
    psd = _opc_geometry_and_moments(opc, cfg.SINGLE_COMPARE_DROP_MIN_UM)
    ql_above_min = np.asarray(psd["lwc_equiv_kg_m3"], dtype=float) / rho_d
    return {
        "time": tobs,
        "rho_d": rho_d,
        # Keep the historical reader-total WELAS LWC and also a quantity
        # reconstructed from exactly the D > SINGLE_COMPARE_DROP_MIN_UM PSD.
        "ql": ql_total,
        "ql_total": ql_total,
        "ql_above_min": ql_above_min,
        # Keep the OPC file's supplied Nd as a separate diagnostic, but derive
        # all PSD-based comparison moments consistently from Conc and Dp.
        "ndrop": np.asarray(opc["ndrop"], dtype=float),
        "nice": np.asarray(opc.get("nice", np.full_like(tobs, np.nan)), dtype=float),
        **psd,
    }


def _model_bin_scheme_flag(model):
    """Return BMM bin_scheme_flag from model output; old files default to 0."""
    if "bin_scheme_flag" not in model:
        return 0
    value = np.asarray(model["bin_scheme_flag"]).reshape(-1)
    if value.size == 0 or not np.isfinite(value[0]):
        return 0
    return int(round(float(value[0])))


def _fixed_warm_mass_edges_for_native_shape(model, warm_shape):
    """Broadcast fixed water-mass lower/upper edges onto native warm-bin shape."""
    if "mbinedges" not in model:
        return None, None
    edges = np.asarray(model["mbinedges"], dtype=float)
    if edges.ndim != 2 or len(warm_shape) != 3:
        return None, None

    a, b = warm_shape[1], warm_shape[2]
    candidates = (edges, edges.T)
    for e in candidates:
        if e.shape == (a, b + 1):
            lower = e[:, :-1][None, :, :]
            upper = e[:, 1:][None, :, :]
            return np.broadcast_to(lower, warm_shape), np.broadcast_to(upper, warm_shape)
        if e.shape == (a + 1, b):
            lower = e[:-1, :][None, :, :]
            upper = e[1:, :][None, :, :]
            return np.broadcast_to(lower, warm_shape), np.broadcast_to(upper, warm_shape)
    return None, None


def _warm_fraction_above_diameter(model, dmin_um):
    """Return the native warm-bin number fraction above a wet-D threshold.

    Full-moving (flag 0) remains a pointwise 0/1 test.  For fixed water-mass
    schemes 1 and 2, Dmin is converted to the water mass required for that wet
    diameter after subtracting the inferred dry-particle volume.  If the
    threshold lies inside a fixed mass bin,

        f = (mupper - mthreshold) / (mupper - mlower),

    clipped to [0,1].
    """
    if "dwet" not in model or "nwat" not in model:
        return None
    d = np.asarray(model["dwet"], dtype=float)
    n = np.asarray(model["nwat"], dtype=float)
    if d.shape != n.shape:
        raise ValueError(f"dwet/nwat shape mismatch: {d.shape} vs {n.shape}")

    dmin = float(dmin_um) * 1.0e-6
    point = np.where(
        np.isfinite(d) & np.isfinite(n) & (d > dmin) & (n > 0.0), 1.0, 0.0
    )

    scheme = _model_bin_scheme_flag(model)
    if scheme == 0:
        return point
    if scheme not in (1, 2):
        raise ValueError(f"Unknown bin_scheme_flag={scheme}")
    if "mwat" not in model or "mbinedges" not in model:
        return point

    mwat = np.asarray(model["mwat"], dtype=float)
    if mwat.shape != d.shape:
        raise ValueError(f"mwat/dwet shape mismatch: {mwat.shape} vs {d.shape}")

    mlower, mupper = _fixed_warm_mass_edges_for_native_shape(model, d.shape)
    if mlower is None:
        return point

    rho_w = 1000.0
    wet_volume = np.pi / 6.0 * np.maximum(d, 0.0)**3
    dry_volume = np.maximum(wet_volume - np.maximum(mwat, 0.0) / rho_w, 0.0)
    threshold_volume = np.pi / 6.0 * dmin**3
    mthreshold = rho_w * np.maximum(threshold_volume - dry_volume, 0.0)

    width = mupper - mlower
    if np.any(~np.isfinite(width)) or np.any(width <= 0.0):
        raise ValueError("Invalid fixed water-mass bin width in model output")

    frac = np.where(
        mthreshold <= mlower, 1.0,
        np.where(
            mthreshold >= mupper, 0.0,
            (mupper - mthreshold) / width,
        ),
    )
    frac = np.clip(frac, 0.0, 1.0)
    return np.where(np.isfinite(d) & np.isfinite(n) & (n > 0.0), frac, 0.0)


def _model_bulk_moments_above(model, dmin_um):
    """Return OPC-like BMM number/size moments above a wet-diameter threshold.

    These moments are calculated from the raw native `(nwat,dwet)` output, not
    from the model activation flag.  For fixed-bin schemes 1 and 2, a cutoff
    inside a native water-mass bin contributes only the corresponding fraction
    of that bin; full-moving particles remain pointwise.
    """
    out = {name: np.full_like(model["time"], np.nan, dtype=float) for name in
           ("number_cm3", "dmean_um", "dvol_um", "deff_um", "rel_disp")}
    if "dwet" not in model or "nwat" not in model:
        return out
    d = np.asarray(model["dwet"], dtype=float)
    n = np.asarray(model["nwat"], dtype=float)
    if d.shape != n.shape:
        raise ValueError(f"dwet/nwat shape mismatch: {d.shape} vs {n.shape}")
    axes = tuple(range(1, n.ndim))
    fraction = _warm_fraction_above_diameter(model, dmin_um)
    w = n * fraction
    m0 = np.sum(w, axis=axes)
    m1 = np.sum(w*d, axis=axes)
    m2 = np.sum(w*d**2, axis=axes)
    m3 = np.sum(w*d**3, axis=axes)
    good0 = m0 > 0.0
    good2 = m2 > 0.0
    dmean = np.full_like(m0, np.nan, dtype=float)
    dvol = np.full_like(m0, np.nan, dtype=float)
    deff = np.full_like(m0, np.nan, dtype=float)
    rel = np.full_like(m0, np.nan, dtype=float)
    dmean[good0] = m1[good0] / m0[good0]
    dvol[good0] = (m3[good0] / m0[good0])**(1.0/3.0)
    deff[good2] = m3[good2] / m2[good2]
    variance = np.zeros_like(m0)
    variance[good0] = np.maximum(m2[good0]/m0[good0] - dmean[good0]**2, 0.0)
    rel[good0] = np.sqrt(variance[good0]) / np.maximum(dmean[good0], np.finfo(float).tiny)
    out["number_cm3"] = m0 * np.asarray(model["rhoa"]) / 1.0e6
    out["dmean_um"] = dmean * 1.0e6
    out["dvol_um"] = dvol * 1.0e6
    out["deff_um"] = deff * 1.0e6
    out["rel_disp"] = rel
    return out


def _model_total_hydrometeor_moments_above(model, dmin_um):
    """Return BMM liquid+ice number/size moments above a diameter threshold.

    The observed OPC/WELAS PSD in mixed-phase experiments is not phase
    resolved.  For a like-for-like size-moment comparison, combine the BMM
    warm distribution (nwat,dwet) and ice distribution (nicem,dmaxice) before
    forming Dmean, Deff and relative dispersion.  Both phases use the same
    optical-diameter lower cutoff used for the observed PSD diagnostics.
    """
    out = {name: np.full_like(model["time"], np.nan, dtype=float) for name in
           ("number_cm3", "dmean_um", "dvol_um", "deff_um", "rel_disp")}

    nt = len(np.asarray(model["time"]))
    m0 = np.zeros(nt, dtype=float)
    m1 = np.zeros(nt, dtype=float)
    m2 = np.zeros(nt, dtype=float)
    m3 = np.zeros(nt, dtype=float)
    have_any = False
    dmin = float(dmin_um) * 1.0e-6

    for dkey, nkey in (("dwet", "nwat"), ("dmaxice", "nicem")):
        if dkey not in model or nkey not in model:
            continue
        d = np.asarray(model[dkey], dtype=float)
        n = np.asarray(model[nkey], dtype=float)
        if d.shape != n.shape:
            raise ValueError(f"{dkey}/{nkey} shape mismatch: {d.shape} vs {n.shape}")
        if d.shape[0] != nt:
            raise ValueError(
                f"{dkey}/{nkey} time dimension {d.shape[0]} does not match time {nt}"
            )
        axes = tuple(range(1, n.ndim))
        if dkey == "dwet" and nkey == "nwat":
            fraction = _warm_fraction_above_diameter(model, dmin_um)
            w = np.where(np.isfinite(n), n * fraction, 0.0)
        else:
            # Ice Dmax is shape-dependent, so keep ice pointwise at Dmin.
            good = np.isfinite(d) & np.isfinite(n) & (d > dmin) & (n > 0.0)
            w = np.where(good, n, 0.0)
        m0 += np.sum(w, axis=axes)
        m1 += np.sum(w * d, axis=axes)
        m2 += np.sum(w * d**2, axis=axes)
        m3 += np.sum(w * d**3, axis=axes)
        have_any = True

    if not have_any:
        return out

    good0 = m0 > 0.0
    good2 = m2 > 0.0
    dmean = np.full(nt, np.nan, dtype=float)
    dvol = np.full(nt, np.nan, dtype=float)
    deff = np.full(nt, np.nan, dtype=float)
    rel = np.full(nt, np.nan, dtype=float)
    dmean[good0] = m1[good0] / m0[good0]
    dvol[good0] = (m3[good0] / m0[good0])**(1.0/3.0)
    deff[good2] = m3[good2] / m2[good2]
    variance = np.zeros(nt, dtype=float)
    variance[good0] = np.maximum(m2[good0]/m0[good0] - dmean[good0]**2, 0.0)
    rel[good0] = np.sqrt(variance[good0]) / np.maximum(
        dmean[good0], np.finfo(float).tiny
    )

    out["number_cm3"] = m0 * np.asarray(model["rhoa"]) / 1.0e6
    out["dmean_um"] = dmean * 1.0e6
    out["dvol_um"] = dvol * 1.0e6
    out["deff_um"] = deff * 1.0e6
    out["rel_disp"] = rel
    return out


def _model_welas_equivalent_ql_above(model, dmin_um):
    """Return WELAS-equivalent BMM liquid water for Dwet > dmin_um.

    WELAS interprets the complete measured wet/optical sphere as liquid water.
    Apply that same virtual-instrument operator to every BMM warm particle
    (nwat+dwet) above the configured lower diameter threshold.  No upper cutoff
    is imposed: this intentionally means *everything* above
    SINGLE_COMPARE_DROP_MIN_UM.

    The returned units are kg liquid-equivalent water per kg dry air because
    nwat is # kg-1 dry air and rho_w*pi/6*D^3 is kg per interpreted drop.
    """
    if "dwet" not in model or "nwat" not in model:
        return np.full_like(model["time"], np.nan, dtype=float)

    d = np.asarray(model["dwet"], dtype=float)
    n = np.asarray(model["nwat"], dtype=float)
    if d.shape != n.shape:
        raise ValueError(f"dwet/nwat shape mismatch: {d.shape} vs {n.shape}")

    fraction = _warm_fraction_above_diameter(model, dmin_um)
    axes = tuple(range(1, n.ndim))
    rho_w = 1000.0
    weighted_number = np.where(
        np.isfinite(d) & np.isfinite(n), n * fraction, 0.0
    )
    return (
        rho_w * np.pi / 6.0
        * np.sum(weighted_number * d**3, axis=axes)
    )


def _model_instrument_number(model, dmin_um):
    """Backward-compatible shorthand for the instrument-like number series."""
    return _model_bulk_moments_above(model, dmin_um)["number_cm3"]


def _common_comparison_psd_edges(obs_edges_um):
    """Return a deliberately coarser common log-D grid for WELAS/BMM plots."""
    obs_edges_um = np.asarray(obs_edges_um, dtype=float)
    if obs_edges_um.ndim != 1 or len(obs_edges_um) < 2:
        raise ValueError("Observed PSD edges must contain at least two values")
    dmin = (
        float(obs_edges_um[0])
        if cfg.SINGLE_COMPARE_PSD_MIN_UM is None
        else float(cfg.SINGLE_COMPARE_PSD_MIN_UM)
    )
    dmax = (
        float(obs_edges_um[-1])
        if cfg.SINGLE_COMPARE_PSD_MAX_UM is None
        else float(cfg.SINGLE_COMPARE_PSD_MAX_UM)
    )
    nbins = int(cfg.SINGLE_COMPARE_PSD_NBINS)
    if nbins < 8:
        raise ValueError("SINGLE_COMPARE_PSD_NBINS must be >= 8")
    if dmin <= 0.0 or dmax <= dmin:
        raise ValueError("Invalid SINGLE_COMPARE_PSD_MIN/MAX_UM")
    # Do not extend beyond WELAS coverage in a direct observation comparison.
    dmin = max(dmin, float(obs_edges_um[0]))
    dmax = min(dmax, float(obs_edges_um[-1]))
    if dmax <= dmin:
        raise ValueError("Requested comparison PSD range does not overlap WELAS")
    return np.logspace(np.log10(dmin), np.log10(dmax), nbins + 1)


def _rebin_dndlog10d(source_psd, source_edges_um, target_edges_um):
    """Conservatively rebin dN/dlog10D between logarithmic diameter grids.

    The source PSD is assumed piecewise constant in log10(D) inside each source
    channel.  Number is integrated over the exact log-diameter overlap with
    every target bin, then divided by the target log width.  Thus the integral
    of the rebinned PSD is conserved over the common diameter range.
    """
    source_psd = np.asarray(source_psd, dtype=float)
    source_edges_um = np.asarray(source_edges_um, dtype=float)
    target_edges_um = np.asarray(target_edges_um, dtype=float)
    if source_psd.ndim != 2:
        raise ValueError("source_psd must have shape (time, diameter_bin)")
    if source_psd.shape[1] != len(source_edges_um) - 1:
        raise ValueError("source PSD/edge size mismatch")
    if np.any(source_edges_um <= 0.0) or np.any(np.diff(source_edges_um) <= 0.0):
        raise ValueError("source PSD edges must be positive and increasing")
    if np.any(target_edges_um <= 0.0) or np.any(np.diff(target_edges_um) <= 0.0):
        raise ValueError("target PSD edges must be positive and increasing")

    sl = np.log10(source_edges_um)
    tl = np.log10(target_edges_um)
    overlap = np.maximum(
        0.0,
        np.minimum(sl[1:, None], tl[None, 1:])
        - np.maximum(sl[:-1, None], tl[None, :-1]),
    )
    # dN/dlogD * overlap width gives number in each source/target overlap.
    number_target = np.nan_to_num(source_psd, nan=0.0) @ overlap
    return number_target / np.diff(tl)[None, :]


def _model_psd_on_fixed_grid(model, edges_um, *, diameter_key="dwet", number_key="nwat"):
    """Conservatively histogram a native moving BMM population onto fixed log-D bins.

    ``number_key`` is a number mixing ratio [kg-1 dry air] carried at the
    instantaneous diameter in ``diameter_key`` [m].  Native BMM warm/ice bins
    are *not* interpreted as fixed wet-size intervals; each native point is
    simply assigned, with its complete number weight, to a diagnostic diameter
    interval.  The result is dN/dlog10D in cm-3.
    """
    if diameter_key not in model or number_key not in model:
        return None
    d = np.asarray(model[diameter_key], dtype=float) * 1.0e6
    n = np.asarray(model[number_key], dtype=float)
    rho = np.asarray(model["rhoa"], dtype=float)
    if d.shape != n.shape or d.shape[0] != len(rho):
        raise ValueError(
            f"Inconsistent model PSD dimensions for {diameter_key}/{number_key}: "
            f"{d.shape} vs {n.shape}; nt={len(rho)}"
        )
    edges_um = np.asarray(edges_um, dtype=float)
    if np.any(edges_um <= 0.0) or np.any(np.diff(edges_um) <= 0.0):
        raise ValueError("PSD diameter edges must be positive and strictly increasing")
    dlog = np.diff(np.log10(edges_um))
    out = np.zeros((len(rho), len(dlog)), dtype=float)
    for it in range(len(rho)):
        diam = d[it].ravel()
        weights = n[it].ravel() * rho[it] / 1.0e6  # # cm-3 represented natively
        good = np.isfinite(diam) & np.isfinite(weights) & (diam > 0.0) & (weights > 0.0)
        if np.any(good):
            per_bin, _ = np.histogram(diam[good], bins=edges_um, weights=weights[good])
            out[it] = per_bin / dlog
    return out


def _model_psd_on_opc_grid(model, edges_um):
    """Return the *complete warm* BMM PSD on the fixed OPC diameter grid.

    This deliberately uses ``nwat+dwet`` rather than the activated-only
    ``nliq`` field.  Unactivated aerosol and haze therefore remain in the same
    model distribution; particles outside the OPC diameter range are simply
    outside this particular instrument comparison.
    """
    return _model_psd_on_fixed_grid(
        model, edges_um, diameter_key="dwet", number_key="nwat"
    )


def _cloud_comparison_window(exp, obs):
    """Choose a cloud comparison window, preferring curated metadata when available."""
    if exp in meta.CLOUD_WINDOWS:
        return np.asarray(meta.CLOUD_WINDOWS[exp], dtype=float)
    onset = meta.CLOUD_ONSET.get(exp, float(obs["time"][0]))
    nd = np.asarray(obs["ndrop"], dtype=float)
    good = np.isfinite(nd) & (obs["time"] >= onset)
    if np.any(good):
        peak = np.nanmax(nd[good])
        active = good & (nd >= max(1.0, 0.05 * peak))
        if np.any(active):
            return np.array([float(obs["time"][active][0]), float(obs["time"][active][-1])])
    return np.array([onset, float(obs["time"][-1])])


def _best_ql_lag(model_time, model_ql, obs_time, obs_ql, window):
    """Diagnostic lag: positive means the observation occurs later than the model."""
    maxlag = float(cfg.SINGLE_COMPARE_MAX_LAG_S)
    step = float(cfg.SINGLE_COMPARE_LAG_STEP_S)
    if maxlag <= 0.0 or step <= 0.0:
        return 0.0, np.nan
    mt = np.asarray(model_time, dtype=float)
    mq = np.asarray(model_ql, dtype=float)
    ot = np.asarray(obs_time, dtype=float)
    oq = np.asarray(obs_ql, dtype=float)
    base = (mt >= window[0]) & (mt <= window[1]) & np.isfinite(mq)
    best_lag, best_corr = 0.0, -np.inf
    for lag in np.arange(-maxlag, maxlag + 0.5 * step, step):
        # If observations are delayed by +lag, obs(t+lag) corresponds to model(t).
        target_obs_time = mt[base] + lag
        inrange = (target_obs_time >= ot[0]) & (target_obs_time <= ot[-1])
        if np.count_nonzero(inrange) < 5:
            continue
        x = mq[base][inrange]
        y = np.interp(target_obs_time[inrange], ot, oq)
        good = np.isfinite(x) & np.isfinite(y)
        if np.count_nonzero(good) < 5 or np.nanstd(x[good]) <= 0.0 or np.nanstd(y[good]) <= 0.0:
            continue
        corr = float(np.corrcoef(x[good], y[good])[0, 1])
        if corr > best_corr:
            best_lag, best_corr = float(lag), corr
    return best_lag, best_corr if np.isfinite(best_corr) else np.nan


def _series_scores(model_time, model_values, obs_time, obs_values, window, lag=0.0):
    mt = np.asarray(model_time, dtype=float)
    mv = np.asarray(model_values, dtype=float)
    ot = np.asarray(obs_time, dtype=float)
    ov = np.asarray(obs_values, dtype=float)
    mask = (mt >= window[0]) & (mt <= window[1]) & np.isfinite(mv)
    query = mt[mask] + lag
    valid = (query >= ot[0]) & (query <= ot[-1])
    if np.count_nonzero(valid) < 2:
        return {"rmse": np.nan, "nrmse": np.nan, "bias": np.nan}
    x = mv[mask][valid]
    y = np.interp(query[valid], ot, ov)
    good = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(good) < 2:
        return {"rmse": np.nan, "nrmse": np.nan, "bias": np.nan}
    err = x[good] - y[good]
    rmse = float(np.sqrt(np.mean(err**2)))
    scale = float(np.nanmax(np.abs(y[good])))
    return {
        "rmse": rmse,
        "nrmse": rmse / scale if scale > 0.0 else np.nan,
        "bias": float(np.mean(err)),
    }


def run_single_experiment(exp, state, data, *, winit=1.3, sce_bins=None,
                          bl_tau_s=None, bl_temp_offset_k=None,
                          bl_wall_water_mode=None, wall_liquid_water_init_kg=None,
                          wall_ice_water_init_kg=None, wall_water_efficiency=None,
                          wall_vapour_transfer_velocity_ms=None,
                          bl_evap_size_exp=None):
    """Generate and run one named experiment, retaining its namelist and output."""
    cfg.OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    outdir = cfg.OUTPUT_ROOT / "single_comparison"
    outdir.mkdir(parents=True, exist_ok=True)
    if not cfg.BMM_EXECUTABLE.exists():
        raise FileNotFoundError(f"BMM executable not found: {cfg.BMM_EXECUTABLE}")
    output_file = outdir / f"output-{exp}.nc"

    # Keep both the temporary SCE override and the main run namelist alive for
    # the duration of main.exe.  The repository SCE namelist is untouched.
    with tempfile.TemporaryDirectory(prefix="iskylab_bmm_single_") as tmpdir:
        tmpdir = Path(tmpdir)
        scefile = _make_temp_sce_namelist(sce_bins, tmpdir)
        namelist_text, _ = make_namelist(
            exp, state, data, output_file, winit=winit, scefile=scefile,
            bl_tau_s=bl_tau_s, bl_temp_offset_k=bl_temp_offset_k,
            bl_wall_water_mode=bl_wall_water_mode,
            wall_liquid_water_init_kg=wall_liquid_water_init_kg,
            wall_ice_water_init_kg=wall_ice_water_init_kg,
            wall_water_efficiency=wall_water_efficiency,
            wall_vapour_transfer_velocity_ms=wall_vapour_transfer_velocity_ms,
            bl_evap_size_exp=bl_evap_size_exp,
        )
        run_namelist = tmpdir / f"namelist-{exp}.in"
        write_text(run_namelist, namelist_text)

        # Retain the generated main namelist for provenance.  Note that when
        # --sce-bins is used its scefile path is intentionally a temporary run
        # path and is not intended for a later manual rerun.
        namelist_file = outdir / f"namelist-{exp}.in"
        write_text(namelist_file, namelist_text)

        print(f"Single run: {exp} -> {output_file}")
        completed = subprocess.run(
            [str(cfg.BMM_EXECUTABLE), str(run_namelist)],
            cwd=cfg.BMM_MODEL_FOLDER,
            check=True,
            text=True,
            capture_output=True,
        )
    if completed.stdout.strip():
        print(completed.stdout, end="" if completed.stdout.endswith("\n") else "\n")
    if completed.stderr.strip():
        print(completed.stderr, end="" if completed.stderr.endswith("\n") else "\n")
    return output_file


def analyse_single_experiment(exp, data, model_file, *, show=True,
                              saturation_time_s=None,
                              shift_cloud_measurements=False):
    """Create time-series and PSD model/observation diagnostics for one experiment."""
    import readModel
    from matplotlib.colors import LogNorm

    model = readModel.readData(model_file, modelStr=exp)[exp]

    # The wall-reservoir diagnostics are newer than some readModel.py versions.
    # Read them directly from the NetCDF output when necessary so the analysis
    # remains compatible with both old and updated model readers.
    reservoir_fields = (
        "chamber_wall_liquid_water",
        "chamber_wall_ice_water",
        "qchamber_wall_liq_evap",
        "qchamber_wall_liq_cond",
        "qchamber_wall_ice_subl",
        "qchamber_wall_ice_dep",
        "qchamber_bl",
        "qchamber_bl_evap",
        "chamber_wall_rh",
        "chamber_wall_vapour_flux",
    )
    missing_reservoir_fields = [name for name in reservoir_fields if name not in model]
    if missing_reservoir_fields:
        from netCDF4 import Dataset
        with Dataset(model_file) as nc:
            for name in missing_reservoir_fields:
                if name in nc.variables:
                    model[name] = np.asarray(nc.variables[name][:]).squeeze()
    obs = _observed_bulk_series(exp, data)
    met = data[f"MeteoCPC-{exp}"]
    window = _cloud_comparison_window(exp, obs)
    onset = meta.CLOUD_ONSET.get(exp, window[0])

    nd_model = np.asarray(model["ndrop"]) * np.asarray(model["rhoa"]) / 1.0e6
    instrument_model = _model_bulk_moments_above(model, cfg.SINGLE_COMPARE_DROP_MIN_UM)
    total_hydrometeor_model = _model_total_hydrometeor_moments_above(
        model, cfg.SINGLE_COMPARE_DROP_MIN_UM
    )
    nd_model_2um = instrument_model["number_cm3"]
    # Liquid-only size moments, retained as useful phase-specific diagnostics.
    deff_model_um = instrument_model["deff_um"]
    rel_model = instrument_model["rel_disp"]
    deff_model_activated_um = np.asarray(model["deff"]) * 1.0e6
    rel_model_activated = np.asarray(model.get("rel_disp_liq", np.full_like(model["time"], np.nan)))
    # Total liquid+ice moments are the direct comparison to the observed merged
    # OPC/WELAS PSD once ice is present.
    deff_model_total_um = total_hydrometeor_model["deff_um"]
    rel_model_total = total_hydrometeor_model["rel_disp"]

    # Ice bulk diagnostics for the main mixed-phase comparison panels.
    qi_model = np.asarray(model.get("qi", np.full_like(model["time"], np.nan)), dtype=float)
    if "nice" in model:
        nice_model_cm3 = np.asarray(model["nice"], dtype=float) * np.asarray(model["rhoa"], dtype=float) / 1.0e6
    elif "nicem" in model:
        nice_native = np.asarray(model["nicem"], dtype=float)
        nice_model_cm3 = np.sum(nice_native, axis=tuple(range(1, nice_native.ndim))) * np.asarray(model["rhoa"], dtype=float) / 1.0e6
    else:
        nice_model_cm3 = np.full_like(model["time"], np.nan, dtype=float)
    nice_obs_cm3 = np.asarray(obs.get("nice", np.full_like(obs["time"], np.nan)), dtype=float)

    # Two deliberately different liquid-water comparisons are retained.
    # 1) total: native BMM ql versus the WELAS reader's total LWC-derived ql.
    # 2) above_min: both sides use a WELAS-like spherical-water retrieval for
    #    every particle with D > SINGLE_COMPARE_DROP_MIN_UM.
    ql_model_total = np.asarray(model["ql"], dtype=float)
    ql_model_above = _model_welas_equivalent_ql_above(
        model, cfg.SINGLE_COMPARE_DROP_MIN_UM
    )
    ql_obs_total = np.asarray(obs["ql_total"], dtype=float)
    ql_obs_above = np.asarray(obs["ql_above_min"], dtype=float)

    ql_mode = str(cfg.SINGLE_COMPARE_QL_MODE).strip().lower().replace("-", "_")
    if ql_mode == "total":
        ql_model_compare = ql_model_total
        ql_obs_compare = ql_obs_total
        ql_mode_label = "total"
    elif ql_mode == "above_min":
        ql_model_compare = ql_model_above
        ql_obs_compare = ql_obs_above
        ql_mode_label = f"WELAS-eq >{cfg.SINGLE_COMPARE_DROP_MIN_UM:g} um"
    else:
        raise ValueError(
            "SINGLE_COMPARE_QL_MODE must be 'total' or 'above_min', got "
            f"{cfg.SINGLE_COMPARE_QL_MODE!r}"
        )

    # Diagnose the WELAS/sample-line lag from the selected liquid-water series.
    # Positive lag means the observation occurs later than the model.  When
    # requested, shift the *single observed time coordinate* by -lag.  Every
    # WELAS-derived quantity (ql, number, Deff, dispersion, ice number and PSD)
    # then follows automatically and consistently.
    diagnosed_lag, lag_corr = _best_ql_lag(
        model["time"], ql_model_compare, obs["time"], ql_obs_compare, window
    )
    applied_obs_shift_s = -diagnosed_lag if shift_cloud_measurements else 0.0
    if shift_cloud_measurements and abs(applied_obs_shift_s) > 0.0:
        obs["time"] = np.asarray(obs["time"], dtype=float) + applied_obs_shift_s
        window = np.asarray(window, dtype=float) + applied_obs_shift_s
        onset = float(onset) + applied_obs_shift_s
        print(
            f"{exp}: shifted all WELAS cloud measurements by "
            f"{applied_obs_shift_s:+.1f} s (diagnosed ql lag "
            f"{diagnosed_lag:+.1f} s)"
        )

    ql_abs = _series_scores(
        model["time"], ql_model_compare, obs["time"], ql_obs_compare, window, lag=0.0
    )
    # Keep the historical lagged score for diagnostics.  If the complete WELAS
    # record has already been shifted, no additional lag is applied here.
    ql_lag = _series_scores(
        model["time"], ql_model_compare, obs["time"], ql_obs_compare, window,
        lag=0.0 if shift_cloud_measurements else diagnosed_lag
    )
    ql_total_abs = _series_scores(
        model["time"], ql_model_total, obs["time"], ql_obs_total, window, lag=0.0
    )
    ql_above_abs = _series_scores(
        model["time"], ql_model_above, obs["time"], ql_obs_above, window, lag=0.0
    )
    nd_score = _series_scores(model["time"], nd_model_2um, obs["time"], obs["ndrop_psd"], window)
    deff_score = _series_scores(model["time"], deff_model_total_um, obs["time"], obs["deff_um"], window)
    rel_score = _series_scores(model["time"], rel_model_total, obs["time"], obs["rel_disp"], window)

    # Integrated liquid-water exposure is much less sensitive to a short sample-line lag than a peak.
    mt_mask = (model["time"] >= window[0]) & (model["time"] <= window[1])
    ot_mask = (obs["time"] >= window[0]) & (obs["time"] <= window[1])
    int_ql_model = float(np.trapezoid(ql_model_compare[mt_mask], np.asarray(model["time"])[mt_mask])) if np.count_nonzero(mt_mask) > 1 else np.nan
    int_ql_obs = float(np.trapezoid(ql_obs_compare[ot_mask], obs["time"][ot_mask])) if np.count_nonzero(ot_mask) > 1 else np.nan
    int_ql_model_total = float(np.trapezoid(ql_model_total[mt_mask], np.asarray(model["time"])[mt_mask])) if np.count_nonzero(mt_mask) > 1 else np.nan
    int_ql_obs_total = float(np.trapezoid(ql_obs_total[ot_mask], obs["time"][ot_mask])) if np.count_nonzero(ot_mask) > 1 else np.nan
    int_ql_model_above = float(np.trapezoid(ql_model_above[mt_mask], np.asarray(model["time"])[mt_mask])) if np.count_nonzero(mt_mask) > 1 else np.nan
    int_ql_obs_above = float(np.trapezoid(ql_obs_above[ot_mask], obs["time"][ot_mask])) if np.count_nonzero(ot_mask) > 1 else np.nan

    # Consistent, colour-blind-friendly styling for the iSKYLAB comparison.
    # Source is primarily encoded by colour and phase/definition by linestyle,
    # so the plots remain readable in colour and when printed in greyscale.
    C_OBS = "#000000"       # observations / WELAS
    C_OBS_ALT = "#7F7F7F"   # alternate observational diagnostic
    C_BMM = "#0072B2"       # BMM liquid / primary model
    C_BMM_TOTAL = "#009E73" # BMM total liquid+ice
    C_ICE = "#D55E00"       # ice
    C_WALL = "#CC79A7"      # chamber wall
    C_FORCING = "#E69F00"   # forcing / auxiliary
    C_ACT = "#56B4E9"       # activated-liquid diagnostic

    fig, axes = plt.subplots(4, 2, figsize=(15, 17), sharex=False)
    ax = axes.ravel()

    # (a) Pressure
    ax[0].plot(
        met["Time"] / 60.0, met["Pressure"],
        color=C_OBS, linewidth=1.8, label="iSKYLAB pressure forcing",
    )
    ax[0].plot(
        model["time"] / 60.0, np.asarray(model["p"]) / 100.0,
        color=C_BMM, linestyle="--", linewidth=1.8, label="BMM pressure",
    )
    ax[0].set_title("(a) Chamber pressure")
    ax[0].set_ylabel("Pressure (hPa)")
    ax[0].legend(fontsize=8, loc="best")
    ax[0].grid(alpha=0.3)

    # (b) Temperature
    ax[1].plot(
        met["Time"] / 60.0, met["Tgw mean"],
        color=C_OBS, linewidth=1.8, label="iSKYLAB gas temperature / forcing",
    )
    if "Tww mean" in met:
        ax[1].plot(
            met["Time"] / 60.0, met["Tww mean"],
            color=C_WALL, linewidth=1.5, label="iSKYLAB wall temperature",
        )
    ax[1].plot(
        model["time"] / 60.0, np.asarray(model["t"]) - 273.15,
        color=C_BMM, linestyle="--", linewidth=1.8, label="BMM gas temperature",
    )
    ax[1].set_title("(b) Chamber temperature")
    ax[1].set_ylabel(r"Temperature ($^\circ$C)")
    ax[1].legend(fontsize=8, loc="best")
    ax[1].grid(alpha=0.3)

    # (c) Water mixing ratios.  Keep total and WELAS-equivalent definitions
    # visibly distinct while using a stable observation/model colour convention.
    total_selected = ql_mode == "total"
    above_selected = ql_mode == "above_min"
    ax[2].plot(
        obs["time"] / 60.0, ql_obs_total * 1.0e3,
        color=C_OBS,
        linestyle="-" if total_selected else "--",
        linewidth=2.0 if total_selected else 1.4,
        alpha=1.0 if total_selected else 0.65,
        label="iSKYLAB/WELAS total liquid water",
    )
    ax[2].plot(
        model["time"] / 60.0, ql_model_total * 1.0e3,
        color=C_BMM,
        linestyle="-" if total_selected else "--",
        linewidth=2.0 if total_selected else 1.4,
        alpha=1.0 if total_selected else 0.65,
        label="BMM total liquid water",
    )
    ax[2].plot(
        obs["time"] / 60.0, ql_obs_above * 1.0e3,
        color=C_OBS_ALT,
        linestyle="-" if above_selected else ":",
        linewidth=2.0 if above_selected else 1.5,
        alpha=1.0 if above_selected else 0.75,
        label=f"iSKYLAB/WELAS liquid water, D > {cfg.SINGLE_COMPARE_DROP_MIN_UM:g} µm",
    )
    ax[2].plot(
        model["time"] / 60.0, ql_model_above * 1.0e3,
        color=C_BMM,
        linestyle="-." if above_selected else ":",
        linewidth=2.0 if above_selected else 1.5,
        alpha=1.0 if above_selected else 0.75,
        label=f"BMM WELAS-equivalent liquid water, D > {cfg.SINGLE_COMPARE_DROP_MIN_UM:g} µm",
    )
    # No shifted comparison curve is drawn by default.  The diagnosed lag is
    # reported numerically only.  When --shift-cloud-measurements is supplied,
    # obs["time"] itself has already been shifted above, so every WELAS-derived
    # quantity is plotted consistently on the shifted time coordinate.
    if np.any(np.isfinite(qi_model)):
        ax[2].plot(
            model["time"] / 60.0, qi_model * 1.0e3,
            color=C_ICE, linestyle="-.", linewidth=1.8,
            label="BMM ice water",
        )
    ax[2].set_title("(c) Hydrometeor water mixing ratio")
    ax[2].set_ylabel(r"Water mixing ratio (g kg$^{-1}$)")
    ax[2].legend(fontsize=7.5, loc="best")
    ax[2].grid(alpha=0.3)

    # (d) Number concentrations.  The axis is stated explicitly in BOTH the
    # ylabel and legend so there is no ambiguity on saved/static figures.
    ax[3].plot(
        obs["time"] / 60.0, obs["ndrop"],
        color=C_OBS, linewidth=1.7,
        label="[LEFT axis] iSKYLAB/WELAS $N_d$ field",
    )
    ax[3].plot(
        obs["time"] / 60.0, obs["ndrop_psd"],
        color=C_OBS_ALT, linestyle=":", linewidth=1.7,
        label=f"[LEFT axis] iSKYLAB/WELAS $N$, D > {cfg.SINGLE_COMPARE_DROP_MIN_UM:g} µm",
    )
    ax[3].plot(
        model["time"] / 60.0, nd_model,
        color=C_BMM, linewidth=1.9,
        label="[LEFT axis] BMM activated liquid $N_d$",
    )
    ax[3].plot(
        model["time"] / 60.0, nd_model_2um,
        color=C_BMM, linestyle="--", linewidth=1.7,
        label=f"[LEFT axis] BMM liquid $N$, D > {cfg.SINGLE_COMPARE_DROP_MIN_UM:g} µm",
    )
    ax[3].set_ylabel(r"LEFT axis: liquid/drop number (cm$^{-3}$)")

    ax3i = ax[3].twinx()
    if np.any(np.isfinite(nice_obs_cm3)):
        ax3i.plot(
            obs["time"] / 60.0, nice_obs_cm3,
            color=C_ICE, linestyle="-.", linewidth=1.8,
            label="[RIGHT axis] iSKYLAB/WELAS ice $N_i$",
        )
    if np.any(np.isfinite(nice_model_cm3)):
        ax3i.plot(
            model["time"] / 60.0, nice_model_cm3,
            color=C_ICE, linestyle=":", linewidth=2.0,
            label="[RIGHT axis] BMM ice $N_i$",
        )
    ax3i.set_ylabel(r"RIGHT axis: ice-crystal number (cm$^{-3}$)", color=C_ICE)
    ax3i.tick_params(axis="y", colors=C_ICE)
    ax3i.spines["right"].set_color(C_ICE)
    ax[3].set_title("(d) Particle number concentration — dual axis")

    h1, l1 = ax[3].get_legend_handles_labels()
    h2, l2 = ax3i.get_legend_handles_labels()
    ax[3].legend(h1 + h2, l1 + l2, fontsize=7.2, loc="best")
    ax[3].grid(alpha=0.3)

    # (e) Effective diameter
    ax[4].plot(
        obs["time"] / 60.0, obs["deff_um"],
        color=C_OBS, linewidth=1.8,
        label="iSKYLAB/WELAS total PSD",
    )
    ax[4].plot(
        model["time"] / 60.0, deff_model_total_um,
        color=C_BMM_TOTAL, linewidth=2.0,
        label=f"BMM total liquid + ice, D > {cfg.SINGLE_COMPARE_DROP_MIN_UM:g} µm",
    )
    ax[4].plot(
        model["time"] / 60.0, deff_model_um,
        color=C_BMM, linestyle="--", linewidth=1.6,
        label=f"BMM liquid only, D > {cfg.SINGLE_COMPARE_DROP_MIN_UM:g} µm",
    )
    ax[4].plot(
        model["time"] / 60.0, deff_model_activated_um,
        color=C_ACT, linestyle=":", linewidth=1.6,
        label="BMM activated liquid",
    )
    ax[4].set_title(r"(e) Effective diameter, $D_{eff}$")
    ax[4].set_ylabel(r"$D_{eff}$ ($\mu$m)")
    ax[4].legend(fontsize=7.5, loc="best")
    ax[4].grid(alpha=0.3)

    # (f) Relative dispersion
    ax[5].plot(
        obs["time"] / 60.0, obs["rel_disp"],
        color=C_OBS, linewidth=1.8,
        label="iSKYLAB/WELAS total PSD",
    )
    ax[5].plot(
        model["time"] / 60.0, rel_model_total,
        color=C_BMM_TOTAL, linewidth=2.0,
        label=f"BMM total liquid + ice, D > {cfg.SINGLE_COMPARE_DROP_MIN_UM:g} µm",
    )
    ax[5].plot(
        model["time"] / 60.0, rel_model,
        color=C_BMM, linestyle="--", linewidth=1.6,
        label=f"BMM liquid only, D > {cfg.SINGLE_COMPARE_DROP_MIN_UM:g} µm",
    )
    ax[5].plot(
        model["time"] / 60.0, rel_model_activated,
        color=C_ACT, linestyle=":", linewidth=1.6,
        label="BMM activated liquid",
    )
    ax[5].set_title(r"(f) Relative diameter dispersion")
    ax[5].set_ylabel(r"Relative dispersion, $\sigma_D/\bar D$")
    ax[5].legend(fontsize=7.5, loc="best")
    ax[5].grid(alpha=0.3)

    # (g) Cumulative chamber/fallout diagnostics
    loss_fields = [
        ("qfan_liq", "Fan loss", C_BMM, "-"),
        ("qwall_liq", "Particle wall loss", C_WALL, "-"),
        ("qfall_liq", "Fallout", C_FORCING, "--"),
        ("qchamber_bl", "BL wall-water net loss (+loss, -source)", C_BMM_TOTAL, "-."),
        ("qchamber_wall_liq_evap", "Wall-liquid evaporation source", C_FORCING, ":"),
        ("qchamber_wall_liq_cond", "Wall-liquid condensation sink", C_WALL, "--"),
        ("qchamber_wall_ice_subl", "Wall-frost sublimation source", C_ACT, ":"),
        ("qchamber_wall_ice_dep", "Wall-frost deposition sink", C_ICE, "--"),
        ("qchamber_bl_evap", "BL particle liquid → vapour", C_ICE, ":"),
    ]
    any_loss = False
    for name, label, colour, linestyle in loss_fields:
        if name in model:
            ax[6].plot(
                model["time"] / 60.0, np.asarray(model[name]) * 1.0e3,
                color=colour, linestyle=linestyle, linewidth=1.8, label=label,
            )
            any_loss = True
    ax[6].set_title("(g) Chamber water exchange and wall reservoirs")
    ax[6].set_ylabel(r"LEFT axis: cumulative water diagnostic (g kg$^{-1}$)")

    # Prognostic wall reservoirs are physical chamber masses [kg], deliberately
    # not mixing ratios.  Plot them in grams on a separate axis: the expected
    # reservoir is normally only a few grams, so a kg-scale axis can hide it.
    ax6r = ax[6].twinx()
    any_reservoir = False
    reservoir_max_g = 0.0
    for name, label, colour, linestyle in (
        ("chamber_wall_liquid_water", "Wall liquid reservoir", C_WALL, "-"),
        ("chamber_wall_ice_water", "Wall ice/frost reservoir", C_ICE, "-."),
    ):
        if name in model:
            reservoir_g = 1.0e3 * np.asarray(model[name], dtype=float)
            ax6r.plot(
                model["time"] / 60.0, reservoir_g,
                color=colour, linestyle=linestyle, linewidth=2.0,
                label=f"[RIGHT axis] {label}",
            )
            finite = reservoir_g[np.isfinite(reservoir_g)]
            if finite.size:
                reservoir_max_g = max(reservoir_max_g, float(np.max(finite)))
            any_reservoir = True
    if any_reservoir:
        ax6r.set_ylabel("RIGHT axis: wall-water reservoir (g)")
        ax6r.set_ylim(0.0, max(1.0, 1.08 * reservoir_max_g))
    h1, l1 = ax[6].get_legend_handles_labels()
    h2, l2 = ax6r.get_legend_handles_labels()
    if any_loss or any_reservoir:
        ax[6].legend(h1 + h2, l1 + l2, fontsize=7.2, loc="best")
    ax[6].grid(alpha=0.3)

    # (h) Summary panel
    ax[7].axis("off")
    sat_time = np.nan
    if cfg.INITIAL_RH_METHOD == "cloud_onset":
        if saturation_time_s is None:
            sat_time = _model_saturation_time(
                exp, meta.CLOUD_ONSET.get(exp, onset), met["Time"]
            )
        else:
            sat_time = float(saturation_time_s)
    sat_label = (
        f"{sat_time/60:.2f} min"
        if np.isfinite(sat_time)
        else "n/a (dewpoint initialisation)"
    )
    summary = (
        f"Cloud comparison: {window[0]/60:.2f}-{window[1]/60:.2f} min\n"
        f"Observed onset marker: {onset/60:.2f} min\n"
        f"Model saturation target: {sat_label}\n"
        f"ql comparison: {ql_mode_label}\n"
        f"Diagnostic ql lag (not automatically applied): {diagnosed_lag:+.0f} s "
        f"(corr={lag_corr:.3f})\n"
        f"Applied WELAS time shift: {applied_obs_shift_s:+.0f} s "
        f"({'ON' if shift_cloud_measurements else 'OFF'})\n"
        f"Applied WELAS shift: {applied_obs_shift_s:+.0f} s\n"
        f"selected ql NRMSE: {ql_abs['nrmse']:.3f}\n"
        f"selected ql NRMSE lagged: {ql_lag['nrmse']:.3f}\n"
        f"total ql NRMSE: {ql_total_abs['nrmse']:.3f}\n"
        f">{cfg.SINGLE_COMPARE_DROP_MIN_UM:g}um ql NRMSE: {ql_above_abs['nrmse']:.3f}\n"
        f"Nd(>{cfg.SINGLE_COMPARE_DROP_MIN_UM:g} um) NRMSE: {nd_score['nrmse']:.3f}\n"
        f"Total Deff NRMSE: {deff_score['nrmse']:.3f}\n"
        f"Total rel. dispersion NRMSE: {rel_score['nrmse']:.3f}\n"
        f"Integral selected ql M/O: {int_ql_model:.4g} / {int_ql_obs:.4g} kg s kg-1\n"
        f"Integral total BMM ql: {int_ql_model_total:.4g} kg s kg-1"
    )
    wall_liq = np.asarray(
        model.get("chamber_wall_liquid_water", np.full_like(model["time"], np.nan)),
        dtype=float,
    )
    wall_ice = np.asarray(
        model.get("chamber_wall_ice_water", np.full_like(model["time"], np.nan)),
        dtype=float,
    )
    wall_total = wall_liq + wall_ice
    if np.any(np.isfinite(wall_total)):
        final_liq_g = (
            1.0e3 * float(wall_liq[np.flatnonzero(np.isfinite(wall_liq))[-1]])
            if np.any(np.isfinite(wall_liq)) else np.nan
        )
        final_ice_g = (
            1.0e3 * float(wall_ice[np.flatnonzero(np.isfinite(wall_ice))[-1]])
            if np.any(np.isfinite(wall_ice)) else np.nan
        )
        max_total_g = 1.0e3 * float(np.nanmax(wall_total))
        summary += (
            f"\nFinal wall liquid/ice: {final_liq_g:.3g} / {final_ice_g:.3g} g"
            f"\nMax total wall reservoir: {max_total_g:.3g} g"
        )
    ax[7].set_title("(h) Comparison summary")
    ax[7].text(0.02, 0.98, summary, va="top", family="monospace", fontsize=9)

    for a in ax[:7]:
        a.axvline(
            onset / 60.0, color="0.45", linestyle=":", linewidth=1.0,
            label=None,
        )
        a.axvspan(
            window[0] / 60.0, window[1] / 60.0,
            color="0.85", alpha=0.18,
        )
        a.set_xlabel("Experiment time (min)")
    ax3i.axvline(onset / 60.0, color="0.45", linestyle=":", linewidth=1.0)

    fig.suptitle(
        f"{exp}: BMM vs iSKYLAB observations\n"
        "Black/grey = observations; blue/green = BMM liquid/total; orange = ice",
        fontsize=14,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.965))

    # Dedicated wall-vapour diagnostic.  This is especially useful for the
    # physical mass-transfer closure (mode 2): RHwall shows the thermodynamic
    # driving force and Jv shows the actual signed area-mean mass flux.
    wall_fig = None
    if "chamber_wall_rh" in model or "chamber_wall_vapour_flux" in model:
        wall_fig, wax = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
        if "chamber_wall_rh" in model:
            wall_rh = np.asarray(model["chamber_wall_rh"], dtype=float)
            wax[0].plot(model["time"] / 60.0, wall_rh, linewidth=1.8)
            wax[0].axhline(1.0, linestyle="--", linewidth=1.0)
            wax[0].set_ylabel(r"$RH_{wall}$")
            wax[0].set_title("Vapour relative humidity with respect to measured wall temperature")
            wax[0].grid(alpha=0.3)
        if "chamber_wall_vapour_flux" in model:
            wall_flux = 1.0e6 * np.asarray(model["chamber_wall_vapour_flux"], dtype=float)
            wax[1].plot(model["time"] / 60.0, wall_flux, linewidth=1.8)
            wax[1].axhline(0.0, linestyle="--", linewidth=1.0)
            wax[1].set_ylabel(r"$J_v$ (mg m$^{-2}$ s$^{-1}$)")
            wax[1].set_title("Signed wall vapour flux (+ wall to air; - air to wall)")
            wax[1].grid(alpha=0.3)
        wax[1].set_xlabel("Experiment time (min)")
        wall_fig.suptitle(f"{exp}: chamber wall-vapour exchange")
        wall_fig.tight_layout()

    # Direct WELAS/model PSD comparison on a deliberately coarser common
    # logarithmic grid.  The native WELAS grid is much finer than the number
    # of independent moving BMM populations, so putting model delta-populations
    # into every native optical channel produces a sparse/striped panel.
    compare_edges = _common_comparison_psd_edges(obs["Dp_edges_um"])
    obs_psd_compare = _rebin_dndlog10d(
        obs["psd"], obs["Dp_edges_um"], compare_edges
    )
    model_psd_liq = _model_psd_on_fixed_grid(
        model, compare_edges, diameter_key="dwet", number_key="nwat"
    )
    model_psd_ice = _model_psd_on_fixed_grid(
        model, compare_edges, diameter_key="dmaxice", number_key="nicem"
    )
    if model_psd_liq is None:
        model_psd = model_psd_ice
    elif model_psd_ice is None:
        model_psd = model_psd_liq
    else:
        model_psd = model_psd_liq + model_psd_ice
    psd_fig = None
    if model_psd is not None:
        psd_fig, pax = plt.subplots(2, 1, figsize=(13, 9), sharex=True, sharey=True)
        # Use the observed WELAS concentration range to define the shared
        # logarithmic colour scale.  The model therefore cannot stretch or
        # compress the colour limits; identical colours in the two panels
        # correspond to identical dN/dlog10D values on the observational scale.
        welas_positive = obs_psd_compare[
            np.isfinite(obs_psd_compare) & (obs_psd_compare > 0.0)
        ]
        if welas_positive.size:
            vmin = max(float(np.nanpercentile(welas_positive, 2.0)), 1.0e-4)
            vmax = max(
                float(np.nanpercentile(welas_positive, 99.5)),
                10.0 * vmin,
            )
        else:
            vmin, vmax = 1.0e-4, 1.0
        norm = LogNorm(vmin=vmin, vmax=vmax)
        pcm0 = pax[0].pcolormesh(
            _time_edges(obs["time"]) / 60.0, compare_edges, obs_psd_compare.T,
            shading="auto", norm=norm,
        )
        pax[0].set_title(
            f"Observed WELAS rebinned to {len(compare_edges)-1} common log-D bins"
        )
        pcm1 = pax[1].pcolormesh(
            _time_edges(model["time"]) / 60.0, compare_edges, model_psd.T,
            shading="auto", norm=norm,
        )
        
        xmin = min(np.nanmin(obs["time"]), np.nanmin(model["time"])) / 60.0
        xmax = max(np.nanmax(obs["time"]),np.nanmax(model["time"])) / 60.0
        
        for a in pax:
            a.set_xlim(xmin, xmax)
            
        pax[1].set_title(
            f"BMM liquid + ice rebinned to the same {len(compare_edges)-1} log-D bins"
        )
        for a in pax:
            a.set_yscale("log")
            a.set_ylabel(r"Wet diameter ($\mu$m)")
            a.axvline(onset / 60.0, color="w", linestyle=":", linewidth=1)
        pax[1].set_xlabel("Time (min)")
        psd_fig.colorbar(pcm1, ax=pax, label=r"dN/dlog$_{10}$D (cm$^{-3}$)")
        psd_fig.suptitle(
            f"{exp}: observed merged and BMM liquid+ice size distributions "
            "(common coarse grid; BMM ice uses Dmax)"
        )

    # Model-only full PSD diagnostic.  The warm panel includes every nwat
    # particle (unactivated aerosol + haze + activated liquid) rather than the
    # activated-only nliq subset.  Ice is shown independently when available.
    n_psd = int(cfg.SINGLE_MODEL_PSD_NBINS)
    if n_psd < 10:
        raise ValueError("SINGLE_MODEL_PSD_NBINS must be >= 10")
    warm_edges = np.logspace(
        np.log10(cfg.SINGLE_MODEL_WARM_PSD_MIN_UM),
        np.log10(cfg.SINGLE_MODEL_WARM_PSD_MAX_UM),
        n_psd + 1,
    )
    warm_psd = _model_psd_on_fixed_grid(
        model, warm_edges, diameter_key="dwet", number_key="nwat"
    )
    ice_edges = np.logspace(
        np.log10(cfg.SINGLE_MODEL_ICE_PSD_MIN_UM),
        np.log10(cfg.SINGLE_MODEL_ICE_PSD_MAX_UM),
        n_psd + 1,
    )
    ice_psd = _model_psd_on_fixed_grid(
        model, ice_edges, diameter_key="dmaxice", number_key="nicem"
    )
    full_psd_fig = None
    if warm_psd is not None:
        have_ice = ice_psd is not None and np.any(np.isfinite(ice_psd) & (ice_psd > 0.0))
        nrows = 2 if have_ice else 1
        full_psd_fig, faxes = plt.subplots(nrows, 1, figsize=(13, 5.0 * nrows), squeeze=False)
        faxes = faxes.ravel()

        def _psd_norm(arr):
            pos = arr[np.isfinite(arr) & (arr > 0.0)]
            if pos.size == 0:
                return LogNorm(vmin=1.0e-4, vmax=1.0)
            vmin = max(float(np.nanpercentile(pos, 1.0)), 1.0e-6)
            vmax = max(float(np.nanpercentile(pos, 99.5)), 10.0 * vmin)
            return LogNorm(vmin=vmin, vmax=vmax)

        wpcm = faxes[0].pcolormesh(
            _time_edges(model["time"]) / 60.0,
            warm_edges,
            warm_psd.T,
            shading="auto",
            norm=_psd_norm(warm_psd),
        )
        faxes[0].set_title("BMM complete warm PSD: aerosol + haze + activated liquid")
        faxes[0].set_yscale("log")
        faxes[0].set_ylabel(r"Wet diameter ($\mu$m)")
        faxes[0].axvline(onset / 60.0, color="w", linestyle=":", linewidth=1)
        full_psd_fig.colorbar(wpcm, ax=faxes[0], label=r"dN/dlog$_{10}$D (cm$^{-3}$)")

        if have_ice:
            ipcm = faxes[1].pcolormesh(
                _time_edges(model["time"]) / 60.0,
                ice_edges,
                ice_psd.T,
                shading="auto",
                norm=_psd_norm(ice_psd),
            )
            faxes[1].set_title(r"BMM ice PSD ($D_{max}$)")
            faxes[1].set_yscale("log")
            faxes[1].set_ylabel(r"$D_{max}$ ($\mu$m)")
            faxes[1].axvline(onset / 60.0, color="w", linestyle=":", linewidth=1)
            full_psd_fig.colorbar(ipcm, ax=faxes[1], label=r"dN/dlog$_{10}$D (cm$^{-3}$)")

        faxes[-1].set_xlabel("Time (min)")
        full_psd_fig.suptitle(f"{exp}: complete BMM particle size distributions")
        full_psd_fig.tight_layout()

    metrics = {
        "experiment": exp,
        "window_start_s": window[0], "window_end_s": window[1],
        "model_saturation_target_s": sat_time,
        "ql_compare_mode": ql_mode,
        "best_ql_lag_s": diagnosed_lag, "best_ql_lag_corr": lag_corr,
        "applied_welas_time_shift_s": applied_obs_shift_s,
        "ql_nrmse_absolute": ql_abs["nrmse"], "ql_nrmse_lagged": ql_lag["nrmse"],
        "ql_total_nrmse_absolute": ql_total_abs["nrmse"],
        "ql_above_min_nrmse_absolute": ql_above_abs["nrmse"],
        "ndrop_psd_nrmse": nd_score["nrmse"], "deff_nrmse": deff_score["nrmse"],
        "rel_disp_nrmse": rel_score["nrmse"],
        "integrated_ql_model_kg_s_kg": int_ql_model,
        "integrated_ql_obs_kg_s_kg": int_ql_obs,
        "integrated_ql_model_total_kg_s_kg": int_ql_model_total,
        "integrated_ql_obs_total_kg_s_kg": int_ql_obs_total,
        "integrated_ql_model_above_min_kg_s_kg": int_ql_model_above,
        "integrated_ql_obs_above_min_kg_s_kg": int_ql_obs_above,
    }
    if cfg.SAVE_SINGLE_COMPARISON:
        outdir = cfg.OUTPUT_ROOT / "single_comparison"
        outdir.mkdir(parents=True, exist_ok=True)
        fig.savefig(outdir / f"comparison-{exp}.png", dpi=180)
        if psd_fig is not None:
            psd_fig.savefig(outdir / f"psd-{exp}.png", dpi=180)
        if full_psd_fig is not None:
            full_psd_fig.savefig(outdir / f"model-full-psd-{exp}.png", dpi=180)
        if wall_fig is not None:
            wall_fig.savefig(outdir / f"wall-vapour-{exp}.png", dpi=180)
        with (outdir / f"metrics-{exp}.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["metric", "value"])
            writer.writerows(metrics.items())
        print(f"Single-experiment diagnostics written to {outdir}")
    if show:
        plt.show()
        if not plt.isinteractive():
            plt.close("all")
    return fig, psd_fig, full_psd_fig, metrics

def main(group=THIS_RUN, run_model=RUN_MODEL, do_analysis=DO_ANALYSIS, do_plot=DO_PLOT,
         experiment=None, winit_single=1.3, saturation_time_min=None, sce_bins=None,
         shift_cloud_measurements=False, bl_tau_s=None, bl_temp_offset_k=None,
         bl_wall_water_mode=None, wall_liquid_water_init_kg=None,
         wall_ice_water_init_kg=None, wall_water_efficiency=None,
         wall_vapour_transfer_velocity_ms=None, bl_evap_size_exp=None,
         interactive_plots=False):
    if interactive_plots:
        plt.ion()
    else:
        plt.ioff()

    if not READ_DATA:
        raise RuntimeError("READ_DATA=False is no longer supported without supplying a cached data dictionary")

    # Resolve the requested target *before* reading data.  Single-experiment
    # mode should never attempt to open files belonging to unrelated cases.
    exp = _normalise_experiment_name(experiment) if experiment is not None else None
    if exp is not None:
        requested_experiments = [exp]
    else:
        if group < 0 or group >= len(meta.BATCH_GROUPS):
            raise ValueError(f"group must be between 0 and {len(meta.BATCH_GROUPS)-1}")
        requested_experiments = meta.BATCH_GROUPS[group]

    data = load_all_data(requested_experiments)
    _write_forcing_smoothing_diagnostics(data)

    if saturation_time_min is not None and exp is None:
        raise ValueError("--cloud-formation-time-min/--saturation-time-min is only valid with --experiment")
    if shift_cloud_measurements and exp is None:
        raise ValueError("--shift-cloud-measurements is currently only valid with --experiment")
    if bl_tau_s is not None and float(bl_tau_s) <= 0.0:
        raise ValueError("--bl-tau-s must be > 0")
    if bl_wall_water_mode is not None and int(bl_wall_water_mode) not in (0, 1, 2):
        raise ValueError("--bl-wall-water-mode must be 0, 1 or 2")
    if wall_liquid_water_init_kg is not None and float(wall_liquid_water_init_kg) < 0.0:
        raise ValueError("--wall-liquid-water-init-kg must be >= 0")
    if wall_ice_water_init_kg is not None and float(wall_ice_water_init_kg) < 0.0:
        raise ValueError("--wall-ice-water-init-kg must be >= 0")
    if wall_water_efficiency is not None and not 0.0 <= float(wall_water_efficiency) <= 1.0:
        raise ValueError("--wall-water-efficiency must be between 0 and 1")
    if wall_vapour_transfer_velocity_ms is not None and float(wall_vapour_transfer_velocity_ms) < 0.0:
        raise ValueError("--wall-vapour-transfer-velocity-ms must be >= 0")
    if bl_evap_size_exp is not None and float(bl_evap_size_exp) < 0.0:
        raise ValueError("--bl-evap-size-exp must be >= 0")
    if saturation_time_min is not None and cfg.INITIAL_RH_METHOD != "cloud_onset":
        raise ValueError(
            "--cloud-formation-time-min/--saturation-time-min requires INITIAL_RH_METHOD='cloud_onset'; "
            "it has no effect with dewpoint initialisation"
        )
    saturation_overrides = {}
    if exp is not None and saturation_time_min is not None:
        saturation_overrides[exp] = float(saturation_time_min) * 60.0
    states = build_initial_state(data, saturation_time_overrides_s=saturation_overrides)

    if experiment is not None:
        if exp not in states:
            raise KeyError(f"No initial-state metadata/data available for {exp}")
        model_file = cfg.OUTPUT_ROOT / "single_comparison" / f"output-{exp}.nc"
        if run_model:
            model_file = run_single_experiment(
                exp, states[exp], data, winit=winit_single, sce_bins=sce_bins,
                bl_tau_s=bl_tau_s, bl_temp_offset_k=bl_temp_offset_k,
                bl_wall_water_mode=bl_wall_water_mode,
                wall_liquid_water_init_kg=wall_liquid_water_init_kg,
                wall_ice_water_init_kg=wall_ice_water_init_kg,
                wall_water_efficiency=wall_water_efficiency,
                wall_vapour_transfer_velocity_ms=wall_vapour_transfer_velocity_ms,
                bl_evap_size_exp=bl_evap_size_exp,
            )
        elif not model_file.exists():
            raise FileNotFoundError(f"Existing single-run output not found: {model_file}")
        if do_analysis:
            return analyse_single_experiment(
                exp, data, model_file, show=do_plot,
                saturation_time_s=states[exp]["saturation_time"],
                shift_cloud_measurements=shift_cloud_measurements,
            )
        return model_file

    batch_sims = meta.BATCH_GROUPS[group]
    winit = meta.GROUP_UPDRAFT[group]

    vals1 = vals2 = None
    if run_model:
        vals1, vals2 = run_batch(
            batch_sims, states, data, winit, sce_bins=sce_bins,
            bl_tau_s=bl_tau_s, bl_temp_offset_k=bl_temp_offset_k,
            bl_wall_water_mode=bl_wall_water_mode,
            wall_liquid_water_init_kg=wall_liquid_water_init_kg,
            wall_ice_water_init_kg=wall_ice_water_init_kg,
            wall_water_efficiency=wall_water_efficiency,
            wall_vapour_transfer_velocity_ms=wall_vapour_transfer_velocity_ms,
        )

    if do_analysis:
        if vals1 is None or vals2 is None:
            vals1 = np.zeros(len(batch_sims))
            vals2 = np.zeros(len(batch_sims))
            for i, exp in enumerate(batch_sims):
                aer = _aerosol_mode_arrays(exp, states[exp], data)
                vals1[i] = np.sum(aer["n1"])
                vals2[i] = np.sum(aer["n2"])
        fig = analyse_batch(batch_sims, states, data, vals1, vals2, meta.GROUP_TYPE[group])
        if do_plot:
            plt.show()
            if not plt.isinteractive():
                plt.close("all")
        return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    target = parser.add_mutually_exclusive_group()
    target.add_argument("--group", type=int, default=None, help="index in experiment_metadata.BATCH_GROUPS")
    target.add_argument("--experiment", help="single experiment, e.g. Exp005 or 5")
    parser.add_argument("--winit", type=float, default=1.3, help="initial/updraft value used for a single run")
    parser.add_argument(
        "--cloud-formation-time-min", "--saturation-time-min",
        dest="saturation_time_min",
        type=float,
        default=None,
        help=(
            "single-run model cloud/saturation target in minutes from experiment "
            "start; sets the initial vapour amount so the measured P/T trajectory "
            "would reach liquid saturation at this time in the absence of water "
            "exchange. The observed WELAS timing is not moved."
        ),
    )
    parser.add_argument(
        "--shift-cloud-measurements",
        action="store_true",
        help=(
            "shift the complete WELAS observation time coordinate by the diagnosed "
            "ql lag before comparison; applies consistently to ql, number, Deff, "
            "dispersion, ice number and the full observed PSD (default: off)"
        ),
    )
    parser.add_argument(
        "--bl-tau-s",
        type=float,
        default=None,
        help=(
            "override chamber_bl_tau (s) in generated BMM namelists for this run; "
            "if omitted use iskylab_config.CHAMBER_BL_TAU"
        ),
    )
    parser.add_argument(
        "--bl-evap-size-exp", "--bl-inhom-size-exp",
        dest="bl_evap_size_exp",
        type=float,
        default=None,
        help=(
            "common chamber-BL evaporation size exponent p. Homogeneous mode "
            "uses p to control size/mass shrinkage; inhomogeneous mode uses p "
            "to control complete-particle selection. p=2 is D2-based in both. "
            "--bl-inhom-size-exp is retained as a deprecated alias."
        ),
    )
    parser.add_argument(
        "--bl-temp-offset-k",
        type=float,
        default=None,
        help=(
            "override chamber_bl_temp_offset (K) in generated BMM namelists for "
            "this run; if omitted use iskylab_config.CHAMBER_BL_TEMP_OFFSET"
        ),
    )
    parser.add_argument(
        "--bl-wall-water-mode",
        type=int,
        choices=(0, 1, 2),
        default=None,
        help=(
            "wall-vapour closure: 0=legacy saturation cap; "
            "1=finite reservoir with fractional relaxation coupled to BL tau; "
            "2=finite reservoir with physical vapour mass-transfer velocity "
            "independent of BL tau"
        ),
    )
    parser.add_argument(
        "--wall-vapour-transfer-velocity-ms",
        type=float,
        default=None,
        help=(
            "physical chamber-wall vapour mass-transfer velocity km in m s-1, "
            "used by --bl-wall-water-mode 2; if omitted use "
            "iskylab_config.CHAMBER_WALL_VAPOUR_TRANSFER_VELOCITY"
        ),
    )
    parser.add_argument(
        "--wall-liquid-water-init-kg",
        type=float,
        default=None,
        help=(
            "initial physical liquid-water mass stored on the chamber wall in kg; "
            "default is 0 (dry wall)"
        ),
    )
    parser.add_argument(
        "--wall-ice-water-init-kg",
        type=float,
        default=None,
        help=(
            "initial physical ice/frost mass stored on the chamber wall in kg; "
            "default is 0 (no initial frost)"
        ),
    )
    parser.add_argument(
        "--wall-water-efficiency",
        type=float,
        default=None,
        help=(
            "fractional equilibration of a BL wall encounter with wall equilibrium, "
            "0=no vapour exchange and 1=full equilibration; evaporation/sublimation "
            "is always capped by the prognostic stored wall-water mass"
        ),
    )
    parser.add_argument(
        "--sce-bins", "--scebins",
        dest="sce_bins",
        type=int,
        default=None,
        help=(
            "set n_binsc in a temporary copy of BMM_MODEL_FOLDER/sce/namelist.in "
            "and point the generated BMM namelist scefile at that copy; "
            "the repository SCE namelist is never modified"
        ),
    )
    parser.add_argument(
        "--plt-ion",
        action="store_true",
        help=(
            "enable matplotlib interactive mode with plt.ion(); default is off. "
            "With the default non-interactive mode, displayed figures are closed "
            "after plt.show() so batch scripts can continue to the next run."
        ),
    )
    parser.add_argument("--no-run", action="store_true", help="analyse existing output without running BMM")
    parser.add_argument("--no-analysis", action="store_true", help="generate/run namelists only")
    parser.add_argument("--no-plot", action="store_true", help="do not display plots (saved diagnostics still written)")
    args = parser.parse_args()
    selected_group = THIS_RUN if args.group is None else args.group
    main(
        group=selected_group,
        experiment=args.experiment,
        winit_single=args.winit,
        saturation_time_min=args.saturation_time_min,
        sce_bins=args.sce_bins,
        shift_cloud_measurements=args.shift_cloud_measurements,
        bl_tau_s=args.bl_tau_s,
        bl_temp_offset_k=args.bl_temp_offset_k,
        bl_evap_size_exp=args.bl_evap_size_exp,
        bl_wall_water_mode=args.bl_wall_water_mode,
        wall_liquid_water_init_kg=args.wall_liquid_water_init_kg,
        wall_ice_water_init_kg=args.wall_ice_water_init_kg,
        wall_water_efficiency=args.wall_water_efficiency,
        wall_vapour_transfer_velocity_ms=args.wall_vapour_transfer_velocity_ms,
        run_model=not args.no_run,
        do_analysis=not args.no_analysis,
        do_plot=not args.no_plot,
        interactive_plots=args.plt_ion,
    )
