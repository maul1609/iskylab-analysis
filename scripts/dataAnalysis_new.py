"""Initialise, batch-run and analyse iSKYLAB experiments with the BMM.

This is the maintained iSKYLAB batch driver.  It reads the supplied chamber
measurements and initial aerosol PSDs, creates one BMM namelist per experiment,
optionally runs ``main.exe``, and compares model/observed cloud-droplet number.

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
* The measured ``Tww_mean`` series can optionally be written as
  ``wall_temp_chamber`` for the revised boundary-layer model.
* Aerosol diameters from the PNSD fitting code are in micrometres internally
  and are converted to metres exactly once when the namelist is written.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter

import experiment_metadata as meta
import iskylab_config as cfg
import readMeteoCPC
import readOPC_Merged
import readPNSD_Mrg_new
import svp
from namelist_utils import read_text, replace_group, set_array, set_value, write_text

R_GAS = 8.314
M_WATER = 18.0e-3
M_DRY_AIR = 28.96e-3
R_D = R_GAS / M_DRY_AIR
R_V = R_GAS / M_WATER
EPSILON = R_D / R_V


def _savgol_time_series(time, values, window_seconds, polyorder, *, name="series"):
    """Return a Savitzky-Golay-smoothed series using a time-based window.

    The chamber records are nominally close to uniform in time, but this helper
    does not assume that every experiment has the same sample interval.  It
    infers the median positive ``dt`` and converts ``window_seconds`` to an odd
    sample count.  If the timestamps are noticeably irregular, the observations
    are first interpolated onto a uniform median-dt grid, filtered there, and
    then interpolated back to the original timestamps.

    End points use SciPy's ``mode="interp"`` so they are not padded with an
    artificial constant/zero value.  Windows that are too short for the chosen
    polynomial order are rejected rather than silently applying a different
    filter.
    """
    time = np.asarray(time, dtype=float)
    values = np.asarray(values, dtype=float)
    if time.ndim != 1 or values.ndim != 1 or len(time) != len(values):
        raise ValueError(f"{name}: time and values must be one-dimensional and equal length")
    if len(time) < polyorder + 3:
        raise ValueError(f"{name}: too few samples for Savitzky-Golay smoothing")
    if window_seconds <= 0.0:
        raise ValueError(f"{name}: smoothing window must be positive")
    if polyorder < 0:
        raise ValueError(f"{name}: polynomial order must be non-negative")

    if not np.all(np.isfinite(time)):
        raise ValueError(f"{name}: time coordinate contains non-finite values")
    dt_all = np.diff(time)
    if np.any(dt_all <= 0.0):
        raise ValueError(f"{name}: time coordinate must be strictly increasing")
    dt_good = dt_all[np.isfinite(dt_all)]
    if dt_good.size == 0:
        raise ValueError(f"{name}: no valid time increments")
    dt = float(np.median(dt_good))

    # Round the requested physical window to the closest sensible odd number
    # of samples.  Require at least polyorder+2 points and then make it odd.
    nwin = int(round(window_seconds / dt))
    nwin = max(nwin, polyorder + 2)
    if nwin % 2 == 0:
        nwin += 1
    if nwin > len(time):
        nwin = len(time) if len(time) % 2 == 1 else len(time) - 1
    if nwin <= polyorder:
        raise ValueError(
            f"{name}: effective window ({nwin} samples) must exceed polyorder={polyorder}"
        )

    # Determine whether direct filtering is justified.  A 5% relative spread
    # in dt is small enough that treating the samples as uniform is harmless;
    # otherwise filter on an explicitly uniform grid.
    rel_irregularity = float(np.max(np.abs(dt_good - dt)) / dt)
    if rel_irregularity <= 0.05:
        smoothed = savgol_filter(values, nwin, polyorder, mode="interp")
    else:
        uniform_time = np.arange(time[0], time[-1] + 0.5 * dt, dt)
        uniform_values = np.interp(uniform_time, time, values)
        nwin_uniform = int(round(window_seconds / dt))
        nwin_uniform = max(nwin_uniform, polyorder + 2)
        if nwin_uniform % 2 == 0:
            nwin_uniform += 1
        if nwin_uniform > len(uniform_time):
            nwin_uniform = len(uniform_time) if len(uniform_time) % 2 == 1 else len(uniform_time) - 1
        if nwin_uniform <= polyorder:
            raise ValueError(f"{name}: effective uniform-grid window is too short")
        uniform_smooth = savgol_filter(
            uniform_values, nwin_uniform, polyorder, mode="interp"
        )
        smoothed = np.interp(time, uniform_time, uniform_smooth)

    return np.asarray(smoothed, dtype=float), {
        "median_dt_s": dt,
        "window_samples": int(nwin),
        "window_seconds_effective": float(nwin * dt),
        "max_dt_irregularity_fraction": rel_irregularity,
    }


def _smooth_chamber_forcing_block(block, *, experiment):
    """Preserve raw forcing observations and replace working P/T fields by smooth data."""
    time = np.asarray(block["Time"], dtype=float)
    settings = [
        ("Tgw mean", cfg.CHAMBER_SMOOTH_TEMP_WINDOW, "gas temperature"),
        ("Tww mean", cfg.CHAMBER_SMOOTH_WALL_TEMP_WINDOW, "wall temperature"),
        ("Pressure", cfg.CHAMBER_SMOOTH_PRESSURE_WINDOW, "pressure"),
    ]
    block.setdefault("forcing_smoothing", {})

    for key, window, label in settings:
        raw_key = f"{key}_raw"
        # Store the originally read observations, including any NaNs, once.
        if raw_key not in block:
            block[raw_key] = np.asarray(block[key], dtype=float).copy()

        # Fill missing observations before filtering.  This filled series is
        # retained too because it makes diagnosing gaps versus filtering easy.
        filled = _fill_nan_linear(time, block[key], name=f"{experiment} {label}")
        block[f"{key}_filled"] = filled.copy()

        if cfg.SMOOTH_CHAMBER_FORCING:
            smooth, info = _savgol_time_series(
                time, filled, window, cfg.CHAMBER_SMOOTH_POLYORDER,
                name=f"{experiment} {label}",
            )
            block[key] = smooth
            correction = smooth - filled
            info.update(
                {
                    "rms_correction": float(np.sqrt(np.mean(correction**2))),
                    "max_abs_correction": float(np.max(np.abs(correction))),
                    "mean_correction": float(np.mean(correction)),
                }
            )
            block["forcing_smoothing"][key] = info
        else:
            block[key] = filled
            block["forcing_smoothing"][key] = {
                "median_dt_s": float(np.median(np.diff(time))),
                "window_samples": 1,
                "window_seconds_effective": 0.0,
                "max_dt_irregularity_fraction": 0.0,
                "rms_correction": 0.0,
                "max_abs_correction": 0.0,
                "mean_correction": 0.0,
            }


def _write_forcing_smoothing_diagnostics(data):
    """Write raw-vs-smoothed forcing plots and a compact CSV summary."""
    if not cfg.SAVE_FORCING_SMOOTHING_DIAGNOSTICS:
        return

    outdir = cfg.OUTPUT_ROOT / "forcing_smoothing"
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []

    for key in readMeteoCPC.metStr:
        if key not in data:
            continue
        met = data[key]
        exp = key.split("-")[-1]
        time = np.asarray(met["Time"], dtype=float)

        for field in ("Tgw mean", "Tww mean", "Pressure"):
            info = met.get("forcing_smoothing", {}).get(field, {})
            rows.append(
                {
                    "experiment": exp,
                    "field": field,
                    "median_dt_s": info.get("median_dt_s", np.nan),
                    "window_samples": info.get("window_samples", np.nan),
                    "window_seconds_effective": info.get("window_seconds_effective", np.nan),
                    "max_dt_irregularity_fraction": info.get("max_dt_irregularity_fraction", np.nan),
                    "rms_correction": info.get("rms_correction", np.nan),
                    "max_abs_correction": info.get("max_abs_correction", np.nan),
                    "mean_correction": info.get("mean_correction", np.nan),
                }
            )

        fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

        axes[0].plot(time, met["Tgw mean_raw"], alpha=0.45, label="gas raw")
        axes[0].plot(time, met["Tgw mean"], label="gas smoothed")
        axes[0].plot(time, met["Tww mean_raw"], alpha=0.45, label="wall raw")
        axes[0].plot(time, met["Tww mean"], label="wall smoothed")
        axes[0].set_ylabel("Temperature (degC)")
        axes[0].legend(ncol=2)
        axes[0].grid()

        delta_raw = np.asarray(met["Tww mean_raw"]) - np.asarray(met["Tgw mean_raw"])
        delta_smooth = np.asarray(met["Tww mean"]) - np.asarray(met["Tgw mean"])
        axes[1].plot(time, delta_raw, alpha=0.45, label="raw")
        axes[1].plot(time, delta_smooth, label="smoothed")
        axes[1].axhline(0.0, linewidth=0.8)
        axes[1].set_ylabel(r"$T_{wall}-T_{gas}$ (K)")
        axes[1].legend()
        axes[1].grid()

        axes[2].plot(time, met["Pressure_raw"], alpha=0.45, label="raw")
        axes[2].plot(time, met["Pressure"], label="smoothed")
        axes[2].set_ylabel("Pressure (hPa)")
        axes[2].set_xlabel("Experiment time (s)")
        axes[2].legend()
        axes[2].grid()

        state = "enabled" if cfg.SMOOTH_CHAMBER_FORCING else "disabled"
        fig.suptitle(f"{exp}: chamber forcing smoothing ({state})")
        fig.tight_layout()
        fig.savefig(outdir / f"forcing-smoothing-{exp}.png", dpi=160)
        plt.close(fig)

    if rows:
        csv_path = outdir / "forcing-smoothing-summary.csv"
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


def load_all_data():
    """Read meteo/CPC, OPC and initial-PNSD files into one dictionary."""
    data = {}

    for i, key in enumerate(readMeteoCPC.metStr):
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
        data[key] = readOPC_Merged.readData(readThis=i, opcStr=key)[key]

    for i, key in enumerate(readPNSD_Mrg_new.npsdStr):
        data[key] = readPNSD_Mrg_new.readData(readThis=i, npsdStr=key)[key]

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
    for key in readMeteoCPC.metStr:
        exp = key.split("-")[-1]
        met = data[key]
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


def build_initial_state(data):
    """Calculate initial thermodynamics and aerosol normalisation metadata."""
    states = {}
    for exp, cloud_time in meta.CLOUD_ONSET.items():
        met_key = f"MeteoCPC-{exp}"
        if met_key not in data:
            continue
        met = data[met_key]
        t0 = _interp_at(met["Time"], met["Tgw mean"], 0.0, name=f"{exp} Tgas") + 273.15
        p0 = _interp_at(met["Time"], met["Pressure"], 0.0, name=f"{exp} pressure") * 100.0
        tc = _interp_at(met["Time"], met["Tgw mean"], cloud_time, name=f"{exp} cloud T") + 273.15
        pc = _interp_at(met["Time"], met["Pressure"], cloud_time, name=f"{exp} cloud p") * 100.0

        if cfg.INITIAL_RH_METHOD == "cloud_onset":
            # Historical initialisation: choose qv so that, in the absence of
            # pre-cloud water exchange, the measured P/T trajectory reaches
            # liquid saturation at the prescribed cloud-onset time.
            es_cloud = float(_svp_liq([tc])[0])
            qv0 = _mixing_ratio_from_vapour_pressure(es_cloud, pc)
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
            aerosol_time = cloud_time
        else:
            raise ValueError(f"Unknown AEROSOL_INIT_TIME={cfg.AEROSOL_INIT_TIME!r}")

        cpc_cm3 = _interp_at(met["Time"], met["CPC_TotBot"], aerosol_time, name=f"{exp} CPC")
        ta = _interp_at(met["Time"], met["Tgw mean"], aerosol_time, name=f"{exp} aerosol T") + 273.15
        pa = _interp_at(met["Time"], met["Pressure"], aerosol_time, name=f"{exp} aerosol p") * 100.0
        rho_d = pa / (R_D * ta)
        aerosol_number_mixing_ratio = cpc_cm3 * 1.0e6 / rho_d  # # kg-1 dry air

        states[exp] = {
            "cloud_time": cloud_time,
            "initT": t0,
            "initP": p0,
            "cloudT": tc,
            "cloudP": pc,
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


def chamber_spec_body(met):
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
    if cfg.CHAMBER_BL_TEMP_MODE == 1:
        blocks.extend(
            [
                "    ! Measured wall temperature (Tww_mean) for BL temp mode 1.",
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

    return {
        "n1": n1,
        "d1_m": d1_out_um * 1.0e-6,
        "sig1": sig1_out,
        "n2": n2,
        "d2_m": d2_out_um * 1.0e-6,
        "sig2": sig2_out,
        "density": density4,
        "kappa": kappa4,
    }


def make_namelist(exp, state, data, output_file, *, winit=1.3):
    """Return a complete BMM namelist for one iSKYLAB experiment."""
    text = read_text(cfg.BMM_TEMPLATE)
    met = data[f"MeteoCPC-{exp}"]
    n = len(met["Time"])
    runtime = min(float(np.nanmax(met["Time"])), 23.0 * 60.0)

    # General initial/run controls.
    text = set_value(text, "outputfile", str(output_file))
    text = set_value(text, "runtime", runtime)
    text = set_value(text, "tinit", state["initT"])
    text = set_value(text, "pinit", state["initP"])
    text = set_value(text, "rhinit", state["initRH"])
    text = set_value(text, "winit", float(winit))
    text = set_value(text, "fallout_flag", bool(cfg.FALLOUT_FLAG))
    text = set_value(text, "residence_depth", float(cfg.RESIDENCE_DEPTH))

    # Chamber observation/physics controls.
    force_p = cfg.FORCE_PRESSURE and not cfg.SYNTHETIC_UPDRAFT
    force_t = cfg.FORCE_TEMPERATURE and not cfg.SYNTHETIC_UPDRAFT
    force_q = cfg.FORCE_QTOT and not cfg.SYNTHETIC_UPDRAFT
    bl_mix = 0 if cfg.SYNTHETIC_UPDRAFT else cfg.CHAMBER_BL_MIX
    fan_loss = 0 if cfg.SYNTHETIC_UPDRAFT else cfg.CHAMBER_FAN_LOSS
    wall_loss = 0 if cfg.SYNTHETIC_UPDRAFT else cfg.CHAMBER_WALL_LOSS

    if force_q and not cfg.WRITE_QTOT_DATA:
        raise ValueError("FORCE_QTOT=True requires WRITE_QTOT_DATA=True")
    if cfg.CHAMBER_BL_TEMP_MODE == 1 and bl_mix and "Tww mean" not in met:
        raise ValueError("BL temp mode 1 requires measured Tww_mean")

    chamber_values = {
        "n_levels_c": n,
        "chamber_force_pressure": force_p,
        "chamber_force_temperature": force_t,
        "chamber_force_qtot": force_q,
        "chamber_bl_mix": bl_mix,
        "chamber_bl_tau": cfg.CHAMBER_BL_TAU,
        "chamber_bl_temp_mode": cfg.CHAMBER_BL_TEMP_MODE,
        "chamber_bl_temp_offset": cfg.CHAMBER_BL_TEMP_OFFSET,
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
    for name, value in chamber_values.items():
        text = set_value(text, name, value)
    text = replace_group(text, "chamber_spec", chamber_spec_body(met))

    # Measured/fitted aerosol initial conditions.
    aer = _aerosol_mode_arrays(exp, state, data)
    for name, values in [
        ("n_aer1(1:3,1:1)", aer["n1"]),
        ("d_aer1(1:3,1:1)", aer["d1_m"]),
        ("sig_aer1(1:3,1:1)", aer["sig1"]),
        ("n_aer1(1:3,2:2)", aer["n2"]),
        ("d_aer1(1:3,2:2)", aer["d2_m"]),
        ("sig_aer1(1:3,2:2)", aer["sig2"]),
        ("density_core1(1:4)", aer["density"]),
        ("kappa_core1(1:4)", aer["kappa"]),
    ]:
        text = set_array(text, name, values)

    return text, aer


def run_batch(batch_sims, states, data, winit):
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
        for nn, exp in enumerate(batch_sims):
            output_file = cfg.OUTPUT_ROOT / f"output{nn:03d}.nc"
            namelist_text, aer = make_namelist(
                exp, states[exp], data, output_file, winit=winit[nn]
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
            dwat = np.where(mwat > 0.0, (mwat / (np.pi / 6.0 * 1000.0)) ** (1.0 / 3.0), 0.0)
            conc = np.sum(np.where(dwat > 2.0e-6, nwat, 0.0), axis=(1, 2))
            time = np.asarray(nc["time"][:])
            p = np.asarray(nc["p"][:])
            t = np.asarray(nc["t"][:])
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
        # Recalculate drop number from bins >2 um using the historical OPC
        # calibration factor; retain the supplied ndrop field separately.
        diam_mask = opc["Dp"] > 2.0
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
    ax.plot(x, ndrop_model / ndrop_model[0], ".-", ms=10)
    ax.plot(x, cdnc_obs / cdnc_obs[0], ".-", ms=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Relative number of drops")
    ax.legend(["model", "data"])
    ax.grid()

    ax = fig.add_subplot(132)
    denom = vals1_cm3 + vals2_cm3
    ax.plot(x, ndrop_model / denom, ".-", ms=10)
    ax.plot(x, cdnc_obs / denom, ".-", ms=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Activated fraction")
    ax.legend(["model", "data"])
    ax.set_ylim((0, 1))
    ax.grid()

    ax = fig.add_subplot(133)
    ax.plot(x, ndrop_model, ".-", ms=10)
    ax.plot(x, cdnc_obs, ".-", ms=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Number of drops (cm$^{-3}$)")
    ax.legend(["model", "data"])
    ax.grid()
    fig.tight_layout()
    return fig


def main(group=THIS_RUN, run_model=RUN_MODEL, do_analysis=DO_ANALYSIS, do_plot=DO_PLOT):
    if group < 0 or group >= len(meta.BATCH_GROUPS):
        raise ValueError(f"group must be between 0 and {len(meta.BATCH_GROUPS)-1}")
    batch_sims = meta.BATCH_GROUPS[group]
    winit = meta.GROUP_UPDRAFT[group]

    if not READ_DATA:
        raise RuntimeError("READ_DATA=False is no longer supported without supplying a cached data dictionary")
    data = load_all_data()
    _write_forcing_smoothing_diagnostics(data)
    states = build_initial_state(data)

    vals1 = vals2 = None
    if run_model:
        vals1, vals2 = run_batch(batch_sims, states, data, winit)

    if do_analysis:
        if vals1 is None or vals2 is None:
            # Reconstruct the initial aerosol totals from the same configuration
            # when analysing existing output files without rerunning the model.
            vals1 = np.zeros(len(batch_sims))
            vals2 = np.zeros(len(batch_sims))
            for i, exp in enumerate(batch_sims):
                aer = _aerosol_mode_arrays(exp, states[exp], data)
                vals1[i] = np.sum(aer["n1"])
                vals2[i] = np.sum(aer["n2"])
        fig = analyse_batch(batch_sims, states, data, vals1, vals2, meta.GROUP_TYPE[group])
        if do_plot:
            plt.show()
        return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", type=int, default=THIS_RUN, help="index in experiment_metadata.BATCH_GROUPS")
    parser.add_argument("--no-run", action="store_true", help="analyse existing /tmp outputs without running BMM")
    parser.add_argument("--no-analysis", action="store_true", help="generate/run namelists only")
    parser.add_argument("--no-plot", action="store_true", help="do not display the summary plot")
    args = parser.parse_args()
    main(
        group=args.group,
        run_model=not args.no_run,
        do_analysis=not args.no_analysis,
        do_plot=not args.no_plot,
    )
