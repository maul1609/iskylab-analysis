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
* The measured ``Tww_mean`` series can optionally be written as
  ``wall_temp_chamber`` for the revised boundary-layer model.
* Aerosol diameters from the PNSD fitting code are in micrometres internally
  and are converted to metres exactly once when the namelist is written.
* ``--experiment Exp005`` runs/analyses one experiment against ql, Nd, Deff,
  relative dispersion and the full OPC size distribution.  WELAS and native
  BMM moving particles are conservatively rebinned onto a configurable common
  coarse log-D grid for a readable like-for-like PSD comparison.
* ``--saturation-time-min`` can move the model saturation target for a single
  run without moving the observed cloud-onset marker or comparison window.
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


def build_initial_state(data, *, saturation_time_overrides_s=None):
    """Calculate initial thermodynamics and aerosol normalisation metadata.

    ``saturation_time_overrides_s`` is intended for one-off single-run CLI
    sensitivity tests.  It changes only the model saturation target used to
    infer the initial vapour amount; observed cloud-onset metadata and plotting
    windows remain fixed.
    """
    states = {}
    saturation_time_overrides_s = saturation_time_overrides_s or {}
    for exp, observed_cloud_time in meta.CLOUD_ONSET.items():
        met_key = f"MeteoCPC-{exp}"
        if met_key not in data:
            continue
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
    runtime = min(float(np.nanmax(met["Time"])), 30.0 * 60.0)

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
        **psd,
    }


def _model_bulk_moments_above(model, dmin_um):
    """Return OPC-like BMM number/size moments above a wet-diameter threshold.

    These moments are calculated from the raw native `(nwat,dwet)` output, not
    from the model activation flag.  They therefore provide an apples-to-apples
    comparison with an optical instrument whose sample definition is diameter
    rather than Koehler activation state.
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
    w = np.where(d > dmin_um * 1.0e-6, n, 0.0)
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

    dmin = float(dmin_um) * 1.0e-6
    good = (
        np.isfinite(d) & np.isfinite(n) &
        (d > dmin) & (n > 0.0)
    )
    axes = tuple(range(1, n.ndim))
    rho_w = 1000.0
    return (
        rho_w * np.pi / 6.0
        * np.sum(np.where(good, n * d**3, 0.0), axis=axes)
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


def run_single_experiment(exp, state, data, *, winit=1.3):
    """Generate and run one named experiment, retaining its namelist and output."""
    cfg.OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    outdir = cfg.OUTPUT_ROOT / "single_comparison"
    outdir.mkdir(parents=True, exist_ok=True)
    if not cfg.BMM_EXECUTABLE.exists():
        raise FileNotFoundError(f"BMM executable not found: {cfg.BMM_EXECUTABLE}")
    output_file = outdir / f"output-{exp}.nc"
    namelist_text, _ = make_namelist(exp, state, data, output_file, winit=winit)
    namelist_file = outdir / f"namelist-{exp}.in"
    write_text(namelist_file, namelist_text)
    print(f"Single run: {exp} -> {output_file}")
    completed = subprocess.run(
        [str(cfg.BMM_EXECUTABLE), str(namelist_file)],
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


def analyse_single_experiment(exp, data, model_file, *, show=True, saturation_time_s=None):
    """Create time-series and PSD model/observation diagnostics for one experiment."""
    import readModel
    from matplotlib.colors import LogNorm

    model = readModel.readData(model_file, modelStr=exp)[exp]
    obs = _observed_bulk_series(exp, data)
    met = data[f"MeteoCPC-{exp}"]
    window = _cloud_comparison_window(exp, obs)
    onset = meta.CLOUD_ONSET.get(exp, window[0])

    nd_model = np.asarray(model["ndrop"]) * np.asarray(model["rhoa"]) / 1.0e6
    instrument_model = _model_bulk_moments_above(model, cfg.SINGLE_COMPARE_DROP_MIN_UM)
    nd_model_2um = instrument_model["number_cm3"]
    deff_model_um = instrument_model["deff_um"]
    rel_model = instrument_model["rel_disp"]
    deff_model_activated_um = np.asarray(model["deff"]) * 1.0e6
    rel_model_activated = np.asarray(model.get("rel_disp_liq", np.full_like(model["time"], np.nan)))

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

    lag, lag_corr = _best_ql_lag(
        model["time"], ql_model_compare, obs["time"], ql_obs_compare, window
    )
    ql_abs = _series_scores(
        model["time"], ql_model_compare, obs["time"], ql_obs_compare, window, lag=0.0
    )
    ql_lag = _series_scores(
        model["time"], ql_model_compare, obs["time"], ql_obs_compare, window, lag=lag
    )
    ql_total_abs = _series_scores(
        model["time"], ql_model_total, obs["time"], ql_obs_total, window, lag=0.0
    )
    ql_above_abs = _series_scores(
        model["time"], ql_model_above, obs["time"], ql_obs_above, window, lag=0.0
    )
    nd_score = _series_scores(model["time"], nd_model_2um, obs["time"], obs["ndrop_psd"], window)
    deff_score = _series_scores(model["time"], deff_model_um, obs["time"], obs["deff_um"], window)
    rel_score = _series_scores(model["time"], rel_model, obs["time"], obs["rel_disp"], window)

    # Integrated liquid-water exposure is much less sensitive to a short sample-line lag than a peak.
    mt_mask = (model["time"] >= window[0]) & (model["time"] <= window[1])
    ot_mask = (obs["time"] >= window[0]) & (obs["time"] <= window[1])
    int_ql_model = float(np.trapezoid(ql_model_compare[mt_mask], np.asarray(model["time"])[mt_mask])) if np.count_nonzero(mt_mask) > 1 else np.nan
    int_ql_obs = float(np.trapezoid(ql_obs_compare[ot_mask], obs["time"][ot_mask])) if np.count_nonzero(ot_mask) > 1 else np.nan
    int_ql_model_total = float(np.trapezoid(ql_model_total[mt_mask], np.asarray(model["time"])[mt_mask])) if np.count_nonzero(mt_mask) > 1 else np.nan
    int_ql_obs_total = float(np.trapezoid(ql_obs_total[ot_mask], obs["time"][ot_mask])) if np.count_nonzero(ot_mask) > 1 else np.nan
    int_ql_model_above = float(np.trapezoid(ql_model_above[mt_mask], np.asarray(model["time"])[mt_mask])) if np.count_nonzero(mt_mask) > 1 else np.nan
    int_ql_obs_above = float(np.trapezoid(ql_obs_above[ot_mask], obs["time"][ot_mask])) if np.count_nonzero(ot_mask) > 1 else np.nan

    fig, axes = plt.subplots(4, 2, figsize=(14, 16), sharex=False)
    ax = axes.ravel()
    # Pressure
    ax[0].plot(met["Time"] / 60.0, met["Pressure"], label="observed/forced")
    ax[0].plot(model["time"] / 60.0, np.asarray(model["p"]) / 100.0, "--", label="BMM")
    ax[0].set_ylabel("Pressure (hPa)"); ax[0].legend(); ax[0].grid()
    # Temperature
    ax[1].plot(met["Time"] / 60.0, met["Tgw mean"], label="gas obs/forcing")
    if "Tww mean" in met:
        ax[1].plot(met["Time"] / 60.0, met["Tww mean"], label="wall")
    ax[1].plot(model["time"] / 60.0, np.asarray(model["t"]) - 273.15, "--", label="BMM gas")
    ax[1].set_ylabel("Temperature (degC)"); ax[1].legend(); ax[1].grid()
    # Liquid water.  Always show both the unrestricted/reader-total quantities
    # and the directly matched WELAS-equivalent D>Dmin virtual-instrument pair.
    total_selected = ql_mode == "total"
    above_selected = ql_mode == "above_min"
    ax[2].plot(
        obs["time"] / 60.0, ql_obs_total * 1.0e3,
        linestyle="-" if total_selected else "--",
        alpha=1.0 if total_selected else 0.55, label="WELAS total LWC",
    )
    ax[2].plot(
        model["time"] / 60.0, ql_model_total * 1.0e3,
        linestyle="-" if total_selected else "--",
        alpha=1.0 if total_selected else 0.55, label="BMM total ql",
    )
    ax[2].plot(
        obs["time"] / 60.0, ql_obs_above * 1.0e3,
        linestyle="-" if above_selected else ":",
        alpha=1.0 if above_selected else 0.55,
        label=f"WELAS PSD >{cfg.SINGLE_COMPARE_DROP_MIN_UM:g} um",
    )
    ax[2].plot(
        model["time"] / 60.0, ql_model_above * 1.0e3,
        linestyle="-" if above_selected else ":",
        alpha=1.0 if above_selected else 0.55,
        label=f"BMM WELAS-eq >{cfg.SINGLE_COMPARE_DROP_MIN_UM:g} um",
    )
    if abs(lag) > 0.0:
        ax[2].plot(
            (obs["time"] - lag) / 60.0, ql_obs_compare * 1.0e3,
            ":", alpha=0.75, label=f"selected obs shifted {-lag:+.0f}s",
        )
    ax[2].set_ylabel(r"$q_l$ (g kg$^{-1}$)"); ax[2].legend(fontsize=8); ax[2].grid()
    # Number
    ax[3].plot(obs["time"] / 60.0, obs["ndrop"], label="OPC Nd field")
    ax[3].plot(obs["time"] / 60.0, obs["ndrop_psd"], ":", label=f"OPC PSD >{cfg.SINGLE_COMPARE_DROP_MIN_UM:g} um")
    ax[3].plot(model["time"] / 60.0, nd_model, label="BMM activated")
    ax[3].plot(model["time"] / 60.0, nd_model_2um, "--", label=f"BMM >{cfg.SINGLE_COMPARE_DROP_MIN_UM:g} um")
    ax[3].set_ylabel(r"$N_d$ (cm$^{-3}$)"); ax[3].legend(fontsize=8); ax[3].grid()
    # Effective diameter
    ax[4].plot(obs["time"] / 60.0, obs["deff_um"], label="OPC")
    ax[4].plot(model["time"] / 60.0, deff_model_um, label=f"BMM >{cfg.SINGLE_COMPARE_DROP_MIN_UM:g} um")
    ax[4].plot(model["time"] / 60.0, deff_model_activated_um, "--", alpha=0.6, label="BMM activated")
    ax[4].set_ylabel(r"$D_{eff}$ ($\mu$m)"); ax[4].legend(); ax[4].grid()
    # Relative dispersion
    ax[5].plot(obs["time"] / 60.0, obs["rel_disp"], label="OPC PSD")
    ax[5].plot(model["time"] / 60.0, rel_model, label=f"BMM >{cfg.SINGLE_COMPARE_DROP_MIN_UM:g} um")
    ax[5].plot(model["time"] / 60.0, rel_model_activated, "--", alpha=0.6, label="BMM activated")
    ax[5].set_ylabel(r"Relative dispersion $\sigma_D/\bar D$"); ax[5].legend(); ax[5].grid()
    # Cumulative chamber/fallout losses
    loss_fields = [
        ("qfan_liq", "fan"), ("qwall_liq", "wall"),
        ("qfall_liq", "fallout"), ("qchamber_bl", "BL wall water"),
    ]
    any_loss = False
    for name, label in loss_fields:
        if name in model:
            ax[6].plot(model["time"] / 60.0, np.asarray(model[name]) * 1.0e3, label=label)
            any_loss = True
    ax[6].set_ylabel(r"Cumulative liquid/water loss (g kg$^{-1}$)")
    if any_loss: ax[6].legend()
    ax[6].grid()
    # Summary panel
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
        f"Diagnostic ql lag: {lag:+.0f} s (corr={lag_corr:.3f})\n"
        f"selected ql NRMSE: {ql_abs['nrmse']:.3f}\n"
        f"selected ql NRMSE lagged: {ql_lag['nrmse']:.3f}\n"
        f"total ql NRMSE: {ql_total_abs['nrmse']:.3f}\n"
        f">{cfg.SINGLE_COMPARE_DROP_MIN_UM:g}um ql NRMSE: {ql_above_abs['nrmse']:.3f}\n"
        f"Nd(>{cfg.SINGLE_COMPARE_DROP_MIN_UM:g} um) NRMSE: {nd_score['nrmse']:.3f}\n"
        f"Deff NRMSE: {deff_score['nrmse']:.3f}\n"
        f"Rel. dispersion NRMSE: {rel_score['nrmse']:.3f}\n"
        f"Integral selected ql M/O: {int_ql_model:.4g} / {int_ql_obs:.4g} kg s kg-1\n"
        f"Integral total BMM ql: {int_ql_model_total:.4g} kg s kg-1"
    )
    ax[7].text(0.02, 0.98, summary, va="top", family="monospace")

    for a in ax[:7]:
        a.axvline(onset / 60.0, color="0.5", linestyle=":", linewidth=1)
        a.axvspan(window[0] / 60.0, window[1] / 60.0, color="0.9", alpha=0.25)
        a.set_xlabel("Time (min)")
    fig.suptitle(f"{exp}: BMM vs iSKYLAB observations")
    fig.tight_layout()

    # Direct WELAS/model PSD comparison on a deliberately coarser common
    # logarithmic grid.  The native WELAS grid is much finer than the number
    # of independent moving BMM populations, so putting model delta-populations
    # into every native optical channel produces a sparse/striped panel.
    compare_edges = _common_comparison_psd_edges(obs["Dp_edges_um"])
    obs_psd_compare = _rebin_dndlog10d(
        obs["psd"], obs["Dp_edges_um"], compare_edges
    )
    model_psd = _model_psd_on_fixed_grid(
        model, compare_edges, diameter_key="dwet", number_key="nwat"
    )
    psd_fig = None
    if model_psd is not None:
        psd_fig, pax = plt.subplots(2, 1, figsize=(13, 9), sharex=True, sharey=True)
        positive = np.concatenate([
            obs_psd_compare[np.isfinite(obs_psd_compare) & (obs_psd_compare > 0.0)],
            model_psd[np.isfinite(model_psd) & (model_psd > 0.0)],
        ])
        if positive.size:
            vmin = max(float(np.nanpercentile(positive, 2.0)), 1.0e-4)
            vmax = max(float(np.nanpercentile(positive, 99.5)), 10.0 * vmin)
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
            f"BMM rebinned to the same {len(compare_edges)-1} log-D bins"
        )
        for a in pax:
            a.set_yscale("log")
            a.set_ylabel(r"Wet diameter ($\mu$m)")
            a.axvline(onset / 60.0, color="w", linestyle=":", linewidth=1)
        pax[1].set_xlabel("Time (min)")
        psd_fig.colorbar(pcm1, ax=pax, label=r"dN/dlog$_{10}$D (cm$^{-3}$)")
        psd_fig.suptitle(
            f"{exp}: observed and model liquid/warm-particle size distributions "
            "(common coarse grid)"
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
        "best_ql_lag_s": lag, "best_ql_lag_corr": lag_corr,
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
        with (outdir / f"metrics-{exp}.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["metric", "value"])
            writer.writerows(metrics.items())
        print(f"Single-experiment diagnostics written to {outdir}")
    if show:
        plt.show()
    return fig, psd_fig, full_psd_fig, metrics

def main(group=THIS_RUN, run_model=RUN_MODEL, do_analysis=DO_ANALYSIS, do_plot=DO_PLOT,
         experiment=None, winit_single=1.3, saturation_time_min=None):
    if not READ_DATA:
        raise RuntimeError("READ_DATA=False is no longer supported without supplying a cached data dictionary")
    data = load_all_data()
    _write_forcing_smoothing_diagnostics(data)

    exp = _normalise_experiment_name(experiment) if experiment is not None else None
    if saturation_time_min is not None and exp is None:
        raise ValueError("--saturation-time-min is only valid with --experiment")
    if saturation_time_min is not None and cfg.INITIAL_RH_METHOD != "cloud_onset":
        raise ValueError(
            "--saturation-time-min requires INITIAL_RH_METHOD='cloud_onset'; "
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
            model_file = run_single_experiment(exp, states[exp], data, winit=winit_single)
        elif not model_file.exists():
            raise FileNotFoundError(f"Existing single-run output not found: {model_file}")
        if do_analysis:
            return analyse_single_experiment(
                exp, data, model_file, show=do_plot,
                saturation_time_s=states[exp]["saturation_time"],
            )
        return model_file

    if group < 0 or group >= len(meta.BATCH_GROUPS):
        raise ValueError(f"group must be between 0 and {len(meta.BATCH_GROUPS)-1}")
    batch_sims = meta.BATCH_GROUPS[group]
    winit = meta.GROUP_UPDRAFT[group]

    vals1 = vals2 = None
    if run_model:
        vals1, vals2 = run_batch(batch_sims, states, data, winit)

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
        return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    target = parser.add_mutually_exclusive_group()
    target.add_argument("--group", type=int, default=None, help="index in experiment_metadata.BATCH_GROUPS")
    target.add_argument("--experiment", help="single experiment, e.g. Exp005 or 5")
    parser.add_argument("--winit", type=float, default=1.3, help="initial/updraft value used for a single run")
    parser.add_argument(
        "--saturation-time-min",
        type=float,
        default=None,
        help=(
            "single-run model saturation target in minutes from experiment start; "
            "overrides the config/default CLOUD_ONSET target without moving the "
            "observed onset marker or comparison window"
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
        run_model=not args.no_run,
        do_analysis=not args.no_analysis,
        do_plot=not args.no_plot,
    )
