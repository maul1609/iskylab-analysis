"""Read commonly used BMM output fields, including optional chamber diagnostics."""

from __future__ import annotations

from netCDF4 import Dataset
import numpy as np

R_D = 8.314 / 28.96e-3

# Variables written only when the corresponding BMM process is enabled.
OPTIONAL_DIAGNOSTICS = [
    # Bulk/PSD diagnostics introduced for direct model-observation comparison.
    "dmean_liq",
    "dvol_liq",
    "rel_disp_liq",
    "dwet",
    "nliq",
    "nwat",
    "mwat",
    "dmaxice",
    "dmean_ice",
    "rel_disp_ice",
    "nicem",
    "mice",
    "qi",
    "nice",
    "qchamber_bl",
    "qchamber_bl_step",
    "qchamber_bl_evap",
    "qchamber_bl_evap_step",
    "qfan_liq",
    "nfan_liq",
    "qfan_ice",
    "nfan_ice",
    "qwall_liq",
    "nwall_liq",
    "qwall_ice",
    "nwall_ice",
    "qfall_liq",
    "nfall_liq",
    "fallrate_liq",
    "qfall_ice",
    "nfall_ice",
    "fallrate_ice",
]


def readData(fileName, modelStr="Model-Exp005"):
    """Return a dictionary of core fields and any available loss diagnostics."""
    with Dataset(fileName) as nc:
        time = np.asarray(nc["time"][:])
        p = np.asarray(nc["p"][:])
        t = np.asarray(nc["t"][:])
        out = {
            "time": time,
            "p": p,
            "t": t,
            "rh": np.asarray(nc["rh"][:]),
            "ndrop": np.asarray(nc["ndrop"][:]),
            "ql": np.asarray(nc["ql"][:]),
            "deff": np.asarray(nc["deff"][:]),
            "rhoa": p / (R_D * t),
        }
        for name in OPTIONAL_DIAGNOSTICS:
            if name in nc.variables:
                out[name] = np.asarray(nc[name][:])
    return {modelStr: out}


if __name__ == "__main__":
    raise SystemExit("Usage: import readModel and call readData(fileName, modelStr)")
