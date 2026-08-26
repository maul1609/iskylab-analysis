"""Diagnose homogeneous/inhomogeneous mixing signatures from iSKYLAB OPC data.

The plot compares the cube of the observed mean-drop diameter relative to an
adiabatic estimate against the retained drop-number fraction.  The original
script contained two unit/algebra errors: air density had misplaced
parentheses and OPC LWC was multiplied by air density a second time when
forming mass per drop.  Both are corrected here.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

import experiment_metadata as meta
import readMeteoCPC
import readOPC_Merged
import svp

R_D = 8.314 / 28.96e-3
R_V = 8.314 / 18.0e-3
EPSILON = R_D / R_V

# Pair used in the historical mixing comparison.
BATCH_SIMS = ["Exp005", "Exp006"]


def _svp_liq(t_k):
    return np.asarray(svp.svp(np.asarray(t_k), "buck2", "liq"), dtype=float)


def analyse_experiment(exp):
    met_key = f"MeteoCPC-{exp}"
    opc_key = f"MergedOPC-{exp}"
    i_met = readMeteoCPC.metStr.index(met_key)
    i_opc = readOPC_Merged.opcStr.index(opc_key)
    met = readMeteoCPC.readData(i_met, met_key)[met_key]
    opc = readOPC_Merged.readData(i_opc, opc_key)[opc_key]

    cloud_time = meta.CLOUD_ONSET[exp]
    cloud_idx = np.flatnonzero(met["Time"] >= cloud_time)
    if cloud_idx.size == 0:
        raise ValueError(f"Cloud-onset time outside meteo series for {exp}")
    i0 = cloud_idx[0]

    p_pa = met["Pressure"] * 100.0
    t_k = met["Tgw mean"] + 273.15
    rho_d = p_pa / (R_D * t_k)

    # Adiabatic liquid-water mixing ratio: saturation mixing ratio at cloud
    # onset minus the saturation mixing ratio along the subsequent measured
    # P/T trajectory, clipped at zero before cloud onset.
    es0 = _svp_liq([t_k[i0]])[0]
    qsat0 = EPSILON * es0 / (p_pa[i0] - es0)
    es = _svp_liq(t_k)
    qsat = EPSILON * es / (p_pa - es)
    ql_ad = np.maximum(qsat0 - qsat, 0.0)

    max_ndrop = np.nanmax(opc["ndrop"])
    nrat = opc["ndrop"] / max_ndrop if max_ndrop > 0.0 else np.full_like(opc["ndrop"], np.nan)

    # Diameter expected if the adiabatic water content were shared equally
    # between the maximum observed number of droplets.
    mass_per_drop_ad = ql_ad * rho_d / (max_ndrop * 1.0e6)
    dmean_ad = (6.0 * mass_per_drop_ad / (np.pi * 1000.0)) ** (1.0 / 3.0) * 1.0e6

    # OPC-derived equal-mass diameter.  LWC is already g m-3, so after the
    # 1e-3 conversion to kg m-3 it must NOT be multiplied by rho_d again.
    with np.errstate(divide="ignore", invalid="ignore"):
        mass_per_drop_obs = opc["lwc"] * 1.0e-3 / (opc["ndrop"] * 1.0e6)
        dmean_lwc = (6.0 * mass_per_drop_obs / (np.pi * 1000.0)) ** (1.0 / 3.0) * 1.0e6

    # Number-weighted mean diameter for measured particles >2 um.
    bins = opc["Dp"] > 2.0
    denom = np.sum(opc["Conc"][:, bins], axis=1)
    dmean = np.full_like(denom, np.nan, dtype=float)
    good = denom > 0.0
    dmean[good] = (
        np.sum(opc["Conc"][:, bins] * opc["Dp"][bins], axis=1)[good] / denom[good]
    )

    # Meteo/OPC products are intended to share the same expansion time grid;
    # interpolate dmean_ad if their lengths/times differ rather than relying on
    # array-position coincidence.
    dmean_ad_opc = np.interp(opc["Time"], met["Time"], dmean_ad)
    with np.errstate(divide="ignore", invalid="ignore"):
        d3_ratio = (dmean / dmean_ad_opc) ** 3

    return opc["Time"], nrat, d3_ratio, dmean_lwc


def main():
    fig, axes = plt.subplots(1, len(BATCH_SIMS), figsize=(6 * len(BATCH_SIMS), 5), squeeze=False)
    for ax, exp in zip(axes[0], BATCH_SIMS):
        time, nrat, d3_ratio, _ = analyse_experiment(exp)
        sc = ax.scatter(nrat, d3_ratio, s=12, c=time)
        ax.set_xlim((0, 1))
        ax.set_ylim((0, 1))
        ax.set_xlabel(r"$N/N_{adia}$")
        ax.set_ylabel(r"$D^3/D_{adia}^3$")
        ax.set_title(exp)
        ax.grid()
        fig.colorbar(sc, ax=ax, label="Time (s)")
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
