"""User-editable configuration for iSKYLAB -> BMM batch simulations.

The analysis/data readers intentionally keep experiment metadata separate from
model-control parameters.  This file contains the latter so that changing a
chamber sensitivity does not require editing the namelist generator itself.

Paths are resolved relative to this repository wherever possible.  Set the
BMM_MODEL_FOLDER environment variable if the BMM source tree lives elsewhere.
"""

from __future__ import annotations

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = REPO_ROOT / "iSKYLAB-data"

# BMM source/build directory.  The historical path is retained as a fallback,
# but an environment variable is preferable when moving between machines.
BMM_MODEL_FOLDER = Path(
    os.environ.get(
        "BMM_MODEL_FOLDER",
        "/Users/mccikpc2/Dropbox/programming/fortran/bmm",
    )
).expanduser()

# A current BMM iSKYLAB namelist is used only as a template.  The batch script
# edits variables by name and replaces the whole &chamber_spec block, so it no
# longer depends on the template containing a particular old value or array
# length.
BMM_TEMPLATE = (
    BMM_MODEL_FOLDER
    / "python"
    / "aida_analysis"
    / "iSKYLAB-namelists"
    / "namelist-aida-exp005.in"
)
BMM_EXECUTABLE = BMM_MODEL_FOLDER / "main.exe"

# Model outputs and generated namelists are kept in the user's /tmp directory.
OUTPUT_ROOT = Path("/tmp") / os.environ.get("USER", "iskylab")
SAVE_GENERATED_NAMELISTS = True

# ---------------------------------------------------------------------------
# Observational forcing
# ---------------------------------------------------------------------------
FORCE_PRESSURE = True
FORCE_TEMPERATURE = True

# Chamber pressure/gas/wall-temperature measurements contain short-period
# instrument noise that can create unrealistically large dT/dt and repeatedly
# move particles across saturation/activation thresholds.  Smooth the forcing
# before it is written to the BMM, while retaining the raw observations for
# diagnostics.  Windows are specified in seconds (not samples).
# Recommended starting configuration for chamber-forcing smoothing.
# Add these to iskylab_config.py.
SMOOTH_CHAMBER_FORCING = True

# "raw", "savgol", or "butterworth"
CHAMBER_SMOOTH_METHOD = "butterworth"

# Robust isolated-spike removal before either smoother.
CHAMBER_DESPIKE = True
CHAMBER_DESPIKE_WINDOW = 10.0       # seconds
CHAMBER_DESPIKE_NSIGMA = 2.0

# Physical smoothing timescales.
# Savitzky-Golay: fitting-window duration.
# Butterworth: -3 dB cutoff period.
CHAMBER_SMOOTH_TEMP_WINDOW = 30.0       # gas temperature, seconds
CHAMBER_SMOOTH_WALL_TEMP_WINDOW = 45.0  # wall-gas delta T, seconds
CHAMBER_SMOOTH_PRESSURE_WINDOW = 15.0   # pressure, seconds

CHAMBER_SMOOTH_POLYORDER = 2
CHAMBER_BUTTERWORTH_ORDER = 2

# Recommended: preserve the physically important wall-minus-gas contrast.
CHAMBER_SMOOTH_WALL_AS_DELTA_T = True

# Write one diagnostic PNG per experiment plus a CSV summary of the smoothing
# corrections under /tmp/$USER/forcing_smoothing.
SAVE_FORCING_SMOOTHING_DIAGNOSTICS = True

# qtot_chamber can be written for diagnostic/sensitivity use without forcing
# the model.  The forcing is deliberately OFF by default because the water
# measurement/interpretation is less certain and because forcing qtot while
# also applying wall/fan sinks can double-count chamber water loss.
WRITE_QTOT_DATA = True
FORCE_QTOT = False

# Source used to construct qtot_chamber when WRITE_QTOT_DATA is True:
#   "vapour_plus_opc" : qv from dew point + OPC liquid-water mixing ratio.
#   "dewpoint_only"   : legacy proxy; only appropriate if the dew-point
#                       instrument is known to represent total sampled water.
QTOT_DATA_MODE = "vapour_plus_opc"

# Initial RH method:
#   "cloud_onset" : infer qv so the model reaches liquid saturation at a
#                   configurable target time along the measured P/T forcing.
#                   With the controls below left at their defaults this is the
#                   historical CLOUD_ONSET time from experiment_metadata.py.
#   "dewpoint"    : use the measured dew point at t=0 directly.
INITIAL_RH_METHOD = "cloud_onset"

# Model saturation-target controls used only when INITIAL_RH_METHOD is
# "cloud_onset".  These alter the initial vapour mixing ratio; they do NOT move
# the observed cloud-onset marker/window used in the comparison plots.
#
# The global shift is useful for systematic sensitivity tests.  Per-experiment
# absolute overrides take precedence.  Values are seconds from the experiment
# time origin, e.g. {"Exp005": 6.6*60.0}.
MODEL_SATURATION_TIME_SHIFT_S = 0.0
MODEL_SATURATION_TIME_OVERRIDES_S = {'Exp005': 8. * 60.0, \
	'Exp006': 12. * 60.0,'Exp026': 10 * 60.0, \
	'Exp007': 12.5 * 60.0,'Exp008': 13. * 60.0,'Exp009': 12. * 60.0, \
	'Exp010': 11. * 60.0,'Exp011': 12.5 * 60.0 , \
	'Exp014': 11. * 60.0,'Exp019': 13. * 60.0}

# Aerosol number normalisation:
#   "t0"          : CPC concentration at t=0.  Recommended when modelling
#                   pre-cloud particle loss (fan, walls or fallout).
#   "cloud_onset" : historical choice; constrains aerosol immediately before
#                   activation and therefore hides any pre-cloud loss.
AEROSOL_INIT_TIME = "cloud_onset"

# ---------------------------------------------------------------------------
# Chamber boundary-layer treatment
# ---------------------------------------------------------------------------
# Coupled chamber BL operator:
#   0 = off
#   1 = on
#
# When enabled, the BL thermodynamics first diagnoses the effective sensible
# temperature
#
#   T_sens = T_gas + CHAMBER_BL_ALPHA_T*(T_wall-T_gas)
#                    + CHAMBER_BL_TEMP_OFFSET
#
# before solving any latent phase change.  The measured wall temperature is
# therefore a boundary condition rather than being assumed to equal the BL-air
# temperature.  Useful limits are:
#   ALPHA_T=0, OFFSET=0 : T_sens follows the measured gas temperature
#   ALPHA_T=1, OFFSET=0 : T_sens follows the measured wall temperature
#   ALPHA_T=0, OFFSET=-0.12 : historical fixed -0.12 K BL sensitivity
CHAMBER_BL_MIX = 1
CHAMBER_BL_TAU = 10.  # s, chamber-scale BL processing/recirculation timescale
CHAMBER_BL_ALPHA_T = 0.0
CHAMBER_BL_TEMP_OFFSET = -0.075  # K, optional unresolved BL sensible-temperature offset
#CHAMBER_BL_TEMP_OFFSET = 0.09  # K, optional unresolved BL sensible-temperature offset

# CHAMBER_BL_TAU = 50.  # s, chamber-scale BL processing/recirculation timescale
# CHAMBER_BL_ALPHA_T = 0.0
#CHAMBER_BL_TEMP_OFFSET=5.0


"""
	No fan: TAU ~ 80; ALPHA_T=1.0 (governed by the wall in quiescent conds); DT=0.0
	Fan: TAU ~ 10; ALPHA_T=0.0 (gas / well mixed), but DT=-0.1
	
	Exp 17: wall is colder!!
CHAMBER_BL_TAU = 10.  # s, chamber-scale BL processing/recirculation timescale
CHAMBER_BL_ALPHA_T = 0.0
CHAMBER_BL_TEMP_OFFSET = -0.6  # K, optional unresolved BL sensible-temperature 

Note for Exp 21 I had to use +0.1 K... maybe something to do with colder / humidy?
wall was warmer in this case
"""

# How thermodynamically required LIQUID evaporation is represented in the PSD:
#   1 = homogeneous diffusional evaporation.  Every activated droplet receives
#       the same finite D^2 decrement; small droplets may naturally evaporate
#       completely while larger droplets shrink.
#   2 = uniform extreme inhomogeneous.  The same fraction of every activated
#       bin evaporates completely; survivors retain their original size.
#   3 = D2-lifetime-weighted extreme inhomogeneous.  Complete-evaporation
#       fractions are biased toward smaller droplets as m_w^(-2/3), capped by
#       the BL-processed fraction; survivors retain their original size.
CHAMBER_BL_EVAP_MODE = 3

# ---------------------------------------------------------------------------
# Drone-fan blade collection
# ---------------------------------------------------------------------------
# 0=off, 1=saturating sigmoid in current particle diameter.
CHAMBER_FAN_LOSS = 0
CHAMBER_FAN_LOSS_KMAX = 7.0e-3 #1.5e-3      # s-1
CHAMBER_FAN_LOSS_D50_REF = 6e-6 #10.0e-6   # m at reference RPM
CHAMBER_FAN_LOSS_EXP = 6.0
CHAMBER_FAN_RPM = 25000.0
CHAMBER_FAN_RPM_REF = 25000.0

# ---------------------------------------------------------------------------
# Non-gravitational wall deposition
# ---------------------------------------------------------------------------
# 0=off, 1=Lai-Nazaroff deposition to cylindrical side walls + ceiling.
# The floor is not included here because gravitational floor loss is handled
# independently by the generic BMM fallout scheme.
CHAMBER_WALL_LOSS = 1
CHAMBER_WALL_USTAR = 0.02 #0.02  # m s-1
CHAMBER_DIAMETER = 1.5    # m
CHAMBER_HEIGHT = 2.15       # m

# ---------------------------------------------------------------------------
# Generic BMM fallout/sedimentation
# ---------------------------------------------------------------------------
# In a chamber run the updated BMM uses CHAMBER_HEIGHT as V/A_floor, so the
# value of residence_depth in the template is not the controlling chamber
# length scale.  It is still set consistently by the generator for clarity.
FALLOUT_FLAG = True
RESIDENCE_DEPTH = CHAMBER_HEIGHT

# Synthetic updraft mode reproduces the old speedFlag behaviour.  It disables
# measured chamber thermodynamic forcing and chamber-specific BL/fan/wall
# processes, while leaving generic BMM options (e.g. fallout) configurable.
SYNTHETIC_UPDRAFT = False

# ---------------------------------------------------------------------------
# Single-experiment model/observation diagnostics
# ---------------------------------------------------------------------------
# Minimum OPC/model wet diameter used for cloud-drop bulk moments.
SINGLE_COMPARE_DROP_MIN_UM = 2.5

# Liquid-water comparison used for the single-experiment ql score/lag:
#   "above_min" : recommended instrument-equivalent comparison.  Reconstruct
#                 liquid-equivalent water from the observed WELAS PSD and the
#                 BMM nwat+dwet population using all particles with
#                 Dwet > SINGLE_COMPARE_DROP_MIN_UM.  WELAS is mimicked by
#                 interpreting the entire measured/model wet sphere as water.
#   "total"     : compare the WELAS reader's total LWC with the BMM's native
#                 total liquid-water mixing ratio (model["ql"]).  This is useful
#                 diagnostically, but the sampled size ranges are not identical.
#
# Both total and >Dmin curves/scores are always calculated and written; this
# option only chooses which pair drives the headline ql lag/NRMSE metrics.
SINGLE_COMPARE_QL_MODE = "above_min"

# Direct WELAS/BMM PSD comparison grid.  WELAS has far more size channels than
# the BMM has independent moving particle populations; putting each BMM point
# onto the native WELAS grid therefore produces a visually sparse/striped model
# panel.  Both observation and model are instead conservatively rebinned onto
# this common logarithmic diagnostic grid for the comparison figure.
#
# Set MIN/MAX to None to use the native WELAS edge range for each experiment.
SINGLE_COMPARE_PSD_NBINS = 48
SINGLE_COMPARE_PSD_MIN_UM = None
SINGLE_COMPARE_PSD_MAX_UM = None

# Diagnostic-only time-lag search.  Absolute-time comparisons are always
# retained; this bounded lag is reported separately to identify likely
# instrument/sample-line delay rather than silently shifting the observations.
SINGLE_COMPARE_MAX_LAG_S = 30.0
SINGLE_COMPARE_LAG_STEP_S = 1.0

# Save the main model/observation panel, the observed/model PSD comparison and
# a compact CSV of scalar scores under OUTPUT_ROOT/single_comparison.
SAVE_SINGLE_COMPARISON = True

# Full model-only PSD diagnostic.  Unlike the direct OPC comparison this uses
# every warm BMM particle (nwat+dwet), including unactivated aerosol and haze.
# Ice uses nicem+dmaxice.  Fixed logarithmic grids are used only for plotting;
# they do not alter the model or imply native BMM wet-bin widths.
SINGLE_MODEL_PSD_NBINS = 180
SINGLE_MODEL_WARM_PSD_MIN_UM = 0.005
SINGLE_MODEL_WARM_PSD_MAX_UM = 1000.0
SINGLE_MODEL_ICE_PSD_MIN_UM = 0.1
SINGLE_MODEL_ICE_PSD_MAX_UM = 10000.0


# ---------------------------------------------------------------------------
# Dust aerosol / ice-active-site settings for AIDAd Exp020-023 and Exp025-027
# ---------------------------------------------------------------------------
# FHH adsorption coefficients are written into the BMM namelist by component.
# SDSA01 uses the Kumar et al. fresh natural-mineral-dust fit as a proxy; ATD03
# uses the dry-generated Arizona Test Dust fit from Kumar et al. (2009).
DUST_FHH = {
    "SDSA01": {"A": 2.25, "B": 1.20},
    "ATD03":  {"A": 0.27, "B": 0.79},
}

# BMM component category strings.  The actual ns(T) fits live in
# bin_microphysics_module.f90, not in this Python repository.
DUST_INP_CATEGORY = {
    "SDSA01": "sdsa01",
    "ATD03": "atd03",
}

# Common cumulative IASD thresholds [deg C], corresponding to 255 ... 240 K
# at exactly 1 K spacing.  There are therefore 16 thresholds.
INP_TEMP_C = [
    -18.15, -19.15, -20.15, -21.15,
    -22.15, -23.15, -24.15, -25.15,
    -26.15, -27.15, -28.15, -29.15,
    -30.15, -31.15, -32.15, -33.15,
]

# Dust experiments are configured as mixed-phase/ice runs and use Koop + INAS.
# DeMott is disabled for these components so the explicit dust IASD reservoir
# is not double counted.
DUST_ENABLE_ICE = True
DUST_ICE_NUCLEATION_MECH = [True, True, False, False]  # Koop, INAS, DeMott, Daily
