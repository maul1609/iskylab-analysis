# iSKYLAB analysis and BMM simulation driver

This repository contains the Python workflow used to initialise, run and analyse
iSKYLAB/AIDAd chamber experiments with the **Bin Microphysics Model (BMM)**.

The maintained entry point is:

```text
scripts/dataAnalysis_new.py
```

The driver can:

- read chamber meteorology, CPC, WELAS/OPC and initial aerosol PSD data;
- smooth measured chamber pressure and temperature forcing;
- construct aerosol initial conditions from fitted measured PSDs;
- generate experiment-specific BMM namelists by variable name;
- run a compiled BMM executable;
- run one experiment or a configured experiment group;
- override important chamber sensitivities from the command line;
- compare BMM output with WELAS/OPC observations;
- construct instrument-like number, liquid-water and size diagnostics;
- compare warm and mixed-phase size distributions on a common grid;
- diagnose wall, boundary-layer, fan and sedimentation/fallout processes;
- configure dust FHH activation and explicit ice-active-site-density (INAS/IASD)
  temperature classes;
- generate variable-length `inp_temp` arrays automatically.

The repository does **not** contain the BMM source tree or the iSKYLAB data.

---

## 1. Repository layout

The main files are:

```text
scripts/
    dataAnalysis_new.py        main BMM driver and analysis
    iskylab_config.py          user-facing model/analysis configuration
    experiment_metadata.py     experiment groups and timing metadata
    batchRuns.sh               example experiment-by-experiment command lines
    namelist_utils.py          robust BMM namelist editing helpers

    readMeteoCPC.py            chamber meteorology/CPC reader
    readOPC_Merged.py          WELAS/OPC reader
    readPNSD_Mrg_new.py        initial aerosol PSD reader/fitter
    readModel.py               BMM NetCDF reader

    svp.py                     saturation-vapour-pressure utilities
    dataAnalysis.py            legacy wrapper
    adiabaticAnalysis.py       additional/legacy diagnostics
    mixingAnalysis.py          additional/legacy diagnostics
    overlayModelPlot.py        additional/legacy plotting

requirements.txt
README.md
```

The configured data location is:

```text
iSKYLAB-data/
```

at the repository root unless changed in `scripts/iskylab_config.py`.

The BMM namelist template is taken from the external BMM tree. The current
configuration expects a template such as:

```text
$BMM_MODEL_FOLDER/python/aida_analysis/iSKYLAB-namelists/namelist-aida-exp005.in
```

The template is only a schema/base configuration. Values are replaced by
variable name, and the complete `&chamber_spec` block is regenerated for each
experiment.

---

## 2. Requirements

Install the Python dependencies with:

```bash
python3 -m pip install -r requirements.txt
```

The current workflow uses NumPy, SciPy, Matplotlib and NetCDF support.

The BMM itself is a separate repository:

```text
https://github.com/UoM-maul1609/bin-microphysics-model
```

Compile BMM using its own build system first. Then point this repository at the
compiled BMM tree, preferably with:

```bash
export BMM_MODEL_FOLDER=/path/to/bin-microphysics-model
```

or edit `BMM_MODEL_FOLDER` in `scripts/iskylab_config.py`.

The driver expects:

```text
$BMM_MODEL_FOLDER/main.exe
```

to exist.

---

## 3. Quick start

Run commands from `scripts/`.

A single experiment:

```bash
cd scripts
python3 dataAnalysis_new.py --experiment Exp006
```

The short numeric form is also accepted:

```bash
python3 dataAnalysis_new.py --experiment 6
```

A mixed-phase/dust experiment:

```bash
python3 dataAnalysis_new.py --experiment Exp021
```

A configured experiment group:

```bash
python3 dataAnalysis_new.py --group 6
```

Single-experiment mode reads only the requested experiment. Missing data from
unrelated experiments therefore do not prevent a single experiment from
running.

---

## 4. Command-line interface

Use:

```bash
python3 dataAnalysis_new.py --help
```

for the current list.

### Run selection

```text
--experiment ExpNNN
--group N
```

These are mutually exclusive.

### Model/run controls

```text
--winit VALUE
--no-run
--no-analysis
--no-plot
```

`--no-run` analyses an existing single-experiment output.

`--no-analysis` generates/runs the model but skips the observational comparison.

`--no-plot` suppresses display of plots; configured saved diagnostics are still
written.

### Cloud/saturation timing

```text
--cloud-formation-time-min X
--saturation-time-min X
```

These are aliases. They are valid for a single experiment when:

```python
INITIAL_RH_METHOD = "cloud_onset"
```

The value changes the **model saturation target used to infer initial vapour**.
It does not move the raw WELAS observations or the curated observational cloud
window.

Example:

```bash
python3 dataAnalysis_new.py \
    --experiment Exp006 \
    --cloud-formation-time-min 12
```

### Chamber boundary-layer overrides

```text
--bl-tau-s X
--bl-temp-offset-k X
--bl-evap-size-exp X
--bl-inhom-size-exp X       deprecated alias of --bl-evap-size-exp
--bl-wall-water-mode {0,1,2}
```

Example:

```bash
python3 dataAnalysis_new.py \
    --experiment Exp006 \
    --cloud-formation-time-min 12 \
    --bl-tau-s 10 \
    --bl-evap-size-exp 2 \
    --bl-wall-water-mode 0 \
    --bl-temp-offset-k -0.09
```

### Wall-water overrides

```text
--wall-vapour-transfer-velocity-ms X
--wall-liquid-water-init-kg X
--wall-ice-water-init-kg X
--wall-water-efficiency X
```

The transfer velocity is used by wall-water mode 2.

### SCE sensitivity

```text
--sce-bins N
--scebins N
```

This creates a **temporary copy** of the BMM SCE namelist, changes only
`n_binsc`, and points the generated main BMM namelist at that temporary file.

The repository/BMM copy of `sce/namelist.in` is not modified.

Example:

```bash
python3 dataAnalysis_new.py \
    --experiment Exp006 \
    --sce-bins 80
```

### Observation-time diagnostic

```text
--shift-cloud-measurements
```

The analysis first diagnoses the liquid-water lag between BMM and WELAS. With
this option, the complete WELAS time coordinate is shifted consistently,
including:

- liquid water;
- number;
- effective diameter;
- relative dispersion;
- ice number;
- the observed PSD.

Without this option the lag remains diagnostic only.

### Matplotlib interactive mode

Interactive plotting is off by default.

Use:

```text
--plt-ion
```

to call `plt.ion()` and keep Matplotlib in interactive mode.

Without `--plt-ion`, the driver uses non-interactive plotting and closes figures
after display, which is useful when calling multiple experiments from a shell
script.

For fully unattended runs, `--no-plot` is still the simplest option.

---

## 5. Example batch commands

`scripts/batchRuns.sh` contains working experiment-by-experiment command lines
for sensitivity studies.

A typical warm-cloud command is:

```bash
python3 dataAnalysis_new.py \
    --experiment Exp006 \
    --cloud-formation-time-min 12 \
    --bl-tau-s 10 \
    --bl-evap-size-exp 2 \
    --bl-wall-water-mode 0 \
    --bl-temp-offset-k -0.09
```

A mixed-phase example is:

```bash
python3 dataAnalysis_new.py \
    --experiment Exp021 \
    --cloud-formation-time-min 12.5 \
    --bl-tau-s 10 \
    --bl-evap-size-exp 2 \
    --bl-wall-water-mode 0 \
    --bl-temp-offset-k 0.4
```

The values in `batchRuns.sh` are experiment sensitivities/examples, not a
replacement for documenting the exact configuration used in a scientific
analysis. For reproducible work, retain the generated namelist.

---

## 6. Output locations

The default output root is:

```text
/tmp/$USER/
```

Single-experiment files are written under:

```text
/tmp/$USER/single_comparison/
```

Typical files are:

```text
output-Exp006.nc
namelist-Exp006.in
comparison-Exp006.png
psd-Exp006.png
model-full-psd-Exp006.png
wall-vapour-Exp006.png
metrics-Exp006.csv
```

The wall-vapour figure is only produced when the corresponding BMM diagnostics
are available.

Forcing-smoothing diagnostics are written under:

```text
/tmp/$USER/forcing_smoothing/
```

including experiment plots and:

```text
forcing-smoothing-method-comparison.csv
```

Generated namelists are retained when:

```python
SAVE_GENERATED_NAMELISTS = True
```

---

## 7. Main configuration

Normal sensitivity studies should be configured in:

```text
scripts/iskylab_config.py
```

rather than by editing the driver.

The main configuration groups are:

- BMM paths and output paths;
- measured chamber forcing;
- forcing smoothing/despiking;
- qtot diagnostics/forcing;
- initial humidity;
- model saturation timing;
- aerosol normalisation;
- chamber BL processing;
- wall-vapour exchange;
- fan-associated particle loss;
- non-gravitational wall deposition;
- sedimentation/fallout;
- chamber geometry;
- WELAS/BMM comparison definitions;
- PSD diagnostic grids;
- dust FHH and INAS/IASD settings.

---

## 8. Chamber meteorological forcing

Pressure, gas temperature and total water are independent model controls:

```python
FORCE_PRESSURE = True
FORCE_TEMPERATURE = True
FORCE_QTOT = False
```

The complete measured time series is written to the generated `&chamber_spec`
block using its actual length.

The driver writes:

```fortran
time_chamber(...)
press_chamber(...)
temp_chamber(...)
```

and, when enabled:

```fortran
qtot_chamber(...)
wall_temp_chamber(...)
```

There is no hard-coded historical chamber-data length.

### Total-water diagnostic

Dew point is converted to **vapour mixing ratio**, not treated as total water.

The recommended diagnostic mode is:

```python
QTOT_DATA_MODE = "vapour_plus_opc"
```

which constructs:

```text
qtot = qv(dew point) + ql(OPC)
```

after interpolation onto the meteorological time grid.

The alternative is:

```python
QTOT_DATA_MODE = "dewpoint_only"
```

`WRITE_QTOT_DATA=True` can write the diagnostic without forcing the model.

Keeping `FORCE_QTOT=False` is useful when explicit chamber water sinks/sources
are active, because imposing observed qtot at the same time can double-count
water exchange.

---

## 9. Chamber forcing smoothing

Measured chamber pressure and temperature contain short-period instrument noise.
The driver can despike and smooth these series before exporting them to BMM.

Current controls include:

```python
SMOOTH_CHAMBER_FORCING = True
CHAMBER_SMOOTH_METHOD = "butterworth"

CHAMBER_DESPIKE = True
CHAMBER_DESPIKE_WINDOW = 10.0
CHAMBER_DESPIKE_NSIGMA = 2.0

CHAMBER_SMOOTH_TEMP_WINDOW = 30.0
CHAMBER_SMOOTH_WALL_TEMP_WINDOW = 45.0
CHAMBER_SMOOTH_PRESSURE_WINDOW = 15.0

CHAMBER_SMOOTH_POLYORDER = 2
CHAMBER_BUTTERWORTH_ORDER = 2
```

Supported smoothing methods are:

```text
raw
savgol
butterworth
```

### Despiking

The pre-filter despiker uses a Hampel-style local median/MAD criterion. It is
designed to replace isolated outliers rather than smooth genuine chamber
structure.

### Savitzky-Golay

The Savitzky-Golay option uses a physical-time window and zero-lag filtering.
For sufficiently irregular sampling the data are temporarily interpolated onto
a median-dt grid.

### Butterworth

The Butterworth option uses a zero-phase `sosfiltfilt` low-pass filter. The
configured timescale is interpreted as the cutoff period.

### Preserving wall-gas temperature contrast

By default:

```python
CHAMBER_SMOOTH_WALL_AS_DELTA_T = True
```

and wall temperature is reconstructed as:

```text
Twall,filtered = Tgas,filtered + filtered(Twall - Tgas)
```

This avoids independently filtering gas and wall temperature in a way that
artificially distorts the wall-minus-gas temperature difference.

Raw, filled, despiked, Savitzky-Golay and Butterworth candidates are retained
for diagnostics.

To generate forcing diagnostics without running BMM:

```bash
python3 dataAnalysis_new.py --no-run --no-analysis --no-plot
```

provided the requested output/data path requirements are satisfied.

---

## 10. Initial humidity and cloud timing

Two initial-humidity methods are available:

```python
INITIAL_RH_METHOD = "cloud_onset"
```

or:

```python
INITIAL_RH_METHOD = "dewpoint"
```

### `cloud_onset`

The initial vapour mixing ratio is chosen so that, in the absence of pre-cloud
water exchange, the measured pressure/temperature trajectory would reach
liquid saturation at the specified model target time.

Persistent controls are:

```python
MODEL_SATURATION_TIME_SHIFT_S = 0.0
MODEL_SATURATION_TIME_OVERRIDES_S = {...}
```

A single run can override the target with:

```bash
--cloud-formation-time-min X
```

The observational onset marker and cloud-comparison window remain separate.

### `dewpoint`

The initial vapour mixing ratio is taken directly from the measured dew point
at `t=0`.

---

## 11. Aerosol initialisation

The initial aerosol PSD is derived from the measured/fitted PNSD and scaled to
the CPC concentration.

Aerosol number can be normalised at:

```python
AEROSOL_INIT_TIME = "t0"
```

or:

```python
AEROSOL_INIT_TIME = "cloud_onset"
```

If fan, wall or fallout losses are represented before cloud formation, `t0`
avoids normalising to a later concentration and then applying the pre-cloud
loss a second time.

The generated namelist includes fitted modal:

```text
number
dry diameter
log-width
component density
kappa
FHH coefficients
INP category
```

as applicable.

Aerosol diameters from the fitting workflow are converted from micrometres to
metres exactly once when written to BMM.

---

## 12. Chamber boundary-layer processing

The current interface separates:

1. whether chamber BL processing is active;
2. the thermodynamic state used for the processed air;
3. the wall-water closure;
4. how required liquid evaporation is represented in the PSD.

### Enable/disable

```python
CHAMBER_BL_MIX = 1
```

with:

```text
0 = off
1 = on
```

`CHAMBER_BL_MIX` is **not** the homogeneous/inhomogeneous evaporation selector.

### BL processing timescale

```python
CHAMBER_BL_TAU = 10.0
```

The BMM uses this as the chamber-scale BL processing/recirculation timescale.

For a timestep `dt`, the corresponding processed-air fraction is based on the
usual exponential relaxation form:

```text
fmix = 1 - exp(-dt/tau)
```

### Effective BL sensible temperature

The processed sensible temperature is based on:

```text
T_sens =
    T_gas
  + CHAMBER_BL_ALPHA_T * (T_wall - T_gas)
  + CHAMBER_BL_TEMP_OFFSET
```

Useful limits are:

```text
alpha = 0, offset = 0   -> follows bulk gas temperature
alpha = 1, offset = 0   -> follows measured wall temperature
```

A small offset can be used as an unresolved chamber/BL thermal sensitivity.

The command-line override is:

```bash
--bl-temp-offset-k X
```

---

## 13. Homogeneous and inhomogeneous evaporation

The way the thermodynamically required **liquid evaporation** is represented in
the PSD is selected with:

```python
CHAMBER_BL_EVAP_MODE = 2
```

Current modes are:

```text
1 = homogeneous evaporation
2 = inhomogeneous complete-particle evaporation
```

The thermodynamic liquid-water target is common to both modes. The difference
is how that target is distributed across the particle population.

### Mode 1: homogeneous

Particle number is retained while liquid mass/wet size decreases.

### Mode 2: inhomogeneous

Selected wet particles completely evaporate their liquid water and return to
aerosol/haze residuals; surviving wet particles retain their size.

### Common size exponent

Both modes use:

```python
CHAMBER_BL_EVAP_SIZE_EXP = p
```

with the command-line override:

```bash
--bl-evap-size-exp p
```

The old:

```bash
--bl-inhom-size-exp
```

name is retained only as a deprecated alias.

Interpretation:

```text
p = 0
    homogeneous: equal fractional liquid-mass shrinkage
    inhomogeneous: uniform complete-particle fraction

p = 2
    homogeneous: common finite D^2-like decrement
    inhomogeneous: inverse-D^2 lifetime weighting

p > 2
    increasingly favours complete evaporation of smaller wet particles
    in the inhomogeneous treatment
```

---

## 14. Wall-water closure modes

The generated BMM namelist supports:

```text
chamber_bl_wall_water_mode = 0, 1 or 2
```

and the driver exposes:

```bash
--bl-wall-water-mode {0,1,2}
```

### Mode 0: legacy saturation-cap/effective closure

Mode 0 retains the historical chamber saturation-adjustment closure used for
reproducibility and chamber sensitivity studies.

It does **not** use the physical wall vapour-transfer velocity.

Therefore:

```bash
--wall-vapour-transfer-velocity-ms
```

has no effect when `--bl-wall-water-mode 0` is selected.

### Mode 1: finite reservoir + fractional relaxation

Mode 1 uses a prognostic finite wall-water reservoir and fractional relaxation
toward wall equilibrium.

The main strength parameter is:

```python
CHAMBER_WALL_WATER_EFFICIENCY
```

The wall cannot evaporate/sublimate more water than is stored in its prognostic
liquid/frost reservoir.

### Mode 2: finite-rate physical vapour mass transfer

Mode 2 uses the measured wall temperature and a vapour mass-transfer velocity:

```python
CHAMBER_WALL_VAPOUR_TRANSFER_VELOCITY
```

with an area-mean flux based on vapour-pressure disequilibrium:

```text
Jv = km * [e_eq(Twall) - e_air] / (Rv * Tgas)
```

In this mode wall-vapour exchange is independent of `CHAMBER_BL_TAU`.

Initial reservoir controls are:

```python
CHAMBER_WALL_LIQUID_WATER_INIT_KG
CHAMBER_WALL_ICE_WATER_INIT_KG
```

The reservoir quantities are physical chamber masses in kg.

---

## 15. Fan-associated particle loss

The optional chamber fan loss is a size-dependent first-order particle sink.

The configured functional form is:

```text
kfan(D) = kmax / [1 + (D50/D)^p]
```

with an RPM-dependent reference size.

Controls are:

```python
CHAMBER_FAN_LOSS
CHAMBER_FAN_LOSS_KMAX
CHAMBER_FAN_LOSS_D50_REF
CHAMBER_FAN_LOSS_EXP
CHAMBER_FAN_RPM
CHAMBER_FAN_RPM_REF
```

The same survival factor is applied to particle number and associated extensive
particle moments.

This is separate from wall deposition and gravitational fallout.

---

## 16. Non-gravitational wall deposition

The wall-loss switch is:

```python
CHAMBER_WALL_LOSS = 1
```

with the smooth-wall deposition treatment controlled by:

```python
CHAMBER_WALL_USTAR
CHAMBER_DIAMETER
CHAMBER_HEIGHT
```

The current configuration uses:

```text
diameter = 1.5 m
height   = 2.15 m
```

The configured non-gravitational wall treatment represents side-wall/ceiling
deposition separately from generic gravitational fallout.

---

## 17. Fallout / sedimentation

Generic BMM sedimentation/fallout is controlled by:

```python
FALLOUT_FLAG = True
RESIDENCE_DEPTH = CHAMBER_HEIGHT
```

For chamber runs, the updated BMM uses the chamber geometry to determine the
well-mixed floor-loss scale, with `CHAMBER_HEIGHT` corresponding to
`V/A_floor` for the cylindrical chamber.

Fan loss, non-gravitational wall loss and fallout are independent switches.

---

## 18. Synthetic-updraft mode

The configuration includes:

```python
SYNTHETIC_UPDRAFT = False
```

When enabled, synthetic-updraft mode disables measured chamber pressure and
temperature forcing and chamber-specific BL/fan/wall processes, while leaving
generic BMM controls such as fallout separately configurable.

---

## 19. Dust activation and ice nucleation

The maintained workflow supports dust-containing AIDAd experiments, including
SDSA01 and ATD03 components configured in `iskylab_config.py`.

### FHH activation

Dust FHH coefficients are selected by component:

```python
DUST_FHH = {
    "SDSA01": {"A": 2.25, "B": 1.20},
    "ATD03":  {"A": 0.27, "B": 0.79},
}
```

and written to BMM through:

```fortran
afhh_core1(...)
bfhh_core1(...)
```

The Python driver supplies the component-specific parameters; the activation
calculation remains in BMM.

### INP category

Dust identity for explicit active-site calculations is passed through:

```fortran
inp_category(...)
```

using configuration such as:

```python
DUST_INP_CATEGORY = {
    "SDSA01": "sdsa01",
    "ATD03":  "atd03",
}
```

The analytical `n_s(T)`/IASD expressions live in BMM, not in this repository.

### Ice nucleation mechanisms

For configured dust experiments the driver can enable ice and write:

```fortran
ice_nucleation_mech(1:4)
```

from:

```python
DUST_ENABLE_ICE
DUST_ICE_NUCLEATION_MECH
```

The current configuration uses Koop + INAS for the maintained dust cases and
does not enable DeMott for those components, avoiding double counting of the
explicit dust IASD reservoir.

---

## 20. Variable-resolution INP temperature classes

`INP_TEMP_C` is no longer restricted to the historical 16-point grid.

The configuration may contain any finite, strictly warm-to-cold sequence, for
example:

```python
import numpy as np

INP_TEMP_C = list(np.linspace(-18.0, -33.0, 50))
```

or:

```python
INP_TEMP_C = list(np.linspace(-18.0, -30.0, 100))
```

The driver automatically writes matching namelist values:

```fortran
n_inp_classes = 50
inp_temp(1:50) = ...
```

or:

```fortran
n_inp_classes = 100
inp_temp(1:100) = ...
```

The old upper bound in the namelist template is detected and replaced
dynamically, so a template containing:

```fortran
inp_temp(1:16)
```

does not restrict the run to 16 classes.

The Python validation requires:

- at least one threshold;
- all values finite;
- strictly decreasing temperatures from warm to cold.

For example:

```text
-18, -18.3, -18.6, ..., -33 degC
```

is valid.

---

## 21. Single-experiment model/observation comparison

Running:

```bash
python3 dataAnalysis_new.py --experiment ExpNNN
```

creates the main model/observation figure plus PSD and metrics diagnostics.

The main comparison contains:

1. chamber pressure;
2. gas/wall temperature;
3. hydrometeor water;
4. particle number;
5. effective diameter;
6. relative dispersion;
7. chamber water-exchange/reservoir diagnostics;
8. a compact metrics summary.

### Pressure and temperature

Measured forcing and BMM output are plotted together. Measured wall temperature
is shown when available.

### Liquid and ice water

Two liquid-water definitions are retained intentionally:

```text
BMM native total ql
WELAS reader total LWC-derived ql
```

and a WELAS-equivalent thresholded quantity reconstructed from the PSD/model
particle distribution.

The selected comparison is controlled by:

```python
SINGLE_COMPARE_QL_MODE = "above_min"
```

with supported values:

```text
above_min
total
```

Ice water `qi` is also plotted when available.

### Number

The analysis keeps distinct:

```text
WELAS supplied Nd
WELAS PSD-integrated N above Dmin
BMM activated liquid Nd
BMM warm-particle N above Dmin
```

Ice number is plotted on a secondary axis when available.

### Effective diameter and dispersion

For mixed-phase experiments, the direct comparison to the merged WELAS PSD uses
combined BMM liquid + ice moments:

```text
Mk,total = sum_liquid(N D^k) + sum_ice(N D^k)
```

with:

```text
Deff = M3/M2
```

and relative dispersion from the combined first and second moments.

Liquid-only and activated-liquid diagnostics remain visible separately.

---

## 22. WELAS-equivalent threshold diagnostics

The default minimum wet diameter is:

```python
SINGLE_COMPARE_DROP_MIN_UM = 2.5
```

WELAS spectra are treated as:

```text
dN/dlog10(D)
```

using reconstructed logarithmic bin edges and their actual `dlogD` widths.

The observational analysis derives:

- number above `Dmin`;
- mean diameter;
- volume-mean diameter;
- effective diameter;
- relative dispersion;
- liquid-equivalent mass above `Dmin`.

The WELAS-equivalent liquid retrieval treats each optical/wet diameter as a
spherical water particle, matching the instrument-style calculation.

For mixed-phase experiments this quantity must be interpreted cautiously:
large observed particles may be ice even though a liquid-equivalent diagnostic
can still be calculated mathematically.

---

## 23. Fixed-bin threshold treatment in BMM diagnostics

For full-moving BMM output (`bin_scheme_flag=0`), the diagnostic wet-diameter
threshold is a pointwise particle test.

For fixed water-mass schemes (`bin_scheme_flag=1` or `2`), if the diagnostic
diameter threshold lies inside a native water-mass bin, the analysis uses the
fraction of that bin above the corresponding water-mass threshold rather than
an all-or-nothing bin-centre decision.

This improves consistency of number and water diagnostics across BMM binning
schemes.

---

## 24. WELAS/BMM PSD comparison

The WELAS has many more optical channels than the BMM has independent particle
populations. Direct comparison on the native WELAS grid would therefore make
the model distribution unnecessarily sparse.

Both distributions are compared on a common coarse logarithmic grid configured
by:

```python
SINGLE_COMPARE_PSD_NBINS = 48
SINGLE_COMPARE_PSD_MIN_UM = None
SINGLE_COMPARE_PSD_MAX_UM = None
```

`None` uses the WELAS bounds for the experiment.

The observed `dN/dlog10D` is conservatively rebinned by exact overlap in
log-diameter space.

The BMM particle populations are histogrammed onto the same diagnostic bins
using their complete number weights.

The colour limits are based on the observed WELAS concentration range so that
matching colours correspond to matching concentrations rather than each panel
being independently rescaled.

### Mixed-phase PSD

When ice is available, the plotted BMM comparison distribution is:

```text
warm/liquid: nwat with dwet
ice:         nicem with dmaxice
```

so the model panel represents liquid + ice hydrometeors.

WELAS provides optical size, whereas BMM ice is currently represented using
`dmaxice`; this distinction should be retained when interpreting the comparison.

---

## 25. Complete model-only PSD

The model-only figure:

```text
model-full-psd-ExpNNN.png
```

covers a much wider size range than the WELAS comparison.

The warm panel uses:

```text
nwat + dwet
```

and therefore includes:

- dry/wet aerosol represented by the warm population;
- haze;
- activated droplets;
- residual warm particles.

It is not restricted to the activated `nliq` subset.

When ice is present, a second panel uses:

```text
nicem + dmaxice
```

The fixed log-D grids are post-processing grids only. They do not change the
native BMM bin representation.

Current controls include:

```python
SINGLE_MODEL_PSD_NBINS = 180

SINGLE_MODEL_WARM_PSD_MIN_UM = 0.005
SINGLE_MODEL_WARM_PSD_MAX_UM = 1000.0

SINGLE_MODEL_ICE_PSD_MIN_UM = 0.1
SINGLE_MODEL_ICE_PSD_MAX_UM = 10000.0
```

---

## 26. Liquid-water lag and optional WELAS time shift

The analysis searches for a bounded ql lag using:

```python
SINGLE_COMPARE_MAX_LAG_S = 30.0
SINGLE_COMPARE_LAG_STEP_S = 1.0
```

A positive diagnosed lag means the observation occurs later than the model.

By default this lag is reported only.

With:

```bash
--shift-cloud-measurements
```

the whole WELAS time coordinate is shifted by the diagnosed amount before
comparison, keeping all cloud measurements internally aligned.

Both absolute and lag-related liquid-water scores are retained.

---

## 27. Comparison metrics

For each single experiment the driver writes:

```text
metrics-ExpNNN.csv
```

including quantities such as:

- comparison window;
- model saturation target;
- selected ql comparison mode;
- diagnosed ql lag and correlation;
- applied WELAS time shift;
- ql NRMSE;
- total-qL NRMSE;
- thresholded-qL NRMSE;
- thresholded number NRMSE;
- total-hydrometeor Deff NRMSE;
- total-hydrometeor relative-dispersion NRMSE;
- integrated model/observed liquid-water exposure.

Where wall-reservoir diagnostics exist, the figure summary also reports
reservoir masses.

---

## 28. Chamber water and wall diagnostics

The main comparison can plot cumulative BMM diagnostics such as:

```text
qfan_liq
qwall_liq
qfall_liq

qchamber_bl
qchamber_bl_evap

qchamber_wall_liq_evap
qchamber_wall_liq_cond
qchamber_wall_ice_subl
qchamber_wall_ice_dep
```

Prognostic wall reservoirs are plotted separately in physical grams:

```text
chamber_wall_liquid_water
chamber_wall_ice_water
```

When available, the dedicated wall-vapour figure also plots:

```text
chamber_wall_rh
chamber_wall_vapour_flux
```

with the vapour flux displayed in `mg m-2 s-1`.

---

## 29. Experiment groups

Batch groups and their updraft/settings metadata are defined in:

```text
scripts/experiment_metadata.py
```

Run a group with:

```bash
python3 dataAnalysis_new.py --group N
```

Group runs use the same aerosol initialisation and namelist-generation
machinery as single runs.

When no `--group` or `--experiment` is supplied, the historical `THIS_RUN`
default in `dataAnalysis_new.py` is used.

---

## 30. Adding a new experiment

A new experiment normally requires:

1. add timing/group metadata to `experiment_metadata.py`;
2. add the meteorology/CPC source to `readMeteoCPC.py`;
3. add the WELAS/OPC source to `readOPC_Merged.py`;
4. add the initial aerosol PSD definition to `readPNSD_Mrg_new.py`;
5. keep parallel reader arrays index-aligned;
6. add new aerosol composition/FHH metadata if required;
7. add a new `inp_category` and corresponding BMM `n_s(T)` implementation for
   a new ice-active dust species;
8. test the experiment in single-run mode before adding it to a group.

Example:

```bash
python3 dataAnalysis_new.py --experiment ExpNNN
```

---

## 31. Recommended workflow

For a model-development or chamber sensitivity run:

```text
1. Build/update BMM.
2. Set BMM_MODEL_FOLDER.
3. Edit iskylab_config.py for persistent settings.
4. Run one representative experiment.
5. Inspect the generated namelist.
6. Inspect forcing-smoothing diagnostics.
7. Inspect the main model/observation figure.
8. Inspect the common-grid PSD comparison.
9. Inspect the complete model PSD.
10. Inspect wall/fan/fallout diagnostics where relevant.
11. Only then run a larger experiment group or batch script.
```

The generated namelist is the best record of the actual model inputs used for
a specific run.

---

## 32. Troubleshooting

### BMM executable not found

Check:

```bash
echo "$BMM_MODEL_FOLDER"
ls "$BMM_MODEL_FOLDER/main.exe"
```

### A single experiment tries to open another experiment

Single mode should only load the requested experiment. Check that the
experiment names and filename/metadata arrays in the reader modules remain
index-aligned.

### Generated `inp_temp` still has the wrong size

Check:

```bash
grep -n "n_inp_classes\|inp_temp" \
    /tmp/$USER/single_comparison/namelist-Exp021.in
```

If the config contains 50 classes, the generated namelist should show:

```text
n_inp_classes = 50
inp_temp(1:50) = ...
```

Also confirm which configuration Python imported:

```bash
python3 - <<'PY'
import iskylab_config as cfg
print(cfg.__file__)
print(len(cfg.INP_TEMP_C))
print(cfg.INP_TEMP_C[:3], cfg.INP_TEMP_C[-3:])
PY
```

### BMM reports `Cannot match namelist object name -34` or another number

This commonly indicates that more `inp_temp` values were written than the
declared namelist slice can accept. Regenerate the namelist with the current
dynamic-INP driver and confirm that `n_inp_classes` and `inp_temp(1:N)` have the
same `N`.

### Wall vapour-transfer velocity appears to do nothing

Check `chamber_bl_wall_water_mode`.

`CHAMBER_WALL_VAPOUR_TRANSFER_VELOCITY` is used by wall-water mode 2, not legacy
mode 0.

### Mixed-phase WELAS-derived liquid-equivalent water becomes very large

The merged WELAS PSD can include ice. A `D^3` spherical-water reconstruction is
then only an instrument-equivalent diagnostic, not a physical liquid-water
measurement. Use the distinct liquid and total-hydrometeor comparisons
appropriately.

### Model output stops early

Inspect the generated:

```text
runtime
```

and the available measured forcing interval. The driver currently caps a single
run at the smaller of the experiment forcing duration and 30 minutes.

### Batch script waits at figures

For unattended operation omit `--plt-ion`, or use:

```bash
--no-plot
```

to suppress displayed figures entirely.

---

## 33. Reproducibility

For a simulation used in a publication, report or model intercomparison, retain:

- the exact BMM commit/version;
- the exact iSKYLAB-analysis commit/version;
- `scripts/iskylab_config.py`;
- the generated BMM namelist;
- the relevant experiment metadata;
- any command-line overrides;
- the BMM NetCDF output;
- the comparison metrics;
- the generated figures.

In particular, do not rely only on the persistent Python configuration when
command-line overrides were used. The generated namelist records the actual BMM
configuration passed to `main.exe`.
