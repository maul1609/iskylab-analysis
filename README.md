# iSKYLAB analysis and BMM simulation driver

This repository contains the analysis and model-driving workflow used to compare
iSKYLAB/AIDAd chamber experiments with the **Bin Microphysics Model (BMM)**.

The maintained entry point is:

```text
scripts/dataAnalysis_new.py
```

The workflow can:

- read the chamber meteorology, CPC, WELAS/OPC and initial aerosol PSD data;
- smooth the measured chamber forcing used by the model;
- fit/construct the initial aerosol size distributions;
- generate a BMM namelist for an individual experiment or experiment group;
- run the externally compiled BMM executable;
- read the BMM NetCDF output;
- compare modelled and observed water, number concentration and size-distribution evolution;
- diagnose chamber wall, fan, boundary-layer and sedimentation losses;
- analyse warm, mixed-phase and ice-containing experiments.

The repository does **not** contain the BMM source tree or the iSKYLAB data.

---

## 1. Repository layout

The important files are:

```text
scripts/
    dataAnalysis_new.py        main driver
    iskylab_config.py          user-facing configuration
    experiment_metadata.py     experiment groups/timing metadata
    readMeteoCPC.py            chamber meteorology/CPC reader
    readOPC_Merged.py          WELAS/OPC reader
    readPNSD_Mrg_new.py        initial aerosol PSD reader/fitter
    readModel.py               BMM NetCDF reader
    namelist-aida-*.in         canonical BMM namelist template

iSKYLAB-data/
    Datasets-V2/
        ...                    experiment data (not distributed here)
```

The analysis code expects `iSKYLAB-data` to sit at the repository root unless
paths are changed in `iskylab_config.py`.

---

## 2. Requirements

### Python

Install the Python requirements from the repository, for example:

```bash
python3 -m pip install -r requirements.txt
```

A recent Python 3 installation with NumPy, SciPy, pandas, matplotlib and NetCDF
support is required by the current workflow.

### BMM

The BMM is a separate repository:

```text
https://github.com/UoM-maul1609/bin-microphysics-model
```

Build the BMM using its own build system before attempting a model run. The
analysis repository only generates the namelist and launches the compiled BMM;
it does not compile the Fortran source automatically.

Tell iSKYLAB where the BMM build is located either in
`scripts/iskylab_config.py` or, preferably, with:

```bash
export BMM_MODEL_FOLDER=/path/to/bin-microphysics-model
```

The path must point to the BMM installation/build location expected by the
configuration in `iskylab_config.py`.

---

## 3. Quick start: run one experiment

Run commands from the `scripts/` directory.

For example, to run Exp006:

```bash
cd scripts
python3 dataAnalysis_new.py --experiment Exp006
```

The short form is also accepted:

```bash
python3 dataAnalysis_new.py --experiment 6
```

For a dust experiment such as Exp021:

```bash
python3 dataAnalysis_new.py --experiment Exp021
```

Single-experiment mode reads only the requested experiment, generates its
namelist, runs BMM, reads the output and makes the comparison diagnostics.

### Useful command-line options

```text
--experiment ExpNNN     run one experiment
--group N               run a group from experiment_metadata.py
--no-run                analyse an existing BMM output without rerunning BMM
--no-analysis           generate/run the namelist but skip analysis
--no-plot               suppress the summary plots
--saturation-time-min X override the model saturation target for one run
```

`--group` and `--experiment` are mutually exclusive. If neither is supplied,
the historical `THIS_RUN` value in the configuration is used.

Examples:

```bash
# run and analyse Exp005
python3 dataAnalysis_new.py --experiment Exp005

# rerun Exp021 but do not produce analysis
python3 dataAnalysis_new.py --experiment Exp021 --no-analysis

# analyse an existing Exp006 output without rerunning BMM
python3 dataAnalysis_new.py --experiment Exp006 --no-run

# one-off saturation timing sensitivity
python3 dataAnalysis_new.py --experiment Exp005 --saturation-time-min 6.6
```

---

## 4. Output locations

Single-experiment outputs are written beneath:

```text
/tmp/$USER/single_comparison/
```

The exact filenames include the experiment identifier. Typical outputs are:

```text
namelist-Exp006.in
comparison-Exp006.png
psd-Exp006.png
model-full-psd-Exp006.png
metrics-Exp006.csv
```

The generated BMM NetCDF file is also retained under the temporary model-output
location used by the driver.

Forcing-smoothing diagnostics are written to:

```text
/tmp/$USER/forcing_smoothing/
```

including one diagnostic plot per experiment and:

```text
forcing-smoothing-summary.csv
```

---

## 5. Main configuration

Most user-adjustable settings are kept in:

```text
scripts/iskylab_config.py
```

This includes:

- BMM location;
- pressure and temperature forcing switches;
- total-water diagnostic/forcing behaviour;
- forcing smoothing;
- model saturation timing;
- aerosol normalisation time;
- chamber boundary-layer treatment;
- fan loss;
- wall deposition;
- sedimentation/fallout;
- chamber geometry;
- model/observation comparison thresholds;
- PSD diagnostic grids;
- lag-diagnostic settings.

This is the file to edit for normal sensitivity studies. Avoid changing
`dataAnalysis_new.py` for settings that are already exposed in the config.

---

## 6. Chamber forcing

The current BMM chamber interface treats the major chamber processes separately.
The driver supports independent forcing of:

- pressure;
- gas temperature;
- total water, if explicitly enabled.

The measured chamber time series are written to the generated `&chamber_spec`
namelist section at the actual experiment length; there is no fixed hard-coded
number of chamber points.

### Pressure and temperature

Pressure and gas temperature forcing can be independently enabled. Forcing is
normally based on the measured chamber time series.

### Total water

Dew point is converted to a **vapour mixing ratio**. It is not silently treated
as total water.

The default diagnostic total-water estimate is:

```text
qtot = qv(dew point) + ql(OPC)
```

with OPC liquid water interpolated onto the meteorological time grid.

`qtot_chamber` can be written diagnostically without being imposed on BMM.
`FORCE_QTOT` is normally kept off because forcing measured total water while
also applying explicit wall/BL/fan losses can double-count chamber water loss.

---

## 7. Chamber forcing smoothing

Short-period noise in measured chamber temperature can produce unrealistic
alternating supersaturation and evaporation tendencies. The model-driving
variables can therefore be smoothed before they are written to the namelist.

Typical settings are:

```python
SMOOTH_CHAMBER_FORCING = True
CHAMBER_SMOOTH_TEMP_WINDOW = 21.0
CHAMBER_SMOOTH_WALL_TEMP_WINDOW = 21.0
CHAMBER_SMOOTH_PRESSURE_WINDOW = 21.0
CHAMBER_SMOOTH_POLYORDER = 2
```

The windows are in **seconds**, not samples.

The following forcing fields are smoothed:

```text
Tgw mean   -> temp_chamber
Tww mean   -> wall_temp_chamber
Pressure   -> press_chamber
```

Dew point and particle observations are not smoothed by this option.

The code retains raw, gap-filled and smoothed versions of the forcing so that
the filtering can be checked explicitly.

To generate only smoothing diagnostics:

```bash
python3 dataAnalysis_new.py --no-run --no-analysis --no-plot
```

---

## 8. Initial humidity and aerosol normalisation

### Initial humidity

Two initial-humidity methods are exposed:

```python
INITIAL_RH_METHOD = "cloud_onset"
```

or:

```python
INITIAL_RH_METHOD = "dewpoint"
```

With `cloud_onset`, the initial vapour amount is chosen so that the measured
pressure/temperature trajectory reaches saturation at the configured model
target time. The target can be adjusted without moving the observed cloud-onset
marker used by the analysis.

The persistent controls are:

```python
MODEL_SATURATION_TIME_SHIFT_S = 0.0
MODEL_SATURATION_TIME_OVERRIDES_S = {}
```

### Aerosol normalisation

Aerosol concentration can be normalised using CPC at `t0` or at cloud onset.
When fan/wall/fallout losses are explicitly represented, `t0` is generally the
cleaner choice because it avoids using a later measured concentration and then
modelling the pre-cloud loss a second time.

---

## 9. Chamber boundary-layer treatment

The chamber boundary-layer process is controlled independently from wall
particle deposition and gravitational fallout.

The current BMM modes are:

```text
chamber_bl_mix = 0   off
chamber_bl_mix = 1   homogeneous boundary-layer processing
chamber_bl_mix = 2   extreme inhomogeneous processing
chamber_bl_mix = 3   D2-weighted extreme inhomogeneous processing
```

Mode 2 removes the same fraction of activated droplets from each size bin by
complete evaporation, leaving survivors at their existing size.

Mode 3 is also an extreme-inhomogeneous treatment but weights complete
evaporation toward smaller droplets using the diffusional lifetime scaling
associated with the D-squared evaporation law.

The boundary-layer air temperature can be based on either:

- a fixed offset from the bulk gas temperature; or
- the measured wall temperature `Tww_mean`.

Processed boundary-layer air mixes back into the chamber on the configured BL
timescale.

If complete evaporation releases residual aerosol, that behaviour is governed
by the BMM `release_aerosol` switch.

---

## 10. Fan, wall and fallout losses

These processes are independent and can be switched on/off separately.

### Fan-associated particle loss

The chamber fan loss uses a saturating sigmoid size dependence:

```text
kfan(D) = kmax / [1 + (D50/D)^p]
```

with optional RPM scaling of `D50`. The main controls are in
`iskylab_config.py`, including:

```text
kmax
D50_ref
p
RPM
reference RPM
```

The same survival fraction is applied to particle number and associated
extensive particle moments.

### Non-gravitational wall deposition

Wall loss can use the Lai-Nazaroff smooth-wall treatment. Chamber geometry and
friction velocity are configurable. Side-wall and ceiling loss are treated as
non-gravitational deposition.

### Gravitational fallout

Generic sedimentation/fallout remains separate from the chamber wall-loss
model. For the AIDAd geometry the chamber height is 2.5 m, which is also the
well-mixed sedimentation depth `V/A_floor` for the cylindrical chamber.

---

## 11. Dust, FHH activation and ice-active-site-density runs

The current analysis supports dust-containing experiments in addition to the
warm NaCl experiments.

The maintained dust set currently includes:

```text
Exp021   SDSA01
Exp022   SDSA01 + NaCl
Exp023   SDSA01 + NaCl
Exp025   ATD03
Exp026   ATD03 + NaCl
Exp027   ATD03 + NaCl
```

Exp020 can be added when the corresponding data files/metadata are available.
Exp024/SDTP02 is not currently part of the maintained simulation set.

### FHH parameters

Adsorption/FHH activation coefficients are written directly into the BMM
namelist using the existing component fields, for example:

```fortran
afhh_core1(...)
bfhh_core1(...)
```

This means the experiment-driving Python code selects the aerosol component and
writes its FHH coefficients, while the activation calculation itself remains in
BMM.

### IASD / n_s parameterisation

For explicit ice-active-site-density calculations the **dust identity** is
specified in the namelist using `inp_category`. The analytical `n_s(T)` fits
belong in `bin_microphysics_module.f90`, not in the Python analysis code.

The current dust categories include:

```text
sdsa01
atd03
```

The active-site temperature classes are also read from the namelist. The current
1 K grid spans 255--240 K and is supplied in **degrees Celsius**, as expected by
BMM:

```fortran
n_inp_classes = 16
inp_temp(1:16) = &
   -18.15, -19.15, -20.15, -21.15, &
   -22.15, -23.15, -24.15, -25.15, &
   -26.15, -27.15, -28.15, -29.15, &
   -30.15, -31.15, -32.15, -33.15
```

When adding a new dust `inp_category`, BMM must map the namelist string to an
internal INP type, recognise it as an explicit INAS/IASD species, and evaluate
its `n_s(T)` expression in `ice_active_site_density()`.

**Important:** dust-specific `n_s(T)` fits should be verified against the
appropriate laboratory dataset before being treated as final. The analysis
repository should not duplicate those fits.

---

## 12. Single-experiment comparison plots

Running:

```bash
python3 dataAnalysis_new.py --experiment ExpNNN
```

produces a main comparison figure containing the principal chamber/model
variables.

### Water

For warm experiments the figure includes observed and modelled liquid-water
quantities.

For mixed-phase/ice experiments the main water panel also includes **BMM ice
water (`qi`)**, so liquid and ice mass evolution can be inspected together.

The observed WELAS total LWC is kept distinct from any liquid-equivalent mass
reconstructed from the measured PSD. In a mixed-phase experiment, the merged
WELAS PSD must not automatically be interpreted as an all-liquid distribution.

If a time-shifted observational ql diagnostic is plotted, it is labelled by the
actual quantity, for example:

```text
WELAS PSD >2 um shifted -30s
```

rather than the ambiguous label `selected obs`.

The lag is diagnostic only and does not modify the BMM simulation or the
unshifted comparison.

### Number concentration

The droplet-number panel retains the liquid/droplet comparison on its main
axis. For ice-containing experiments a secondary axis also shows:

```text
WELAS Ni
BMM Ni
```

so ice number can be compared without compressing it onto the much larger
liquid-number scale.

### Effective diameter and relative dispersion

The analysis keeps the liquid diagnostics but also computes **total
hydrometeor** moments by combining the BMM liquid and ice populations.

For a total moment of order `k`:

```text
Mk,total = sum_liquid(N D^k) + sum_ice(N D^k)
```

The total effective diameter is then:

```text
Deff,total = M3,total / M2,total
```

and the total relative dispersion is calculated from the combined first and
second diameter moments.

These total quantities are useful for comparison with the merged WELAS/OPC
size distribution after ice forms. Liquid-only model quantities remain plotted
as additional diagnostics.

---

## 13. PSD comparison with WELAS/OPC

The WELAS has many more optical channels than there are independent moving BMM
particle populations. Direct comparison is therefore performed on a common
coarse logarithmic diameter grid.

The grid is configured with:

```python
SINGLE_COMPARE_PSD_NBINS = 48
SINGLE_COMPARE_PSD_MIN_UM = None
SINGLE_COMPARE_PSD_MAX_UM = None
```

`None` uses the measured WELAS lower/upper bounds for the experiment.

The observed PSD is conservatively rebinned in log-diameter space. BMM particle
populations are histogrammed onto the same diagnostic bins, conserving particle
number before division by `Delta log10(D)`.

### Mixed-phase PSDs

For ice-containing experiments the BMM distribution compared with the merged
WELAS PSD is the sum of:

```text
liquid/warm particles: nwat with dwet
ice particles:         nicem with dmaxice
```

so the plotted BMM PSD represents the **total liquid + ice hydrometeor
population**, rather than liquid alone.

The current model/observation diameter mapping should still be interpreted with
care: WELAS provides an optical particle size whereas BMM ice is currently
placed using `dmaxice`.

---

## 14. Complete model PSD

The model-only figure:

```text
model-full-psd-ExpNNN.png
```

is intentionally broader than the direct WELAS comparison.

Its warm-particle panel uses all `nwat` together with `dwet`, including:

- unactivated aerosol;
- haze particles;
- activated droplets;
- residual warm particles.

`nliq` is retained as an activated-liquid diagnostic and is not used to define
the complete warm PSD.

A separate ice panel uses:

```text
nicem + dmaxice
```

when ice is present.

The model PSD grids are post-processing grids only; they do not imply that BMM
has fixed wet-diameter bins.

---

## 15. WELAS/OPC moments and threshold diagnostics

OPC/WELAS spectra are treated as `dN/dlog10(D)`. Their actual logarithmic bin
widths are reconstructed from the reported diameter centres rather than using a
hard-coded constant width.

The same bin widths are used to calculate observational diagnostics such as:

- mean diameter;
- volume-mean diameter;
- `Deff = M3/M2`;
- relative dispersion;
- number above the selected wet-diameter threshold.

The default comparison threshold is controlled by:

```python
SINGLE_COMPARE_DROP_MIN_UM = 2.0
```

For warm experiments, the same threshold can be used to construct an
instrument-like BMM liquid-water/number comparison from `nwat,dwet`.

For mixed-phase experiments, a PSD-derived `D^3` mass integral over all WELAS
particles must **not** be interpreted as liquid water because the measured PSD
can contain ice.

---

## 16. Batch runs

Experiment groups are defined in:

```text
scripts/experiment_metadata.py
```

Run a group with:

```bash
python3 dataAnalysis_new.py --group N
```

Group mode uses the same namelist-generation and BMM-driving machinery as
single-experiment mode.

When adding a new experiment, keep the experiment ordering consistent across
all reader metadata. In particular, any index-based filename and metadata
arrays in the reader modules must refer to exactly the same experiment at a
given index.

---

## 17. Adding a new experiment

A typical new experiment requires updates in the following places:

1. Add the experiment metadata/timing to `experiment_metadata.py`.
2. Add the corresponding meteorology/CPC file to `readMeteoCPC.py`.
3. Add the WELAS/OPC file to `readOPC_Merged.py`.
4. Add the initial aerosol PSD information to `readPNSD_Mrg_new.py`.
5. Ensure all parallel arrays in each reader remain index-aligned.
6. Add any new aerosol composition/FHH metadata required for namelist generation.
7. If it is a new ice-active dust, add the corresponding `inp_category` and
   `n_s(T)` implementation to BMM.
8. Test first with single-experiment mode:

```bash
python3 dataAnalysis_new.py --experiment ExpNNN
```

Only after the single run works should the experiment be added to a batch group.

---

## 18. Recommended development workflow

For a new or modified chamber configuration:

```text
1. Compile/update BMM.
2. Set BMM_MODEL_FOLDER.
3. Edit iskylab_config.py for the desired chamber physics.
4. Run one representative experiment.
5. Inspect the generated namelist.
6. Inspect forcing-smoothing diagnostics.
7. Check the main model/observation comparison.
8. Check the common-grid PSD comparison.
9. Check the complete model PSD.
10. Only then run the full experiment group.
```

For example:

```bash
cd scripts
python3 dataAnalysis_new.py --experiment Exp006
```

followed by a batch run when satisfied:

```bash
python3 dataAnalysis_new.py --group 6
```

---

## 19. Notes on model/observation interpretation

Several plotted quantities are intentionally kept separate rather than forced
into a single definition:

- BMM physical `ql` versus an instrument-like mass reconstructed from
  `nwat,dwet`;
- WELAS total LWC versus a PSD-derived mass integral;
- activated BMM liquid number versus an OPC-like count above a wet-diameter
  threshold;
- liquid-only moments versus total liquid+ice hydrometeor moments;
- WELAS optical diameter versus BMM ice `dmaxice`.

This is deliberate. The comparison code is intended to expose instrument and
phase-definition differences rather than hide them inside a single metric.

---

## 20. Current chamber geometry

The default AIDAd chamber geometry used by the driver is a cylinder with:

```text
diameter = 2.0 m
height   = 2.5 m
```

The geometry is used by the chamber wall/fallout treatments. Keep these values
consistent between `iskylab_config.py` and the generated BMM namelist when
performing chamber simulations.

---

## 21. Troubleshooting

### A single experiment tries to open another experiment's file

Single-experiment mode should only load the requested experiment. If an Exp021
run attempts to open Exp020, check that the filename and metadata arrays in the
reader modules are aligned. Do not leave a filename for an omitted experiment
in one parallel array while removing its label from another.

### BMM does not run

Check:

```bash
echo "$BMM_MODEL_FOLDER"
```

and confirm the compiled BMM executable is present where
`iskylab_config.py` expects it.

Also inspect the generated namelist under `/tmp/$USER/single_comparison/`.

### Model output ends before the observations

Check the runtime written into the generated namelist and any runtime cap in the
analysis configuration/code. The run duration should normally cover the full
requested experimental time series.

### Mixed-phase WELAS-derived liquid water looks far too large

Do not interpret the entire merged WELAS PSD as spherical liquid once ice is
present. Large ice particles make a `D^3` liquid-equivalent integral enormous.
Use the supplied WELAS LWC for liquid-water comparison and use the merged PSD
for total hydrometeor size/moment comparisons.

---

## 22. Reproducibility

For a simulation intended for publication or comparison, retain:

- the exact BMM commit/version;
- the exact iSKYLAB analysis commit/version;
- the generated BMM namelist;
- the relevant `iskylab_config.py` settings;
- the experiment metadata;
- the BMM NetCDF output;
- the comparison metrics and figures.

The generated namelist is especially important because it records the actual
model configuration used for the run.
