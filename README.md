# iSKYLAB analysis and BMM batch driver

This repository reads the iSKYLAB experiment files in `./iSKYLAB-data`, fits the
initial aerosol PSDs, generates BMM namelists, runs experiment batches, and
compares model output with the chamber measurements.

## Current BMM chamber interface

The maintained entry point is `scripts/dataAnalysis_new.py`.  Model/chamber
controls are collected in `scripts/iskylab_config.py`; the analysis script no
longer edits exact copied strings from one historical namelist.

The batch driver supports the current BMM chamber controls:

- independent pressure, gas-temperature and total-water forcing;
- chamber BL mixing `0=off`, `1=homogeneous`, `2=extreme inhomogeneous`;
- BL temperature from either a fixed gas-wall offset or the measured
  `Tww_mean` wall-temperature time series;
- sigmoid drone-fan particle loss (`kmax`, `D50`, exponent, RPM);
- Lai-Nazaroff non-gravitational wall deposition (`u*` plus chamber geometry);
- independent generic gravitational fallout/sedimentation.

The chamber geometry defaults to a 2.0 m diameter × 2.5 m high cylinder.  The
same 2.5 m height is written as the generic residence depth for clarity; in a
chamber run the updated BMM uses `chamber_height` as `V/A_floor` for fallout.

## Running

Set the BMM location either by editing `scripts/iskylab_config.py` or, preferably,
with an environment variable:

```bash
export BMM_MODEL_FOLDER=/path/to/bmm
```

Then from `scripts/` run, for example:

```bash
python dataAnalysis_new.py --group 6
```

Useful options are:

```text
--group N       select an experiment group from experiment_metadata.py
--no-run        analyse existing /tmp outputs without rerunning BMM
--no-analysis   generate/run namelists only
--no-plot       suppress the summary plot
```

Generated namelists and NetCDF files are written to `/tmp/$USER` by default.


## Chamber forcing smoothing

Short-period noise in the measured gas temperature can create unrealistically
large alternating heating/cooling tendencies in the BMM and repeatedly move
aerosol/droplets across saturation and activation thresholds.  The batch driver
therefore supports Savitzky-Golay smoothing of the **forcing variables only**.
It is enabled by default in `scripts/iskylab_config.py`:

```python
SMOOTH_CHAMBER_FORCING = True
CHAMBER_SMOOTH_TEMP_WINDOW = 21.0
CHAMBER_SMOOTH_WALL_TEMP_WINDOW = 21.0
CHAMBER_SMOOTH_PRESSURE_WINDOW = 21.0
CHAMBER_SMOOTH_POLYORDER = 2
```

The windows are in seconds, not samples.  The code infers the median sampling
interval for each experiment; if timestamps are noticeably irregular it
interpolates to a uniform median-dt grid for filtering and then maps the result
back onto the original time coordinate.  SciPy's `mode="interp"` treatment is
used at the ends to avoid artificial padding.

The following working fields are smoothed before namelist generation:

- `Tgw mean` -> `temp_chamber`;
- `Tww mean` -> `wall_temp_chamber` when BL temperature mode 1 is used;
- `Pressure` -> `press_chamber`.

Dew point and particle observations are **not** smoothed by this option.  The
raw P/T observations are retained in memory as `*_raw`, with the gap-filled
pre-filter signal retained as `*_filled`.

With `SAVE_FORCING_SMOOTHING_DIAGNOSTICS=True`, the driver writes one diagnostic
plot per experiment plus `forcing-smoothing-summary.csv` under
`/tmp/$USER/forcing_smoothing/`.  Each plot compares raw and smoothed gas/wall
temperatures, explicitly shows `Twall-Tgas`, and compares raw/smoothed pressure.
The CSV reports the effective window, sampling interval, timestamp irregularity,
and RMS/maximum smoothing corrections.  Check these diagnostics when changing
the window; 11, 21 and 31 s are sensible sensitivity values for nominal 1 Hz
records.

To generate only these diagnostics without running or analysing the BMM, use:

```bash
python dataAnalysis_new.py --no-run --no-analysis --no-plot
```

## Water-variable convention

`TDew` is converted to a **vapour mixing ratio**.  It is not silently labelled
as total water.  The default diagnostic `qtot_chamber` is

```text
qtot = qv(dew point) + ql(OPC)
```

with OPC LWC interpolated onto the chamber meteo time grid.  If the dew-point
instrument is known to measure an evaporated total-water sample instead, set
`QTOT_DATA_MODE = "dewpoint_only"` explicitly.

`qtot_chamber` can be written without being used as forcing.  `FORCE_QTOT` is
off by default because the total-water interpretation is less certain and
because forcing it while also applying chamber sinks can double-count water
loss.

## Initial RH and aerosol normalisation

Two choices are exposed in `iskylab_config.py`:

- `INITIAL_RH_METHOD="cloud_onset"` chooses initial vapour so the measured P/T
  trajectory reaches saturation at a configurable model target time.  With
  `MODEL_SATURATION_TIME_SHIFT_S=0` and an empty
  `MODEL_SATURATION_TIME_OVERRIDES_S`, this is the historical observed cloud
  onset from `experiment_metadata.py`.  The model target can be shifted without
  moving the observational onset marker/window;
- `INITIAL_RH_METHOD="dewpoint"` uses the measured t=0 dew point directly.

Aerosol normalisation can use CPC at t=0 or at cloud onset.  `t0` is recommended
when fan/wall/fallout losses are active because using cloud-onset CPC and then
modelling pre-cloud loss would double represent that loss.

## Important code fixes in this revision

- Removed the obsolete `chamber_override` workflow.
- Removed all hard-coded `n_levels_c=1853` and giant chamber-array replacement
  strings.  `&chamber_spec` is regenerated to the actual experiment length.
- Namelist variables are edited by name and missing variables now raise errors
  instead of silently leaving template values unchanged.
- Fixed fallback PNSD diameter units: fitted/fallback diameters are held in µm
  and converted to metres exactly once.
- Fixed NaN filling for chamber water/dew-point series; missing values are
  interpolated rather than all being assigned the value before the first NaN.
- Corrected the air-density parentheses and double-density factor in
  `mixingAnalysis.py`, and corrected the plot-axis labels.
- Made OPC effective diameter safe for empty spectra.
- Updated `readModel.py` to expose the new BL/fan/wall/fallout diagnostics when
  present.
- Fixed dormant `svp.py` errors in the WMO power terms, Clausius latent heats,
  Python-3 range use and invalid-method error handling.
- Data paths are resolved from the repository rather than from the current
  working directory.

## Dependencies

See `requirements.txt`.  The repository does not include the experiment data
or the BMM executable/build dependencies.

## Single-experiment model/observation comparison

The driver now supports a dedicated development/validation mode, for example:

```bash
python dataAnalysis_new.py --experiment Exp005
```

or simply `--experiment 5`.  A one-off saturation-timing sensitivity can be
requested without editing metadata, for example:

```bash
python dataAnalysis_new.py --experiment Exp005 --saturation-time-min 6.6
```

This changes the **model saturation target** used to infer the initial vapour
amount.  It does not move the observed cloud-onset marker or curated comparison
window.  For persistent or batch sensitivities use
`MODEL_SATURATION_TIME_SHIFT_S` and `MODEL_SATURATION_TIME_OVERRIDES_S` in
`iskylab_config.py`.

The single-run workflow generates one named namelist/output under
`/tmp/$USER/single_comparison/` and writes:

- `comparison-Exp005.png`: pressure and temperature forcing, liquid-water
  mixing ratio, droplet number, effective diameter, relative dispersion and
  cumulative chamber/fallout loss diagnostics;
- `psd-Exp005.png`: observed WELAS and BMM wet-particle size distributions
  conservatively rebinned onto the **same configurable coarse log-D grid and
  colour scale**;
- `metrics-Exp005.csv`: cloud-window scores including ql, number, effective
  diameter and relative dispersion, plus a bounded diagnostic ql lag.

The absolute-time comparison is always retained.  The lag calculation is a
separate diagnostic only (`SINGLE_COMPARE_MAX_LAG_S`, default +/-30 s) intended
to expose likely sample-line/instrument delay rather than silently align the
model and observations.  Positive reported lag means the observation occurs
later than the model.

OPC spectra are treated as `dN/dlog10D`.  Their actual logarithmic bin widths
are reconstructed from the diameter centres rather than using the historical
hard-coded `0.009839` factor.  The same widths are used to calculate OPC
`Dmean`, volume-mean diameter, `Deff=M3/M2`, relative dispersion and an
instrument-like droplet number above `SINGLE_COMPARE_DROP_MIN_UM` (2 um by
default).

The updated BMM writes each native warm-bin number together with its current
wet diameter.  Those native bins are **not** interpreted as fixed wet-size
intervals.  For the direct WELAS/BMM figure, the optical PSD is first
number-conservatively integrated onto a configurable coarse common log-D grid,
and `(nwat,dwet)` is histogrammed onto those same edges before division by
`Delta log10D`.  This avoids giving the BMM artificial WELAS-scale resolution
while retaining an apples-to-apples diagnostic grid.

For number concentration, the single-run plot deliberately shows both the BMM
activation diagnostic and an OPC-like BMM count above the selected wet-diameter
threshold.  This makes instrument-threshold effects visible instead of folding
them into an activation comparison.

### Complete model PSD (including unactivated particles)

The direct OPC comparison remains instrument-limited, but the single-run mode
also writes `model-full-psd-EXP.png`.  Its warm panel uses **all** `nwat` with
`dwet`, so unactivated aerosol, haze and activated droplets occupy one
continuous model distribution.  `nliq` is not used for this full PSD; it is
retained only as an activated-only diagnostic.  A second panel shows the ice
population from `nicem` with `dmaxice` when ice is present.

The plotting grids are fixed logarithmic diagnostic grids configured by
`SINGLE_MODEL_*_PSD_*` in `iskylab_config.py`.  They are post-processing grids
only and do not imply that the native BMM wet diameters have fixed bin widths.

The historical batch interface is retained unchanged:

```bash
python dataAnalysis_new.py --group 6
```

`--group` and `--experiment` are mutually exclusive.  If neither is supplied,
the historical `THIS_RUN` group is used.

### Coarse common WELAS/BMM comparison grid

The WELAS has many more optical size channels than there are independent BMM
moving particle populations.  The direct comparison therefore no longer puts
each BMM population into the native WELAS channels.  Instead both datasets are
number-conservatively rebinned to a common logarithmic grid controlled by:

```python
SINGLE_COMPARE_PSD_NBINS = 48
SINGLE_COMPARE_PSD_MIN_UM = None
SINGLE_COMPARE_PSD_MAX_UM = None
```

`None` uses the measured WELAS lower/upper edge for that experiment.  The
model-only full PSD plot remains separate and still uses the broader
`SINGLE_MODEL_*_PSD_*` grid.
