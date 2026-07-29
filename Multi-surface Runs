# User Note: Bailey (2026) Three-Surface Runoff Model

## 1. Purpose

The Python script **`Bailey2026 RunoffModel_3surface.py`** implements a nonlinear-reservoir runoff model in which each subcatchment can contain up to three independently parameterised surface types:

* `PAVED`
* `ROOF`
* `ADDITIONAL`

Each active surface is assigned its own area, runoff coefficient, slope, α parameter and β parameter. Rainfall is converted to inflow separately for each surface, routed through the nonlinear-reservoir model, and then combined to obtain:

1. the total outflow from each subcatchment; and
2. the total outflow from the complete modelled catchment.

The script can be run either as a single baseline simulation or as a parallel α–β parameter sweep.

---

## 2. Required files

The supplied files are:

```text
Bailey2026 RunoffModel_3surface.py
subcatchments_3_surfaces.csv
Raingauges(1).csv
```

Place all three files in the same folder.

### Important filename check

The script currently refers to:

```python
CATCHMENT_PARAMETERS_FILE = Path(
    "./CatchmentA_subcatchments_3_surfaces 2 only-fixedparams.csv"
)
RAINFALL_CSV_PATH = Path("./Raingauges.csv")
```

These names do not match the supplied CSV filenames. Either rename the CSV files or change the script settings to:

```python
CATCHMENT_PARAMETERS_FILE = Path("./subcatchments_3_surfaces.csv")
RAINFALL_CSV_PATH = Path("./Raingauges(1).csv")
```

The script must be run from the folder containing these files unless full file paths are provided.

---

## 3. Python requirements

The script requires Python 3 and the following packages:

```text
numpy
pandas
scipy
```

Install them from a command prompt or terminal using:

```bash
pip install numpy pandas scipy
```

The Python standard-library modules used for file handling, timing and parallel processing do not require separate installation.

---

## 4. Model structure

For every active surface, the script calculates:

```text
k = √slope / α
m = 5/3 + β
```

The surface is then represented by a nonlinear storage–discharge relationship using the calculated values of `k` and `m`.

The exponent convention is explicitly fixed in the script as:

```python
m = (5.0 / 3.0) + beta
```

The sign of β should therefore not be reversed when preparing the parameter CSV or defining parameter sweeps.

A surface is active only when its area is greater than zero. For example:

```text
ADDITIONAL_area_ha = 0
```

causes the `ADDITIONAL` surface to be omitted from routing for that subcatchment.

Although the model supports three surface types, the supplied subcatchment CSV activates only the `PAVED` and `ROOF` surfaces. All supplied `ADDITIONAL_area_ha` values are zero.

---

## 5. Subcatchment input file

The file **`subcatchments_3_surfaces.csv`** defines the characteristics of each subcatchment and its three potential surface types.

It must contain the following general columns:

| Column            | Description                                        |
| ----------------- | -------------------------------------------------- |
| `subcatchment_id` | Unique subcatchment identifier                     |
| `rain_gauge`      | Rainfall-gauge column assigned to the subcatchment |

It must also contain the following five columns for each surface type:

```text
[SURFACE]_area_ha
[SURFACE]_runoff_coeff
[SURFACE]_slope
[SURFACE]_alpha
[SURFACE]_beta
```

The complete required surface columns are therefore:

```text
PAVED_area_ha
PAVED_runoff_coeff
PAVED_slope
PAVED_alpha
PAVED_beta

ROOF_area_ha
ROOF_runoff_coeff
ROOF_slope
ROOF_alpha
ROOF_beta

ADDITIONAL_area_ha
ADDITIONAL_runoff_coeff
ADDITIONAL_slope
ADDITIONAL_alpha
ADDITIONAL_beta
```

Additional columns, such as `FM Ref` and `Comments`, may be included and are ignored by the model.

### Supplied subcatchment file

The supplied file contains six subcatchments:

```text
A
B
C
D
E
F
```

All six subcatchments:

* contain active paved and roof surfaces;
* assign rainfall gauge `RG2`;
* set the additional-surface area to zero.

The total active areas in the supplied file are approximately:

```text
Paved area: 0.458 ha
Roof area:  0.316 ha
```

### Input requirements

* Each `subcatchment_id` should be unique.
* Surface areas must be supplied in hectares.
* A surface with an area greater than zero must have valid values for all its other parameters.
* Runoff coefficients would normally lie between 0 and 1.
* Slopes must not be negative.
* α must be greater than zero.
* The value in `rain_gauge` must exactly match a rainfall-column heading.
* The model strips leading and trailing spaces from subcatchment and rainfall-gauge identifiers.
* Rows without a subcatchment identifier or rainfall-gauge assignment are removed.
* A subcatchment with no positive surface areas is skipped with a warning.

Even when the `ADDITIONAL` surface is not used, its required columns must remain present in the CSV.

---

## 6. Rainfall input file

The rainfall file must contain:

```text
time_minutes
```

followed by one or more rainfall-gauge columns, for example:

```text
time_minutes,RG1,RG2,RG3
0,0,0,0
2,0,1.2,0
4,0,2.0,0
```

Rainfall values must be supplied as **intensities in millimetres per hour**.

The script converts rainfall intensity to inflow using:

* the area assigned to the individual surface;
* the surface runoff coefficient; and
* the rainfall gauge assigned to the subcatchment.

Rainfall is linearly interpolated between the supplied time points. Rainfall before the first value and after the final value is set to zero.

### Supplied rainfall file

The supplied `Raingauges(1).csv` file contains:

* `time_minutes`;
* rainfall gauges `RG1` to `RG13`;
* values at two-minute intervals; and
* a rainfall record extending from 0 to 48 minutes.

The final supplied rainfall time is therefore:

```text
48 minutes = 0.8 hours
```

All supplied subcatchments use `RG2`, which is present in the rainfall file.

### Missing rainfall data

The script applies:

```python
df.dropna()
```

to the complete rainfall table. A missing value in any rainfall-gauge column therefore removes that entire time row, even where the affected gauge is not being used.

Missing rainfall values should be checked and corrected before running the model.

---

## 7. Simulation time settings

The principal time settings are:

```python
S0 = 0.0
START_TIME = 0.0
END_TIME = 1488.0
TIME_STEP = 1.0 / 30.0
```

These represent:

| Setting      | Meaning                        |
| ------------ | ------------------------------ |
| `S0`         | Initial storage                |
| `START_TIME` | Simulation start time in hours |
| `END_TIME`   | Simulation end time in hours   |
| `TIME_STEP`  | Output interval in hours       |

The supplied time step is:

```text
1/30 hour = 2 minutes
```

### Important duration check

The script is currently configured to run for 1,488 hours, whereas the supplied rainfall file covers only 0.8 hours. Rainfall will consequently be set to zero after 0.8 hours.

For a simulation covering only the supplied rainfall record, change:

```python
END_TIME = 0.8
```

For another rainfall dataset, `END_TIME` would normally be set using:

```text
final time_minutes value ÷ 60
```

A longer period may be retained where a post-event recession period is deliberately required.

---

## 8. Numerical solver settings

The model attempts the following numerical solvers in order:

```python
SOLVER_METHODS_TO_TRY = ["LSODA", "BDF", "Radau", "RK45"]
```

`LSODA` is used first. If it fails for a particular surface and parameter combination, the model attempts the remaining methods.

The default tolerances are:

```python
SOLVER_RTOL = 1e-5
SOLVER_ATOL = 1e-8
SOLVER_MAX_STEP = TIME_STEP
```

Negative storage values arising from numerical error are clipped to zero. The model does not introduce an artificial positive storage floor or associated dry-weather outflow.

These solver settings should normally be retained unless a numerical sensitivity test is being undertaken.

---

## 9. Running a single baseline simulation

Set:

```python
MODE = "single_run"
```

The model will use the α and β values contained in the subcatchment CSV for every active surface.

Run the script from a terminal opened in the model folder:

```bash
python "Bailey2026 RunoffModel_3surface.py"
```

On systems where Python is invoked using `python3`, use:

```bash
python3 "Bailey2026 RunoffModel_3surface.py"
```

The output is saved as:

```text
routing_results_3surface/baseline_csv_alpha_beta.csv
```

The standard output file contains:

| Column                           | Description                                                    |
| -------------------------------- | -------------------------------------------------------------- |
| `Time_hr`                        | Elapsed simulation time in hours                               |
| One column for each subcatchment | Combined outflow from all active surfaces in that subcatchment |
| `Total_Outflow_m3s`              | Sum of all subcatchment outflows                               |

The subcatchment columns are labelled using their `subcatchment_id` values.

---

## 10. Saving surface-specific results

By default:

```python
SAVE_SURFACE_DETAIL_COLUMNS = False
```

The saved file therefore contains only subcatchment totals and the overall catchment total.

To save outflow from each individual surface, change this to:

```python
SAVE_SURFACE_DETAIL_COLUMNS = True
```

Additional columns will then be added using names such as:

```text
A_PAVED_Outflow_m3s
A_ROOF_Outflow_m3s
A_ADDITIONAL_Outflow_m3s
```

A surface-specific column is created only for an active surface.

Saving surface details increases output-file size, particularly for long simulations and large parameter sweeps.

---

## 11. Output precision and compression

Saved numerical values are rounded using:

```python
OUTPUT_DECIMALS = 5
```

This setting affects only the saved CSV files. Internal model calculations retain their normal numerical precision.

To change the number of saved decimal places, edit this value, for example:

```python
OUTPUT_DECIMALS = 8
```

By default:

```python
COMPRESS_CSV_OUTPUTS = False
```

To save compressed CSV files, use:

```python
COMPRESS_CSV_OUTPUTS = True
```

The output filenames will then end in:

```text
.csv.gz
```

Compressed files can be read directly by pandas and many spreadsheet or data-analysis programs.

---

## 12. Running a parallel parameter sweep

Set:

```python
MODE = "parallel_sweep"
```

The script provides two different ways to define the sweep:

```python
SWEEP_MODE = "explicit_incomplete_runs"
```

or:

```python
SWEEP_MODE = "normal_arrays"
```

Only the settings associated with the selected sweep mode are used.

---

## 13. Explicit incomplete-run mode

The supplied script currently uses:

```python
SWEEP_MODE = "explicit_incomplete_runs"
```

In this mode, the script ignores the general α and β arrays and runs only the parameter combinations listed in:

```python
INCOMPLETE_RUNS
```

Each row must contain four values in the following order:

```text
PAVED α, PAVED β, ROOF α, ROOF β
```

For example:

```python
INCOMPLETE_RUNS = np.array([
    [0.011, -0.3, 0.004, -0.9],
    [0.011, -0.3, 0.006, -0.9],
    [0.011, -0.3, 0.010, -0.9],
], dtype=float)
```

The supplied settings therefore run three parameter combinations.

This mode is intended for running a known list of outstanding combinations without recreating an entire Cartesian parameter grid.

Duplicate parameter combinations are detected and cause the sweep to stop with an error.

---

## 14. Normal-array sweep mode

To construct a sweep from the α and β arrays, set:

```python
SWEEP_MODE = "normal_arrays"
```

The surface-specific arrays are:

```python
PAVED_ALPHA_VALUES
PAVED_BETA_VALUES

ROOF_ALPHA_VALUES
ROOF_BETA_VALUES

ADDITIONAL_ALPHA_VALUES
ADDITIONAL_BETA_VALUES
```

The script forms all α–β pairs within each active surface list.

The surfaces varied by the sweep are controlled by:

```python
SURFACES_TO_SWEEP = ["PAVED", "ROOF"]
```

Surfaces not included in this list retain the α and β values from the subcatchment CSV.

To include the third surface, use:

```python
SURFACES_TO_SWEEP = ["PAVED", "ROOF", "ADDITIONAL"]
```

However, varying the additional-surface parameters has no effect where all `ADDITIONAL_area_ha` values are zero.

---

## 15. Sweep combination modes

When `SWEEP_MODE` is set to `normal_arrays`, the surface parameter pairs can be combined in three ways.

### Cartesian mode

```python
SWEEP_COMBINATION_MODE = "cartesian"
```

Every parameter pair for one surface is combined with every parameter pair for the other active surfaces.

This can generate a large number of runs. The expected number should be checked before execution.

### Linked-by-index mode

```python
SWEEP_COMBINATION_MODE = "linked_by_index"
```

The first paved pair is combined with the first roof pair, the second paved pair with the second roof pair, and so on.

All active surface-pair lists must contain the same number of entries.

### One-surface-at-a-time mode

```python
SWEEP_COMBINATION_MODE = "one_surface_at_a_time"
```

One surface type is varied while the other surface types retain the parameter values contained in the subcatchment CSV.

This mode is useful for separately examining paved-, roof- and additional-surface sensitivity.

### Interaction with explicit-run mode

`SWEEP_COMBINATION_MODE` is ignored when:

```python
SWEEP_MODE = "explicit_incomplete_runs"
```

---

## 16. Restricting a sweep to selected subcatchments

By default:

```python
TARGET_SUBCATCHMENTS = None
```

The sweep parameters are applied to every active surface of the selected type in all subcatchments.

To change parameters only in selected subcatchments, provide their identifiers:

```python
TARGET_SUBCATCHMENTS = ["A", "B"]
```

Non-targeted subcatchments retain their CSV α and β values.

Identifiers must match the values in the `subcatchment_id` column.

---

## 17. Parallel-processing settings

The number of worker processes is currently defined by:

```python
MAX_WORKERS = os.cpu_count() - 2
```

This attempts to leave two processor cores unused.

A fixed value can be used to reduce memory demand or maintain computer responsiveness:

```python
MAX_WORKERS = 4
```

For safer operation on computers with a small number of processor cores, the setting may be written as:

```python
MAX_WORKERS = max(1, (os.cpu_count() or 2) - 2)
```

Progress is reported according to:

```python
REPORT_EVERY = 5
```

The parameter sweep is processed in small batches to improve stability, particularly on Windows.

The script should be executed as a complete Python file rather than by running isolated sections in an interactive console. The required multiprocessing protection is already included under:

```python
if __name__ == "__main__":
```

---

## 18. Parameter-sweep outputs

Individual parameter-run files are saved in:

```text
3surface_alpha_beta_sweep_outputs2/
```

Filenames identify the surface parameters used, for example:

```text
paved_a0.011_b-0.30__roof_a0.004_b-0.90.csv
```

Each successful file contains:

* simulation time;
* total outflow from each subcatchment;
* optional surface-specific outflows; and
* total catchment outflow.

A sweep summary is saved as:

```text
3surface_alpha_beta_sweep_outputs2/
3surface_alpha_beta_sweep_summary.csv
```

The summary records:

* run name;
* success or failure status;
* output filename or error details;
* paved α and β where varied;
* roof α and β where varied; and
* additional-surface α and β where varied.

The summary is written even where a sweep is interrupted or some runs fail.

---

## 19. Output folders and overwriting

The script automatically creates:

```text
routing_results_3surface
3surface_alpha_beta_sweep_outputs2
```

Existing output files with identical names may be overwritten without a separate confirmation. Results that need to be retained should be copied, renamed or moved before the model is rerun.

The output-folder names can be changed in the user settings:

```python
RESULTS_FOLDER = Path("./routing_results_3surface")
SWEEP_OUTPUT_FOLDER = Path("./3surface_alpha_beta_sweep_outputs2")
```

---

## 20. Common errors

### Catchment file not found

Check that:

* the CSV is in the working folder;
* the filename in `CATCHMENT_PARAMETERS_FILE` is correct; and
* the script is being run from the intended directory.

The filename in the supplied script does not currently match the supplied subcatchment CSV.

### Rainfall file not found

Check the value assigned to:

```python
RAINFALL_CSV_PATH
```

The supplied file is named `Raingauges(1).csv`, whereas the script currently refers to `Raingauges.csv`.

### Missing columns in catchment parameter CSV

The error message lists the required columns that could not be found. Check spelling, capitalisation and underscores.

Surface names must use the prefixes:

```text
PAVED_
ROOF_
ADDITIONAL_
```

### Rainfall CSV must contain `time_minutes`

Confirm that the time column is named exactly:

```text
time_minutes
```

### Gauge not found in rainfall data

A value in the subcatchment `rain_gauge` column does not match a rainfall-column heading.

Check for:

* spelling differences;
* additional spaces;
* missing gauges; or
* changes in capitalisation.

### Invalid α

Every active surface must have:

```text
α > 0
```

### Negative slope

Surface slopes must not be negative. A zero slope is accepted by the script but produces `k = 0` and therefore no modelled discharge.

### Invalid exponent

The script requires:

```text
5/3 + β > 0
```

A sufficiently negative β value causes the run to fail.

### No active surfaces

At least one surface in each included subcatchment must have an area greater than zero.

### Broken process pool

This can occur where a worker process is terminated during a large parallel sweep. Possible responses include:

* reducing `MAX_WORKERS`;
* shortening the simulation period;
* reducing the number of simultaneous parameter combinations;
* disabling surface-detail columns;
* enabling compressed output; or
* dividing the sweep into smaller explicit run lists.

### Very large output files

Output size increases with:

* simulation duration;
* number of subcatchments;
* number of parameter combinations;
* output frequency;
* saved decimal precision; and
* inclusion of surface-detail columns.

Before a large sweep, check `END_TIME`, `TIME_STEP`, `OUTPUT_DECIMALS` and `SAVE_SURFACE_DETAIL_COLUMNS`.

---

## 21. Recommended workflow

1. Place the Python script and both CSV files in the same folder.
2. Correct the two input filenames in the script.
3. Check that the rainfall duration and `END_TIME` are consistent.
4. Confirm that every subcatchment rainfall-gauge assignment exists in the rainfall CSV.
5. Set `MODE = "single_run"`.
6. Run the baseline model and inspect the combined hydrograph.
7. Enable `SAVE_SURFACE_DETAIL_COLUMNS` where separate paved and roof responses need to be checked.
8. Confirm the exponent convention `m = 5/3 + β`.
9. Define the required parameter combinations.
10. Set `MODE = "parallel_sweep"` only after the baseline simulation completes successfully.
11. Check the sweep summary for failed runs.
12. Retain copies of the script, input CSVs and settings with the outputs so that each simulation can be reproduced.
