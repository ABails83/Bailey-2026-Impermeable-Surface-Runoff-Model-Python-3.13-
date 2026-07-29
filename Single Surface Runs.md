# Bailey-2026-Impermeable-Surface-Runoff-Model-Python-3.13-
# User Note: Bailey (2026) Single-Surface Runoff Model

## 1. Purpose

The Python script **`Bailey2026 RunoffModel - single surface.py`** implements a nonlinear-reservoir runoff model for one or more urban subcatchments. Rainfall inputs are converted to subcatchment inflows and routed using the specified surface slope and calibrated α–β parameters.

The script can be used in two modes:

1. **Single run:** uses the α and β values assigned to each subcatchment in `Subcatchments.csv`.
2. **Parallel parameter sweep:** runs multiple combinations of α and β and saves a separate combined hydrograph for each parameter pair.

---

## 2. Required files

Place the following three files in the same folder:

```text
Bailey2026 RunoffModel - single surface.py
Subcatchments.csv
Raingauges.csv
```

The script uses relative file paths:

```python
CATCHMENT_PARAMETERS_FILE = Path("./Subcatchments.csv")
RAINFALL_CSV_PATH = Path("./Raingauges.csv")
```

The code should therefore be run from the folder containing these files. Alternatively, these paths may be replaced with full file paths.

---

## 3. Python requirements

The script requires Python 3 and the following packages:

```text
numpy
pandas
scipy
matplotlib
```

These packages can be installed from a terminal or command prompt using:

```bash
pip install numpy pandas scipy matplotlib
```

`matplotlib` is required even when plot generation is disabled because it is imported when the script starts.

---

## 4. Subcatchment input file

The file **`Subcatchments.csv`** defines the properties and model parameters of each subcatchment. It must contain the following columns:

| Column            | Description                                              |
| ----------------- | -------------------------------------------------------- |
| `subcatchment_id` | Unique name or identifier for the subcatchment           |
| `area_ha`         | Subcatchment contributing area in hectares               |
| `runoff_coeff`    | Dimensionless runoff coefficient                         |
| `slope`           | Surface slope expressed as a decimal                     |
| `alpha`           | Calibrated α parameter used in a single run              |
| `beta`            | Calibrated β parameter used in a single run              |
| `rain_gauge`      | Name of the rainfall column assigned to the subcatchment |

An example row is:

```text
subcatchment_id,area_ha,runoff_coeff,slope,alpha,beta,rain_gauge
A,0.100,1,0.04,0.01,-0.5,RG2
```

Additional columns, such as flow-monitor references or comments, may be included and are ignored by the model.

### Input requirements

* Each `subcatchment_id` should be unique.
* `area_ha` should be greater than zero.
* `runoff_coeff` would normally be between 0 and 1.
* `slope` must not be negative.
* `alpha` must be greater than zero.
* The value in `rain_gauge` must exactly match a column heading in `Raingauges.csv`, including spelling and capitalisation.
* Rows containing missing values in any required column are removed during import.

---

## 5. Rainfall input file

The file **`Raingauges.csv`** contains the rainfall time series. The first column must be:

```text
time_minutes
```

All other columns represent rainfall gauges, for example:

```text
time_minutes,RG1,RG2,RG3,RG4
0,0,0,0,0
2,0,0,0,0
4,0,1.2,0,0
```

Rainfall values must be supplied as **intensities in millimetres per hour**. The script internally converts these values to metres per second and then to subcatchment inflow in cubic metres per second using:

* the subcatchment area;
* the runoff coefficient; and
* the selected rainfall-gauge column.

### Rainfall-file requirements

* `time_minutes` must contain numerical elapsed times in minutes.
* Times should increase sequentially.
* Rainfall-gauge names must agree with those specified in `Subcatchments.csv`.
* Rainfall cells should be numerical.
* Missing data should be filled or removed before running the model. The script applies `dropna()` to the complete rainfall table, meaning that a missing value in any gauge column causes the entire time row to be removed.
* Rainfall outside the period covered by the CSV is set to zero.

The supplied rainfall example extends to **48 minutes, or 0.8 hours**. The current script is configured to simulate to 1,488 hours. Consequently, rainfall is zero after 0.8 hours when the supplied example file is used. For a simulation covering only the example rainfall period, change:

```python
END_TIME = 0.8
```

For other datasets, `END_TIME` should normally be set to the final `time_minutes` value divided by 60.

---

## 6. Main user settings

The main settings are located near the beginning of the script under:

```python
# USER SETTINGS
```

### Run mode

Select one of the following:

```python
MODE = "single_run"
```

or:

```python
MODE = "parallel_sweep"
```

### Simulation period

The model time settings are:

```python
S0 = 0.0
START_TIME = 0.0
END_TIME = 1488.0
TIME_STEP = 1.0 / 30.0
```

These settings mean:

* initial storage is zero;
* the simulation starts at 0 hours;
* the simulation ends at 1,488 hours; and
* results are reported every 1/30 hour, equivalent to 2 minutes.

The output interval does not need to be identical to the rainfall interval because rainfall is linearly interpolated between the supplied values.

The model should be used with the unit convention under which it was calibrated. Rainfall is supplied in millimetres per hour, the model time vector is expressed in hours, and the generated flow columns are labelled in cubic metres per second. Changes to these conventions should only be made after checking the governing equation and parameter units.

---

## 7. Running a single simulation

Set:

```python
MODE = "single_run"
```

The model then uses the individual `alpha` and `beta` values contained in each row of `Subcatchments.csv`.

From a terminal opened in the model folder, run:

```bash
python "Bailey2026 RunoffModel - single surface.py"
```

On some systems, the command may instead be:

```bash
python3 "Bailey2026 RunoffModel - single surface.py"
```

The script reports each subcatchment as it is processed and creates the following output:

```text
routing_results_TSR/combined_outflow_hydrographs.csv
```

This file contains:

| Column                      | Description                                     |
| --------------------------- | ----------------------------------------------- |
| `Time_hr`                   | Elapsed simulation time in hours                |
| One column per subcatchment | Modelled outflow from that subcatchment in m³/s |
| `Total_Outflow_m3s`         | Sum of all modelled subcatchment outflows       |

### Optional individual outputs

To save inflow, outflow and storage results separately for every subcatchment, change:

```python
SAVE_INDIVIDUAL_RESULTS = True
```

The files are saved in:

```text
routing_results_TSR/
```

Each individual file contains:

```text
Time_hr
Inflow_m3s
Outflow_m3s
Storage_m3
```

### Optional plots

To generate an inflow–outflow plot for each subcatchment, change:

```python
MAKE_PLOTS = True
```

Plots are saved as PNG files in:

```text
routing_plots_TSR/
```

---

## 8. Running an α–β parameter sweep

Set:

```python
MODE = "parallel_sweep"
```

The values tested are defined by:

```python
alpha_values = ...
beta_values = ...
```

The supplied settings test 14 α values and 9 β values, giving:

```text
14 × 9 = 126 simulations
```

By default:

```python
target_subcatchments = None
```

This applies each α–β pair to every subcatchment.

To apply the sweep only to selected subcatchments, specify their identifiers, for example:

```python
target_subcatchments = ["A", "B"]
```

When:

```python
preserve_non_target_values = True
```

subcatchments not included in the target list retain the α and β values given in `Subcatchments.csv`.

The identifiers in `target_subcatchments` must use the same data type and spelling as those in the CSV.

### Parallel processing

The number of worker processes is set by:

```python
MAX_WORKERS = max(1, (os.cpu_count() or 2) - 1)
```

This normally uses all but one available processor core. A smaller number may be specified manually where computer memory or responsiveness is a concern, for example:

```python
MAX_WORKERS = 4
```

The sweep should be run as a complete Python script rather than by executing isolated sections in an interactive environment. The script already contains the required:

```python
if __name__ == "__main__":
```

guard for multiprocessing.

### Sweep outputs

Each successful parameter combination produces a file in:

```text
alpha_beta_sweep_outputs/
```

Files are named using the tested parameters, for example:

```text
alpha_0.010_beta_0.40.csv
```

Each file contains the outflow hydrographs for all subcatchments and their total.

A summary file is also generated:

```text
alpha_beta_sweep_outputs/alpha_beta_sweep_summary.csv
```

The summary records:

* run name;
* α value;
* β value;
* success or failure status; and
* output filename or error message.

Parameter sweeps can generate many large files. The simulation duration, output interval, number of parameter combinations and available disk space should therefore be checked before starting a sweep.

---

## 9. Output folders

The following folders are created automatically if they do not already exist:

```text
routing_results_TSR
routing_plots_TSR
alpha_beta_sweep_outputs
```

Existing files with the same names may be overwritten without an additional warning. Previous results should therefore be copied or renamed before rerunning the model where they need to be retained.

---

## 10. Common errors

### “Catchment parameter file not found”

Check that:

* `Subcatchments.csv` is in the working folder;
* its filename agrees with the path in the script; and
* the script is being run from the correct directory.

### “Rainfall data could not be loaded”

Check that:

* `Raingauges.csv` is present;
* the file contains `time_minutes`; and
* all required values are numerical and correctly formatted.

### “Gauge not found in rainfall file”

A `rain_gauge` entry in `Subcatchments.csv` does not exactly match a rainfall-column heading. Check for differences in spelling, capitalisation or spaces.

### “No subcatchments loaded”

One or more required columns may be absent, or all rows may contain missing required values. Confirm that the required column names have not been changed.

### Invalid or unstable results

Check that:

* α is greater than zero;
* slope is not negative;
* area and runoff coefficient are reasonable;
* the simulation period agrees with the rainfall record;
* α and β use the same calibration and unit conventions as the model; and
* the rainfall values are intensities in millimetres per hour rather than rainfall depths per interval.

---

## 11. Recommended workflow

Before a full model run:

1. Confirm the two CSV filenames and column headings.
2. Check that every subcatchment has a valid rainfall-gauge assignment.
3. Set `END_TIME` to the required simulation duration.
4. Run the model in `single_run` mode.
5. Inspect the combined hydrograph and, where necessary, enable individual results and plots.
6. Run `parallel_sweep` only after the baseline model has completed successfully.
7. Retain a copy of the input CSV files and script settings alongside the generated results to provide a reproducible record of each model run.
