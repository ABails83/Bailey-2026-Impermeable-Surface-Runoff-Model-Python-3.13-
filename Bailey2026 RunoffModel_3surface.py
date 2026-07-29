import itertools
import os
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

# =========================================================
# USER SETTINGS
# =========================================================

# Modes:
#   "single_run"       = run the alpha/beta values already in the catchment CSV.
#   "parallel_sweep"  = run a parameter sweep for the three surface types.
MODE = "single_run"

CATCHMENT_PARAMETERS_FILE = Path("./CatchmentA_subcatchments_3_surfaces 2 only-fixedparams.csv")
RAINFALL_CSV_PATH = Path("./Raingauges.csv")

RESULTS_FOLDER = Path("./routing_results_3surface")
SWEEP_OUTPUT_FOLDER = Path("./3surface_alpha_beta_sweep_outputs2")

RESULTS_FOLDER.mkdir(exist_ok=True)
SWEEP_OUTPUT_FOLDER.mkdir(exist_ok=True)

# Time settings
S0 = 0.0
START_TIME = 0.0
END_TIME = 1488.0
TIME_STEP = 1.0 / 30.0   # hours; 1/30 hr = 2 minutes

# Solver settings
# LSODA is usually more robust than RK45 for these final sweep runs because it
# automatically switches between non-stiff and stiff integration when the
# storage-discharge equation becomes numerically awkward near zero storage.
# The fallback methods are only used if the primary method fails for a run.
SOLVER_METHODS_TO_TRY = ["LSODA", "BDF", "Radau", "RK45"]
SOLVER_RTOL = 1e-5
SOLVER_ATOL = 1e-8
SOLVER_MAX_STEP = TIME_STEP

times = np.arange(START_TIME, END_TIME + TIME_STEP / 2, TIME_STEP)
t_span = (START_TIME, END_TIME)

# ---------------------------------------------------------
# Model exponent
# ---------------------------------------------------------
# Fixed thesis convention: m = 5/3 + beta.
# This is hard-coded below to avoid sign-convention ambiguity.

# ---------------------------------------------------------
# Sweep settings
# ---------------------------------------------------------
# The sweep is defined independently for each surface type. Each list contains
# explicit (alpha, beta) pairs. This mirrors a specific-pairs workflow and avoids
# accidental runs from broad alpha x beta grids unless you deliberately create them.

# Sweep settings
# Options:
#   "normal_arrays"           = use the separate alpha/beta arrays below.
#   "explicit_incomplete_runs" = use INCOMPLETE_RUNS directly and avoid rerunning completed cartesian combinations.
SWEEP_MODE = "explicit_incomplete_runs"

# These arrays are retained for normal/cartesian testing. They are ignored when
# SWEEP_MODE = "explicit_incomplete_runs".
PAVED_ALPHA_VALUES = np.array([0.008, 0.011, 0.014, 0.018])
PAVED_BETA_VALUES  = np.array([-0.7, -0.5, -0.3])

ROOF_ALPHA_VALUES  = np.array([0.004, 0.006, 0.008, 0.010])
ROOF_BETA_VALUES   = np.array([-0.9, -0.7, -0.5, -0.3])

# Dummy values because ADDITIONAL is inactive if ADDITIONAL_area_ha = 0.
ADDITIONAL_ALPHA_VALUES = np.array([0.01])
ADDITIONAL_BETA_VALUES = np.array([0.0])

# Exact incomplete runs only. Columns are:
# [PAVED_alpha, PAVED_beta, ROOF_alpha, ROOF_beta]
INCOMPLETE_RUNS = np.array([
    [0.011, -0.3, 0.004, -0.9],
    [0.011, -0.3, 0.006, -0.9],
    [0.011, -0.3, 0.010, -0.9],

  
], dtype=float)

PAVED_ALPHA_BETA_PAIRS = [
    (round(a, 6), round(b, 6))
    for a in PAVED_ALPHA_VALUES
    for b in PAVED_BETA_VALUES
]

ROOF_ALPHA_BETA_PAIRS = [
    (round(a, 6), round(b, 6))
    for a in ROOF_ALPHA_VALUES
    for b in ROOF_BETA_VALUES
]

ADDITIONAL_ALPHA_BETA_PAIRS = [
    (round(a, 6), round(b, 6))
    for a in ADDITIONAL_ALPHA_VALUES
    for b in ADDITIONAL_BETA_VALUES
]

# Sweep combination mode:
#   "cartesian"       = all PAVED x ROOF x ADDITIONAL combinations.
#                       This can become very large: 6 x 6 x 6 = 216 runs.
#   "linked_by_index" = run the nth PAVED pair with the nth ROOF pair and nth ADDITIONAL pair.
#                       Lists must have equal lengths; 6 pairs = 6 runs.
#   "one_surface_at_a_time" = vary one surface at a time and keep the other surfaces at their CSV values.
#                       Useful for isolating paved/roof/additional sensitivity.
SWEEP_COMBINATION_MODE = "cartesian"

# Surface types to include in sweep. Non-swept surface types keep CSV alpha/beta.
SURFACES_TO_SWEEP = ["PAVED", "ROOF"]

# Optional: restrict changes to selected subcatchments only.
# None = apply sweep pair to every active surface of that type in every subcatchment.
# Example: TARGET_SUBCATCHMENTS = ["SUB001", "SUB002"]
TARGET_SUBCATCHMENTS = None

MAX_WORKERS = os.cpu_count() - 2
REPORT_EVERY = 5

# Save detailed per-surface columns as well as subcatchment totals.
# False keeps files smaller: Time_hr + each subcatchment total + Total_Outflow_m3s.
SAVE_SURFACE_DETAIL_COLUMNS = False

# Output precision settings. These affect saved CSV files only, not internal calculations.
OUTPUT_DECIMALS = 5
COMPRESS_CSV_OUTPUTS = False


# =========================================================
# PRECISION-CONTROLLED OUTPUT
# =========================================================

def prepare_output_dataframe(df, decimals=OUTPUT_DECIMALS):
    """Round numeric columns before writing to disk to reduce CSV file size."""
    out = df.copy()
    numeric_cols = out.select_dtypes(include=[np.number]).columns
    out[numeric_cols] = out[numeric_cols].round(decimals)
    return out


def write_csv_reduced_precision(df, filepath, decimals=OUTPUT_DECIMALS):
    """Write a DataFrame with controlled floating-point precision."""
    filepath = Path(filepath)
    if COMPRESS_CSV_OUTPUTS and filepath.suffix != ".gz":
        filepath = filepath.with_suffix(filepath.suffix + ".gz")
    output_df = prepare_output_dataframe(df, decimals=decimals)
    output_df.to_csv(filepath, index=False, float_format=f"%.{decimals}f")
    return filepath


# =========================================================
# CORE MODEL PHYSICS
# =========================================================

def calculate_outflow(S, k, m):
    """Calculate reservoir outflow without adding artificial storage.

    Negative numerical storage values are clipped to zero, but no positive
    storage floor is applied. This avoids introducing artificial dry-weather
    outflow when runoff depths/storage are very small.
    """
    S = np.asarray(S, dtype=float)
    S_eff = np.maximum(S, 0.0)
    Q = k * np.power(S_eff, m)
    return Q.item() if Q.shape == () else Q


def non_linear_reservoir_ode(t, S, I_func, k, m):
    """ODE with a physical zero-storage constraint.

    The discharge is calculated from max(S, 0), with no artificial positive
    storage floor. If the solver steps slightly below zero during dry periods,
    the derivative is constrained so storage does not continue decreasing.
    """
    I = float(I_func(t))
    S_eff = max(float(S[0]), 0.0)
    Q = float(calculate_outflow(S_eff, k, m))
    dSdt = I - Q

    if S_eff <= 0.0 and dSdt < 0.0:
        dSdt = 0.0

    return np.array([dSdt])


def calculate_exponent(beta):
    """Return the storage-discharge exponent using the fixed thesis convention.

    The model uses m = 5/3 + beta for all surfaces and all runs.
    """
    return (5.0 / 3.0) + float(beta)


# =========================================================
# DATA IMPORT
# =========================================================

def import_subcatchment_parameters(filepath, verbose=True):
    df = pd.read_csv(filepath)

    required_columns = {
        "subcatchment_id",
        "rain_gauge",
        "PAVED_area_ha",
        "PAVED_runoff_coeff",
        "PAVED_slope",
        "PAVED_alpha",
        "PAVED_beta",
        "ROOF_area_ha",
        "ROOF_runoff_coeff",
        "ROOF_slope",
        "ROOF_alpha",
        "ROOF_beta",
        "ADDITIONAL_area_ha",
        "ADDITIONAL_runoff_coeff",
        "ADDITIONAL_slope",
        "ADDITIONAL_alpha",
        "ADDITIONAL_beta",
    }

    missing = required_columns - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in catchment parameter CSV: {sorted(missing)}")

    df = df.dropna(subset=["subcatchment_id", "rain_gauge"]).copy()

    subcatchments = {}
    for _, row in df.iterrows():
        sub_id = str(row["subcatchment_id"]).strip()
        gauge = str(row["rain_gauge"]).strip()

        surfaces = []
        for surface_type in ["PAVED", "ROOF", "ADDITIONAL"]:
            area_col = f"{surface_type}_area_ha"
            if pd.notna(row[area_col]) and float(row[area_col]) > 0:
                surfaces.append({
                    "surface_type": surface_type,
                    "area_ha": float(row[f"{surface_type}_area_ha"]),
                    "runoff_coeff": float(row[f"{surface_type}_runoff_coeff"]),
                    "slope": float(row[f"{surface_type}_slope"]),
                    "alpha": float(row[f"{surface_type}_alpha"]),
                    "beta": float(row[f"{surface_type}_beta"]),
                    "rain_gauge": gauge,
                })

        if not surfaces:
            print(f"WARNING: {sub_id} has no active surfaces")
            continue

        subcatchments[sub_id] = surfaces

    if verbose:
        print(f"Imported {len(subcatchments)} subcatchments")
    return subcatchments


def load_rainfall_dataframe(filepath):
    df = pd.read_csv(filepath)
    if "time_minutes" not in df.columns:
        raise ValueError("Rainfall CSV must contain 'time_minutes' column")
    return df.dropna().sort_values("time_minutes")


def build_rainfall_dict(rainfall_df):
    rainfall_dict = {"time_hr": rainfall_df["time_minutes"].to_numpy() / 60.0}
    for col in rainfall_df.columns:
        if col != "time_minutes":
            rainfall_dict[col] = rainfall_df[col].to_numpy()
    return rainfall_dict


# =========================================================
# INFLOW GENERATION AND ROUTING
# =========================================================

def create_inflow_from_arrays(time_hr, intensity_mm_hr, area_ha, runoff_coeff):
    intensity_m_s = intensity_mm_hr / 3.6e6
    area_m2 = area_ha * 10000.0
    inflow = runoff_coeff * intensity_m_s * area_m2

    return interp1d(
        time_hr,
        inflow,
        kind="linear",
        bounds_error=False,
        fill_value=(0.0, 0.0),
    )


def route_non_linear_reservoir(S0, t_span, times, inflow_func, k, m):
    if inflow_func is None:
        return None

    if len(times) < 2:
        raise ValueError("times must contain at least 2 values")

    last_error = None
    last_message = None

    # Try the most stable solver first. If a particular parameter combination is
    # numerically awkward, fall back to other implicit solvers before finally
    # trying the original explicit RK45 method. This prevents one bad surface
    # combination from crashing the whole sweep.
    for method in SOLVER_METHODS_TO_TRY:
        try:
            sol = solve_ivp(
                fun=non_linear_reservoir_ode,
                t_span=t_span,
                y0=[S0],
                t_eval=times,
                args=(inflow_func, k, m),
                method=method,
                rtol=SOLVER_RTOL,
                atol=SOLVER_ATOL,
                max_step=SOLVER_MAX_STEP,
            )

            if sol.success:
                storage = np.maximum(sol.y[0], 0.0)
                outflow = calculate_outflow(storage, k, m)
                inflow = inflow_func(times)

                return pd.DataFrame({
                    "Time_hr": times,
                    "Inflow_m3s": inflow,
                    "Outflow_m3s": outflow,
                    "Storage_m3": storage,
                })

            last_message = f"{method}: {sol.message}"

        except Exception as e:
            last_error = e
            last_message = f"{method}: {e}"

    if last_error is not None:
        raise RuntimeError(
            "All solver methods failed. Last error was: " + str(last_message)
        ) from last_error

    raise RuntimeError("All solver methods failed. Last message was: " + str(last_message))


# =========================================================
# PARAMETER OVERRIDE HELPERS
# =========================================================

def normalise_surface_name(surface_type):
    return str(surface_type).strip().upper()


def should_override_surface(subcatchment_id, surface_type, surface_params):
    surface_type = normalise_surface_name(surface_type)

    if surface_type not in SURFACES_TO_SWEEP:
        return False

    if surface_type not in surface_params:
        return False

    if TARGET_SUBCATCHMENTS is None:
        return True

    return str(subcatchment_id) in {str(x) for x in TARGET_SUBCATCHMENTS}


def get_alpha_beta_for_surface(subcatchment_id, surface, surface_params):
    surface_type = normalise_surface_name(surface["surface_type"])

    if should_override_surface(subcatchment_id, surface_type, surface_params):
        alpha_use, beta_use = surface_params[surface_type]
    else:
        alpha_use = surface["alpha"]
        beta_use = surface["beta"]

    return float(alpha_use), float(beta_use)


def make_run_name(surface_params):
    parts = []
    for surface_type in ["PAVED", "ROOF", "ADDITIONAL"]:
        if surface_type in surface_params:
            alpha, beta = surface_params[surface_type]
            parts.append(f"{surface_type.lower()}_a{alpha:.3f}_b{beta:.2f}")
    return "__".join(parts)


def build_sweep_parameter_sets():
    """Build the list of parameter sets to run.

    In normal_arrays mode, this keeps the original cartesian/linked/one-at-a-time
    behaviour using the separate alpha and beta arrays.

    In explicit_incomplete_runs mode, each row of INCOMPLETE_RUNS is converted
    directly into one PAVED + ROOF parameter set. This avoids recreating the full
    cartesian product and therefore avoids rerunning completed simulations.
    """

    if SWEEP_MODE == "explicit_incomplete_runs":
        incomplete = np.asarray(INCOMPLETE_RUNS, dtype=float)

        if incomplete.ndim != 2 or incomplete.shape[1] != 4:
            raise ValueError(
                "INCOMPLETE_RUNS must be a 2D array with four columns: "
                "[PAVED_alpha, PAVED_beta, ROOF_alpha, ROOF_beta]"
            )

        parameter_sets = []
        seen_run_names = set()

        for row in incomplete:
            paved_alpha, paved_beta, roof_alpha, roof_beta = row

            surface_params = {
                "PAVED": (round(float(paved_alpha), 6), round(float(paved_beta), 6)),
                "ROOF": (round(float(roof_alpha), 6), round(float(roof_beta), 6)),
            }

            run_name = make_run_name(surface_params)
            if run_name in seen_run_names:
                raise ValueError(f"Duplicate run in INCOMPLETE_RUNS: {run_name}")

            seen_run_names.add(run_name)
            parameter_sets.append(surface_params)

        return parameter_sets

    if SWEEP_MODE != "normal_arrays":
        raise ValueError(
            "SWEEP_MODE must be 'normal_arrays' or 'explicit_incomplete_runs'"
        )

    pair_lists = {
        "PAVED": PAVED_ALPHA_BETA_PAIRS,
        "ROOF": ROOF_ALPHA_BETA_PAIRS,
        "ADDITIONAL": ADDITIONAL_ALPHA_BETA_PAIRS,
    }

    active_surfaces = [normalise_surface_name(s) for s in SURFACES_TO_SWEEP]

    for surface_type in active_surfaces:
        if surface_type not in pair_lists:
            raise ValueError(f"Unknown surface type in SURFACES_TO_SWEEP: {surface_type}")
        if not pair_lists[surface_type]:
            raise ValueError(f"No alpha/beta pairs supplied for {surface_type}")

    if SWEEP_COMBINATION_MODE == "cartesian":
        ordered_lists = [pair_lists[s] for s in active_surfaces]
        parameter_sets = []
        for combo in itertools.product(*ordered_lists):
            parameter_sets.append(dict(zip(active_surfaces, combo)))
        return parameter_sets

    if SWEEP_COMBINATION_MODE == "linked_by_index":
        lengths = [len(pair_lists[s]) for s in active_surfaces]
        if len(set(lengths)) != 1:
            raise ValueError(
                "For linked_by_index mode, all active surface pair lists must have equal length"
            )
        parameter_sets = []
        for i in range(lengths[0]):
            parameter_sets.append({s: pair_lists[s][i] for s in active_surfaces})
        return parameter_sets

    if SWEEP_COMBINATION_MODE == "one_surface_at_a_time":
        parameter_sets = []
        for surface_type in active_surfaces:
            for pair in pair_lists[surface_type]:
                parameter_sets.append({surface_type: pair})
        return parameter_sets

    raise ValueError(
        "SWEEP_COMBINATION_MODE must be 'cartesian', 'linked_by_index', or 'one_surface_at_a_time'"
    )


# =========================================================
# ROUTING FOR ONE PARAMETER SET
# =========================================================

def run_model_for_parameter_set(task):
    run_name = task["run_name"]
    surface_params = task["surface_params"]
    catchment_file = task["catchment_file"]
    rainfall_dict = task["rainfall_dict"]
    times = task["times"]
    t_span = task["t_span"]
    S0 = task["S0"]
    output_folder = Path(task["output_folder"])

    try:
        subcatchments = import_subcatchment_parameters(catchment_file, verbose=False)

        combined_outflows = {}
        detailed_columns = {}
        rainfall_time_hr = rainfall_dict["time_hr"]

        for subcatchment_id, surfaces in subcatchments.items():
            subcatchment_total_outflow = None

            for surface in surfaces:
                surface_type = normalise_surface_name(surface["surface_type"])
                gauge = surface["rain_gauge"]

                if gauge not in rainfall_dict:
                    raise ValueError(
                        f"Gauge '{gauge}' not found in rainfall data for {subcatchment_id} {surface_type}"
                    )

                alpha_use, beta_use = get_alpha_beta_for_surface(
                    subcatchment_id, surface, surface_params
                )

                if alpha_use <= 0:
                    raise ValueError(
                        f"Invalid alpha={alpha_use} for {subcatchment_id} {surface_type}"
                    )

                slope = float(surface["slope"])
                if slope < 0:
                    raise ValueError(
                        f"Negative slope={slope} for {subcatchment_id} {surface_type}"
                    )

                k = (slope ** 0.5) / alpha_use
                m = calculate_exponent(beta_use)

                if m <= 0:
                    raise ValueError(
                        f"Invalid exponent m={m:.6f} for {subcatchment_id} {surface_type}; "
                        f"using m = 5/3 + beta with beta={beta_use}"
                    )

                inflow_func = create_inflow_from_arrays(
                    rainfall_time_hr,
                    rainfall_dict[gauge],
                    float(surface["area_ha"]),
                    float(surface["runoff_coeff"]),
                )

                results = route_non_linear_reservoir(
                    S0=S0,
                    t_span=t_span,
                    times=times,
                    inflow_func=inflow_func,
                    k=k,
                    m=m,
                )

                outflow = results["Outflow_m3s"].to_numpy()

                if subcatchment_total_outflow is None:
                    subcatchment_total_outflow = outflow.copy()
                else:
                    subcatchment_total_outflow += outflow

                if SAVE_SURFACE_DETAIL_COLUMNS:
                    detailed_columns[f"{subcatchment_id}_{surface_type}_Outflow_m3s"] = outflow

            if subcatchment_total_outflow is None:
                raise RuntimeError(f"No active routed surfaces for {subcatchment_id}")

            combined_outflows[subcatchment_id] = subcatchment_total_outflow

        if not combined_outflows:
            raise RuntimeError("No valid outputs produced")

        combined_df = pd.DataFrame({
            sub_id: combined_outflows[sub_id]
            for sub_id in sorted(combined_outflows.keys())
        })
        combined_df.insert(0, "Time_hr", times)

        if SAVE_SURFACE_DETAIL_COLUMNS:
            detail_df = pd.DataFrame(detailed_columns)
            combined_df = pd.concat([combined_df, detail_df], axis=1)

        subcatchment_cols = [c for c in combined_outflows.keys()]
        combined_df["Total_Outflow_m3s"] = combined_df[subcatchment_cols].sum(axis=1)

        out_file = output_folder / f"{run_name}.csv"
        out_file = write_csv_reduced_precision(combined_df, out_file)

        return {
            "run_name": run_name,
            "surface_params": surface_params,
            "status": "success",
            "message": str(out_file),
        }

    except Exception as e:
        return {
            "run_name": run_name,
            "surface_params": surface_params,
            "status": "failed",
            "message": f"{e}\n{traceback.format_exc()}",
        }


# =========================================================
# SINGLE RUN
# =========================================================

def run_single_baseline_model():
    rainfall_df = load_rainfall_dataframe(RAINFALL_CSV_PATH)
    rainfall_dict = build_rainfall_dict(rainfall_df)

    task = {
        "run_name": "baseline_csv_alpha_beta",
        "surface_params": {},
        "catchment_file": str(CATCHMENT_PARAMETERS_FILE),
        "rainfall_dict": rainfall_dict,
        "times": times,
        "t_span": t_span,
        "S0": S0,
        "output_folder": str(RESULTS_FOLDER),
    }

    result = run_model_for_parameter_set(task)
    if result["status"] != "success":
        raise RuntimeError(result["message"])

    print(f"\nSingle baseline simulation complete: {result['message']}")


# =========================================================
# PARALLEL SWEEP DRIVER
# =========================================================

def run_parallel_sweep():
    rainfall_df = load_rainfall_dataframe(RAINFALL_CSV_PATH)
    rainfall_dict = build_rainfall_dict(rainfall_df)

    parameter_sets = build_sweep_parameter_sets()
    total_runs = len(parameter_sets)

    print("=" * 80)
    print("3-SURFACE ALPHA-BETA PARAMETER SWEEP")
    print("=" * 80)
    print(f"Catchment file          : {CATCHMENT_PARAMETERS_FILE}")
    print(f"Rainfall file           : {RAINFALL_CSV_PATH}")
    print(f"Output folder           : {SWEEP_OUTPUT_FOLDER}")
    print(f"Sweep mode              : {SWEEP_MODE}")
    print(f"Combination mode        : {SWEEP_COMBINATION_MODE}")
    print(f"Surfaces swept          : {SURFACES_TO_SWEEP}")
    print("Exponent convention     : m = 5/3 + beta")
    print(f"Initial storage S0      : {S0}")
    print(f"Output decimals         : {OUTPUT_DECIMALS}")
    print(f"Compressed CSV outputs  : {COMPRESS_CSV_OUTPUTS}")
    print(f"Workers                 : {MAX_WORKERS}")
    print(f"Total runs              : {total_runs}")
    print("=" * 80)

    tasks = []
    for surface_params in parameter_sets:
        run_name = make_run_name(surface_params)
        tasks.append({
            "run_name": run_name,
            "surface_params": surface_params,
            "catchment_file": str(CATCHMENT_PARAMETERS_FILE),
            "rainfall_dict": rainfall_dict,
            "times": times,
            "t_span": t_span,
            "S0": S0,
            "output_folder": str(SWEEP_OUTPUT_FOLDER),
        })

    start_time = time.perf_counter()
    results_summary = []

    n_complete = 0
    n_success = 0
    n_failed = 0

    # Submit only a small number of jobs at once. This is much more stable on
    # Windows than submitting the full sweep to the process pool in one go,
    # because each process receives rainfall arrays, solves many reservoirs,
    # creates long output arrays, and writes a CSV file.
    batch_size = max(1, MAX_WORKERS * 2)

    def append_result(result):
        nonlocal n_success, n_failed

        flat_result = {
            "run_name": result["run_name"],
            "status": result["status"],
            "message": result["message"],
        }

        for surface_type in ["PAVED", "ROOF", "ADDITIONAL"]:
            pair = result.get("surface_params", {}).get(surface_type)
            if pair is not None:
                flat_result[f"{surface_type}_alpha"] = pair[0]
                flat_result[f"{surface_type}_beta"] = pair[1]

        results_summary.append(flat_result)

        if result["status"] == "success":
            n_success += 1
        else:
            n_failed += 1

    try:
        for batch_start in range(0, total_runs, batch_size):
            batch_tasks = tasks[batch_start:batch_start + batch_size]
            batch_number = (batch_start // batch_size) + 1
            total_batches = int(np.ceil(total_runs / batch_size))

            print(
                f"\n[batch {batch_number}/{total_batches}] "
                f"submitting {len(batch_tasks)} jobs"
            )

            with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
                futures = {
                    executor.submit(run_model_for_parameter_set, task): task
                    for task in batch_tasks
                }

                for future in as_completed(futures):
                    task = futures[future]
                    n_complete += 1

                    try:
                        result = future.result()
                    except BrokenProcessPool as e:
                        result = {
                            "run_name": task["run_name"],
                            "surface_params": task.get("surface_params", {}),
                            "status": "failed",
                            "message": (
                                "Broken process pool. A worker process was terminated "
                                f"abruptly while this run was active or pending: {e}"
                            ),
                        }
                    except Exception as e:
                        result = {
                            "run_name": task["run_name"],
                            "surface_params": task.get("surface_params", {}),
                            "status": "failed",
                            "message": f"{e}\n{traceback.format_exc()}",
                        }

                    append_result(result)

                    elapsed = time.perf_counter() - start_time
                    avg_time = elapsed / n_complete if n_complete else 0.0
                    remaining = total_runs - n_complete
                    eta_seconds = avg_time * remaining

                    if (n_complete % REPORT_EVERY == 0) or (n_complete == total_runs):
                        print(
                            f"[progress {n_complete:>4}/{total_runs}] "
                            f"success={n_success} failed={n_failed} | "
                            f"elapsed={elapsed/60:.2f} min | "
                            f"eta={eta_seconds/60:.2f} min | "
                            f"last={result['run_name']} ({result['status']})"
                        )

                        if result["status"] == "failed":
                            print(f"    reason: {result['message']}")

    finally:
        total_elapsed = time.perf_counter() - start_time

        summary_df = pd.DataFrame(results_summary)
        if not summary_df.empty:
            summary_df = summary_df.sort_values(["status", "run_name"]).reset_index(drop=True)

        summary_file = SWEEP_OUTPUT_FOLDER / "3surface_alpha_beta_sweep_summary.csv"
        summary_file = write_csv_reduced_precision(summary_df, summary_file)

        print("\n" + "=" * 80)
        print("3-SURFACE PARAMETER SWEEP FINISHED")
        print("=" * 80)
        print(f"Successful runs : {n_success}")
        print(f"Failed runs     : {n_failed}")
        print(f"Summary file    : {summary_file}")
        print(f"Total elapsed   : {total_elapsed/60:.2f} minutes")
        print("=" * 80)


# =========================================================
# MAIN
# =========================================================

if __name__ == "__main__":
    if MODE == "single_run":
        run_single_baseline_model()
    elif MODE == "parallel_sweep":
        run_parallel_sweep()
    else:
        raise ValueError("MODE must be 'single_run' or 'parallel_sweep'")
