import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import matplotlib.pyplot as plt
import os
import time

# =========================================================
# USER SETTINGS
# =========================================================

MODE = "single_run"   # "single_run" or "parallel_sweep"

CATCHMENT_PARAMETERS_FILE = Path("./Subcatchments.csv")
RAINFALL_CSV_PATH = Path("./Raingauges.csv")

RESULTS_FOLDER = Path("./routing_results_TSR")
PLOT_FOLDER = Path("./routing_plots_TSR")
SWEEP_OUTPUT_FOLDER = Path("./alpha_beta_sweep_outputs")

RESULTS_FOLDER.mkdir(exist_ok=True)
PLOT_FOLDER.mkdir(exist_ok=True)
SWEEP_OUTPUT_FOLDER.mkdir(exist_ok=True)

# Time settings
S0 = 0.0
START_TIME = 0.0
END_TIME = 1488.0
TIME_STEP = 1.0 / 30.0   # hours

times = np.arange(START_TIME, END_TIME + TIME_STEP / 2, TIME_STEP)
t_span = (START_TIME, END_TIME)

# Sweep settings
alpha_values = np.concatenate(([0.025, 0.03, 0.04, 0.05], np.arange(0.002, 0.022, 0.002)))
beta_values = np.arange(0.1, 1.0, 0.1)

# Set to None to apply each alpha/beta pair to all subcatchments
target_subcatchments = None

# If target_subcatchments is not None, preserve original alpha/beta
# values for non-targeted subcatchments
preserve_non_target_values = True

MAX_WORKERS = max(1, (os.cpu_count() or 2) - 1)
REPORT_EVERY = 1

SAVE_INDIVIDUAL_RESULTS = False
MAKE_PLOTS = False


# =========================================================
# CORE MODEL PHYSICS
# =========================================================

def calculate_outflow(S, k, m):
    S = np.asarray(S)
    S_eff = np.maximum(S, 0.0)
    Q = k * np.power(S_eff, m)
    return Q.item() if Q.shape == () else Q


def non_linear_reservoir_ode(t, S, I_func, k, m):
    Q = calculate_outflow(S, k, m)
    I = float(I_func(t))
    dSdt = np.array([I]) - Q
    return dSdt


# =========================================================
# DATA IMPORT
# =========================================================

def import_subcatchment_parameters(filepath):
    try:
        df = pd.read_csv(filepath)

        required_columns = {
            "subcatchment_id",
            "area_ha",
            "runoff_coeff",
            "slope",
            "alpha",
            "beta",
            "rain_gauge",
        }

        if not required_columns.issubset(df.columns):
            raise ValueError(f"CSV must contain {required_columns}")

        df = df.dropna(subset=required_columns).sort_values("subcatchment_id")

        subcatchments = {}
        for _, row in df.iterrows():
            sub_id = row["subcatchment_id"]
            subcatchments[sub_id] = {
                "area_ha": float(row["area_ha"]),
                "runoff_coeff": float(row["runoff_coeff"]),
                "slope": float(row["slope"]),
                "alpha": float(row["alpha"]),
                "beta": float(row["beta"]),
                "rain_gauge": row["rain_gauge"],
            }

        print(f"Imported {len(subcatchments)} subcatchments")
        return subcatchments

    except FileNotFoundError:
        print("Catchment parameter file not found")
        return {}
    except Exception as e:
        print(f"Error importing subcatchment parameters: {e}")
        return {}


def load_rainfall_dataframe(filepath):
    try:
        df = pd.read_csv(filepath)

        if "time_minutes" not in df.columns:
            raise ValueError("Rainfall CSV must contain 'time_minutes' column")

        df = df.dropna().sort_values("time_minutes")
        return df

    except Exception as e:
        print(f"Error loading rainfall file: {e}")
        return None


def create_inflow_from_rainfall_df(df, gauge_column, area_ha, runoff_coeff):
    try:
        if gauge_column not in df.columns:
            raise ValueError(f"Gauge '{gauge_column}' not found in rainfall file")

        intensity_mm_hr = df[gauge_column].to_numpy()

        # Convert mm/hr -> m/s
        intensity_m_s = intensity_mm_hr / 3.6e6

        area_m2 = area_ha * 10000
        inflow = runoff_coeff * intensity_m_s * area_m2

        time_hr = df["time_minutes"].to_numpy() / 60

        inflow_func = interp1d(
            time_hr,
            inflow,
            kind="linear",
            bounds_error=False,
            fill_value=(0.0, 0.0),
        )

        return inflow_func

    except Exception as e:
        print(f"Error creating inflow for gauge {gauge_column}: {e}")
        return None


# =========================================================
# ROUTING ENGINE
# =========================================================

def route_non_linear_reservoir(S0, t_span, times, inflow_func, k, m):
    if inflow_func is None:
        return None

    sol = solve_ivp(
        fun=non_linear_reservoir_ode,
        t_span=t_span,
        y0=[S0],
        t_eval=times,
        args=(inflow_func, k, m),
        method="RK45",
        rtol=1e-6,
        atol=1e-9,
        max_step=(times[1] - times[0]),
    )

    Storage = np.maximum(sol.y[0], 0.0)
    Outflow = calculate_outflow(Storage, k, m)
    Inflow = inflow_func(times)

    results = pd.DataFrame(
        {
            "Time_hr": times,
            "Inflow_m3s": Inflow,
            "Outflow_m3s": Outflow,
            "Storage_m3": Storage,
        }
    )

    return results


# =========================================================
# SINGLE BASELINE RUN
# =========================================================

def run_single_baseline_model():
    subcatchments = import_subcatchment_parameters(CATCHMENT_PARAMETERS_FILE)
    rainfall_df = load_rainfall_dataframe(RAINFALL_CSV_PATH)

    if rainfall_df is None:
        raise RuntimeError("Rainfall data could not be loaded")

    if not subcatchments:
        raise RuntimeError("No subcatchments loaded")

    combined_outflows = {}
    time_vector = None

    for subcatchment_id, subcatchment in subcatchments.items():
        print(f"\nRunning subcatchment {subcatchment_id}")

        area = subcatchment["area_ha"]
        runoff_coeff = subcatchment["runoff_coeff"]
        slope = subcatchment["slope"]
        alpha = subcatchment["alpha"]
        beta = subcatchment["beta"]
        gauge = subcatchment["rain_gauge"]

        k = (slope ** 0.5) / alpha
        m = (5 / 3) - beta

        inflow_func = create_inflow_from_rainfall_df(
            rainfall_df,
            gauge,
            area,
            runoff_coeff
        )

        results = route_non_linear_reservoir(
            S0=S0,
            t_span=t_span,
            times=times,
            inflow_func=inflow_func,
            k=k,
            m=m
        )

        if results is None:
            print(f"Skipping {subcatchment_id} due to error")
            continue

        if time_vector is None:
            time_vector = results["Time_hr"].values

        combined_outflows[subcatchment_id] = results["Outflow_m3s"].values

        if SAVE_INDIVIDUAL_RESULTS:
            results_file = RESULTS_FOLDER / f"routing_results_{subcatchment_id}.csv"
            results.to_csv(results_file, index=False)

        if MAKE_PLOTS:
            plt.figure(figsize=(10, 6))
            plt.plot(results["Time_hr"], results["Inflow_m3s"], label="Inflow")
            plt.plot(results["Time_hr"], results["Outflow_m3s"], label="Outflow")
            plt.xlabel("Time (hr)")
            plt.ylabel("Flow (m3/s)")
            plt.title(f"Subcatchment {subcatchment_id} (Gauge: {gauge})")
            plt.legend()
            plot_file = PLOT_FOLDER / f"routing_plot_{subcatchment_id}.png"
            plt.savefig(plot_file, dpi=300, bbox_inches="tight")
            plt.close()

    if len(combined_outflows) == 0:
        raise RuntimeError("No valid outputs produced")

    combined_df = pd.DataFrame(dict(sorted(combined_outflows.items())))
    combined_df.insert(0, "Time_hr", time_vector)
    combined_df["Total_Outflow_m3s"] = combined_df.drop(columns=["Time_hr"]).sum(axis=1)

    combined_file = RESULTS_FOLDER / "combined_outflow_hydrographs.csv"
    combined_df.to_csv(combined_file, index=False)

    print(f"\nCombined outflow file saved to {combined_file}")
    print("\nSingle baseline simulation complete.")


# =========================================================
# PARALLEL SWEEP WORKER
# =========================================================

def run_single_alpha_beta(
    alpha,
    beta,
    catchment_file,
    rainfall_file,
    output_folder,
    S0,
    t_span,
    times,
    target_subcatchments=None,
    preserve_non_target_values=True
):
    run_name = f"alpha_{alpha:.3f}_beta_{beta:.2f}"
    out_file = Path(output_folder) / f"{run_name}.csv"

    try:
        subcatchments = import_subcatchment_parameters(catchment_file)
        rainfall_df = load_rainfall_dataframe(rainfall_file)

        if rainfall_df is None:
            return {
                "run_name": run_name,
                "alpha": alpha,
                "beta": beta,
                "status": "failed",
                "message": "Rainfall data could not be loaded"
            }

        if not subcatchments:
            return {
                "run_name": run_name,
                "alpha": alpha,
                "beta": beta,
                "status": "failed",
                "message": "No subcatchments loaded"
            }

        combined_outflows = {}
        time_vector = None

        for subcatchment_id, subcatchment in subcatchments.items():
            area = subcatchment["area_ha"]
            runoff_coeff = subcatchment["runoff_coeff"]
            slope = subcatchment["slope"]
            gauge = subcatchment["rain_gauge"]

            if target_subcatchments is None:
                alpha_use = alpha
                beta_use = beta
            else:
                if subcatchment_id in target_subcatchments:
                    alpha_use = alpha
                    beta_use = beta
                else:
                    if preserve_non_target_values:
                        alpha_use = subcatchment["alpha"]
                        beta_use = subcatchment["beta"]
                    else:
                        alpha_use = alpha
                        beta_use = beta

            if alpha_use <= 0:
                return {
                    "run_name": run_name,
                    "alpha": alpha,
                    "beta": beta,
                    "status": "failed",
                    "message": f"Invalid alpha={alpha_use} for subcatchment {subcatchment_id}"
                }

            if slope < 0:
                return {
                    "run_name": run_name,
                    "alpha": alpha,
                    "beta": beta,
                    "status": "failed",
                    "message": f"Negative slope={slope} for subcatchment {subcatchment_id}"
                }

            k = (slope ** 0.5) / alpha_use
            m = (5.0 / 3.0) - beta_use

            inflow_func = create_inflow_from_rainfall_df(
                df=rainfall_df,
                gauge_column=gauge,
                area_ha=area,
                runoff_coeff=runoff_coeff
            )

            if inflow_func is None:
                return {
                    "run_name": run_name,
                    "alpha": alpha,
                    "beta": beta,
                    "status": "failed",
                    "message": f"Could not create inflow for gauge '{gauge}' in subcatchment {subcatchment_id}"
                }

            results = route_non_linear_reservoir(
                S0=S0,
                t_span=t_span,
                times=times,
                inflow_func=inflow_func,
                k=k,
                m=m
            )

            if results is None:
                return {
                    "run_name": run_name,
                    "alpha": alpha,
                    "beta": beta,
                    "status": "failed",
                    "message": f"Routing failed for subcatchment {subcatchment_id}"
                }

            if time_vector is None:
                time_vector = results["Time_hr"].values

            combined_outflows[subcatchment_id] = results["Outflow_m3s"].values

        if len(combined_outflows) == 0:
            return {
                "run_name": run_name,
                "alpha": alpha,
                "beta": beta,
                "status": "failed",
                "message": "No valid outputs produced"
            }

        combined_df = pd.DataFrame(dict(sorted(combined_outflows.items())))
        combined_df.insert(0, "Time_hr", time_vector)

        sub_cols = [c for c in combined_df.columns if c != "Time_hr"]
        combined_df["Total_Outflow_m3s"] = combined_df[sub_cols].sum(axis=1)

        combined_df.to_csv(out_file, index=False)

        return {
            "run_name": run_name,
            "alpha": alpha,
            "beta": beta,
            "status": "success",
            "message": str(out_file)
        }

    except Exception as e:
        return {
            "run_name": run_name,
            "alpha": alpha,
            "beta": beta,
            "status": "failed",
            "message": str(e)
        }


# =========================================================
# PARALLEL SWEEP DRIVER
# =========================================================

def run_parallel_sweep():
    alpha_beta_pairs = [(a, b) for a in alpha_values for b in beta_values]
    total_runs = len(alpha_beta_pairs)

    print("=" * 70)
    print("ALPHA-BETA PARALLEL SWEEP")
    print("=" * 70)
    print(f"Catchment file : {CATCHMENT_PARAMETERS_FILE}")
    print(f"Rainfall file  : {RAINFALL_CSV_PATH}")
    print(f"Output folder  : {SWEEP_OUTPUT_FOLDER}")
    print(f"Workers        : {MAX_WORKERS}")
    print(f"Total runs     : {total_runs}")
    print(f"Report every   : {REPORT_EVERY} completions")
    print("=" * 70)

    start_time = time.perf_counter()
    results_summary = []

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {}

        for idx, (alpha, beta) in enumerate(alpha_beta_pairs, start=1):
            run_name = f"alpha_{alpha:.3f}_beta_{beta:.2f}"

            future = executor.submit(
                run_single_alpha_beta,
                alpha,
                beta,
                str(CATCHMENT_PARAMETERS_FILE),
                str(RAINFALL_CSV_PATH),
                str(SWEEP_OUTPUT_FOLDER),
                S0,
                t_span,
                times,
                target_subcatchments,
                preserve_non_target_values
            )

            futures[future] = (alpha, beta, run_name)
            print(f"[submitted {idx:>4}/{total_runs}] {run_name}")

        print("\nAll jobs submitted. Processing...\n")

        n_complete = 0
        n_success = 0
        n_failed = 0

        for future in as_completed(futures):
            alpha, beta, run_name = futures[future]
            n_complete += 1

            try:
                result = future.result()
            except Exception as e:
                result = {
                    "run_name": run_name,
                    "alpha": alpha,
                    "beta": beta,
                    "status": "failed",
                    "message": str(e)
                }

            results_summary.append(result)

            if result["status"] == "success":
                n_success += 1
            else:
                n_failed += 1

            elapsed = time.perf_counter() - start_time
            avg_time_per_completed = elapsed / n_complete if n_complete > 0 else 0.0
            remaining = total_runs - n_complete
            eta_seconds = avg_time_per_completed * remaining

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

    total_elapsed = time.perf_counter() - start_time

    summary_df = pd.DataFrame(results_summary)
    summary_df = summary_df.sort_values(["status", "alpha", "beta"]).reset_index(drop=True)

    summary_file = SWEEP_OUTPUT_FOLDER / "alpha_beta_sweep_summary.csv"
    summary_df.to_csv(summary_file, index=False)

    print("\n" + "=" * 70)
    print("PARALLEL SWEEP COMPLETE")
    print("=" * 70)
    print(f"Successful runs : {n_success}")
    print(f"Failed runs     : {n_failed}")
    print(f"Summary file    : {summary_file}")
    print(f"Total elapsed   : {total_elapsed/60:.2f} minutes")
    print("=" * 70)


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