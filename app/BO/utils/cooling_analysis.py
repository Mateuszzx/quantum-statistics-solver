import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm.auto import tqdm

from qss.solver_options import sa_defaults, CoolingSchedule
from .sa_helpers import run_single_experiment
from .config_loader import load_config

sns.set_theme(style="whitegrid")

def run_cooling_analysis_experiment():
    """
    Runs the cooling analysis experiment based on the 'cooling.yaml' config.
    Returns a list of result dictionaries (from run_single_experiment).
    """
    config = load_config("cooling")

    # Physics parameters
    n = config["physics"]["n"]
    beta = config["physics"]["beta"]
    H_a = config["physics"]["H_a"]
    dist_type = config["physics"]["dist_type"]

    # Experiment parameters
    schedules_str = config["experiment"]["schedules"]
    start_temps = config["experiment"]["start_temps"]

    # Convert schedule strings to Enums
    schedules = [getattr(CoolingSchedule, s) for s in schedules_str]

    results = []

    print("Starting Cooling Analysis Simulations...")

    # Iterate over schedules and start temperatures
    for sched in tqdm(schedules, desc="Schedules"):
        for T0 in tqdm(start_temps, desc="Start Temps", leave=False):

            # Create specific options for this run
            current_solver_config = config["solver"].copy()
            current_solver_config["sched"] = sched
            current_solver_config["T0"] = T0

            opts = sa_defaults(**current_solver_config)

            # Define a label for the plot
            label = f"{sched.name}, T0={T0}"

            # Run experiment
            res = run_single_experiment(
                n, beta, opts,
                H_a=H_a,
                dist_type=dist_type,
                log_prefix=f"sa_traj_{sched.name}_{T0}"
            )

            if not res["log_df"].empty:
                res["log_df"]["label"] = label
                res["schedule"] = sched.name
                res["T0"] = T0
                results.append(res)

    print("Simulations complete.")
    return results


def plot_cooling_analysis(results):
    """
    Plots the optimization trajectories from the experiment results.
    """
    if not results:
        print("No results to plot.")
        return

    schedules_found = sorted(list(set(r["schedule"] for r in results)))
    t0_found = sorted(list(set(r["T0"] for r in results)))

    n_rows = len(t0_found)
    n_cols = len(schedules_found)

    if n_rows == 0 or n_cols == 0:
        print("No valid results found.")
        return

    fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(4 * n_cols, 3 * n_rows), sharex=True, sharey=True)
    fig.suptitle("Simulated Annealing Trajectory: Cooling Schedule & Start Temp")

    # Handle case where n_rows or n_cols is 1 (axes is not 2D array)
    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)

    for i, t0 in enumerate(t0_found):
        for j, sched_name in enumerate(schedules_found):
            ax = axes[i, j]

            # Find the matching result
            matching_dfs = [r["log_df"] for r in results if r["schedule"] == sched_name and r["T0"] == t0]

            for df in matching_dfs:
                # Plot Energy
                ax.plot(df['iteration'], df['energy'] + 1e-16, label="Energy", alpha=0.7)

                # Plot Temperature on twin axis
                ax2 = ax.twinx()
                ax2.plot(df['iteration'], df['temperature'], color='red', linestyle='--', label='SA Temp', alpha=0.5)

                ax2.tick_params(axis='y', labelcolor='red')

                # Only show right y-label on the last column
                if j == n_cols - 1:
                    ax2.set_ylabel("SA Temperature")
                else:
                    ax2.set_yticklabels([])

            # Labels and Titles
            if i == 0:
                ax.set_title(sched_name)
            if j == 0:
                ax.set_ylabel(f"T0={t0}\nCost (Energy)")
            if i == n_rows - 1:
                ax.set_xlabel("Iteration")


    plt.tight_layout()
    plt.show()

