import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm.auto import tqdm
from pathlib import Path
from qss.solver_options import sa_defaults
from .sa_helpers import run_single_experiment
from .config_loader import load_config

def run_convergence_speed_experiment():
    config = load_config("convergence_speed")
    output_csv = Path("../data") / "convergence_speed_stats.csv"
    n_fixed = config["physics"]["n_fixed"]
    beta_fixed = config["physics"]["beta_fixed"]

    alphas = np.round(np.linspace(
        config["ranges"]["alpha_start"],
        config["ranges"]["alpha_end"],
        config["ranges"]["alpha_steps"]
    ), 3)

    T0s = np.round(np.logspace(
        config["ranges"]["T0_start_exp"],
        config["ranges"]["T0_end_exp"],
        config["ranges"]["T0_steps"]
    ))

    speed_stats = []

    print("Running Convergence Speed Experiment...")
    for alpha in tqdm(alphas, desc="alpha loop"):
        for T0 in tqdm(T0s, desc="T0 loop", leave=False):
            opts_dict = config["solver"].copy()
            opts_dict["T0"] = float(T0)
            opts_dict["alpha"] = float(alpha)
            opts_dict["step_decay"] = float(alpha)

            opts = sa_defaults(**opts_dict)

            res = run_single_experiment(
                n_fixed, beta_fixed, opts,
                H_a=config["physics"]["H_a"],
                dist_type=config["physics"]["dist_type"],
                log_prefix="speed_stat"
            )

            # Add parameters to result
            res["alpha"] = alpha
            res["T0"] = T0
            speed_stats.append(res)

    df_speed = pd.DataFrame(speed_stats)

    # Ensure output directory exists
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    df_speed.drop(columns=["log_df"]).to_csv(output_csv, index=False)
    print(f"Saved convergence speed statistics to {Path(output_csv).resolve()}")

    return df_speed

def plot_convergence_speed(df_speed: pd.DataFrame) -> None:
    n_fixed = 0.8
    beta_fixed = 5.0

    # Plot Heatmap
    pivot_table = df_speed.pivot(index="alpha", columns="T0", values="conv_iter")

    plt.figure(figsize=(10, 8))
    sns.heatmap(pivot_table, annot=True, fmt=".0f", cmap="viridis_r")
    plt.title(f"Convergence Speed (Iterations to < 1e-6)\nn={n_fixed}, T={1/beta_fixed}")
    plt.xlabel("Start Temperature (T0)")
    plt.ylabel("Cooling Rate (alpha)")
    plt.show()
