import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm.auto import tqdm
from pathlib import Path
from qss.solver_options import sa_defaults
from .sa_helpers import run_single_experiment
from .config_loader import load_config

def run_step_sigma_experiment():
    config = load_config("step_sigma")

    # Fix n and T
    n_fixed = config["physics"]["n_fixed"]
    beta_fixed = config["physics"]["beta_fixed"]

    step_sigmas = np.linspace(
        config["ranges"]["sigma_start"],
        config["ranges"]["sigma_end"],
        config["ranges"]["sigma_steps"]
    )
    iter_stats = []

    print("Running Step Sigma Experiment...")
    for sigma in tqdm(step_sigmas, desc="sigma loop"):
        # Modify options
        opts_dict = config["solver"].copy()
        opts_dict["step_sigma"] = float(sigma)
        opts = sa_defaults(**opts_dict)

        res = run_single_experiment(
            n_fixed, beta_fixed, opts,
            H_a=config["physics"]["H_a"],
            dist_type=config["physics"]["dist_type"],
            log_prefix="sigma_stat"
        )
        res["step_sigma"] = sigma
        iter_stats.append(res)

    df_iter = pd.DataFrame(iter_stats)

    # Save summary to CSV
    output_csv = "data/sa_step_sigma_stats.csv"
    df_iter.drop(columns=["log_df"]).to_csv(output_csv, index=False)
    print(f"Saved step sigma statistics to {Path(output_csv).resolve()}")

    return df_iter

def plot_step_sigma_statistics(df_iter):
    n_fixed = 0.8
    beta_fixed = 5.0

    plt.figure(figsize=(8, 6))
    sns.lineplot(data=df_iter, x="step_sigma", y="best_iter", marker="o")
    plt.title(f"Iteration of Best Solution vs Step Sigma (n={n_fixed}, T={1/beta_fixed})")
    plt.xlabel("Step Sigma")
    plt.ylabel("Iteration Found")
    plt.show()
