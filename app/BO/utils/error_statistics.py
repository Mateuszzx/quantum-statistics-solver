import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm.auto import tqdm
from pathlib import Path
from qss.solver_options import sa_defaults
from .sa_helpers import run_single_experiment
from .config_loader import load_config

def run_error_statistics_experiment():
    config = load_config("error_stats")

    # Grid of n and T
    n_values = np.linspace(
        config["ranges"]["n_start"],
        config["ranges"]["n_end"],
        config["ranges"]["n_steps"]
    )
    beta_values = np.logspace(
        config["ranges"]["beta_start_exp"],
        config["ranges"]["beta_end_exp"],
        config["ranges"]["beta_steps"]
    )

    error_stats = []

    # Base options
    base_opts = sa_defaults(**config["solver"])

    print("Running Error Statistics Experiment...")
    for n in tqdm(n_values, desc="n loop"):
        for beta in tqdm(beta_values, desc="beta loop", leave=False):
            res = run_single_experiment(
                n, beta, base_opts,
                H_a=config["physics"]["H_a"],
                dist_type=config["physics"]["dist_type"],
                log_prefix="error_stat"
            )
            error_stats.append(res)

    return pd.DataFrame(error_stats)

def plot_error_statistics(df_error):
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # 1. Mean Error vs Temperature (averaged over n)
    mean_err_T = df_error.groupby("T")["final_error"].mean().reset_index()
    sns.lineplot(data=mean_err_T, x="T", y="final_error", marker="o", ax=axes[0])
    axes[0].set_xscale("log")
    axes[0].set_title("Mean Error vs Temperature (Averaged over n)")
    axes[0].set_ylabel("Mean Absolute Error")

    # 2. Mean Error vs Band Filling (averaged over T)
    mean_err_n = df_error.groupby("n")["final_error"].mean().reset_index()
    sns.lineplot(data=mean_err_n, x="n", y="final_error", marker="o", ax=axes[1])
    axes[1].set_title("Mean Error vs Band Filling (Averaged over T)")
    axes[1].set_ylabel("Mean Absolute Error")

    plt.tight_layout()
    plt.show()
