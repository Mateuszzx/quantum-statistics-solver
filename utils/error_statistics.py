import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from qss.solver_options import sa_defaults

from .config_loader import load_config
from .sa_helpers import run_single_experiment


def run_error_statistics_experiment() -> pd.DataFrame:
    config = load_config("error_stats")

    # Grid of n and T
    n_values = np.linspace(
        config["ranges"]["n_start"], config["ranges"]["n_end"], config["ranges"]["n_steps"]
    )
    beta_values = np.logspace(
        config["ranges"]["beta_start_exp"],
        config["ranges"]["beta_end_exp"],
        config["ranges"]["beta_steps"],
    )

    error_stats = []

    # Base options
    base_opts = sa_defaults(**config["solver"])

    print("Running Error Statistics Experiment...")
    for n in tqdm(n_values, desc="n loop"):
        for beta in tqdm(beta_values, desc="beta loop", leave=False):
            res = run_single_experiment(
                n,
                beta,
                base_opts,
                h_a=config["physics"]["H_a"],
                dist_type=config["physics"]["dist_type"],
                log_prefix="error_stat",
            )
            error_stats.append(res)

    return pd.DataFrame(error_stats)


def plot_error_statistics(df_error: pd.DataFrame) -> None:
    _, axes = plt.subplots(1, 2, figsize=(16, 6))

    # 1. Mean Error vs Temperature (averaged over n)
    stats_t = df_error.groupby("T")["final_error"].agg(["mean", "std"]).reset_index()
    axes[0].errorbar(
        stats_t["T"],
        stats_t["mean"],
        yerr=stats_t["std"],
        fmt="o",
        capsize=5,
        ecolor="gray",
        color="blue",
    )
    axes[0].set_xscale("log")
    axes[0].set_title("Mean Error vs Temperature (Averaged over n)")
    axes[0].set_xlabel("Temperature (T)")
    axes[0].set_ylabel("Mean Absolute Error")
    axes[0].grid(True, which="both", ls="--", alpha=0.5)

    # 2. Mean Error vs Band Filling (averaged over T)
    stats_n = df_error.groupby("n")["final_error"].agg(["mean", "std"]).reset_index()
    axes[1].errorbar(
        stats_n["n"],
        stats_n["mean"],
        yerr=stats_n["std"],
        fmt="o",
        capsize=5,
        ecolor="gray",
        color="green",
    )
    axes[1].set_title("Mean Error vs Band Filling (Averaged over T)")
    axes[1].set_xlabel("Band Filling (n)")
    axes[1].set_ylabel("Mean Absolute Error")
    axes[1].grid(True, which="both", ls="--", alpha=0.5)

    plt.tight_layout()
    plt.show()
