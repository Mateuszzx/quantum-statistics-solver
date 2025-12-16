import uuid
import pandas as pd
from pathlib import Path
from tqdm.auto import tqdm
import sys
import os


project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../"))
if project_root not in sys.path:
    sys.path.append(project_root)

from qss.solver_options import sa_defaults, CoolingSchedule
from qss.chemical_potential import chemical_potential


def run_sa_experiment(
    n_vals,
    beta_vals,
    H_a_vals,
    distros,
    e_k,
    N=100,
    solver_opts=None,
    logs_dir="logs",
    output_file="chemical_potential_experiment_results.csv",
    x0=0.0,
    a0=-10.0,
    b0=10.0
):
    """
    Runs the Simulated Annealing experiment for chemical potential.

    Parameters:
    -----------
    n_vals : array-like
        Values of band filling (n).
    beta_vals : array-like
        Values of inverse temperature (beta).
    H_a_vals : array-like
        Values of magnetic field (H_a).
    distros : list
        List of DistType distributions to test.
    e_k : array-like
        Energy dispersion matrix.
    N : int
        System size (linear dimension).
    solver_opts : CombinedSolverOptions, optional
        Solver configuration. If None, default "patient" SA settings are used.
    logs_dir : str or Path
        Directory to save individual simulation logs.
    output_file : str
        Filename for the aggregated results CSV.
    x0 : float
        Initial guess for chemical potential.
    a0 : float
        Lower bound for search range.
    b0 : float
        Upper bound for search range.

    Returns:
    --------
    pd.DataFrame
        Aggregated results.
    """

    # --- Default Solver Configuration if not provided ---
    if solver_opts is None:
        print("Using default patient SA settings...")
        py_opts = sa_defaults(
            T0=1e3,
            T_min=1e-9,
            alpha=0.999,       # Slow cooling
            max_iters=50000,   # High total iters
            min_iters=2000,    # Force exploration
            stall_iters=1000,  # High stall tolerance
            step_decay=0.999,  # Slow step decay
            step_sigma=5.0,    # High initial sigma
            bounds=[(-10., 10.)],
            scale_sigma_by_T=True,
            reflect_bounds=True,
            sched=CoolingSchedule.Exponential
        )
        native_opts = py_opts.to_qss()
    else:
        native_opts = solver_opts

    logs_path = Path(logs_dir)
    logs_path.mkdir(parents=True, exist_ok=True)
    print(f"Logs will be saved to: {logs_path.absolute()}")

    results = []

    # Calculate total iterations for progress bar
    total_sims = len(distros) * len(beta_vals) * len(H_a_vals) * len(n_vals)

    with tqdm(total=total_sims, desc="Simulations") as pbar:
        for d in distros:
            for b in beta_vals:
                for h in H_a_vals:
                    for n_val in n_vals:
                        try:
                            # Generate unique log filename
                            rnd_suffix = uuid.uuid4().hex[:6]
                            dist_name = str(d).split('.')[-1] # e.g. "FD"
                            log_filename = logs_path / f"log_sa_{dist_name}_b{b:.2f}_n{n_val:.2f}_H{h:.2f}_{rnd_suffix}.csv"

                            # Set log file in options
                            if hasattr(native_opts, 'log_file'):
                                native_opts.log_file = str(log_filename)
                            elif hasattr(native_opts, 'sa'):
                                native_opts.sa.log_file = str(log_filename)

                            # Run calculation
                            res_df = chemical_potential(
                                distros=[d],
                                n_vals=[n_val],
                                beta_vals=[b],
                                H_a_vals=[h],
                                N=N,
                                E=e_k,
                                solver_opts=native_opts,
                                x0=x0,
                                a0=a0,
                                b0=b0,
                            )

                            results.append(res_df)

                        except Exception as e:
                            print(f"\nError for: {d}, n={n_val}, beta={b}, H={h}")
                            print(f"Details: {e}")

                        pbar.update(1)

    # Aggregate and Save
    if results:
        final_df = pd.concat(results)
        final_df.to_csv(output_file)
        print(f"\nExperiment finished. Aggregated results saved to: {output_file}")
        return final_df
    else:
        print("\nNo results collected.")
        return pd.DataFrame()
