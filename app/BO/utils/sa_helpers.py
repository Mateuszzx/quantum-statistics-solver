import numpy as np
import pandas as pd
from pathlib import Path
import uuid
import copy
from qss.chemical_potential import chemical_potential
from qss.distributions import DistType, prepare_2D_lattice_data
from qss.solver_options import CombinedSolverOptions

N = 100
t = -1
t_prim = 0.25
kx = np.linspace(-np.pi, np.pi, N, endpoint=False)
ky = np.linspace(-np.pi, np.pi, N, endpoint=False)
e_k, _, _ = prepare_2D_lattice_data(kx, ky, t, t_prim)

H_a = 0.0
dist_type = DistType.FD

def run_single_experiment(n_val,
                          beta_val,
                          opts: CombinedSolverOptions,
                          H_a=0.0,
                          dist_type=DistType.FD,
                          log_prefix="exp"):

    """Runs a single experiment and returns the result dataframe (from log) and final stats."""

    log_dir = Path("logs")
    log_dir.mkdir(exist_ok=True)

    run_id = uuid.uuid4().hex[:8]
    log_file_path = log_dir / f"{log_prefix}_{run_id}.csv"
    log_file_abs = str(log_file_path.resolve())

    opts_copy = copy.deepcopy(opts)
    opts_copy.sa.log_file = log_file_abs

    native_opts = opts_copy.to_qss()

    res_df = chemical_potential(
        [dist_type],
        [n_val],
        [beta_val],
        [H_a],
        e_k,
        solver_opts=native_opts
    )

    final_residual = res_df.iloc[0]['residual']

    try:
        log_df = pd.read_csv(log_file_abs)
        # Find iteration of best solution (minimum energy)
        min_energy = log_df['energy'].min()
        best_iter = log_df[log_df['energy'] == min_energy]['iteration'].iloc[0]

        # Convergence speed: Iterations to reach 1e-6
        threshold = 1e-6
        conv_df = log_df[log_df['energy'].abs() < threshold]
        if not conv_df.empty:
            conv_iter = conv_df['iteration'].iloc[0]
        else:
            conv_iter = log_df['iteration'].max()

    except Exception as e:
        print(f"Error reading log {log_file_abs}: {e}")
        best_iter = np.nan
        conv_iter = np.nan
        log_df = pd.DataFrame()

    return {
        "n": n_val,
        "beta": beta_val,
        "T": 1.0/beta_val,
        "final_error": final_residual,
        "best_iter": best_iter,
        "conv_iter": conv_iter,
        "log_df": log_df,
        "log_file": log_file_abs # Return the path so we can map results to files
    }
