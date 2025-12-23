from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

from app.BO.utils.experiment_utils import run_sa_experiment
from qss.distributions import DistType, prepare_2D_lattice_data
from qss.solver_options import CoolingSchedule, sa_defaults
import tools.scientific_plot as sp

# Set page config
st.set_page_config(page_title="Quantum Statistics Solver", layout="wide")

st.title("Quantum Statistics Solver - Chemical Potential Analysis")

# --- Sidebar Parameters ---
st.sidebar.header("Parameters")

# Lattice
st.sidebar.subheader("Lattice")
N = st.sidebar.number_input("N (Lattice size)", min_value=10, max_value=10000, value=100, step=10)
t = st.sidebar.number_input("t (Hopping parameter)", value=-1.0)
t_prim = st.sidebar.number_input("t' (Next-nearest hopping)", value=0.25)

# Band Filling
st.sidebar.subheader("Band Filling (n)")
col1, col2 = st.sidebar.columns(2)
n_min = col1.number_input("n min", min_value=0.0, max_value=1.0, value=0.01)
n_max = col2.number_input("n max", min_value=0.0, max_value=1.0, value=0.99)
n_count = st.sidebar.number_input("n steps", min_value=2, max_value=100, value=10)

# Physics
st.sidebar.subheader("Physics")
beta_val = st.sidebar.number_input("Beta (1/kBT)", min_value=0.1, value=1.0)
H_a_val = st.sidebar.number_input("H_a (Magnetic Field)", value=0.0)

# Distributions
st.sidebar.subheader("Distributions")
dist_map = {
    "FD": DistType.FD,
    "FLB": DistType.FLB,
    "NFL": DistType.NFL,
    "NFLB": DistType.NFLB,
    "SSL": DistType.SSL,
}
selected_dist_names = st.sidebar.multiselect(
    "Select Distributions", options=list(dist_map.keys()), default=list(dist_map.keys())
)
distros = [dist_map[name] for name in selected_dist_names]

# Solver Options
with st.sidebar.expander("Solver Options"):
    max_iters = st.number_input("Max Iterations", value=50000)
    min_iters = st.number_input("Min Iterations", value=2000)
    stall_iters = st.number_input("Stall Iterations", value=1000)

    T0 = st.number_input("Initial Temp (T0)", value=1000.0)
    T_min = st.number_input("Min Temp (T_min)", value=1e-9, format="%.1e")

    alpha = st.number_input(
        "Cooling Rate (alpha)", value=0.999, min_value=0.0, max_value=1.0, format="%.4f"
    )
    step_decay = st.number_input(
        "Step Decay", value=0.999, min_value=0.0, max_value=1.0, format="%.4f"
    )
    step_sigma = st.number_input("Step Sigma", value=5.0)

    sched_options = {
        "Exponential": CoolingSchedule.Exponential,
        "Linear": CoolingSchedule.Linear,
        "Logarithmic": CoolingSchedule.Logarithmic,
        "Sawtooth": CoolingSchedule.Sawtooth,
    }
    sched_name = st.selectbox("Cooling Schedule", options=list(sched_options.keys()))
    sched = sched_options[sched_name]

    scale_sigma_by_t = st.checkbox("Scale Sigma by T", value=True)
    reflect_bounds = st.checkbox("Reflect Bounds", value=True)

# Search Configuration
st.sidebar.subheader("Search Configuration")
search_min = st.sidebar.number_input("Search Min (a0)", value=-10.0)
search_max = st.sidebar.number_input("Search Max (b0)", value=10.0)

use_random_x0 = st.sidebar.checkbox("Random Initial Guess", value=False)
if not use_random_x0:
    x0_val = st.sidebar.number_input("Initial Guess (x0)", value=0.0)
else:
    # Generate a random value now, but note it will be fixed for the run unless we move it inside the loop
    # For now, we pass a single value.
    x0_val = float(np.random.uniform(search_min, search_max))
    st.sidebar.info(f"Random x0: {x0_val:.4f}")

# --- Main Execution ---

if st.sidebar.button("Run Simulation", type="primary"):
    if not distros:
        st.error("Please select at least one distribution.")
    else:
        with st.spinner("Running simulation..."):
            try:
                # Prepare Data
                kx = np.linspace(-np.pi, np.pi, N, endpoint=False)
                ky = np.linspace(-np.pi, np.pi, N, endpoint=False)
                e_k, _, _ = prepare_2D_lattice_data(kx, ky, t, t_prim)

                n_vals = np.linspace(n_min, n_max, n_count, dtype=np.float64)
                beta_vals = np.array([beta_val], dtype=np.float64)
                H_a_vals = np.array([H_a_val], dtype=np.float64)

                # Solver Config
                py_opts = sa_defaults(
                    T0=T0,
                    T_min=T_min,
                    alpha=alpha,
                    max_iters=int(max_iters),
                    min_iters=int(min_iters),
                    stall_iters=int(stall_iters),
                    step_decay=step_decay,
                    step_sigma=step_sigma,
                    bounds=[(float(search_min), float(search_max))],
                    scale_sigma_by_T=scale_sigma_by_t,
                    reflect_bounds=reflect_bounds,
                    sched=sched,
                )
                native_opts = py_opts.to_qss()

                # Output paths
                logs_dir = Path("logs")
                logs_dir.mkdir(exist_ok=True)
                output_file = Path("data/gui_results.csv")
                output_file.parent.mkdir(exist_ok=True)

                # Run Experiment
                mu_df = run_sa_experiment(
                    n_vals=n_vals,
                    beta_vals=beta_vals,
                    h_a_vals=H_a_vals,
                    distros=distros,
                    e_k=e_k,
                    N=N,
                    solver_opts=native_opts,
                    logs_dir=str(logs_dir),
                    output_file=str(output_file),
                    x0=x0_val,
                    a0=search_min,
                    b0=search_max,
                )

                st.success("Simulation Complete!")

                # Layout: Plot on left (or top), Data below
                st.subheader("Chemical Potential vs Band Filling")

                # Plotting
                plt.close("all")

                # Call the plotting function from tools
                sp.plot_result_n(
                    results=mu_df,
                    n=n_vals,
                    betas=beta_vals,
                    y_label=r"chemical potential, $\mu / k_B N$",
                    distros=distros,
                    H_a=H_a_val,
                    y_col="mu",
                    legend_position="center left",
                    figsize=(10, 6),
                    cols=1,
                )

                # Display the plot
                st.pyplot(plt.gcf())

                # Display Data
                with st.expander("Show Raw Data"):
                    st.dataframe(mu_df)

                    csv = mu_df.to_csv(index=False).encode("utf-8")
                    st.download_button(
                        "Download CSV",
                        csv,
                        "chemical_potential_results.csv",
                        "text/csv",
                        key="download-csv",
                    )

            except Exception as e:
                st.error(f"An error occurred: {e}")
                # Print traceback for debugging
                import traceback

                st.text(traceback.format_exc())
else:
    st.info("Adjust parameters in the sidebar and click 'Run Simulation' to start.")
