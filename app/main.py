import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(style="whitegrid")

from qss.solver_options import brent_defaults, sa_defaults
from qss.chemical_potential import solve_chemical_potential
from qss.distributions import DistType, prepare_2D_lattice_data

if __name__ == "__main__":
    print("Hello_world")
    MU_B  = 5.788381806e-5  # eV/T
    N = 100
    t = -1
    t_prim = 0.25

    kx = np.linspace(-np.pi, np.pi, N, endpoint=False)
    ky = np.linspace(-np.pi, np.pi, N, endpoint=False)
    e_k, _, _ = prepare_2D_lattice_data(kx, ky, t, t_prim)


    # Band filling values
    n_min = 0.1
    n_max = 0.99
    n_count = 40

    # Magnetic field values
    Ha_min = 0.0
    Ha_max = 10.0
    Ha_count = 5

    # Beta values
    beta_min = 1
    beta_max = 2
    beta_count = 2

    distros = [DistType.FD, DistType.FLB, DistType.NFLB]  # ND
    # n  = np.array([0.5, 0.8, 0.9], dtype=np.float64)  # NN
    n  = np.linspace(n_min, n_max, n_count, dtype=np.float64)  # NN
    # beta = np.array([0.1, 1.0], dtype=np.float64)  # NB
    beta = np.logspace(beta_min, beta_max, beta_count, dtype=np.float64)  # NB
    H_a = np.array([0.0], dtype=np.float64) # NH
    # H_a = np.linspace(Ha_min, Ha_max, Ha_count, dtype=np.float64) # NH

    py_opts = brent_defaults(abs_tol=1e-9, rel_tol=1e-6, expand_factor=1.42, initial_step=11, max_expands=1, max_iters=5000)
    native_opts= py_opts.to_qss()

    mu_df = solve_chemical_potential(
        distros=distros,
        n_assumed=n,
        beta=beta,
        H_a=H_a,
        N=N,
        E=e_k,
        solver_opts=native_opts,
        x0=10.,
        a0=-5.0,
        b0=5.0,
    )
    mu_df.to_csv("chemical_potential_results.csv")

    data = mu_df.loc[(["FD","FLB", "NFLB"], n, beta[-1], H_a[0]), :]
    fig,ax = plt.subplots()
    sns.lineplot(data= data, x="n", y="mu", hue="dist", markers=True, ax=ax)
    ax.set_ylim(-3,3)
    plt.show()
