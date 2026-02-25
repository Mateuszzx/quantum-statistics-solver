// Example executable: compute chemical potential mu(n) for several distributions
// Writes CSV: dist,n,mu,iterations,converged,error
// Then prints a short Python snippet for plotting.

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <iomanip>

#include "qss/chemical_potential.hpp"
#include "qss/solver.hpp"

static std::vector<double> linspace(double a, double b, std::size_t N)
{
    std::vector<double> v; v.reserve(N);
    if (N == 0) return v;
    if (N == 1) { v.push_back(a); return v; }
    double step = (b - a) / static_cast<double>(N - 1);
    for (std::size_t i = 0; i < N; ++i) v.push_back(a + step * static_cast<double>(i));
    return v;
}

int main()
{
    // Lattice / energy parameters (simple tight-binding demo)
    const int N = 100;                // lattice size per dimension
    const double t = -1.0;            // hopping
    const double t_prime = 0.25;      // next-nearest hopping
    const double beta = 10.0;         // inverse temperature
    const double H_a_T = 0.0;         // magnetic field (Tesla)
    const double g = 2.0;             // g-factor
    const double muB = qss::QSS_MU_B_EV_PER_T; // Bohr magneton (eV/T)

    // Fillings range
    double n_min = 0.10;
    double n_max = 0.99;
    std::size_t n_count = 60; // resolution
    auto n_values = linspace(n_min, n_max, n_count);

    // Build 2D energy lattice E(kx,ky)
    // k in [-pi,pi)
    auto k = linspace(-M_PI, M_PI, static_cast<std::size_t>(N));
    qss::Mat E(N, std::vector<double>(N));
    for (int ix = 0; ix < N; ++ix)
    {
        for (int iy = 0; iy < N; ++iy)
        {
            double kx = k[ix];
            double ky = k[iy];
            // Tight-binding dispersion with next-nearest term
            double eps = -2.0 * t * (std::cos(kx) + std::cos(ky))
                         -4.0 * t_prime * std::cos(kx) * std::cos(ky);
            E[ix][iy] = eps;
        }
    }

    // Combined solver options setup
    qss::CombinedSolverOptions opts;
    opts.bracketer_type = qss::BracketerType::Expansion;
    opts.root.abs_tol = 1e-9;
    opts.root.rel_tol = 1e-6;
    opts.root.max_iters = 5000;
    opts.bracketer_opts.expand_factor = 1.4;
    opts.bracketer_opts.initial_step = 10.0;
    opts.bracketer_opts.max_expands = 2;

    // Simulated annealing parameters (only used if solver switched below)
    opts.sa.T0 = 1.0;          // initial temperature
    opts.sa.T_min = 1e-8;      // stopping temperature
    opts.sa.alpha = 0.80;      // decay factor for exponential schedule
    opts.sa.sched = qss::CoolingSchedule::Exponential;
    opts.sa.max_iters = 50000; // global cap
    opts.sa.min_iters = 70;   // ensure minimum exploration
    opts.sa.stall_iters = 50; // terminate if no improvement
    opts.sa.step_sigma = 12.0; // proposal stddev
    opts.sa.step_decay = 35;   // iterations between sigma decay
    opts.sa.scale_sigma_by_T = true; // adapt step size with temperature
    opts.sa.bounds = {{-5.0, 5.0}};  // 1D search bounds for mu
    opts.sa.reflect_bounds = true;   // reflect when out of bounds
    opts.sa.seed = 0;                // deterministic runs if non-zero

    // Choose solver: default Brent, override via env QSS_SOLVER=SA
    const char* solver_env = "SA";
    if (solver_env && std::string(solver_env) == "SA") {
        opts.solver = qss::SolverType::SimulatedAnnealing;
        std::cout << "[qss] Using Simulated Annealing solver configuration\n";
    } else {
        opts.solver = qss::SolverType::Brent;
        std::cout << "[qss] Using Brent root-finding solver configuration\n";
    }

    // Bracket window and initial guess for chemical potential root-finding
    double a0 = -5.0;
    double b0 = 5.0;
    double x0 = 0.0;

    struct DistLabel { qss::DistType dist; const char* label; };
    std::vector<DistLabel> distros = {
        { qss::DistType::FD,   "FD" },
        { qss::DistType::FLB,  "FLB" },
        { qss::DistType::NFLB, "NFLB" }
    };

    // Prepare output CSV
    const std::string out_path = "mu_vs_n.csv";
    std::ofstream ofs(out_path);
    if (!ofs.is_open())
    {
        std::cerr << "Failed to open output file: " << out_path << std::endl;
        return 1;
    }
    ofs << "dist,n,mu,iterations,converged,error\n";
    ofs << std::setprecision(12);

    // Loop over distributions and fillings
    for (const auto &dl : distros)
    {
        for (double n_fill : n_values)
        {
            auto R = qss::chemical_potential_B(
                dl.dist,
                beta,
                n_fill,
                static_cast<double>(N),
                E,
                H_a_T,
                0.0, // U (not used for these distros)
                g,
                muB,
                opts,
                a0,
                b0,
                x0);
            ofs << dl.label << ','
                << n_fill << ','
                << R.mu << ','
                << R.iterations << ','
                << (R.converged ? 1 : 0) << ','
                << R.error << '\n';
        }
    }
    ofs.close();

    std::cout << "Wrote chemical potential data to " << out_path << "\n";
    std::cout << "Plot with Python (after building bindings):\n";
    std::cout << "\n";
    std::cout << "import pandas as pd\n"
              << "import seaborn as sns, matplotlib.pyplot as plt\n"
              << "df = pd.read_csv('" << out_path << "')\n"
              << "sns.lineplot(data=df, x='n', y='mu', hue='dist')\n"
              << "plt.show()\n";

    return 0;
}
