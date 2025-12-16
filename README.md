# Quantum Statistics Solver

[![CI](https://github.com/YOUR_USERNAME/quantum-statistics-solver/workflows/CI/badge.svg)](https://github.com/YOUR_USERNAME/quantum-statistics-solver/actions)
[![codecov](https://codecov.io/gh/YOUR_USERNAME/quantum-statistics-solver/branch/main/graph/badge.svg)](https://codecov.io/gh/YOUR_USERNAME/quantum-statistics-solver)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

High-performance library for quantum statistics calculations, focusing on strongly correlated electron systems.

## Features

- **Distributions**: Fermi-Dirac (FD), Fermi Liquid with Boltzmann (FLB), Non-Fermi Liquid with Boltzmann (NFLB).
- **Thermodynamics**: Efficient calculation of Chemical Potential ($\mu$), Entropy ($S$), and Energy ($E$).
- **Solvers**: Robust root-finding (Brent's method) and global optimization (Simulated Annealing) for solving $\mu(n)$.
- **Performance**: Core logic implemented in C++17 for speed, exposed to Python via `nanobind`.
- **Visualization**: Built-in plotting tools for scientific publication-quality figures.

## Installation

### Prerequisites
- C++17 compatible compiler
- CMake >= 3.15
- Python >= 3.10

### Install from Source

```bash
git clone https://github.com/Mateuszzx/quantum-statistics-solver.git
cd quantum-statistics-solver
pip install -e .
```

## Usage

### Python API

Calculate chemical potential for different distributions:

```python
import numpy as np
from qss import chemical_potential, DistType
from qss.distributions import prepare_2D_lattice_data

# 1. Setup Lattice
N = 100
kx = np.linspace(-np.pi, np.pi, N)
ky = np.linspace(-np.pi, np.pi, N)
e_k, _, _ = prepare_2D_lattice_data(kx, ky, t=-1.0, t_prim=0.25)

# 2. Define Parameters
n_vals = np.linspace(0.1, 0.9, 20)
beta_vals = [10.0, 20.0] # 1/kBT
H_a_vals = [0.0]         # Magnetic field

# 3. Calculate Chemical Potential
df = chemical_potential(
    distros=[DistType.FD, DistType.FLB],
    n_vals=n_vals,
    beta_vals=beta_vals,
    H_a_vals=H_a_vals,
    E=e_k
)

print(df.head())
```

### Plotting

Use the included scientific plotting tools:

```python
import tools.scientific_plot as sp
import matplotlib.pyplot as plt

# Plot Entropy vs Temperature
sp.plot_results_temperature(
    results=df,
    n=[0.5, 0.8],
    betas=beta_vals,
    y_label=r"$\mu$",
    distros=["FD", "FLB"]
)
plt.show()
```

## Development

This project uses a `Makefile` to manage development tasks.

### Setup

```bash
make venv
source .venv/bin/activate
make dev-install
```

### Testing & Coverage

```bash
# Run Python tests
make test-python

# Run C++ tests
make test-cpp

# Generate coverage reports (HTML)
make coverage-all
```

### Code Quality

```bash
make check-all  # Runs linting, formatting, and type checking
```

## Project Structure

```txt
quantum-statistics-solver/
├── src/                    # C++ source code
├── include/                # C++ headers
├── qss/                    # Python package
├── tests/
│   ├── cpp/               # C++ tests (Google Test)
│   └── python/            # Python tests (pytest)
├── docs/                  # Documentation
├── tools/                 # Development scripts
├── .github/workflows/     # CI/CD pipelines
├── Makefile               # Build automation
├── pyproject.toml         # Python config & dependencies
├── codecov.yml            # Coverage reporting config
└── .pre-commit-config.yaml # Git hooks config
```

## CI/CD Pipeline

GitHub Actions workflows (`.github/workflows/`):

- **ci.yml**: Main CI pipeline (Linux)
- **coverage.yml**: Code coverage checks
- **pre-commit.yml**: Code quality checks
- **documentation.yml**: Build and deploy docs

All workflows upload coverage to Codecov for tracking.
