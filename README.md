# Quantum Statistics Solver

High-performance library for quantum statistics calculations, focusing on strongly correlated electron systems.

## Features

- **Distributions**: Fermi-Dirac (FD), Fermi Liquid with Boltzmann (FLB), Non-Fermi Liquid with Boltzmann (NFLB).
- **Thermodynamics**: Efficient calculation of Chemical Potential ($\mu$), Entropy ($S$), and Energy ($E$).
- **Solvers**: Robust root-finding (Brent's method) and global optimization (Simulated Annealing) for solving $\mu(n)$.
- **Performance**: Core logic implemented in C++17 for speed, exposed to Python via `nanobind`.
- **Visualization**: Built-in plotting tools for scientific publication-quality figures.

## GUI (Streamlit)

The repository includes an interactive GUI at `gui/web_app.py` that runs the same “notebook-style” sweeps over:

- band filling \(n\)
- temperature \(k_B T/|t|\)
- magnetic field \(H_a\)

and computes \(\mu\), \(S\), \(E\), \(F\), \(M\), and \(n_{kd}\).

### Run

From the `quantum-statistics-solver/` directory:

```bash
streamlit run gui/web_app.py
```

If you use the project `Makefile`, there may also be a helper target (e.g. `make gui-launch`).

### What the GUI shows

- **Panels tab**: pick the **quantity**, **x-axis variable**, **panel sweep variable**, and **fixed parameters**; plot styling supports:
  - configurable panels-per-row and figure size
  - optional x/y axis limits
  - per-distribution colors/linestyles matching `tools/scientific_plot_v2.py`
  - linear/log x-scale toggle
- **Distributions tab**:
  - notebook-style distribution cuts (when `tools/scientific_plot_v2.py` is available)
  - a **2D occupancy heatmap** (2D lattice) computed using the **solved \(\mu\)** from results
- **Solver parameters**: solver configuration is available in a collapsible section in the main window.

### Output data

The GUI writes the chemical potential solve grid to:

- `data/gui_results.csv`

and also provides download buttons for derived quantities in the Panels tab.

## Installation

### Prerequisites

- C++17 compatible compiler
- CMake >= 3.15
- Python >= 3.11
- [uv](https://docs.astral.sh/uv/)

### Install from Source

```bash
git clone https://github.com/Mateuszzx/quantum-statistics-solver.git
cd quantum-statistics-solver
uv venv
uv pip install --python .venv/bin/python -e .
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

Use the included scientific plotting tools (see `tools/scientific_plot_v2.py` in the parent `qss-results` checkout):

```python
import matplotlib.pyplot as plt

plt.show()
```

## Development

This project uses a `Makefile` to manage development tasks.

### Setup

Requires [uv](https://docs.astral.sh/uv/).

```bash
make dev-install
source .venv/bin/activate
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

### Documentation

#### Prerequisites (system packages)

The documentation build uses:

- **Sphinx** (installed into the project venv by `make docs`)
- **Doxygen** + **Graphviz** (optional, for C++ API pages via Breathe)

On Debian/Ubuntu, you can install the system dependencies with:

```bash
make install
```

or manually:

```bash
sudo apt-get install -y doxygen graphviz python3-venv python3-dev
```

LaTeX is **not required** (Doxygen formulas are rendered via MathJax in HTML).

Build the Sphinx (and optional Doxygen) docs:

```bash
make docs
```

Serve the generated HTML on `http://localhost:8000`:

```bash
make docs-serve
```

Docs are written to `build/docs/html/index.html`.

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
