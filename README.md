# Quantum Statistics Solver

High-performance library for quantum statistics calculations, focusing on strongly correlated electron systems.

## Features

- **Distributions**: FD, FLB, NFL, NFLB, SSL, HUBBARD.
- **Thermodynamics**: Chemical potential ($\mu$), entropy ($S$), heat capacity ($C$), energy ($E$), free energy ($F$), magnetization ($M$), $n_{kd}$.
- **Solvers**: Brent root-finding and simulated annealing for $\mu(n)$; optional **prescribed $\mu$** (skip the root solve).
- **Performance**: Core logic in C++17 via `nanobind`.
- **GUI**: Streamlit app for sweeps, panels, and BZ occupancy plots.

## GUI (Streamlit)

Interactive app under `gui/` (entry point: `gui/web_app.py`).

Sweeps over $n$ (or prescribed $\mu$), $k_B T/|t|$, $H_a$, and $U/|t|$, then plots registered derived quantities ($\mu$, $S$, $C$, $E$, $F$, $M$, $n_{kd}$).

### Run

From `quantum-statistics-solver/`:

```bash
make gui-launch
```

Equivalent:

```bash
uv run streamlit run gui/web_app.py
```

Widget defaults live in `gui/defaults.yaml` (override with `QSS_GUI_DEFAULTS=/path/to.yaml`).

### What the GUI shows

- **Panels**: quantity / X axis / panel sweep / fixed parameters; optional **Prescribe $\mu$** (skips $\mu(n)$ solve); solver options in a sidebar expander when solving from $n$.
- **Distributions**: standalone occupancy map / 1D cuts; $\mu$ solved from $n$ or prescribed. Uses `gui/scientific_plot.py` for multi-curve 1D BZ cuts (does **not** require Panel experiment results).
- **Plot style**: panels per row, figure size, optional axis limits, linear/log X scale.

### Output data

- `data/gui_results.csv` — chemical-potential grid from Panels
- Download buttons for the active quantity in the Panels data expander
- Logs under `logs/`

More detail: `docs/gui.rst`.

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
# or full dev env:
make dev-install
```

## Usage

### Python API

```python
import numpy as np
from qss import chemical_potential, DistType
from qss.distributions import prepare_2D_lattice_data

N = 100
kx = np.linspace(-np.pi, np.pi, N)
ky = np.linspace(-np.pi, np.pi, N)
e_k, _, _ = prepare_2D_lattice_data(kx, ky, t=-1.0, t_prim=0.25)

# Solve μ from band filling n
df = chemical_potential(
    distros=[DistType.FD, DistType.FLB],
    n_vals=np.linspace(0.1, 0.9, 20),
    kBT_vals=[0.1, 0.05],
    H_a_vals=[0.0],
    E=e_k,
    U_vals=[0.0],
)

# Or prescribe μ and skip the root solve (n is computed from occupancy)
df_mu = chemical_potential(
    distros=[DistType.FD],
    mu_vals=[-1.0, 0.0, 1.0],
    kBT_vals=[0.1],
    H_a_vals=[0.0],
    E=e_k,
)

print(df.head())
```

### Plotting

Panel / publication-style helpers used by the GUI live in `gui/panel_plot.py` and `gui/scientific_plot.py`.

## Development

### Setup

```bash
make dev-install
source .venv/bin/activate
```

### Testing & Coverage

```bash
make test-python
make test-cpp
make coverage-all
```

### Code Quality

```bash
make check          # Python checks + clang-tidy
make precommit-run  # pre-commit hooks
```

### Documentation

System packages (Debian/Ubuntu):

```bash
make install
# or: sudo apt-get install -y doxygen graphviz python3-venv python3-dev
```

```bash
make docs
make docs-serve   # http://localhost:8000
```

HTML output: `build/docs/html/index.html`.

## Project Structure

```txt
quantum-statistics-solver/
├── src/                     # C++ source
├── include/                 # C++ headers
├── qss/                     # Python package
├── gui/                     # Streamlit app (web_app.py, defaults.yaml, …)
├── utils/                   # Experiment helpers (e.g. run_experiment)
├── scripts/                 # Thin wrappers used by CI (→ Makefile)
├── tests/
│   ├── cpp/
│   └── python/
├── docs/
├── .github/workflows/
├── Makefile
├── pyproject.toml
├── uv.lock
└── .pre-commit-config.yaml
```

## CI/CD Pipeline

GitHub Actions (`.github/workflows/`):

- **ci.yml** — build and test
- **coverage.yml** — coverage
- **pre-commit.yml** — hooks
- **documentation.yml** — docs build / deploy

## Citing

**If you use this software in research, please cite the accompanying article** (and, where appropriate, this repository).

- Machine-readable citation: [`CITATION.cff`](CITATION.cff) (GitHub → *Cite this repository*)
- BibTeX: [`docs/references.bib`](docs/references.bib)

```bibtex
@misc{kuterba2026nonstandardquantumstatisticsstrongly,
      title={Nonstandard quantum statistics for strongly correlated systems: Superexclusion in momentum space and comparison with exact results},
      author={Piotr Kuterba and Mateusz Wójcik and Józef Spałek},
      year={2026},
      eprint={2611.26489},
      archivePrefix={arXiv},
      primaryClass={cond-mat.str-el},
      url={https://arxiv.org/abs/2611.26489},
}
```
