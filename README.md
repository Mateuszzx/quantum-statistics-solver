# Quantum Statistics Solver

[![CI](https://github.com/YOUR_USERNAME/quantum-statistics-solver/workflows/CI/badge.svg)](https://github.com/YOUR_USERNAME/quantum-statistics-solver/actions)
[![codecov](https://codecov.io/gh/YOUR_USERNAME/quantum-statistics-solver/branch/main/graph/badge.svg)](https://codecov.io/gh/YOUR_USERNAME/quantum-statistics-solver)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

High-performance library for quantum statistics calculations.

## Quick Start

```bash
# Clone repository
git clone https://github.com/your-org/quantum-statistics-solver.git
cd quantum-statistics-solver

# Set up development environment
make dev-install

# Activate virtual environment
source .venv/bin/activate

# Build and test
make build
make test

# Run quality checks
make check-all
```

## Development Tools

This project uses modern development tools configured in `pyproject.toml`:

| Tool | Purpose | Command | Config |
|------|---------|---------|--------|
| **pytest** | Test runner | `make test-python` | `[tool.pytest]` |
| **coverage** | Code coverage | `make coverage-python` | `[tool.coverage]` |
| **ruff** | Fast linter | `make lint-python` | `[tool.ruff]` |
| **black** | Code formatter | `make format-python` | `[tool.black]` |
| **mypy** | Type checker | `make typecheck-python` | `[tool.mypy]` |
| **codespell** | Spell checker | `make codespell` | `[tool.codespell]` |
| **pre-commit** | Git hooks | `make precommit-install` | `.pre-commit-config.yaml` |

### Code Coverage

Coverage reports are generated for both C++ and Python code:

```bash
# Generate coverage reports
make coverage-all

# View reports
open build/coverage/cpp/index.html      # C++ coverage
open build/coverage/python/index.html   # Python coverage
```

Coverage is automatically tracked in CI and uploaded to [Codecov](https://codecov.io):
- Configuration: `codecov.yml`
- CI integration: `.github/workflows/ci.yml`
- Target: 80% project coverage, 70% patch coverage

## Project Structure

```
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

- **ubuntu.yml**: Main CI pipeline (Linux, macOS)
- **windows.yml**: Windows builds (MSVC, MinGW)
- **pre-commit.yml**: Code quality checks
- **documentation.yml**: Build and deploy docs

All workflows upload coverage to Codecov for tracking.
