# ==============================================================================
# Quantum Statistics Solver - Project Makefile
# ==============================================================================
# This Makefile provides convenient targets for building, testing, and
# maintaining the quantum-statistics-solver project.
#
# Quick Start:
#   make help                # Show all available targets
#   make dev-install         # Set up development environment
#   make build              # Build C++ library and Python bindings
#   make test               # Run all tests (C++ and Python)
#   make check-all          # Run all quality checks
#
# Documentation: https://github.com/your-org/quantum-statistics-solver
# ==============================================================================

SHELL := /bin/bash -o pipefail
.DEFAULT_GOAL := help

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

# Build directories
BUILD_DIR ?= build
SRC_DIR   ?= src
DOCS_DIR  ?= docs

# Python configuration
VENV_DIR  ?= .venv
PYTHON    ?= $(VENV_DIR)/bin/python
PIP       ?= $(VENV_DIR)/bin/pip

# Documentation output
DOCS_HTML := $(BUILD_DIR)/docs/html
DOXY_OUT  := $(BUILD_DIR)/doxygen

# Build configuration
CMAKE_GENERATOR ?=
CMAKE_BUILD_TYPE ?= Release

# Tool configuration
CLANG_TIDY_PY ?= python tools/run-clang-tidy.py
CLANG_TIDY_JOBS ?= $(shell nproc 2>/dev/null || echo 4)

# Source files for linting
SRC_FILES := $(shell find $(SRC_DIR) -type f \( -name '*.c' -o -name '*.cpp' -o -name '*.h' -o -name '*.hpp' \) 2>/dev/null)

# ------------------------------------------------------------------------------
# Phony Targets Declaration
# ------------------------------------------------------------------------------
.PHONY: help all clean install \
        venv dev-install \
        build configure test test-cpp test-python \
        coverage coverage-cpp coverage-python coverage-all coverage-clean \
        lint lint-cpp lint-python lint-all \
        format format-cpp format-python format-all \
        typecheck-python codespell codespell-fix \
        check-python check-all fix-python \
        docs docs-clean docs-serve \
        precommit-install precommit-run precommit-ci \
        prepare compile-commands

# ==============================================================================
# Help & Information
# ==============================================================================

help: ## Show this help message with all available targets
	@echo "Quantum Statistics Solver - Makefile"
	@echo "===================================="
	@echo ""
	@awk 'BEGIN{FS=":.*##"; section=""} \
		/^## /{section=substr($$0,4); printf "\n\033[1m%s\033[0m\n", section; next} \
		/^[a-zA-Z0-9_-]+:.*##/{printf "  \033[36m%-22s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)
	@echo ""

all: build test ## Default: build and test everything

# ==============================================================================
## Environment Setup
# ==============================================================================

install: ## Install system dependencies (Debian/Ubuntu only)
	@echo "[info] Installing system packages (requires sudo)"
	sudo apt-get update -y
	sudo apt-get install -y --no-install-recommends \
		gcc g++ cmake make git python3-dev python3-venv \
		clang-format clang-tidy doxygen graphviz \
		lcov gcovr

venv: ## Create Python virtual environment
	@echo "[info] Creating virtual environment: $(VENV_DIR)"
	@if [ ! -x "$(PYTHON)" ]; then \
		python3 -m venv "$(VENV_DIR)"; \
		$(PIP) install --upgrade pip wheel setuptools; \
	fi

dev-install: venv ## Install package in development mode with all dependencies
	@echo "[info] Installing development environment"
	$(PIP) install -e ".[dev,test,docs]"
	@echo "[info] Development environment ready!"
	@echo "      Activate with: source $(VENV_DIR)/bin/activate"

# ==============================================================================
## Build System
# ==============================================================================

prepare: ## Create/refresh build directory
	@echo "[info] Preparing build directory: $(BUILD_DIR)"
	@mkdir -p "$(BUILD_DIR)"

configure: prepare ## Configure CMake build system
	@echo "[info] Configuring CMake ($(CMAKE_BUILD_TYPE))"
	@rm -f "$(BUILD_DIR)/CMakeCache.txt"
	cmake -S . -B "$(BUILD_DIR)" \
		$(if $(CMAKE_GENERATOR),-G "$(CMAKE_GENERATOR)",) \
		-DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) \
		-DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
		-DENABLE_PY_MODULE=ON \
		-DENABLE_TESTING=ON

reconfigure: ## Force reconfigure (cleans CMake cache)
	@echo "[info] Force reconfiguring (cleaning cache)"
	@rm -rf "$(BUILD_DIR)/CMakeCache.txt" "$(BUILD_DIR)/CMakeFiles"
	@$(MAKE) configure

build: configure ## Build C++ library and Python bindings
	@echo "[info] Building project"
	cmake --build "$(BUILD_DIR)" -j$(CLANG_TIDY_JOBS)

compile-commands: ## Generate compile_commands.json for IDE/tools
	@echo "[info] Generating compile_commands.json"
	@$(MAKE) configure CMAKE_BUILD_TYPE=Debug

clean: ## Remove all build artifacts and caches
	@echo "[info] Cleaning build artifacts"
	@rm -rf "$(BUILD_DIR)" dist .coverage* .pytest_cache .mypy_cache .ruff_cache
	@find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	@find . -type f -name "*.pyc" -delete 2>/dev/null || true

# ==============================================================================
## Testing
# ==============================================================================

test: test-cpp test-python ## Run all tests (C++ and Python)

test-cpp: build ## Run C++ tests with Google Test
	@echo "[info] Running C++ tests"
	@cd "$(BUILD_DIR)" && ctest --output-on-failure -j$(CLANG_TIDY_JOBS)

test-python: venv ## Run Python tests with pytest
	@echo "[info] Running Python tests"
	$(PYTHON) -m pytest tests/python/ -v

# ==============================================================================
## Code Coverage
# ==============================================================================
# Coverage reports track which code is exercised by tests.
#
# Local usage:
#   make coverage-python     # Generate Python coverage report
#   make coverage-cpp        # Generate C++ coverage report
#   make coverage-all        # Generate both reports
#   open build/coverage/python/index.html  # View in browser
#
# CI Integration (automatic):
#   - Coverage reports uploaded to Codecov.io in CI workflows
#   - Codecov posts comments on pull requests
#   - View dashboard: https://codecov.io/gh/<your-username>/<repo>
#
# Configuration:
#   - Python coverage: pyproject.toml [tool.coverage.*]
#   - C++ coverage: CMakeLists.txt (ENABLE_COVERAGE option)
#   - Codecov behavior: codecov.yml
#
# Manual upload (if needed):
#   codecov -f build/coverage/coverage.xml -t YOUR_TOKEN

coverage-python: venv ## Generate Python coverage report (pytest-cov)
	@echo "[info] Running Python tests with coverage"
	@mkdir -p "$(BUILD_DIR)/coverage/python"
	@PYTHONPATH="$(CURDIR)" $(PYTHON) -m coverage run \
        --source="$(CURDIR)/qss" \
        -m pytest tests/python/ -v
	@$(PYTHON) -m coverage report -m
	@$(PYTHON) -m coverage html -d "$(BUILD_DIR)/coverage/python"
	@$(PYTHON) -m coverage xml -o "$(BUILD_DIR)/coverage/coverage.xml"
	@echo ""
	@echo "=================================="
	@echo "Python Coverage Report Generated"
	@echo "=================================="
	@echo "HTML: $(BUILD_DIR)/coverage/python/index.html"
	@echo "XML:  $(BUILD_DIR)/coverage/coverage.xml"
	@echo ""
	@echo "To upload to Codecov manually:"
	@echo "  codecov -f $(BUILD_DIR)/coverage/coverage.xml -t YOUR_TOKEN"
	@echo ""

coverage-cpp: ## Generate C++ coverage report (requires lcov/gcovr)
	@echo "[info] Building with coverage instrumentation"
	@cmake -S . -B "$(BUILD_DIR)" \
		-DCMAKE_BUILD_TYPE=Debug \
		-DENABLE_COVERAGE=ON \
		-DCMAKE_CXX_FLAGS="--coverage -O0 -g" \
		-DCMAKE_EXE_LINKER_FLAGS="--coverage"
	@cmake --build "$(BUILD_DIR)" -j$(CLANG_TIDY_JOBS)
	@echo "[info] Running tests to generate coverage data"
	@cd "$(BUILD_DIR)" && ctest --output-on-failure || true
	@echo "[info] Generating coverage report"
	@mkdir -p "$(BUILD_DIR)/coverage/cpp"
	@if command -v gcovr >/dev/null 2>&1; then \
		gcovr -r . --object-directory "$(BUILD_DIR)" \
			--html --html-details \
			-o "$(BUILD_DIR)/coverage/cpp/index.html" \
			--exclude 'tests/.*' --exclude 'external/.*'; \
	else \
		lcov --capture --directory "$(BUILD_DIR)" \
			--output-file "$(BUILD_DIR)/coverage/cpp/coverage.info" \
			--rc geninfo_unexecuted_blocks=1; \
		lcov --remove "$(BUILD_DIR)/coverage/cpp/coverage.info" \
			'/usr/*' '*/tests/*' '*/external/*' \
			--output-file "$(BUILD_DIR)/coverage/cpp/coverage.cleaned.info"; \
		genhtml "$(BUILD_DIR)/coverage/cpp/coverage.cleaned.info" \
			--output-directory "$(BUILD_DIR)/coverage/cpp"; \
	fi
	@echo "[info] C++ coverage report: $(BUILD_DIR)/coverage/cpp/index.html"

coverage-all: coverage-cpp coverage-python ## Generate both C++ and Python coverage
	@echo ""
	@echo "=================================="
	@echo "Coverage Reports Generated"
	@echo "=================================="
	@echo "C++:    $(BUILD_DIR)/coverage/cpp/index.html"
	@echo "Python: $(BUILD_DIR)/coverage/python/index.html"
	@echo ""
	@echo "Upload to Codecov (CI only):"
	@echo "  codecov -f $(BUILD_DIR)/coverage/coverage.xml"
	@echo ""

coverage-clean: ## Clean coverage artifacts
	@echo "[info] Cleaning coverage data"
	@rm -rf "$(BUILD_DIR)/coverage" .coverage* coverage.xml
	@find . -name "*.gcda" -delete
	@find . -name "*.gcno" -delete

# ==============================================================================
## Code Quality - C++
# ==============================================================================

lint-cpp: compile-commands ## Run clang-tidy on C++ code (non-fatal)
	@echo "[info] Running clang-tidy (non-fatal)"
	@$(CLANG_TIDY_PY) \
		-p "$(BUILD_DIR)" -j $(CLANG_TIDY_JOBS) \
		-header-filter='(include|src|tests)/' -quiet \
		$(SRC_FILES) || true

lint-cpp-ci: compile-commands ## Run clang-tidy with -Werror (CI mode)
	@echo "[info] Running clang-tidy (CI mode: fail on warnings)"
	@$(CLANG_TIDY_PY) \
		-p "$(BUILD_DIR)" -j $(CLANG_TIDY_JOBS) \
		-header-filter='(include|src|tests)/' \
		-warnings-as-errors='*' \
		$(SRC_FILES)

format-cpp: ## Format C++ code with clang-format
	@echo "[info] Formatting C++ code"
	@python tools/run-clang-format.py -r $(SRC_DIR) include/ || true

# ==============================================================================
## Code Quality - Python
# ==============================================================================
# Python quality checks follow modern Python standards:
# - ruff: Fast linter (replaces flake8, isort, pyupgrade, etc.)
# - black: Opinionated code formatter
# - mypy: Static type checker
# - codespell: Spell checker
# Configuration: pyproject.toml

lint-python: venv ## Run ruff linter on Python code
	@echo "[info] Running ruff linter"
	$(PYTHON) -m ruff check qss/ tests/python/ || true

lint-python-fix: venv ## Run ruff with auto-fix
	@echo "[info] Running ruff with auto-fix"
	$(PYTHON) -m ruff check --fix qss/ tests/python/

format-python: venv ## Format Python code with black
	@echo "[info] Formatting Python code with black"
	$(PYTHON) -m black qss/ tests/python/

format-python-check: venv ## Check Python code formatting (no changes)
	@echo "[info] Checking Python formatting"
	$(PYTHON) -m black --check qss/ tests/python/

typecheck-python: venv ## Run mypy type checker
	@echo "[info] Running mypy type checker"
	$(PYTHON) -m mypy qss/ || true

codespell: venv ## Run spell checker on code and docs
	@echo "[info] Running codespell"
	$(PYTHON) -m codespell qss/ tests/python/ README.md docs/ || true

codespell-fix: venv ## Fix spelling errors automatically
	@echo "[info] Fixing spelling errors"
	$(PYTHON) -m codespell -w qss/ tests/python/ README.md docs/

check-python: lint-python format-python-check typecheck-python codespell ## Run all Python checks
	@echo "[info] All Python quality checks complete"

fix-python: lint-python-fix format-python codespell-fix ## Auto-fix Python issues
	@echo "[info] All Python auto-fixes applied"

# ==============================================================================
## Code Quality - Combined
# ==============================================================================

lint-all: lint-cpp lint-python ## Run all linters (C++ and Python)

format-all: format-cpp format-python ## Format all code

check-all: check-python lint-cpp ## Run all quality checks
	@echo "[info] All quality checks complete"

# ==============================================================================
## Pre-commit Hooks
# ==============================================================================
# Pre-commit runs quality checks before each commit.
# Configuration: .pre-commit-config.yaml

precommit-install: venv ## Install pre-commit hooks
	@echo "[info] Installing pre-commit hooks"
	$(PYTHON) -m pre_commit install

precommit-run: venv ## Run pre-commit on all files (non-fatal)
	@echo "[info] Running pre-commit hooks"
	$(PYTHON) -m pre_commit run --all-files || true

precommit-ci: venv ## Run pre-commit in CI mode (fail on errors)
	@echo "[info] Running pre-commit (CI mode)"
	$(PYTHON) -m pre_commit run --all-files

# ==============================================================================
## Documentation
# ==============================================================================
# Documentation is built with Sphinx (Python) and Doxygen (C++).
# Output: build/docs/html/index.html
# Configuration: docs/conf.py, docs/Doxyfile

docs: venv ## Build all documentation (Sphinx + Doxygen)
	@echo "[info] Building documentation"
	@mkdir -p "$(BUILD_DIR)/docs"
	@if [ -f docs/Doxyfile ]; then \
		doxygen docs/Doxyfile; \
	fi
	$(PYTHON) -m sphinx.cmd.build -b html docs "$(DOCS_HTML)"
	@echo "[info] Documentation built: $(DOCS_HTML)/index.html"

docs-serve: docs ## Serve documentation locally (http://localhost:8000)
	@echo "[info] Serving documentation at http://localhost:8000"
	@cd "$(DOCS_HTML)" && $(PYTHON) -m http.server 8000

docs-clean: ## Clean documentation build artifacts
	@echo "[info] Cleaning documentation"
	@rm -rf "$(DOCS_HTML)" "$(DOXY_OUT)"

# ==============================================================================
# End of Makefile
# ==============================================================================
