GUI (Streamlit)
===============

If you use QSS (library or GUI) in published work, please cite
Kuterba, Wójcik \& Spałek, arXiv:2611.26489 — see ``CITATION.cff``,
``docs/references.bib``, and the *Citing* section of the project README.

The Streamlit application entry point is ``gui/web_app.py``. Views and widgets live in
separate modules under ``gui/`` (see :ref:`gui-layout` below).

Run
---

From the ``quantum-statistics-solver/`` directory::

   streamlit run gui/web_app.py

Or, if you use the project Makefile::

   make gui-launch

Defaults
--------

Widget defaults (lattice ``N``, hopping ``t`` / ``t'``, sweep grids, solver Brent/SA,
plot sizes, distribution state point, …) live in ``gui/defaults.yaml``. Edit that file
and rerun the app — no Python changes needed. Override the path with
``QSS_GUI_DEFAULTS=/path/to/my_defaults.yaml`` for alternate presets.

After changing YAML, use a fresh browser session or clear Streamlit cache if widgets
still show old values (Streamlit keeps per-widget state in ``session_state``).

Overview
--------

The sidebar switches between two views: **Panels** and **Distributions**.

**Panels** runs a Cartesian grid sweep over:

- band filling :math:`n`
- temperature :math:`k_B T/|t|`
- magnetic field :math:`H_a`
- on-site interaction :math:`U/|t|`

For each chosen distribution (FD, FLB, NFL, NFLB, SSL, HUBBARD) it solves
:math:`\mu` and then computes derived columns registered in ``gui/derived_quantities.py``:

- chemical potential :math:`\mu`
- entropy :math:`S`
- heat capacity :math:`C`
- energy :math:`E` and free energy :math:`F`
- magnetization :math:`M`
- :math:`n_{kd}`

To add another post-processed quantity, implement the function in ``qss/`` and register
one ``DerivedQuantitySpec`` in ``derived_quantities.py`` (no change to ``web_app.py``).

Panels view
-----------

Workflow:

1. Configure **Lattice**, **Distributions**, and **Experiment parameters** in the
   sidebar (each sweep variable is in a collapsed expander by default).
2. Optionally enable **Prescribe μ** to skip the μ(n) root solve and sweep μ
   directly (band filling *n* is then computed from occupancy). Otherwise open
   **Solver parameters** (collapsed expander): Brent or simulated annealing.
3. Click **Run experiment**.
4. Choose **Quantity**, **X axis variable**, **X scale** (linear / log), and **Sweep in
   panels**; fix the remaining sweep variables (includes μ when prescribed).
5. Adjust **Plot style** (figure size, panels per row, optional axis limits).

Panel line colors and styles follow the same palette as ``tools/scientific_plot_v2.py``
(see ``gui/panel_plot.py``).

Distributions view
----------------

Standalone workflow — does **not** use Panels experiment results:

- Set a single state point :math:`(n, k_B T, H_a, U)` in the sidebar.
- **Occupancy map / 1D cut**: solves :math:`\mu` from :math:`n` for the selected
  distribution, then plots a 2D map or a 1D cut (horizontal / diagonal).
- **1D BZ cuts** (2D lattice): uses ``gui/scientific_plot.py`` (``ScientificPlotter``)
  to overlay several :math:`n` and :math:`k_B T` values; :math:`\mu` is solved per curve.

Output data
-----------

The Panels experiment writes the chemical-potential grid to::

   quantum-statistics-solver/data/gui_results.csv

Logs go under ``logs/``. The **Data (table + download)** expander exports the
currently selected quantity as CSV.

.. _gui-layout:

Package layout
--------------

+---------------------------+------------------------------------------+
| Module                    | Role                                     |
+===========================+==========================================+
| ``web_app.py``            | Entry point, view switch                 |
| ``view_panels.py``        | Experiment run and sweep panels          |
| ``view_distributions.py`` | Occupancy maps and 1D BZ cuts            |
| ``defaults.yaml`` /       | Editable widget defaults                 |
| ``defaults.py``           |                                          |
| ``derived_quantities.py`` | Registry of plottable quantities         |
| ``sweep_input.py``        | Sweep UI (scalar / list / lin / log)     |
| ``solver_controls.py``    | Brent / SA sidebar widgets               |
| ``sidebar_controls.py``   | Lattice and distribution multiselect     |
| ``panel_plot.py``         | Multi-panel matplotlib plots             |
| ``scientific_plot.py``    | 1D BZ cut plots (Distributions)        |
| ``occupancy_plot.py``     | 2D occupancy / cut plots                 |
+---------------------------+------------------------------------------+
