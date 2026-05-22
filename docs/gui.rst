GUI (Streamlit)
===============

The Streamlit application lives in ``quantum-statistics-solver/gui/web_app.py``.

Run it from the ``quantum-statistics-solver/`` directory::

   streamlit run gui/web_app.py

Defaults
--------

Widget defaults (lattice ``N``, hopping ``t``, sweep grids, solver Brent/SA, plot sizes, …)
live in ``gui/defaults.yaml``. Edit that file and rerun the app — no Python changes needed.
Override the path with ``QSS_GUI_DEFAULTS=/path/to/my_defaults.yaml`` if you want several presets.

Overview
--------

The GUI runs a Cartesian grid sweep over:

- band filling :math:`n`
- temperature :math:`k_B T/|t|`
- magnetic field :math:`H_a`

For each chosen distribution (FD, FLB, NFL, NFLB, SSL, HUBBARD) it solves :math:`\mu` and then computes:

- entropy :math:`S`
- energy :math:`E` and free energy :math:`F`
- magnetization :math:`M`
- :math:`n_{kd}`

Panels tab
----------

The **Panels** tab is the main plotting workflow:

- choose **Quantity**
- choose **X axis variable** and **X scale** (linear / log)
- choose **Sweep in panels** (facet variable) and set the remaining variables as **Fixed**
- adjust plot style (panels-per-row, per-panel figure size, optional x/y limits)

The distribution colors and line styles match ``tools/scientific_plot_v2.py``.

Distributions tab
-----------------

The **Distributions** tab provides:

- Notebook-style distribution plots (when the repository layout contains ``tools/scientific_plot_v2.py``).
- A **2D occupancy heatmap** (for 2D lattice) computed using the **solved :math:`\mu`** from the results grid.

Output data
-----------

The GUI writes the chemical-potential grid to::

   quantum-statistics-solver/data/gui_results.csv

and provides download buttons for derived quantities from the active panel selection.
