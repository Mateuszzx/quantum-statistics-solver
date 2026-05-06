from pathlib import Path
import sys

project = "qss"
author = "qss authors"

# Try to import the built extension from .venv or build dir
root = Path(__file__).resolve().parents[1]
venv_site = next((root.glob(".venv/lib/python*/site-packages")), None)
if venv_site:
    sys.path.insert(0, str(venv_site))
build_py = root / "build" / "Release" / "python"
if build_py.is_dir():
    sys.path.insert(0, str(build_py))

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.ifconfig",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
    "myst_parser",
    "breathe",
]

autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}
napoleon_numpy_docstring = True
napoleon_google_docstring = True

html_theme = "furo"

# Breathe setup (optional if Doxygen XML is present)
doxy_xml = root / "build" / "doxygen" / "xml"
breathe_projects = {}
if doxy_xml.is_dir():
    breathe_projects["qss"] = str(doxy_xml)
breathe_default_project = "qss" if "qss" in breathe_projects else ""

# Keep build noise low
nitpicky = False
