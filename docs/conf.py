# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import shutil


def copy_code(old_dir, new_dir):
    if os.path.exists(new_dir):
        shutil.rmtree(new_dir)
    os.makedirs(new_dir)

    folders = ["common", "mbsolve", "results", "setups", "solver"]
    # Copy folders
    for f in folders:
        old_folder = os.path.join(old_dir, f)
        new_folder = os.path.join(new_dir, f)
        shutil.copytree(old_folder, new_folder)
    return


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "monacoQC"
copyright = "2025, TUEICPH"
author = "TUEICPH"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinxcontrib.matlab",
    "sphinxcontrib.katex",
    "myst_parser",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

matlab_original_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)))
matlab_src_dir = os.path.join(os.path.dirname(__file__), "_build/matlab")
copy_code(matlab_original_dir, matlab_src_dir)


primary_domain = "mat"
# autoclass_content = 'class'
matlab_short_links = True
matlab_auto_link = "basic"
autodoc_member_order = "bysource"
matlab_show_property_specs = False
matlab_show_property_default_value = False
matlab_short_links = True
autodoc_default_options = {"private-members": True, "show-inheritance": True}

# Napoleon
napoleon_numpy_docstring = False
napoleon_google_docstring = True
napoleon_include_init_with_doc = True
napoleon_custom_sections = [
    ("Syntax", "params_style"),
    ("Input Arguments", "params_style"),
    ("Output Arguments", "params_style"),
    ("Name Value Arguments", "params_style"),
    ("Properties", "params_style"),
]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
# html_static_path = ["_static"]
