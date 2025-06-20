## 1. Project Overview

This repository contains the **Evolutionary Sparse Learning with Paired Species Contrast (ESL-PSC)** toolkit. The project consists of two main parts:
1.  A core **Command-Line Interface (CLI)** for performing scientific analysis.
2.  A **Graphical User Interface (GUI)** built with PyQt6 that acts as a user-friendly wrapper around the CLI.

The primary goal is to analyze molecular convergence signatures from multiple sequence alignments.

## 2. Directory Structure

The repository is organized as follows:
.
├── .github/          # Configuration for GitHub Actions (CI).
├── bin/              # Pre-compiled helper binaries (OS-specific). 
├── esl_psc_cli/      # Source code for the core command-line tool.
│   ├── init.py
│   ├── esl_multimatrix.py # Main logic for the multi-matrix CLI entry point.
│   ├── esl_integrator.py  # Core integration and model-building logic.
│   ├── esl_psc_functions.py # Helper functions, including subprocess calls.
│   └── ...
├── gui/              # Source code for the PyQt6 GUI application.
│   ├── init.py
│   ├── main.py       # Main GUI entry point.
│   ├── core/         # Core GUI logic (e.g., worker threads).
│   └── ui/           # UI components (windows, pages, widgets).
├── images/           # Images used in the README.
├── photosynthesis_alignments/ # Demo data: alignment files.
├── tests/            # Automated tests for the project.
│   └── test_cli_smoke.py  # A high-level smoke test for the CLI.
├── .gitignore
├── AGENTS.MD         # This instruction file.
├── README.md         # General documentation for human users.
├── esl_multimatrix.py # A stub for the main CLI entry point
└── requirements-gui.txt # Python dependencies for the GUI.

## 3. Verification Workflow

Before submitting a pull request, please run the following checks to ensure code quality and functionality. Ensure the development dependencies are installed:

```bash
# 1. Core + GUI dependencies
pip install -r requirements-gui.txt

# 2. Dev-only tools
pip install pytest pytest-qt flake8
```


### Step 3.1: Code Linting – Logical Errors Only

We lint only for `F`-codes (potential bugs).  The research helpers in `additional_code/` are *not* part of the main library and are excluded.

```bash
# Lint for F-codes, skip stylistic issues and the research helpers
flake8 --select=F --exclude additional_code .
```

### Step 3.2: Automated Testing (Check for Functionality)

The project contains smoke tests to verify that both the CLI and GUI components can run without fatal errors.

* `test_cli_smoke.py`: Runs the command-line tool with demo data to ensure it produces the correct output files.
* `test_gui_smoke.py`: Instantiates the GUI application in an "offscreen" mode to ensure the UI can be built without crashing.
Run all tests:

```bash
pytest -q
```

The `tests/conftest.py` file automatically adds the project root to `sys.path` and provides a compatibility fixture so **no environment variables need to be set**.

All tests must pass.

