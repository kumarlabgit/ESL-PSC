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

Before submitting a pull request, please run the following checks to ensure code quality and functionality. It is assumed the necessary development dependencies (`pytest`, `pytest-qt`, `flake8`) are already installed in the environment.

### Step 3.1: Code Linting (Check for Errors)

This project uses `flake8` to check for logical errors in the code, such as undefined names or unused imports. The following command is configured to **ignore purely stylistic warnings** (like line length or whitespace) and focus only on potential bugs.

```bash
# Run the linter to find logical code errors (F-codes) and ignore stylistic warnings (E/W-codes)
flake8 . --ignore=E,W
```

### Step 3.2: Automated Testing (Check for Functionality)

The project contains smoke tests to verify that both the CLI and GUI components can run without fatal errors.

* `test_cli_smoke.py`: Runs the command-line tool with demo data to ensure it produces the correct output files.
* `test_gui_smoke.py`: Instantiates the GUI application in an "offscreen" mode to ensure the UI can be built without crashing.
Run all tests using the pytest command:

```bash
# Run all tests in the 'tests/' directory
pytest
```
All tests must pass.