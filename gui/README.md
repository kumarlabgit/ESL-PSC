# ESL-PSC GUI Wizard

A graphical user interface for running ESL-PSC analyses with an intuitive wizard interface.

## Features

- Step-by-step wizard for configuring ESL-PSC analyses
- Intuitive file selection
- Parameter configuration with sensible defaults
- Real-time progress tracking
- Cross-platform support (Windows, macOS, Linux)
 - Interactive tree viewer with phenotype coloring:
   - Binary (-1/1) phenotypes color labels red/blue
   - Continuous float phenotypes render via red→blue gradient (low→high)
   - Pair-based coloring overrides phenotype colors; missing phenotypes appear black

## Installation

1. Ensure you have Python 3.8 or later installed
2. Install the required dependencies:
   ```bash
   pip install -r requirements-gui.txt
   ```

## Running the GUI

From the project root directory, run:

```bash
python -m gui.main
```

## Development

### Project Structure

- `gui/main.py` - Application entry point
- `gui/ui/` - User interface components
- `gui/ui/pages/` - Wizard pages
- `gui/ui/widgets/` - Reusable UI widgets
- `gui/core/` - Business logic and data handling
- `gui/resources/` - Icons, styles, and other resources

