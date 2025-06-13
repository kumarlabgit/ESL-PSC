# ESL-PSC GUI Wizard

A graphical user interface for running ESL-PSC analyses with an intuitive wizard interface.

## Features

- Step-by-step wizard for configuring ESL-PSC analyses
- Intuitive file selection with drag-and-drop support
- Parameter configuration with sensible defaults
- Real-time progress tracking
- Cross-platform support (Windows, macOS, Linux)

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
  - `main_window.py` - Main window and wizard implementation
- `gui/core/` - Business logic and data handling
- `gui/resources/` - Icons, styles, and other resources

