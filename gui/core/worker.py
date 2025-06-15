"""
Worker thread for running ESL-PSC commands.
"""
import sys
import shlex
import contextlib
import io
from esl_psc_cli import esl_multimatrix
from PyQt6.QtCore import QObject, pyqtSignal, QRunnable, pyqtSlot, QThreadPool

class WorkerSignals(QObject):
    """Defines the signals available from a running worker thread."""
    output = pyqtSignal(str)
    error = pyqtSignal(str)
    finished = pyqtSignal(int)  # exit code
    progress = pyqtSignal(int)  # progress percentage


class ESLWorker(QRunnable):
    """Worker thread for running ESL-PSC commands."""
    
    def __init__(self, command):
        """Initialize the worker with the command to run."""
        super().__init__()
        self.command = command
        self.signals = WorkerSignals()
        self.is_running = False
    
    @pyqtSlot()
    def run(self):
        """Execute esl_multimatrix in-process and relay its output."""
        self.is_running = True
        try:
            # Split the raw command string exactly as a shell would
            args = shlex.split(self.command)

            stdout_buf = io.StringIO()
            stderr_buf = io.StringIO()
            exit_code = 0

            # Capture everything the script prints
            with contextlib.redirect_stdout(stdout_buf), contextlib.redirect_stderr(stderr_buf):
                try:
                    esl_multimatrix.main(args)
                except SystemExit as se:
                    # esl_multimatrix calls sys.exit() → convert to int
                    exit_code = se.code if isinstance(se.code, int) else 1

            # Stream captured stdout
            for line in stdout_buf.getvalue().splitlines():
                self.signals.output.emit(line)

            # Stream captured stderr
            for line in stderr_buf.getvalue().splitlines():
                self.signals.error.emit(line)

            self.signals.finished.emit(exit_code)

        except Exception as e:
            self.signals.error.emit(f"Error running ESL-PSC: {e}")
            self.signals.finished.emit(1)
        finally:
            self.is_running = False
    
    def stop(self):
        """Stop the running process."""
        self.is_running = False
