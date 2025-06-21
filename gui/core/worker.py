"""
Worker thread for running ESL-PSC commands.
"""
import contextlib
import io
import re
from esl_psc_cli import esl_multimatrix
from PyQt6.QtCore import QObject, pyqtSignal, QRunnable, pyqtSlot

class WorkerSignals(QObject):
    """Defines the signals available from a running worker thread."""
    output = pyqtSignal(str)
    error = pyqtSignal(str)
    finished = pyqtSignal(int)  # exit code

    overall_progress = pyqtSignal(int) # Overall progress (0-100)
    step_progress = pyqtSignal(int)    # Progress of current step (0-100)
    step_status = pyqtSignal(str)      # Text description of current step


class ESLWorker(QRunnable):
    """Worker thread for running ESL-PSC commands."""
    
    def __init__(self, command_args):
        """Initialize the worker with the command to run."""
        super().__init__()
        self.command_args = command_args
        self.signals = WorkerSignals()
        self.is_running = False
        self.was_stopped = False # Flag to indicate if stop() was called
    
    @pyqtSlot()
    def run(self):
        """Execute esl_multimatrix in-process and stream its output."""
        self.is_running = True
        exit_code = 0

        class StreamEmitter(io.TextIOBase):
            def __init__(self, worker, stream_type='stdout'):
                super().__init__()
                self.worker = worker
                self.signals = worker.signals
                self._buf = ""
                self.stream_type = stream_type # 'stdout' or 'stderr'

            def write(self, s):
                if not self.worker.is_running: return 0
                # Normalize carriage returns that may come from the CLI
                self._buf += s.replace("\r", "\n")
                while "\n" in self._buf:
                    line, self._buf = self._buf.split("\n", 1)
                    if self.stream_type == 'stdout':
                        # Only parse stdout for progress updates
                        is_progress = self._parse_progress(line)
                        if not is_progress:
                            # Emit the line even if it's empty to preserve original spacing
                            self.signals.output.emit(line)
                    else:  # stderr
                        # Always forward stderr lines to preserve spacing
                        self.signals.error.emit(line)
                return len(s)

            def flush(self):
                pass # Writes are handled immediately

            def _parse_progress(self, line):
                # Overall combo progress: "--- Processing combo 1 of 16 (combo_0) ---"
                m_combo = re.search(r"Processing combo (\d+) of (\d+)", line)
                if m_combo:
                    current, total = map(int, m_combo.groups())
                    if total > 0:
                        self.signals.overall_progress.emit(int(current / total * 100))
                    self.signals.step_status.emit(line.strip().replace("---", "").strip())
                    self.signals.step_progress.emit(0)
                    return True

                # Step progress: "run 123 of 400"
                m_step = re.search(r"run (\d+) of (\d+)", line, re.IGNORECASE)
                if m_step:
                    current, total = map(int, m_step.groups())
                    if total > 0:
                        self.signals.step_progress.emit(int(current / total * 100))
                    # This is a progress update, but we still want to see it in the log
                    return False

                # Step status for major phases
                status_keywords = ["Building models...", "Calculating predictions and/or weights...", "Running ESL preprocess...", "Generating alignments for"]
                for keyword in status_keywords:
                    if keyword in line:
                        self.signals.step_status.emit(line)
                        self.signals.step_progress.emit(0)
                        return False # Also show this status in the log

                return False # Not a progress line

        try:
            # Create two separate streams
            out_stream = StreamEmitter(self, stream_type='stdout')
            err_stream = StreamEmitter(self, stream_type='stderr')

            # Redirect streams and run the main function
            with contextlib.redirect_stdout(out_stream), contextlib.redirect_stderr(err_stream):
                if self.is_running:
                    esl_multimatrix.main(self.command_args)

        except SystemExit as se:
            exit_code = se.code if isinstance(se.code, int) else 1
        except Exception as e:
            import traceback
            self.signals.error.emit(f"An unexpected worker error occurred: {e}\n{traceback.format_exc()}")
            exit_code = 1
        finally:
            self.is_running = False
            # Only emit finished signal if it wasn't stopped by the user
            if not self.was_stopped:
                self.signals.finished.emit(exit_code)
    
    def stop(self):
        """Flags the worker to stop and emits the finished signal."""
        if self.is_running:
            self.is_running = False
            self.was_stopped = True
            self.signals.finished.emit(-1) # Emit a special code for user stop