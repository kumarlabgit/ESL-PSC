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
        """Execute esl_multimatrix in‑process, stream its output, and update progress."""
        self.is_running = True
        exit_code = 0

        class StreamEmitter(io.TextIOBase):
            """A minimal text IO object that forwards every complete line via Qt signals."""
            def __init__(self, out_sig, err_sig, prog_sig):
                super().__init__()
                self.out_sig = out_sig
                self.err_sig = err_sig
                self.prog_sig = prog_sig
                self._buf = ""

            # stdout handler ----------------------------------------------------
            def write(self, s):
                self._buf += s
                while "\n" in self._buf:
                    line, self._buf = self._buf.split("\n", 1)
                    self.out_sig.emit(line)
                    self._parse_progress(line)
                return len(s)

            def flush(self):
                if self._buf:
                    self.out_sig.emit(self._buf)
                    self._parse_progress(self._buf)
                    self._buf = ""

            # stderr handler ----------------------------------------------------
            def _parse_progress(self, line):
                m = re.search(r'run\s+(\d+)\s+of\s+(\d+)', line, re.IGNORECASE)
                if m:
                    cur, total = map(int, m.groups())
                    if total > 0:
                        pct = int(cur / total * 100)
                        self.prog_sig.emit(pct)

        try:
            args = self.command

            out_stream = StreamEmitter(self.signals.output,
                                    self.signals.error,
                                    self.signals.progress)
            err_stream = StreamEmitter(self.signals.error,
                                    self.signals.error,
                                    self.signals.progress)

            # Redirect stdout/stderr so every write() is pushed immediately
            with contextlib.redirect_stdout(out_stream), contextlib.redirect_stderr(err_stream):
                try:
                    esl_multimatrix.main(args)
                except SystemExit as se:
                    exit_code = se.code if isinstance(se.code, int) else 1

            if exit_code == 0:
                self.signals.progress.emit(100)  # force 100 % on clean exit

        except Exception as e:
            self.signals.error.emit(f"Error running ESL-PSC: {e}")
            exit_code = 1
        finally:
            self.signals.finished.emit(exit_code)
            self.is_running = False
    
    def stop(self):
        """Stop the running process."""
        self.is_running = False
