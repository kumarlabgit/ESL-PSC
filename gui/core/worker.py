"""
Worker thread for running ESL-PSC commands.
"""
import contextlib
import os
import io
import re
import threading
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
        # Track deletion-canceler progress counts for step progress bar
        self.del_total_combos: int | None = None  # total combos printed by deletion_canceler
        self.del_current_combo: int = 0           # number of combos completed so far
        # Track current combo index for prettier status messages
        self.current_combo: int | None = None
        self.total_combos: int | None = None
    
    @pyqtSlot()
    def run(self):
        """Execute esl_multimatrix in-process and stream its output."""
        self.is_running = True
        original_cwd = os.getcwd()  # Preserve GUI's working directory
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
                    # Save for later status messages
                    try:
                        self.worker.current_combo = current
                        self.worker.total_combos = total
                    except Exception:
                        pass
                    self.signals.step_status.emit(line.strip().replace("---", "").strip())
                    self.signals.step_progress.emit(0)
                    return True

                # Deletion-canceler overall progress – first record total # combos
                m_del_total = re.search(r"Generated (\d+) species combinations", line)
                if m_del_total:
                    # Save expected total combinations for later progress calc
                    try:
                        self.worker.del_total_combos = int(m_del_total.group(1))
                    except Exception:
                        self.worker.del_total_combos = None
                    # Not a direct progress update we want to show in bar yet
                    return False

                # Deletion-canceler per-combo progress: "Generating alignments for: <species...>"
                if "Generating alignments for" in line:
                    # Increment completed combo counter and compute % if total known
                    if self.worker.del_total_combos:
                        self.worker.del_current_combo += 1
                        pct = int(self.worker.del_current_combo / self.worker.del_total_combos * 100)
                        self.signals.step_progress.emit(pct)
                    # Update status text regardless
                    self.signals.step_status.emit(line.strip())
                    return False

                # ESL preprocess step indicator
                if "Running ESL preprocess" in line or "preprocess_" in line or "preprocess_mac" in line:
                    # Reset deletion-canceler counters so they don't affect later steps
                    self.worker.del_total_combos = None
                    self.worker.del_current_combo = 0
                    # Compose friendly status
                    combo_msg = ""
                    if self.worker.current_combo and self.worker.total_combos:
                        combo_msg = f" for combo {self.worker.current_combo} of {self.worker.total_combos}"
                    friendly = f"Running ESL preprocess{combo_msg}..."
                    # Reset counters for per-file progress
                    self.pre_total_files: int | None = None
                    self.pre_current_file: int = 0
                    self.signals.step_status.emit(friendly)
                    self.signals.output.emit(friendly)
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
                # Capture paths.txt to know file count
                m_paths = re.search(r"paths\.txt", line) and "preprocess" in line
                if m_paths:
                    # Attempt to extract the paths.txt token (last arg ending with paths.txt)
                    try:
                        parts = line.split()
                        for token in reversed(parts):
                            if token.endswith("paths.txt"):
                                paths_file = token
                                break
                        else:
                            paths_file = ""
                        if paths_file and os.path.isfile(paths_file):
                            with open(paths_file, "r", encoding="utf-8", errors="ignore") as pf:
                                self.pre_total_files = sum(1 for _ in pf)
                    except Exception:
                        self.pre_total_files = None
                    return True  # suppress command line output

                # Preprocess per-file progress
                if line.startswith("Processing FASTA file"):
                    if getattr(self, "pre_total_files", None):
                        self.pre_current_file += 1
                        pct = int(self.pre_current_file / self.pre_total_files * 100)
                        self.signals.step_progress.emit(pct)
                        if self.pre_total_files >= 5 and self.pre_current_file % max(1, self.pre_total_files // 10) == 0:
                            self.signals.step_status.emit(f"Preprocessing alignments ({self.pre_current_file}/{self.pre_total_files})")
                    return True  # hide individual file lines

                # Filter noisy full path lines produced by preprocess
                if "-alignments" in line and "/combo_" in line:
                    return True

                status_keywords = ["Building models...", "Calculating predictions and/or weights...", "Running ESL preprocess...", "Generating alignments for"]
                for keyword in status_keywords:
                    if keyword in line:
                        self.signals.step_status.emit(line)
                        self.signals.step_progress.emit(0)
                        return False # Also show this status in the log

                return False # Not a progress line

        try:
            # Create two separate StreamEmitters (Python-level redirection)
            out_stream = StreamEmitter(self, stream_type='stdout')
            err_stream = StreamEmitter(self, stream_type='stderr')

            # ------------------------------------------------------------------
            # Low-level (C / subprocess) redirection: duplicate fd 1 & 2 to pipes
            # ------------------------------------------------------------------
            orig_stdout_fd = os.dup(1)
            orig_stderr_fd = os.dup(2)

            out_r, out_w = os.pipe()
            err_r, err_w = os.pipe()

            # Point process-wide stdout/stderr to the write ends of our pipes
            os.dup2(out_w, 1)
            os.dup2(err_w, 2)

            # Close duplicate write ends – the fds 1/2 now refer to them
            os.close(out_w)
            os.close(err_w)

            def _pipe_reader(fd: int, emitter: StreamEmitter):
                with os.fdopen(fd, "r", buffering=1, errors="replace") as f:
                    for line in f:
                        emitter.write(line)

            # Background threads to read from the pipes
            t_out = threading.Thread(target=_pipe_reader, args=(out_r, out_stream), daemon=True)
            t_err = threading.Thread(target=_pipe_reader, args=(err_r, err_stream), daemon=True)
            t_out.start()
            t_err.start()

            # Redirect *Python-level* stdout/stderr and run the CLI
            with contextlib.redirect_stdout(out_stream), contextlib.redirect_stderr(err_stream):
                if self.is_running:
                    esl_multimatrix.main(self.command_args)

            # Ensure readers flush remaining data
            t_out.join(timeout=0.1)
            t_err.join(timeout=0.1)

        except SystemExit as se:
            exit_code = se.code if isinstance(se.code, int) else 1
        except Exception as e:
            import traceback
            self.signals.error.emit(f"An unexpected worker error occurred: {e}\n{traceback.format_exc()}")
            exit_code = 1
        finally:
            # Restore original file descriptors so that future output is not captured
            try:
                if 'orig_stdout_fd' in locals():
                    os.dup2(orig_stdout_fd, 1)
                if 'orig_stderr_fd' in locals():
                    os.dup2(orig_stderr_fd, 2)
            except Exception:
                pass
            for _fd_name in ('orig_stdout_fd', 'orig_stderr_fd'):
                try:
                    if _fd_name in locals():
                        os.close(locals()[_fd_name])
                except Exception:
                    pass

            self.is_running = False
            # Restore the GUI's original working directory in case esl_multimatrix changed it
            try:
                os.chdir(original_cwd)
            except Exception:
                # Silently ignore – worst case, the old path no longer exists
                pass
            # Only emit finished signal if it wasn't stopped by the user
            if not self.was_stopped:
                self.signals.finished.emit(exit_code)
    
    def stop(self):
        """Flags the worker to stop and emits the finished signal."""
        if self.is_running:
            self.is_running = False
            self.was_stopped = True
            self.signals.finished.emit(-1) # Emit a special code for user stop