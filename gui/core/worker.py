"""
Worker thread for running ESL-PSC commands.
"""
import os
import io
import re
import threading
import subprocess
import sys
from pathlib import Path
from PySide6.QtCore import QObject, Signal, QRunnable, Slot

class WorkerSignals(QObject):
    """Defines the signals available from a running worker thread."""
    output = Signal(str)
    error = Signal(str)
    finished = Signal(int)  # exit code

    overall_progress = Signal(int) # Overall progress (0-100)
    step_progress = Signal(int)    # Progress of current step (0-100)
    step_status = Signal(str)      # Text description of current step


class ESLWorker(QRunnable):
    """Worker thread for running ESL-PSC commands."""
    
    def __init__(self, command_args):
        """Initialize the worker with the command to run."""
        super().__init__()
        self.command_args = command_args
        self.signals = WorkerSignals()
        self.is_running = False
        self.was_stopped = False # Flag to indicate if stop() was called
        self.process = None
        # Track deletion-canceler progress counts for step progress bar
        self.del_total_combos: int | None = None  # total combos printed by deletion_canceler
        self.del_current_combo: int = 0           # number of combos completed so far
        # Track current combo index for prettier status messages
        self.current_combo: int | None = None
        self.total_combos: int | None = None
        # Keep the alignments directory for preprocessing progress fallback
        self.alignments_dir: str | None = None
        if '--alignments_dir' in self.command_args:
            try:
                idx = self.command_args.index('--alignments_dir') + 1
                if idx < len(self.command_args):
                    self.alignments_dir = self.command_args[idx]
            except Exception:
                self.alignments_dir = None
    
    @Slot()
    def run(self):
        """Execute esl_multimatrix in a subprocess and stream its output."""
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
                    if self.worker.alignments_dir and os.path.isdir(self.worker.alignments_dir):
                        try:
                            self.pre_total_files = len([
                                f for f in os.listdir(self.worker.alignments_dir)
                                if f.endswith('.fas')
                            ])
                        except Exception:
                            self.pre_total_files = None
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
                        elif self.worker.alignments_dir and os.path.isdir(self.worker.alignments_dir):
                            self.pre_total_files = len([
                                f for f in os.listdir(self.worker.alignments_dir)
                                if f.endswith('.fas')
                            ])
                    except Exception:
                        self.pre_total_files = None
                    return True  # suppress command line output

                # Preprocess per-file progress
                if line.startswith("Processing FASTA file"):
                    if self.pre_total_files is None and self.worker.alignments_dir and os.path.isdir(self.worker.alignments_dir):
                        try:
                            self.pre_total_files = len([
                                f for f in os.listdir(self.worker.alignments_dir)
                                if f.endswith('.fas')
                            ])
                        except Exception:
                            self.pre_total_files = None
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
            out_stream = StreamEmitter(self, stream_type='stdout')
            err_stream = StreamEmitter(self, stream_type='stderr')

            def _build_command() -> list[str]:
                """
                • If we’re running from source (argv[0] ends with .py), use
                    the current Python to launch the module with -m.
                • Otherwise we’re inside the packaged bundle: call the helper
                  binary `esl_multimatrix(.exe)` that lives next to the GUI
                  launcher (`main` on macOS, `ESL-PSC.exe` on Windows).
                """
                launcher_path = Path(os.path.realpath(sys.argv[0]))
                running_from_source = launcher_path.suffix.lower() == ".py"

                if running_from_source:
                    return [
                        sys.executable, "-u", "-m",
                        "esl_psc_cli.esl_multimatrix",
                        *self.command_args,
                    ]

                # -------- packaged path --------
                exe_name = "esl_multimatrix" + (".exe" if os.name == "nt" else "")
                cli_helper = launcher_path.with_name(exe_name)

                # Pass the bundle’s MacOS/ (or Windows dir) so the helper
                # can find bin/preprocess and bin/sg_lasso
                bundle_dir = launcher_path.parent  # Contents/MacOS or the exe dir

                return [
                    str(cli_helper),
                    "--esl_main_dir", str(bundle_dir),
                    *self.command_args,
                ]

            command = _build_command()
            self.process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,
            )

            def _reader(pipe, emitter: StreamEmitter):
                for line in iter(pipe.readline, ''):
                    if not self.is_running:
                        break
                    emitter.write(line)
                pipe.close()

            t_out = threading.Thread(target=_reader, args=(self.process.stdout, out_stream), daemon=True)
            t_err = threading.Thread(target=_reader, args=(self.process.stderr, err_stream), daemon=True)
            t_out.start()
            t_err.start()

            self.process.wait()
            t_out.join()
            t_err.join()
            exit_code = self.process.returncode

        except Exception as e:
            import traceback
            self.signals.error.emit(
                f"An unexpected worker error occurred: {e}\n{traceback.format_exc()}"
            )
            exit_code = 1
        finally:
            self.is_running = False
            if self.process and self.process.poll() is None:
                self.process.kill()
            self.process = None
            try:
                os.chdir(original_cwd)
            except Exception:
                pass
            if not self.was_stopped:
                self.signals.finished.emit(exit_code)
    
    def stop(self):
        """Flags the worker to stop and emits the finished signal."""
        if self.is_running:
            self.is_running = False
            self.was_stopped = True
            if self.process and self.process.poll() is None:
                try:
                    self.process.kill()
                except Exception:
                    pass
            self.signals.finished.emit(-1) # Emit a special code for user stop