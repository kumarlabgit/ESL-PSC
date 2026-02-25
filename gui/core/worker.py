"""
Worker thread for running ESL-PSC commands.
"""
from __future__ import annotations
import os
import io
import re
import threading
import subprocess
import shutil
import sys
from pathlib import Path
from PySide6.QtCore import QObject, Signal, QRunnable, Slot
from contextlib import redirect_stdout, redirect_stderr # noqa: F401
from esl_psc_cli import esl_psc_functions as ecf

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
        self.original_cwd: str | None = None
        # Track deletion-canceler progress counts for step progress bar
        self.del_total_combos: int | None = None  # total combos printed by deletion_canceler
        self.del_current_combo: int = 0           # number of combos completed so far
        # Flag indicating the global deletion-canceler phase has completed
        self.del_cancel_done: bool = False
        # Track current combo index for prettier status messages
        self.current_combo: int | None = None
        self.total_combos: int | None = None
        # Progress phase tracking for per-combo "sub-bar" so that the bar fills
        # only once per combo instead of resetting multiple times.
        # 0 = not started, 1 = deletion-canceler, 2 = preprocess, 3 = build models, 4 = predictions
        self.phase: int = 0
        self.phase_weights: list[int] = [10, 10, 60, 20]  # must sum to 100
        self.total_phases: int = len(self.phase_weights)

        # Keep the alignments directory for preprocessing progress fallback
        self.alignments_dir: str | None = None
        if '--alignments_dir' in self.command_args:
            try:
                idx = self.command_args.index('--alignments_dir') + 1
                if idx < len(self.command_args):
                    self.alignments_dir = self.command_args[idx]
            except Exception:
                self.alignments_dir = None

    @staticmethod
    def _resolve_unified_rust_binary() -> Path | None:
        if os.name == "nt":
            exe_names = [
                "esl-psc.exe",
                "esl-psc",
            ]
        else:
            exe_names = ["esl-psc"]

        candidates: list[Path] = []
        seen: set[str] = set()

        # Source-tree location
        try:
            repo_root = Path(__file__).resolve().parents[2]
            for exe_name in exe_names:
                candidates.append(repo_root / "esl_psc_rs" / "target" / "release" / exe_name)
                candidates.append(repo_root / "bin" / exe_name)
        except Exception:
            pass

        # Packaged app sibling location
        try:
            launcher = Path(os.path.realpath(sys.argv[0]))
            for exe_name in exe_names:
                candidates.append(launcher.with_name(exe_name))
                candidates.append(launcher.parent / exe_name)
                candidates.append(launcher.parent / "bin" / exe_name)
        except Exception:
            pass

        # PATH fallback
        for exe_name in exe_names:
            which = shutil.which(exe_name)
            if which:
                candidates.append(Path(which))

        for cand in candidates:
            key = str(cand)
            if key in seen:
                continue
            seen.add(key)
            if cand.is_file() and os.access(cand, os.X_OK):
                return cand
        return None

    @staticmethod
    def get_command_preview_prefix() -> str:
        rust_bin = ESLWorker._resolve_unified_rust_binary()
        if rust_bin is not None:
            return str(rust_bin)
        return f"{sys.executable} -u -m esl_psc_cli.esl_multimatrix"

    @staticmethod
    def _split_plot_flags(command_args: list[str]) -> tuple[str | None, list[str]]:
        """Return (plot_mode, args_without_plot_flags)."""
        has_cont = "--make_continuous_plot" in command_args or "--make-continuous-plot" in command_args
        has_violin = "--make_sps_plot" in command_args or "--make-sps-plot" in command_args
        has_kde = "--make_sps_kde_plot" in command_args or "--make-sps-kde-plot" in command_args

        mode: str | None = None
        if has_cont:
            mode = "continuous"
        elif has_violin:
            mode = "violin"
        elif has_kde:
            mode = "kde"

        plot_flags = {
            "--make_continuous_plot",
            "--make-continuous-plot",
            "--make_sps_plot",
            "--make-sps-plot",
            "--make_sps_kde_plot",
            "--make-sps-kde-plot",
        }
        filtered = [arg for arg in command_args if arg not in plot_flags]
        return mode, filtered

    @staticmethod
    def _arg_value(command_args: list[str], *names: str, default: str | None = None) -> str | None:
        for i, token in enumerate(command_args):
            if token in names and i + 1 < len(command_args):
                return command_args[i + 1]
        return default

    @staticmethod
    def _arg_pair(command_args: list[str], *names: str) -> tuple[str, str] | None:
        for i, token in enumerate(command_args):
            if token in names and i + 2 < len(command_args):
                return command_args[i + 1], command_args[i + 2]
        return None

    def _run_inprocess_plot(self, mode: str, command_args: list[str]) -> bool:
        """Generate plots using the GUI's Python runtime (no extra bundled runtime)."""
        if "--no_pred_output" in command_args or "--no-pred-output" in command_args:
            self.signals.output.emit("[INFO] Plot requested, but --no_pred_output is set; skipping plot generation.")
            return True

        output_dir = self._arg_value(command_args, "--output_dir", "--output-dir")
        output_base = self._arg_value(
            command_args,
            "--output_file_base_name",
            "--output-file-base-name",
        )
        if not output_dir or not output_base:
            self.signals.error.emit(
                "Unable to generate plots: missing --output_dir or --output_file_base_name."
            )
            return False

        min_genes_raw = self._arg_value(command_args, "--min_genes", "--min-genes", default="0")
        try:
            min_genes = int(min_genes_raw or "0")
        except ValueError:
            min_genes = 0

        pred_csv = Path(output_dir) / f"{output_base}_species_predictions.csv"
        if not pred_csv.is_file():
            self.signals.error.emit(
                f"Unable to generate plots: predictions CSV not found: {pred_csv}"
            )
            return False

        plot_args = [
            "--mode", mode,
            "--pred_csv", str(pred_csv),
            "--title", output_base,
            "--min_genes", str(min_genes),
        ]
        if mode != "continuous":
            pheno_names = self._arg_pair(command_args, "--pheno_names", "--pheno-names")
            if pheno_names:
                plot_args.extend([
                    "--pheno_name1", pheno_names[0],
                    "--pheno_name2", pheno_names[1],
                ])

        try:
            from esl_psc_cli.plot_cli import main as plot_cli_main

            self.signals.output.emit(
                f"[INFO] Generating {mode} prediction plot with bundled Python runtime..."
            )
            rc = int(plot_cli_main(plot_args))
            if rc != 0:
                self.signals.error.emit(
                    f"Plot generation failed with exit code {rc}."
                )
                return False
            self.signals.output.emit("[INFO] Plot generation completed.")
            return True
        except Exception as exc:
            self.signals.error.emit(f"Plot generation failed: {exc}")
            return False
    
    @Slot()
    def run(self):
        """Execute esl_multimatrix in a subprocess and stream its output."""
        self.is_running = True
        self.original_cwd = os.getcwd()  # Preserve GUI's working directory
        exit_code = 0

        class StreamEmitter(io.TextIOBase):
            def __init__(self, worker, stream_type='stdout'):
                super().__init__()
                self.worker = worker
                self.signals = worker.signals
                self._buf = ""
                self.stream_type = stream_type # 'stdout' or 'stderr'

                # Helper to scale step progress across the three phases so the
                # GUI sub-progress bar goes from 0→100 once per combo.
                def _emit_step_progress(pct: int):
                    """
                    Scale the raw percentage of the *current* phase into the
                    overall 0–100 % range for the combo.

                    pct: 0-100 value representing progress within the phase.
                    """
                    phase = max(1, self.worker.phase)
                    # Sum of weights of completed phases
                    completed = sum(self.worker.phase_weights[:phase-1])
                    current_span = self.worker.phase_weights[phase-1]
                    scaled = int(completed + pct * current_span / 100)
                    # Clamp to 100 and detect combo completion
                    end_of_combo = False
                    if pct >= 100 and phase >= self.worker.total_phases:
                        scaled = 100
                        end_of_combo = True
                    self.signals.step_progress.emit(min(scaled, 100))
                    # Update overall progress once the combo has fully completed
                    if end_of_combo and self.worker.current_combo and self.worker.total_combos:
                        self.signals.overall_progress.emit(int(self.worker.current_combo / self.worker.total_combos * 100))

                # Store on self for reuse
                self._emit_step_progress = _emit_step_progress

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
                            # Drop whitespace-only lines during preprocess to reduce noise
                            if not (self.worker.phase == 2 and not line.strip()):
                                self.signals.output.emit(line)
                    else:  # stderr
                        # Always forward stderr lines to preserve spacing
                        self.signals.error.emit(line)
                return len(s)

            def flush(self):
                pass # Writes are handled immediately

            def flush_buffer(self):
                """Process any remaining data in the buffer."""
                if self._buf:
                    if self.stream_type == 'stdout':
                        self.signals.output.emit(self._buf)
                    else:
                        self.signals.error.emit(self._buf)
                    self._buf = ""

            def _parse_progress(self, line):
                # Overall combo progress: "--- Processing combo 1 of 16 (combo_0) ---"
                m_combo = re.search(r"Processing combo (\d+) of (\d+)", line)
                if m_combo:
                    current, total = map(int, m_combo.groups())
                    # Update overall progress immediately to reflect combos already completed
                    if total > 0 and current > 1:
                        self.signals.overall_progress.emit(int((current - 1) / total * 100))
                    # Save for later status messages and end-of-combo updates
                    try:
                        self.worker.current_combo = current
                        self.worker.total_combos = total
                        # Mark deletion-canceler as completed; subsequent "Generating alignments" lines belong to combo-specific work
                        self.worker.del_cancel_done = True
                    except Exception:
                        pass
                    # New combo → reset phase tracker
                    self.worker.phase = 0
                    self.signals.step_status.emit(line.strip().replace("---", "").strip())
                    self._emit_step_progress(0)
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

                # Deletion-canceler per-combo progress: original line lists every species,
                # which overflows the GUI. Replace it with a concise message.
                if "Generating alignments for" in line:
                    # Treat as deletion-canceler progress *only* before the first combo starts
                    if not self.worker.del_cancel_done:
                        self.worker.phase = 1
                        # Increment completed combo counter and compute % if total known
                        if self.worker.del_total_combos:
                            self.worker.del_current_combo += 1
                            pct = int(self.worker.del_current_combo / self.worker.del_total_combos * 100)
                            # During global deletion-canceler phase use the entire bar (0–100)
                            self.signals.step_progress.emit(pct)
                        # Build a concise status message
                        combo_msg = f"combo {self.worker.del_current_combo}"
                        if self.worker.del_total_combos:
                            combo_msg += f" of {self.worker.del_total_combos}"
                        friendly = f"Generating alignments for {combo_msg}..."
                        # Emit to status bar and to log (output). Suppress original long line.
                        self.signals.step_status.emit(friendly)
                        self.signals.output.emit(friendly)
                        return True
                    # After deletion-canceler finished, per-combo alignment generation is noisy – suppress in progress bars and logs
                    return True

                # ESL preprocess step indicator
                if "Running ESL preprocess" in line or "preprocess_" in line or "preprocess_mac" in line:
                    # We’re now in phase 2 (preprocess)
                    self.worker.phase = 2
                    # Reset deletion-canceler counters so they don't affect later steps
                    self.worker.del_total_combos = None
                    self.worker.del_current_combo = 0
                    # Compose friendly status
                    combo_msg = ""
                    if self.worker.current_combo and self.worker.total_combos:
                        combo_msg = f" for combo {self.worker.current_combo} of {self.worker.total_combos}"
                    friendly = f"Running ESL preprocess{combo_msg}..."
                    # Reset counters for per-file progress
                    self.pre_total_files: int | None = 0
                    self.pre_current_file: int = 0
                    if self.worker.alignments_dir and os.path.isdir(self.worker.alignments_dir):
                        try:
                            self.pre_total_files = len([
                                f for f in os.listdir(self.worker.alignments_dir)
                                if ecf.is_fasta(f)
                            ])
                        except Exception:
                            self.pre_total_files = None
                    self.signals.step_status.emit(friendly)
                    self.signals.output.emit(friendly)
                    self._emit_step_progress(0)
                    return True

                # Step progress: "run 123 of 400"
                m_step = re.search(r"run (\d+) of (\d+)", line, re.IGNORECASE)
                if m_step:
                    # If we are still in earlier phases, assume this is model building (phase 3).
                    if self.worker.phase < 3:
                        self.worker.phase = 3
                    current, total = map(int, m_step.groups())
                    if total > 0:
                        self._emit_step_progress(int(current / total * 100))
                    # This is a progress update, but we still want to see it in the log
                    return False

                # Group penalties adjustment status (CLI prints this right before replacement)
                if "adjusting group penalties" in line.lower():
                    # Treat as entering the model-building phase soon; reset phase progress if advancing
                    if self.worker.phase < 3:
                        self.worker.phase = 3
                        self._emit_step_progress(0)
                    self.signals.step_status.emit("adjusting group penalties...")
                    return False

                # Median penalty calculation status (printed before computing median GP)
                if "calculating median group penalty" in line.lower():
                    if self.worker.phase < 3:
                        self.worker.phase = 3
                        self._emit_step_progress(0)
                    self.signals.step_status.emit("calculating median group penalty...")
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
                                if ecf.is_fasta(f)
                            ])
                    except Exception:
                        self.pre_total_files = None
                    if getattr(self, "pre_total_files", None) is None and self.worker.alignments_dir and os.path.isdir(self.worker.alignments_dir):
                        try:
                            self.pre_total_files = len([
                                f for f in os.listdir(self.worker.alignments_dir)
                                if ecf.is_fasta(f)
                            ])
                        except Exception:
                            self.pre_total_files = None
                    if getattr(self, "pre_total_files", None):
                        self.pre_current_file += 1
                        pct = int(self.pre_current_file / self.pre_total_files * 100)
                        self._emit_step_progress(pct)
                        if self.pre_total_files >= 5 and self.pre_current_file % max(1, self.pre_total_files // 10) == 0:
                            self.signals.step_status.emit(f"Preprocessing alignments ({self.pre_current_file}/{self.pre_total_files})")
                    return True  # hide individual file lines

                status_keywords = [
                    "Building models...",
                    "Calculating predictions and/or weights...",
                    "Running ESL preprocess...",
                    "Generating alignments for"
                ]
                for keyword in status_keywords:
                    if keyword in line:
                        # Determine phase implied by keyword
                        new_phase = self.worker.phase
                        if "Generating alignments" in keyword:
                            new_phase = 1
                        elif "Running ESL preprocess" in keyword:
                            new_phase = 2
                        elif "Building models" in keyword:
                            new_phase = 3
                        elif "Calculating predictions" in keyword:
                            new_phase = 4

                        # If we have moved to a *new* phase, update tracker and reset sub-progress
                        if new_phase > self.worker.phase:
                            self.worker.phase = new_phase
                            # Reset within new phase
                            self._emit_step_progress(0)
                        # Otherwise (same phase), keep progress where it is.
                        self.signals.step_status.emit(line)
                        return False  # Also show this status in the log

                return False # Not a progress line

        def _run_subprocess(command: list[str]) -> int:
            out_stream = StreamEmitter(self, stream_type='stdout')
            err_stream = StreamEmitter(self, stream_type='stderr')

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
            return self.process.returncode

        # Prefer unified Rust runner when available.
        rust_bin = self._resolve_unified_rust_binary()
        # Detect a packaged run either by sys.frozen or by the launcher not ending in '.py'
        packaged = getattr(sys, 'frozen', False) or Path(sys.argv[0]).suffix.lower() != ".py"

        if rust_bin is not None:
            try:
                plot_mode, rust_args = self._split_plot_flags(self.command_args)
                command = [str(rust_bin), *rust_args]
                self.signals.output.emit(f"[INFO] Running unified Rust CLI: {rust_bin}")
                if plot_mode is not None:
                    self.signals.output.emit(
                        "[INFO] Plot flags will be handled in-process by the GUI Python runtime."
                    )
                exit_code = _run_subprocess(command)
                if exit_code == 0 and plot_mode is not None:
                    if not self._run_inprocess_plot(plot_mode, rust_args):
                        exit_code = 1
            except Exception as e:
                import traceback
                if not self.was_stopped:
                    self.signals.error.emit(
                        f"An unexpected worker error occurred: {e}\n{traceback.format_exc()}"
                    )
                    exit_code = 1
                else:
                    exit_code = -1
            finally:
                self.is_running = False
                if self.process and self.process.poll() is None:
                    self.process.kill()
                self.process = None
                try:
                    os.chdir(self.original_cwd)
                except Exception:
                    pass
                if not self.was_stopped:
                    self.signals.finished.emit(exit_code)
        elif packaged:
            # Fallback for packaged builds when Rust CLI isn't present.
            try:
                out_stream = StreamEmitter(self, stream_type='stdout')
                err_stream = StreamEmitter(self, stream_type='stderr')
                # Import only in fallback path to avoid eagerly pulling the full
                # Python CLI stack into packaged apps that run Rust by default.
                from esl_psc_cli.esl_multimatrix import main as esl_main

                with redirect_stdout(out_stream), redirect_stderr(err_stream):
                    esl_main(self.command_args)

                # Flush any remaining output
                out_stream.flush_buffer()
                err_stream.flush_buffer()
                exit_code = 0

            except Exception as e:
                import traceback
                if not self.was_stopped:
                    self.signals.error.emit(
                        f"An unexpected worker error occurred: {e}\n{traceback.format_exc()}"
                    )
                    exit_code = 1
                else:
                    # Treat exceptions raised after a user-initiated stop as cancellation
                    exit_code = -1
            finally:
                self.is_running = False
                try:
                    os.chdir(self.original_cwd)
                except Exception:
                    pass
                if not self.was_stopped:
                    self.signals.finished.emit(exit_code)
        else:
            # Source fallback when Rust binary is unavailable.
            try:
                command = [
                    sys.executable, "-u", "-m",
                    "esl_psc_cli.esl_multimatrix",
                    *self.command_args,
                ]
                self.signals.output.emit("[INFO] Unified Rust CLI not found; falling back to Python CLI.")
                exit_code = _run_subprocess(command)
            except Exception as e:
                import traceback
                if not self.was_stopped:
                    self.signals.error.emit(
                        f"An unexpected worker error occurred: {e}\n{traceback.format_exc()}"
                    )
                    exit_code = 1
                else:
                    # Suppress errors that arise due to user-initiated stop and mark as cancelled
                    exit_code = -1
            finally:
                self.is_running = False
                if self.process and self.process.poll() is None:
                    self.process.kill()
                self.process = None
                try:
                    os.chdir(self.original_cwd)
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
            # Attempt to restore the original working directory immediately on stop
            try:
                if self.original_cwd and os.getcwd() != self.original_cwd:
                    os.chdir(self.original_cwd)
            except Exception:
                pass
            # Emit special code for user stop
            self.signals.finished.emit(-1)
