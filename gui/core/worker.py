"""
Worker thread for running ESL-PSC commands.
"""
import subprocess
import sys
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
        """Run the command and emit signals for output/errors."""
        self.is_running = True
        
        try:
            # Create process
            process = subprocess.Popen(
                self.command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                shell=True,
                bufsize=1
            )
            
            # Read output in real-time
            while True:
                output = process.stdout.readline()
                if output == '' and process.poll() is not None:
                    break
                if output:
                    self.signals.output.emit(output.strip())
            
            # Get any remaining output/error
            stdout, stderr = process.communicate()
            
            if stdout:
                self.signals.output.emit(stdout.strip())
            if stderr:
                self.signals.error.emit(stderr.strip())
            
            # Emit finished signal with exit code
            self.signals.finished.emit(process.returncode)
            
        except Exception as e:
            self.signals.error.emit(f"Error running command: {str(e)}")
            self.signals.finished.emit(1)
        finally:
            self.is_running = False
    
    def stop(self):
        """Stop the running process."""
        self.is_running = False
