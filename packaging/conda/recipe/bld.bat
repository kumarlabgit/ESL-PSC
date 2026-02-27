@echo on
setlocal

cargo build --release --manifest-path esl_psc_rs\Cargo.toml
if errorlevel 1 exit /b 1

if not exist "%LIBRARY_BIN%" mkdir "%LIBRARY_BIN%"
set "BIN_PATH=esl_psc_rs\target\release\esl-psc.exe"
if not exist "%BIN_PATH%" set "BIN_PATH=target\release\esl-psc.exe"
if not exist "%BIN_PATH%" (
  echo could not find built binary esl-psc.exe in expected target directories
  exit /b 1
)
copy /Y "%BIN_PATH%" "%LIBRARY_BIN%\esl-psc.exe"
if errorlevel 1 exit /b 1

if not exist "%SP_DIR%\esl_psc_cli" mkdir "%SP_DIR%\esl_psc_cli"
xcopy /E /I /Y esl_psc_cli "%SP_DIR%\esl_psc_cli"
if errorlevel 1 exit /b 1

if not exist "%SP_DIR%\gui" mkdir "%SP_DIR%\gui"
if not exist "%SP_DIR%\gui\core" mkdir "%SP_DIR%\gui\core"
copy /Y gui\__init__.py "%SP_DIR%\gui\__init__.py"
if errorlevel 1 exit /b 1
copy /Y gui\core\fast_scan.py "%SP_DIR%\gui\core\fast_scan.py"
if errorlevel 1 exit /b 1
copy /Y gui\core\fasta_io.py "%SP_DIR%\gui\core\fasta_io.py"
if errorlevel 1 exit /b 1
copy /Y gui\core\ancestral_reconstruction.py "%SP_DIR%\gui\core\ancestral_reconstruction.py"
if errorlevel 1 exit /b 1
