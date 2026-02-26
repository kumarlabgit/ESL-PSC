@echo off
setlocal

set "SCRIPT_DIR=%~dp0"
set "TOOLKIT_ROOT=%SCRIPT_DIR%.."
set "PYTHONPATH=%TOOLKIT_ROOT%\python;%PYTHONPATH%"
set "PATH=%TOOLKIT_ROOT%\bin;%PATH%"

if defined ESL_PSC_PYTHON (
  "%ESL_PSC_PYTHON%" -m esl_psc_cli.plot_cli %*
  goto :eof
)

where py >nul 2>nul
if %errorlevel%==0 (
  py -3 -m esl_psc_cli.plot_cli %*
) else (
  python -m esl_psc_cli.plot_cli %*
)
