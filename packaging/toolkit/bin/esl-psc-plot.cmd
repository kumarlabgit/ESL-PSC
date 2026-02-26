@echo off
setlocal

set "SCRIPT_DIR=%~dp0"
set "TOOLKIT_ROOT=%SCRIPT_DIR%.."
"%TOOLKIT_ROOT%\bin\esl-psc.exe" plot %*
