@echo off
setlocal
cd /d "%~dp0"

where pythonw >nul 2>nul
if errorlevel 1 (
    python thicken_hex_mesh_gui.pyw
) else (
    start "" pythonw thicken_hex_mesh_gui.pyw
)
