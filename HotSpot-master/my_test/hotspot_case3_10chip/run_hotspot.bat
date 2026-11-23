@echo off
set ROOT=%~dp0..\HotSpot-master\HotSpot-master
if not exist outputs mkdir outputs
if exist "%ROOT%\hotspot.exe" (
  set HS=%ROOT%\hotspot.exe
) else if exist "%ROOT%\hotspot" (
  set HS=%ROOT%\hotspot
) else (
  echo [error] compile HotSpot first: cd HotSpot-master\HotSpot-master && make
  exit /b 1
)
"%HS%" -c case3_10chip.config -f case3_10chip.flp -p case3_10chip.ptrace -materials_file case3_10chip.materials -model_type grid -grid_rows 64 -grid_cols 64 -steady_file outputs/case3_10chip.steady -grid_steady_file outputs/case3_10chip.grid.steady
