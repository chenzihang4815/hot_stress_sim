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
"%HS%" -c lbfgs_L2__basic.config -f lbfgs_L2__basic.flp -p lbfgs_L2__basic.ptrace -materials_file lbfgs_L2__basic.materials -model_type grid -grid_rows 64 -grid_cols 64 -steady_file outputs/lbfgs_L2__basic.steady -grid_steady_file outputs/lbfgs_L2__basic.grid.steady
