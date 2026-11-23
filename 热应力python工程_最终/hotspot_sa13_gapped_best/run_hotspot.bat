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
"%HS%" -c sa13_gapped_best.config -f sa13_gapped_best.flp -p sa13_gapped_best.ptrace -materials_file sa13_gapped_best.materials -model_type grid -grid_rows 64 -grid_cols 64 -steady_file outputs/sa13_gapped_best.steady -grid_steady_file outputs/sa13_gapped_best.grid.steady
