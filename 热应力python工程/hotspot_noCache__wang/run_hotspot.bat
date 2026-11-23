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
"%HS%" -c cpu_gpu_cim_noCache__wang.config -f cpu_gpu_cim_noCache__wang.flp -p cpu_gpu_cim_noCache__wang.ptrace -materials_file cpu_gpu_cim_noCache__wang.materials -model_type grid -grid_rows 64 -grid_cols 64 -steady_file outputs/cpu_gpu_cim_noCache__wang.steady -grid_steady_file outputs/cpu_gpu_cim_noCache__wang.grid.steady
