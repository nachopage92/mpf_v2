:: parametros generales
set test=patchtest-0
set npts=91-91
set E=1000.0
set v=0.3
set plot=False
set show=False
set npmin=9
set npmax=12

:: preproceso
python3 preproceso2d.py --func Preproceso --test %test% --npts %npts% --young %E% --poisson %v%

:: generacion de nubes
genclouds.exe %test% %npts% %show% %npmin% %npmax% false Windows

:: ejecuta test
echo ./DATOS/%test%_%npts% > NUBESPUNT.DAT
punto.exe

:::: graficar solucion exacta
::python3 postproceso2d.py --func PLOT_SOL_EXACT --test %test% --young %E% --poisson %v%

:: graficar solucion numerica
python3 postproceso2d.py --func PLOT_SOL_APROX --test %test% --npts %npts%

:: graficar nubes
python3 postproceso2d.py --func PLOT_CLOUD --test %test% --npts %npts%
