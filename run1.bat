:: parametros generales
set test=patchtest-0
set E=1000.0
set v=0.3
set plot=False
set show=False
set npmin=9
set npmax=12

::  preproceso
python3 preproceso2d.py --func Preproceso --test %test% --npts 3-3 --young %E% --poisson %v% 
python3 preproceso2d.py --func Preproceso --test %test% --npts 5-5 --young %E% --poisson %v% 
python3 preproceso2d.py --func Preproceso --test %test% --npts 11-11 --young %E% --poisson %v% 
python3 preproceso2d.py --func Preproceso --test %test% --npts 51-51 --young %E% --poisson %v% 
::python3 preproceso2d.py --func Preproceso --test %test% --npts 91-91 --young %E% --poisson %v% 


:: generacion de nubes
genclouds.exe %test% 3-3 %show% %npmin% %npmax% false Windows
genclouds.exe %test% 5-5 %show% %npmin% %npmax% false Windows
genclouds.exe %test% 11-11 %show% %npmin% %npmax% false Windows
genclouds.exe %test% 51-51 %show% %npmin% %npmax% false Windows
::genclouds.exe %test% 91-91 %show% %npmin% %npmax% false Windows

:: ejecutar test
echo ./DATOS/%test%_3-3 > NUBESPUNT.DAT
punto.exe
echo ./DATOS/%test%_5-5 > NUBESPUNT.DAT
punto.exe
echo ./DATOS/%test%_11-11 > NUBESPUNT.DAT
punto.exe
echo ./DATOS/%test%_51-51 > NUBESPUNT.DAT
punto.exe
::echo ./DATOS/%test%_91-91 > NUBESPUNT.DAT
::punto.exe

:: graficar solucion exacta
thon3 postproceso2d.py --func PLOT_SOL_EXACT --test %test% --young %E% --poisson %v% 

:: graficar solucion numerica
python3 postproceso2d.py --func PLOT_SOL_APROX --test %test% --npts 3-3
python3 postproceso2d.py --func PLOT_SOL_APROX --test %test% --npts 5-5
python3 postproceso2d.py --func PLOT_SOL_APROX --test %test% --npts 11-11
python3 postproceso2d.py --func PLOT_SOL_APROX --test %test% --npts 51-51
::python3 postproceso2d.py --func PLOT_SOL_APROX --test %test% --npts 91-91
