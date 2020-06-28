#!/bin/bash

# parametros generales
test='patchtest-0'
npts='3-3'
E='1000.0'
v='0.3'
plot='False'
show='False'
npmin='9'
npmax='12'

# preproceso
args='--func Preproceso --test '$test' --npts '$npts' --young '$E' --poisson '$v
echo python3 preproceso2d.py $args
python3 ./preproceso2d.py $args

# generacion de nubes
args=$test' '$npts' '$show' '$npmin' '$npmax' false Linux'
echo ./genclouds $args
./genclouds $args

# ejecuta test
input_filename='./DATOS/'$test'_'$npts
echo $input_filename > NUBESPUNT.DAT
#./punto
./punto2

## graficar solucion exacta
#python3 postproceso2d.py --func PLOT_SOL_EXACT --test %test% --young %E% --poisson %v%

## graficar solucion numerica
#python3 postproceso2d.py --func PLOT_SOL_APROX --test $test --npts $npts

## graficar nubes
#python3 postproceso2d.py --func PLOT_CLOUD --test %test% --npts %npts%

