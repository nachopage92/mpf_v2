#!/bin/bash

# parametros generales

all_test='patchtest-0 patchtest-1 patchtest-2 infinite-plate'
all_test='infinite-plate'
all_npts='3-3 5-5 11-11 21-21 31-31 41-41 51-51'

#all_test='cantilever'
#all_npts='17-3 25-4 33-5 41-6 49-7 57-8 65-9 72-10 81-11'


E='1000.0'
v='0.3'
#plot='True'
plot='False'
#show='True'
show='False'
npmin='9'
npmax='12'


# :::::::::::::::::::::::::

preproceso () {
for test in $all_test; do
    for npts in $all_npts; do    
        # preproceso
        args='--func Preproceso --test '$test' --npts '$npts' --young '$E' --poisson '$v' --plot '$plot' --show '$show
        echo python3 preproceso2d.py $args
        python3 ./preproceso2d.py $args
    done
done
}

run_gencld () {
for test in $all_test; do
    for npts in $all_npts; do    
        # generacion de nubes
        args=$test' '$npts' '$show' '$npmin' '$npmax' false Linux'
        echo ./genclouds $args
        ./genclouds $args
    done
done
}
    
run_test2d () {
for test in $all_test; do
    for npts in $all_npts; do    
        # ejecuta test
        input_filename='./DATOS/'$test'_'$npts
        echo $input_filename > NUBESPUNT.DAT
        #./punto
        ./punto2
    done
done
}
    
plot_appsol () {
for test in $all_test; do
    for npts in $all_npts; do    
        # graficar solucion numerica
        args='--func PLOT_SOL_APROX --test '$test' --npts '$npts
        echo python3 postproceso2d.py $args
        python3 postproceso2d.py $args
    done
done
}

plot_clouds () {
for test in $all_test; do
    for npts in $all_npts; do    
        # graficar nubes
        args='--func PLOT_CLOUD --test '$test' --npts '$npts
        echo python3 postproceso2d.py $args
        python3 postproceso2d.py $args
    done
done
}

plot_exactsol () {
for test in $all_test; do
    # graficar solucion exacta
    args='--func PLOT_SOL_EXACT --test '$test' --young '$E' --poisson '$v
    echo python3 postproceso2d.py $args
    python3 postproceso2d.py $args
done
}


run_all () {
for test in $all_test; do
    for npts in $all_npts; do    

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
        
        if [ "$plot" == 'True' ] ; then
            # graficar solucion numerica
            python3 postproceso2d.py --func PLOT_SOL_APROX --test $test --npts $npts
        fi
        
    done
done
}