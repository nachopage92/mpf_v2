#Start of the Makefile.
#Defining Variables.
obj_punto = puntofin2d-fl.f90
obj_punto2 = mefpm.f90
obj_cloud = mrgrnk.f90 NUBES.f90
obj_maxent = priorweightfunction.f90 lbfgs.f maxent.f90
f90comp = gfortran
FLAG = -cpp -ggdb -fcheck=all -fbacktrace -ffpe-trap=zero,invalid,overflow
# if blas lapack is installed, then use:
#DEFS = -DUSELAPBLAS

# FLAGS PARA GFORTRAN Y GDB (recomendado utilizar zero, invalid y overflow)
#	-‘invalid’ (invalid floating point operation, such as SQRT(-1.0))
#	-‘zero’ (division by zero)
#	-‘overflow’ (overflow in a floating point operation)
#	-‘underflow’ (underflow in a floating point operation)
#	-‘inexact’ (loss of precision during operation)
#	-‘denormal’ (operation performed on a denormal value)

all: punto punto2 genclouds

punto: $(obj_punto)
	$(f90comp) -o punto $(FLAG) $(obj_punto)

punto2: $(obj_punto2)
	$(f90comp) -o punto2 $(FLAG) $(obj_punto2) $(obj_maxent) $(DEFS)

genclouds: $(obj_cloud)
	$(f90comp) -o genclouds $(FLAG) $(obj_cloud)
