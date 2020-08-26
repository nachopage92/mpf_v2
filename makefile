#Start of the Makefile.
#Defining Variables.
obj_punto = puntofin2d-fl.f90
obj_punto2 = bcg2.f mefpm_le2d.f90
obj_punto3 = bcg2.f mefpm_poisson2d.f90
obj_punto4 = lsqrblas.f lsqr.f mefpm_le2d_2.f90
obj_punto5 = lsqrblas.f lsqr.f mefpm_poisson2d_2.f90
obj_cloud = mrgrnk.f90 NUBES.f90
obj_maxent = lbfgs.f maxent.f90
obj_maxent_v2 = lme_simplificado.f90
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

#all: punto punto2 punto3 genclouds
all: punto punto2 genclouds

punto: $(obj_punto)
	$(f90comp) -o punto $(FLAG) $(obj_punto)

punto2: $(obj_punto2) $(obj_maxent) $(obj_maxent_v2)
	$(f90comp) -o punto2 $(FLAG) $(obj_maxent) $(obj_maxent_v2) $(obj_punto2) $(DEFS)

punto3: $(obj_punto3) $(obj_maxent) $(obj_maxent_v2)
	$(f90comp) -o punto3 $(FLAG) $(obj_maxent) $(obj_maxent_v2) $(obj_punto3) $(DEFS)

punto4: $(obj_punto4) $(obj_maxent) $(obj_maxent_v2)
	$(f90comp) -o punto4 $(FLAG) $(obj_maxent) $(obj_maxent_v2) $(obj_punto4) $(DEFS)

genclouds: $(obj_cloud)
	$(f90comp) -o genclouds $(FLAG) $(obj_cloud)
