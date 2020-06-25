#Start of the Makefile.
#Defining Variables.
obj_punto = puntofin2d-fl.f90
obj_cloud = mrgrnk.f90 NUBES.f90
f90comp = gfortran
FLAG = -ggdb -fcheck=all -fbacktrace -ffpe-trap=zero,invalid,overflow

punto: $(obj_punto)
	$(f90comp) -o punto $(FLAG) $(obj_punto)

genclouds: $(obj_cloud)
	$(f90comp) -o genclouds $(FLAG) $(obj_cloud)
