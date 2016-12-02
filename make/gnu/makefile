stap : stap.o datain.o elcal.o assem.o globals.o truss.o memalloc.o 4Q.o 9Q.o mathkernel.o hexahedral.o BEAM.o
	gfortran -o stap stap.o datain.o elcal.o assem.o globals.o truss.o memalloc.o 4Q.o 9Q.o mathkernel.o hexahedral.o BEAM.o

stap.o : ../../src/stap.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/stap.f90

datain.o : ../../src/datain.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/datain.f90

assem.o : ../../src/assem.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/assem.f90

elcal.o : ../../src/elcal.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/elcal.f90

truss.o : ../../src/truss.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/truss.f90
	
4Q.o : ../../src/4Q.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/4Q.f90
	
9Q.o : ../../src/9Q.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/9Q.f90
	
hexahedral.o : ../../src/hexahedral.f90 ../../src/globals.f90 globals.mod memallocate.mod mathkernel.mod
	gfortran -c  ../../src/hexahedral.f90

BEAM.o : ../../src/BEAM.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/BEAM.f90

globals.o : ../../src/globals.f90
	gfortran -c  ../../src/globals.f90 

globals.mod : ../../src/globals.f90 
	gfortran -c  ../../src/globals.f90 

memalloc.o : ../../src/memalloc.f90
	gfortran -c  ../../src/memalloc.f90

memallocate.mod : ../../src/memalloc.f90
	gfortran -c  ../../src/memalloc.f90
	
mathkernel.o : ../../src/mathkernel.f90
	gfortran -c  ../../src/mathkernel.f90

mathkernel.mod : ../../src/mathkernel.f90
	gfortran -c  ../../src/mathkernel.f90

clean : 
	rm stap stap.o datain.o elcal.o assem.o globals.o globals.mod truss.o memalloc.o memallocate.mod 4Q.o 9Q.o hexahedral.o BEAM.o mathkernel.o mathkernel.mod LOAD.TMP ELMNT.TMP
