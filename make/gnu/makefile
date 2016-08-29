stap : stap.o datain.o elcal.o assem.o  globals.o truss.o memalloc.o
	gfortran -o stap stap.o datain.o elcal.o assem.o  globals.o truss.o memalloc.o

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

globals.o : ../../src/globals.f90
	gfortran -c  ../../src/globals.f90 

globals.mod : ../../src/globals.f90 
	gfortran -c  ../../src/globals.f90 

memalloc.o : ../../src/memalloc.f90
	gfortran -c  ../../src/memalloc.f90

memallocate.mod : ../../src/memalloc.f90
	gfortran -c  ../../src/memalloc.f90

clean : 
	rm stap stap.o datain.o elcal.o assem.o globals.o globals.mod truss.o memalloc.o memallocate.mod
