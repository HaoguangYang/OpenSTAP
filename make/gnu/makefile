stap : stap.o node.o list.o solver_mode.o datain.o elcal.o assem.o globals.o truss.o memalloc.o 3T.o 4Q.o 9Q.o mathkernel.o 8H.o 27H.o BEAM.o plate.o shell.o plate8Q.o shell8Q.o postprocessor.o modalanlsis.o LANCZOS.o mathops.o ORTHMGS.o jacobi.o
	gfortran -o stap  stap.o node.o list.o solver_mode.o datain.o elcal.o assem.o globals.o truss.o memalloc.o 3T.o 4Q.o 9Q.o mathkernel.o 8H.o 27H.o BEAM.o plate.o shell.o plate8Q.o shell8Q.o postprocessor.o modalanlsis.o LANCZOS.o mathops.o ORTHMGS.o jacobi.o

stap.o : ../../src/stap.f90 ../../src/globals.f90 globals.mod memallocate.mod postprocessor.o
	gfortran -c  ../../src/stap.f90

node.o : ../../src/module/node.f90
	gfortran -c ../../src/module/node.f90
    
node.mod : ../../src/module/node.f90
	gfortran -c ../../src/module/node.f90
    
list.o : ../../src/module/list.f90 ../../src/module/node.f90
	gfortran -c ../../src/module/list.f90
    
list.mod : ../../src/module/list.f90 ../../src/module/node.f90
	gfortran -c ../../src/module/list.f90
	
solver_mode.o : ../../src/module/solver_mode.f90 ../../src/globals.f90 list.mod node.mod
	gfortran -c ../../src/module/solver_mode.f90
    
datain.o : ../../src/datain.f90 ../../src/globals.f90 globals.mod memallocate.mod node.mod solver_mode.o
	gfortran -c  ../../src/datain.f90

assem.o : ../../src/assem.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/assem.f90

elcal.o : ../../src/elcal.f90 ../../src/globals.f90 globals.mod memallocate.mod postprocessor.o
	gfortran -c  ../../src/elcal.f90

truss.o : ../../src/elements/truss.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/elements/truss.f90
	
3T.o : ../../src/elements/3T.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/elements/3T.f90
	
4Q.o : ../../src/elements/4Q.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/elements/4Q.f90
	
9Q.o : ../../src/elements/9Q.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/elements/9Q.f90
	
8H.o : ../../src/elements/8H.f90 ../../src/globals.f90 globals.mod memallocate.mod mathkernel.mod
	gfortran -c  ../../src/elements/8H.f90
	
27H.o : ../../src/elements/8H.f90 ../../src/globals.f90 globals.mod memallocate.mod mathkernel.mod
	gfortran -c  ../../src/elements/27H.f90

BEAM.o : ../../src/elements/BEAM.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/elements/BEAM.f90
	
shell.o : ../../src/elements/shell.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/elements/shell.f90
	
plate.o : ../../src/elements/plate.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/elements/plate.f90
	
shell8Q.o : ../../src/elements/shell8Q.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/elements/shell8Q.f90
	
plate8Q.o : ../../src/elements/plate8Q.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/elements/plate8Q.f90

globals.o : ../../src/globals.f90
	gfortran -c  ../../src/globals.f90 

globals.mod : ../../src/globals.f90 
	gfortran -c  ../../src/globals.f90 

memalloc.o : ../../src/memalloc.f90
	gfortran -c  ../../src/memalloc.f90

memallocate.mod : ../../src/memalloc.f90
	gfortran -c  ../../src/memalloc.f90
	
mathkernel.o : ../../src/module/mathkernel.f90
	gfortran -c  ../../src/module/mathkernel.f90

mathkernel.mod : ../../src/module/mathkernel.f90
	gfortran -c  ../../src/module/mathkernel.f90
	
postprocessor.o : ../../src/postprocessor.f90 mathkernel.o memallocate.mod
	gfortran -c  ../../src/postprocessor.f90
	
modalanlsis.o : ../../src/modalanlsis.f90 globals.mod
	gfortran -c  ../../src/modalanlsis.f90
	
LANCZOS.o : ../../src/modal/LANCZOS.f90 mathops.o ORTHMGS.o jacobi.o
	gfortran -c ../../src/modal/LANCZOS.f90
	
mathops.o : ../../src/modal/mathops.f90
	gfortran -c ../../src/modal/mathops.f90
	
ORTHMGS.o : ../../src/modal/ORTHMGS.f90 mathops.o
	gfortran -c ../../src/modal/ORTHMGS.f90
	
jacobi.o : ../../src/modal/jacobi.f90
	gfortran -c ../../src/modal/jacobi.f90

clean : 
	rm stap stap.o node.o node_class.mod list.o list_class.mod solver_mode.o datain.o elcal.o assem.o globals.o globals.mod truss.o memalloc.o memallocate.mod 3T.o 4Q.o 9Q.o 8H.o 27H.o BEAM.o plate.o shell.o plate8Q.o shell8Q.o mathkernel.o mathkernel.mod postprocessor.o modalanlsis.o LANCZOS.o mathops.o ORTHMGS.o jacobi.o
