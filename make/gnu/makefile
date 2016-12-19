stap : stap.o node.o list.o solver_mode.o datain.o elcal.o assem.o globals.o truss.o memalloc.o 3T.o 4Q.o 9Q.o mathkernel.o hexahedral.o BEAM.o plate.o shell.o plate8Q.o shell8Q.o postprocessor.o
	gfortran -o stap stap.o node.o list.o solver_mode.o datain.o elcal.o assem.o globals.o truss.o memalloc.o 3T.o 4Q.o 9Q.o mathkernel.o hexahedral.o BEAM.o plate.o shell.o plate8Q.o shell8Q.o postprocessor.o

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
	
solver_mode.o : ../../src/solver_mode.f90 ../../src/globals.f90 list.mod node.mod
	gfortran -c ../../src/solver_mode.f90
    
datain.o : ../../src/datain.f90 ../../src/globals.f90 globals.mod memallocate.mod node.mod solver_mode.o
	gfortran -c  ../../src/datain.f90

assem.o : ../../src/assem.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/assem.f90

elcal.o : ../../src/elcal.f90 ../../src/globals.f90 globals.mod memallocate.mod postprocessor.o
	gfortran -c  ../../src/elcal.f90

truss.o : ../../src/truss.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/truss.f90
	
3T.o : ../../src/3T.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/3T.f90
	
4Q.o : ../../src/4Q.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/4Q.f90
	
9Q.o : ../../src/9Q.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/9Q.f90
	
hexahedral.o : ../../src/hexahedral.f90 ../../src/globals.f90 globals.mod memallocate.mod mathkernel.mod
	gfortran -c  ../../src/hexahedral.f90

BEAM.o : ../../src/BEAM.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/BEAM.f90
	
shell.o : ../../src/shell.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/shell.f90
	
plate.o : ../../src/plate.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/plate.f90
	
shell8Q.o : ../../src/shell8Q.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/shell8Q.f90
	
plate8Q.o : ../../src/plate8Q.f90 ../../src/globals.f90 globals.mod memallocate.mod
	gfortran -c  ../../src/plate8Q.f90

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
	
postprocessor.o : ../../src/postprocessor.f90 mathkernel.o memallocate.mod
	gfortran -c  ../../src/postprocessor.f90

clean : 
	rm stap stap.o node.o node_class.mod list.o list_class.mod solver_mode.o datain.o elcal.o assem.o globals.o globals.mod truss.o memalloc.o memallocate.mod 3T.o 4Q.o 9Q.o hexahedral.o BEAM.o plate.o shell.o plate8Q.o shell8Q.o mathkernel.o mathkernel.mod postprocessor.o
