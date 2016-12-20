subsp : subsp.o LANCZOS.o sub.o ORTHMGS.o jacobi.o
	gfortran -o subsp  subsp.o LANCZOS.o sub.o ORTHMGS.o jacobi.o

subsp.o : ../subsp.f90 LANCZOS.o
	gfortran -c ../subsp.f90
	
LANCZOS.o : ../LANCZOS.f90 sub.o ORTHMGS.o jacobi.o
	gfortran -c ../LANCZOS.f90
	
sub.o : ../sub.f90
	gfortran -c ../sub.f90
	
ORTHMGS.o : ../ORTHMGS.f90 sub.o
	gfortran -c ../ORTHMGS.f90
	
jacobi.o : ../jacobi.f90
	gfortran -c ../jacobi.f90

clean :
	rm fort.16 jacobi.o LANCZOS.o ORTHMGS.o sub.o subsp subsp.o
