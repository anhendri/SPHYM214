OPTS   = -fdollar-ok
LAPACK = -llapack -lblas
PGPLOT = -I%GWMODDIR% -Wl,--subsystem,console -lpgplot -lGrWin -mwindows
OPENMP = -fopenmp
MPI    = -I C:\MPICH2\include -L C:\MPICH2\lib -lmpi -lfmpich2g -lmpicxx
DEBUG  = -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace 
PROFIL = -pg

all:STRUCTURE.exe

STRUCTURE.exe:STRUCTURE.o

	gfortran -o STRUCTURE.exe *.o $(LAPACK)


STRUCTURE.o: 

	gfortran STRUCTURE.f90 -c
	
clean:

	del *.o

