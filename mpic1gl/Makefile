#Makefile for 1D Electrostatic Gridless OpenMP PIC codes

# Makefile gfortran compiler with MacOS X

#OpenMP
#MPFC = gfortran -fopenmp
#MPCC = gcc -fopenmp
#NoOpenMP
#MPFC = gfortran
#MPCC = gcc

#FC90 = gfortran
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fcheck=bounds -fdefault-real-8 -Wall -std=f95

#CCOPTS = -O3 -Wall -std=c99

# Makefile Intel compiler with Mac OS X

# OpenMP
#MPFC = ifort -qopenmp
#MPCC = gcc -fopenmp
#NoOpenMP
#MPFC = ifort
#MPCC = gcc

#FC90 = ifort
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -std=c99
#LEGACY = -nofor_main

# Makefile Intel compiler with Linux

# OpenMP
#MPFC = ifort -qopenmp
#MPCC = icc -qopenmp
#NoOpenMP
#MPFC = ifort
#MPCC = icc

#FC90 = ifort
#CC = icc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -std=c99
#LEGACY = -nofor_main

# Makefile gfortran compiler with Linux

# OpenMP
MPFC = gfortran -fopenmp
MPCC = gcc -fopenmp
#NoOpenMP
#MPFC = gfortran
#MPCC = gcc

FC90 = gfortran
CC = gcc

#OPTS90 = -O3
OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -Wall -std=f95

CCOPTS = -O3 -Wall -std=c99

# Makefile PGI compiler with Linux

# OpenMP
#MPFC = pgf90
#MPCC = gcc -fopenmp

#FC90 = pgf90
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard

#CCOPTS = -O3 -Wall -std=c99

# Makefile Nag compiler with Linux

#FC90 = nagfor
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -default_kinds:r=64
#OPTS90 = -O3 -default_kinds:r=64 -C=array -nan -w=all -f95

#CCOPTS = -O3 -Wall -std=c99

# Makefile Nag compiler with Linux

# OpenMP
#MPFC = nagfor -openmp
#MPCC = gcc -fopenmp

#FC90 = nagfor
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -default_kinds:r=64
#OPTS90 = -O3 -default_kinds:r=64 -C=array -nan -w=all -f95

#CCOPTS = -O3 -Wall -std=c99

# Makefile Cray compiler with Linux

# OpenMP
#MPFC = ftn -h omp
#MPCC = cc -h omp
#NoOpenMP
#MPFC = ftn -h noomp
#MPCC = cc -h noomp

#FC90 = ftn
#CC = cc

#OPTS90 = -O 3
#OPTS90 = -O 3 -s real64
#OPTS90 = -O 3 -R b -s real64 -en

#CCOPTS = -O 3 -h c99 -h conform

#

# Linkage rules

all : fmpic1gl cmpic1gl_f

fmpic1gl : fmpic1gl.o fmpush1gl.o fomplib.o dtimer.o
	$(MPFC) $(OPTS90) -o fmpic1gl fmpic1gl.o fmpush1gl.o fomplib.o \
	mpush1gl_h.o omplib_h.o dtimer.o

cmpic1gl : cmpic1gl.o cmpush1gl.o complib.o dtimer.o
	$(MPCC) $(CCOPTS) -o cmpic1gl cmpic1gl.o cmpush1gl.o  complib.o \
	dtimer.o -lm

fmpic1gl_c : fmpic1gl_c.o cmpush1gl.o complib.o dtimer.o
	$(MPFC) $(OPTS90) -o fmpic1gl_c fmpic1gl_c.o cmpush1gl.o complib.o \
	dtimer.o

cmpic1gl_f : cmpic1gl.o cmpush1gl_f.o complib_f.o fmpush1gl.o fomplib.o \
              dtimer.o
	$(MPFC) $(OPTS90) $(LEGACY) -o cmpic1gl_f cmpic1gl.o cmpush1gl_f.o \
	complib_f.o fmpush1gl.o fomplib.o dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fomplib.o : omplib.f
	$(MPFC) $(OPTS90) -o fomplib.o -c omplib.f

omplib_h.o : omplib_h.f90
	$(FC90) $(OPTS90) -o omplib_h.o -c omplib_h.f90

complib.o : omplib.c
	$(MPCC) $(CCOPTS) -o complib.o -c omplib.c

complib_f.o : omplib_f.c
	$(MPCC) $(CCOPTS) -o complib_f.o -c omplib_f.c

fmpush1gl.o : mpush1gl.f
	$(MPFC) $(OPTS90) -o fmpush1gl.o -c mpush1gl.f

mpush1gl_h.o : mpush1gl_h.f90
	$(FC90) $(OPTS90) -o mpush1gl_h.o -c mpush1gl_h.f90

cmpush1gl.o : mpush1gl.c
	$(CC) $(CCOPTS) -o cmpush1gl.o -c mpush1gl.c

cmpush1gl_f.o : mpush1gl_f.c
	$(MPCC) $(CCOPTS) -o cmpush1gl_f.o -c mpush1gl_f.c

fmpic1gl.o : mpic1gl.f90 mpush1gl_h.o omplib_h.o
	$(FC90) $(OPTS90) -o fmpic1gl.o -c mpic1gl.f90

cmpic1gl.o : mpic1gl.c
	$(CC) $(CCOPTS) -o cmpic1gl.o -c mpic1gl.c

fmpic1gl_c.o : mpic1gl_c.f90
	$(FC90) $(OPTS90) -o fmpic1gl_c.o -c mpic1gl_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fmpic1gl cmpic1gl fmpic1gl_c cmpic1gl_f
