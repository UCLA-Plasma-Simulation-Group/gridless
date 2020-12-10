#Makefile for 2D Electrostatic Gridless PIC codes

# Makefile gfortran compiler with MacOS X

#FC90 = gfortran
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fcheck=bounds -fdefault-real-8 -Wall -std=f95

#CCOPTS = -O3 -Wall -std=c99

#LEGACY =

# Makefile Intel compiler with Mac OS X

#FC90 = ifort
#CC = icc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -std=c99

# Makefile Intel compiler with Linux

#FC90 = ifort
#CC = icc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -std=c99

#LEGACY = -nofor_main

# Makefile gfortran compiler with Linux

FC90 = gfortran
CC = gcc

#OPTS90 = -O3
OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -Wall -std=f95

CCOPTS = -O3 -Wall -std=c99

# Makefile PGI compiler with Linux

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

# Makefile Cray compiler with Linux

#FC90 = ftn
#CC = cc

#OPTS90 = -O 3
#OPTS90 = -O 3 -s real64
#OPTS90 = -O 3 -R b -s real64 -en

#CCOPTS = -O 3 -h c99 -h conform

#

# Linkage rules

all : fpic2gl cpic2gl_f

fpic2gl : fpic2gl.o fpush2gl.o dtimer.o
	$(FC90) $(OPTS90) -o fpic2gl fpic2gl.o fpush2gl.o push2gl_h.o \
        dtimer.o

cpic2gl : cpic2gl.o cpush2gl.o dtimer.o
	$(CC) $(CCOPTS) -o cpic2gl cpic2gl.o cpush2gl.o dtimer.o -lm

fpic2gl_c : fpic2gl_c.o cpush2gl.o dtimer.o
	$(FC90) $(OPTS90) -o fpic2gl_c fpic2gl_c.o cpush2gl.o dtimer.o

cpic2gl_f : cpic2gl.o cpush2gl_f.o fpush2gl.o dtimer.o
	$(FC90) $(OPTS90) $(LEGACY) -o cpic2gl_f cpic2gl.o cpush2gl_f.o fpush2gl.o \
        dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fpush2gl.o : push2gl.f
	$(FC90) $(OPTS90) -o fpush2gl.o -c push2gl.f

push2gl_h.o : push2gl_h.f90
	$(FC90) $(OPTS90) -o push2gl_h.o -c push2gl_h.f90

cpush2gl.o : push2gl.c
	$(CC) $(CCOPTS) -o cpush2gl.o -c push2gl.c

cpush2gl_f.o : push2gl_f.c
	$(CC) $(CCOPTS) -o cpush2gl_f.o -c push2gl_f.c

fpic2gl.o : pic2gl.f90 push2gl_h.o
	$(FC90) $(OPTS90) -o fpic2gl.o -c pic2gl.f90

cpic2gl.o : pic2gl.c
	$(CC) $(CCOPTS) -o cpic2gl.o -c pic2gl.c

fpic2gl_c.o : pic2gl_c.f90
	$(FC90) $(OPTS90) -o fpic2gl_c.o -c pic2gl_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fpic2gl cpic2gl fpic2gl_c cpic2gl_f
