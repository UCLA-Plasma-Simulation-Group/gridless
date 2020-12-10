#Makefile for 1D Electrostatic Gridless PIC codes

# Makefile gfortran compiler with MacOS X

#FC90 = gfortran
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fcheck=bounds -fdefault-real-8 -Wall -std=f95

#CCOPTS = -O3 -Wall -std=c99

# Makefile Intel compiler with Mac OS X

#FC90 = ifort
#CC = icc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -std=c99
#LEGACY = -nofor_main

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

all : fpic1gl cpic1gl_f

fpic1gl : fpic1gl.o fpush1gl.o dtimer.o
	$(FC90) $(OPTS90) -o fpic1gl fpic1gl.o fpush1gl.o push1gl_h.o \
        dtimer.o

cpic1gl : cpic1gl.o cpush1gl.o dtimer.o
	$(CC) $(CCOPTS) -o cpic1gl cpic1gl.o cpush1gl.o dtimer.o -lm

fpic1gl_c : fpic1gl_c.o cpush1gl.o dtimer.o
	$(FC90) $(OPTS90) -o fpic1gl_c fpic1gl_c.o cpush1gl.o dtimer.o

cpic1gl_f : cpic1gl.o cpush1gl_f.o fpush1gl.o dtimer.o
	$(FC90) $(OPTS90) $(LEGACY) -o cpic1gl_f cpic1gl.o cpush1gl_f.o fpush1gl.o \
        dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fpush1gl.o : push1gl.f
	$(FC90) $(OPTS90) -o fpush1gl.o -c push1gl.f

push1gl_h.o : push1gl_h.f90
	$(FC90) $(OPTS90) -o push1gl_h.o -c push1gl_h.f90

cpush1gl.o : push1gl.c
	$(CC) $(CCOPTS) -o cpush1gl.o -c push1gl.c

cpush1gl_f.o : push1gl_f.c
	$(CC) $(CCOPTS) -o cpush1gl_f.o -c push1gl_f.c

fpic1gl.o : pic1gl.f90 push1gl_h.o
	$(FC90) $(OPTS90) -o fpic1gl.o -c pic1gl.f90

cpic1gl.o : pic1gl.c
	$(CC) $(CCOPTS) -o cpic1gl.o -c pic1gl.c

fpic1gl_c.o : pic1gl_c.f90
	$(FC90) $(OPTS90) -o fpic1gl_c.o -c pic1gl_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fpic1gl cpic1gl fpic1gl_c cpic1gl_f
