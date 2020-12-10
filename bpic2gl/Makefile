#Makefile for 2-1/2D Electromagnetic Gridless PIC codes

# Makefile gfortran compiler with MacOS X

#FC90 = gfortran
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fcheck=bounds -fdefault-real-8 -Wall -std=f95

#CCOPTS = -O3 -Wall -std=c99

# Makefile Intel compiler with Mac OS X

#FC90 = ifort
#CC = gcc

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

all : fbpic2gl cbpic2gl_f

fbpic2gl : fbpic2gl.o fbpush2gl.o dtimer.o
	$(FC90) $(OPTS90) -o fbpic2gl fbpic2gl.o fbpush2gl.o bpush2gl_h.o \
        dtimer.o

cbpic2gl : cbpic2gl.o cbpush2gl.o dtimer.o
	$(CC) $(CCOPTS) -o cbpic2gl cbpic2gl.o cbpush2gl.o dtimer.o -lm

fbpic2gl_c : fbpic2gl_c.o cbpush2gl.o dtimer.o
	$(FC90) $(OPTS90) -o fbpic2gl_c fbpic2gl_c.o cbpush2gl.o dtimer.o

cbpic2gl_f : cbpic2gl.o cbpush2gl_f.o fbpush2gl.o dtimer.o
	$(FC90) $(OPTS90) $(LEGACY) -o cbpic2gl_f cbpic2gl.o cbpush2gl_f.o \
	    fbpush2gl.o dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fbpush2gl.o : bpush2gl.f
	$(FC90) $(OPTS90) -o fbpush2gl.o -c bpush2gl.f

bpush2gl_h.o : bpush2gl_h.f90
	$(FC90) $(OPTS90) -o bpush2gl_h.o -c bpush2gl_h.f90

cbpush2gl.o : bpush2gl.c
	$(CC) $(CCOPTS) -o cbpush2gl.o -c bpush2gl.c

cbpush2gl_f.o : bpush2gl_f.c
	$(CC) $(CCOPTS) -o cbpush2gl_f.o -c bpush2gl_f.c

fbpic2gl.o : bpic2gl.f90 bpush2gl_h.o
	$(FC90) $(OPTS90) -o fbpic2gl.o -c bpic2gl.f90

cbpic2gl.o : bpic2gl.c
	$(CC) $(CCOPTS) -o cbpic2gl.o -c bpic2gl.c

fbpic2gl_c.o : bpic2gl_c.f90
	$(FC90) $(OPTS90) -o fbpic2gl_c.o -c bpic2gl_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fbpic2gl cbpic2gl fbpic2gl_c cbpic2gl_f
