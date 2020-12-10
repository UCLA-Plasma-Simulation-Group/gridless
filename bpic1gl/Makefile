#Makefile for 1-2/2D Electromagnetic Gridless PIC codes

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

# Makefile Cray compiler with Linux

#FC90 = ftn
#CC = cc

#OPTS90 = -O 3
#OPTS90 = -O 3 -s real64
#OPTS90 = -O 3 -R b -s real64 -en

#CCOPTS = -O 3 -h c99 -h conform

#

# Linkage rules

all : fbpic1gl cbpic1gl_f

fbpic1gl : fbpic1gl.o fbpush1gl.o dtimer.o
	$(FC90) $(OPTS90) -o fbpic1gl fbpic1gl.o fbpush1gl.o bpush1gl_h.o \
        dtimer.o

cbpic1gl : cbpic1gl.o cbpush1gl.o dtimer.o
	$(CC) $(CCOPTS) -o cbpic1gl cbpic1gl.o cbpush1gl.o dtimer.o -lm

fbpic1gl_c : fbpic1gl_c.o cbpush1gl.o dtimer.o
	$(FC90) $(OPTS90) -o fbpic1gl_c fbpic1gl_c.o cbpush1gl.o dtimer.o

cbpic1gl_f : cbpic1gl.o cbpush1gl_f.o fbpush1gl.o dtimer.o
	$(FC90) $(OPTS90) $(LEGACY) -o cbpic1gl_f cbpic1gl.o cbpush1gl_f.o \
	fbpush1gl.o dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fbpush1gl.o : bpush1gl.f
	$(FC90) $(OPTS90) -o fbpush1gl.o -c bpush1gl.f

bpush1gl_h.o : bpush1gl_h.f90
	$(FC90) $(OPTS90) -o bpush1gl_h.o -c bpush1gl_h.f90

cbpush1gl.o : bpush1gl.c
	$(CC) $(CCOPTS) -o cbpush1gl.o -c bpush1gl.c

cbpush1gl_f.o : bpush1gl_f.c
	$(CC) $(CCOPTS) -o cbpush1gl_f.o -c bpush1gl_f.c

fbpic1gl.o : bpic1gl.f90 bpush1gl_h.o
	$(FC90) $(OPTS90) -o fbpic1gl.o -c bpic1gl.f90

cbpic1gl.o : bpic1gl.c
	$(CC) $(CCOPTS) -o cbpic1gl.o -c bpic1gl.c

fbpic1gl_c.o : bpic1gl_c.f90
	$(FC90) $(OPTS90) -o fbpic1gl_c.o -c bpic1gl_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fbpic1gl cbpic1gl fbpic1gl_c cbpic1gl_f
