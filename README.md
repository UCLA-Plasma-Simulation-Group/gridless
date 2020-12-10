# Skeleton Gridless Codes

Skeleton codes are bare-bones but fully functional PIC codes containing all the crucial elements but not the diagnostics and initial conditions typical of production codes.  The only diagnostic are the energy values.  The procedures are written in Fortran, but there are wrapper functions so enable to main code to run from either Fortran or C.

# Serial gridless codes

The basic serial codes do not make use of any parallelism, and are the base codes for students or researchers who are unfamiliar with girdless particle codes.

## Electrostatic:

1.  1D Serial Electrostatic Spectral code:  pic1gl
2.  2D Serial Electrostatic Spectral code:  pic2gl

## Electromagnetic:

3.  1-2/2D Serial Electromagnetic Spectral code:  bpic1gl
4.  2-1/2D Serial Electromagnetic Spectral code:  bpic2gl


# Parallel gridless codes

These codes illustrate how to use shared memory parallelism with OpenMP. The push procedures are parallelized over particles, the deposit procedures are parallelized over Fourier modes.

## Electrostatic:

1. 1D Parallel Electrostatic Spectral code:  mpic1gl
2. 2D Parallel Electrostatic Spectral code:  mpic2gl

## Electromagnetic:

3. 1-2/2D Parallel Electromagnetic Spectral code:  mbpic1gl
4. 2-1/2D Parallel Electromagnetic Spectral code:  mbpic2gl
