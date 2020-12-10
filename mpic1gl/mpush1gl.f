c Fortran Library for Skeleton 1D Electrostatic Gridless OpenMP PIC Code
c written by Viktor K. Decyk, UCLA
c-----------------------------------------------------------------------
      subroutine DISTR1(part,vtx,vdx,npx,idimp,nop,nx,ipbc)
c for 1d code, this subroutine calculates initial particle co-ordinate
c and velocity, with uniform density and maxwellian velocity with drift
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c vtx = thermal velocity of particles in x direction
c vdx = drift velocity of particles x direction
c npx = number of particles distributed in x direction
c idimp = size of phase space = 2
c nop = number of particles
c nx = system length in x direction
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer npx, idimp, nop, nx, ipbc
      real part, vtx, vdx
      dimension part(idimp,nop)
c local data
      integer j
      real edgelx, at1, sum1
      double precision dsum1
      double precision ranorm
c set boundary values
      edgelx = 0.0
      at1 = real(nx)/real(npx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         at1 = real(nx-2)/real(npx)
      endif
c uniform density profile
      do 10 j = 1, npx
      part(1,j) = edgelx + at1*(real(j) - .5)
   10 continue
c maxwellian velocity distribution
      do 20 j = 1, npx
      part(2,j) = vtx*ranorm()
   20 continue
c add correct drift
      dsum1 = 0.0d0
      do 30 j = 1, npx
      dsum1 = dsum1 + part(2,j)
   30 continue
      sum1 = dsum1
      sum1 = sum1/real(npx) - vdx
      do 40 j = 1, npx
      part(2,j) = part(2,j) - sum1
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPUSH1GL(part,fx,qbm,dt,ek,idimp,nop,nx,nxhm,nxvh)
c for 1d code, this subroutine updates particle co-ordinate and
c velocity using leap-frog scheme in time using gridless spectral
c version, with periodic boundaries
c OpenMP version, parallelization over particles
c (16 flops + tan + div)*nxhm plus 8 flops/particle
c 2*nxhm loads, 2 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t))*dt, where q/m is charge/mass,
c and x(t+dt) = x(t) + vx(t+dt/2)*dt
c fx(x(t)) is calculated from the expression
c fx(x) = sum(fx(n)*exp(sqrt(-1)*2*n*pi*x/nx))
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c fx(n) = x component of force/charge at fourier mode (n-1)
c if nxhm=nx/2, nx/2 modes are included, in a format compatible with the
c real to complex FFT used in traditional PIC
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2)
c idimp = size of phase space = 2
c nop = number of particles
c nx = system length in x direction
c nxhm = number of fourier modes kept
c nxvh = dimension of field array, must be >= nxhm
      implicit none
      integer nop, idimp, nx, nxhm, nxvh
      real part, qbm, dt, ek
      complex fx
      dimension part(idimp,nop), fx(nxvh)
c local data
      integer i, j
      real zero, anx, dnx, qtm, tn, at1, at2, dkx, x, vx, dx
      double precision sum1, ex
      complex zt1
      zero = 0.0
      anx = real(nx)
      dnx = 6.28318530717959/anx
      qtm = qbm*dt
      sum1 = 0.0d0
c loop over particles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,x,vx,at1,dkx,tn,at2,dx,ex,zt1) REDUCTION(+:sum1)
      do 20 i = 1, nop
      x = part(1,i)
      vx = part(2,i)
c find electric field
      ex = 0.0d0
c mode numbers 0 < kx < nx/2
      at1 = dnx*x
      do 10 j = 2, nxhm
      dkx = at1*real(j-1)
c zt1 = cmplx(cos(dkx),sin(dkx))
      tn = tan(0.5*dkx)
      at2 = (tn + tn)/(1 + tn*tn)
      zt1 = cmplx(1.0-tn*at2,at2)
      ex = ex + real(fx(j)*zt1)
   10 continue
      ex = 2.0d0*ex + real(fx(1))
c special case to match conventional PIC code
      if (nxhm.eq.(nx/2)) then
         dkx = at1*real(nxhm)
         dx = ex + aimag(fx(1))*cos(dkx)
      else
         dx = ex
      endif
c new velocity
      dx = vx + qtm*dx
c average kinetic energy
      sum1 = sum1 + (dx + vx)**2
      part(2,i) = dx
c new position
      dx = x + dx*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
c set new position
      part(1,i) = dx
   20 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine MDPOST1GL(part,q,qm,nop,idimp,nx,nxhm,nxvh)
c for 1d code, this subroutine calculates particle charge density
c using gridless spectral version, periodic boundaries
c OpenMP version, parallelization over fourier modes
c (11 flops + tan + div)*nxhm plus 1 flop/particle
c nxhm loads and stores
c input: all, output: q
c charge density is calculated from the expression:
c q(n) = sum(qm*exp(-sqrt(-1)*2*n*pi*x/nx))
c part(1,n) = position x of particle n
c q(n) = charge density at fourier grid point n
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nx = system length in x direction
c nxhm = number of fourier modes kept
c nxvh = dimension of charge array, must be >= nxhm
      implicit none
      integer nop, idimp, nx, nxhm, nxvh
      real part, qm
      complex q
      dimension part(idimp,nop), q(nxvh)
c local data
      integer i, j
      real dnx, qmn, at1, dkx, tn, at2
      dnx = 6.28318530717959/real(nx)
      qmn = qm/real(nx)
c find fourier components
      do 20 i = 1, nop
c mode numbers 0 < kx < nxhm
      at1 = dnx*part(1,i)
!$OMP PARALLEL DO PRIVATE(j,dkx,tn,at2)
      do 10 j = 2, nxhm
      dkx = at1*real(j-1)
c cmplx(cos(dkx),-sin(dkx))
      tn = tan(0.5*dkx)
      at2 = (tn + tn)/(1 + tn*tn)
      q(j) = q(j) + qmn*cmplx(1.0-tn*at2,-at2)
   10 continue
!$OMP END PARALLEL DO
c special case to match conventional PIC code
      if (nxhm.eq.(nx/2)) then
         dkx = at1*real(nxhm)
         q(1) = cmplx(real(q(1))+qmn,aimag(q(1))+qmn*cos(dkx))
      else
         q(1) = cmplx(real(q(1))+qmn,0.0)
      endif
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FFC1INITGL(ffc,ax,affp,nx,nxhm)
c calculates tables needed by a one dimensional gridless field solvers.
c Truncated fourier series gives a sinc function instead of a delta 
c function for particle.  The smoothing function, the convolution of a
c box car function with a gaussian, integrates the sinc shape to obtain
c a smoother shape with minimized oscillations.
c input: ax,affp,nx,nxhm, output: ffc
c g(kx) = (affp/kx**2)*s(kx)
c s(kx) = sin(kx*w)/kx*w)*exp(-(kx*ax)**2/2)
c where kx = 2pi*j/nx, w = 1.9269/(dnx*real(nxhm)),
c Si(1.9269) = pi/2 and j-1 = fourier mode number
c aimag(ffc(j)) = finite-size particle shape factor s
c real(ffc(j)) = potential green's function g
c ax = half-width of particle in x direction
c affp = normalization constant = nx/np, where np=number of particles
c nx = system length in x direction
c nxhm = number of fourier modes kept in x
      implicit none
      integer nx, nxhm
      real ax, affp
      complex ffc
      dimension ffc(nxhm)
c local data
      integer j
      real dnx, dkx, w, at1, at2
      dnx = 6.28318530717959/real(nx)
      w = 1.9269/(dnx*real(nxhm))
c prepare form factor array
      do 10 j = 2, nxhm
      dkx = dnx*real(j - 1)
      at1 = dkx*dkx
      at2 = dkx*w
      at2 = (sin(at2)/at2)*exp(-.5*(dkx*ax)**2)
      ffc(j) = cmplx(affp*at2/at1,at2)
   10 continue
      ffc(1) = cmplx(affp,1.0)
      return
      end
c-----------------------------------------------------------------------
      subroutine POIS1GL(q,fx,isign,ffc,ax,affp,we,nx,nxhm,nxvh)
c this subroutine solves 1d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions, for gridless spectral code
c for isign = 0, input: isign,ax,affp,nx,nxhm,nxvh output: ffc
c for isign /= 0, input: q,ffc,isign,nx,nxhm,nxvh output: fx,we
c equation used is:
c fx(kx) = -sqrt(-1)*kx*g(kx)*s(kx)*q(kx),
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c g(k) = (affp/k**2)*s(k), and s(k) = exp(-(k*ax)**2/2)
c q(j) = complex charge density for fourier mode (j-1)
c fx(j) = x component of complex force/charge, or fourier mode (j-1)
c if isign = 0, form factor array is prepared
c aimag(ffc(j)) = finite-size particle shape factor s
c real(ffc(j)) = potential green's function g
c for fourier mode (j-1)
c ax = half-width of particle in x direction
c affp = normalization constant = nx/np, where np = number of particles
c electric field energy is also calculated, using
c we = nx*sum((affp/k**2)*|q(k)*s(k)|**2)
c nx = system length in x direction
c nxhm = number of fourier modes kept
c nxvh = dimension of field array, must be >= nxhm
      implicit none
      integer isign, nx, nxhm, nxvh
      real ax, affp, we
      complex q, fx, ffc
      dimension q(nxvh), fx(nxvh), ffc(nxhm)
c local data
      integer j
      real dnx, dkx, at1, at2
      complex zt1, zero
      double precision wp
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.0,0.0)
      if (isign.ne.0) go to 20
c prepare form factor array
      do 10 j = 2, nxhm
      dkx = dnx*real(j - 1)
      at1 = exp(-.5*(dkx*ax)**2)
      ffc(j) = cmplx(affp*at1/(dkx*dkx),at1)
   10 continue
      ffc(1) = cmplx(affp,1.0)
      return
      if (isign.ne.0) go to 20
c calculate force/charge and sum field energy
   20 wp = 0.0d0
c mode numbers 0 < kx < nxhm
      do 30 j = 2, nxhm
      at1 = real(ffc(j))*aimag(ffc(j))
      at2 = dnx*real(j - 1)*at1
      zt1 = cmplx(aimag(q(j)),-real(q(j)))
      fx(j) = at2*zt1
      wp = wp + at1*q(j)*conjg(q(j))
   30 continue
c mode number kx = 0
      fx(1) = zero
      we = real(nx)*wp
      return
      return
      end
c-----------------------------------------------------------------------
      function ranorm()
c this program calculates a random number y from a gaussian distribution
c with zero mean and unit variance, according to the method of
c mueller and box:
c    y(k) = (-2*ln(x(k)))**1/2*sin(2*pi*x(k+1))
c    y(k+1) = (-2*ln(x(k)))**1/2*cos(2*pi*x(k+1)),
c where x is a random number uniformly distributed on (0,1).
c written for the ibm by viktor k. decyk, ucla
      implicit none
      integer iflg,isc,i1,r1,r2,r4,r5
      double precision ranorm,h1l,h1u,h2l,r0,r3,asc,bsc,temp
      save iflg,r1,r2,r4,r5,h1l,h1u,h2l,r0
      data r1,r2,r4,r5 /885098780,1824280461,1396483093,55318673/
      data h1l,h1u,h2l /65531.0d0,32767.0d0,65525.0d0/
      data iflg,r0 /0,0.0d0/
      if (iflg.eq.0) go to 10
      ranorm = r0
      r0 = 0.0d0
      iflg = 0
      return
   10 isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1 - (r1/isc)*isc
      r3 = h1l*dble(r1) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2/isc
      isc = r2 - i1*isc
      r0 = h1l*dble(r2) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1 = r3 - dble(isc)*bsc
      temp = dsqrt(-2.0d0*dlog((dble(r1) + dble(r2)*asc)*asc))
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r4 - (r4/isc)*isc
      r3 = h2l*dble(r4) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r5/isc
      isc = r5 - i1*isc
      r0 = h2l*dble(r5) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r5 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r4 = r3 - dble(isc)*bsc
      r0 = 6.28318530717959d0*((dble(r4) + dble(r5)*asc)*asc)
      ranorm = temp*dsin(r0)
      r0 = temp*dcos(r0)
      iflg = 1
      return
      end
c-----------------------------------------------------------------------
      function randum()
c this is a version of the random number generator dprandom due to
c c. bingham and the yale computer center, producing numbers
c in the interval (0,1).  written for the sun by viktor k. decyk, ucla
      implicit none
      integer isc,i1,r1,r2
      double precision randum,h1l,h1u,r0,r3,asc,bsc
      save r1,r2,h1l,h1u
      data r1,r2 /1271199957,1013501921/
      data h1l,h1u /65533.0d0,32767.0d0/
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1 - (r1/isc)*isc
      r3 = h1l*dble(r1) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2/isc
      isc = r2 - i1*isc
      r0 = h1l*dble(r2) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1 = r3 - dble(isc)*bsc
      randum = (dble(r1) + dble(r2)*asc)*asc
      return
      end
c The following extra procedure is not currently used:
c-----------------------------------------------------------------------
      subroutine MEFIELD1GL(fx,gx,xp,nx,nxhm,nxvh)
c for 1d code, this subroutine calculates real space 1d field gx at
c location xp, from complex fourier co-efficients fx
c OpenMP version, parallelization over fourier modes
c (16 flops + tan + div)*nxhm, nxhm loads, 1 store
c input: all, output: gx
c equations used are:
c gx = sum(fx(n)*exp(sqrt(-1)*2*n*pi*xp/nx)
c fx(n) = x component of force/charge at fourier mode n-1
c if nxhm=nx/2, nx/2 mode is included, in a format compatible with the
c real to complex FFT used in traditional PIC
c gx = x component of force/charge at location xp
c xp = position to evaluate field
c nx = system length in x direction
c nxhm = number of fourier modes kept in x
c nxvh = second dimension of vector field arrays, must be >= nxhm
      implicit none
      integer nx, nxhm, nxvh
      real gx, xp
      complex fx
      dimension fx(nxvh)
c local data
      integer j
      real dnx, at1, dkx, tn, at2
      double precision ex
      complex zt1
      dnx = 6.28318530717959/real(nx)
c find field
      ex = 0.0d0
c mode numbers 0 < kx < nxhm
      at1 = dnx*xp
!$OMP PARALLEL DO PRIVATE(j,dkx,tn,at2,zt1) REDUCTION(+:ex)
      do 10 j = 2, nxhm
      dkx = at1*real(j-1)
c zt1 = cmplx(cos(dkx),sin(dkx))
      tn = tan(0.5*dkx)
      at2 = (tn + tn)/(1 + tn*tn)
      zt1 = cmplx(1.0-tn*at2,at2)
      ex = ex + real(fx(j)*zt1)
   10 continue
!$OMP END PARALLEL DO
      ex = 2.0d0*ex + real(fx(1))
c special case to match conventional PIC code
      if (nxhm.eq.(nx/2)) then
         dkx = at1*real(nxhm)
         at1 = cos(dkx)
         ex = ex + aimag(fx(1))*at1
      endif
      gx = real(ex)
      return
      end
