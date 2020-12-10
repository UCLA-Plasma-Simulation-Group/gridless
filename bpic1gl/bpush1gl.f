c Fortran Library for Skeleton 1-2/2D Electromagnetic Gridless PIC Code
c written by Viktor K. Decyk, UCLA
c-----------------------------------------------------------------------
      subroutine DISTR1H(part,vtx,vty,vtz,vdx,vdy,vdz,npx,idimp,nop,nx, 
     1ipbc)
c for 1-2/2d code, this subroutine calculates initial particle
c co-ordinate and velocity, with uniform density and maxwellian
c velocity with drift
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c part(3,n) = velocity vy of particle n
c part(4,n) = velocity vz of particle n
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c npx = number of particles distributed in x direction
c idimp = size of phase space = 4
c nop = number of particles
c nx = system length in x direction
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer npx, idimp, nop, nx, ipbc
      real part, vtx, vty, vtz, vdx, vdy,v dz
      dimension part(idimp,nop)
c local data
      integer j
      real edgelx, at1, sum1, sum2, sum3
      double precision dsum1, dsum2, dsum3
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
      part(3,j) = vty*ranorm()
      part(4,j) = vtz*ranorm()
   20 continue
c add correct drift
      dsum1 = 0.0d0
      dsum2 = 0.0d0
      dsum3 = 0.0d0
      do 30 j = 1, npx
      dsum1 = dsum1 + part(2,j)
      dsum2 = dsum2 + part(3,j)
      dsum3 = dsum3 + part(4,j)
   30 continue
      sum1 = dsum1
      sum2 = dsum2
      sum3 = dsum3
      at1 = 1.0/real(npx)
      sum1 = at1*sum1 - vdx
      sum2 = at1*sum2 - vdy
      sum3 = at1*sum3 - vdz
      do 40 j = 1, npx
      part(2,j) = part(2,j) - sum1
      part(3,j) = part(3,j) - sum2
      part(4,j) = part(4,j) - sum3
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine BPUSH13GL(part,fxyz,byz,sctx,omx,qbm,dt,dtc,ek,idimp,  
     1nop,nx,nxhm,nxvh)
c for 1-2/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time using gridless spectral
c version, with magnetic field and periodic boundaries.
c Using the Boris Mover.
c baseline scalar version
c 48*nxhm + 68 flops/particle, 5*nxhm + 4 loads, 4 stores, 1 divide
c input: all, output: part, sctx, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fx(x(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fy(x(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fz(x(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omy = (q/m)*by(x(t)), and omz = (q/m)*bz(x(t)).
c position equation used is:
c x(t+dt) = x(t) + vx(t+dt/2)*dt
c fx(x(t)) is calculated from the expression
c fx(x) = sum(fxyz(1,n)*exp(sqrt(-1)*2*n*pi*x/nx))
c similarly for fy(x), fz(x), by(x), bz(x)
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c part(3,n) = velocity vy of particle n
c part(4,n) = velocity vz of particle n
c fxyz(i,n) = i component of force/charge at fourier grid point n
c byz(1,n) = y component of magnetic field at fourier grid point n
c byz(2,n) = z component of magnetic field at fourier grid point n
c sctx = scratch array for sines and cosines
c omx = magnetic field electron cyclotron frequency in x
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt)**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx = system length in x direction
c nxhm = number of fourier modes kept
c nxvh = dimension of field array, must be >= nxhm
      implicit none
      integer idimp, nop, nx, nxhm, nxvh
      real part, omx, qbm, dt, dtc, ek
      complex fxyz, byz, sctx
      dimension part(idimp,nop), fxyz(3,nxvh), byz(2,nxvh)
      dimension sctx(nxvh)
c local data
      integer i, j
      real zero, anx, dnx, dkx, at1
      real qtmh, dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt
      real omt, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real anorm
      double precision sum1, ex, ey, ez, by, bz
      complex zt1
      zero = 0.0
      anx = real(nx)
      dnx = 6.28318530717959/anx
      qtmh = 0.5*qbm*dt
      sum1 = 0.0d0
      do 30 i = 1, nop
c find electric and magnetic field
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
      by = 0.0d0
      bz = 0.0d0
c mode numbers 0 < kx < nxhm
      do 20 j = 2, nxhm
      zt1 = sctx(j-1)
      ex = ex + real(fxyz(1,j)*zt1)
      ey = ey + real(fxyz(2,j)*zt1)
      ez = ez + real(fxyz(3,j)*zt1)
      by = by + real(byz(1,j)*zt1)
      bz = bz + real(byz(2,j)*zt1)
   20 continue
      ex = 2.0d0*ex + real(fxyz(1,1))
      ey = 2.0d0*ey + real(fxyz(2,1))
      ez = 2.0d0*ez + real(fxyz(3,1))
      by = 2.0d0*by + real(byz(1,1))
      bz = 2.0d0*bz + real(byz(2,1))
c special case to match conventional PIC code
      if (nxhm.eq.(nx/2)) then
         at1 = real(sctx(nxhm))
         dx = ex + aimag(fxyz(1,1))*at1
         dy = ey + aimag(fxyz(2,1))*at1
         dz = ez + aimag(fxyz(3,1))*at1
         oy = by + aimag(byz(1,1))*at1
         oz = bz + aimag(byz(2,1))*at1
      else
         dx = ex
         dy = ey
         dz = ez
         oy = by
         oz = bz
      endif
      ox = omx
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(2,i) + dx
      acy = part(3,i) + dy
      acz = part(4,i) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(2,i) = dx
      part(3,i) = dy
      part(4,i) = dz
c new position
      dx = part(1,i) + dx*dtc
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
c set new position
      part(1,i) = dx
   30 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine RBPUSH13GL(part,fxyz,byz,sctx,omx,qbm,dt,dtc,ci,ek,    
     1idimp,nop,nx,nxhm,nxvh)
c for 1-2/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time using gridless spectral
c version, for relativistic particles with magnetic field and periodic
c boundaries
c Using the Boris Mover.
c scalar version using guard cells
c 48*nxhm + 78 flops/particle, 5*nxhm + 4 loads, 4 stores
c 4 divides, 2 sqrts
c input: all, output: part, sctx, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),)*dt) +
c    .5*(q/m)*fx(x(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t)))*dt) +
c    .5*(q/m)*fy(x(t),)*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) +
c    .5*(q/m)*fz(x(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omy = (q/m)*by(x(t))*gami and omz = (q/m)*bz(x(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equation used is:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c fx(x(t)) is calculated from the expression
c fx(x) = sum(fxyz(1,n)*exp(sqrt(-1)*2*n*pi*x/nx))
c similarly for fy(x), fz(x), by(x), bz(x)
c part(1,n) = position x of particle n
c part(2,n) = momentum px of particle n
c part(3,n) = momentum py of particle n
c part(4,n) = momentum pz of particle n
c fxyz(i,n) = i component of force/charge at fourier grid point n
c byz(1,n) = y component of magnetic field at fourier grid point n
c byz(2,n) = z component of magnetic field at fourier grid point n
c sctx = scratch array for sines and cosines
c omx = magnetic field electron cyclotron frequency in x
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t))*dt)**2)/(1. + gami)
c idimp = size of phase space = 4
c nop = number of particles
c nx = system length in x direction
c nxhm = number of fourier modes kept
c nxvh = dimension of field array, must be >= nxhm
      implicit none
      integer idimp, nop, nx, nxhm, nxvh
      real part, omx, qbm, dt, dtc, ci, ek
      complex fxyz, byz, sctx
      dimension part(idimp,nop), fxyz(3,nxvh), byz(2,nxvh)
      dimension sctx(nxvh)
c local data
      integer i, j
      real zero, anx, dnx, dkx, at1
      real qtmh, dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt
      real omt, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real anorm, ci2, p2, gami, qtmg, dtg
      double precision sum1, ex, ey, ez, by, bz
      complex zt1
      zero = 0.0
      anx = real(nx)
      dnx = 6.28318530717959/anx
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      sum1 = 0.0d0
      do 30 i = 1, nop
c find electric and magnetic field
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
      by = 0.0d0
      bz = 0.0d0
c mode numbers 0 < kx < nxhm
      do 20 j = 2, nxhm
      zt1 = sctx(j-1)
      ex = ex + real(fxyz(1,j)*zt1)
      ey = ey + real(fxyz(2,j)*zt1)
      ez = ez + real(fxyz(3,j)*zt1)
      by = by + real(byz(1,j)*zt1)
      bz = bz + real(byz(2,j)*zt1)
   20 continue
      ex = 2.0d0*ex + real(fxyz(1,1))
      ey = 2.0d0*ey + real(fxyz(2,1))
      ez = 2.0d0*ez + real(fxyz(3,1))
      by = 2.0d0*by + real(byz(1,1))
      bz = 2.0d0*bz + real(byz(2,1))
c special case to match conventional PIC code
      if (nxhm.eq.(nx/2)) then
         at1 = real(sctx(nxhm))
         dx = ex + aimag(fxyz(1,1))*at1
         dy = ey + aimag(fxyz(2,1))*at1
         dz = ez + aimag(fxyz(3,1))*at1
         oy = by + aimag(byz(1,1))*at1
         oz = bz + aimag(byz(2,1))*at1
      else
         dx = ex
         dy = ey
         dz = ez
         oy = by
         oz = bz
      endif
      ox = omx
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(2,i) + dx
      acy = part(3,i) + dy
      acz = part(4,i) + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c renormalize magnetic field
      qtmg = qtmh*gami
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(2,i) = dx
      part(3,i) = dy
      part(4,i) = dz
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,i) + dx*dtg
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
c set new position
      part(1,i) = dx
   30 continue
c normalize kinetic energy
      ek = ek + sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine ABPUSH13GL(part,fxyz,byz,sctx,omx,qbm,dt,dtc,ek,idimp, 
     1nop,nx,nxhm,nxvh)
c for 1-2/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time using gridless spectral
c version, with magnetic field and periodic boundaries.
c Using the Analytic Boris Mover,
c assumes constant E, B fields during a time step
c baseline scalar version
c 48*nxhm + 110 flops/particle, 5*nxhm + 4 loads, 4 stores
c 2 divides, 1 sqrt, 1 tangent
c input: all, output: part, sctx, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + (q/m)*gx(x(t))) +
c    rot(2)*(vy(t-dt/2) + (q/m)*gy(x(t))) +
c    rot(3)*(vz(t-dt/2) + (q/m)*gz(x(t)))) + (q/m)*gx(x(t)))
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + (q/m)*gx(x(t))) +
c    rot(5)*(vy(t-dt/2) + (q/m)*gy(x(t))) +
c    rot(6)*(vz(t-dt/2) + (q/m)*gz(x(t)))) + (q/m)*gy(x(t)))
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + (q/m)*gx(x(t))) +
c    rot(8)*(vy(t-dt/2) + (q/m)*gy(x(t))t) +
c    rot(9)*(vz(t-dt/2) + (q/m)*gz(x(t)))) + (q/m)*gz(x(t)))
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - tan(om*dt/2)**2 + 2*((omx/om)*tan(om*dt/2))**2)/norm
c    rot(2) = 2*((omz/om)*tan(om*dt/2) 
c              + (omx*omy/om**2)*tan(om*dt/2)**2)/norm
c    rot(3) = 2*(-(omy/om)*tan(om*dt/2)
c              + (omx*omz/om**2)*tan(om*dt/2)**2)/norm
c    rot(4) = 2*(-(omz/om)*tan(om*dt/2) 
c              + (omx*omy/om**2)*tan(om*dt/2)**2)/norm
c    rot(5) = (1 - tan(om*dt/2)**2 + 2*((omy/om)*tan(om*dt/2))**2)/norm
c    rot(6) = 2*((omx/om)*tan(om*dt/2) 
c              + (omy*omz/om**2)*tan(om*dt/2)**2)/norm
c    rot(7) = 2*((omy/om)*tan(om*dt/2)
c              + (omx*omz/om**2)*tan(om*dt/2)**2)/norm
c    rot(8) = 2*(-(omx/om)*tan(om*dt/2)
c              + (omy*omz/om**2)*tan(om*dt/2)**2)/norm
c    rot(9) = (1 - tan(om*dt/2)**2 + 2*((omz/om)*tan(om*dt/2))**2)/norm
c norm = 1 + tan(om*dt/2))**2 and om = sqrt(omx**2 + omy**2 + omz**2)
c gx(x) = 0.5*flx(x)*dt + (fx(x)-flx(x))*tan(0.5*om*dt)/om
c gy(x) = 0.5*fly(x)*dt + (fy(x)-fly(x))*tan(0.5*om*dt)/om
c gz(x) = 0.5*flz(x)*dt + (fz(x)-flz(x))*tan(0.5*om*dt)/om
c where flx(x) = fpl(x)*omx/om**2, fly(x) = fpl(x)*omy/om**2,
c flz(x) = fpl(x)*omz/om**2, and fpl(x) = fx(x)*omx+fy(x)*omy+fz(x)*omz
c the rotation matrix is determined by:
c omy = (q/m)*by(x(t)), and omz = (q/m)*bz(x(t)).
c position equation used is:
c x(t+dt) = x(t) + vx(t+dt/2)*dt
c fx(x(t)) is calculated from the expression
c fx(x) = sum(fxyz(1,n)*exp(sqrt(-1)*2*n*pi*x/nx))
c similarly for fy(x), fz(x), by(x), bz(x)
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c part(3,n) = velocity vy of particle n
c part(4,n) = velocity vz of particle n
c fxyz(i,n) = i component of force/charge at fourier grid point n
c byz(1,n) = y component of magnetic field at fourier grid point n
c byz(2,n) = z component of magnetic field at fourier grid point n
c sctx = scratch array for sines and cosines
c omx = magnetic field electron cyclotron frequency in x
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt)**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx = system length in x direction
c nxhm = number of fourier modes kept
c nxvh = dimension of field array, must be >= nxhm
      implicit none
      integer idimp, nop, nx, nxhm, nxvh
      real part, omx, qbm, dt, dtc, ek
      complex fxyz, byz, sctx
      dimension part(idimp,nop), fxyz(3,nxvh), byz(2,nxvh)
      dimension sctx(nxvh)
c local data
      integer i, j
      real zero, anx, dth, dnx, dkx, at1, omt, omti, epl
      real qtmh, dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1, ex, ey, ez, by, bz
      complex zt1
      zero = 0.0
      anx = real(nx)
      dth = 0.5*dt
      dnx = 6.28318530717959/anx
      sum1 = 0.0d0
      do 30 i = 1, nop
c find electric and magnetic field
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
      by = 0.0d0
      bz = 0.0d0
c mode numbers 0 < kx < nxhm
      do 20 j = 2, nxhm
      zt1 = sctx(j-1)
      ex = ex + real(fxyz(1,j)*zt1)
      ey = ey + real(fxyz(2,j)*zt1)
      ez = ez + real(fxyz(3,j)*zt1)
      by = by + real(byz(1,j)*zt1)
      bz = bz + real(byz(2,j)*zt1)
   20 continue
      ex = 2.0d0*ex + real(fxyz(1,1))
      ey = 2.0d0*ey + real(fxyz(2,1))
      ez = 2.0d0*ez + real(fxyz(3,1))
      by = 2.0d0*by + real(byz(1,1))
      bz = 2.0d0*bz + real(byz(2,1))
c special case to match conventional PIC code
      if (nxhm.eq.(nx/2)) then
         at1 = real(sctx(nxhm))
         dx = ex + aimag(fxyz(1,1))*at1
         dy = ey + aimag(fxyz(2,1))*at1
         dz = ez + aimag(fxyz(3,1))*at1
         oy = by + aimag(byz(1,1))*at1
         oz = bz + aimag(byz(2,1))*at1
      else
         dx = ex
         dy = ey
         dz = ez
         oy = by
         oz = bz
      endif
      ox = omx
c normalize electric field
      dx = qbm*dx
      dy = qbm*dy
      dz = qbm*dz
c half acceleration
      acx = part(2,i)
      acy = part(3,i)
      acz = part(4,i)
      omxt = acx + dx*dth
      omyt = acy + dy*dth
      omzt = acz + dz*dth
c time-centered kinetic energy
      sum1 = sum1 + (omxt*omxt + omyt*omyt + omzt*omzt)
c normalize magnetic field
      ox = qbm*ox
      oy = qbm*oy
      oz = qbm*oz
      qtmh = dth
c correct the half-acceleration by decomposing E field into components
c parallel and perpendicular to the B field
      omt = sqrt(ox*ox + oy*oy + oz*oz)
      omti = 0.0
      epl = dx*ox + dy*oy + dz*oz
      if (omt.gt.0.0) then
         omti = 1.0/omt
         qtmh = omti*tan(dth*omt)
      endif
      epl = epl*omti*omti
      omxt = epl*ox
      omyt = epl*oy
      omzt = epl*oz
      dx = omxt*dth + (dx - omxt)*qtmh
      dy = omyt*dth + (dy - omyt)*qtmh
      dz = omzt*dth + (dz - omzt)*qtmh
      acx = acx + dx
      acy = acy + dy
      acz = acz + dz
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(2,i) = dx
      part(3,i) = dy
      part(4,i) = dz
c new position
      dx = part(1,i) + dx*dtc
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
c set new position
      part(1,i) = dx
   30 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine ARBPUSH13GL(part,fxyz,byz,sctx,omx,qbm,dt,dtc,ci,ek,   
     1idimp,nop,nx,nxhm,nxvh)
c for 1-2/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time using gridless spectral
c version, for relativistic particles with magnetic field and periodic
c boundaries
c Using the Analytic Boris Mover,
c assumes constant E, B fields, and gamma during a time step
c scalar version using guard cells
c 48*nxhm + 194 flops/particle, 5*nxhm + 4 loads, 4 stores
c 6 divides, 3 sqrts, 1 tangent
c input: all, output: part, sctx, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + (q/m)*gx(x(t))) +
c    rot(2)*(py(t-dt/2) + (q/m)*gy(x(t))) +
c    rot(3)*(pz(t-dt/2) + (q/m)*gz(x(t)))) + (q/m)*gx(x(t)))
c py(t+dt/2) = rot(4)*(px(t-dt/2) + (q/m)*gx(x(t))) +
c    rot(5)*(py(t-dt/2) + (q/m)*gy(x(t))) +
c    rot(6)*(pz(t-dt/2) + (q/m)*gz(x(t)))) + (q/m)*gy(x(t)))
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + (q/m)*gx(x(t))) +
c    rot(8)*(py(t-dt/2) + (q/m)*gy(x(t))t) +
c    rot(9)*(pz(t-dt/2) + (q/m)*gz(x(t)))) + (q/m)*gz(x(t)))
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - tan(om*dt/2)**2 + 2*((omx/om)*tan(om*dt/2))**2)/norm
c    rot(2) = 2*((omz/om)*tan(om*dt/2) 
c              + (omx*omy/om**2)*tan(om*dt/2)**2)/norm
c    rot(3) = 2*(-(omy/om)*tan(om*dt/2)
c              + (omx*omz/om**2)*tan(om*dt/2)**2)/norm
c    rot(4) = 2*(-(omz/om)*tan(om*dt/2) 
c              + (omx*omy/om**2)*tan(om*dt/2)**2)/norm
c    rot(5) = (1 - tan(om*dt/2)**2 + 2*((omy/om)*tan(om*dt/2))**2)/norm
c    rot(6) = 2*((omx/om)*tan(om*dt/2) 
c              + (omy*omz/om**2)*tan(om*dt/2)**2)/norm
c    rot(7) = 2*((omy/om)*tan(om*dt/2)
c              + (omx*omz/om**2)*tan(om*dt/2)**2)/norm
c    rot(8) = 2*(-(omx/om)*tan(om*dt/2)
c              + (omy*omz/om**2)*tan(om*dt/2)**2)/norm
c    rot(9) = (1 - tan(om*dt/2)**2 + 2*((omz/om)*tan(om*dt/2))**2)/norm
c norm = 1 + tan(om*dt/2))**2 and om = sqrt(omx**2 + omy**2 + omz**2)
c gx(x) = 0.5*flx(x)*dt + (fx(x)-flx(x))*tan(0.5*om*dt)/om
c gy(x) = 0.5*fly(x)*dt + (fy(x)-fly(x))*tan(0.5*om*dt)/om
c gz(x) = 0.5*flz(x)*dt + (fz(x)-flz(x))*tan(0.5*om*dt)/om
c where flx(x) = fpl(x)*omx/om**2, fly(x) = fpl(x)*omy/om**2,
c flz(x) = fpl(x)*omz/om**2, and fpl(x) = fx(x)*omx+fy(x)*omy+fz(x)*omz
c the rotation matrix is determined by:
c omy = (q/m)*by(x(t))*gami and omz = (q/m)*bz(x(t))*gami, where
c we approximate analytic gamma function with 4th order taylor series:
c gam = gam(0) + dg0*tau+d2g0*tau**2/2+d3g0*tau**3/6+d4g0*tau**3/24
c where gam(0) = sqrt(1.0 + p2(0)*ci*ci)
c p2(0) = px(t-dt/2)**2+py(t-dt/2)**2+pz(t-dt/2)**2
c dg0 = ug0/(ci*ci) = dgamma/dtau
c d2g0 = u2g0/(ci*ci) = d2gamma/dtau2
c d3g0 = u3g0/(ci*ci) = -omt*omt*(dgamma_perp/dtau)
c d4g0 = u4g0/(ci*ci) = -omt*omt*(d2gamma_perp/dtau2)
c then using the result t = integral of gamma(tau)*dtau, we can
c approximate tau(t) using a one pass Newton method and set gami = tau/t
c position equation used is:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c fx(x(t)) is calculated from the expression
c fx(x) = sum(fxyz(1,n)*exp(sqrt(-1)*2*n*pi*x/nx))
c similarly for fy(x), fz(x), by(x), bz(x)
c part(1,n) = position x of particle n
c part(2,n) = momentum px of particle n
c part(3,n) = momentum py of particle n
c part(4,n) = momentum pz of particle n
c fxyz(i,n) = i component of force/charge at fourier grid point n
c byz(1,n) = y component of magnetic field at fourier grid point n
c byz(2,n) = z component of magnetic field at fourier grid point n
c sctx = scratch array for sines and cosines
c omx = magnetic field electron cyclotron frequency in x
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = p2(0)/(1.+sqrt(1.+p2(0)*ci*ci))+ug0*th+u2g0*th**2/2+u3g0*th**3/6
c where th = tau/2
c idimp = size of phase space = 4
c nop = number of particles
c nx = system length in x direction
c nxhm = number of fourier modes kept
c nxvh = dimension of field array, must be >= nxhm
      implicit none
      integer idimp, nop, nx, nxhm, nxvh
      real part, omx, qbm, dt, dtc, ci, ek
      complex fxyz, byz, sctx
      dimension part(idimp,nop), fxyz(3,nxvh), byz(2,nxvh)
      dimension sctx(nxvh)
c local data
      integer i, j
      real zero, sixth, anx, dth, dnx, dkx, ci2, at1
      real dx, dy, dz, ox, oy, oz, px, py, pz, p2, gam0, wk, ug0, u2g0
      real acx, acy, acz, omt, omti, epl, omxt, omyt, omzt, ugt0, u2gt0
      real gami, tau, f, fp, qtmg, dtg
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1, ex, ey, ez, by, bz
      complex zt1
      zero = 0.0
      sixth = 1.0/6.0
      anx = real(nx)
      dth = 0.5*dt
      dnx = 6.28318530717959/anx
      ci2 = ci*ci
      sum1 = 0.0d0
      do 30 i = 1, nop
c find electric and magnetic field
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
      by = 0.0d0
      bz = 0.0d0
c mode numbers 0 < kx < nxhm
      do 20 j = 2, nxhm
      zt1 = sctx(j-1)
      ex = ex + real(fxyz(1,j)*zt1)
      ey = ey + real(fxyz(2,j)*zt1)
      ez = ez + real(fxyz(3,j)*zt1)
      by = by + real(byz(1,j)*zt1)
      bz = bz + real(byz(2,j)*zt1)
   20 continue
      ex = 2.0d0*ex + real(fxyz(1,1))
      ey = 2.0d0*ey + real(fxyz(2,1))
      ez = 2.0d0*ez + real(fxyz(3,1))
      by = 2.0d0*by + real(byz(1,1))
      bz = 2.0d0*bz + real(byz(2,1))
c special case to match conventional PIC code
      if (nxhm.eq.(nx/2)) then
         at1 = real(sctx(nxhm))
         dx = ex + aimag(fxyz(1,1))*at1
         dy = ey + aimag(fxyz(2,1))*at1
         dz = ez + aimag(fxyz(3,1))*at1
         oy = by + aimag(byz(1,1))*at1
         oz = bz + aimag(byz(2,1))*at1
      else
         dx = ex
         dy = ey
         dz = ez
         oy = by
         oz = bz
      endif
      ox = omx
c normalize electric field
      dx = qbm*dx
      dy = qbm*dy
      dz = qbm*dz
c normalize magnetic field
      ox = qbm*ox
      oy = qbm*oy
      oz = qbm*oz
c read momentum
      px = part(2,i)
      py = part(3,i)
      pz = part(4,i)
c find initial gamma
      p2 = px*px + py*py + pz*pz
      gam0 = sqrt(1.0 + p2*ci2)
c initial kinetic energy
      wk = p2/(1.0 + gam0)
c dgamma/dtau
      ug0 = dx*px + dy*py + dz*pz
c d2gamma/dtau2
      acx = py*oz - pz*oy
      acy = pz*ox - px*oz
      acz = px*oy - py*ox
      u2g0 = gam0*(dx*dx + dy*dy + dz*dz) + (dx*acx + dy*acy + dz*acz)
c correct the half-acceleration by decomposing E field into components
c parallel and perpendicular to the B field
      omt = sqrt(ox*ox + oy*oy + oz*oz)
      omti = 0.0
      if (omt.gt.0.0) omti = 1.0/omt
      epl = dx*ox + dy*oy + dz*oz
      epl = epl*omti*omti
c E parallel
      omxt = epl*ox
      omyt = epl*oy
      omzt = epl*oz
c E perp
      dx = dx - omxt
      dy = dy - omyt
      dz = dz - omzt
c dgamma_perp/dtau
      ugt0 = dx*px + dy*py + dz*pz
c d2gamma_perp/dtau2
      u2gt0 = gam0*(dx*dx + dy*dy + dz*dz) + (dx*acx + dy*acy + dz*acz)
c find gami with one pass Newton method and fourth order taylor series
      gami = 1.0/gam0
      tau = dt*gami
      at1 = omt*omt
      f = (0.5*ug0 + sixth*(u2g0 - 0.25*at1*(ugt0 + 0.2*(u2gt0*tau))*tau
     1)*tau)*tau
      fp = (ug0 + (0.5*u2g0 - at1*sixth*(ugt0 + 0.25*u2gt0*tau)*tau)*tau
     1)*tau
      gami = gami*(1.0 - f*ci2/(gam0 + fp*ci2))
c time-centered kinetic energy
      at1 = -at1*ugt0
      wk = wk + 0.5*(ug0 + 0.25*(u2g0 + sixth*at1*tau)*tau)*tau
      sum1 = sum1 + wk
c set proper time
      qtmg = dth*gami
      if (omt.gt.0.0) qtmg = omti*tan(qtmg*omt)
      at1 = qtmg/gami
c modified half acceleration
      dx = omxt*dth + dx*at1
      dy = omyt*dth + dy*at1
      dz = omzt*dth + dz*at1
      acx = px + dx
      acy = py + dy
      acz = pz + dz
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      px = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      py = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      pz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(2,i) = px
      part(3,i) = py
      part(4,i) = pz
c update inverse gamma
      p2 = px*px + py*py + pz*pz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,i) + px*dtg
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
c set new position
      part(1,i) = dx
   30 continue
c normalize kinetic energy
      ek = ek + sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine EARBPUSH13GL(part,fxyz,byz,sctx,omx,qbm,dt,dtc,ci,ek,  
     1idimp,nop,nx,nxhm,nxvh)
c for 1-2/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time using gridless spectral
c version, for relativistic particles with magnetic field and periodic
c boundaries
c Using the Exact Analytic mover, a variant of the algorithm developed
c by Fei Li et al., where E and B fields are constant and gamma varies
c during a time step.
c scalar version using guard cells
c (48 FLOPs+sin+cos)*nxhm + 186 FLOPs, 11 divides, 5 sqrts per particle
c plus 51 FLOPs, 2 divides, 1 tan, 1 tanh per Newton iteration/particle
c input: all, output: part, sctx, ek
c The particle momenta are advanced in natural coordinates, defined as:
c up is in the direction ep = E parallel to the B field,
c ut is in the direction of et = E perpendicular to the B field,
c ud is in the ExB direction.  Momenta are transformed from cartesian
c coordinates to natural coordinates before the advance, then
c transformed back to cartesian coordinates.
c momenta equations used are:
c up = up + ep*dt
c ut = ut + etp*dt + (om/w)*(F*(cos(w*tau)-1.0) + G*sin(w*tau))
c ud = ved*gam(tau) + F*sin(w*tau) - H*cos(w*tau), where
c F = (w/omega)*ut(0) - om2t*(omega*w)*ci2*ep*up(0)
c G = gam0*(et - al2*om2t) - ud(0)*omega
c H = om2t*omega*gam0 - ud(0)
c om = sqrt(qbm*omx)**2 + (qbm*by)**2 + (qbm*bz)**2)
c etp = om2t*al*al, ved = om*om2t, om2t = et/(om*om + al*al
c the gyration frequency w and al are given by:
c w = sqrt(0.5*(sqrt(om*om-e*e)**2 + 4.0*(ci*ep*om)**2))+om*om-e*e))
c al = sqrt(0.5*(sqrt(om*om-e*e)**2 + 4.0*(ci*ep*om)**2))-om*om+e*e))
c where e = sqrt((qbm*ex)**2+(qbm*ey)**2+(qbm*ez)**2)
c gam(tau) is the relativistic factor given by:
c gam = A*sin(w*tau) - B*cos(w*tau) + C*sinh(al*tau) + D*cos(al*tau),
c where the constants A, B, C, D are defined in terms of initial
c derivatives of the gamma function, as follows:
c  w*A = (w*w*dg0 - om*om*dgp0)/(w*w + al*al)
c  w*B = (w*w*d2g0 - om*om*d2gp0)/(w*w + al*al)
c  al*C = (al*al*dg0 + om*om*dgp0)/(w*w + al*al)
c  al*C = (al*al*d2g0 + om*om*d2gp0)/(w*w + al*al)
c and where the initial derivatives of the gamma function are:
c  dgp0 = ci*ci*ep*up(0)
c  d2gp0 = ci*ci*ep*ep*gam(0)
c  dg0 = dgp0 + ci*ci*et*ut(0)
c  d2g0 = d2gp0 + ci*ci*et*(et*gam(0) - ud(0)*om)
c the proper time tau can be determined from the relativistic factor
c using the equation t = integral gam(tau)*dtau.  This involves finding
c zeros of the function f via an iteration scheme, where
c f = gam(0)*tau - dt - (A*(cos(w*tau)-1.0) + B*(sin(w*tau)-w*tau))/w
c                + (C*(cosh(al*tau)-1.0) + C*(sinh(al*tau)-w*tau))/al
c once we know tau(dt), we can evaluate the momenta up, ut, and ud
c and transform back to cartesian coordinates
c position equation used is:
c x(t+dt) = x(t) + px(t+dt/2)*dtg, where
c dtg = dtc/sqrt(1.0+(px(t+dt/2)**2+py(t+dt/2)**2+pz(t+dt/2)**2))*ci*ci)
c fx(x(t)) is calculated from the expression
c fx(x) = sum(fxyz(1,n)*exp(sqrt(-1)*2*n*pi*x/nx))
c similarly for fy(x), fz(x), by(x), bz(x)
c part(1,n) = position x of particle n
c part(2,n) = momentum px of particle n
c part(3,n) = momentum py of particle n
c part(4,n) = momentum pz of particle n
c fxyz(i,n) = i component of force/charge at fourier grid point n
c byz(1,n) = y component of magnetic field at fourier grid point n
c byz(2,n) = z component of magnetic field at fourier grid point n
c sctx = scratch array for sines and cosines
c omx = magnetic field electron cyclotron frequency in x
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using the average
c value of energies at the beginning and end of the time step:
c ek = 0.5*((px(t-dt/2)*2+py(t-dt/2)**2+pz(t-dt/2)**2)/(1.0+gam(0))
c    +      (px(t+dt/2)*2+py(t+dt/2)**2+pz(t+dt/2)**2)/(1.0+gam(tau)))
c idimp = size of phase space = 4
c nop = number of particles
c nx = system length in x direction
c nxhm = number of fourier modes kept
c nxvh = dimension of field array, must be >= nxhm
      implicit none
      integer idimp, nop, nx, nxhm, nxvh
      real part, omx, qbm, dt, dtc, ci, ek
      complex fxyz, byz, sctx
      dimension part(idimp,nop), fxyz(3,nxvh), byz(2,nxvh)
      dimension sctx(nxvh)
c local data
      integer imax
      real half, sixth, s1, s2, s3, s4, c1, c2, c3, c4, weps, teps, deps
      real erps
c imax = maximum iteration count
      parameter (imax=25)
      parameter (half=0.5,sixth=1.0/6.0)
      parameter (s1=1.0/6.0,s2=1.0/20.0,s3=1.0/42.0,s4=1.0/72.0)
      parameter (c1=0.5,c2=1.0/12.0,c3=1.0/30.0,c4=1.0/56.0)
c weps = taylor series expansion criteria
      parameter (weps=1.0e-1)
c deps = taylor series expansion criteria for differences
      parameter (deps=1.0e-3)
      integer i, j, ii
      real zero, anx, dnx, dkx, ci2, at1, small, prec, vresult, erm
      real dx, dy, dz, ox, oy, oz, tx, ty, tz, om2, e2, omega, di, ep
      real et, px, py, pz, up, ut, ud, w2, wl2, al2, w, al, wi, ali, wi2
      real ali2, om2t, ws, p2, gam0, dgp0, dgt0, d2gp0, d2gt0, dg0, d2g0
      real wt2, ea, eb, ec, ed, fa, fb, tau, cs, sn, csh, snh, fpp, fp
      real ft, f, t2, t3, wt, wd, tn, tn2, tn2i, csd, snd, cshd, snhd
      real csc, snc, cshc, snhc, fc, fpi, gam, gami, wk, dtg
      double precision sum1, ex, ey, ez, by, bz, dt1, err
      complex zt1
      data small /1.0e-12/
      prec = 1.0 + small
c teps = tau precision; erps = orthogonality precision
c detect autodouble precision
      if (vresult(prec).gt.1.0) then
         teps = 1.0e-14
         erps = 1.0e-10
c default single precision
      else
         teps = 4.0e-8
         erps = 4.0e-3
      endif
      zero = 0.0
      anx = real(nx)
      dnx = 6.28318530717959/anx
      ci2 = ci*ci
      erm = 0.0
      sum1 = 0.0d0
      do 50 i = 1, nop
c find electric and magnetic field
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
      by = 0.0d0
      bz = 0.0d0
c mode numbers 0 < kx < nxhm
      do 20 j = 2, nxhm
      zt1 = sctx(j-1)
      ex = ex + real(fxyz(1,j)*zt1)
      ey = ey + real(fxyz(2,j)*zt1)
      ez = ez + real(fxyz(3,j)*zt1)
      by = by + real(byz(1,j)*zt1)
      bz = bz + real(byz(2,j)*zt1)
   20 continue
      ex = 2.0d0*ex + real(fxyz(1,1))
      ey = 2.0d0*ey + real(fxyz(2,1))
      ez = 2.0d0*ez + real(fxyz(3,1))
      by = 2.0d0*by + real(byz(1,1))
      bz = 2.0d0*bz + real(byz(2,1))
c special case to match conventional PIC code
      if (nxhm.eq.(nx/2)) then
         at1 = real(sctx(nxhm))
         dx = ex + aimag(fxyz(1,1))*at1
         dy = ey + aimag(fxyz(2,1))*at1
         dz = ez + aimag(fxyz(3,1))*at1
         oy = by + aimag(byz(1,1))*at1
         oz = bz + aimag(byz(2,1))*at1
      else
         dx = ex
         dy = ey
         dz = ez
         oy = by
         oz = bz
      endif
c normalize magnetic field
      ox = qbm*omx
      oy = qbm*oy
      oz = qbm*oz
c normalize electric field
      tx = qbm*dx
      ty = qbm*dy
      tz = qbm*dz
c create direction cosines to translate from/to cartesian coordinates
c first find the direction along the magnetic field B
      om2 = ox*ox + oy*oy + oz*oz
      omega = sqrt(om2)
      if (om2.gt.0.0) then
         di = 1.0/omega
         ox = ox*di
         oy = oy*di
         oz = oz*di
c if omega = 0, then use the y direction
      else
         ox = 0.0
         oy = 1.0
         oz = 0.0
      endif
c then find the direction along the electric field Et perpendicular to B
      dt1 = dble(tx*ox) + dble(ty*oy) + dble(tz*oz)
      tx = tx - dt1*ox
      ty = ty - dt1*oy
      tz = tz - dt1*oz
      ep = dt1
      e2 = tx*tx + ty*ty + tz*tz
      et = sqrt(e2)
      if (et > 0.0) then
         di = 1.0/et
         tx = tx*di
         ty = ty*di
         tz = tz*di
c then find the direction along Et x B
         dx = ty*oz - tz*oy
         dy = tz*ox - tx*oz
         dz = tx*oy - ty*ox
c check for roundoff error
         err = dble(tx*ox) + dble(ty*oy) + dble(tz*oz)
         if (err > erps) then
            write (*,*) 'Error: Et not normal to omega=',et,err
            et = 0.0
         else
            err = dble(tx*dx) + dble(ty*dy) + dble(tz*dz)
            if (err > erps) then
               write (*,*) 'Error: Et d=',et,err
            endif
         endif
      endif
c special case Et = 0, or round off error detected
      if (et==0.0) then      
c first find direction with smallest component of B
         ii = 1
         dx = abs(ox)
         dy = abs(oy)
         dz = abs(oz)
         di = dx
         if (dy <= dx) then
            ii = 2
            di = dy
            if ((dx.gt.dy).and.(dy.eq.dz)) ii = 3
         endif
         if (dz.lt.di) ii = 3 
c then the cross product of that direction with B
         if (ii.eq.1) then
            dz = 1.0/sqrt(oy*oy + oz*oz)
            dx = 0.0
            dy = -oz*dz
            dz = oy*dz
         else if (ii.eq.2) then
            dz = 1.0/sqrt(ox*ox + oz*oz)
            dx = oz*dz
            dy = 0.0
            dz = -ox*dz
         else if (ii.eq.3) then
            dz = 1.0/sqrt(ox*ox + oy*oy)
            dx = -oy*dz
            dy = ox*dz
            dz = 0.0
         endif
c then find the direction along minus d x B
         tx = dz*oy - dy*oz
         ty = dx*oz - dz*ox
         tz = dy*ox - dx*oy
      endif
c calculate frequencies
      e2 = (e2 + ep*ep)*ci2
      w2 = om2 - e2
      wl2 = sqrt(w2*w2 + 4.0*ci2*(ep*omega)**2)
      al2 = 0.5*(wl2 - w2)
      w2 = 0.5*(wl2 + w2)
      w = sqrt(w2)
      al = sqrt(al2)
      wi = 0.0
      if (w > 0.0) wi = 1.0/w
      ali = 0.0
      if (al > 0.0) ali = 1.0/al
      wi2 = wi*wi
      ali2 = ali*ali
c calculate weights
      om2t = om2 + al2
      if (om2t > 0.0) om2t = et/om2t
      ws = 0.0
      if (omega /= 0.0) ws = w/omega
c translate momenta from cartesian coordinates
      px = part(2,i)
      py = part(3,i)
      pz = part(4,i)
      up = px*ox + py*oy + pz*oz
      ut = px*tx + py*ty + pz*tz
      ud = px*dx + py*dy + pz*dz
c find initial gamma
      p2 = up*up + ut*ut + ud*ud
      gam0 = sqrt(1.0 + p2*ci2)
c calculate initial kinetic energy
      wk = p2/(1.0 + gam0)
c partial derivatives of gamma
      dgp0 = ci2*ep*up
      dgt0 = ci2*et*ut
      d2gp0 = ci2*ep*ep*gam0
      d2gt0 = ci2*et*(et*gam0 - ud*omega)
      dg0 = dgp0 + dgt0
      d2g0 = d2gp0 + d2gt0
c calculate trigonometeric and hyperbolic coefficients
      wt2 = 0.0
      if (wl2 > 0.0) wt2 = 1.0/wl2
      ea = wt2*(w2*dg0 - om2*dgp0)
      eb = wt2*(w2*d2g0 - om2*d2gp0)
      ec = wt2*(al2*dg0 + om2*dgp0)
      ed = wt2*(al2*d2g0 + om2*d2gp0)
      if (wl2==0.0) then
         ea = dgt0
         eb = d2gt0
      endif
      fa = ws*ut - om2t*omega*wi*ci2*ep*up
      fb = gam0*(et - al2*om2t) - ud*omega
      fc = om2t*omega*gam0 - ud
c zeroth order guess for tau
      tau = dt/gam0
c
c iteration loop for finding tau(t)
      cs = 1.0
      sn = 0.0
      csh = 1.0
      snh = 0.0
      fpp = 0.0
      fp = 0.0
      ft = -tau
      do 30 ii = 1, imax
      t2 = tau*tau; t3 = t2*tau
c calculate trigonometric functions
      wt = w*tau
      wd = w*ft
      if (wt > weps) then
         if (abs(wd) > deps) then
            tn = tan(0.5*wt)
            tn2 = tn*tn
            tn2i = 1.0/(1.0 + tn2)
            cs = (1.0 - tn2)*tn2i
            sn = (tn + tn)*tn2i
c second order taylor series approximation for wd
         else
            wt2 = wd*wd
            wd = -wd*(1.0 - s1*wt2)
            wt2 = 1.0 - c1*wt2
c third order taylor series approximation for wd
c           wd = -wd*(1.0 - s1*wt2*(1.0 - s2*wt2))
c           wt2 = 1.0 - c1*wt2*(1.0 - c2*wt2)
            f = cs*wt2 - sn*wd
            sn = sn*wt2 + cs*wd
            cs = f
         endif
c calculate special functions
         csc = cs - 1.0
         snc = sn*wi
         csd = csc*wi2
         snd = (snc - tau)*wi2
c fourth order taylor series approximation for wt
      else
         wt2 = wt*wt
         csd = -c1*(1.0 - c2*wt2*(1.0 - c3*wt2))*t2
         snd = -s1*(1.0 - s2*wt2*(1.0 - s3*wt2))*t3
c fifth order taylor series approximation for wt
c        csd = -c1*(1.0 - c2*wt2*(1.0 - c3*wt2*(1.0 - c4*wt2)))*t2
c        snd = -s1*(1.0 - s2*wt2*(1.0 - s3*wt2*(1.0 - s4*wt2)))*t3
         csc = w2*csd
         snc = tau + w2*snd
         cs = csc + 1.0
         sn = w*snc
      endif
c calculate hyperbolic functions
      wt = al*tau
      wd = al*ft
      if (wt > weps) then
         if (abs(wd) > deps) then
            tn = tanh(0.5*wt)
            tn2 = tn*tn;
            tn2i = 1.0/(1.0 - tn2)
            csh = (1.0 + tn2)*tn2i
            snh = (tn + tn)*tn2i
c second order taylor series approximation for wd
         else
            wt2 = wd*wd
            wd = -wd*(1.0 + s1*wt2)
            wt2 = 1.0 + c1*wt2
c third order taylor series approximation for wd
c           wd = -wd*(1.0 + s1*wt2*(1.0 + s2*wt2))
c           wt2 = 1.0d+ c1*wt2*(1.0 + c2*wt2)
            f = csh*wt2 + snh*wd
            snh = snh*wt2 + csh*wd
            csh = f
         endif
c calculate special functions
         cshc = csh - 1.0
         snhc = snh*ali
         cshd = cshc*ali2
         snhd = (snhc - tau)*ali2
c fourth order taylor series approximation for wt
      else
         wt2 = wt*wt
         cshd = c1*(1.0 + c2*wt2*(1.0 + c3*wt2))*t2
         snhd = s1*(1.0 + s2*wt2*(1.0 + s3*wt2))*t3
c fifth order taylor series approximation for wt
c        cshd = c1*(1.0 + c2*wt2*(1.0 + c3*wt2*(1.0 + c4*wt2)))*t2
c        snhd = s1*(1.0 + s2*wt2*(1.0 + s3*wt2*(1.0 + s4*wt2)))*t3
         cshc = al2*cshd
         snhc = tau + al2*snhd
         csh = cshc + 1.0
         snh = al*snhc
      endif
c gam = gamma(tau)
      gam = gam0 + (ea*snc - eb*csd) + (ec*snhc + ed*cshd)
      fpi = 1.0/gam
c calculate time expression whose root we seek
      f = gam0*tau - (ea*csd + eb*snd) + (ec*cshd + ed*snhd) - dt
c newton's quadratic method
      ft = f*fpi
c either add Halley's optional cubic correction
c fpp = dgamma/dtau
c     fpp = (ea*cs + eb*snc) + (ec*csh + ed*snhc)
c     ft = ft/(1.0d0 - 0.5d0*ft*fpp*fpi)
c or add Householder's optional quartic correction
c fpp = dgamma/dtau
c     fpp = (ea*cs + eb*snc) + (ec*csh + ed*snhc)
c fp = d2gamma/dtau2
c     fp = (eb*cs - ea*w*sn) + (ec*al*snh + ed*csh)
c     wt2 = ft*fpp*fpi
c     ft = ft*(1.0d0 - 0.5d0*wt2)/(1.0d0 - wt2 + sixth*(ft*ft)*fp*fpi)
c update tau: ft = -delta tau
      tau = tau - ft
      if (abs(f) < teps) go to 40
   30 continue
c
c convergence failure
   40 continue
      if (ii.gt.imax) then
         write (*,*) i,'tau error:f,ft=',f,ft
      endif
c update gamma
c fpp = dgamma/dt
      if (fpp.eq.0.0) fpp = (ea*cs + eb*snc) + (ec*csh + ed*snhc)
c first order taylor series
      gam = gam - fpp*ft
c second order taylor series
c fpp = dgamma/dtau
c     if (fp.eq.0.0) fpp = (ea*cs + eb*snc) + (ec*csh + ed*snhc)
c fp = d2gamma/dtau2
c     if (fp.eq.0.0) fp = (eb*cs - ea*w*sn) + (ec*al*snh + ed*csh)
c     gam = gam - (fpp - half*fp*ft)*ft
c
c update sine/cosine functions
c first order taylor series
      wd = w*ft
      f = sn*wd
      csc = csc + f
      snc = snc - cs*ft
      f = cs + f
      sn = sn - cs*wd
      cs = f
c second order taylor series
c     wd = w*ft
c     wt2 = wd*wd
c     f = -ft*(1.0d0 - sixth*wt2)
c     wd = half*wt2
c     wt2 = 1.0d0 - wd
c     sc = csc*wt2 - wd
c     wd = w*f
c     sc = csc - sn*wd
c     snc = snc*wt2 + cs*f
c     f = cs*wt2 - sn*wd
c     sn = sn*wt2 + cs*wd
c     cs = f
c update momenta
      up = up + ep*dt
      ut = ut + om2t*al2*dt + (omega*wi*fa*csc + fb*snc)
      ud = om2t*omega*gam + fa*sn - fc*cs
c calculate inverse gamma and average kinetic energy
      gami = 1.0/gam
      p2 = up*up + ut*ut + ud*ud
      sum1 = sum1 + dble(0.5*(wk + gami*p2/(1.0 + gami)))
c sanity check
c     erm = max(erm,sqrt(1.0d0 + p2*ci2))
c translate momenta to cartesian coordinates
      px = up*ox + ut*tx + ud*dx
      py = up*oy + ut*ty + ud*dy
      pz = up*oz + ut*tz + ud*dz
      part(2,i) = px
      part(3,i) = py
      part(4,i) = pz
c new position
      dtg = dtc*gami
      dx = part(1,i) + px*dtg
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
c set new position
      part(1,i) = dx
   50 continue
c normalize kinetic energy
      ek = ek + sum1
c sanity check
c     write (*,*) 'gamma sanity check=',erm
      return
      end
c-----------------------------------------------------------------------
      subroutine DPOST1GL(part,q,sctx,qm,nop,idimp,nx,nxhm,nxvh)
c for 1d code, this subroutine calculates particle charge density
c using gridless spectral version, periodic boundaries
c baseline scalar version
c 4*nxhm flops/particle
c input: all, output: q, sctx
c charge density is calculated from the expression:
c q(n) = sum(qm*exp(-sqrt(-1)*2*n*pi*x/nx))
c part(1,n) = position x of particle n
c q(n) = charge density at fourier grid point n
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nx = system length in x direction
c nxhm = number of fourier modes kept
c nxvh = dimension of charge array, must be >= nxhm
      implicit none
      integer nop, idimp, nx, nxhm, nxvh
      real part, qm
      complex q, sctx
      dimension part(idimp,nop), q(nxvh), sctx(nxvh)
c local data
      integer i, j
      real qmn, dnx, dkx
      dnx = 6.28318530717959/real(nx)
      qmn = qm/real(nx)
c find fourier components
      do 30 i = 1, nop
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c mode numbers 0 < kx < nxhm
      do 20 j = 2, nxhm
      q(j) = q(j) + qmn*sctx(j-1)
   20 continue
c special case to match conventional PIC code
      if (nxhm.eq.(nx/2)) then
         q(1) = cmplx(real(q(1))+qmn,aimag(q(1))+qmn*real(sctx(nxhm)))
      else
         q(1) = cmplx(real(q(1))+qmn,0.0)
      endif
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DJPOST1GL(part,cu,sctx,cux0,qm,dt,nop,idimp,nx,nxhm,   
     1nxvh)
c for 1-2/2d code, this subroutine calculates particle current density
c using gridless spectral version, periodic boundaries
c in addition, particle positions are advanced a half time-step
c baseline scalar version
c 8*nxhm flops/particle + 2 flops/particle, 2*nxhm loads and stores
c input: all, output: part, cu
c current density is calculated from the expression:
c cu(i,n) = sum(qmi*exp(-sqrt(-1)*2*n*pi*x/nx))
c where qmi = qm*vi, where i = x,y,z
c part(1,n) = position x of particle n
c part(2,n) = x velocity of particle n
c part(3,n) = y velocity of particle n
c part(4,n) = z velocity of particle n
c cu(i,n) = charge density at fourier grid point n
c sctx = scratch array for sines and cosines
c cux0 = net current in x direction
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 4
c nx = system length in x direction
c nxhm = number of fourier modes kept
c nxvh = second dimension of current array, must be >= nxhm
      implicit none
      integer nop, idimp, nx, nxhm, nxvh
      real part, cux0, qm, dt
      complex cu, sctx
      dimension part(idimp,nop), cu(2,nxvh), sctx(nxvh)
c local data
      integer i, j
      real zero, anx, qmn, dnx, dkx, at1
      real vx, vy, vz, dx
      complex zt1
      double precision sum1
      zero = 0.0
      anx = real(nx)
      dnx = 6.28318530717959/anx
      qmn = qm/anx
      sum1 = 0.0d0
c find fourier components
      do 30 i = 1, nop
      vx = qmn*part(2,i)
      vy = qmn*part(3,i)
      vz = qmn*part(4,i)
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c mode numbers 0 < kx < nxhm
      do 20 j = 2, nxhm
      zt1 = sctx(j-1)
      cu(1,j) = cu(1,j) + vy*zt1
      cu(2,j) = cu(2,j) + vz*zt1
   20 continue
c special case to match conventional PIC code
      if (nxhm.eq.(nx/2)) then
         at1 = real(sctx(nxhm))
         cu(1,1) = cmplx(real(cu(1,1))+vy,aimag(cu(1,1))+vy*at1)
         cu(2,1) = cmplx(real(cu(2,1))+vz,aimag(cu(2,1))+vz*at1)
      else
         cu(1,1) = cmplx(real(cu(1,1))+vy,0.0)
         cu(2,1) = cmplx(real(cu(2,1))+vz,0.0)
      endif
c sum current in x direction
      sum1 = sum1 + vx
c advance position half a time-step
      dx = part(1,i) + part(2,i)*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
c set new position
      part(1,i) = dx
   30 continue
      cux0 = cux0 + sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine RDJPOST1GL(part,cu,sctx,cux0,qm,dt,ci,nop,idimp,nx,nxhm
     1,nxvh)
c for 1-2/2d code, this subroutine calculates particle current density
c using gridless spectral version for relativistic particles,
c periodic boundaries
c in addition, particle positions are advanced a half time-step
c baseline scalar version
c 8*nxhm flops/particle + 12 flops/particle, 2*nxhm loads and stores
c input: all, output: part, sctx, cu
c current density is calculated from the expression:
c cu(i,n) = sum(qci*exp(-sqrt(-1)*2*n*pi*x/nx))
c and qci = qm*pi*gami, where i = y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n
c part(2,n) = x momentum of particle n
c part(3,n) = y momentum of particle n
c part(4,n) = z momentum of particle n
c cu(i,n) = charge density at fourier grid point n
c sctx = scratch array for sines and cosines
c cux0 = net current in x direction
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 4
c nx = system length in x direction
c nxhm = number of fourier modes kept
c nxvh = first dimension of current array, must be >= nxhm
      implicit none
      integer nop, idimp, nx, nxhm, nxvh
      real part, cux0, qm, dt, ci
      complex cu, sctx
      dimension part(idimp,nop), cu(2,nxvh), sctx(nxvh)
c local data
      integer i, j
      real zero, anx, qmn, dnx, dkx, at1
      real vx, vy, vz, dx
      real ci2, p2, gami, px
      complex zt1
      double precision sum1
      zero = 0.0
      anx = real(nx)
      dnx = 6.28318530717959/anx
      qmn = qm/anx
      ci2 = ci*ci
      sum1 = 0.0d0
c find fourier components
      do 30 i = 1, nop
c find inverse gamma
      vx = part(2,i)
      vy = part(3,i)
      vz = part(4,i)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
      px = vx*gami
      vx = qmn*vx*gami
      vy = qmn*vy*gami
      vz = qmn*vz*gami
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c mode numbers 0 < kx < nxhm
      do 20 j = 2, nxhm
      zt1 = sctx(j-1)
      cu(1,j) = cu(1,j) + vy*zt1
      cu(2,j) = cu(2,j) + vz*zt1
   20 continue
c special case to match conventional PIC code
      if (nxhm.eq.(nx/2)) then
         at1 = real(sctx(nxhm))
         cu(1,1) = cmplx(real(cu(1,1))+vy,aimag(cu(1,1))+vy*at1)
         cu(2,1) = cmplx(real(cu(2,1))+vz,aimag(cu(2,1))+vz*at1)
      else
         cu(1,1) = cmplx(real(cu(1,1))+vy,0.0)
         cu(2,1) = cmplx(real(cu(2,1))+vz,0.0)
      endif
c sum current in x direction
      sum1 = sum1 + vx
c advance position half a time-step
      dx = part(1,i) + px*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
c set new position
      part(1,i) = dx
   30 continue
      cux0 = cux0 + sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine D2JPOST1GL(part,dcu,sctx,qm,nop,idimp,nx,nxhm,nxvh)
c for 1-2/2d code, this subroutine calculates the time derivative of the
c particle current density for force-free particles,
c used to initialize the darwin tranverse electric field
c using gridless spectral version, periodic boundaries
c baseline scalar version
c 2 + 21*nxhm flops per particle + 4 + 2*nxhm loads and stores
c plus nxhm sines and cosines/particle
c input: all, output: part, dcu, sctx
c derivative of current density is calculated from the expression:
c dcu(i,n) = sum(qdi*exp(-sqrt(-1)*kn*x))
c where qdi = qm*vi*d, where i = y,z,
c kn = 2*n*pi/nx, k.vi = kn*vx, d = -sqrt(-1)*k.vi
c part(1,n) = position x of particle n
c part(2,n) = x velocity of particle n
c part(3,n) = y velocity of particle n
c part(4,n) = z velocity of particle n
c cu(i,n) = current density at fourier grid point n
c dcu(i,n) = derivative of current at fourier grid point n
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nx = system length in x direction
c nxhm = number of fourier modes kept in x
c nxvh = first dimension of current array, must be >= nxhm
      implicit none
      integer nop, idimp, nx, nxhm, nxvh
      real part, qm
      complex dcu, sctx
      dimension part(idimp,nop), dcu(2,nxvh)
      dimension sctx(nxvh)
c local data
      integer i, j
      real qmn, dnx, dkx
      real x, ux, vy, vz, dp
      complex zt1
      dnx = 6.28318530717959/real(nx)
      qmn = qm/real(nx)
c find fourier components
      do 30 i = 1, nop
      x = part(1,i)
      ux = part(2,i)
      vy = qmn*part(3,i)
      vz = qmn*part(4,i)
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*x
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c mode numbers 0 < kx < nxhm
      do 20 j = 2, nxhm
      dkx = dnx*real(j-1)
      dp = dkx*ux
c deposit current derivative
      zt1 = sctx(j-1)*cmplx(0.0,-dp)
      dcu(1,j) = dcu(1,j) + vy*zt1
      dcu(2,j) = dcu(2,j) + vz*zt1
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine RD2JPOST1GL(part,dcu,sctx,qm,ci,nop,idimp,nx,nxhm,nxvh)
c for 1-2/2d code, this subroutine calculates the time derivative of the
c particle current density for force-free particles,
c used to initialize the darwin tranverse electric field
c using gridless spectral version for relativistic particles,
c periodic boundaries
c baseline scalar version
c 11 + 23*nxhm flops per particle + 4 + 2*nxhm loads and stores
c plus nxhm sines and cosines divides and sqrts/particle
c input: all, output: part, dcu, sctx
c derivative of current density is calculated from the expression:
c dcu(i,n,m) = sum(qdi*exp(-sqrt(-1)*kn*x))
c and qdi = qm*vi*d, where i = y,z,
c where vi = pi*gami, gami = 1./sqrt(1.+sum(pi**2)*ci*ci) 
c kn = 2*n*pi/nx, k.vi = kn*vx and d = -sqrt(-1)*k.vi
c part(1,n) = position x of particle n
c part(2,n) = x momentum of particle n
c part(3,n) = y momentum of particle n
c part(4,n) = z momentum of particle n
c cu(i,n) = current density at fourier grid point n
c dcu(i,n) = derivative of current at fourier grid point n
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 4
c nx = system length in x direction
c nxhm = number of fourier modes kept in x
c nxvh = first dimension of current array, must be >= nxhm
      implicit none
      integer nop, idimp, nx, nxhm, nxvh
      real part, qm, ci
      complex dcu, sctx
      dimension part(idimp,nop), dcu(2,nxvh)
      dimension sctx(nxvh)
c local data
      integer i, j
      real qmn, dnx, dkx
      real x, ux, vy, vz, ci2, p2, gami, dp
      complex zt1
      dnx = 6.28318530717959/real(nx)
      qmn = qm/real(nx)
      ci2 = ci*ci
c find fourier components
      do 30 i = 1, nop
      x = part(1,i)
      ux = part(2,i)
c find inverse gamma
      vy = part(3,i)
      vz = part(4,i)
      p2 = ux*ux + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
      vy = qmn*(vy*gami)
      vz = qmn*(vz*gami)
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*x
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c mode numbers 0 < kx < nxhm
      do 20 j = 2, nxhm
      dkx = dnx*real(j-1)
      dp = dkx*ux*gami
c deposit current derivative
      zt1 = sctx(j-1)*cmplx(0.0,-dp)
      dcu(1,j) = dcu(1,j) + vy*zt1
      dcu(2,j) = dcu(2,j) + vz*zt1
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DSJPOST1GL(part,cu,dcu,sctx,qm,ci,nop,idimp,nx,nxhm,   
     1nxvh)
c for 1-2/2d code, this subroutine calculates a rescaled particle
c current density and its time derivative for force-free particles,
c used to initialize the electromagnetic fields
c using gridless spectral version for non-relativistic particles,
c periodic boundaries
c baseline scalar version
c 33*nxhm flops and 2*nxhm divides per particle
c 4*nxhm + 4 loads and 4*nxhm stores
c plus nxhm sines and cosines/particle
c input: all, output: part, cu, dcu, sctx
c current density is calculated from the expression:
c cu(i,n) = sum(qci*exp(-sqrt(-1)*kn*x))
c derivative of current density is calculated from the expression:
c dcu(i,n) = sum(qdi*exp(-sqrt(-1)*kn*x))
c where qci = qm*vi*s, qdi = qci*d, where i = y,z,
c s = 1.0/(1.0 - (k.vi*ci)**2)/(k*k))
c kn = 2*n*pi/nx, k*k = kn*kn
c k.vi = kn*vx and d = -sqrt(-1)*k.vi
c if nxhm=nx/2, nx/2 mode is included, in a format compatible with the
c real to complex FFT used in traditional PIC
c part(1,n) = position x of particle n
c part(2,n) = x velocity of particle n
c part(3,n) = y velocity of particle n
c part(4,n) = z velocity of particle n
c cu(i,n) = current density at fourier grid point n
c dcu(i,n) = derivative of current at fourier grid point n
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 5
c nx = system length in x direction
c nxhm = number of fourier modes kept in x
c nxvh = first dimension of current array, must be >= nxhm
      implicit none
      integer nop, idimp, nx, nxhm, nxvh
      real part, qm, ci
      complex cu, dcu, sctx
      dimension part(idimp,nop), cu(2,nxvh), dcu(2,nxvh)
      dimension sctx(nxvh)
c local data
      integer i, j
      real qmn, dnx, dkx, dk2, at1
      real x, ux, vy, vz, ci2, sp, dp
      complex zt1
      dnx = 6.28318530717959/real(nx)
      qmn = qm/real(nx)
      ci2 = ci*ci
c find fourier components
      do 30 i = 1, nop
      x = part(1,i)
      ux = part(2,i)
      vy = qmn*part(3,i)
      vz = qmn*part(4,i)
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*x
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c mode numbers 0 < kx < nxhm
      do 20 j = 2, nxhm
      dkx = dnx*real(j-1)
      dk2 = dkx*dkx
      dp = dkx*ux
      sp = 1.0/(1.0 - dp*dp*ci2/dk2)
c deposit current
      zt1 = sctx(j-1)*sp
      cu(1,j) = cu(1,j) + vy*zt1
      cu(2,j) = cu(2,j) + vz*zt1
c deposit current derivative
      zt1 = zt1*cmplx(0.0,-dp)
      dcu(1,j) = dcu(1,j) + vy*zt1
      dcu(2,j) = dcu(2,j) + vz*zt1
   20 continue
c special case to match conventional PIC code
c deposit current
      if (nxhm.eq.(nx/2)) then
         dkx = dnx*real(nxhm)
         dk2 = dkx*dkx
         dp = dkx*ux
         sp = 1.0/(1.0 - dp*dp*ci2/dk2)
         at1 = real(sctx(nxhm))*sp
         cu(1,1) = cmplx(real(cu(1,1))+vy,aimag(cu(1,1))+vy*at1)
         cu(2,1) = cmplx(real(cu(2,1))+vz,aimag(cu(2,1))+vz*at1)
      else
         cu(1,1) = cmplx(real(cu(1,1))+vy,0.0)
         cu(2,1) = cmplx(real(cu(2,1))+vz,0.0)
      endif
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine RDSJPOST1GL(part,cu,dcu,sctx,qm,ci,nop,idimp,nx,nxhm,  
     1nxvh)
c for 1-2/2d code, this subroutine calculates a rescaled particle
c current density and its time derivative for force-free particles,
c used to initialize the electromagnetic fields
c using gridless spectral version for relativistic particles,
c periodic boundaries
c baseline scalar version
c 44*nxhm flops, 3*nxhm divides, and nxhm sqrts per particle
c 4*nxhm + 4 loads and 4*nxhm stores
c plus nxhm sines and cosines/particle
c input: all, output: part, cu, dcu, sctx
c current density is calculated from the expression:
c cu(i,n,m) = sum(qci*exp(-sqrt(-1)*kn*x))
c derivative of current density is calculated from the expression:
c dcu(i,n,m) = sum(qdi*exp(-sqrt(-1)*kn*x))
c where qci = qm*vi*s, qdi = qci*d, where i = y,z,
c gam2 = 1.+sum(pi**2)*ci*ci), gami = 1./sqrt(gam2), vi = pi*gami, 
c s = gam2/(1.0 + ((sum(pi**2) - (k.pi)**2/(k*k))*ci*ci
c kn = 2*n*pi/nx, k*k = kn*kn,
c k.pi = kn*px and d = -sqrt(-1)*k.vi
c if nxhm=nx/2, nx/2 mode is included, in a format compatible with the
c real to complex FFT used in traditional PIC
c part(1,n) = position x of particle n
c part(2,n) = x momentum of particle n
c part(3,n) = y momentum of particle n
c part(4,n) = z momentum of particle n
c cu(i,n) = current density at fourier grid point n
c dcu(i,n) = derivative of current at fourier grid point n
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 4
c nx = system length in x direction
c nxhm = number of fourier modes kept in x
c nxvh = first dimension of current array, must be >= nxhm
      implicit none
      integer nop, idimp, nx, nxhm, nxvh
      real part, qm, ci
      complex cu, dcu, sctx
      dimension part(idimp,nop), cu(2,nxvh), dcu(2,nxvh)
      dimension sctx(nxvh)
c local data
      integer i, j
      real qmn, dnx, dkx, dk2, at1
      real x, ux, vy, vz, ci2, p2, gam2, gami
      real sp, dp
      complex zt1
      dnx = 6.28318530717959/real(nx)
      qmn = qm/real(nx)
      ci2 = ci*ci
c find fourier components
      do 30 i = 1, nop
      x = part(1,i)
      ux = part(2,i)
c find inverse gamma
      vy = part(3,i)
      vz = part(4,i)
      p2 = ux*ux + vy*vy + vz*vz
      gam2 = 1.0 + p2*ci2
      gami = 1.0/sqrt(gam2)
      vy = qmn*(vy*gami)
      vz = qmn*(vz*gami)
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*x
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c mode numbers 0 < kx < nxhm
      do 20 j = 2, nxhm
      dkx = dnx*real(j-1)
      dk2 = dkx*dkx
      dp = dkx*ux
      sp = dp*dp/dk2
      dp = dp*gami
      sp = gam2/(1.0 + (p2 - sp)*ci2)
c deposit current
      zt1 = sctx(j-1)*sp
      cu(1,j) = cu(1,j) + vy*zt1
      cu(2,j) = cu(2,j) + vz*zt1
c deposit current derivative
      zt1 = zt1*cmplx(0.0,-dp)
      dcu(1,j) = dcu(1,j) + vy*zt1
      dcu(2,j) = dcu(2,j) + vz*zt1
   20 continue
c special case to match conventional PIC code
c deposit current
      if (nxhm.eq.(nx/2)) then
         dkx = dnx*real(nxhm)
         dk2 = dkx*dkx
         dp = dkx*ux
         sp = dp*dp/dk2
         dp = dp*gami
         sp = gam2/(1.0 + (p2 - sp)*ci2)
         at1 = real(sctx(nxhm))*sp
         cu(1,1) = cmplx(real(cu(1,1))+vy,aimag(cu(1,1))+vy*at1)
         cu(2,1) = cmplx(real(cu(2,1))+vz,aimag(cu(2,1))+vz*at1)
      else
         cu(1,1) = cmplx(real(cu(1,1))+vy,0.0)
         cu(2,1) = cmplx(real(cu(2,1))+vz,0.0)
      endif
   30 continue
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
      subroutine IBPOIS13GL(cu,byz,ffc,ci,wm,nx,nxhm,nxvh)
c this subroutine solves 1-2/2d poisson's equation in fourier space for
c magnetic field, with periodic boundary conditions, 
c for gridless spectral code
c input: cu,ffc,ci,nx,nxv, output: byz,wm
c approximate flop count is: 29*nxhm
c the magnetic field is calculated using the equations:
c by(kx) = -ci*ci*sqrt(-1)*g(kx)*kx*cuz(kx),
c bz(kx) = ci*ci*sqrt(-1)*g(kx)*kx*cuy(kx),
c where kx = 2pi*j/nx, and j = fourier mode number,
c g(kx) = (affp/kx**2)*s(kx),
c s(kx) = exp(-((kx*ax)**2)/2)
c cu(i,j) = complex current density for fourier mode (j-1)
c byz(i,j) = i component of complex magnetic field
c all for fourier mode (j-1)
c aimag(ffc(j)) = finite-size particle shape factor s
c for fourier mode (j-1)
c real(ffc(j)) = potential green's function g
c for fourier mode (j-1)
c ci = reciprocal of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*sum((affp/kx**2)*ci*ci*|cu(kx)*s(kx)|**2), where
c affp = normalization constant = nx/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx = system length in x direction
c nxhm = number of fourier modes kept
c nxvh = second dimension of field arrays, must be >= nxhm
      implicit none
      integer nx, nxhm, nxvh
      real ci, wm
      complex cu, byz, ffc
      dimension cu(2,nxvh), byz(2,nxvh), ffc(nxhm)
c local data
      integer j
      real dnx, ci2, at1, at2
      complex zero, zt1, zt2
      double precision wp
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nxhm
      do 10 j = 2, nxhm
      at1 = ci2*real(ffc(j))
      at2 = dnx*real(j - 1)*at1
      at1 = at1*aimag(ffc(j))
      zt1 = cmplx(-aimag(cu(2,j)),real(cu(2,j)))
      zt2 = cmplx(-aimag(cu(1,j)),real(cu(1,j)))
      byz(1,j) = -at2*zt1
      byz(2,j) = at2*zt2
      wp = wp + at1*(cu(1,j)*conjg(cu(1,j)) + cu(2,j)*conjg(cu(2,j)))
   10 continue
c mode number kx = 0
      byz(1,1) = zero
      byz(2,1) = zero
      wm = real(nx)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine AMAXWEL1GL(eyz,byz,cu,ffc,ci,dt,wf,wm,nx,nxhm,nxvh)
c this subroutine solves 1d maxwell's equation in fourier space for
c transverse electric and magnetic fields with periodic boundary
c conditions, using an analytic scheme due to Irving Haber
c for gridless spectral code.
c input: all, output: wf, wm, eyz, byz
c approximate flop count is: 82*nxhm
c equations being solved are:
c (c*Bn)/dt = -ick X En and
c En/dt = ick X (c*Bn) - JTn
c Note that the normalization of the E and B fields differ:
c B is normalized to the dimensionless cyclotron frequency.
c solutions are given by:
c En+1 = C*En + iS*(k X (c*Bn)) - S*JTn+1/2/c
c c*Bn+1 = C*(c*Bn) - iS*(k X En)/c + iS*T*(k X JTn+1/2)
c where En = input eyz, En+1 = output eyz
c Bn = input byz, Bn+1 = output byz, JTn+1/2 = affp*cu*s(kx)
c C = cos(k*c*dt),  S = sin(k*c*dt)/k, T = tan(k*c*dt/2)/kc
c kx = 2pi*j/nx, k = sqrt(dkx*dkx), c = 1.0/ci
c and s(kx) = exp(-(kx*ax)**2/2)
c cu(i,j) = complex current density
c eyz(i,j) = complex transverse electric field
c byz(i,j) = complex magnetic field
c for component i, all for fourier mode (j-1)
c real(ffc(1)) = affp = normalization constant = nx/np,
c where np=number of particles
c aimag(ffc(j) = finite-size particle shape factor s,
c s(kx) = exp(-(kx*ax)**2/2), for fourier mode (j-1)
c ci = reciprocal of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*sum((1/affp)*|eyz(kx)|**2)
c magnetic field energy is also calculated, using
c wm = nx*sum((c2/affp)*|byz(kx)|**2), where c2 = 1./(ci*ci)
c nx = system length in x direction
c nxhm = number of fourier modes kept in x
c nxvh = first dimension of field arrays, must be >= nxhm
      implicit none
      integer nx, nxhm, nxvh
      real ci, dt, wf, wm
      complex eyz, byz, cu, ffc
      dimension eyz(2,nxvh), byz(2,nxvh), cu(2,nxvh), ffc(nxhm)
c local data
      integer j
      real dnx, dth, cc, cdth, affp, anorm, dkx, t2
      real t, c, s, sc, afs, aft
      complex zero, zt1, zt2, zt5, zt6, zt7, zt8
      double precision wp, ws
      if (ci.le.0.0) return
      dnx = 6.28318530717959/real(nx)
      dth = 0.5*dt
      cc = 1.0/ci
      cdth = cc*dth
      affp = real(ffc(1))
      zero = cmplx(0.0,0.0)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
c calculate the electromagnetic fields
      do 10 j = 2, nxhm
      dkx = dnx*real(j - 1)
      aft = affp*aimag(ffc(j))*ci
      t2 = sqrt(dkx*dkx)
      t = t2*cdth
      t2 = 1.0/t2
      s = t + t
      c = cos(s)
      s = sin(s)*t2
      t = tan(t)*t2
      afs = s*aft
      sc = s*cc
      aft = aft*t
      s = s*ci
c calculate iB
      zt1 = cmplx(-aimag(byz(2,j)),real(byz(2,j)))
      zt2 = cmplx(-aimag(byz(1,j)),real(byz(1,j)))
c update electric field
      zt5 = c*eyz(1,j) - sc*(dkx*zt1) - afs*cu(1,j)
      zt6 = c*eyz(2,j) + sc*(dkx*zt2) - afs*cu(2,j)
c calculate iE
      zt1 = cmplx(-aimag(eyz(2,j)),real(eyz(2,j)))
      zt2 = cmplx(-aimag(eyz(1,j)),real(eyz(1,j)))
c store electric field and calculate energy
      eyz(1,j) = zt5
      eyz(2,j) = zt6
      ws = ws + anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
c calculate ijperp
      zt7 = cmplx(-aimag(cu(2,j)),real(cu(2,j)))
      zt8 = cmplx(-aimag(cu(1,j)),real(cu(1,j)))
c update magnetic field 
      zt5 = c*byz(1,j) + s*dkx*(zt1 - aft*zt7)
      zt6 = c*byz(2,j) - s*dkx*(zt2 - aft*zt8)
c store magnetic field and calculate energy
      byz(1,j) = zt5
      byz(2,j) = zt6
      wp = wp + anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
   10 continue
c mode number kx = 0
c     afs = affp*aimag(ffc(1))*dt
c     eyz(1,1) = eyz(1,1) - afs*cu(1,1)
c     eyz(2,1) = eyz(2,1) - afs*cu(2,1)
      byz(1,1) = zero
      byz(2,1) = zero
      eyz(1,1) = zero
      eyz(2,1) = zero
      wf = real(nx)*ws
      wm = real(nx)*(cc*cc)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine EMFIELD1GL(fxyz,fx,eyz,ffc,nxhm,nxvh)
c this subroutine merges complex vector fields
c includes additional smoothing, for gridless spectral code
      implicit none
      integer nxhm, nxvh
      complex fxyz, fx, eyz, ffc
      dimension fxyz(3,nxvh), fx(nxvh), eyz(2,nxvh)
      dimension ffc(nxhm)
c local data
      integer j
      real at1
c add the fields
      do 10 j = 1, nxhm
      at1 = aimag(ffc(j))
      fxyz(1,j) = fx(j)
      fxyz(2,j) = eyz(1,j)*at1
      fxyz(3,j) = eyz(2,j)*at1
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine BMFIELD1GL(fyz,eyz,ffc,nxhm,nxvh)
c this subroutine copies complex vector fields
c includes additional smoothing, for gridless spectral code
      implicit none
      integer nxhm, nxvh
      complex fyz, eyz, ffc
      dimension fyz(2,nxvh), eyz(2,nxvh)
      dimension ffc(nxhm)
c local data
      integer j
      real at1
      do 10 j = 1, nxhm
      at1 = aimag(ffc(j))
      fyz(1,j) = eyz(1,j)*at1
      fyz(2,j) = eyz(2,j)*at1
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine EPOIS13GL(dcu,eyz,isign,ffe,ax,affp,wp0,ci,wf,nx,nxhm, 
     1nxvh)
c this subroutine solves 1-2/2d poisson's equation in fourier space for
c transverse electric field (or convolution of transverse electric field
c over particle shape), with periodic boundary conditions,
c for gridless spectral code
c using algorithm described in J. Busnardo-Neto, P. L. Pritchett,
c A. T. Lin, and J. M. Dawson, J. Computational Phys. 23, 300 (1977).
c for isign = 0, input: isign,ax,affp,wp0,nx,nxvh, output:ffe
c for isign =/ 0, input: dcu,ffe,isign,ci,nx,nxvh,nxhd, output: eyz,wf
c approximate flop count is: 25*nxhm
c if isign = 0, form factor array is prepared
c if isign = -1, smoothed transverse electric field is calculated
c using the equation:
c ey(kx) = -ci*ci*g(kx)*dcuy(kx)*s(kx)
c ez(kx) = -ci*ci*g(kx)*dcuz(kx)*s(kx)
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c g(kx) = (affp/(kx**2+wp0*ci2*s(kx)**2))*s(kx),
c s(kx) = exp(-((kx*ax)**2+)/2)
c if isign = 1, unsmoothed transverse electric field is calculated
c using the equation:
c ey(kx) = -ci*ci*g(kx)*dcuy(kx)
c ez(kx) = -ci*ci*g(kx)*dcuz(kx)
c dcu(i,j) = transverse part of complex derivative of current for
c fourier mode (j-1)
c eyz(1,j) = y component of complex transverse electric field
c eyz(2,j) = z component of complex transverse electric field
c all for fourier mode (j-1)
c aimag(ffe(j)) = finite-size particle shape factor s
c for fourier mode (j-1)
c real(ffe(j)) = potential green's function g
c for fourier mode (j-1)
c ax = half-width of particle in x direction
c affp = normalization constant = nx/np, where np=number of particles
c wp0 = normalized total plasma frequency squared
c ci = reciprocal of velocity of light
c transverse electric field energy is also calculated, using
c wf = nx*sum((affp/(kx**2*ci*ci)**2)*|dcu(kx)*s(kx)|**2)
c this expression is valid only if the derivative of current is
c divergence-free
c nx = system length in x direction
c nxhm = number of fourier modes kept
c nxvh = second dimension of field arrays, must be >= nxhm
      implicit none
      integer isign, nx, nxhm, nxvh
      real ax, affp, wp0, ci, wf
      complex dcu, eyz, ffe
      dimension dcu(2,nxvh), eyz(2,nxvh)
      dimension ffe(nxvh)
c local data
      integer j
      real dnx, ci2, wpc, dkx, at1, at2
      complex zero
      double precision wp
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
      if (isign.ne.0) go to 20
      wpc = wp0*ci2
c prepare form factor array
      do 10 j = 2, nxhm
      dkx = dnx*real(j - 1)
      at2 = exp(-.5*(dkx*ax)**2)
      ffe(j) = cmplx(affp*at2/(dkx*dkx+ wpc*at2*at2),at2)
   10 continue
      ffe(1) = cmplx(affp,1.0)
      return
c calculate smoothed transverse electric field and sum field energy
   20 if (isign.gt.0) go to 40
      wp = 0.0d0
c mode numbers 0 < kx < nxhm
      do 30 j = 2, nxhm
      at2 = -ci2*real(ffe(j))
      at1 = at2*aimag(ffe(j))
      at2 = at2*at2
      eyz(1,j) = at1*dcu(1,j)
      eyz(2,j) = at1*dcu(2,j)
      wp = wp + at2*(dcu(1,j)*conjg(dcu(1,j)) + dcu(2,j)*conjg(dcu(2,j))
     1)
   30 continue
c mode number kx = 0
      eyz(1,1) = zero
      eyz(2,1) = zero
      wf = real(nx)*wp/real(ffe(1))
      return
c calculate unsmoothed transverse electric field and sum field energy
   40 wp = 0.0d0
c mode numbers 0 < kx < nxhm
      do 50 j = 2, nxhm
      at2 = -ci2*real(ffe(j))
      at1 = at2*at2
      eyz(1,j) = at2*dcu(1,j)
      eyz(2,j) = at2*dcu(2,j)
      wp = wp + at1*(dcu(1,j)*conjg(dcu(1,j)) + dcu(2,j)*conjg(dcu(2,j))
     1)
   50 continue
c mode number kx = 0
      eyz(1,1) = zero
      eyz(2,1) = zero
      wf = real(nx)*wp/real(ffe(1))
      return
      end
c-----------------------------------------------------------------------
      subroutine EVFIELD13GL(fxy,gxy,sctx,xp,nx,nxhm,nxvh)
c for 1-2/2d code, this subroutine calculates real space 1d field gxy at
c location xp, from complex fourier co-efficients fxy
c input: all, output: gxy, sctx
c equations used are:
c gxy(1) = sum(fxy(1,n)*exp(sqrt(-1)*2*n*pi*xp/nx)
c gxy(2) = sum(fxy(2,n)*exp(sqrt(-1)*2*n*pi*xp/nx)
c gxy(3) = sum(fxy(3,n)*exp(sqrt(-1)*2*n*pi*xp/nx)
c fxy(1,n) = x component of force/charge at fourier mode n-1
c fxy(2,n) = y component of force/charge at fourier mode n-1
c fxy(3,n) = z component of force/charge at fourier mode n-1
c if nxhm=nx/2, nx/2 mode is included, in a format compatible with the
c real to complex FFT used in traditional PIC
c gxy(1) = x component of force/charge at location xp
c gxy(2) = y component of force/charge at location xp
c gxy(3) = z component of force/charge at location xp
c sctx = scratch array for sines and cosines
c xp = position to evaluate field
c nx = system length in x direction
c nxhm = number of fourier modes kept in x
c nxvh = second dimension of vector field arrays, must be >= nxhm
      implicit none
      integer nx, nxhm, nxvh
      real xp
      real gxy
      complex fxy, sctx
      dimension gxy(3)
      dimension fxy(3,nxvh), sctx(nxvh)
c local data
      integer j
      real dnx, dkx, at1
      double precision ex, ey, ez
      complex zt1
      dnx = 6.28318530717959/real(nx)
c find field
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*xp
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
c mode numbers 0 < kx < nxhm
      do 20 j = 2, nxhm
      zt1 = sctx(j-1)
      ex = ex + real(fxy(1,j)*zt1)
      ey = ey + real(fxy(2,j)*zt1)
      ez = ez + real(fxy(3,j)*zt1)
   20 continue
      ex = 2.0d0*ex + real(fxy(1,1))
      ey = 2.0d0*ey + real(fxy(2,1))
      ez = 2.0d0*ez + real(fxy(3,1))
c special case to match conventional PIC code
      if (nxhm.eq.(nx/2)) then
         at1 = real(sctx(nxhm))
         ex = ex + aimag(fxy(1,1))*at1
         ey = ey + aimag(fxy(2,1))*at1
         ez = ez + aimag(fxy(3,1))*at1
      endif
      gxy(1) = real(ex)
      gxy(2) = real(ey)
      gxy(3) = real(ez)
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
c-----------------------------------------------------------------------
      function vresult(prec)
      implicit none
      real prec, vresult
      vresult = prec
      return
      end
c The following extra procedures are not currently used:
c-----------------------------------------------------------------------
      subroutine MJPOST1GL(part,amu,sctx,qm,nop,idimp,nx,nxhm,nxvh)
c for 1-2/2d code, this subroutine calculates particle momentum flux
c using gridless spectral version, periodic boundaries
c baseline scalar version
c 8*nxhm flops/particle + 3 flops/particle, 2*nxhm loads and stores
c input: all, output: amu
c momentum flux is calculated from the expression:
c amu(i,n) = sum(qci*exp(-sqrt(-1)*2*n*pi*x/nx))
c and qci = qm*vj*vk, where jk = xy,xz for i = 1, 2
c where vj = vj(t-dt/2) and vk = vk(t-dt/2)
c part(1,n) = position x of particle n at t
c part(2,n) = x velocity of particle n at t - dt/2
c part(3,n) = y velocity of particle n at t - dt/2
c part(4,n) = z velocity of particle n at t - dt/2
c amu(i,j) = ith component of momentum flux
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nx = system length in x direction
c nxhm = number of fourier modes kept
c nxvh = second dimension of flux array, must be >= nxhm
      implicit none
      integer nop, idimp, nx, nxhm, nxvh
      real part, qm
      complex amu, sctx
      dimension part(idimp,nop), amu(2,nxvh), sctx(nxvh)
c local data
      integer i, j
      real qmn, dnx, dkx, at1
      real vx, v1, v2
      complex zt1
      dnx = 6.28318530717959/real(nx)
      qmn = qm/real(nx)
c find fourier components
      do 30 i = 1, nop
      vx = qmn*part(2,i)
      v1 = vx*part(3,i)
      v2 = vx*part(4,i)
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c mode numbers 0 < kx < nxhm
      do 20 j = 2, nxhm
      zt1 = sctx(j-1)
      amu(1,j) = amu(1,j) + v1*zt1
      amu(2,j) = amu(2,j) + v2*zt1
   20 continue
c special case to match conventional PIC code
      if (nxhm.eq.(nx/2)) then
         at1 = real(sctx(nxhm))
         amu(1,1) = cmplx(real(amu(1,1))+v1,aimag(amu(1,1))+v1*at1)
         amu(2,1) = cmplx(real(amu(2,1))+v2,aimag(amu(2,1))+v2*at1)
      else
         amu(1,1) = cmplx(real(amu(1,1))+v1,0.0)
         amu(2,1) = cmplx(real(amu(2,1))+v2,0.0)
      endif
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine RMJPOST1GL(part,amu,sctx,qm,ci,nop,idimp,nx,nxhm,nxvh)
c for 1-2/2d code, this subroutine calculates particle momentum flux
c using gridless spectral version for relativistic particles,
c and periodic boundaries
c baseline scalar version
c 8*nxhm flops/particle + 11 flops/particle, 2*nxhm loads and stores
c input: all, output: amu
c momentum flux is calculated from the expression:
c amu(i,n) = sum(qci*exp(-sqrt(-1)*2*n*pi*x/nx)
c and qci = qm*pj*pk*gami2, where jk = xy,xz, for i = 1, 2
c where pj = pj(t-dt/2) and pk = pk(t-dt/2)
c where gami2 = 1./(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n at t
c part(2,n) = x momentum of particle n at t - dt/2
c part(3,n) = y momentum of particle n at t - dt/2
c part(4,n) = z momentum of particle n at t - dt/2
c amu(i,j) = ith component of momentum flux
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 4
c nx = system length in x direction
c nxhm = number of fourier modes kept
c nxvh = first dimension of flux array, must be >= nxhm
      implicit none
      integer nop, idimp, nx, nxhm, nxvh
      real part, qm, ci
      complex amu, sctx
      dimension part(idimp,nop), amu(2,nxvh), sctx(nxvh)
c local data
      integer i, j
      real qmn, ci2, gami2, dnx, dkx, at1
      real vx, vy, vz, p2, v1, v2
      complex zt1
      dnx = 6.28318530717959/real(nx)
      qmn = qm/real(nx)
      ci2 = ci*ci
c find fourier components
      do 30 i = 1, nop
      vx = part(2,i)
      vy = part(3,i)
      vz = part(4,i)
      p2 = vx*vx + vy*vy + vz*vz
      gami2 = qmn/(1.0 + p2*ci2)
      v1 = (vx*vy)*gami2
      v2 = (vx*vz)*gami2
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c mode numbers 0 < kx < nxhm
      do 20 j = 2, nxhm
      zt1 = sctx(j-1)
      amu(1,j) = amu(1,j) + v1*zt1
      amu(2,j) = amu(2,j) + v2*zt1
   20 continue
c special case to match conventional PIC code
      if (nxhm.eq.(nx/2)) then
         at1 = real(sctx(nxhm))
         amu(1,1) = cmplx(real(amu(1,1))+v1,aimag(amu(1,1))+v1*at1)
         amu(2,1) = cmplx(real(amu(2,1))+v2,aimag(amu(2,1))+v2*at1)
      else
         amu(1,1) = cmplx(real(amu(1,1))+v1,0.0)
         amu(2,1) = cmplx(real(amu(2,1))+v2,0.0)
      endif
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DCUPERP13GL(dcu,amu,nx,nxhm,nxvh)
c this subroutine calculates transverse part of the derivative of
c the current density from the momentum flux
c in 1-2/2d with periodic boundary conditions,
c for gridless spectral code
c the transverse part of the derivative of the current is calculated
c using the equations:
c dcu(1,kx) = -sqrt(-1)*kx*vx*vy
c dcu(2,kx) = -sqrt(-1)*kx*vx*vz
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c except for dcu(i,kx=pi) = dcu(i,kx=0) = 0.
c amu(1,j) = xy component of complex momentum flux
c amu(2,j) = xz component of complex momentum flux
c all for fourier mode (j-1)
c nx = system length in x direction
c nxhm = number of fourier modes kept
c nxvh = second dimension of field arrays, must be >= nxhm
      implicit none
      integer nx, nxhm, nxvh
      complex dcu, amu
      dimension dcu(2,nxvh), amu(2,nxvh)
c local data
      integer j
      real dnx, dkx
      complex zero, zt1, zt2
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.0,0.0)
c mode numbers 0 < kx < nxhm
      do 10 j = 2, nxhm
      dkx = dnx*real(j - 1)
      zt2 = cmplx(aimag(amu(1,j)),-real(amu(1,j)))
      dcu(1,j) = dkx*zt2
      zt1 = cmplx(aimag(amu(2,j)),-real(amu(2,j)))
      dcu(2,j) = dkx*zt1
   10 continue
c mode number kx = 0
      dcu(1,1) = zero
      dcu(2,1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine MAXWEL1GL(eyz,byz,cu,ffc,ci,dt,wf,wm,nx,nxhm,nxvh)
c this subroutine solves 1d maxwell's equation in fourier space for
c transverse electric and magnetic fields with periodic boundary
c conditions, for gridless spectral code
c input: all, output: wf, wm, eyz, byz
c approximate flop count is: 87*nxhm
c the magnetic field is first updated half a step using the equations:
c by(kx) = by(kx) + .5*dt*sqrt(-1)*kx*ez(kx)
c bz(kx) = bz(kx) - .5*dt*sqrt(-1)*kx*ey(kx)
c the electric field is then updated a whole step using the equations:
c ey(kx) = ey(kx) - c2*dt*sqrt(-1)*kx*bz(kx) - affp*dt*cuy(kx)*s(kx)
c ez(kx) = ez(kx) + c2*dt*sqrt(-1)*kx*by(kx) - affp*dt*cuz(kx)*s(kx)
c the magnetic field is finally updated the remaining half step with
c the new electric field and the previous magnetic field equations.
c where kx = 2pi*j/nx, c2 = 1./(ci*ci)
c and s(kx) = exp(-((kx*ax)**2/2)
c j = fourier mode numbers
c and similarly for by, bz.
c cu(i,j) = complex current density
c eyz(i,j) = complex transverse electric field
c byz(i,j) = complex magnetic field
c for component i, all for fourier mode (j-1)
c real(ffc(1)) = affp = normalization constant = nx/np,
c where np=number of particles
c aimag(ffc(j)) = finite-size particle shape factor s,
c s(kx) = exp(-((kx*ax)**2)/2)
c for fourier mode (j-1)
c ci = reciprocal of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*sum((1/affp)*|eyz(kx)|**2)
c magnetic field energy is also calculated, using
c wm = nx*sum((c2/affp)*|byz(kx)|**2)
c nx = system length in x direction
c nxhm = number of fourier modes kept
c nxvh = first dimension of field arrays, must be >= nxhm
      implicit none
      integer nx, nxhm, nxvh
      real ci, dt, wf, wm
      complex eyz, byz, cu, ffc
      dimension eyz(2,nxvh), byz(2,nxvh), cu(2,nxvh), ffc(nxhm)
c local data
      integer j
      real dnx, dth, c2, cdt, affp, adt, anorm, dkx, afdt
      complex zero, zt1, zt2, zt5, zt6, zt8, zt9
      double precision wp, ws
      if (ci.le.0.0) return
      dnx = 6.28318530717959/real(nx)
      dth = 0.5*dt
      c2 = 1.0/(ci*ci)
      cdt = c2*dt
      affp = real(ffc(1))
      adt = affp*dt
      zero = cmplx(0.0,0.0)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
c calculate the electromagnetic fields
      do 10 j = 2, nxhm
      dkx = dnx*real(j - 1)
      afdt = adt*aimag(ffc(j))
c update magnetic field half time step
      zt1 = cmplx(-aimag(eyz(2,j)),real(eyz(2,j)))
      zt2 = cmplx(-aimag(eyz(1,j)),real(eyz(1,j)))
      zt5 = byz(1,j) + dth*(dkx*zt1)
      zt6 = byz(2,j) - dth*(dkx*zt2)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt8 = eyz(1,j) - cdt*(dkx*zt1) - afdt*cu(1,j)
      zt9 = eyz(2,j) + cdt*(dkx*zt2) - afdt*cu(2,j)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      eyz(1,j) = zt8
      eyz(2,j) = zt9
      ws = ws + anorm*(zt8*conjg(zt8) + zt9*conjg(zt9))
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2)
      byz(1,j) = zt5
      byz(2,j) = zt6
      wp = wp + anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
   10 continue
c mode number kx = 0
      eyz(1,1) = zero
      eyz(2,1) = zero
      byz(1,1) = zero
      byz(2,1) = zero
      wf = real(nx)*ws
      wm = real(nx)*c2*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine FIELD13GL(fxy,gxy,sctx,nx,nxhm,nxvh)
c for 1-2/2d code, this subroutine calculates real space 1d field gxy at
c integer grid co-ordinates from complex fourier co-efficients fxy
c input: all, output: gxy, sctx
c equations used are:
c gxy(1,x) = sum(fxy(1,n)*exp(sqrt(-1)*2*n*pi*x/nx))
c gxy(2,x) = sum(fxy(2,n)*exp(sqrt(-1)*2*n*pi*x/nx))
c gxy(3,x) = sum(fxy(3,n)*exp(sqrt(-1)*2*n*pi*x/nx))
c fxy(1,n) = x component of force/charge at fourier mode n-1
c fxy(2,n) = y component of force/charge at fourier mode n-1
c fxy(3,n) = z component of force/charge at fourier mode n-1
c if nxhm=nx/2, nx/2 mode is included, in a format compatible with the
c real to complex FFT used in traditional PIC
c gxy(1,j) = x component of force/charge at grid j
c gxy(2,j) = y component of force/charge at grid j
c gxy(3,j) = z component of force/charge at grid j
c sctx = scratch array for sines and cosines
c nx = system length in x direction
c nxhm = number of fourier modes kept in x
c nxvh = second dimension of field arrays, must be >= nxhm
      implicit none
      integer nx, nxhm, nxvh
      real gxy
      complex fxy, sctx
      dimension gxy(3,2*nxvh)
      dimension fxy(3,nxvh), sctx(nxvh)
c local data
      integer j, n
      real dnx, dkx, at1
      double precision ex, ey, ez
      complex zt1
      dnx = 6.28318530717959/real(nx)
      do 30 n = 1, nx
c find field
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*real(n-1)
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
c mode numbers 0 < kx < nxhm
      do 20 j = 2, nxhm
      zt1 = sctx(j-1)
      ex = ex + (real(fxy(1,j)*zt1))
      ey = ey + (real(fxy(2,j)*zt1))
      ez = ez + (real(fxy(3,j)*zt1))
   20 continue
      ex = 2.0d0*ex + real(fxy(1,1))
      ey = 2.0d0*ey + real(fxy(2,1))
      ez = 2.0d0*ez + real(fxy(3,1))
c special case to match conventional PIC code
      if (nxhm.eq.(nx/2)) then
         at1 = real(sctx(nxhm))
         ex = ex + aimag(fxy(1,1))*at1
         ey = ey + aimag(fxy(2,1))*at1
         ez = ez + aimag(fxy(3,1))*at1
      endif
      gxy(1,n) = real(ex)
      gxy(2,n) = real(ey)
      gxy(3,n) = real(ez)
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine ELFIELD1GL(q,fx,ffc,we,nx,nxhm,nxvh)
c this subroutine solves 1d poisson's equation in fourier space for
c unsmoothed longitudinal electric field with periodic boundary
c conditions, for gridless spectral code.
c input: q,ffc,nx,nxhm,nxvh output: fx,we
c equation used is:
c fx(kx) = -sqrt(-1)*kx*g(kx)*q(kx),
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c g(k) = (affp/k**2)*s(k), and s(k) = exp(-(k*ax)**2/2)
c q(j) = complex charge density for fourier mode (j-1)
c fx(j) = x component of complex electric field, for fourier mode (j-1)
c aimag(ffc(j)) = finite-size particle shape factor s
c real(ffc(j) = potential green's function g
c for fourier mode (j-1)
c electric field energy is also calculated, using
c we = nx*sum((affp/k**2)*|q(k)*s(k)|**2)
c where affp = normalization constant = nx/np, where np=number of
c particles
c nx = system length in x direction
c nxhm = number of fourier modes kept
c nxvh = dimension of field array, must be >= nxhm
      implicit none
      integer nx, nxhm, nxvh
      real we
      complex q, fx, ffc
      dimension q(nxvh), fx(nxvh), ffc(nxhm)
c local data
      integer j
      real dnx, at1, at2
      complex zt1, zero
      double precision wp
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.0,0.0)
c calculate unsmoothed longitudinal electric field and sum field energy
      wp = 0.0
c mode numbers 0 < kx < nxhm
      do 10 j = 2, nxhm
      at1 = real(ffc(j))
      at2 = dnx*real(j - 1)*at1
      at1 = at1*aimag(ffc(j))
      zt1 = cmplx(aimag(q(j)),-real(q(j)))
      fx(j) = at2*zt1
      wp = wp + at1*(q(j)*conjg(q(j)))
   10 continue
c mode number kx = 0
      fx(1) = zero
      we = real(nx)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PSMOOTH1GL(q,qs,ffc,nxhm,nxvh)
c this subroutine provides a scalar smoothing function
c with periodic boundary conditions, for gridless spectral code.
c input: q,ffc,nxhm,nxvh, output: qs
c smoothing is calculated using the equation:
c qs(kx) = q(kx)*s(kx)
c where s(kx) = exp(-((kx*ax)**2)/2), and
c where kx = 2pi*j/nx, and j = fourier mode numbers
c q(j) = complex charge density
c aimag(ffc(j)) = finite-size particle shape factor s
c for fourier mode j-1
c ny = system length in y direction
c nxhm = number of fourier modes kept in x
c nxvh = first dimension of field arrays, must be >= nxhm
      implicit none
      integer nxhm, nxvh
      complex q, qs, ffc
      dimension q(nxvh), qs(nxvh)
      dimension ffc(nxhm)
c local data
      integer j
      real at1
c calculate smoothing
c mode numbers 0 <= kx < nxhm
      do 10 j = 1, nxhm
      at1 = aimag(ffc(j))
      qs(j) = at1*q(j)
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine POYNT1GL(q,eyz,byz,ffc,ex0,sx,sy,sz,nx,nxhm,nxvh)
c this subroutine calculates the momentum in the electromagnetic field
c given by the poynting flux, for gridless spectral code.
c inputs are the charge density, transverse electric field, and magnetic
c field.  outputs are sx, sy, sz
c equation used is:
c sx = sum((eyz(1,j))*conjg(byz(2,j))-eyz(2,j)*conjg(byz(1,j)))
c sy = sum(-fx(j)*conjg(byz(2,j)))
c sz = sum(fx(j)*conjg(byz(1,j))), where
c fx(kx) = -sqrt(-1)*kx*g(kx)*q(kx),
c kx = 2pi*j/nx, and j = fourier mode numbers,
c g(kx) = (affp/kx**2)*s(kx), and s(k) = exp(-(k*ax)**2/2)
c q(j) = complex charge density for fourier mode (j-1)
c eyz(1,j) = y component of transverse electric field,
c eyz(2,j) = z component of transverse electric field,
c byz(1,j) = y component of magnetic field,
c byz(2,j) = z component of magnetic field,
c all for fourier mode (j-1)
c real(ffc(j)) = potential green's function g for fourier mode (j-1)
c sx/sy/sz = x/y/z components of field momentum
c nx = system length in x direction
c nxhm = number of fourier modes kept
c nxvh = second dimension of field arrays, must be >= nxhm
      implicit none
      integer nx, nxhm, nxvh
      real ex0, sx, sy, sz
      complex q, eyz, byz, ffc
      dimension q(nxvh), eyz(2,nxvh), byz(2,nxvh)
      dimension ffc(nxhm)
c local data
      integer j
      real dnx, affp, at1
      complex zt1
      double precision wx, wy, wz
      dnx = 6.28318530717959/real(nx)
      affp = real(ffc(1))
c calculate total electromagnetic field momentum
      wx = 0.0d0
      wy = 0.0d0
      wz = 0.0d0
c mode numbers 0 < kx < nxhm
      do 10 j = 2, nxhm
      at1 = dnx*real(j - 1)*real(ffc(j))
      zt1 = cmplx(aimag(q(j)),-real(q(j)))
      zt1 = at1*zt1
      wx = wx + (eyz(1,j)*conjg(byz(2,j))-eyz(2,j)*conjg(byz(1,j)))
      wy = wy - zt1*conjg(byz(2,j))
      wz = wz + zt1*conjg(byz(1,j))
   10 continue
c mode number kx = 0
      wx = 2.0*wx + (eyz(1,1)*conjg(byz(2,1))-eyz(2,1)*conjg(byz(1,1)))
      wy = 2.0*wy - ex0*real(byz(2,1))
      wz = 2.0*wz + ex0*real(byz(1,1))
c normalize
      at1 = real(nx)/affp
      sx = at1*wx
      sy = at1*wy
      sz = at1*wz
      return
      end
c-----------------------------------------------------------------------
      subroutine A0RBPUSH13GL(part,fxyz,byz,sctx,omx,qbm,dt,dtc,ci,ek,  
     1idimp,nop,nx,nxhm,nxvh)
c for 1-2/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time using gridless spectral
c version, for relativistic particles with magnetic field and periodic
c boundaries
c Using the Analytic Boris Mover,
c assumes constant E, B fields, and gamma during a time step
c gamma and energy calculated with Boris half push with electric field
c scalar version using guard cells
c 48*nxhm + 118 flops/particle, 5*nxhm + 4 loads, 4 stores
c 5 divides, 3 sqrts, 1 tangent
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + (q/m)*gx(x(t))) +
c    rot(2)*(py(t-dt/2) + (q/m)*gy(x(t))) +
c    rot(3)*(pz(t-dt/2) + (q/m)*gz(x(t)))) + (q/m)*gx(x(t)))
c py(t+dt/2) = rot(4)*(px(t-dt/2) + (q/m)*gx(x(t))) +
c    rot(5)*(py(t-dt/2) + (q/m)*gy(x(t))) +
c    rot(6)*(pz(t-dt/2) + (q/m)*gz(x(t)))) + (q/m)*gy(x(t)))
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + (q/m)*gx(x(t))) +
c    rot(8)*(py(t-dt/2) + (q/m)*gy(x(t))t) +
c    rot(9)*(pz(t-dt/2) + (q/m)*gz(x(t)))) + (q/m)*gz(x(t)))
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - tan(om*dt/2)**2 + 2*((omx/om)*tan(om*dt/2))**2)/norm
c    rot(2) = 2*((omz/om)*tan(om*dt/2) 
c              + (omx*omy/om**2)*tan(om*dt/2)**2)/norm
c    rot(3) = 2*(-(omy/om)*tan(om*dt/2)
c              + (omx*omz/om**2)*tan(om*dt/2)**2)/norm
c    rot(4) = 2*(-(omz/om)*tan(om*dt/2) 
c              + (omx*omy/om**2)*tan(om*dt/2)**2)/norm
c    rot(5) = (1 - tan(om*dt/2)**2 + 2*((omy/om)*tan(om*dt/2))**2)/norm
c    rot(6) = 2*((omx/om)*tan(om*dt/2) 
c              + (omy*omz/om**2)*tan(om*dt/2)**2)/norm
c    rot(7) = 2*((omy/om)*tan(om*dt/2)
c              + (omx*omz/om**2)*tan(om*dt/2)**2)/norm
c    rot(8) = 2*(-(omx/om)*tan(om*dt/2)
c              + (omy*omz/om**2)*tan(om*dt/2)**2)/norm
c    rot(9) = (1 - tan(om*dt/2)**2 + 2*((omz/om)*tan(om*dt/2))**2)/norm
c norm = 1 + tan(om*dt/2))**2 and om = sqrt(omx**2 + omy**2 + omz**2)
c gx(x) = 0.5*flx(x)*dt + (fx(x)-flx(x))*tan(0.5*om*dt)/om
c gy(x) = 0.5*fly(x)*dt + (fy(x)-fly(x))*tan(0.5*om*dt)/om
c gz(x) = 0.5*flz(x)*dt + (fz(x)-flz(x))*tan(0.5*om*dt)/om
c where flx(x) = fpl(x)*omx/om**2, fly(x) = fpl(x)*omy/om**2,
c flz(x) = fpl(x)*omz/om**2, and fpl(x) = fx(x)*omx+fy(x)*omy+fz(x)*omz
c the rotation matrix is determined by:
c omy = (q/m)*by(x(t))*gami and omz = (q/m)*bz(x(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equation used is:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c fx(x(t)) is calculated from the expression
c fx(x) = sum(fxyz(1,n)*exp(sqrt(-1)*2*n*pi*x/nx))
c similarly for fy(x), fz(x), by(x), bz(x)
c part(1,n) = position x of particle n
c part(2,n) = momentum px of particle n
c part(3,n) = momentum py of particle n
c part(4,n) = momentum pz of particle n
c fxyz(i,n) = i component of force/charge at fourier grid point n
c byz(1,n) = y component of magnetic field at fourier grid point n
c byz(2,n) = z component of magnetic field at fourier grid point n
c sctx = scratch array for sines and cosines
c omx = magnetic field electron cyclotron frequency in x
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t))*dt)**2)/(1. + gami)
c idimp = size of phase space = 4
c nop = number of particles
c nx = system length in x direction
c nxhm = number of fourier modes kept
c nxvh = dimension of field array, must be >= nxhm
      implicit none
      integer idimp, nop, nx, nxhm, nxvh
      real part, omx, qbm, dt, dtc, ci, ek
      complex fxyz, byz, sctx
      dimension part(idimp,nop), fxyz(3,nxvh), byz(2,nxvh)
      dimension sctx(nxvh)
c local data
      integer i, j
      real zero, anx, dth, dnx, dkx, ci2, at1
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt
      real p2, gam, gami, qtmg, dtg, omt, omti, epl
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1, ex, ey, ez, by, bz
      complex zt1
      zero = 0.0
      anx = real(nx)
      dth = 0.5*dt
      dnx = 6.28318530717959/anx
      ci2 = ci*ci
      sum1 = 0.0d0
      do 30 i = 1, nop
c find electric and magnetic field
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
      by = 0.0d0
      bz = 0.0d0
c mode numbers 0 < kx < nxhm
      do 20 j = 2, nxhm
      zt1 = sctx(j-1)
      ex = ex + real(fxyz(1,j)*zt1)
      ey = ey + real(fxyz(2,j)*zt1)
      ez = ez + real(fxyz(3,j)*zt1)
      by = by + real(byz(1,j)*zt1)
      bz = bz + real(byz(2,j)*zt1)
   20 continue
      ex = 2.0d0*ex + real(fxyz(1,1))
      ey = 2.0d0*ey + real(fxyz(2,1))
      ez = 2.0d0*ez + real(fxyz(3,1))
      by = 2.0d0*by + real(byz(1,1))
      bz = 2.0d0*bz + real(byz(2,1))
c special case to match conventional PIC code
      if (nxhm.eq.(nx/2)) then
         at1 = real(sctx(nxhm))
         dx = ex + aimag(fxyz(1,1))*at1
         dy = ey + aimag(fxyz(2,1))*at1
         dz = ez + aimag(fxyz(3,1))*at1
         oy = by + aimag(byz(1,1))*at1
         oz = bz + aimag(byz(2,1))*at1
      else
         dx = ex
         dy = ey
         dz = ez
         oy = by
         oz = bz
      endif
      ox = omx
c normalize electric field
      dx = qbm*dx
      dy = qbm*dy
      dz = qbm*dz
c half acceleration to calculate gamma
      acx = part(2,i)
      acy = part(3,i)
      acz = part(4,i)
      omxt = acx + dx*dth
      omyt = acy + dy*dth
      omzt = acz + dz*dth
c find gamma and inverse gamma
      p2 = omxt*omxt + omyt*omyt + omzt*omzt
      gam = sqrt(1.0 + p2*ci2)
      gami = 1.0/gam
c normalize magnetic field
      ox = qbm*ox
      oy = qbm*oy
      oz = qbm*oz
      qtmg = dth*gami
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c correct the half-acceleration by decomposing E field into components
c parallel and perpendicular to the B field
      omt = sqrt(ox*ox + oy*oy + oz*oz)
      omti = 0.0
      epl = dx*ox + dy*oy + dz*oz
      if (omt.gt.0.0) then
         omti = 1.0/omt
         qtmg = omti*tan(qtmg*omt)
      endif
      epl = epl*omti*omti
      gam = gam*qtmg
      omxt = epl*ox
      omyt = epl*oy
      omzt = epl*oz
      dx = omxt*dth + (dx - omxt)*gam
      dy = omyt*dth + (dy - omyt)*gam
      dz = omzt*dth + (dz - omzt)*gam
      acx = acx + dx
      acy = acy + dy
      acz = acz + dz
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(2,i) = dx
      part(3,i) = dy
      part(4,i) = dz
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,i) + dx*dtg
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
c set new position
      part(1,i) = dx
   30 continue
c normalize kinetic energy
      ek = ek + sum1
      return
      end
