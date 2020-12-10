c Fortran Library for Skeleton 2-1/2D Electromagnetic Gridless PIC Code
c written by Viktor K. Decyk, UCLA
c-----------------------------------------------------------------------
      subroutine DISTR2H(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,idimp,nop,
     1nx,ny,ipbc)
c for 2-1/2d code, this subroutine calculates initial particle
c co-ordinates and velocities with uniform density and maxwellian
c velocity with drift
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c npx/npy = initial number of particles distributed in x/y direction
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer npx, npy, idimp, nop, nx, ny, ipbc
      real vtx, vty, vtz, vdx, vdy, vdz
      real part
      dimension part(idimp,nop)
c local data
      integer j, k, k1, npxy
      real edgelx, edgely, at1, at2, at3, sum1, sum2, sum3
      double precision dsum1, dsum2, dsum3
      double precision ranorm
      npxy = npx*npy
c set boundary values
      edgelx = 0.0
      edgely = 0.0
      at1 = real(nx)/real(npx)
      at2 = real(ny)/real(npy)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         at1 = real(nx-2)/real(npx)
         at2 = real(ny-2)/real(npy)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         at1 = real(nx-2)/real(npx)
      endif
c uniform density profile
      do 20 k = 1, npy
      k1 = npx*(k - 1)
      at3 = edgely + at2*(real(k) - 0.5)
      do 10 j = 1, npx
      part(1,j+k1) = edgelx + at1*(real(j) - 0.5)
      part(2,j+k1) = at3
   10 continue
   20 continue
c maxwellian velocity distribution
      do 30 j = 1, npxy
      part(3,j) = vtx*ranorm()
      part(4,j) = vty*ranorm()
      part(5,j) = vtz*ranorm()
   30 continue
c add correct drift
      dsum1 = 0.0d0
      dsum2 = 0.0d0
      dsum3 = 0.0d0
      do 40 j = 1, npxy
      dsum1 = dsum1 + part(3,j)
      dsum2 = dsum2 + part(4,j)
      dsum3 = dsum3 + part(5,j)
   40 continue
      sum1 = dsum1
      sum2 = dsum2
      sum3 = dsum3
      at1 = 1./real(npxy)
      sum1 = at1*sum1 - vdx
      sum2 = at1*sum2 - vdy
      sum3 = at1*sum3 - vdz
      do 50 j = 1, npxy
      part(3,j) = part(3,j) - sum1
      part(4,j) = part(4,j) - sum2
      part(5,j) = part(5,j) - sum3
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine BPUSH23GL(part,fxy,bxy,sctx,qbm,dt,dtc,ek,idimp,nop,nx,
     1ny,nxhm,nyhm,nxvh,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time using gridless spectral
c version, with magnetic field, periodic boundaries and various particle
c boundary conditions.  Using the Boris Mover.
c 60*nxhm*nyhm + 70 flops/particle, 1 divide, 3*nxhm*nyhm + 5 loads,
c 5 stores
c plus (nxhm + nyhm) sines and cosines/particle
c input: all, output: part, sctx, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
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
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c fx(x(t),y(t)) and fy(x(t),y(t)) are calculated from the expression
c fx(x,y) = sum(fxy(1,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c fxy(i,n,m) = i component of force/charge at fourier grid point n,m
c bxy(i,n,m) = i component of magnetic field at fourier grid point n,m
c if nxhm=nx/2 and nyhm=ny/2, nx/2 and ny/2 modes are included, in a
c format compatible with the real to complex FFT used in traditional PIC
c sctx = scratch array for sines and cosines
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of field arrays, must be >= nxhm
c nyv = second dimension of field arrays, must be >= 2*nyhm
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxhm, nyhm, nxvh, nyv, ipbc
      real part, qbm, dt, dtc, ek
      complex fxy, bxy, sctx
      dimension part(idimp,nop), fxy(3,nxvh,nyv), bxy(3,nxvh,nyv)
      dimension sctx(nxvh)
c local data
      integer i, j, k, k1, nxh, nyh, ny2
      real edgelx, edgely, edgerx, edgery, dnx, dny, dkx, dky, at1, at3
      real qtmh, dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt
      real omt, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real anorm, x, y
      double precision sum1, exl, eyl, ezl, ex, ey, ez
      double precision bxl, byl, bzl, bx, by, bz
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = 2*nyhm + 2
      qtmh = 0.5*qbm*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.0
         edgely = 0.0
         edgerx = real(nx)
         edgery = real(ny)
      else if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 0.0
         edgerx = real(nx-1)
         edgery = real(ny)
      endif
      do 50 i = 1, nop
      x = part(1,i)
      y = part(2,i)
c find electric and magnetic field
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*x
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
      bx = 0.0d0
      by = 0.0d0
      bz = 0.0d0
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k-1)*y
      zt1 = cmplx(cos(dky),sin(dky))
      zt2 = conjg(zt1)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      bxl = 0.0d0
      byl = 0.0d0
      bzl = 0.0d0
      do 20 j = 2, nxhm
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      exl = exl + (real(fxy(1,j,k)*zt3) + real(fxy(1,j,k1)*zt4))
      eyl = eyl + (real(fxy(2,j,k)*zt3) + real(fxy(2,j,k1)*zt4))
      ezl = ezl + (real(fxy(3,j,k)*zt3) + real(fxy(3,j,k1)*zt4))
      bxl = bxl + (real(bxy(1,j,k)*zt3) + real(bxy(1,j,k1)*zt4))
      byl = byl + (real(bxy(2,j,k)*zt3) + real(bxy(2,j,k1)*zt4))
      bzl = bzl + (real(bxy(3,j,k)*zt3) + real(bxy(3,j,k1)*zt4))
   20 continue
c mode number kx = 0
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         zt2 = zt2*real(sctx(nxhm))
         exl = exl + (real(fxy(1,1,k)*zt1) + real(fxy(1,1,k1)*zt2))
         eyl = eyl + (real(fxy(2,1,k)*zt1) + real(fxy(2,1,k1)*zt2))
         ezl = ezl + (real(fxy(3,1,k)*zt1) + real(fxy(3,1,k1)*zt2))
         bxl = bxl + (real(bxy(1,1,k)*zt1) + real(bxy(1,1,k1)*zt2))
         byl = byl + (real(bxy(2,1,k)*zt1) + real(bxy(2,1,k1)*zt2))
         bzl = bzl + (real(bxy(3,1,k)*zt1) + real(bxy(3,1,k1)*zt2))
      else
         exl = exl + real(fxy(1,1,k)*zt1)
         eyl = eyl + real(fxy(2,1,k)*zt1)
         ezl = ezl + real(fxy(3,1,k)*zt1)
         bxl = bxl + real(bxy(1,1,k)*zt1)
         byl = byl + real(bxy(2,1,k)*zt1)
         bzl = bzl + real(bxy(3,1,k)*zt1)
      endif
      ex = ex + exl
      ey = ey + eyl
      ez = ez + ezl
      bx = bx + bxl
      by = by + byl
      bz = bz + bzl
   30 continue
c mode number ky = 0
      k1 = nyhm + 1
      dky = dny*real(nyh)*y
      at1 = cos(dky)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      bxl = 0.0d0
      byl = 0.0d0
      bzl = 0.0d0
      do 40 j = 2, nxhm
      zt3 = sctx(j-1)
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         zt1 = at1*zt3
         exl = exl + (real(fxy(1,j,1)*zt3) + real(fxy(1,j,k1)*zt1))
         eyl = eyl + (real(fxy(2,j,1)*zt3) + real(fxy(2,j,k1)*zt1))
         ezl = ezl + (real(fxy(3,j,1)*zt3) + real(fxy(3,j,k1)*zt1))
         bxl = bxl + (real(bxy(1,j,1)*zt3) + real(bxy(1,j,k1)*zt1))
         byl = byl + (real(bxy(2,j,1)*zt3) + real(bxy(2,j,k1)*zt1))
         bzl = bzl + (real(bxy(3,j,1)*zt3) + real(bxy(3,j,k1)*zt1))
      else
         exl = exl + real(fxy(1,j,1)*zt3)
         eyl = eyl + real(fxy(2,j,1)*zt3)
         ezl = ezl + real(fxy(3,j,1)*zt3)
         bxl = bxl + real(bxy(1,j,1)*zt3)
         byl = byl + real(bxy(2,j,1)*zt3)
         bzl = bzl + real(bxy(3,j,1)*zt3)
      endif
   40 continue
      ex = 2.0d0*(ex + exl)
      ey = 2.0d0*(ey + eyl)
      ez = 2.0d0*(ez + ezl)
      bx = 2.0d0*(bx + bxl)
      by = 2.0d0*(by + byl)
      bz = 2.0d0*(bz + bzl)
      at3 = real(sctx(nxhm))
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         ex = ex + (real(fxy(1,1,1)) + aimag(fxy(1,1,1))*at3)
         ey = ey + (real(fxy(2,1,1)) + aimag(fxy(2,1,1))*at3)
         ez = ez + (real(fxy(3,1,1)) + aimag(fxy(3,1,1))*at3)
         bx = bx + (real(bxy(1,1,1)) + aimag(bxy(1,1,1))*at3)
         by = by + (real(bxy(2,1,1)) + aimag(bxy(2,1,1))*at3)
         bz = bz + (real(bxy(3,1,1)) + aimag(bxy(3,1,1))*at3)
      else
         ex = ex + real(fxy(1,1,1))
         ey = ey + real(fxy(2,1,1))
         ez = ez + real(fxy(3,1,1))
         bx = bx + real(bxy(1,1,1))
         by = by + real(bxy(2,1,1))
         bz = bz + real(bxy(3,1,1))
      endif
c special case to match conventional PIC code
      if (nxhm.eq.nyh) then
         if (nxhm.eq.nxh) then
            ex = ex + (real(fxy(1,1,k1)) + aimag(fxy(1,1,k1))*at3)*at1
            ey = ey + (real(fxy(2,1,k1)) + aimag(fxy(2,1,k1))*at3)*at1
            ez = ez + (real(fxy(3,1,k1)) + aimag(fxy(3,1,k1))*at3)*at1
            bx = bx + (real(bxy(1,1,k1)) + aimag(bxy(1,1,k1))*at3)*at1
            by = by + (real(bxy(2,1,k1)) + aimag(bxy(2,1,k1))*at3)*at1
            bz = bz + (real(bxy(3,1,k1)) + aimag(bxy(3,1,k1))*at3)*at1
         else
            ex = ex + real(fxy(1,1,k1))*at1
            ey = ey + real(fxy(2,1,k1))*at1
            ez = ez + real(fxy(3,1,k1))*at1
            bx = bx + real(bxy(1,1,k1))*at1
            by = by + real(bxy(2,1,k1))*at1
            bz = bz + real(bxy(3,1,k1))*at1
         endif
      endif
      dx = ex
      dy = ey
      dz = ez
      ox = bx
      oy = by
      oz = bz
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,i) + dx
      acy = part(4,i) + dy
      acz = part(5,i) + dz
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
      part(3,i) = dx
      part(4,i) = dy
      part(5,i) = dz
c new position
      dx = x + dx*dtc
      dy = y + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i)
            part(3,i) = -part(3,i)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,i)
            part(4,i) = -part(4,i)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i)
            part(3,i) = -part(3,i)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,i) = dx
      part(2,i) = dy
   50 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine RBPUSH23GL(part,fxy,bxy,sctx,qbm,dt,dtc,ci,ek,idimp,nop
     1,nx,ny,nxhm,nyhm,nxvh,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time using gridless spectral
c version, for relativistic particles with magnetic field, periodic
c boundaries and various particle boundary conditions.
c Using the Boris Mover.
c 60*nxhm*nyhm + 82 flops/particle, 4 divides, 2 sqrts,
c 3*nxhm*nyhm + 5 loads, 5 stores
c plus (nxhm + nyhm) sines and cosines/particle
c input: all, output: part, sctx, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
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
c omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
c omz = (q/m)*bz(x(t),y(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c fx(x(t),y(t)) and fy(x(t),y(t)) are calculated from the expression
c fx(x,y) = sum(fxy(1,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c part(5,n) = momentum pz of particle n
c fxy(i,n,m) = i component of force/charge at fourier grid point n,m
c bxy(i,n,m) = i component of magnetic field at fourier grid point n,m
c if nxhm=nx/2 and nyhm=ny/2, nx/2 and ny/2 modes are included, in a
c format compatible with the real to complex FFT used in traditional PIC
c sctx = scratch array for sines and cosines
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of field arrays, must be >= nxhm
c nyv = second dimension of field arrays, must be >= 2*nyhm
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxhm, nyhm, nxvh, nyv, ipbc
      real part, qbm, dt, dtc, ci, ek
      complex fxy, bxy, sctx
      dimension part(idimp,nop), fxy(3,nxvh,nyv), bxy(3,nxvh,nyv)
      dimension sctx(nxvh)
c local data
      integer i, j, k, k1, nxh, nyh, ny2
      real edgelx, edgely, edgerx, edgery, dnx, dny, dkx, dky, at1, at3
      real qtmh, dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt
      real omt, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real anorm, ci2, p2, gami, qtmg, dtg, x, y
      double precision sum1, exl, eyl, ezl, ex, ey, ez
      double precision bxl, byl, bzl, bx, by, bz
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = 2*nyhm + 2
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.0
         edgely = 0.0
         edgerx = real(nx)
         edgery = real(ny)
      else if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 0.0
         edgerx = real(nx-1)
         edgery = real(ny)
      endif
      do 50 i = 1, nop
      x = part(1,i)
      y = part(2,i)
c find electric and magnetic field
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*x
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
      bx = 0.0d0
      by = 0.0d0
      bz = 0.0d0
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k-1)*y
      zt1 = cmplx(cos(dky),sin(dky))
      zt2 = conjg(zt1)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      bxl = 0.0d0
      byl = 0.0d0
      bzl = 0.0d0
      do 20 j = 2, nxhm
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      exl = exl + (real(fxy(1,j,k)*zt3) + real(fxy(1,j,k1)*zt4))
      eyl = eyl + (real(fxy(2,j,k)*zt3) + real(fxy(2,j,k1)*zt4))
      ezl = ezl + (real(fxy(3,j,k)*zt3) + real(fxy(3,j,k1)*zt4))
      bxl = bxl + (real(bxy(1,j,k)*zt3) + real(bxy(1,j,k1)*zt4))
      byl = byl + (real(bxy(2,j,k)*zt3) + real(bxy(2,j,k1)*zt4))
      bzl = bzl + (real(bxy(3,j,k)*zt3) + real(bxy(3,j,k1)*zt4))
   20 continue
c mode number kx = 0
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         zt2 = zt2*real(sctx(nxhm))
         exl = exl + (real(fxy(1,1,k)*zt1) + real(fxy(1,1,k1)*zt2))
         eyl = eyl + (real(fxy(2,1,k)*zt1) + real(fxy(2,1,k1)*zt2))
         ezl = ezl + (real(fxy(3,1,k)*zt1) + real(fxy(3,1,k1)*zt2))
         bxl = bxl + (real(bxy(1,1,k)*zt1) + real(bxy(1,1,k1)*zt2))
         byl = byl + (real(bxy(2,1,k)*zt1) + real(bxy(2,1,k1)*zt2))
         bzl = bzl + (real(bxy(3,1,k)*zt1) + real(bxy(3,1,k1)*zt2))
      else
         exl = exl + real(fxy(1,1,k)*zt1)
         eyl = eyl + real(fxy(2,1,k)*zt1)
         ezl = ezl + real(fxy(3,1,k)*zt1)
         bxl = bxl + real(bxy(1,1,k)*zt1)
         byl = byl + real(bxy(2,1,k)*zt1)
         bzl = bzl + real(bxy(3,1,k)*zt1)
      endif
      ex = ex + exl
      ey = ey + eyl
      ez = ez + ezl
      bx = bx + bxl
      by = by + byl
      bz = bz + bzl
   30 continue
c mode number ky = 0
      k1 = nyhm + 1
      dky = dny*real(nyh)*y
      at1 = cos(dky)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      bxl = 0.0d0
      byl = 0.0d0
      bzl = 0.0d0
      do 40 j = 2, nxhm
      zt3 = sctx(j-1)
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         zt1 = at1*zt3
         exl = exl + (real(fxy(1,j,1)*zt3) + real(fxy(1,j,k1)*zt1))
         eyl = eyl + (real(fxy(2,j,1)*zt3) + real(fxy(2,j,k1)*zt1))
         ezl = ezl + (real(fxy(3,j,1)*zt3) + real(fxy(3,j,k1)*zt1))
         bxl = bxl + (real(bxy(1,j,1)*zt3) + real(bxy(1,j,k1)*zt1))
         byl = byl + (real(bxy(2,j,1)*zt3) + real(bxy(2,j,k1)*zt1))
         bzl = bzl + (real(bxy(3,j,1)*zt3) + real(bxy(3,j,k1)*zt1))
      else
         exl = exl + real(fxy(1,j,1)*zt3)
         eyl = eyl + real(fxy(2,j,1)*zt3)
         ezl = ezl + real(fxy(3,j,1)*zt3)
         bxl = bxl + real(bxy(1,j,1)*zt3)
         byl = byl + real(bxy(2,j,1)*zt3)
         bzl = bzl + real(bxy(3,j,1)*zt3)
      endif
   40 continue
      ex = 2.0d0*(ex + exl)
      ey = 2.0d0*(ey + eyl)
      ez = 2.0d0*(ez + ezl)
      bx = 2.0d0*(bx + bxl)
      by = 2.0d0*(by + byl)
      bz = 2.0d0*(bz + bzl)
      at3 = real(sctx(nxhm))
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         ex = ex + (real(fxy(1,1,1)) + aimag(fxy(1,1,1))*at3)
         ey = ey + (real(fxy(2,1,1)) + aimag(fxy(2,1,1))*at3)
         ez = ez + (real(fxy(3,1,1)) + aimag(fxy(3,1,1))*at3)
         bx = bx + (real(bxy(1,1,1)) + aimag(bxy(1,1,1))*at3)
         by = by + (real(bxy(2,1,1)) + aimag(bxy(2,1,1))*at3)
         bz = bz + (real(bxy(3,1,1)) + aimag(bxy(3,1,1))*at3)
      else
         ex = ex + real(fxy(1,1,1))
         ey = ey + real(fxy(2,1,1))
         ez = ez + real(fxy(3,1,1))
         bx = bx + real(bxy(1,1,1))
         by = by + real(bxy(2,1,1))
         bz = bz + real(bxy(3,1,1))
      endif
c special case to match conventional PIC code
      if (nxhm.eq.nyh) then
         if (nxhm.eq.nxh) then
            ex = ex + (real(fxy(1,1,k1)) + aimag(fxy(1,1,k1))*at3)*at1
            ey = ey + (real(fxy(2,1,k1)) + aimag(fxy(2,1,k1))*at3)*at1
            ez = ez + (real(fxy(3,1,k1)) + aimag(fxy(3,1,k1))*at3)*at1
            bx = bx + (real(bxy(1,1,k1)) + aimag(bxy(1,1,k1))*at3)*at1
            by = by + (real(bxy(2,1,k1)) + aimag(bxy(2,1,k1))*at3)*at1
            bz = bz + (real(bxy(3,1,k1)) + aimag(bxy(3,1,k1))*at3)*at1
         else
            ex = ex + real(fxy(1,1,k1))*at1
            ey = ey + real(fxy(2,1,k1))*at1
            ez = ez + real(fxy(3,1,k1))*at1
            bx = bx + real(bxy(1,1,k1))*at1
            by = by + real(bxy(2,1,k1))*at1
            bz = bz + real(bxy(3,1,k1))*at1
         endif
      endif
      dx = ex
      dy = ey
      dz = ez
      ox = bx
      oy = by
      oz = bz
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,i) + dx
      acy = part(4,i) + dy
      acz = part(5,i) + dz
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
      part(3,i) = dx
      part(4,i) = dy
      part(5,i) = dz
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = x + dx*dtg
      dy = y + dy*dtg
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i)
            part(3,i) = -part(3,i)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,i)
            part(4,i) = -part(4,i)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i)
            part(3,i) = -part(3,i)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,i) = dx
      part(2,i) = dy
   50 continue
c normalize kinetic energy
      ek = ek + sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine ABPUSH23GL(part,fxy,bxy,sctx,qbm,dt,dtc,ek,idimp,nop,nx
     1,ny,nxhm,nyhm,nxvh,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time using gridless spectral
c version, with magnetic field, periodic boundaries and various particle
c boundary conditions.
c Using the Analytic Boris Mover,
c assumes constant E, B fields during a time step
c 60*nxhm*nyhm + 114 flops/particle, 3*nxhm*nyhm + 5 loads, 5 stores
c 2 divides, 1 sqrt, 1 tangent
c plus (nxhm + nyhm) sines and cosines/particle
c input: all, output: part, sctx, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + (q/m)*gx(x(t),y(t))) +
c    rot(2)*(vy(t-dt/2) + (q/m)*gy(x(t),y(t)))) +
c    rot(3)*(vz(t-dt/2) + (q/m)*gz(x(t),y(t))))) + (q/m)*gx(x(t),y(t))))
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + (q/m)*gx(x(t),y(t)))) +
c    rot(5)*(vy(t-dt/2) + (q/m)*gy(x(t),y(t)))) +
c    rot(6)*(vz(t-dt/2) + (q/m)*gz(x(t),y(t))))) + (q/m)*gy(x(t),y(t))))
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + (q/m)*gx(x(t),y(t)))) +
c    rot(8)*(vy(t-dt/2) + (q/m)*gy(x(t),y(t)))t) +
c    rot(9)*(vz(t-dt/2) + (q/m)*gz(x(t),y(t))))) + (q/m)*gz(x(t),y(t))))
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
c gx(x,y) = 0.5*flx(x,y)*dt + (fx(x,y)-flx(x,y))*tan(0.5*om*dt)/om
c gy(x,y) = 0.5*fly(x,y)*dt + (fy(x,y)-fly(x,y))*tan(0.5*om*dt)/om
c gz(x,y) = 0.5*flz(x,y)*dt + (fz(x,y)-flz(x,y))*tan(0.5*om*dt)/om
c where flx(x,y) = fpl(x,y)*omx/om**2, fly(x,y) = fpl(x,y)*omy/om**2,
c flz(x,y) = fpl(x,y)*omz/om**2,
c and fpl(x,y) = fx(x,y)*omx+fy(x,y)*omy+fz(x,y)*omz
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c fx(x(t),y(t)) and fy(x(t),y(t)) are calculated from the expression
c fx(x,y) = sum(fxy(1,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c fxy(i,n,m) = i component of force/charge at fourier grid point n,m
c bxy(i,n,m) = i component of magnetic field at fourier grid point n,m
c if nxhm=nx/2 and nyhm=ny/2, nx/2 and ny/2 modes are included, in a
c format compatible with the real to complex FFT used in traditional PIC
c sctx = scratch array for sines and cosines
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of field arrays, must be >= nxhm
c nyv = second dimension of field arrays, must be >= 2*nyhm
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxhm, nyhm, nxvh, nyv, ipbc
      real part, qbm, dt, dtc, ek
      complex fxy, bxy, sctx
      dimension part(idimp,nop), fxy(3,nxvh,nyv), bxy(3,nxvh,nyv)
      dimension sctx(nxvh)
c local data
      integer i, j, k, k1, nxh, nyh, ny2
      real edgelx, edgely, edgerx, edgery, dnx, dny, dkx, dky, at1, at3
      real qtmh, dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt
      real omt, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real anorm, dth, omti, epl, x, y
      double precision sum1, exl, eyl, ezl, ex, ey, ez
      double precision bxl, byl, bzl, bx, by, bz
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = 2*nyhm + 2
      dth = 0.5*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.0
         edgely = 0.0
         edgerx = real(nx)
         edgery = real(ny)
      else if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 0.0
         edgerx = real(nx-1)
         edgery = real(ny)
      endif
      do 50 i = 1, nop
      x = part(1,i)
      y = part(2,i)
c find electric and magnetic field
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*x
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
      bx = 0.0d0
      by = 0.0d0
      bz = 0.0d0
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k-1)*y
      zt1 = cmplx(cos(dky),sin(dky))
      zt2 = conjg(zt1)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      bxl = 0.0d0
      byl = 0.0d0
      bzl = 0.0d0
      do 20 j = 2, nxhm
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      exl = exl + (real(fxy(1,j,k)*zt3) + real(fxy(1,j,k1)*zt4))
      eyl = eyl + (real(fxy(2,j,k)*zt3) + real(fxy(2,j,k1)*zt4))
      ezl = ezl + (real(fxy(3,j,k)*zt3) + real(fxy(3,j,k1)*zt4))
      bxl = bxl + (real(bxy(1,j,k)*zt3) + real(bxy(1,j,k1)*zt4))
      byl = byl + (real(bxy(2,j,k)*zt3) + real(bxy(2,j,k1)*zt4))
      bzl = bzl + (real(bxy(3,j,k)*zt3) + real(bxy(3,j,k1)*zt4))
   20 continue
c mode number kx = 0
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         zt2 = zt2*real(sctx(nxhm))
         exl = exl + (real(fxy(1,1,k)*zt1) + real(fxy(1,1,k1)*zt2))
         eyl = eyl + (real(fxy(2,1,k)*zt1) + real(fxy(2,1,k1)*zt2))
         ezl = ezl + (real(fxy(3,1,k)*zt1) + real(fxy(3,1,k1)*zt2))
         bxl = bxl + (real(bxy(1,1,k)*zt1) + real(bxy(1,1,k1)*zt2))
         byl = byl + (real(bxy(2,1,k)*zt1) + real(bxy(2,1,k1)*zt2))
         bzl = bzl + (real(bxy(3,1,k)*zt1) + real(bxy(3,1,k1)*zt2))
      else
         exl = exl + real(fxy(1,1,k)*zt1)
         eyl = eyl + real(fxy(2,1,k)*zt1)
         ezl = ezl + real(fxy(3,1,k)*zt1)
         bxl = bxl + real(bxy(1,1,k)*zt1)
         byl = byl + real(bxy(2,1,k)*zt1)
         bzl = bzl + real(bxy(3,1,k)*zt1)
      endif
      ex = ex + exl
      ey = ey + eyl
      ez = ez + ezl
      bx = bx + bxl
      by = by + byl
      bz = bz + bzl
   30 continue
c mode number ky = 0
      k1 = nyhm + 1
      dky = dny*real(nyh)*y
      at1 = cos(dky)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      bxl = 0.0d0
      byl = 0.0d0
      bzl = 0.0d0
      do 40 j = 2, nxhm
      zt3 = sctx(j-1)
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         zt1 = at1*zt3
         exl = exl + (real(fxy(1,j,1)*zt3) + real(fxy(1,j,k1)*zt1))
         eyl = eyl + (real(fxy(2,j,1)*zt3) + real(fxy(2,j,k1)*zt1))
         ezl = ezl + (real(fxy(3,j,1)*zt3) + real(fxy(3,j,k1)*zt1))
         bxl = bxl + (real(bxy(1,j,1)*zt3) + real(bxy(1,j,k1)*zt1))
         byl = byl + (real(bxy(2,j,1)*zt3) + real(bxy(2,j,k1)*zt1))
         bzl = bzl + (real(bxy(3,j,1)*zt3) + real(bxy(3,j,k1)*zt1))
      else
         exl = exl + real(fxy(1,j,1)*zt3)
         eyl = eyl + real(fxy(2,j,1)*zt3)
         ezl = ezl + real(fxy(3,j,1)*zt3)
         bxl = bxl + real(bxy(1,j,1)*zt3)
         byl = byl + real(bxy(2,j,1)*zt3)
         bzl = bzl + real(bxy(3,j,1)*zt3)
      endif
   40 continue
      ex = 2.0d0*(ex + exl)
      ey = 2.0d0*(ey + eyl)
      ez = 2.0d0*(ez + ezl)
      bx = 2.0d0*(bx + bxl)
      by = 2.0d0*(by + byl)
      bz = 2.0d0*(bz + bzl)
      at3 = real(sctx(nxhm))
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         ex = ex + (real(fxy(1,1,1)) + aimag(fxy(1,1,1))*at3)
         ey = ey + (real(fxy(2,1,1)) + aimag(fxy(2,1,1))*at3)
         ez = ez + (real(fxy(3,1,1)) + aimag(fxy(3,1,1))*at3)
         bx = bx + (real(bxy(1,1,1)) + aimag(bxy(1,1,1))*at3)
         by = by + (real(bxy(2,1,1)) + aimag(bxy(2,1,1))*at3)
         bz = bz + (real(bxy(3,1,1)) + aimag(bxy(3,1,1))*at3)
      else
         ex = ex + real(fxy(1,1,1))
         ey = ey + real(fxy(2,1,1))
         ez = ez + real(fxy(3,1,1))
         bx = bx + real(bxy(1,1,1))
         by = by + real(bxy(2,1,1))
         bz = bz + real(bxy(3,1,1))
      endif
c special case to match conventional PIC code
      if (nxhm.eq.nyh) then
         if (nxhm.eq.nxh) then
            ex = ex + (real(fxy(1,1,k1)) + aimag(fxy(1,1,k1))*at3)*at1
            ey = ey + (real(fxy(2,1,k1)) + aimag(fxy(2,1,k1))*at3)*at1
            ez = ez + (real(fxy(3,1,k1)) + aimag(fxy(3,1,k1))*at3)*at1
            bx = bx + (real(bxy(1,1,k1)) + aimag(bxy(1,1,k1))*at3)*at1
            by = by + (real(bxy(2,1,k1)) + aimag(bxy(2,1,k1))*at3)*at1
            bz = bz + (real(bxy(3,1,k1)) + aimag(bxy(3,1,k1))*at3)*at1
         else
            ex = ex + real(fxy(1,1,k1))*at1
            ey = ey + real(fxy(2,1,k1))*at1
            ez = ez + real(fxy(3,1,k1))*at1
            bx = bx + real(bxy(1,1,k1))*at1
            by = by + real(bxy(2,1,k1))*at1
            bz = bz + real(bxy(3,1,k1))*at1
         endif
      endif
      dx = ex
      dy = ey
      dz = ez
      ox = bx
      oy = by
      oz = bz
c normalize electric field
      dx = qbm*dx
      dy = qbm*dy
      dz = qbm*dz
c half acceleration
      acx = part(3,i)
      acy = part(4,i)
      acz = part(5,i)
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
      part(3,i) = dx
      part(4,i) = dy
      part(5,i) = dz
c new position
      dx = x + dx*dtc
      dy = y + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i)
            part(3,i) = -part(3,i)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,i)
            part(4,i) = -part(4,i)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i)
            part(3,i) = -part(3,i)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,i) = dx
      part(2,i) = dy
   50 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine ARBPUSH23GL(part,fxy,bxy,sctx,qbm,dt,dtc,ci,ek,idimp,  
     1nop,nx,ny,nxhm,nyhm,nxvh,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time using gridless spectral
c version, for relativistic particles with magnetic field, periodic
c boundaries and various particle boundary conditions.
c Using the Analytic Boris Mover,
c assumes constant E, B fields, and gamma during a time step
c scalar version using guard cells
c 60*nxhm*nyhm + 194 flops/particle, 3*nxhm*nyhm + 5 loads, 5 stores
c 6 divides, 3 sqrts, 1 tangent
c plus (nxhm + nyhm) sines and cosines/particle
c input: all, output: part, sctx, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + (q/m)*gx(x(t),y(t))) +
c    rot(2)*(py(t-dt/2) + (q/m)*gy(x(t),y(t))) +
c    rot(3)*(pz(t-dt/2) + (q/m)*gz(x(t),y(t)))) + (q/m)*gx(x(t),y(t)))
c py(t+dt/2) = rot(4)*(px(t-dt/2) + (q/m)*gx(x(t)),y(t)) +
c    rot(5)*(py(t-dt/2) + (q/m)*gy(x(t),y(t))) +
c    rot(6)*(pz(t-dt/2) + (q/m)*gz(x(t),y(t)))) + (q/m)*gy(x(t),y(t)))
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + (q/m)*gx(x(t),y(t))) +
c    rot(8)*(py(t-dt/2) + (q/m)*gy(x(t),y(t))t) +
c    rot(9)*(pz(t-dt/2) + (q/m)*gz(x(t),y(t)))) + (q/m)*gz(x(t),y(t)))
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
c gx(x,y) = 0.5*flx(x,y)*dt + (fx(x,y)-flx(x,y))*tan(0.5*om*dt)/om
c gy(x,y) = 0.5*fly(x,y)*dt + (fy(x,y)-fly(x,y))*tan(0.5*om*dt)/om
c gz(x,y) = 0.5*flz(x,y)*dt + (fz(x,y)-flz(x,y))*tan(0.5*om*dt)/om
c where flx(x,y) = fpl(x,y)*omx/om**2, fly(x,y) = fpl(x,y)*omy/om**2,
c flz(x,y) = fpl(x,y)*omz/om**2,
c and fpl(x,y) = fx(x,y)*omx+fy(x,y)*omy+fz(x,y)*omz
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
c omz = (q/m)*bz(x(t),y(t))*gami, where
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
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c fx(x(t),y(t)) and fy(x(t),y(t)) are calculated from the expression
c fx(x,y) = sum(fxy(1,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c part(5,n) = momentum pz of particle n
c fxy(i,n,m) = i component of force/charge at fourier grid point n,m
c bxy(i,n,m) = i component of magnetic field at fourier grid point n,m
c if nxhm=nx/2 and nyhm=ny/2, nx/2 and ny/2 modes are included, in a
c format compatible with the real to complex FFT used in traditional PIC
c sctx = scratch array for sines and cosines
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = p2(0)/(1.+sqrt(1.+p2(0)*ci*ci))+ug0*th+u2g0*th**2/2+u3g0*th**3/6
c where th = tau/2
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of field arrays, must be >= nxhm
c nyv = second dimension of field arrays, must be >= 2*nyhm
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxhm, nyhm, nxvh, nyv, ipbc
      real part, qbm, dt, dtc, ci, ek
      complex fxy, bxy, sctx
      dimension part(idimp,nop), fxy(3,nxvh,nyv), bxy(3,nxvh,nyv)
      dimension sctx(nxvh)
c local data
      integer i, j, k, k1, nxh, nyh, ny2
      real edgelx, edgely, edgerx, edgery
      real sixth, dnx, dny, dth, ci2, dkx, dky, at1, at3
      real dx, dy, dz, ox, oy, oz, px, py, pz, p2, gam0, wk, ug0, u2g0
      real acx, acy, acz, omt, omti, epl, omxt, omyt, omzt, ugt0, u2gt0
      real gami, tau, f, fp, qtmg, dtg, x, y
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1, exl, eyl, ezl, ex, ey, ez
      double precision bxl, byl, bzl, bx, by, bz
      complex zt1, zt2, zt3, zt4
      sixth = 1.0/6.0
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = 2*nyhm + 2
      dth = 0.5*dt
      ci2 = ci*ci
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.0
         edgely = 0.0
         edgerx = real(nx)
         edgery = real(ny)
      else if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 0.0
         edgerx = real(nx-1)
         edgery = real(ny)
      endif
      do 50 i = 1, nop
      x = part(1,i)
      y = part(2,i)
c find electric and magnetic field
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*x
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
      bx = 0.0d0
      by = 0.0d0
      bz = 0.0d0
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k-1)*y
      zt1 = cmplx(cos(dky),sin(dky))
      zt2 = conjg(zt1)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      bxl = 0.0d0
      byl = 0.0d0
      bzl = 0.0d0
      do 20 j = 2, nxhm
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      exl = exl + (real(fxy(1,j,k)*zt3) + real(fxy(1,j,k1)*zt4))
      eyl = eyl + (real(fxy(2,j,k)*zt3) + real(fxy(2,j,k1)*zt4))
      ezl = ezl + (real(fxy(3,j,k)*zt3) + real(fxy(3,j,k1)*zt4))
      bxl = bxl + (real(bxy(1,j,k)*zt3) + real(bxy(1,j,k1)*zt4))
      byl = byl + (real(bxy(2,j,k)*zt3) + real(bxy(2,j,k1)*zt4))
      bzl = bzl + (real(bxy(3,j,k)*zt3) + real(bxy(3,j,k1)*zt4))
   20 continue
c mode number kx = 0
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         zt2 = zt2*real(sctx(nxhm))
         exl = exl + (real(fxy(1,1,k)*zt1) + real(fxy(1,1,k1)*zt2))
         eyl = eyl + (real(fxy(2,1,k)*zt1) + real(fxy(2,1,k1)*zt2))
         ezl = ezl + (real(fxy(3,1,k)*zt1) + real(fxy(3,1,k1)*zt2))
         bxl = bxl + (real(bxy(1,1,k)*zt1) + real(bxy(1,1,k1)*zt2))
         byl = byl + (real(bxy(2,1,k)*zt1) + real(bxy(2,1,k1)*zt2))
         bzl = bzl + (real(bxy(3,1,k)*zt1) + real(bxy(3,1,k1)*zt2))
      else
         exl = exl + real(fxy(1,1,k)*zt1)
         eyl = eyl + real(fxy(2,1,k)*zt1)
         ezl = ezl + real(fxy(3,1,k)*zt1)
         bxl = bxl + real(bxy(1,1,k)*zt1)
         byl = byl + real(bxy(2,1,k)*zt1)
         bzl = bzl + real(bxy(3,1,k)*zt1)
      endif
      ex = ex + exl
      ey = ey + eyl
      ez = ez + ezl
      bx = bx + bxl
      by = by + byl
      bz = bz + bzl
   30 continue
c mode number ky = 0
      k1 = nyhm + 1
      dky = dny*real(nyh)*y
      at1 = cos(dky)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      bxl = 0.0d0
      byl = 0.0d0
      bzl = 0.0d0
      do 40 j = 2, nxhm
      zt3 = sctx(j-1)
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         zt1 = at1*zt3
         exl = exl + (real(fxy(1,j,1)*zt3) + real(fxy(1,j,k1)*zt1))
         eyl = eyl + (real(fxy(2,j,1)*zt3) + real(fxy(2,j,k1)*zt1))
         ezl = ezl + (real(fxy(3,j,1)*zt3) + real(fxy(3,j,k1)*zt1))
         bxl = bxl + (real(bxy(1,j,1)*zt3) + real(bxy(1,j,k1)*zt1))
         byl = byl + (real(bxy(2,j,1)*zt3) + real(bxy(2,j,k1)*zt1))
         bzl = bzl + (real(bxy(3,j,1)*zt3) + real(bxy(3,j,k1)*zt1))
      else
         exl = exl + real(fxy(1,j,1)*zt3)
         eyl = eyl + real(fxy(2,j,1)*zt3)
         ezl = ezl + real(fxy(3,j,1)*zt3)
         bxl = bxl + real(bxy(1,j,1)*zt3)
         byl = byl + real(bxy(2,j,1)*zt3)
         bzl = bzl + real(bxy(3,j,1)*zt3)
      endif
   40 continue
      ex = 2.0d0*(ex + exl)
      ey = 2.0d0*(ey + eyl)
      ez = 2.0d0*(ez + ezl)
      bx = 2.0d0*(bx + bxl)
      by = 2.0d0*(by + byl)
      bz = 2.0d0*(bz + bzl)
      at3 = real(sctx(nxhm))
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         ex = ex + (real(fxy(1,1,1)) + aimag(fxy(1,1,1))*at3)
         ey = ey + (real(fxy(2,1,1)) + aimag(fxy(2,1,1))*at3)
         ez = ez + (real(fxy(3,1,1)) + aimag(fxy(3,1,1))*at3)
         bx = bx + (real(bxy(1,1,1)) + aimag(bxy(1,1,1))*at3)
         by = by + (real(bxy(2,1,1)) + aimag(bxy(2,1,1))*at3)
         bz = bz + (real(bxy(3,1,1)) + aimag(bxy(3,1,1))*at3)
      else
         ex = ex + real(fxy(1,1,1))
         ey = ey + real(fxy(2,1,1))
         ez = ez + real(fxy(3,1,1))
         bx = bx + real(bxy(1,1,1))
         by = by + real(bxy(2,1,1))
         bz = bz + real(bxy(3,1,1))
      endif
c special case to match conventional PIC code
      if (nxhm.eq.nyh) then
         if (nxhm.eq.nxh) then
            ex = ex + (real(fxy(1,1,k1)) + aimag(fxy(1,1,k1))*at3)*at1
            ey = ey + (real(fxy(2,1,k1)) + aimag(fxy(2,1,k1))*at3)*at1
            ez = ez + (real(fxy(3,1,k1)) + aimag(fxy(3,1,k1))*at3)*at1
            bx = bx + (real(bxy(1,1,k1)) + aimag(bxy(1,1,k1))*at3)*at1
            by = by + (real(bxy(2,1,k1)) + aimag(bxy(2,1,k1))*at3)*at1
            bz = bz + (real(bxy(3,1,k1)) + aimag(bxy(3,1,k1))*at3)*at1
         else
            ex = ex + real(fxy(1,1,k1))*at1
            ey = ey + real(fxy(2,1,k1))*at1
            ez = ez + real(fxy(3,1,k1))*at1
            bx = bx + real(bxy(1,1,k1))*at1
            by = by + real(bxy(2,1,k1))*at1
            bz = bz + real(bxy(3,1,k1))*at1
         endif
      endif
      dx = ex
      dy = ey
      dz = ez
      ox = bx
      oy = by
      oz = bz
c normalize electric field
      dx = qbm*dx
      dy = qbm*dy
      dz = qbm*dz
c normalize magnetic field
      ox = qbm*ox
      oy = qbm*oy
      oz = qbm*oz
c read momentum
      px = part(3,i)
      py = part(4,i)
      pz = part(5,i)
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
      part(3,i) = px
      part(4,i) = py
      part(5,i) = pz
c update inverse gamma
      p2 = px*px + py*py + pz*pz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = x + px*dtg
      dy = y + py*dtg
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i)
            part(3,i) = -part(3,i)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,i)
            part(4,i) = -part(4,i)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i)
            part(3,i) = -part(3,i)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,i) = dx
      part(2,i) = dy
   50 continue
c normalize kinetic energy
      ek = ek + sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine EARBPUSH23GL(part,fxy,bxy,sctx,qbm,dt,dtc,ci,ek,idimp, 
     1nop,nx,ny,nxhm,nyhm,nxvh,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time using gridless spectral
c version, for relativistic particles with magnetic field, periodic
c boundaries and various particle boundary conditions.
c Using the Exact Analytic mover, a variant of the algorithm developed
c by Fei Li et al., where E and B fields are constant and gamma varies
c during a time step.
c scalar version using guard cells
c (48 FLOPs+sin+cos)*nxhm + 197 FLOPs, 12 divides, 6 sqrts per particle
c plus 48 FLOPs, 3 divides, 1 tan, 1 tanh per Newton iteration/particle
c plus (nxhm + nyhm) sines and cosines/particle
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
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c fx(x(t),y(t)) and fy(x(t),y(t)) are calculated from the expression
c fx(x,y) = sum(fxy(1,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c part(5,n) = momentum pz of particle n
c fxy(i,n,m) = i component of force/charge at fourier grid point n,m
c bxy(i,n,m) = i component of magnetic field at fourier grid point n,m
c if nxhm=nx/2 and nyhm=ny/2, nx/2 and ny/2 modes are included, in a
c format compatible with the real to complex FFT used in traditional PIC
c sctx = scratch array for sines and cosines
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using the average
c value of energies at the beginning and end of the time step:
c ek = 0.5*((px(t-dt/2)*2+py(t-dt/2)**2+pz(t-dt/2)**2)/(1.0+gam(0))
c    +      (px(t+dt/2)*2+py(t+dt/2)**2+pz(t+dt/2)**2)/(1.0+gam(tau)))
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of field arrays, must be >= nxhm
c nyv = second dimension of field arrays, must be >= 2*nyhm
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxhm, nyhm, nxvh, nyv, ipbc
      real part, qbm, dt, dtc, ci, ek
      complex fxy, bxy, sctx
      dimension part(idimp,nop), fxy(3,nxvh,nyv), bxy(3,nxvh,nyv)
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
      integer i, j, k, k1, nxh, nyh, ny2, ii
      real edgelx, edgely, edgerx, edgery
      real dnx, dny, dkx, dky, ci2, at1, at3, small, prec, vresult, erm
      real dx, dy, dz, ox, oy, oz, tx, ty, tz, om2, e2, omega, di, ep
      real et, px, py, pz, up, ut, ud, w2, wl2, al2, w, al, wi, ali, wi2
      real ali2, om2t, ws, p2, gam0, dgp0, dgt0, d2gp0, d2gt0, dg0, d2g0
      real wt2, ea, eb, ec, ed, fa, fb, tau, cs, sn, csh, snh, fpp, fp
      real ft, f, t2, t3, wt, wd, tn, tn2, tn2i, csd, snd, cshd, snhd
      real csc, snc, cshc, snhc, fc, fpi, gam, gami, wk, dtg, x, y
      double precision sum1, exl, eyl, ezl, ex, ey, ez
      double precision bxl, byl, bzl, bx, by, bz, dt1, err
      complex zt1, zt2, zt3, zt4
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
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = 2*nyhm + 2
      ci2 = ci*ci
      erm = 0.0
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.0
         edgely = 0.0
         edgerx = real(nx)
         edgery = real(ny)
      else if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 0.0
         edgerx = real(nx-1)
         edgery = real(ny)
      endif
      do 70 i = 1, nop
      x = part(1,i)
      y = part(2,i)
c find electric and magnetic field
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*x
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
      bx = 0.0d0
      by = 0.0d0
      bz = 0.0d0
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k-1)*y
      zt1 = cmplx(cos(dky),sin(dky))
      zt2 = conjg(zt1)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      bxl = 0.0d0
      byl = 0.0d0
      bzl = 0.0d0
      do 20 j = 2, nxhm
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      exl = exl + (real(fxy(1,j,k)*zt3) + real(fxy(1,j,k1)*zt4))
      eyl = eyl + (real(fxy(2,j,k)*zt3) + real(fxy(2,j,k1)*zt4))
      ezl = ezl + (real(fxy(3,j,k)*zt3) + real(fxy(3,j,k1)*zt4))
      bxl = bxl + (real(bxy(1,j,k)*zt3) + real(bxy(1,j,k1)*zt4))
      byl = byl + (real(bxy(2,j,k)*zt3) + real(bxy(2,j,k1)*zt4))
      bzl = bzl + (real(bxy(3,j,k)*zt3) + real(bxy(3,j,k1)*zt4))
   20 continue
c mode number kx = 0
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         zt2 = zt2*real(sctx(nxhm))
         exl = exl + (real(fxy(1,1,k)*zt1) + real(fxy(1,1,k1)*zt2))
         eyl = eyl + (real(fxy(2,1,k)*zt1) + real(fxy(2,1,k1)*zt2))
         ezl = ezl + (real(fxy(3,1,k)*zt1) + real(fxy(3,1,k1)*zt2))
         bxl = bxl + (real(bxy(1,1,k)*zt1) + real(bxy(1,1,k1)*zt2))
         byl = byl + (real(bxy(2,1,k)*zt1) + real(bxy(2,1,k1)*zt2))
         bzl = bzl + (real(bxy(3,1,k)*zt1) + real(bxy(3,1,k1)*zt2))
      else
         exl = exl + real(fxy(1,1,k)*zt1)
         eyl = eyl + real(fxy(2,1,k)*zt1)
         ezl = ezl + real(fxy(3,1,k)*zt1)
         bxl = bxl + real(bxy(1,1,k)*zt1)
         byl = byl + real(bxy(2,1,k)*zt1)
         bzl = bzl + real(bxy(3,1,k)*zt1)
      endif
      ex = ex + exl
      ey = ey + eyl
      ez = ez + ezl
      bx = bx + bxl
      by = by + byl
      bz = bz + bzl
   30 continue
c mode number ky = 0
      k1 = nyhm + 1
      dky = dny*real(nyh)*y
      at1 = cos(dky)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      bxl = 0.0d0
      byl = 0.0d0
      bzl = 0.0d0
      do 40 j = 2, nxhm
      zt3 = sctx(j-1)
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         zt1 = at1*zt3
         exl = exl + (real(fxy(1,j,1)*zt3) + real(fxy(1,j,k1)*zt1))
         eyl = eyl + (real(fxy(2,j,1)*zt3) + real(fxy(2,j,k1)*zt1))
         ezl = ezl + (real(fxy(3,j,1)*zt3) + real(fxy(3,j,k1)*zt1))
         bxl = bxl + (real(bxy(1,j,1)*zt3) + real(bxy(1,j,k1)*zt1))
         byl = byl + (real(bxy(2,j,1)*zt3) + real(bxy(2,j,k1)*zt1))
         bzl = bzl + (real(bxy(3,j,1)*zt3) + real(bxy(3,j,k1)*zt1))
      else
         exl = exl + real(fxy(1,j,1)*zt3)
         eyl = eyl + real(fxy(2,j,1)*zt3)
         ezl = ezl + real(fxy(3,j,1)*zt3)
         bxl = bxl + real(bxy(1,j,1)*zt3)
         byl = byl + real(bxy(2,j,1)*zt3)
         bzl = bzl + real(bxy(3,j,1)*zt3)
      endif
   40 continue
      ex = 2.0d0*(ex + exl)
      ey = 2.0d0*(ey + eyl)
      ez = 2.0d0*(ez + ezl)
      bx = 2.0d0*(bx + bxl)
      by = 2.0d0*(by + byl)
      bz = 2.0d0*(bz + bzl)
      at3 = real(sctx(nxhm))
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         ex = ex + (real(fxy(1,1,1)) + aimag(fxy(1,1,1))*at3)
         ey = ey + (real(fxy(2,1,1)) + aimag(fxy(2,1,1))*at3)
         ez = ez + (real(fxy(3,1,1)) + aimag(fxy(3,1,1))*at3)
         bx = bx + (real(bxy(1,1,1)) + aimag(bxy(1,1,1))*at3)
         by = by + (real(bxy(2,1,1)) + aimag(bxy(2,1,1))*at3)
         bz = bz + (real(bxy(3,1,1)) + aimag(bxy(3,1,1))*at3)
      else
         ex = ex + real(fxy(1,1,1))
         ey = ey + real(fxy(2,1,1))
         ez = ez + real(fxy(3,1,1))
         bx = bx + real(bxy(1,1,1))
         by = by + real(bxy(2,1,1))
         bz = bz + real(bxy(3,1,1))
      endif
c special case to match conventional PIC code
      if (nxhm.eq.nyh) then
         if (nxhm.eq.nxh) then
            ex = ex + (real(fxy(1,1,k1)) + aimag(fxy(1,1,k1))*at3)*at1
            ey = ey + (real(fxy(2,1,k1)) + aimag(fxy(2,1,k1))*at3)*at1
            ez = ez + (real(fxy(3,1,k1)) + aimag(fxy(3,1,k1))*at3)*at1
            bx = bx + (real(bxy(1,1,k1)) + aimag(bxy(1,1,k1))*at3)*at1
            by = by + (real(bxy(2,1,k1)) + aimag(bxy(2,1,k1))*at3)*at1
            bz = bz + (real(bxy(3,1,k1)) + aimag(bxy(3,1,k1))*at3)*at1
         else
            ex = ex + real(fxy(1,1,k1))*at1
            ey = ey + real(fxy(2,1,k1))*at1
            ez = ez + real(fxy(3,1,k1))*at1
            bx = bx + real(bxy(1,1,k1))*at1
            by = by + real(bxy(2,1,k1))*at1
            bz = bz + real(bxy(3,1,k1))*at1
         endif
      endif
      dx = ex
      dy = ey
      dz = ez
      ox = bx
      oy = by
      oz = bz
c normalize magnetic field
      ox = qbm*ox
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
      px = part(3,i)
      py = part(4,i)
      pz = part(5,i)
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
      do 50 ii = 1, imax
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
      if (abs(f) < teps) go to 60
   50 continue
c
c convergence failure
   60 continue
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
      part(3,i) = px
      part(4,i) = py
      part(5,i) = pz
c new position
      dtg = dtc/sqrt(1.0 + p2*ci2)
      dx = x + px*dtg
      dy = y + py*dtg
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i)
            part(3,i) = -part(3,i)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,i)
            part(4,i) = -part(4,i)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i)
            part(3,i) = -part(3,i)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,i) = dx
      part(2,i) = dy
   70 continue
c normalize kinetic energy
      ek = ek + sum1
c sanity check
c     write (*,*) 'gamma sanity check=',erm
      return
      end
c-----------------------------------------------------------------------
      subroutine DPOST2GL(part,q,sctx,qm,nop,idimp,nx,ny,nxhm,nyhm,nxvh,
     1nyv)
c for 2d code, this subroutine calculates particle charge density
c using gridless spectral version, periodic boundaries
c baseline scalar version
c 16*nxhm*nyhm flops/particle
c plus (nxhm + nyhm) sines and cosines/particle
c input: all, output: q, sctx
c charge density is calculated from the expression:
c q(n,m) = sum(qm*exp(-sqrt(-1)*2*n*pi*x/nx)*exp(-sqrt(-1)*2*m*pi*y/ny))
c if nxhm=nx/2 and nyhm=ny/2, nx/2 and ny/2 modes are calculated, in a
c format compatible with the real to complex FFT used in traditional PIC
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c q(n,m) = charge density at fourier grid point n,m
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of field arrays, must be >= nxhm
c nyv = second dimension of field arrays, must be >= 2*nyhm
      implicit none
      integer nop, idimp, nx, ny, nxhm, nyhm, nxvh, nyv
      real part, qm
      complex q, sctx
      dimension part(idimp,nop), q(nxvh,nyv), sctx(nxvh)
c local data
      integer i, j, k, k1, nxh, nyh, ny2
      real qmn, dnx, dny, dkx, dky, at1, at3, x, y
      complex zt1, zt2, zt3
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = 2*nyhm + 2
      qmn = qm/real(nx*ny)
c find fourier components
      do 50 i = 1, nop
      x = part(1,i)
      y = part(2,i)
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*x
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k-1)*y
      zt1 = qmn*cmplx(cos(dky),-sin(dky))
      zt2 = conjg(zt1)
      do 20 j = 2, nxhm
      zt3 = sctx(j-1)
      q(j,k) = q(j,k) + zt1*zt3
      q(j,k1) = q(j,k1) + zt2*zt3
   20 continue
c mode number kx = 0
      q(1,k) = q(1,k) + zt1
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         q(1,k1) = q(1,k1) + zt2*real(sctx(nxhm))
      else
         q(1,k1) = q(1,k1) + zt2
      endif
   30 continue
c mode number ky = 0
      k1 = nyhm + 1
      dky = dny*real(nyh)*y
      at1 = qmn*cos(dky)
      do 40 j = 2, nxhm
      zt3 = sctx(j-1)
      q(j,1) = q(j,1) + qmn*zt3
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         q(j,k1) = q(j,k1) + at1*zt3
      endif
   40 continue
      at3 = real(sctx(nxhm))
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         q(1,1) = cmplx(real(q(1,1))+qmn,aimag(q(1,1))+qmn*at3)
      else
         q(1,1) = cmplx(real(q(1,1))+qmn,0.0)
      endif
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         if (nxhm.eq.nxh) then
            q(1,k1) = cmplx(real(q(1,k1))+at1,aimag(q(1,k1))+at1*at3)
         else
            q(1,k1) = cmplx(real(q(1,k1))+at1,0.0)
         endif
      endif
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DJPOST2GL(part,cu,sctx,qm,dt,nop,idimp,nx,ny,nxhm,nyhm,
     1nxvh,nyv,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using gridless spectral version, periodic boundaries
c in addition, particle positions are advanced a half time-step
c baseline scalar version
c 36*nxhm*nyhm + 7 flops/particle, 6*nxhm*nyhm + 5 loads and stores
c plus (nxhm + nyhm) sines and cosines/particle
c input: all, output: part, cu, sctx
c charge density is calculated from the expression:
c cu(i,n,m) = sum(qmi*exp(-sqrt(-1)*2*n*pi*x/nx)*
c                exp(-sqrt(-1)*2*m*pi*y/ny))
c where qmi = qm*vi, where i = x,y,z
c if nxhm=nx/2 and nyhm=ny/2, nx/2 and ny/2 modes are calculated, in a
c format compatible with the real to complex FFT used in traditional PIC
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c part(5,n) = z velocity of particle n
c cu(i,n,m) = current density at fourier grid point n,m
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of current array, must be >= nxhm
c nyv = second dimension of current array, must be >= 2*nyhm
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nop, idimp, nx, ny, nxhm, nyhm, nxvh, nyv, ipbc
      real part, qm, dt
      complex cu, sctx
      dimension part(idimp,nop), cu(3,nxvh,nyv), sctx(nxvh)
c local data
      integer i, j, k, k1, nxh, nyh, ny2
      real qmn, dnx, dny, dkx, dky, at1, at3
      real x, y, vx, vy, vz, edgelx, edgely, edgerx, edgery, dx, dy
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = 2*nyhm + 2
      qmn = qm/real(nx*ny)
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.0
         edgely = 0.0
         edgerx = real(nx)
         edgery = real(ny)
      else if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 0.0
         edgerx = real(nx-1)
         edgery = real(ny)
      endif
c find fourier components
      do 50 i = 1, nop
      x = part(1,i)
      y = part(2,i)
      vx = qmn*part(3,i)
      vy = qmn*part(4,i)
      vz = qmn*part(5,i)
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*x
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k-1)*y
      zt1 = cmplx(cos(dky),-sin(dky))
      zt2 = conjg(zt1)
      do 20 j = 2, nxhm
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      cu(1,j,k) = cu(1,j,k) + vx*zt3
      cu(2,j,k) = cu(2,j,k) + vy*zt3
      cu(3,j,k) = cu(3,j,k) + vz*zt3
      cu(1,j,k1) = cu(1,j,k1) + vx*zt4
      cu(2,j,k1) = cu(2,j,k1) + vy*zt4
      cu(3,j,k1) = cu(3,j,k1) + vz*zt4
   20 continue
c mode number kx = 0
      cu(1,1,k) = cu(1,1,k) + vx*zt1
      cu(2,1,k) = cu(2,1,k) + vy*zt1
      cu(3,1,k) = cu(3,1,k) + vz*zt1
c special case to match conventional PIC code
      if (nxhm.eq.nxh) zt2 = zt2*real(sctx(nxhm))
      cu(1,1,k1) = cu(1,1,k1) + vx*zt2
      cu(2,1,k1) = cu(2,1,k1) + vy*zt2
      cu(3,1,k1) = cu(3,1,k1) + vz*zt2
   30 continue
c mode number ky = 0
      k1 = nyhm + 1
      dky = dny*real(nyh)*y
      at1 = cos(dky)
      do 40 j = 2, nxhm
      zt3 = sctx(j-1)
      cu(1,j,1) = cu(1,j,1) + vx*zt3
      cu(2,j,1) = cu(2,j,1) + vy*zt3
      cu(3,j,1) = cu(3,j,1) + vz*zt3
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         zt1 = at1*zt3
         cu(1,j,k1) = cu(1,j,k1) + vx*zt1
         cu(2,j,k1) = cu(2,j,k1) + vy*zt1
         cu(3,j,k1) = cu(3,j,k1) + vz*zt1
      endif
   40 continue
      at3 = real(sctx(nxhm))
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         cu(1,1,1) = cmplx(real(cu(1,1,1))+vx,aimag(cu(1,1,1))+vx*at3)
         cu(2,1,1) = cmplx(real(cu(2,1,1))+vy,aimag(cu(2,1,1))+vy*at3)
         cu(3,1,1) = cmplx(real(cu(3,1,1))+vz,aimag(cu(3,1,1))+vz*at3)
      else
         cu(1,1,1) = cmplx(real(cu(1,1,1))+vx,0.0)
         cu(2,1,1) = cmplx(real(cu(2,1,1))+vy,0.0)
         cu(3,1,1) = cmplx(real(cu(3,1,1))+vz,0.0)
      endif
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         if (nxhm.eq.nxh) then
            at3 = at1*at3
            cu(1,1,k1) = cmplx(real(cu(1,1,k1))+vx*at1,aimag(cu(1,1,k1))
     1+vx*at3)
            cu(2,1,k1) = cmplx(real(cu(2,1,k1))+vy*at1,aimag(cu(2,1,k1))
     1+vy*at3)
            cu(3,1,k1) = cmplx(real(cu(3,1,k1))+vz*at1,aimag(cu(3,1,k1))
     1+vz*at3)
         else
            cu(1,1,k1) = cmplx(real(cu(1,1,k1))+vx*at1,0.0)
            cu(2,1,k1) = cmplx(real(cu(2,1,k1))+vy*at1,0.0)
            cu(3,1,k1) = cmplx(real(cu(3,1,k1))+vz*at1,0.0)
         endif
      endif
c advance position half a time-step
      dx = x + part(3,i)*dt
      dy = y + part(4,i)*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i)
            part(3,i) = -part(3,i)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,i)
            part(4,i) = -part(4,i)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i)
            part(3,i) = -part(3,i)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,i) = dx
      part(2,i) = dy
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine RDJPOST2GL(part,cu,sctx,qm,dt,ci,nop,idimp,nx,ny,nxhm, 
     1nyhm,nxvh,nyv,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using gridless spectral version for relativistic particles,
c periodic boundaries
c in addition, particle positions are advanced a half time-step
c baseline scalar version
c 36*nxhm*nyhm + 17 flops/particle, 1 divide, 1, sqrt,
c 6*nxhm*nyhm + 5 loads and stores
c plus (nxhm + nyhm) sines and cosines/particle
c input: all, output: part, cu, sctx
c charge density is calculated from the expression:
c cu(i,n,m) = sum(qci*exp(-sqrt(-1)*2*n*pi*x/nx)*
c                exp(-sqrt(-1)*2*m*pi*y/ny))
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c if nxhm=nx/2 and nyhm=ny/2, nx/2 and ny/2 modes are calculated, in a
c format compatible with the real to complex FFT used in traditional PIC
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x momentum of particle n
c part(4,n) = y momentum of particle n
c part(5,n) = z momentum of particle n
c cu(i,n,m) = current density at fourier grid point n,m
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of current array, must be >= nxhm
c nyv = second dimension of current array, must be >= 2*nyhm
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nop, idimp, nx, ny, nxhm, nyhm, nxvh, nyv, ipbc
      real part, qm, dt, ci
      complex cu, sctx
      dimension part(idimp,nop), cu(3,nxvh,nyv), sctx(nxvh)
c local data
      integer i, j, k, k1, nxh, nyh, ny2
      real qmn, dnx, dny, dkx, dky, at1, at3
      real x, y, vx, vy, vz, edgelx, edgely, edgerx, edgery, dx, dy
      real ci2, p2, gami, px, py, pz
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = 2*nyhm + 2
      qmn = qm/real(nx*ny)
      ci2 = ci*ci
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.0
         edgely = 0.0
         edgerx = real(nx)
         edgery = real(ny)
      else if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 0.0
         edgerx = real(nx-1)
         edgery = real(ny)
      endif
c find fourier components
      do 50 i = 1, nop
      x = part(1,i)
      y = part(2,i)
c find inverse gamma
      vx = part(3,i)
      vy = part(4,i)
      vz = part(5,i)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
      px = vx*gami
      py = vy*gami
      pz = vz*gami
      vx = qmn*px
      vy = qmn*py
      vz = qmn*pz
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*x
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k-1)*y
      zt1 = cmplx(cos(dky),-sin(dky))
      zt2 = conjg(zt1)
      do 20 j = 2, nxhm
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      cu(1,j,k) = cu(1,j,k) + vx*zt3
      cu(2,j,k) = cu(2,j,k) + vy*zt3
      cu(3,j,k) = cu(3,j,k) + vz*zt3
      cu(1,j,k1) = cu(1,j,k1) + vx*zt4
      cu(2,j,k1) = cu(2,j,k1) + vy*zt4
      cu(3,j,k1) = cu(3,j,k1) + vz*zt4
   20 continue
c mode number kx = 0
      cu(1,1,k) = cu(1,1,k) + vx*zt1
      cu(2,1,k) = cu(2,1,k) + vy*zt1
      cu(3,1,k) = cu(3,1,k) + vz*zt1
c special case to match conventional PIC code
      if (nxhm.eq.nxh) zt2 = zt2*real(sctx(nxhm))
      cu(1,1,k1) = cu(1,1,k1) + vx*zt2
      cu(2,1,k1) = cu(2,1,k1) + vy*zt2
      cu(3,1,k1) = cu(3,1,k1) + vz*zt2
   30 continue
c mode number ky = 0
      k1 = nyhm + 1
      dky = dny*real(nyh)*y
      at1 = cos(dky)
      do 40 j = 2, nxhm
      zt3 = sctx(j-1)
      cu(1,j,1) = cu(1,j,1) + vx*zt3
      cu(2,j,1) = cu(2,j,1) + vy*zt3
      cu(3,j,1) = cu(3,j,1) + vz*zt3
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         zt1 = at1*zt3
         cu(1,j,k1) = cu(1,j,k1) + vx*zt1
         cu(2,j,k1) = cu(2,j,k1) + vy*zt1
         cu(3,j,k1) = cu(3,j,k1) + vz*zt1
      endif
   40 continue
      at3 = real(sctx(nxhm))
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         cu(1,1,1) = cmplx(real(cu(1,1,1))+vx,aimag(cu(1,1,1))+vx*at3)
         cu(2,1,1) = cmplx(real(cu(2,1,1))+vy,aimag(cu(2,1,1))+vy*at3)
         cu(3,1,1) = cmplx(real(cu(3,1,1))+vz,aimag(cu(3,1,1))+vz*at3)
      else
         cu(1,1,1) = cmplx(real(cu(1,1,1))+vx,0.0)
         cu(2,1,1) = cmplx(real(cu(2,1,1))+vy,0.0)
         cu(3,1,1) = cmplx(real(cu(3,1,1))+vz,0.0)
      endif
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         if (nxhm.eq.nxh) then
            at3 = at1*at3
            cu(1,1,k1) = cmplx(real(cu(1,1,k1))+vx*at1,aimag(cu(1,1,k1))
     1+vx*at3)
            cu(2,1,k1) = cmplx(real(cu(2,1,k1))+vy*at1,aimag(cu(2,1,k1))
     1+vy*at3)
            cu(3,1,k1) = cmplx(real(cu(3,1,k1))+vz*at1,aimag(cu(3,1,k1))
     1+vz*at3)
         else
            cu(1,1,k1) = cmplx(real(cu(1,1,k1))+vx*at1,0.0)
            cu(2,1,k1) = cmplx(real(cu(2,1,k1))+vy*at1,0.0)
            cu(3,1,k1) = cmplx(real(cu(3,1,k1))+vz*at1,0.0)
         endif
      endif
c advance position half a time-step
      dx = x + px*dt
      dy = y + py*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i)
            part(3,i) = -part(3,i)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,i)
            part(4,i) = -part(4,i)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i)
            part(3,i) = -part(3,i)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,i) = dx
      part(2,i) = dy
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine D2JPOST2GL(part,dcu,sctx,qm,nop,idimp,nx,ny,nxhm,nyhm, 
     1nxvh,nyv)
c for 2-1/2d code, this subroutine calculates the time derivative of the
c particle current density for force-free particles,
c used to initialize the darwin tranverse electric field
c using gridless spectral version, periodic boundaries
c baseline scalar version
c 43*nxhm*nyhm flops per particle
c 12*nxhm*nyhm + 5 loads and stores
c plus (nxhm + nyhm) sines and cosines/particle
c input: all, output: part, dcu, sctx
c derivative of current density is calculated from the expression:
c dcu(i,n,m) = sum(qdi*exp(-sqrt(-1)*kn*x)*exp(-sqrt(-1)*km*y))
c where qdi = qm*vi*d, where i = x,y,z
c kn = 2*n*pi/nx, km = 2*m*pi/ny, k.vi = kn*vx+km*vy, d = -sqrt(-1)*k.vi
c if nxhm=nx/2 and nyhm=ny/2, nx/2 and ny/2 modes are calculated, in a
c format compatible with the real to complex FFT used in traditional PIC
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c part(5,n) = z velocity of particle n
c cu(i,n,m) = charge density at fourier grid point n,m
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of current array, must be >= nxhm
c nyv = second dimension of current array, must be >= 2*nyhm
      implicit none
      integer nop, idimp, nx, ny, nxhm, nyhm, nxvh, nyv
      real part, qm
      complex dcu, sctx
      dimension part(idimp,nop), dcu(3,nxvh,nyv)
      dimension sctx(nxvh)
c local data
      integer i, j, k, k1, nxh, nyh, ny2
      real qmn, dnx, dny, dkx, dky, dkv, at1
      real x, y, ux, uy, vx, vy, vz, dp, dm
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = 2*nyhm + 2
      qmn = qm/real(nx*ny)
c find fourier components
      do 50 i = 1, nop
      x = part(1,i)
      y = part(2,i)
      ux = part(3,i)
      uy = part(4,i)
      vx = qmn*ux
      vy = qmn*uy
      vz = qmn*part(5,i)
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*x
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k-1)
      dkv = dky*uy
      dky = dky*y
      zt1 = cmplx(cos(dky),-sin(dky))
      zt2 = conjg(zt1)
      do 20 j = 2, nxhm
      dkx = dnx*real(j-1)
      dm = dkx*ux
      dp = dm + dkv
      dm = dm - dkv
      zt3 = sctx(j-1)
c deposit current derivative
      zt4 = zt2*zt3*cmplx(0.0,-dm)
      zt3 = zt1*zt3*cmplx(0.0,-dp)
      dcu(1,j,k) = dcu(1,j,k) + vx*zt3
      dcu(2,j,k) = dcu(2,j,k) + vy*zt3
      dcu(3,j,k) = dcu(3,j,k) + vz*zt3
      dcu(1,j,k1) = dcu(1,j,k1) + vx*zt4
      dcu(2,j,k1) = dcu(2,j,k1) + vy*zt4
      dcu(3,j,k1) = dcu(3,j,k1) + vz*zt4
   20 continue
c mode number kx = 0
c deposit current derivative
      dp = dkv
      zt1 = zt1*cmplx(0.0,-dp)
      dcu(1,1,k) = dcu(1,1,k) + vx*zt1
      dcu(2,1,k) = dcu(2,1,k) + vy*zt1
      dcu(3,1,k) = dcu(3,1,k) + vz*zt1
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         dkx = dnx*real(nxh)
         dm = dkx*ux - dkv
         zt2 = zt2*real(sctx(nxhm))*cmplx(0.0,-dm)
         dcu(1,1,k1) = dcu(1,1,k1) + vx*zt2
         dcu(2,1,k1) = dcu(2,1,k1) + vy*zt2
         dcu(3,1,k1) = dcu(3,1,k1) + vz*zt2
      else
         zt2 = zt2*cmplx(0.0,dp)
         dcu(1,1,k1) = dcu(1,1,k1) + vx*zt2
         dcu(2,1,k1) = dcu(2,1,k1) + vy*zt2
         dcu(3,1,k1) = dcu(3,1,k1) + vz*zt2
      endif
   30 continue
c mode number ky = 0
      k1 = nyhm + 1
      dky = dny*real(nyh)
      dkv = dky*uy
      dky = dky*y
      at1 = cos(dky)
      do 40 j = 2, nxhm
      dkx = dnx*real(j-1)
      dp = dkx*ux
      zt3 = sctx(j-1)
c deposit current derivative
      zt3 = zt3*cmplx(0.0,-dp)
      dcu(1,j,1) = dcu(1,j,1) + vx*zt3
      dcu(2,j,1) = dcu(2,j,1) + vy*zt3
      dcu(3,j,1) = dcu(3,j,1) + vz*zt3
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         dp = dkx*ux + dkv
         zt1 = at1*sctx(j-1)*cmplx(0.0,-dp)
         dcu(1,j,k1) = dcu(1,j,k1) + vx*zt1
         dcu(2,j,k1) = dcu(2,j,k1) + vy*zt1
         dcu(3,j,k1) = dcu(3,j,k1) + vz*zt1
      endif
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine RD2JPOST2GL(part,dcu,sctx,qm,ci,nop,idimp,nx,ny,nxhm,  
     1nyhm,nxvh,nyv)
c for 2-1/2d code, this subroutine calculates the time derivative of the
c particle current density for force-free particles,
c used to initialize the darwin tranverse electric field
c using gridless spectral version for relativistic particles,
c periodic boundaries
c baseline scalar version
c 62*nxhm*nyhm flops, 1 divides, 1 sqrt per particle
c 12*nxhm*nyhm + 5 loads and stores
c plus (nxhm + nyhm) sines and cosines/particle
c input: all, output: part, dcu, sctx
c derivative of current density is calculated from the expression:
c dcu(i,n,m) = sum(qdi*exp(-sqrt(-1)*kn*x)*exp(-sqrt(-1)*km*y))
c and qdi = qm*vi*d, where i = x,y,z,
c where vi = pi*gami, gami = 1./sqrt(1.+sum(pi**2)*ci*ci) 
c kn = 2*n*pi/nx, km = 2*m*pi/ny, k.vi = kn*vx+km*vy
c and d = -sqrt(-1)*k.vi
c if nxhm=nx/2 and nyhm=ny/2, nx/2 and ny/2 modes are calculated, in a
c format compatible with the real to complex FFT used in traditional PIC
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x momentum of particle n
c part(4,n) = y momentum of particle n
c part(5,n) = z momentum of particle n
c cu(i,n,m) = charge density at fourier grid point n,m
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of current array, must be >= nxhm
c nyv = second dimension of current array, must be >= 2*nyhm
      implicit none
      integer nop, idimp, nx, ny, nxhm, nyhm, nxvh, nyv
      real part, qm, ci
      complex dcu, sctx
      dimension part(idimp,nop), dcu(3,nxvh,nyv)
      dimension sctx(nxvh)
c local data
      integer i, j, k, k1, nxh, nyh, ny2
      real qmn, dnx, dny, dkx, dky, dkv, at1
      real x, y, ux, uy, vx, vy, vz, ci2, p2, gami, px, py, pz, dp, dm
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = 2*nyhm + 2
      qmn = qm/real(nx*ny)
      ci2 = ci*ci
c find fourier components
      do 50 i = 1, nop
      x = part(1,i)
      y = part(2,i)
      ux = part(3,i)
      uy = part(4,i)
c find inverse gamma
      vx = ux
      vy = uy
      vz = part(5,i)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
      px = vx*gami
      py = vy*gami
      pz = vz*gami
      vx = qmn*px
      vy = qmn*py
      vz = qmn*pz
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*x
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k-1)
      dkv = dky*uy
      dky = dky*y
      zt1 = cmplx(cos(dky),-sin(dky))
      zt2 = conjg(zt1)
      do 20 j = 2, nxhm
      dkx = dnx*real(j-1)
      dm = dkx*ux
      dp = (dm + dkv)*gami
      dm = (dm - dkv)*gami
      zt3 = sctx(j-1)
c deposit current derivative
      zt4 = zt2*zt3*cmplx(0.0,-dm)
      zt3 = zt1*zt3*cmplx(0.0,-dp)
      dcu(1,j,k) = dcu(1,j,k) + vx*zt3
      dcu(2,j,k) = dcu(2,j,k) + vy*zt3
      dcu(3,j,k) = dcu(3,j,k) + vz*zt3
      dcu(1,j,k1) = dcu(1,j,k1) + vx*zt4
      dcu(2,j,k1) = dcu(2,j,k1) + vy*zt4
      dcu(3,j,k1) = dcu(3,j,k1) + vz*zt4
   20 continue
c mode number kx = 0
c deposit current derivative
      dp = dkv*gami
      zt1 = zt1*cmplx(0.0,-dp)
      dcu(1,1,k) = dcu(1,1,k) + vx*zt1
      dcu(2,1,k) = dcu(2,1,k) + vy*zt1
      dcu(3,1,k) = dcu(3,1,k) + vz*zt1
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         dkx = dnx*real(nxh)
         dm = dkx*ux
         dm = (dm - dkv)*gami
         zt2 = zt2*real(sctx(nxhm))*cmplx(0.0,-dm)
         dcu(1,1,k1) = dcu(1,1,k1) + vx*zt2
         dcu(2,1,k1) = dcu(2,1,k1) + vy*zt2
         dcu(3,1,k1) = dcu(3,1,k1) + vz*zt2
      else
         zt2 = zt2*cmplx(0.0,dp)
         dcu(1,1,k1) = dcu(1,1,k1) + vx*zt2
         dcu(2,1,k1) = dcu(2,1,k1) + vy*zt2
         dcu(3,1,k1) = dcu(3,1,k1) + vz*zt2
      endif
   30 continue
c mode number ky = 0
      k1 = nyhm + 1
      dky = dny*real(nyh)
      dkv = dky*uy
      dky = dky*y
      at1 = cos(dky)
      do 40 j = 2, nxhm
      dkx = dnx*real(j-1)
      dp = dkx*ux*gami
      zt3 = sctx(j-1)
c deposit current derivative
      zt3 = zt3*cmplx(0.0,-dp)
      dcu(1,j,1) = dcu(1,j,1) + vx*zt3
      dcu(2,j,1) = dcu(2,j,1) + vy*zt3
      dcu(3,j,1) = dcu(3,j,1) + vz*zt3
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         dp = (dkx*ux + dkv)*gami
         zt1 = at1*sctx(j-1)*cmplx(0.0,-dp)
         dcu(1,j,k1) = dcu(1,j,k1) + vx*zt1
         dcu(2,j,k1) = dcu(2,j,k1) + vy*zt1
         dcu(3,j,k1) = dcu(3,j,k1) + vz*zt1
      endif
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DSJPOST2GL(part,cu,dcu,sctx,qm,ci,nop,idimp,nx,ny,nxhm,
     1nyhm,nxvh,nyv)
c for 2-1/2d code, this subroutine calculates a rescaled particle
c current density and its time derivative for force-free particles,
c used to initialize the electromagnetic fields
c using gridless spectral version for non-relativistic particles,
c periodic boundaries
c baseline scalar version
c 55*nxhm*nyhm flops and 4*nxhm*nyhm divides per particle
c 12*nxhm*nyhm + 5 loads and stores
c plus (nxhm + nyhm) sines and cosines/particle
c input: all, output: part, cu, sctx
c current density is calculated from the expression:
c cu(i,n,m) = sum(qci*exp(-sqrt(-1)*kn*x)*exp(-sqrt(-1)*km*y))
c derivative of current density is calculated from the expression:
c dcu(i,n,m) = sum(qdi*exp(-sqrt(-1)*kn*x)*exp(-sqrt(-1)*km*y))
c where qci = qm*vi*s, qdi = qci*d, where i = x,y,z,
c s = 1.0/(1.0 - (k.vi*ci)**2)/(k*k))
c kn = 2*n*pi/nx, km = 2*m*pi/ny, k*k = kn*kn + km*km,
c k.vi = kn*vx+km*vy and d = -sqrt(-1)*k.vi
c if nxhm=nx/2 and nyhm=ny/2, nx/2 and ny/2 modes are calculated, in a
c format compatible with the real to complex FFT used in traditional PIC
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c part(5,n) = z velocity of particle n
c cu(i,n,m) = charge density at fourier grid point n,m
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of current array, must be >= nxhm
c nyv = second dimension of current array, must be >= 2*nyhm
      implicit none
      integer nop, idimp, nx, ny, nxhm, nyhm, nxvh, nyv
      real part, qm, ci
      complex cu, dcu, sctx
      dimension part(idimp,nop), cu(3,nxvh,nyv), dcu(3,nxvh,nyv)
      dimension sctx(nxvh)
c local data
      integer i, j, k, k1, nxh, nyh, ny2
      real qmn, dnx, dny, dkx, dky, dky2, dk2, dkv, at1, at3
      real x, y, ux, uy, vx, vy, vz, ci2, sp, sm, dp, dm
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = 2*nyhm + 2
      qmn = qm/real(nx*ny)
      ci2 = ci*ci
c find fourier components
      do 50 i = 1, nop
      x = part(1,i)
      y = part(2,i)
      ux = part(3,i)
      uy = part(4,i)
      vx = qmn*ux
      vy = qmn*uy
      vz = qmn*part(5,i)
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*x
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k-1)
      dky2 = dky*dky
      dkv = dky*uy
      dky = dky*y
      zt1 = cmplx(cos(dky),-sin(dky))
      zt2 = conjg(zt1)
      do 20 j = 2, nxhm
      dkx = dnx*real(j-1)
      dk2 = dkx*dkx + dky2
      dm = dkx*ux
      dp = dm + dkv
      dm = dm - dkv
      sp = 1.0/(1.0 - dp*dp*ci2/dk2)
      sm = 1.0/(1.0 - dm*dm*ci2/dk2)
      zt3 = sctx(j-1)
c deposit current
      zt4 = zt2*zt3*sm
      zt3 = zt1*zt3*sp
      cu(1,j,k) = cu(1,j,k) + vx*zt3
      cu(2,j,k) = cu(2,j,k) + vy*zt3
      cu(3,j,k) = cu(3,j,k) + vz*zt3
      cu(1,j,k1) = cu(1,j,k1) + vx*zt4
      cu(2,j,k1) = cu(2,j,k1) + vy*zt4
      cu(3,j,k1) = cu(3,j,k1) + vz*zt4
c deposit current derivative
      zt3 = zt3*cmplx(0.0,-dp)
      zt4 = zt4*cmplx(0.0,-dm)
      dcu(1,j,k) = dcu(1,j,k) + vx*zt3
      dcu(2,j,k) = dcu(2,j,k) + vy*zt3
      dcu(3,j,k) = dcu(3,j,k) + vz*zt3
      dcu(1,j,k1) = dcu(1,j,k1) + vx*zt4
      dcu(2,j,k1) = dcu(2,j,k1) + vy*zt4
      dcu(3,j,k1) = dcu(3,j,k1) + vz*zt4
   20 continue
c mode number kx = 0
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         dkx = dnx*real(nxh)
         dk2 = dkx*dkx + dky2
         dm = dkx*ux
         dp = dm + dkv
         dm = dm - dkv
         sp = 1.0/(1.0 - dp*dp*ci2/dk2)
         sm = 1.0/(1.0 - dm*dm*ci2/dk2)
c deposit current
         zt1 = zt1*sp
         zt2 = zt2*real(sctx(nxhm))*sm
         cu(1,1,k) = cu(1,1,k) + vx*zt1
         cu(2,1,k) = cu(2,1,k) + vy*zt1
         cu(3,1,k) = cu(3,1,k) + vz*zt1
         cu(1,1,k1) = cu(1,1,k1) + vx*zt2
         cu(2,1,k1) = cu(2,1,k1) + vy*zt2
         cu(3,1,k1) = cu(3,1,k1) + vz*zt2
c deposit current derivative
         zt1 = zt1*cmplx(0.0,-dp)
         zt2 = zt2*cmplx(0.0,-dm)
         dcu(1,1,k) = dcu(1,1,k) + vx*zt1
         dcu(2,1,k) = dcu(2,1,k) + vy*zt1
         dcu(3,1,k) = dcu(3,1,k) + vz*zt1
         dcu(1,1,k1) = dcu(1,1,k1) + vx*zt2
         dcu(2,1,k1) = dcu(2,1,k1) + vy*zt2
         dcu(3,1,k1) = dcu(3,1,k1) + vz*zt2
      else
         dp = dkv
         sp = 1.0/(1.0 - dp*dp*ci2/dky2)
c deposit current
         zt1 = zt1*sp
         zt2 = zt2*sp
         cu(1,1,k) = cu(1,1,k) + vx*zt1
         cu(2,1,k) = cu(2,1,k) + vy*zt1
         cu(3,1,k) = cu(3,1,k) + vz*zt1
         cu(1,1,k1) = cu(1,1,k1) + vx*zt2
         cu(2,1,k1) = cu(2,1,k1) + vy*zt2
         cu(3,1,k1) = cu(3,1,k1) + vz*zt2
c deposit current derivative
         zt1 = zt1*cmplx(0.0,-dp)
         zt2 = zt2*cmplx(0.0,dp)
         dcu(1,1,k) = dcu(1,1,k) + vx*zt1
         dcu(2,1,k) = dcu(2,1,k) + vy*zt1
         dcu(3,1,k) = dcu(3,1,k) + vz*zt1
         dcu(1,1,k1) = dcu(1,1,k1) + vx*zt2
         dcu(2,1,k1) = dcu(2,1,k1) + vy*zt2
         dcu(3,1,k1) = dcu(3,1,k1) + vz*zt2
      endif
   30 continue
c mode number ky = 0
      k1 = nyhm + 1
      dky = dny*real(nyh)
      dky2 = dky*dky
      dkv = dky*uy
      dky = dky*y
      at1 = cos(dky)
      do 40 j = 2, nxhm
      dkx = dnx*real(j-1)
      dk2 = dkx*dkx
      dp = dkx*ux
      sp = 1.0/(1.0 - dp*dp*ci2/dk2)
      zt3 = sctx(j-1)
c deposit current
      zt3 = zt3*sp
      cu(1,j,1) = cu(1,j,1) + vx*zt3
      cu(2,j,1) = cu(2,j,1) + vy*zt3
      cu(3,j,1) = cu(3,j,1) + vz*zt3
c deposit current derivative
      zt3 = zt3*cmplx(0.0,-dp)
      dcu(1,j,1) = dcu(1,j,1) + vx*zt3
      dcu(2,j,1) = dcu(2,j,1) + vy*zt3
      dcu(3,j,1) = dcu(3,j,1) + vz*zt3
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         dk2 = dkx*dkx + dky2
         dp = dkx*ux + dkv
         sp = 1.0/(1.0 - dp*dp*ci2/dk2)
c deposit current
         zt1 = at1*sp*sctx(j-1)
         cu(1,j,k1) = cu(1,j,k1) + vx*zt1
         cu(2,j,k1) = cu(2,j,k1) + vy*zt1
         cu(3,j,k1) = cu(3,j,k1) + vz*zt1
c deposit current derivative
         zt1 = zt1*cmplx(0.0,-dp)
         dcu(1,j,k1) = dcu(1,j,k1) + vx*zt1
         dcu(2,j,k1) = dcu(2,j,k1) + vy*zt1
         dcu(3,j,k1) = dcu(3,j,k1) + vz*zt1
      endif
   40 continue
      dkx = dnx*real(nxh)
      dk2 = dkx*dkx
      dp = dkx*ux
      sp = 1.0/(1.0 - dp*dp*ci2/dk2)
      at3 = real(sctx(nxhm))*sp
c special case to match conventional PIC code
c deposit current
      if (nxhm.eq.nxh) then
         cu(1,1,1) = cmplx(real(cu(1,1,1))+vx,aimag(cu(1,1,1))+vx*at3)
         cu(2,1,1) = cmplx(real(cu(2,1,1))+vy,aimag(cu(2,1,1))+vy*at3)
         cu(3,1,1) = cmplx(real(cu(3,1,1))+vz,aimag(cu(3,1,1))+vz*at3)
      else
         cu(1,1,1) = cmplx(real(cu(1,1,1))+vx,0.0)
         cu(2,1,1) = cmplx(real(cu(2,1,1))+vy,0.0)
         cu(3,1,1) = cmplx(real(cu(3,1,1))+vz,0.0)
      endif
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         if (nxhm.eq.nxh) then
            dp = dp + dkv
            sp = 1.0/(1.0 - dp*dp*ci2/dk2)
            at3 = real(sctx(nxhm))
            at1 = at1*sp
            at3 = at1*at3
            cu(1,1,k1) = cmplx(real(cu(1,1,k1))+vx*at1,aimag(cu(1,1,k1))
     1+vx*at3)
            cu(2,1,k1) = cmplx(real(cu(2,1,k1))+vy*at1,aimag(cu(2,1,k1))
     1+vy*at3)
            cu(3,1,k1) = cmplx(real(cu(3,1,k1))+vz*at1,aimag(cu(3,1,k1))
     1+vz*at3)
         else
            cu(1,1,k1) = cmplx(real(cu(1,1,k1))+vx*at1,0.0)
            cu(2,1,k1) = cmplx(real(cu(2,1,k1))+vy*at1,0.0)
            cu(3,1,k1) = cmplx(real(cu(3,1,k1))+vz*at1,0.0)
         endif
      endif
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine RDSJPOST2GL(part,cu,dcu,sctx,qm,ci,nop,idimp,nx,ny,nxhm
     1,nyhm,nxvh,nyv)
c for 2-1/2d code, this subroutine calculates a rescaled particle
c current density and its time derivative for force-free particles,
c used to initialize the electromagnetic fields
c using gridless spectral version for relativistic particles,
c periodic boundaries
c baseline scalar version
c 97*nxhm*nyhm flops and 4*nxhm*nyhm divides per particle
c 12*nxhm*nyhm + 5 loads and stores
c plus (nxhm + nyhm) sines and cosines/particle
c input: all, output: part, cu, sctx
c current density is calculated from the expression:
c cu(i,n,m) = sum(qci*exp(-sqrt(-1)*kn*x)*exp(-sqrt(-1)*km*y))
c derivative of current density is calculated from the expression:
c dcu(i,n,m) = sum(qdi*exp(-sqrt(-1)*kn*x)*exp(-sqrt(-1)*km*y))
c where qci = qm*vi*s, qdi = qci*d, where i = x,y,z,
c gam2 = 1.+sum(pi**2)*ci*ci), gami = 1./sqrt(gam2), vi = pi*gami, 
c s = gam2/(1.0 + ((sum(pi**2) - (k.pi)**2/(k*k))*ci*ci
c kn = 2*n*pi/nx, km = 2*m*pi/ny, k*k = kn*kn + km*km,
c k.pi = kn*px+km*py and d = -sqrt(-1)*k.vi
c if nxhm=nx/2 and nyhm=ny/2, nx/2 and ny/2 modes are calculated, in a
c format compatible with the real to complex FFT used in traditional PIC
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x momentum of particle n
c part(4,n) = y momentum of particle n
c part(5,n) = z momentum of particle n
c cu(i,n,m) = charge density at fourier grid point n,m
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of current array, must be >= nxhm
c nyv = second dimension of current array, must be >= 2*nyhm
      implicit none
      integer nop, idimp, nx, ny, nxhm, nyhm, nxvh, nyv
      real part, qm, ci
      complex cu, dcu, sctx
      dimension part(idimp,nop), cu(3,nxvh,nyv), dcu(3,nxvh,nyv)
      dimension sctx(nxvh)
c local data
      integer i, j, k, k1, nxh, nyh, ny2
      real qmn, dnx, dny, dkx, dky, dky2, dk2, dkv, at1, at3
      real x, y, ux, uy, vx, vy, vz, ci2, p2, gam2, gami, px, py, pz
      real sp, sm, dp, dm
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = 2*nyhm + 2
      qmn = qm/real(nx*ny)
      ci2 = ci*ci
c find fourier components
      do 50 i = 1, nop
      x = part(1,i)
      y = part(2,i)
      ux = part(3,i)
      uy = part(4,i)
c find inverse gamma
      vx = ux
      vy = uy
      vz = part(5,i)
      p2 = vx*vx + vy*vy + vz*vz
      gam2 = 1.0 + p2*ci2
      gami = 1.0/sqrt(gam2)
      px = vx*gami
      py = vy*gami
      pz = vz*gami
      vx = qmn*px
      vy = qmn*py
      vz = qmn*pz
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*x
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k-1)
      dky2 = dky*dky
      dkv = dky*uy
      dky = dky*y
      zt1 = cmplx(cos(dky),-sin(dky))
      zt2 = conjg(zt1)
      do 20 j = 2, nxhm
      dkx = dnx*real(j-1)
      dk2 = dkx*dkx + dky2
      dm = dkx*ux
      dp = dm + dkv
      dm = dm - dkv
      sp = dp*dp/dk2
      sm = dm*dm/dk2
      dp = dp*gami
      dm = dm*gami
      sp = gam2/(1.0 + (p2 - sp)*ci2)
      sm = gam2/(1.0 + (p2 - sm)*ci2)
      zt3 = sctx(j-1)
c deposit current
      zt4 = zt2*zt3*sm
      zt3 = zt1*zt3*sp
      cu(1,j,k) = cu(1,j,k) + vx*zt3
      cu(2,j,k) = cu(2,j,k) + vy*zt3
      cu(3,j,k) = cu(3,j,k) + vz*zt3
      cu(1,j,k1) = cu(1,j,k1) + vx*zt4
      cu(2,j,k1) = cu(2,j,k1) + vy*zt4
      cu(3,j,k1) = cu(3,j,k1) + vz*zt4
c deposit current derivative
      zt3 = zt3*cmplx(0.0,-dp)
      zt4 = zt4*cmplx(0.0,-dm)
      dcu(1,j,k) = dcu(1,j,k) + vx*zt3
      dcu(2,j,k) = dcu(2,j,k) + vy*zt3
      dcu(3,j,k) = dcu(3,j,k) + vz*zt3
      dcu(1,j,k1) = dcu(1,j,k1) + vx*zt4
      dcu(2,j,k1) = dcu(2,j,k1) + vy*zt4
      dcu(3,j,k1) = dcu(3,j,k1) + vz*zt4
   20 continue
c mode number kx = 0
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         dkx = dnx*real(nxh)
         dk2 = dkx*dkx + dky2
         dm = dkx*ux
         dp = dm + dkv
         dm = dm - dkv
         sp = dp*dp/dk2
         sm = dm*dm/dk2
         dp = dp*gami
         dm = dm*gami
         sp = gam2/(1.0 + (p2 - sp)*ci2)
         sm = gam2/(1.0 + (p2 - sm)*ci2)
c deposit current
         zt1 = zt1*sp
         zt2 = zt2*real(sctx(nxhm))*sm
         cu(1,1,k) = cu(1,1,k) + vx*zt1
         cu(2,1,k) = cu(2,1,k) + vy*zt1
         cu(3,1,k) = cu(3,1,k) + vz*zt1
         cu(1,1,k1) = cu(1,1,k1) + vx*zt2
         cu(2,1,k1) = cu(2,1,k1) + vy*zt2
         cu(3,1,k1) = cu(3,1,k1) + vz*zt2
c deposit current derivative
         zt1 = zt1*cmplx(0.0,-dp)
         zt2 = zt2*cmplx(0.0,-dm)
         dcu(1,1,k) = dcu(1,1,k) + vx*zt1
         dcu(2,1,k) = dcu(2,1,k) + vy*zt1
         dcu(3,1,k) = dcu(3,1,k) + vz*zt1
         dcu(1,1,k1) = dcu(1,1,k1) + vx*zt2
         dcu(2,1,k1) = dcu(2,1,k1) + vy*zt2
         dcu(3,1,k1) = dcu(3,1,k1) + vz*zt2
      else
         dp = dkv
         sp = dp*dp/dky2
         dp = dp*gami
         sp = gam2/(1.0 + (p2 - sp)*ci2)
c deposit current
         zt1 = zt1*sp
         zt2 = zt2*sp
         cu(1,1,k) = cu(1,1,k) + vx*zt1
         cu(2,1,k) = cu(2,1,k) + vy*zt1
         cu(3,1,k) = cu(3,1,k) + vz*zt1
         cu(1,1,k1) = cu(1,1,k1) + vx*zt2
         cu(2,1,k1) = cu(2,1,k1) + vy*zt2
         cu(3,1,k1) = cu(3,1,k1) + vz*zt2
c deposit current derivative
         zt1 = zt1*cmplx(0.0,-dp)
         zt2 = zt2*cmplx(0.0,dp)
         dcu(1,1,k) = dcu(1,1,k) + vx*zt1
         dcu(2,1,k) = dcu(2,1,k) + vy*zt1
         dcu(3,1,k) = dcu(3,1,k) + vz*zt1
         dcu(1,1,k1) = dcu(1,1,k1) + vx*zt2
         dcu(2,1,k1) = dcu(2,1,k1) + vy*zt2
         dcu(3,1,k1) = dcu(3,1,k1) + vz*zt2
      endif
   30 continue
c mode number ky = 0
      k1 = nyhm + 1
      dky = dny*real(nyh)
      dky2 = dky*dky
      dkv = dky*uy
      dky = dky*y
      at1 = cos(dky)
      do 40 j = 2, nxhm
      dkx = dnx*real(j-1)
      dk2 = dkx*dkx
      dp = dkx*ux
      sp = dp*dp/dk2
      dp = dp*gami
      sp = gam2/(1.0 + (p2 - sp)*ci2)
      zt3 = sctx(j-1)
c deposit current
      zt3 = zt3*sp
      cu(1,j,1) = cu(1,j,1) + vx*zt3
      cu(2,j,1) = cu(2,j,1) + vy*zt3
      cu(3,j,1) = cu(3,j,1) + vz*zt3
c deposit current derivative
      zt3 = zt3*cmplx(0.0,-dp)
      dcu(1,j,1) = dcu(1,j,1) + vx*zt3
      dcu(2,j,1) = dcu(2,j,1) + vy*zt3
      dcu(3,j,1) = dcu(3,j,1) + vz*zt3
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         dk2 = dkx*dkx + dky2
         dp = dkx*ux + dkv
         sp = dp*dp/dk2
         dp = dp*gami
         sp = gam2/(1.0 + (p2 - sp)*ci2)
c deposit current
         zt1 = at1*sp*sctx(j-1)
         cu(1,j,k1) = cu(1,j,k1) + vx*zt1
         cu(2,j,k1) = cu(2,j,k1) + vy*zt1
         cu(3,j,k1) = cu(3,j,k1) + vz*zt1
c deposit current derivative
         zt1 = zt1*cmplx(0.0,-dp)
         dcu(1,j,k1) = dcu(1,j,k1) + vx*zt1
         dcu(2,j,k1) = dcu(2,j,k1) + vy*zt1
         dcu(3,j,k1) = dcu(3,j,k1) + vz*zt1
      endif
   40 continue
      dkx = dnx*real(nxh)
      dk2 = dkx*dkx
      dp = dkx*ux
      sp = dp*dp/dk2
      dp = dp*gami
      sp = gam2/(1.0 + (p2 - sp)*ci2)
      at3 = real(sctx(nxhm))*sp
c special case to match conventional PIC code
c deposit current
      if (nxhm.eq.nxh) then
         cu(1,1,1) = cmplx(real(cu(1,1,1))+vx,aimag(cu(1,1,1))+vx*at3)
         cu(2,1,1) = cmplx(real(cu(2,1,1))+vy,aimag(cu(2,1,1))+vy*at3)
         cu(3,1,1) = cmplx(real(cu(3,1,1))+vz,aimag(cu(3,1,1))+vz*at3)
      else
         cu(1,1,1) = cmplx(real(cu(1,1,1))+vx,0.0)
         cu(2,1,1) = cmplx(real(cu(2,1,1))+vy,0.0)
         cu(3,1,1) = cmplx(real(cu(3,1,1))+vz,0.0)
      endif
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         if (nxhm.eq.nxh) then
            dp = dp + dkv
            sp = dp*dp/dk2
            dp = dp*gami
            sp = gam2/(1.0 + (p2 - sp)*ci2)
            at3 = real(sctx(nxhm))
            at1 = at1*sp
            at3 = at1*at3
            cu(1,1,k1) = cmplx(real(cu(1,1,k1))+vx*at1,aimag(cu(1,1,k1))
     1+vx*at3)
            cu(2,1,k1) = cmplx(real(cu(2,1,k1))+vy*at1,aimag(cu(2,1,k1))
     1+vy*at3)
            cu(3,1,k1) = cmplx(real(cu(3,1,k1))+vz*at1,aimag(cu(3,1,k1))
     1+vz*at3)
         else
            cu(1,1,k1) = cmplx(real(cu(1,1,k1))+vx*at1,0.0)
            cu(2,1,k1) = cmplx(real(cu(2,1,k1))+vy*at1,0.0)
            cu(3,1,k1) = cmplx(real(cu(3,1,k1))+vz*at1,0.0)
         endif
      endif
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FFC2INITGL(ffc,ax,ay,affp,nx,ny,nxhm,nyhm)
c calculates tables needed by a two dimensional gridless field solvers
c Truncated fourier series gives a sinc function instead of a delta 
c function for particle.  The smoothing function, the convolution of a
c box car function with a gaussian, integrates the sinc shape to obtain
c a smoother shape with minimized oscillations.
c input: ax,ay,affp,nx,ny,nxhm,nyhm, output: ffc
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = sin(kx*wx)/kx*wx)*sin(ky*wy)/ky*wy)
c            *exp(-((kx*ax)**2+(ky*ay)**2)/2)
c where kx = 2pi*j/nx and ky = 2pi*k/ny, j-1,k-1 = fourier mode numbers,
c wx = 1.9269/(dnx*real(nxhm)), wy = 1.9269/(dny*real(nyhm)), and
c Si(1.9269) = pi/2
c aimag(ffc(j,k)) = finite-size particle shape factor s
c real(ffc(j,k)) = potential green's function g
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
      implicit none
      integer nx, ny, nxhm, nyhm
      real ax, ay, affp
      complex ffc
      dimension ffc(nxhm,nyhm)
c local data
      integer j, k
      real dnx, dny, dkx, dky, wx, wy, at1, at2, at3, at4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      wx = 1.9269/(dnx*real(nxhm))
      wy = 1.9269/(dny*real(nyhm))
c prepare form factor array
      do 20 k = 1, nyhm
      dky = dny*real(k - 1)
      at1 = dky*dky
      if (k.eq.1) then
         at2 = 1.0
      else
         at2 = dky*wy
         at2 = (sin(at2)/at2)*exp(-.5*(dky*ay)**2)
      endif
      do 10 j = 1, nxhm
      dkx = dnx*real(j - 1)
      at3 = dkx*dkx + at1
      if (j.eq.1) then
         at4 = 1.0
      else
         at4 = dkx*wx
         at4 = (sin(at4)/at4)*exp(-.5*(dkx*ax)**2)
      endif
      at4 = at2*at4
      if (at3.eq.0.0) then
         ffc(j,k) = cmplx(affp,1.0)
      else
         ffc(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine POIS23GL(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxhm,nyhm,
     1nxvh,nyv)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions, for gridless spectral code.
c Zeros out z component
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxhm,nyhm,nxvh,nyv
c output: ffc
c for isign /= 0, input: q,ffc,isign,nx,ny,nxhm,nyhm,nxvh,nyv
c output: fxy,we
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2),
c q(j,k) = complex charge density for fourier mode (j-1,k-1)
c fxy(1,j,k) = x component of complex force/charge,
c fxy(2,j,k) = y component of complex force/charge,
c fxy(3,j,k) = zero,
c all for fourier mode (j-1,k-1)
c zero row at k = nyhm + 1 is for compatibility with FFT solvers
c when nxhm = nx/2 and nyhm = ny/2
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated
c aimag(ffc(j,k)) = finite-size particle shape factor s
c real(ffc(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of field arrays, must be >= nxhm
c nyv = second dimension of field arrays, must be >= 2*nyhm
      implicit none
      integer isign, nx, ny, nxhm, nyhm, nxvh, nyv
      real ax, ay, affp, we
      complex q, fxy, ffc
      dimension q(nxvh,nyv), fxy(3,nxvh,nyv)
      dimension ffc(nxhm,nyhm)
c local data
      integer j, k, k1, ny2
      real dnx, dny, dkx, dky, at1, at2, at3, at4
      complex zt1, zt2, zero
      double precision wp, swp
      ny2 = 2*nyhm + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyhm
      dky = dny*real(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nxhm
      dkx = dnx*real(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.0) then
         ffc(j,k) = cmplx(affp,1.0)
      else
         ffc(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
c calculate force/charge and sum field energy
   30 swp = 0.0d0
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 50 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k - 1)
      wp = 0.0
      do 40 j = 2, nxhm
      at1 = real(ffc(j,k))*aimag(ffc(j,k))
      at2 = dnx*real(j - 1)*at1
      at3 = dky*at1
      zt1 = cmplx(aimag(q(j,k)),-real(q(j,k)))
      zt2 = cmplx(aimag(q(j,k1)),-real(q(j,k1)))
c ky > 0
      fxy(1,j,k) = at2*zt1
      fxy(2,j,k) = at3*zt1
      fxy(3,j,k) = zero
c ky < 0
      fxy(1,j,k1) = at2*zt2
      fxy(2,j,k1) = -at3*zt2
      fxy(3,j,k1) = zero
      wp = wp + at1*(q(j,k)*conjg(q(j,k)) + q(j,k1)*conjg(q(j,k1)))
   40 continue
      swp = swp + wp
   50 continue
c mode number kx = 0
      wp = 0.0
      do 60 k = 2, nyhm
      k1 = ny2 - k
      at1 = real(ffc(1,k))*aimag(ffc(1,k))
      at3 = dny*real(k - 1)*at1
c ky > 0
      zt1 = cmplx(aimag(q(1,k)),-real(q(1,k)))
      fxy(1,1,k) = zero
      fxy(2,1,k) = at3*zt1
      fxy(3,1,k) = zero
c ky < 0
      fxy(1,1,k1) = zero
      fxy(3,1,k1) = zero
      if (nyhm.eq.(ny/2)) then
         fxy(2,1,k1) = zero
      else
         zt1 = cmplx(aimag(q(1,k1)),-real(q(1,k1)))
         fxy(2,1,k1) = -at3*zt1
      endif
      wp = wp + at1*(q(1,k)*conjg(q(1,k)))
   60 continue
      swp = swp + wp
c mode number ky = 0
      k1 = nyhm + 1
      wp = 0.0
      do 70 j = 2, nxhm
      at1 = real(ffc(j,1))*aimag(ffc(j,1))
      at2 = dnx*real(j - 1)*at1
      zt1 = cmplx(aimag(q(j,1)),-real(q(j,1)))
      fxy(1,j,1) = at2*zt1
      fxy(2,j,1) = zero
      fxy(3,j,1) = zero
c unused extra row
      fxy(1,j,k1) = zero
      fxy(2,j,k1) = zero
      fxy(3,j,k1) = zero
      wp = wp + at1*(q(j,1)*conjg(q(j,1)))
   70 continue
      swp = swp + wp
c mode number kx = 0, ky = 0
      fxy(1,1,1) = zero
      fxy(2,1,1) = zero
      fxy(3,1,1) = zero
c unused extra element
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      fxy(3,1,k1) = zero
      we = real(nx*ny)*swp
      return
      end
c-----------------------------------------------------------------------
      subroutine CUPERP2GL(cu,nx,ny,nxhm,nyhm,nxvh,nyv)
c this subroutine calculates the transverse current in fourier space
c for gridless spectral code.
c input: all, output: cu
c approximate flop count is: 36*nxhm*nyhm
c and nxhm*nyhm divides
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for cux(kx=pi) = cuy(kx=pi) = 0, cux(ky=pi) = cuy(ky=pi) = 0,
c and cux(kx=0,ky=0) = cuy(kx=0,ky=0) = 0.
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c zero row at k = nyhm + 1 is for compatibility with FFT solvers
c when nxhm = nx/2 and nyhm = ny/2
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of current array, must be >= nxhm
c nyv = second dimension of current array, must be >= 2*nyhm
      implicit none
      integer nx, ny, nxhm, nyhm, nxvh, nyv
      complex cu
      dimension cu(3,nxvh,nyv)
c local data
      integer ny2, j, k, k1
      real dnx, dny, dkx, dky, dky2, at1
      complex zero, zt1
      ny2 = 2*nyhm + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
c calculate transverse part of current
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 20 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nxhm
      dkx = dnx*real(j - 1)
      at1 = 1.0/(dkx*dkx + dky2)
      zt1 = at1*(dkx*cu(1,j,k) + dky*cu(2,j,k))
c ky > 0
      cu(1,j,k) = cu(1,j,k) - dkx*zt1
      cu(2,j,k) = cu(2,j,k) - dky*zt1
c ky < 0
      zt1 = at1*(dkx*cu(1,j,k1) - dky*cu(2,j,k1))
      cu(1,j,k1) = cu(1,j,k1) - dkx*zt1
      cu(2,j,k1) = cu(2,j,k1) + dky*zt1
   10 continue
   20 continue
c mode number kx = 0
      do 30 k = 2, nyhm
      k1 = ny2 - k
c ky > 0
      cu(2,1,k) = zero
c ky < 0
      if (nyhm.eq.(ny/2)) cu(1,1,k1) = zero
      cu(2,1,k1) = zero
   30 continue
c mode number ky = 0
      k1 = nyhm + 1
      do 40 j = 2, nxhm
      cu(1,j,1) = zero
c unused extra row
      cu(1,j,k1) = zero
      cu(2,j,k1) = zero
   40 continue
c unused extra element
      cu(1,1,k1) = zero
      cu(2,1,k1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine IBPOIS23GL(cu,bxy,ffc,ci,wm,nx,ny,nxhm,nyhm,nxvh,nyv)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field, with periodic boundary conditions,
c for gridless spectral code.
c input: cu,ffc,ci,nx,ny,nxhm,nyhm,nxvh, output: bxy,wm
c approximate flop count is: 90*nxhm*nyhm + 40*(nxhm + nyhm)
c the magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky)),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bx(kx=pi) = by(kx=pi) = bz(kx=pi) = bx(ky=pi) = by(ky=pi) = bz(ky=pi) 
c = 0, and bx(kx=0,ky=0) = by(kx=0,ky=0) = bz(kx=0,ky=0) = 0.
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c bxy(i,j,k) = i component of complex magnetic field
c zero row at k = nyhm + 1 is for compatibility with FFT solvers
c when nxhm = nx/2 and nyhm = ny/2
c all for fourier mode (j-1,k-1)
c aimag(ffc(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffc(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ci = reciprocal of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci*
c    |cu(kx,ky)*s(kx,ky)|**2), where
c affp = normalization constant = nx*ny/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of field arrays, must be >= nxhm
c nyv = second dimension of field arrays, must be >= 2*nyhm
      implicit none
      integer nx, ny, nxhm, nyhm, nxvh, nyv
      real ci, wm
      complex cu, bxy, ffc
      dimension cu(3,nxvh,nyv), bxy(3,nxvh,nyv)
      dimension ffc(nxhm,nyhm)
c local data
      integer ny2, j, k, k1
      real dnx, dny, dky, ci2, at1, at2, at3
      complex zero, zt1, zt2, zt3
      double precision wp, swp
      ny2 = 2*nyhm + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      swp = 0.0d0
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 20 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k - 1)
      wp = 0.0
      do 10 j = 2, nxhm
      at1 = ci2*real(ffc(j,k))
      at2 = dnx*real(j - 1)*at1
      at3 = dky*at1
      at1 = at1*aimag(ffc(j,k))
      zt1 = cmplx(-aimag(cu(3,j,k)),real(cu(3,j,k)))
      zt2 = cmplx(-aimag(cu(2,j,k)),real(cu(2,j,k)))
      zt3 = cmplx(-aimag(cu(1,j,k)),real(cu(1,j,k)))
c ky > 0
      bxy(1,j,k) = at3*zt1
      bxy(2,j,k) = -at2*zt1
      bxy(3,j,k) = at2*zt2 - at3*zt3
c ky < 0
      zt1 = cmplx(-aimag(cu(3,j,k1)),real(cu(3,j,k1)))
      zt2 = cmplx(-aimag(cu(2,j,k1)),real(cu(2,j,k1)))
      zt3 = cmplx(-aimag(cu(1,j,k1)),real(cu(1,j,k1)))
      bxy(1,j,k1) = -at3*zt1
      bxy(2,j,k1) = -at2*zt1
      bxy(3,j,k1) = at2*zt2 + at3*zt3
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k))                         
     1   + cu(2,j,k)*conjg(cu(2,j,k)) + cu(3,j,k)*conjg(cu(3,j,k))      
     2   + cu(1,j,k1)*conjg(cu(1,j,k1)) + cu(2,j,k1)*conjg(cu(2,j,k1))  
     3   + cu(3,j,k1)*conjg(cu(3,j,k1)))
   10 continue
      swp = swp + wp
   20 continue
c mode number kx = 0
      wp = 0.0
      do 30 k = 2, nyhm
      k1 = ny2 - k
      at1 = ci2*real(ffc(1,k))
      at3 = dny*real(k - 1)*at1
      at1 = at1*aimag(ffc(1,k))
c ky > 0
      zt1 = cmplx(-aimag(cu(3,1,k)),real(cu(3,1,k)))
      zt3 = cmplx(-aimag(cu(1,1,k)),real(cu(1,1,k)))
      bxy(1,1,k) = at3*zt1
      bxy(2,1,k) = zero
      bxy(3,1,k) = -at3*zt3
c ky < 0
      bxy(2,1,k1) = zero
      if (nyhm.eq.(ny/2)) then
         bxy(1,1,k1) = zero
         bxy(3,1,k1) = zero
      else
         zt1 = cmplx(-aimag(cu(3,1,k1)),real(cu(3,1,k1)))
         zt3 = cmplx(-aimag(cu(1,1,k1)),real(cu(1,1,k1)))
         bxy(1,1,k1) = -at3*zt1
         bxy(3,1,k1) = at3*zt3
      endif
      wp = wp + at1*(cu(1,1,k)*conjg(cu(1,1,k))                         
     1   + cu(2,1,k)*conjg(cu(2,1,k)) + cu(3,1,k)*conjg(cu(3,1,k)))
   30 continue
      swp = swp + wp
c mode number ky = 0
      k1 = nyhm + 1
      wp = 0.0
      do 40 j = 2, nxhm
      at1 = ci2*real(ffc(j,1))
      at2 = dnx*real(j - 1)*at1
      at1 = at1*aimag(ffc(j,1))
      zt1 = cmplx(-aimag(cu(3,j,1)),real(cu(3,j,1)))
      zt2 = cmplx(-aimag(cu(2,j,1)),real(cu(2,j,1)))
      bxy(1,j,1) = zero
      bxy(2,j,1) = -at2*zt1
      bxy(3,j,1) = at2*zt2
c unused extra row
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      wp = wp + at1*(cu(1,j,1)*conjg(cu(1,j,1))                         
     1   + cu(2,j,1)*conjg(cu(2,j,1)) + cu(3,j,1)*conjg(cu(3,j,1)))
   40 continue
      swp = swp + wp
c mode number kx = 0, ky = 0
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
c unused extra element
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wm = real(nx*ny)*swp
      return
      end
c-----------------------------------------------------------------------
      subroutine AMAXWEL2GL(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxhm,nyhm, 
     1nxvh,nyv)
c this subroutine solves 2-1/2d maxwell's equation in fourier space for
c transverse electric and magnetic fields with periodic boundary
c conditions, using an analytic scheme due to Irving Haber
c for gridless spectral code.
c input: all, output: wf, wm, exy, bxy
c approximate flop count is: 245*nxhm*nyhm + 84*(nxhm + nyhm)
c equations being solved are:
c (c*Bn)/dt = -ick X En and
c En/dt = ick X (c*Bn) - JTn
c Note that the normalization of the E and B fields differ:
c B is normalized to the dimensionless cyclotron frequency.
c solutions are given by:
c En+1 = C*En + iS*(k X (c*Bn)) - S*JTn+1/2/c
c c*Bn+1 = C*(c*Bn) - iS*(k X En)/c + iS*T*(k X JTn+1/2)
c where En = input exy, En+1 = output exy
c Bn = input bxy, Bn+1 = output bxy, JTn+1/2 =  affp*cu*s(kx,ky)
c C = cos(k*c*dt),  S = sin(k*c*dt)/k, T = tan(k*c*dt/2)/kc
c kx = 2pi*j/nx, ky = 2pi*k/ny, k = sqrt(kx*kx + ky*ky),
c c = 1.0/ci and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c cu(i,j,k) = complex current density
c exy(i,j,k) = complex transverse electric field
c bxy(i,j,k) = complex magnetic field
c for component i, all for fourier mode (j-1,k-1)
c zero row at k = nyhm + 1 is for compatibility with FFT solvers
c when nxhm = nx/2 and nyhm = ny/2
c real(ffc(1,1)) = affp = normalization constant = nx*ny/np,
c where np=number of particles
c aimag(ffc(j,k)) = finite-size particle shape factor s,
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c for fourier mode (j-1,k-1)
c ci = reciprocal of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny**sum((1/affp)*|exy(kx,ky)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny**sum((c2/affp)*|bxy(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of field arrays, must be >= nxhm
c nyv = second dimension of field arrays, must be >= 2*nyhm
      implicit none
      integer nx, ny, nxhm, nyhm, nxvh, nyv
      real ci, dt, wf, wm
      complex exy, bxy, cu, ffc
      dimension exy(3,nxvh,nyv), bxy(3,nxvh,nyv), cu(3,nxvh,nyv)
      dimension ffc(nxhm,nyhm)
c local data
      integer ny2, j, k, k1
      real dnx, dny, dth, cc, cdth, affp, anorm, dkx, dky, dky2, t2
      real t, c, s, sc, afs, aft
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9
      double precision wp, ws, swp, sws
      if (ci.le.0.0) return
      ny2 = 2*nyhm + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dth = 0.5*dt
      cc = 1.0/ci
      cdth = cc*dth
      affp = real(ffc(1,1))
      zero = cmplx(0.0,0.0)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      sws = 0.0d0
      swp = 0.0d0
c calculate the electromagnetic fields
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 20 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k - 1)
      dky2 = dky*dky
      ws = 0.0d0
      wp = 0.0d0
      do 10 j = 2, nxhm
      dkx = dnx*real(j - 1)
      aft = affp*aimag(ffc(j,k))*ci
      t2 = sqrt(dkx*dkx + dky2)
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
c update fields for ky > 0
c calculate iB
      zt1 = cmplx(-aimag(bxy(3,j,k)),real(bxy(3,j,k)))
      zt2 = cmplx(-aimag(bxy(2,j,k)),real(bxy(2,j,k)))
      zt3 = cmplx(-aimag(bxy(1,j,k)),real(bxy(1,j,k)))
c update electric field
      zt4 = c*exy(1,j,k) + sc*(dky*zt1) - afs*cu(1,j,k)
      zt5 = c*exy(2,j,k) - sc*(dkx*zt1) - afs*cu(2,j,k)
      zt6 = c*exy(3,j,k) + sc*(dkx*zt2 - dky*zt3) - afs*cu(3,j,k)
c calculate iE
      zt1 = cmplx(-aimag(exy(3,j,k)),real(exy(3,j,k)))
      zt2 = cmplx(-aimag(exy(2,j,k)),real(exy(2,j,k)))
      zt3 = cmplx(-aimag(exy(1,j,k)),real(exy(1,j,k)))
c store electric field and calculate energy
      exy(1,j,k) = zt4
      exy(2,j,k) = zt5
      exy(3,j,k) = zt6
      ws = ws + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
c calculate ijperp
      zt7 = cmplx(-aimag(cu(3,j,k)),real(cu(3,j,k)))
      zt8 = cmplx(-aimag(cu(2,j,k)),real(cu(2,j,k)))
      zt9 = cmplx(-aimag(cu(1,j,k)),real(cu(1,j,k)))
c update magnetic field 
      zt4 = c*bxy(1,j,k) - s*dky*(zt1 - aft*zt7)
      zt5 = c*bxy(2,j,k) + s*dkx*(zt1 - aft*zt7)
      zt6 = c*bxy(3,j,k) - s*(dkx*(zt2 - aft*zt8) - dky*(zt3 - aft*zt9))
c store magnetic field and calculate energy
      bxy(1,j,k) = zt4
      bxy(2,j,k) = zt5
      bxy(3,j,k) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
c update fields for ky < 0
c calculate iB
      zt1 = cmplx(-aimag(bxy(3,j,k1)),real(bxy(3,j,k1)))
      zt2 = cmplx(-aimag(bxy(2,j,k1)),real(bxy(2,j,k1)))
      zt3 = cmplx(-aimag(bxy(1,j,k1)),real(bxy(1,j,k1)))
c update electric field
      zt4 = c*exy(1,j,k1) - sc*(dky*zt1) - afs*cu(1,j,k1)
      zt5 = c*exy(2,j,k1) - sc*(dkx*zt1) - afs*cu(2,j,k1)
      zt6 = c*exy(3,j,k1) + sc*(dkx*zt2 + dky*zt3) - afs*cu(3,j,k1)
c calculate iE
      zt1 = cmplx(-aimag(exy(3,j,k1)),real(exy(3,j,k1)))
      zt2 = cmplx(-aimag(exy(2,j,k1)),real(exy(2,j,k1)))
      zt3 = cmplx(-aimag(exy(1,j,k1)),real(exy(1,j,k1)))
c store electric field and calculate energy
      exy(1,j,k1) = zt4
      exy(2,j,k1) = zt5
      exy(3,j,k1) = zt6
      ws = ws + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
c calculate ijperp
      zt7 = cmplx(-aimag(cu(3,j,k1)),real(cu(3,j,k1)))
      zt8 = cmplx(-aimag(cu(2,j,k1)),real(cu(2,j,k1)))
      zt9 = cmplx(-aimag(cu(1,j,k1)),real(cu(1,j,k1)))
c update magnetic field 
      zt4 = c*bxy(1,j,k1) + s*dky*(zt1 - aft*zt7)
      zt5 = c*bxy(2,j,k1) + s*dkx*(zt1 - aft*zt7)
      zt6 = c*bxy(3,j,k1) - s*(dkx*(zt2 - aft*zt8) + dky*(zt3 - aft*zt9)
     1)
c store magnetic field and calculate energy
      bxy(1,j,k1) = zt4
      bxy(2,j,k1) = zt5
      bxy(3,j,k1) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
   10 continue
      swp = swp + wp
      sws = sws + ws
   20 continue
c mode number kx = 0
      ws = 0.0d0
      wp = 0.0d0
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k - 1)
      aft = affp*aimag(ffc(1,k))*ci
      t2 = sqrt(dky*dky)
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
c update fields for ky > 0
c calculate iB
      zt1 = cmplx(-aimag(bxy(3,1,k)),real(bxy(3,1,k)))
      zt3 = cmplx(-aimag(bxy(1,1,k)),real(bxy(1,1,k)))
c update electric field
      zt4 = c*exy(1,1,k) + sc*(dky*zt1) - afs*cu(1,1,k)
      zt6 = c*exy(3,1,k) - sc*(dky*zt3) - afs*cu(3,1,k)
c calculate iE
      zt1 = cmplx(-aimag(exy(3,1,k)),real(exy(3,1,k)))
      zt3 = cmplx(-aimag(exy(1,1,k)),real(exy(1,1,k)))
c store electric field and calculate energy
      exy(1,1,k) = zt4
      exy(2,1,k) = zero
      exy(3,1,k) = zt6
      ws = ws + anorm*(zt4*conjg(zt4) + zt6*conjg(zt6))
c calculate ijperp
      zt7 = cmplx(-aimag(cu(3,1,k)),real(cu(3,1,k)))
      zt9 = cmplx(-aimag(cu(1,1,k)),real(cu(1,1,k)))
c update magnetic field 
      zt4 = c*bxy(1,1,k) - s*dky*(zt1 - aft*zt7)
      zt6 = c*bxy(3,1,k) + s*dky*(zt3 - aft*zt9)
c store magnetic field and calculate energy
      bxy(1,1,k) = zt4
      bxy(2,1,k) = zero
      bxy(3,1,k) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt6*conjg(zt6))
c update fields for ky < 0
      bxy(2,1,k1) = zero
      exy(2,1,k1) = zero
      if (nyhm.eq.(ny/2)) then
         bxy(1,1,k1) = zero
         bxy(3,1,k1) = zero
         exy(1,1,k1) = zero
         exy(3,1,k1) = zero
      else
c calculate iB
         zt1 = cmplx(-aimag(bxy(3,1,k1)),real(bxy(3,1,k1)))
         zt3 = cmplx(-aimag(bxy(1,1,k1)),real(bxy(1,1,k1)))
c update electric field
         zt4 = c*exy(1,1,k1) - sc*(dky*zt1) - afs*cu(1,1,k1)
         zt6 = c*exy(3,1,k1) + sc*(dky*zt3) - afs*cu(3,1,k1)
c calculate iE
         zt1 = cmplx(-aimag(exy(3,1,k1)),real(exy(3,1,k1)))
         zt3 = cmplx(-aimag(exy(1,1,k1)),real(exy(1,1,k1)))
c store electric field and calculate energy
         exy(1,1,k1) = zt4
         exy(3,1,k1) = zt6
c calculate ijperp
         zt7 = cmplx(-aimag(cu(3,1,k1)),real(cu(3,1,k1)))
         zt9 = cmplx(-aimag(cu(1,1,k1)),real(cu(1,1,k1)))
c update magnetic field 
         zt4 = c*bxy(1,1,k1) + s*dky*(zt1 - aft*zt7)
         zt6 = c*bxy(3,1,k1) - s*dky*(zt3 - aft*zt9)
c store magnetic field and calculate energy
         bxy(1,1,k1) = zt4
         bxy(3,1,k1) = zt6
      endif
   30 continue
      swp = swp + wp
      sws = sws + ws
c mode number ky = 0
      k1 = nyhm + 1
      ws = 0.0d0
      wp = 0.0d0
      do 40 j = 2, nxhm
      dkx = dnx*real(j - 1)
      aft = affp*aimag(ffc(j,1))*ci
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
      zt1 = cmplx(-aimag(bxy(3,j,1)),real(bxy(3,j,1)))
      zt2 = cmplx(-aimag(bxy(2,j,1)),real(bxy(2,j,1)))
c update electric field
      zt5 = c*exy(2,j,1) - sc*(dkx*zt1) - afs*cu(2,j,1)
      zt6 = c*exy(3,j,1) + sc*(dkx*zt2) - afs*cu(3,j,1)
c calculate iE
      zt1 = cmplx(-aimag(exy(3,j,1)),real(exy(3,j,1)))
      zt2 = cmplx(-aimag(exy(2,j,1)),real(exy(2,j,1)))
c store electric field and calculate energy
      exy(1,j,1) = zero
      exy(2,j,1) = zt5
      exy(3,j,1) = zt6
      ws = ws + anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
c calculate ijperp
      zt7 = cmplx(-aimag(cu(3,j,1)),real(cu(3,j,1)))
      zt8 = cmplx(-aimag(cu(2,j,1)),real(cu(2,j,1)))
c update magnetic field 
      zt5 = c*bxy(2,j,1) + s*dkx*(zt1 - aft*zt7)
      zt6 = c*bxy(3,j,1) - s*dkx*(zt2 - aft*zt8)
c store magnetic field and calculate energy
      bxy(1,j,1) = zero
      bxy(2,j,1) = zt5
      bxy(3,j,1) = zt6
      wp = wp + anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
c unused extra row
      exy(1,j,k1) = zero
      exy(2,j,k1) = zero
      exy(3,j,k1) = zero
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
   40 continue
      swp = swp + wp
      sws = sws + ws
c mode number kx = 0, ky = 0
c     afs = affp*aimag(ffc(1,1))*dt
c     exy(1,1,1) = exy(1,1,1) - afs*cu(1,1,1)
c     exy(2,1,1) = exy(2,1,1) - afs*cu(2,1,1)
c     exy(3,1,1) = exy(3,1,1) - afs*cu(3,1,1)
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      exy(1,1,1) = zero
      exy(2,1,1) = zero
      exy(3,1,1) = zero
c unused extra element
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      wf = real(nx*ny)*sws
      wm = real(nx*ny)*(cc*cc)*swp
      return
      end
c-----------------------------------------------------------------------
      subroutine EMFIELD2GL(fxy,exy,ffc,isign,nxhm,nyhm,nxvh,nyv)
c this subroutine either adds complex vector fields if isign > 0
c or copies complex vector fields if isign < 0
c includes additional smoothing, for gridless spectral code.
      implicit none
      integer isign, nxhm, nyhm, nxvh, nyv
      complex fxy, exy, ffc
      dimension fxy(3,nxvh,nyv), exy(3,nxvh,nyv)
      dimension ffc(nxhm,nyhm)
c local data
      integer i, j, k, ny2, k1
      real at1
      ny2 = 2*nyhm + 2
c add the fields
      if (isign.gt.0) then
         do 30 k = 2, nyhm
         k1 = ny2 - k
         do 20 j = 1, nxhm
         at1 = aimag(ffc(j,k))
         do 10 i = 1, 3
         fxy(i,j,k) = fxy(i,j,k) + exy(i,j,k)*at1
         fxy(i,j,k1) = fxy(i,j,k1) + exy(i,j,k1)*at1
   10    continue
   20    continue
   30    continue
         k1 = nyhm + 1
         do 50 j = 1, nxhm
         at1 = aimag(ffc(j,1))
         do 40 i = 1, 3
         fxy(i,j,1) = fxy(i,j,1) + exy(i,j,1)*at1
         fxy(i,j,k1) = fxy(i,j,k1) + exy(i,j,k1)*at1
   40    continue
   50    continue
c copy the fields
      else if (isign.lt.0) then
         do 80 k = 2, nyhm
         k1 = ny2 - k
         do 70 j = 1, nxhm
         at1 = aimag(ffc(j,k))
         do 60 i = 1, 3
         fxy(i,j,k) = exy(i,j,k)*at1
         fxy(i,j,k1) = exy(i,j,k1)*at1
   60    continue
   70    continue
   80    continue
         k1 = nyhm + 1
         do 100 j = 1, nxhm
         at1 = aimag(ffc(j,1))
         do 90 i = 1, 3
         fxy(i,j,1) = exy(i,j,1)*at1
         fxy(i,j,k1) = exy(i,j,k1)*at1
   90    continue
  100    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine EPOIS23GL(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,ny,
     1nxhm,nyhm,nxvh,nyv)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c transverse electric field (or convolution of transverse electric field
c over particle shape), with periodic boundary conditions,
c for gridless spectral code.
c using algorithm described in J. Busnardo-Neto, P. L. Pritchett,
c A. T. Lin, and J. M. Dawson, J. Computational Phys. 23, 300 (1977).
c for isign = 0, input: isign,ax,ay,affp,wp0,nx,ny,nxhm,nyhm, output:ffe
c for isign /= 0, input: dcu,ffe,isign,ci,nx,ny,nxvh,nyv,nxhm,nyhm,
c output: exy,wf
c approximate flop count is: 68*nxhm*nyhm + 33*(nxhm + nyhm)
c if isign = 0, form factor array is prepared
c if isign = -1, smoothed transverse electric field is calculated
c using the equation:
c ex(kx,ky) = -ci*ci*g(kx,ky)*dcux(kx,ky)*s(kx,ky)
c ey(kx,ky) = -ci*ci*g(kx,ky)*dcuy(kx,ky)*s(kx,ky)
c ez(kx,ky) = -ci*ci*g(kx,ky)*dcuz(kx,ky)*s(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2+wp0*ci2*s(kx,ky)**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = ex(ky=pi) = ey(ky=pi) = ez(ky=pi) 
c = 0, and ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
c if isign = 1, unsmoothed transverse electric field is calculated
c using the equation:
c ex(kx,ky) = -ci*ci*g(kx,ky)*dcux(kx,ky)
c ey(kx,ky) = -ci*ci*g(kx,ky)*dcuy(kx,ky)
c ez(kx,ky) = -ci*ci*g(kx,ky)*dcuz(kx,ky)
c dcu(i,j,k) = transverse part of complex derivative of current for
c fourier mode (j-1,k-1)
c exy(1,j,k) = x component of complex transverse electric field
c exy(2,j,k) = y component of complex transverse electric field
c exy(3,j,k) = z component of complex transverse electric field
c all for fourier mode (j-1,k-1)
c zero row at k = nyhm + 1 is for compatibility with FFT solvers
c when nxhm = nx/2 and nyhm = ny/2
c aimag(ffe(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffe(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c wp0 = normalized total plasma frequency squared
c ci = reciprical of velocity of light
c transverse electric field energy is also calculated, using
c wf = nx*ny*sum((affp/((kx**2+ky**2)*ci*ci)**2)
c    |dcu(kx,ky)*s(kx,ky)|**2)
c this expression is valid only if the derivative of current is
c divergence-free
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of field arrays, must be >= nxhm
c nyv = second dimension of field arrays, must be >= 2*nyhm
      implicit none
      integer isign, nx, ny, nxhm, nyhm, nxvh, nyv
      real ax, ay, affp, wp0, ci, wf
      complex dcu, exy, ffe
      dimension dcu(3,nxvh,nyv), exy(3,nxvh,nyv)
      dimension ffe(nxhm,nyhm)
c local data
      integer ny2, j, k, k1
      real dnx, dny, ci2, wpc, dkx, dky, at1, at2, at3, at4
      complex zero
      double precision wp, swp
      ny2 = 2*nyhm + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
      if (isign.ne.0) go to 30
      wpc = wp0*ci2
c prepare form factor array
      do 20 k = 1, nyhm
      dky = dny*real(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nxhm
      dkx = dnx*real(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.0) then
         ffe(j,k) = cmplx(affp,1.)
      else
         ffe(j,k) = cmplx(affp*at4/(at3 + wpc*at4*at4),at4)
      endif
   10 continue
   20 continue
      return
c calculate smoothed transverse electric field and sum field energy
   30 if (isign.gt.0) go to 80
      swp = 0.0d0
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 50 k = 2, nyhm
      k1 = ny2 - k
      wp = 0.0
      do 40 j = 2, nxhm
      at2 = -ci2*real(ffe(j,k))
      at1 = at2*aimag(ffe(j,k))
      at2 = at2*at2
c ky > 0
      exy(1,j,k) = at1*dcu(1,j,k)
      exy(2,j,k) = at1*dcu(2,j,k)
      exy(3,j,k) = at1*dcu(3,j,k)
c ky < 0
      exy(1,j,k1) = at1*dcu(1,j,k1)
      exy(2,j,k1) = at1*dcu(2,j,k1)
      exy(3,j,k1) = at1*dcu(3,j,k1)
      wp = wp + at2*(dcu(1,j,k)*conjg(dcu(1,j,k))                       
     1 + dcu(2,j,k)*conjg(dcu(2,j,k)) + dcu(3,j,k)*conjg(dcu(3,j,k))    
     2 + dcu(1,j,k1)*conjg(dcu(1,j,k1)) + dcu(2,j,k1)*conjg(dcu(2,j,k1))
     3 + dcu(3,j,k1)*conjg(dcu(3,j,k1)))
   40 continue
      swp = swp + wp
   50 continue
c mode number kx = 0
      wp = 0.0
      do 60 k = 2, nyhm
      k1 = ny2 - k
      at2 = -ci2*real(ffe(1,k))
      at1 = at2*aimag(ffe(1,k))
      at2 = at2*at2
c ky > 0
      exy(1,1,k) = at1*dcu(1,1,k)
      exy(2,1,k) = at1*dcu(2,1,k)
      exy(3,1,k) = at1*dcu(3,1,k)
c ky < 0
      if (nyhm.eq.(ny/2)) then
         exy(1,1,k1) = zero
         exy(2,1,k1) = zero
         exy(3,1,k1) = zero
      else
         exy(1,1,k1) = at1*dcu(1,1,k1)
         exy(2,1,k1) = at1*dcu(2,1,k1)
         exy(3,1,k1) = at1*dcu(3,1,k1)
      endif
      wp = wp + at2*(dcu(1,1,k)*conjg(dcu(1,1,k))                       
     1+ dcu(2,1,k)*conjg(dcu(2,1,k)) + dcu(3,1,k)*conjg(dcu(3,1,k)))
   60 continue
      swp = swp + wp
c mode number ky = 0
      k1 = nyhm + 1
      wp = 0.0
      do 70 j = 2, nxhm
      at2 = -ci2*real(ffe(j,1))
      at1 = at2*aimag(ffe(j,1))
      at2 = at2*at2
      exy(1,j,1) = at1*dcu(1,j,1)
      exy(2,j,1) = at1*dcu(2,j,1)
      exy(3,j,1) = at1*dcu(3,j,1)
c unused extra row
      exy(1,j,k1) = zero
      exy(2,j,k1) = zero
      exy(3,j,k1) = zero
      wp = wp + at2*(dcu(1,j,1)*conjg(dcu(1,j,1))                       
     1 + dcu(2,j,1)*conjg(dcu(2,j,1)) + dcu(3,j,1)*conjg(dcu(3,j,1)))
   70 continue
      swp = swp + wp
      exy(1,1,1) = zero
      exy(2,1,1) = zero
      exy(3,1,1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      wf = real(nx*ny)*swp/real(ffe(1,1))
      return
c calculate unsmoothed transverse electric field and sum field energy
   80 swp = 0.0d0
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 100 k = 2, nyhm
      k1 = ny2 - k
      wp = 0.0
      do 90 j = 2, nxhm
      at2 = -ci2*real(ffe(j,k))
      at1 = at2*at2
c ky > 0
      exy(1,j,k) = at2*dcu(1,j,k)
      exy(2,j,k) = at2*dcu(2,j,k)
      exy(3,j,k) = at2*dcu(3,j,k)
c ky < 0
      exy(1,j,k1) = at2*dcu(1,j,k1)
      exy(2,j,k1) = at2*dcu(2,j,k1)
      exy(3,j,k1) = at2*dcu(3,j,k1)
      wp = wp + at1*(dcu(1,j,k)*conjg(dcu(1,j,k))                       
     1+ dcu(2,j,k)*conjg(dcu(2,j,k)) + dcu(3,j,k)*conjg(dcu(3,j,k))     
     2 + dcu(1,j,k1)*conjg(dcu(1,j,k1)) + dcu(2,j,k1)*conjg(dcu(2,j,k1))
     3 + dcu(3,j,k1)*conjg(dcu(3,j,k1)))
   90 continue
      swp = swp + wp
  100 continue
c mode number kx = 0
      wp = 0.0
      do 110 k = 2, nyhm
      k1 = ny2 - k
      at2 = -ci2*real(ffe(1,k))
      at1 = at2*at2
c ky < 0
      exy(1,1,k) = at2*dcu(1,1,k)
      exy(2,1,k) = at2*dcu(2,1,k)
      exy(3,1,k) = at2*dcu(3,1,k)
c ky < 0
      if (nyhm.eq.(ny/2)) then
         exy(1,1,k1) = zero
         exy(2,1,k1) = zero
         exy(3,1,k1) = zero
      else
         exy(1,1,k1) = at2*dcu(1,1,k1)
         exy(2,1,k1) = at2*dcu(2,1,k1)
         exy(3,1,k1) = at2*dcu(3,1,k1)
      endif
      wp = wp + at1*(dcu(1,1,k)*conjg(dcu(1,1,k))                       
     1+ dcu(2,1,k)*conjg(dcu(2,1,k)) + dcu(3,1,k)*conjg(dcu(3,1,k)))
  110 continue
      swp = swp + wp
c mode number ky = 0
      k1 = nyhm + 1
      wp = 0.0
      do 120 j = 2, nxhm
      at2 = -ci2*real(ffe(j,1))
      at1 = at2*at2
      exy(1,j,1) = at2*dcu(1,j,1)
      exy(2,j,1) = at2*dcu(2,j,1)
      exy(3,j,1) = at2*dcu(3,j,1)
c unused extra row
      exy(1,j,k1) = zero
      exy(2,j,k1) = zero
      exy(3,j,k1) = zero
      wp = wp + at1*(dcu(1,j,1)*conjg(dcu(1,j,1))                       
     1 + dcu(2,j,1)*conjg(dcu(2,j,1)) + dcu(3,j,1)*conjg(dcu(3,j,1)))
  120 continue
      swp = swp + wp
c mode number kx = 0, ky = 0
      exy(1,1,1) = zero
      exy(2,1,1) = zero
      exy(3,1,1) = zero
c unused extra element
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      wf = real(nx*ny)*swp/real(ffe(1,1))
      return
      end
c-----------------------------------------------------------------------
      subroutine EVFIELD23GL(fxy,gxy,sctx,xp,yp,nx,ny,nxhm,nyhm,nxvh,nyv
     1)
c for 2-1/2d code, this subroutine calculates real space 2d field gxy at
c location xp, yp, from complex fourier co-efficients fxy
c input: all, output: gxy, sctx
c equations used are:
c gxy(1) = sum(fxy(1,n,m)*exp(sqrt(-1)*2*n*pi*xp/nx)*
c                             exp(sqrt(-1)*2*m*pi*yp/ny))
c gxy(2) = sum(fxy(2,n,m)*exp(sqrt(-1)*2*n*pi*xp/nx)*
c                             exp(sqrt(-1)*2*m*pi*yp/ny))
c gxy(3) = sum(fxy(3,n,m)*exp(sqrt(-1)*2*n*pi*xp/nx)*
c                             exp(sqrt(-1)*2*m*pi*y/ny))
c fxy(1,n,m) = x component of field at fourier mode (n-1,m-1)
c fxy(2,n,m) = y component of field at fourier mode (n-1,m-1)
c fxy(3,n,m) = z component of field at fourier mode (n-1,m-1)
c if nxhm=nx/2 and nyhm=ny/2, nx/2 and ny/2 modes are included, in a
c format compatible with the real to complex FFT used in traditional PIC
c gxy(1) = x component of field at location (xp,yp)
c gxy(2) = y component of field at location (xp,yp)
c gxy(3) = z component of field at location (xp,yp)
c sctx = scratch array for sines and cosines
c xp/yp = position to evaluate field
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of field arrays, must be >= nxhm
c nyv = second dimension of field arrays, must be >= 2*nyhm
      implicit none
      integer nx, ny, nxhm, nyhm, nxvh, nyv
      real xp, yp
      real gxy
      complex fxy, sctx
      dimension gxy(3)
      dimension fxy(3,nxvh,nyv), sctx(nxvh)
c local data
      integer j, k, k1, nxh, nyh, ny2
      real dnx, dny, dkx, dky, at1, at3
      double precision exl, eyl, ezl, ex, ey, ez
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = 2*nyhm + 2
c find electric field
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*xp
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k-1)*yp
      zt1 = cmplx(cos(dky),sin(dky))
      zt2 = conjg(zt1)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      do 20 j = 2, nxhm
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      exl = exl + (real(fxy(1,j,k)*zt3) + real(fxy(1,j,k1)*zt4))
      eyl = eyl + (real(fxy(2,j,k)*zt3) + real(fxy(2,j,k1)*zt4))
      ezl = ezl + (real(fxy(3,j,k)*zt3) + real(fxy(3,j,k1)*zt4))
   20 continue
c mode numbers kx = 0
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         zt2 = zt2*real(sctx(nxhm))
         exl = exl + (real(fxy(1,1,k)*zt1) + real(fxy(1,1,k1)*zt2))
         eyl = eyl + (real(fxy(2,1,k)*zt1) + real(fxy(2,1,k1)*zt2))
         ezl = ezl + (real(fxy(3,1,k)*zt1) + real(fxy(3,1,k1)*zt2))
      else
         exl = exl + real(fxy(1,1,k)*zt1)
         eyl = eyl + real(fxy(2,1,k)*zt1)
         ezl = ezl + real(fxy(3,1,k)*zt1)
      endif
      ex = ex + exl
      ey = ey + eyl
      ez = ez + ezl
   30 continue
c mode numbers ky = 0
      k1 = nyhm + 1
      dky = dny*real(nyh)*yp
      at1 = cos(dky)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      do 40 j = 2, nxhm
      zt3 = sctx(j-1)
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         zt1 = at1*zt3
         exl = exl + (real(fxy(1,j,1)*zt3) + real(fxy(1,j,k1)*zt1))
         eyl = eyl + (real(fxy(2,j,1)*zt3) + real(fxy(2,j,k1)*zt1))
         ezl = ezl + (real(fxy(3,j,1)*zt3) + real(fxy(3,j,k1)*zt1))
      else
         exl = exl + (real(fxy(1,j,1)*zt3))
         eyl = eyl + (real(fxy(2,j,1)*zt3))
         ezl = ezl + (real(fxy(3,j,1)*zt3))
      endif
   40 continue
      ex = 2.0d0*(ex + exl)
      ey = 2.0d0*(ey + eyl)
      ez = 2.0d0*(ez + ezl)
      at3 = real(sctx(nxhm))
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         ex = ex + (real(fxy(1,1,1)) + aimag(fxy(1,1,1))*at3)
         ey = ey + (real(fxy(2,1,1)) + aimag(fxy(2,1,1))*at3)
         ez = ez + (real(fxy(3,1,1)) + aimag(fxy(3,1,1))*at3)
      else
         ex = ex + real(fxy(1,1,1))
         ey = ey + real(fxy(2,1,1))
         ez = ez + real(fxy(3,1,1))
      endif
c special case to match conventional PIC code
      if (nxhm.eq.nyh) then
         if (nxhm.eq.nxh) then
            ex = ex + (real(fxy(1,1,k1)) + aimag(fxy(1,1,k1))*at3)*at1
            ey = ey + (real(fxy(2,1,k1)) + aimag(fxy(2,1,k1))*at3)*at1
            ez = ez + (real(fxy(3,1,k1)) + aimag(fxy(3,1,k1))*at3)*at1
         else
            ex = ex + real(fxy(1,1,k1))*at1
            ey = ey + real(fxy(2,1,k1))*at1
            ez = ez + real(fxy(3,1,k1))*at1
         endif
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
      subroutine MJPOST2GL(part,amu,sctx,qm,nop,idimp,nx,ny,nxhm,nyhm,  
     1nxvh,nyv)
c for 2-1/2d code, this subroutine calculates particle momentum flux
c using gridless spectral version, periodic boundaries
c baseline scalar version
c 44*nxhm*nyhm + 10 flops/particle, 2*NX*NY loads and stores
c plus (nxhm + nyhm) sines and cosines/particle
c input: all, output: amu, sctx
c momentum flux is calculated from the expression:
c amu(i,n,m) = sum(qci*exp(-sqrt(-1)*2*n*pi*x/nx)*
c                  exp(-sqrt(-1)*2*m*pi*y/ny))
c where qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c and vj = vj(t-dt/2), and vk = vk(t-dt/2)
c if nxhm=nx/2 and nyhm=ny/2, nx/2 and ny/2 modes are calculated, in a
c format compatible with the real to complex FFT used in traditional PIC
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = x velocity of particle n at t - dt/2
c part(4,n) = y velocity of particle n at t - dt/2
c part(5,n) = z velocity of particle n at t - dt/2
c amu(i,j,k) = ith component of momentum flux
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of flux array, must be >= nxhm
c nyv = second dimension of flux array, must be >= 2*nyhm
      implicit none
      integer nop, idimp, nx, ny, nxhm, nyhm, nxvh, nyv
      real part, qm
      complex amu, sctx
      dimension part(idimp,nop), amu(4,nxvh,nyv), sctx(nxvh)
c local data
      integer i, j, k, k1, nxh, nyh, ny2
      real qmn, dnx, dny, dkx, dky, at1, at3
      real vx, vy, vz, v1, v2, v3, v4
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = 2*nyhm + 2
      qmn = qm/real(nx*ny)
c find fourier components
      do 50 i = 1, nop
      vx = part(3,i)
      vy = part(4,i)
      vz = part(5,i)
      v1 = qmn*(vx*vx - vy*vy)
      v2 = qmn*vx*vy
      v3 = qmn*vz*vx
      v4 = qmn*vz*vy
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k-1)*part(2,i)
      zt1 = cmplx(cos(dky),-sin(dky))
      zt2 = conjg(zt1)
      do 20 j = 2, nxhm
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      amu(1,j,k) = amu(1,j,k) + v1*zt3
      amu(2,j,k) = amu(2,j,k) + v2*zt3
      amu(3,j,k) = amu(3,j,k) + v3*zt3
      amu(4,j,k) = amu(4,j,k) + v4*zt3
      amu(1,j,k1) = amu(1,j,k1) + v1*zt4
      amu(2,j,k1) = amu(2,j,k1) + v2*zt4
      amu(3,j,k1) = amu(3,j,k1) + v3*zt4
      amu(4,j,k1) = amu(4,j,k1) + v4*zt4
   20 continue
c mode number kx = 0
      amu(1,1,k) = amu(1,1,k) + v1*zt1
      amu(2,1,k) = amu(2,1,k) + v2*zt1
      amu(3,1,k) = amu(3,1,k) + v3*zt1
      amu(4,1,k) = amu(4,1,k) + v4*zt1
c special case to match conventional PIC code
      if (nxhm.eq.nxh) zt2 = zt2*real(sctx(nxhm))
      amu(1,1,k1) = amu(1,1,k1) + v1*zt2
      amu(2,1,k1) = amu(2,1,k1) + v2*zt2
      amu(3,1,k1) = amu(3,1,k1) + v3*zt2
      amu(4,1,k1) = amu(4,1,k1) + v4*zt2
   30 continue
c mode number ky = 0
      k1 = nyhm + 1
      dky = dny*real(nyh)*part(2,i)
      at1 = cos(dky)
      do 40 j = 2, nxhm
      zt3 = sctx(j-1)
      amu(1,j,1) = amu(1,j,1) + v1*zt3
      amu(2,j,1) = amu(2,j,1) + v2*zt3
      amu(3,j,1) = amu(3,j,1) + v3*zt3
      amu(4,j,1) = amu(4,j,1) + v4*zt3
      zt1 = at1*zt3
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         zt1 = at1*zt3
         amu(1,j,k1) = amu(1,j,k1) + v1*zt1
         amu(2,j,k1) = amu(2,j,k1) + v2*zt1
         amu(3,j,k1) = amu(3,j,k1) + v3*zt1
         amu(4,j,k1) = amu(4,j,k1) + v4*zt1
      endif
   40 continue
      at3 = real(sctx(nxh))
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         amu(1,1,1) = cmplx(real(amu(1,1,1))+v1,aimag(amu(1,1,1))+v1*at3
     1)
         amu(2,1,1) = cmplx(real(amu(2,1,1))+v2,aimag(amu(2,1,1))+v2*at3
     1)
         amu(3,1,1) = cmplx(real(amu(3,1,1))+v3,aimag(amu(3,1,1))+v3*at3
     1)
         amu(4,1,1) = cmplx(real(amu(4,1,1))+v4,aimag(amu(4,1,1))+v4*at3
     1)
      else
         amu(1,1,1) = cmplx(real(amu(1,1,1))+v1,0.0)
         amu(2,1,1) = cmplx(real(amu(2,1,1))+v2,0.0)
         amu(3,1,1) = cmplx(real(amu(3,1,1))+v3,0.0)
         amu(4,1,1) = cmplx(real(amu(4,1,1))+v4,0.0)
      endif
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         if (nxhm.eq.nxh) then
            at3 = at1*at3
            amu(1,1,k1) = cmplx(real(amu(1,1,k1))+v1*at1,               
     1aimag(amu(1,1,k1))+v1*at3)
            amu(2,1,k1) = cmplx(real(amu(2,1,k1))+v2*at1,               
     1aimag(amu(2,1,k1))+v2*at3)
            amu(3,1,k1) = cmplx(real(amu(3,1,k1))+v3*at1,               
     1aimag(amu(3,1,k1))+v3*at3)
            amu(4,1,k1) = cmplx(real(amu(4,1,k1))+v4*at1,               
     1aimag(amu(4,1,k1))+v4*at3)
         else
            amu(1,1,k1) = cmplx(real(amu(1,1,k1))+v1*at1,0.0)
            amu(2,1,k1) = cmplx(real(amu(2,1,k1))+v2*at1,0.0)
            amu(3,1,k1) = cmplx(real(amu(3,1,k1))+v3*at1,0.0)
            amu(4,1,k1) = cmplx(real(amu(4,1,k1))+v4*at1,0.0)
         endif
      endif
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine RMJPOST2GL(part,amu,sctx,qm,ci,nop,idimp,nx,ny,nxhm,   
     1nyhm,nxvh,nyv)
c for 2-1/2d code, this subroutine calculates particle momentum flux
c using gridless spectral version for relativistic particles,
c and periodic boundaries
c baseline scalar version
c 44*nxhm*nyhm + 10 flops/particle, 2*NX*NY loads and stores
c plus (nxhm + nyhm) sines and cosines/particle
c input: all, output: amu, sctx
c momentum flux is calculated from the expression:
c amu(i,n,m) = sum(qci*exp(-sqrt(-1)*2*n*pi*x/nx)*
c                  exp(-sqrt(-1)*2*m*pi*y/ny))
c and qci = qm*pj*pk*gami2, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where pj = pj(t-dt/2) and pk = pk(t-dt/2)
c where gami2 = 1./(1.+sum(pi**2)*ci*ci)
c if nxhm=nx/2 and nyhm=ny/2, nx/2 and ny/2 modes are calculated, in a
c format compatible with the real to complex FFT used in traditional PIC
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = x momentum of particle n at t - dt/2
c part(4,n) = y momentum of particle n at t - dt/2
c part(5,n) = z momentum of particle n at t - dt/2
c amu(i,j,k) = ith component of momentum flux
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of flux array, must be >= nxhm
c nyv = second dimension of flux array, must be >= 2*nyhm
      implicit none
      integer nop, idimp, nx, ny, nxhm, nyhm, nxvh, nyv
      real part, qm, ci
      complex amu, sctx
      dimension part(idimp,nop), amu(4,nxvh,nyv), sctx(nxvh)
c local data
      integer i, j, k, k1, nxh, nyh, ny2
      real qmn, ci2, gami2, dnx, dny, dkx, dky, at1, at3
      real vx, vy, vz, p2, v1, v2, v3, v4
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = 2*nyhm + 2
      qmn = qm/real(nx*ny)
      ci2 = ci*ci
c find fourier components
      do 50 i = 1, nop
      vx = part(3,i)
      vy = part(4,i)
      vz = part(5,i)
      p2 = vx*vx + vy*vy + vz*vz
      gami2 = qmn/(1.0 + p2*ci2)
      v1 = (vx*vx - vy*vy)*gami2
      v2 = (vx*vy)*gami2
      v3 = (vz*vx)*gami2
      v4 = (vz*vy)*gami2
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k-1)*part(2,i)
      zt1 = cmplx(cos(dky),-sin(dky))
      zt2 = conjg(zt1)
      do 20 j = 2, nxhm
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      amu(1,j,k) = amu(1,j,k) + v1*zt3
      amu(2,j,k) = amu(2,j,k) + v2*zt3
      amu(3,j,k) = amu(3,j,k) + v3*zt3
      amu(4,j,k) = amu(4,j,k) + v4*zt3
      amu(1,j,k1) = amu(1,j,k1) + v1*zt4
      amu(2,j,k1) = amu(2,j,k1) + v2*zt4
      amu(3,j,k1) = amu(3,j,k1) + v3*zt4
      amu(4,j,k1) = amu(4,j,k1) + v4*zt4
   20 continue
c mode number kx = 0
      amu(1,1,k) = amu(1,1,k) + v1*zt1
      amu(2,1,k) = amu(2,1,k) + v2*zt1
      amu(3,1,k) = amu(3,1,k) + v3*zt1
      amu(4,1,k) = amu(4,1,k) + v4*zt1
c special case to match conventional PIC code
      if (nxhm.eq.nxh) zt2 = zt2*real(sctx(nxhm))
      amu(1,1,k1) = amu(1,1,k1) + v1*zt2
      amu(2,1,k1) = amu(2,1,k1) + v2*zt2
      amu(3,1,k1) = amu(3,1,k1) + v3*zt2
      amu(4,1,k1) = amu(4,1,k1) + v4*zt2
   30 continue
c mode number ky = 0
      k1 = nyhm + 1
      dky = dny*real(nyh)*part(2,i)
      at1 = cos(dky)
      do 40 j = 2, nxhm
      zt3 = sctx(j-1)
      amu(1,j,1) = amu(1,j,1) + v1*zt3
      amu(2,j,1) = amu(2,j,1) + v2*zt3
      amu(3,j,1) = amu(3,j,1) + v3*zt3
      amu(4,j,1) = amu(4,j,1) + v4*zt3
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         zt1 = at1*zt3
         amu(1,j,k1) = amu(1,j,k1) + v1*zt1
         amu(2,j,k1) = amu(2,j,k1) + v2*zt1
         amu(3,j,k1) = amu(3,j,k1) + v3*zt1
         amu(4,j,k1) = amu(4,j,k1) + v4*zt1
      endif
   40 continue
      at3 = real(sctx(nxh))
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         amu(1,1,1) = cmplx(real(amu(1,1,1))+v1,aimag(amu(1,1,1))+v1*at3
     1)
         amu(2,1,1) = cmplx(real(amu(2,1,1))+v2,aimag(amu(2,1,1))+v2*at3
     1)
         amu(3,1,1) = cmplx(real(amu(3,1,1))+v3,aimag(amu(3,1,1))+v3*at3
     1)
         amu(4,1,1) = cmplx(real(amu(4,1,1))+v4,aimag(amu(4,1,1))+v4*at3
     1)
      else
         amu(1,1,1) = cmplx(real(amu(1,1,1))+v1,0.0)
         amu(2,1,1) = cmplx(real(amu(2,1,1))+v2,0.0)
         amu(3,1,1) = cmplx(real(amu(3,1,1))+v3,0.0)
         amu(4,1,1) = cmplx(real(amu(4,1,1))+v4,0.0)
      endif
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         if (nxhm.eq.nxh) then
            at3 = at1*at3
            amu(1,1,k1) = cmplx(real(amu(1,1,k1))+v1*at1,               
     1aimag(amu(1,1,k1))+v1*at3)
            amu(2,1,k1) = cmplx(real(amu(2,1,k1))+v2*at1,               
     1aimag(amu(2,1,k1))+v2*at3)
            amu(3,1,k1) = cmplx(real(amu(3,1,k1))+v3*at1,               
     1aimag(amu(3,1,k1))+v3*at3)
            amu(4,1,k1) = cmplx(real(amu(4,1,k1))+v4*at1,               
     1aimag(amu(4,1,k1))+v4*at3)
         else
            amu(1,1,k1) = cmplx(real(amu(1,1,k1))+v1*at1,0.0)
            amu(2,1,k1) = cmplx(real(amu(2,1,k1))+v2*at1,0.0)
            amu(3,1,k1) = cmplx(real(amu(3,1,k1))+v3*at1,0.0)
            amu(4,1,k1) = cmplx(real(amu(4,1,k1))+v4*at1,0.0)
         endif
      endif
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DCUPERP23GL(dcu,amu,nx,ny,nxhm,nyhm,nxvh,nyv)
c this subroutine calculates transverse part of the derivative of
c the current density from the momentum flux
c in 2-1/2d with periodic boundary conditions,
c for gridless spectral code.
c the derivative of the current is calculated using the equations:
c dcu(1,kx,ky) = -sqrt(-1)*(kx*vx*vx+ky*vx*vy)
c dcu(2,kx,ky) = -sqrt(-1)*(kx*vx*vy+ky*vy*vy)
c dcu(3,kx,ky) = -sqrt(-1)*(kx*vx*vz+ky*vy*vz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for dcu(i,kx=pi) = dcu(i,ky=pi) =  dcu(i,kx=0,ky=0) = 0.
c the transverse part is calculated using the equation:
c dcu(1,kx,ky) = dcu(1,kx,ky)-kx*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c dcu(2,kx,ky) = dcu(2,kx,ky)-ky*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c on output:
c dcu(i,j,k) = transverse part of complex derivative of current for
c fourier mode (j-1,k-1)
c amu(1,j,k) = xx component of complex momentum flux
c amu(2,j,k) = xy component of complex momentum flux
c amu(3,j,k) = zx component of complex momentum flux
c amu(4,j,k) = zy component of complex momentum flux
c all for fourier mode (j-1,k-1)
c zero row at k = nyhm + 1 is for compatibility with FFT solvers
c when nxhm = nx/2 and nyhm = ny/2
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = second dimension of field arrays, must be >= nxhm
c nyv = third dimension of field arrays, must be >= 2*nyhm
      implicit none
      integer nx, ny, nxhm, nyhm, nxvh, nyv
      complex dcu, amu
      dimension dcu(3,nxvh,nyv), amu(4,nxvh,nyv)
c local data
      integer ny2, j, k, k1
      real dnx, dny, dky, dky2, dkx, dkx2, dkxy, dkxy2, at1
      complex zero, zt1, zt2, zt3
      ny2 = 2*nyhm + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 20 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nxhm
      dkx = dnx*real(j - 1)
      dkx2 = dkx*dkx
      dkxy = dkx*dky
      dkxy2 = dky2 - dkx2
      at1 = 1.0/(dkx2 + dky2)
c ky > 0
      zt1 = cmplx(aimag(amu(1,j,k)),-real(amu(1,j,k)))
      zt2 = cmplx(aimag(amu(2,j,k)),-real(amu(2,j,k)))
      zt3 = at1*(dkxy*zt1 + dkxy2*zt2)
      dcu(1,j,k) = dky*zt3
      dcu(2,j,k) = -dkx*zt3
      zt1 = cmplx(aimag(amu(3,j,k)),-real(amu(3,j,k)))
      zt2 = cmplx(aimag(amu(4,j,k)),-real(amu(4,j,k)))
      dcu(3,j,k) = dkx*zt1 + dky*zt2
c ky < 0
      zt1 = cmplx(aimag(amu(1,j,k1)),-real(amu(1,j,k1)))
      zt2 = cmplx(aimag(amu(2,j,k1)),-real(amu(2,j,k1)))
      zt3 = at1*(dkxy*zt1 - dkxy2*zt2)
      dcu(1,j,k1) = dky*zt3
      dcu(2,j,k1) = dkx*zt3
      zt1 = cmplx(aimag(amu(3,j,k1)),-real(amu(3,j,k1)))
      zt2 = cmplx(aimag(amu(4,j,k1)),-real(amu(4,j,k1)))
      dcu(3,j,k1) = dkx*zt1 - dky*zt2
   10 continue
   20 continue
c mode number kx = 0
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k - 1)
c ky > 0
      zt2 = cmplx(aimag(amu(2,1,k)),-real(amu(2,1,k)))
      dcu(1,1,k) = dky*zt2
      dcu(2,1,k) = zero
      zt2 = cmplx(aimag(amu(4,1,k)),-real(amu(4,1,k)))
      dcu(3,1,k) = dky*zt2
c ky < 0
      dcu(2,1,k1) = zero
      if (nyhm.eq.(ny/2)) then
         dcu(1,1,k1) = zero
         dcu(3,1,k1) = zero
      else
         zt2 = cmplx(aimag(amu(2,1,k1)),-real(amu(2,1,k1)))
         dcu(1,1,k1) = -dky*zt2
         zt2 = cmplx(aimag(amu(4,1,k1)),-real(amu(4,1,k1)))
         dcu(3,1,k1) = -dky*zt2
      endif
   30 continue
c mode number ky = 0
      k1 = nyhm + 1
      do 40 j = 2, nxhm
      dkx = dnx*real(j - 1)
      zt2 = cmplx(aimag(amu(2,j,1)),-real(amu(2,j,1)))
      dcu(1,j,1) = zero
      dcu(2,j,1) = dkx*zt2
      zt1 = cmplx(aimag(amu(3,j,1)),-real(amu(3,j,1)))
      dcu(3,j,1) = dkx*zt1
c unused extra row
      dcu(1,j,k1) = zero
      dcu(2,j,k1) = zero
      dcu(3,j,k1) = zero
   40 continue
c mode numbers kx = 0, ky = 0
      dcu(1,1,1) = zero
      dcu(2,1,1) = zero
      dcu(3,1,1) = zero
c unused extra element
      dcu(1,1,k1) = zero
      dcu(2,1,k1) = zero
      dcu(3,1,k1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine ADCUPERP23GL(dcu,amu,nx,ny,nxhm,nyhm,nxvh,nyv)
c this subroutine calculates transverse part of the derivative of
c the current density from the momentum flux and acceleration density
c in 2-1/2d with periodic boundary condition,
c for gridless spectral code.
c the derivative of the current is calculated using the equations:
c dcu(1,kx,ky) = dcu(1,kx,ky)-sqrt(-1)*(kx*vx*vx+ky*vx*vy)
c dcu(2,kx,ky) = dcu(2,kx,ky)-sqrt(-1)*(kx*vx*vy+ky*vy*vy)
c dcu(3,kx,ky) = dcu(3,kx,ky)-sqrt(-1)*(kx*vx*vz+ky*vy*vz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for dcu(i,kx=pi) = dcu(i,ky=pi) =  dcu(i,kx=0,ky=0) = 0.
c the transverse part is calculated using the equation:
c dcu(1,kx,ky) = dcu(1,kx,ky)-kx*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c dcu(2,kx,ky) = dcu(2,kx,ky)-ky*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c on input:
c dcu(i,j,k) = complex acceleration density for fourier mode (j-1,k-1)
c on output:
c dcu(i,j,k) = transverse part of complex derivative of current for
c fourier mode (j-1,k-1)
c amu(1,j,k) = xx component of complex momentum flux
c amu(2,j,k) = xy component of complex momentum flux
c amu(3,j,k) = zx component of complex momentum flux
c amu(4,j,k) = zy component of complex momentum flux
c all for fourier mode (j-1,k-1)
c zero row at k = nyhm + 1 is for compatibility with FFT solvers
c when nxhm = nx/2 and nyhm = ny/2
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = second dimension of field arrays, must be >= nxhm
c nyv = third dimension of field arrays, must be >= 2*nyhm
      implicit none
      integer nx, ny, nxhm, nyhm, nxvh, nyv
      complex dcu, amu
      dimension dcu(3,nxvh,nyv), amu(4,nxvh,nyv)
c local data
      integer ny2, j, k, k1
      real dnx, dny, dky, dky2, dkx, dkx2, dkxy, dkxy2, at1
      complex zero, zt1, zt2, zt3
      ny2 = 2*nyhm + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 20 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nxhm
      dkx = dnx*real(j - 1)
      dkx2 = dkx*dkx
      dkxy = dkx*dky
      dkxy2 = dky2 - dkx2
      at1 = 1.0/(dkx2 + dky2)
c ky > 0
      zt1 = cmplx(aimag(amu(1,j,k)),-real(amu(1,j,k)))
      zt2 = cmplx(aimag(amu(2,j,k)),-real(amu(2,j,k)))
      zt3 = at1*(dky*dcu(1,j,k) - dkx*dcu(2,j,k) + dkxy*zt1 + dkxy2*zt2)
      dcu(1,j,k) = dky*zt3
      dcu(2,j,k) = -dkx*zt3
      zt1 = cmplx(aimag(amu(3,j,k)),-real(amu(3,j,k)))
      zt2 = cmplx(aimag(amu(4,j,k)),-real(amu(4,j,k)))
      dcu(3,j,k) = dcu(3,j,k) + dkx*zt1 + dky*zt2
c ky < 0
      zt1 = cmplx(aimag(amu(1,j,k1)),-real(amu(1,j,k1)))
      zt2 = cmplx(aimag(amu(2,j,k1)),-real(amu(2,j,k1)))
      zt3 = at1*(dky*dcu(1,j,k1) + dkx*dcu(2,j,k1) + dkxy*zt1 - dkxy2*zt
     12)
      dcu(1,j,k1) = dky*zt3
      dcu(2,j,k1) = dkx*zt3
      zt1 = cmplx(aimag(amu(3,j,k1)),-real(amu(3,j,k1)))
      zt2 = cmplx(aimag(amu(4,j,k1)),-real(amu(4,j,k1)))
      dcu(3,j,k1) = dcu(3,j,k1) + dkx*zt1 - dky*zt2
   10 continue
   20 continue
c mode number kx = 0
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k - 1)
c ky > 0
      zt2 = cmplx(aimag(amu(2,1,k)),-real(amu(2,1,k)))
      dcu(1,1,k) = dcu(1,1,k) + dky*zt2
      dcu(2,1,k) = zero
      zt2 = cmplx(aimag(amu(4,1,k)),-real(amu(4,1,k)))
      dcu(3,1,k) = dcu(3,1,k) + dky*zt2
c ky < 0
      dcu(2,1,k1) = zero
      if (nyhm.eq.(ny/2)) then
         dcu(1,1,k1) = zero
         dcu(3,1,k1) = zero
      else
         zt2 = cmplx(aimag(amu(2,1,k1)),-real(amu(2,1,k1)))
         dcu(1,1,k1) = dcu(1,1,k1) - dky*zt2
         zt2 = cmplx(aimag(amu(4,1,k1)),-real(amu(4,1,k1)))
         dcu(3,1,k1) = dcu(3,1,k1) - dky*zt2
      endif
   30 continue
c mode number ky = 0
      k1 = nyhm + 1
      do 40 j = 2, nxhm
      dkx = dnx*real(j - 1)
      zt2 = cmplx(aimag(amu(2,j,1)),-real(amu(2,j,1)))
      dcu(1,j,1) = zero
      dcu(2,j,1) = dcu(2,j,1) + dkx*zt2
      zt1 = cmplx(aimag(amu(3,j,1)),-real(amu(3,j,1)))
      dcu(3,j,1) = dcu(3,j,1) + dkx*zt1
c unused extra row
      dcu(1,j,k1) = zero
      dcu(2,j,k1) = zero
      dcu(3,j,k1) = zero
   40 continue
c mode numbers kx = 0, ky = 0
      dcu(1,1,1) = zero
      dcu(2,1,1) = zero
c unused extra element
      dcu(1,1,k1) = zero
      dcu(2,1,k1) = zero
      dcu(3,1,k1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine MAXWEL2GL(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxhm,nyhm,  
     1nxvh,nyv)
c this subroutine solves 2-1/2d maxwell's equation in fourier space for
c transverse electric and magnetic fields with periodic boundary
c conditions, for gridless spectral code.
c input: all, output: wf, wm, exy, bxy
c approximate flop count is: 286*nxhm*nyhm + 84*(nxhm + nyhm)
c the magnetic field is first updated half a step using the equations:
c bx(kx,ky) = bx(kx,ky) - .5*dt*sqrt(-1)*ky*ez(kx,ky)
c by(kx,ky) = by(kx,ky) + .5*dt*sqrt(-1)*kx*ez(kx,ky)
c bz(kx,ky) = bz(kx,ky) - .5*dt*sqrt(-1)*(kx*ey(kx,ky)-ky*ex(kx,ky))
c the electric field is then updated a whole step using the equations:
c ex(kx,ky) = ex(kx,ky) + c2*dt*sqrt(-1)*ky*bz(kx,ky)
c                       - affp*dt*cux(kx,ky)*s(kx,ky)
c ey(kx,ky) = ey(kx,ky) - c2*dt*sqrt(-1)*kx*bz(kx,ky)
c                       - affp*dt*cuy(kx,ky)*s(kx,ky)
c ez(kx,ky) = ez(kx,ky) + c2*dt*sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))
c                       - affp*dt*cuz(kx,ky)*s(kx,ky)
c the magnetic field is finally updated the remaining half step with
c the new electric field and the previous magnetic field equations.
c Note that the normalization of the E and B fields differ:
c B is normalized to the dimensionless cyclotron frequency.
c where kx = 2pi*j/nx, ky = 2pi*k/ny, c2 = 1./(ci*ci)
c and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c j,k = fourier mode numbers, except for
c ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
c ex(ky=pi) = ey(ky=pi) = ex(ky=pi) = 0,
c ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
c and similarly for bx, by, bz.
c cu(i,j,k) = complex current density
c exy(i,j,k) = complex transverse electric field
c bxy(i,j,k) = complex magnetic field
c for component i, all for fourier mode (j-1,k-1)
c zero row at k = nyhm + 1 is for compatibility with FFT solvers
c when nxhm = nx/2 and nyhm = ny/2
c real(ffc(1,1)) = affp = normalization constant = nx*ny/np,
c where np=number of particles
c aimag(ffc(j,k)) = finite-size particle shape factor s,
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c for fourier mode (j-1,k-1)
c ci = reciprocal of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny**sum((1/affp)*|exy(kx,ky)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny**sum((c2/affp)*|bxy(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of field arrays, must be >= nxhm
c nyv = second dimension of field arrays, must be >= 2*nyhm
      implicit none
      integer nx, ny, nxhm, nyhm, nxvh, nyv
      real ci, dt, wf, wm
      complex exy, bxy, cu, ffc
      dimension exy(3,nxvh,nyv), bxy(3,nxvh,nyv), cu(3,nxvh,nyv)
      dimension ffc(nxhm,nyhm)
c local data
      integer ny2, j, k, k1
      real dnx, dny, dth, c2, cdt, affp, anorm, dkx, dky, afdt, adt
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9
      double precision wp, ws
      if (ci.le.0.0) return
      ny2 = 2*nyhm + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dth = 0.5*dt
      c2 = 1.0/(ci*ci)
      cdt = c2*dt
      affp = real(ffc(1,1))
      adt = affp*dt
      zero = cmplx(0.0,0.0)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
c calculate the electromagnetic fields
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 20 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k - 1)
      do 10 j = 2, nxhm
      dkx = dnx*real(j - 1)
      afdt = adt*aimag(ffc(j,k))
c update magnetic field half time step, ky > 0
      zt1 = cmplx(-aimag(exy(3,j,k)),real(exy(3,j,k)))
      zt2 = cmplx(-aimag(exy(2,j,k)),real(exy(2,j,k)))
      zt3 = cmplx(-aimag(exy(1,j,k)),real(exy(1,j,k)))
      zt4 = bxy(1,j,k) - dth*(dky*zt1)
      zt5 = bxy(2,j,k) + dth*(dkx*zt1)
      zt6 = bxy(3,j,k) - dth*(dkx*zt2 - dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,j,k) + cdt*(dky*zt1) - afdt*cu(1,j,k)
      zt8 = exy(2,j,k) - cdt*(dkx*zt1) - afdt*cu(2,j,k)
      zt9 = exy(3,j,k) + cdt*(dkx*zt2 - dky*zt3) - afdt*cu(3,j,k)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exy(1,j,k) = zt7
      exy(2,j,k) = zt8
      exy(3,j,k) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      zt4 = zt4 - dth*(dky*zt1)
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2 - dky*zt3)
      bxy(1,j,k) = zt4
      bxy(2,j,k) = zt5
      bxy(3,j,k) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
c update magnetic field half time step, ky < 0
      zt1 = cmplx(-aimag(exy(3,j,k1)),real(exy(3,j,k1)))
      zt2 = cmplx(-aimag(exy(2,j,k1)),real(exy(2,j,k1)))
      zt3 = cmplx(-aimag(exy(1,j,k1)),real(exy(1,j,k1)))
      zt4 = bxy(1,j,k1) + dth*(dky*zt1)
      zt5 = bxy(2,j,k1) + dth*(dkx*zt1)
      zt6 = bxy(3,j,k1) - dth*(dkx*zt2 + dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,j,k1) - cdt*(dky*zt1) - afdt*cu(1,j,k1)
      zt8 = exy(2,j,k1) - cdt*(dkx*zt1) - afdt*cu(2,j,k1)
      zt9 = exy(3,j,k1) + cdt*(dkx*zt2 + dky*zt3) - afdt*cu(3,j,k1)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exy(1,j,k1) = zt7
      exy(2,j,k1) = zt8
      exy(3,j,k1) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      zt4 = zt4 + dth*(dky*zt1)
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2 + dky*zt3)
      bxy(1,j,k1) = zt4
      bxy(2,j,k1) = zt5
      bxy(3,j,k1) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
   10 continue
   20 continue
c mode number kx = 0
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k - 1)
      afdt = adt*aimag(ffc(1,k))
c update magnetic field half time step, ky > 0
      zt1 = cmplx(-aimag(exy(3,1,k)),real(exy(3,1,k)))
      zt3 = cmplx(-aimag(exy(1,1,k)),real(exy(1,1,k)))
      zt4 = bxy(1,1,k) - dth*(dky*zt1)
      zt6 = bxy(3,1,k) + dth*(dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,1,k) + cdt*(dky*zt1) - afdt*cu(1,1,k)
      zt9 = exy(3,1,k) - cdt*(dky*zt3) - afdt*cu(3,1,k)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exy(1,1,k) = zt7
      exy(2,1,k) = zero
      exy(3,1,k) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt9*conjg(zt9))
      zt4 = zt4 - dth*(dky*zt1)
      zt6 = zt6 + dth*(dky*zt3)
      bxy(1,1,k) = zt4
      bxy(2,1,k) = zero
      bxy(3,1,k) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt6*conjg(zt6))
c update magnetic field half time step, ky < 0
      bxy(2,1,k1) = zero
      exy(2,1,k1) = zero
      if (nyhm.eq.(ny/2)) then
         bxy(1,1,k1) = zero
         bxy(3,1,k1) = zero
         exy(1,1,k1) = zero
         exy(3,1,k1) = zero
      else
         zt1 = cmplx(-aimag(exy(3,1,k1)),real(exy(3,1,k1)))
         zt3 = cmplx(-aimag(exy(1,1,k1)),real(exy(1,1,k1)))
         zt4 = bxy(1,1,k1) + dth*(dky*zt1)
         zt6 = bxy(3,1,k1) - dth*(dky*zt3)
c update electric field whole time step
         zt1 = cmplx(-aimag(zt6),real(zt6))
         zt3 = cmplx(-aimag(zt4),real(zt4))
         zt7 = exy(1,1,k1) - cdt*(dky*zt1) - afdt*cu(1,1,k1)
         zt9 = exy(3,1,k1) + cdt*(dky*zt3) - afdt*cu(3,1,k1)
c update magnetic field half time step and store electric field
         zt1 = cmplx(-aimag(zt9),real(zt9))
         zt3 = cmplx(-aimag(zt7),real(zt7))
         exy(1,1,k1) = zt7
         exy(3,1,k1) = zt9
         zt4 = zt4 + dth*(dky*zt1)
         zt6 = zt6 - dth*(dky*zt3)
         bxy(1,1,k1) = zt4
         bxy(3,1,k1) = zt6
      endif
   30 continue
c mode number ky = 0
      k1 = nyhm + 1
      do 40 j = 2, nxhm
      dkx = dnx*real(j - 1)
      afdt = adt*aimag(ffc(j,1))
c update magnetic field half time step
      zt1 = cmplx(-aimag(exy(3,j,1)),real(exy(3,j,1)))
      zt2 = cmplx(-aimag(exy(2,j,1)),real(exy(2,j,1)))
      zt5 = bxy(2,j,1) + dth*(dkx*zt1)
      zt6 = bxy(3,j,1) - dth*(dkx*zt2)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt8 = exy(2,j,1) - cdt*(dkx*zt1) - afdt*cu(2,j,1)
      zt9 = exy(3,j,1) + cdt*(dkx*zt2) - afdt*cu(3,j,1)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      exy(1,j,1) = zero
      exy(2,j,1) = zt8
      exy(3,j,1) = zt9
      ws = ws + anorm*(zt8*conjg(zt8) + zt9*conjg(zt9))
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2)
      bxy(1,j,1) = zero
      bxy(2,j,1) = zt5
      bxy(3,j,1) = zt6
      wp = wp + anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
c unused extra row
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      exy(1,j,k1) = zero
      exy(2,j,k1) = zero
      exy(3,j,k1) = zero
   40 continue
c mode number kx = 0, ky = 0
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      exy(1,1,1) = zero
      exy(2,1,1) = zero
      exy(3,1,1) = zero
c unused extra element
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      wf = real(nx*ny)*ws
      wm = real(nx*ny)*c2*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine FIELD23GL(fxy,gxy,sctx,nx,ny,nxhm,nyhm,nxvh,nyv)
c for 2-1/2d code, this subroutine calculates real space 2d field gxy at
c integer grid co-ordinates from complex fourier co-efficients fxy
c input: all, output: gxy, sctx
c equations used are:
c gxy(1,x,y) = sum(fxy(1,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c                             exp(sqrt(-1)*2*m*pi*y/ny))
c gxy(2,x,y) = sum(fxy(2,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c                             exp(sqrt(-1)*2*m*pi*y/ny))
c gxy(3,x,y) = sum(fxy(3,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c                             exp(sqrt(-1)*2*m*pi*y/ny))
c fxy(1,n,m) = x component of force/charge at fourier mode (n-1,m-1)
c fxy(2,n,m) = y component of force/charge at fourier mode (n-1,m-1)
c fxy(3,n,m) = z component of force/charge at fourier mode (n-1,m-1)
c if nxhm=nx/2 and nyhm=ny/2, nx/2 and ny/2 modes are included, in a
c format compatible with the real to complex FFT used in traditional PIC
c gxy(1,j,k) = x component of force/charge at grid (j,k)
c gxy(2,j,k) = y component of force/charge at grid (j,k)
c gxy(3,j,k) = z component of force/charge at grid (j,k)
c sctx = scratch array for sines and cosines
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of field arrays, must be >= nxhm
c nyv = second dimension of field arrays, must be >= 2*nyhm
      implicit none
      integer nx, ny, nxhm, nyhm, nxvh, nyv
      real gxy
      complex fxy, sctx
      dimension gxy(3,2*nxvh,nyv)
      dimension fxy(3,nxvh,nyv), sctx(nxvh)
c local data
      integer j, k, k1, nxh, nyh, ny2, n, m
      real dnx, dny, dkx, dky, at1, at3
      double precision exl, eyl, ezl, ex, ey, ez
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = 2*nyhm + 2
      do 60 m = 1, ny
      do 50 n = 1, nx
c find electric field
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*real(n-1)
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k-1)*real(m-1)
      zt1 = cmplx(cos(dky),sin(dky))
      zt2 = conjg(zt1)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      do 20 j = 2, nxhm
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      exl = exl + (real(fxy(1,j,k)*zt3) + real(fxy(1,j,k1)*zt4))
      eyl = eyl + (real(fxy(2,j,k)*zt3) + real(fxy(2,j,k1)*zt4))
      ezl = ezl + (real(fxy(3,j,k)*zt3) + real(fxy(3,j,k1)*zt4))
   20 continue
c mode numbers kx = 0
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         zt2 = zt2*real(sctx(nxhm))
         exl = exl + (real(fxy(1,1,k)*zt1) + real(fxy(1,1,k1)*zt2))
         eyl = eyl + (real(fxy(2,1,k)*zt1) + real(fxy(2,1,k1)*zt2))
         ezl = ezl + (real(fxy(3,1,k)*zt1) + real(fxy(3,1,k1)*zt2))
      else
         exl = exl + real(fxy(1,1,k)*zt1)
         eyl = eyl + real(fxy(2,1,k)*zt1)
         ezl = ezl + real(fxy(3,1,k)*zt1)
      endif
      ex = ex + exl
      ey = ey + eyl
      ez = ez + ezl
   30 continue
c mode numbers ky = 0
      k1 = nyhm + 1
      dky = dny*real(nyh)*real(m-1)
      at1 = cos(dky)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      do 40 j = 2, nxhm
      zt3 = sctx(j-1)
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         zt1 = at1*zt3
         exl = exl + (real(fxy(1,j,1)*zt3) + real(fxy(1,j,k1)*zt1))
         eyl = eyl + (real(fxy(2,j,1)*zt3) + real(fxy(2,j,k1)*zt1))
         ezl = ezl + (real(fxy(3,j,1)*zt3) + real(fxy(3,j,k1)*zt1))
      else
         exl = exl + (real(fxy(1,j,1)*zt3))
         eyl = eyl + (real(fxy(2,j,1)*zt3))
         ezl = ezl + (real(fxy(3,j,1)*zt3))
      endif
   40 continue
      ex = 2.0d0*(ex + exl)
      ey = 2.0d0*(ey + eyl)
      ez = 2.0d0*(ez + ezl)
      at3 = real(sctx(nxhm))
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         ex = ex + (real(fxy(1,1,1)) + aimag(fxy(1,1,1))*at3)
         ey = ey + (real(fxy(2,1,1)) + aimag(fxy(2,1,1))*at3)
         ez = ez + (real(fxy(3,1,1)) + aimag(fxy(3,1,1))*at3)
      else
         ex = ex + real(fxy(1,1,1))
         ey = ey + real(fxy(2,1,1))
         ez = ez + real(fxy(3,1,1))
      endif
c special case to match conventional PIC code
      if (nxhm.eq.nyh) then
         if (nxhm.eq.nxh) then
            ex = ex + (real(fxy(1,1,k1)) + aimag(fxy(1,1,k1))*at3)*at1
            ey = ey + (real(fxy(2,1,k1)) + aimag(fxy(2,1,k1))*at3)*at1
            ez = ez + (real(fxy(3,1,k1)) + aimag(fxy(3,1,k1))*at3)*at1
         else
            ex = ex + real(fxy(1,1,k1))*at1
            ey = ey + real(fxy(2,1,k1))*at1
            ez = ez + real(fxy(3,1,k1))*at1
         endif
      endif
      gxy(1,n,m) = real(ex)
      gxy(2,n,m) = real(ey)
      gxy(3,n,m) = real(ez)
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine ELFIELD23GL(q,fxy,ffc,we,nx,ny,nxhm,nyhm,nxvh,nyv)
c this subroutine solves 2d poisson's equation in fourier space for
c unsmoothed longitudinal electric field with periodic boundary
c conditions, for gridless spectral code.
c Zeros out z component
c input: q,ffc,nx,ny,nxhm,nyhm,nxvh,nyv, output: fxy,we
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2),
c q(j,k) = complex charge density for fourier mode (j-1,k-1)
c fxy(1,j,k) = x component of complex electric field,
c fxy(2,j,k) = y component of complex electric field,
c fxy(3,j,k) = zero,
c all for fourier mode (j-1,k-1)
c zero row at k = nyhm + 1 is for compatibility with FFT solvers
c when nxhm = nx/2 and nyhm = ny/2
c aimag(ffc(j,k)) = finite-size particle shape factor s
c real(ffc(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c where affp = normalization constant = nx*ny/np, where np=number of
c particles
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of field arrays, must be >= nxhm
c nyv = second dimension of field arrays, must be >= 2*nyhm
      implicit none
      integer nx, ny, nxhm, nyhm, nxvh, nyv
      real we
      complex q, fxy, ffc
      dimension q(nxvh,nyv), fxy(3,nxvh,nyv)
      dimension ffc(nxhm,nyhm)
c local data
      integer j, k, k1, ny2
      real dnx, dny, dky, at1, at2, at3
      complex zt1, zt2, zero
      double precision wp, swp
      ny2 = 2*nyhm + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
c calculate unsmoothed longitudinal electric field and sum field energy
      swp = 0.0d0
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 20 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k - 1)
      wp = 0.0
      do 10 j = 2, nxhm
      at1 = real(ffc(j,k))
      at2 = dnx*real(j - 1)*at1
      at3 = dky*at1
      at1 = at1*aimag(ffc(j,k))
      zt1 = cmplx(aimag(q(j,k)),-real(q(j,k)))
      zt2 = cmplx(aimag(q(j,k1)),-real(q(j,k1)))
c ky > 0
      fxy(1,j,k) = at2*zt1
      fxy(2,j,k) = at3*zt1
      fxy(3,j,k) = zero
c ky < 0
      fxy(1,j,k1) = at2*zt2
      fxy(2,j,k1) = -at3*zt2
      fxy(3,j,k1) = zero
      wp = wp + at1*(q(j,k)*conjg(q(j,k)) + q(j,k1)*conjg(q(j,k1)))
   10 continue
      swp = swp + wp
   20 continue
c mode number kx = 0
      wp = 0.0
      do 30 k = 2, nyhm
      k1 = ny2 - k
      at1 = real(ffc(1,k))
      at3 = dny*real(k - 1)*at1
      at1 = at1*aimag(ffc(1,k))
c ky > 0
      zt1 = cmplx(aimag(q(1,k)),-real(q(1,k)))
      fxy(1,1,k) = zero
      fxy(2,1,k) = at3*zt1
      fxy(3,1,k) = zero
c ky < 0
      fxy(1,1,k1) = zero
      fxy(3,1,k1) = zero
      if (nyhm.eq.(ny/2)) then
         fxy(2,1,k1) = zero
      else
         zt1 = cmplx(aimag(q(1,k1)),-real(q(1,k1)))
         fxy(2,1,k1) = -at3*zt1
      endif
      wp = wp + at1*(q(1,k)*conjg(q(1,k)))
   30 continue
      swp = swp + wp
c mode number ky = 0
      k1 = nyhm + 1
      wp = 0.0
      do 40 j = 2, nxhm
      at1 = real(ffc(j,1))
      at2 = dnx*real(j - 1)*at1
      at1 = at1*aimag(ffc(j,1))
      zt1 = cmplx(aimag(q(j,1)),-real(q(j,1)))
      fxy(1,j,1) = at2*zt1
      fxy(2,j,1) = zero
      fxy(3,j,1) = zero
c unused extra row
      fxy(1,j,k1) = zero
      fxy(2,j,k1) = zero
      fxy(3,j,k1) = zero
      wp = wp + at1*(q(j,1)*conjg(q(j,1)))
   40 continue
      swp = swp + wp
c mode number kx = 0, ky = 0
      fxy(1,1,1) = zero
      fxy(2,1,1) = zero
      fxy(3,1,1) = zero
c unused extra element
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      fxy(3,1,k1) = zero
      we = real(nx*ny)*swp
      return
      end
c-----------------------------------------------------------------------
      subroutine PSMOOTH2GL(q,qs,ffc,ny,nxhm,nyhm,nxvh,nyv)
c this subroutine provides a scalar smoothing function
c with periodic boundary conditions, for gridless spectral code.
c input: q,ffc,nx,ny,nxhm,nyhm,nxvh,nyv, output: qs
c approximate flop count is: 4*nxhm*nyhm + 2*(nxhm + nyhm)
c where nxc = nxhm - 1, nyc = nyhm - 1
c smoothing is calculated using the equation:
c qs(kx,ky) = q(kx,ky)*s(kx,ky)
c where s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), and
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers
c q(j,k) = complex charge density
c aimag(ffc(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c zero row at k = nyhm + 1 is for compatibility with FFT solvers
c when nxhm = nx/2 and nyhm = ny/2
c ny = system length in y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of field arrays, must be >= nxhm
c nyv = second dimension of field arrays, must be >= 2*nyhm
      implicit none
      integer ny, nxhm, nyhm, nxvh, nyv
      complex q, qs, ffc
      dimension q(nxvh,nyv), qs(nxvh,nyv)
      dimension ffc(nxhm,nyhm)
c local data
      integer j, k, k1, ny2
      real at1
      complex zero
      ny2 = 2*nyhm + 2
      zero = cmplx(0.0,0.0)
c calculate smoothing
c mode numbers 0 <= kx < nxhm and 0 < ky < nyhm
      do 20 k = 2, nyhm
      k1 = ny2 - k
      do 10 j = 2, nxhm
      at1 = aimag(ffc(j,k))
      qs(j,k) = at1*q(j,k)
      qs(j,k1) = at1*q(j,k1)
   10 continue
   20 continue
c mode number kx = 0
      do 30 k = 2, nyhm
      k1 = ny2 - k
      at1 = aimag(ffc(1,k))
c ky > 0
      qs(1,k) = at1*q(1,k)
c ky < 0
      if (nyhm.eq.(ny/2)) then
         qs(1,k1) = zero
      else
         qs(1,k1) = at1*q(1,k1)
      endif
   30 continue
c mode number ky = 0
      k1 = nyhm + 1
      do 40 j = 1, nxhm
      at1 = aimag(ffc(j,1))
      qs(j,1) = at1*q(j,1)
c unused extra row
      qs(j,k1) = zero
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine POYNT2GL(q,exy,bxy,ffc,sx,sy,sz,nx,ny,nxhm,nyhm,nxvh,  
     1nyv)
c this subroutine calculates the total momentum in the electromagnetic
c field given by the poynting flux, for gridless spectral code.  inputs
c are the charge density, transverse electric field, and magnetic field.
c outputs are sx, sy, sz
c equation used is:
c sx = sum((fy(j,k)+exy(2,j,k))*conjg(bxy(3,j,k))-exy(3,j,k)*
c conjg(bxy(2,j,k)))
c sy = sum(exy(3,j,k)*conjg(bxy(1,j,k))-(fx(j,k)+exy(1,j,k))*
c conjg(bxy(3,j,k)))
c sz = sum((fx(j,k)+exy(1,j,k))*conjg(bxy(2,j,k))-(fy(j,k)+exy(2,j,k))*
c conjg(bxy(1,j,k))), where
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*q(kx,ky),
c kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky)
c q(j,k) = complex charge density for fourier mode (j-1,k-1)
c exy(1,j,k) = x component of transverse electric field,
c exy(2,j,k) = y component of transverse electric field,
c exy(3,j,k) = z component of transverse electric field,
c bxy(1,j,k) = x component of magnetic field,
c bxy(2,j,k) = y component of magnetic field,
c bxy(3,j,k) = z component of magnetic field,
c all for fourier mode (j-1,k-1)
c zero row at k = nyhm + 1 is for compatibility with FFT solvers
c when nxhm = nx/2 and nyhm = ny/2
c real(ffc(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c sx/sy/sz = x/y/z components of field momentum
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of field arrays, must be >= nxhm
c nyv = second dimension of field arrays, must be >= 2*nyhm
      implicit none
      integer nx, ny, nxhm, nyhm, nxvh, nyv
      real sx, sy, sz
      complex q, exy, bxy, ffc
      dimension q(nxvh,nyv), exy(3,nxvh,nyv), bxy(3,nxvh,nyv)
      dimension ffc(nxhm,nyhm)
c local data
      integer ny2, j, k, k1
      real dnx, dny, affp, dky, at1, at2, at3
      complex zt1, zt2
      double precision wx, wy, wz, wxl, wyl, wzl
      ny2 = 2*nyhm + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      affp = real(ffc(1,1))
c calculate force/charge and sum field energy
      wx = 0.0d0
      wy = 0.0d0
      wz = 0.0d0
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 20 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k - 1)
      wxl = 0.0d0
      wyl = 0.0d0
      wzl = 0.0d0
      do 10 j = 2, nxhm
      at1 = real(ffc(j,k))
      at2 = dnx*real(j - 1)*at1
      at3 = dky*at1
c ky > 0
      zt1 = cmplx(aimag(q(j,k)),-real(q(j,k)))
      zt2 = at3*zt1 + exy(2,j,k)
      zt1 = at2*zt1 + exy(1,j,k)
      wxl = wxl + (zt2*conjg(bxy(3,j,k))-exy(3,j,k)*conjg(bxy(2,j,k)))
      wyl = wyl + (exy(3,j,k)*conjg(bxy(1,j,k))-zt1*conjg(bxy(3,j,k)))
      wzl = wzl + (zt1*conjg(bxy(2,j,k))-zt2*conjg(bxy(1,j,k)))
c ky < 0
      zt2 = cmplx(aimag(q(j,k1)),-real(q(j,k1)))
      zt1 = at2*zt2 + exy(1,j,k1)
      zt2 = -at3*zt2 + exy(2,j,k1)
      wxl = wxl + (zt2*conjg(bxy(3,j,k1))-exy(3,j,k1)*conjg(bxy(2,j,k1))
     1)
      wyl = wyl + (exy(3,j,k1)*conjg(bxy(1,j,k1))-zt1*conjg(bxy(3,j,k1))
     1)
      wzl = wzl + (zt1*conjg(bxy(2,j,k1))-zt2*conjg(bxy(1,j,k1)))
   10 continue
      wx = wx + wxl
      wy = wy + wyl
      wz = wz + wzl
   20 continue
c mode numbers kx = 0
      wxl = 0.0d0
      wyl = 0.0d0
      wzl = 0.0d0
      do 30 k = 2, nyhm
      at1 = real(ffc(1,k))
      at3 = dny*real(k - 1)*at1
c ky > 0
      zt1 = cmplx(aimag(q(1,k)),-real(q(1,k)))
      zt2 = at3*zt1 + exy(2,1,k)
      zt1 = exy(1,1,k)
      wxl = wxl + (zt2*conjg(bxy(3,1,k))-exy(3,1,k)*conjg(bxy(2,1,k)))
      wyl = wyl + (exy(3,1,k)*conjg(bxy(1,1,k))-zt1*conjg(bxy(3,1,k)))
      wzl = wzl + (zt1*conjg(bxy(2,1,k))-zt2*conjg(bxy(1,1,k)))
c ky < 0
      if (nyhm.ne.(ny/2)) then
         zt1 = cmplx(aimag(q(1,k1)),-real(q(1,k1)))
         zt2 = -at3*zt1 + exy(2,1,k1)
         zt1 = exy(1,1,k1)
         wxl = wxl + (zt2*conjg(bxy(3,1,k1))                            
     1 - exy(3,1,k1)*conjg(bxy(2,1,k1)))
         wyl = wyl + (exy(3,1,k1)*conjg(bxy(1,1,k1))                    
     1 - zt1*conjg(bxy(3,1,k1)))
         wzl = wzl + (zt1*conjg(bxy(2,1,k1))-zt2*conjg(bxy(1,1,k1)))
      endif
   30 continue
      wx = wx + wxl
      wy = wy + wyl
      wz = wz + wzl
c mode numbers ky = 0
      wxl = 0.0d0
      wyl = 0.0d0
      wzl = 0.0d0
      do 40 j = 2, nxhm
      at1 = real(ffc(j,1))
      at2 = dnx*real(j - 1)*at1
      zt1 = cmplx(aimag(q(j,1)),-real(q(j,1)))
      zt2 = exy(2,j,1)
      zt1 = at2*zt1 + exy(1,j,1)
      wxl = wxl + (zt2*conjg(bxy(3,j,1))-exy(3,j,1)*conjg(bxy(2,j,1)))
      wyl = wyl + (exy(3,j,1)*conjg(bxy(1,j,1))-zt1*conjg(bxy(3,j,1)))
      wzl = wzl + (zt1*conjg(bxy(2,j,1))-zt2*conjg(bxy(1,j,1)))
   40 continue
      wx = wx + wxl
      wy = wy + wyl
      wz = wz + wzl
      at1 = 2.0*real(nx*ny)/affp
      sx = at1*wx
      sy = at1*wy
      sz = at1*wz
      return
      end
c-----------------------------------------------------------------------
      subroutine A0RBPUSH23GL(part,fxy,bxy,sctx,qbm,dt,dtc,ci,ek,idimp, 
     1nop,nx,ny,nxhm,nyhm,nxvh,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time using gridless spectral
c version, for relativistic particles with magnetic field, periodic
c boundaries and various particle boundary conditions.
c Using the Analytic Boris Mover,
c assumes constant E, B fields, and gamma during a time step
c gamma and energy calculated with Boris half push with electric field
c scalar version using guard cells
c 60*nxhm*nyhm + 122 flops/particle, 3*nxhm*nyhm+ 5 loads, 5 stores
c 5 divides, 3 sqrts, 1 tangent
c plus (nxhm + nyhm) sines and cosines/particle
c input: all, output: part, sctx, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + (q/m)*gx(x(t),y(t))) +
c    rot(2)*(py(t-dt/2) + (q/m)*gy(x(t),y(t))) +
c    rot(3)*(pz(t-dt/2) + (q/m)*gz(x(t),y(t)))) + (q/m)*gx(x(t),y(t)))
c py(t+dt/2) = rot(4)*(px(t-dt/2) + (q/m)*gx(x(t)),y(t)) +
c    rot(5)*(py(t-dt/2) + (q/m)*gy(x(t),y(t))) +
c    rot(6)*(pz(t-dt/2) + (q/m)*gz(x(t),y(t)))) + (q/m)*gy(x(t),y(t)))
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + (q/m)*gx(x(t),y(t))) +
c    rot(8)*(py(t-dt/2) + (q/m)*gy(x(t),y(t))t) +
c    rot(9)*(pz(t-dt/2) + (q/m)*gz(x(t),y(t)))) + (q/m)*gz(x(t),y(t)))
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
c gx(x,y) = 0.5*flx(x,y)*dt + (fx(x,y)-flx(x,y))*tan(0.5*om*dt)/om
c gy(x,y) = 0.5*fly(x,y)*dt + (fy(x,y)-fly(x,y))*tan(0.5*om*dt)/om
c gz(x,y) = 0.5*flz(x,y)*dt + (fz(x,y)-flz(x,y))*tan(0.5*om*dt)/om
c where flx(x,y) = fpl(x,y)*omx/om**2, fly(x,y) = fpl(x,y)*omy/om**2,
c flz(x,y) = fpl(x,y)*omz/om**2,
c and fpl(x,y) = fx(x,y)*omx+fy(x,y)*omy+fz(x,y)*omz
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
c omz = (q/m)*bz(x(t),y(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c fx(x(t),y(t)) and fy(x(t),y(t)) are calculated from the expression
c fx(x,y) = sum(fxy(1,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c part(5,n) = momentum pz of particle n
c fxy(i,n,m) = i component of force/charge at fourier grid point n,m
c bxy(i,n,m) = i component of magnetic field at fourier grid point n,m
c if nxhm=nx/2 and nyhm=ny/2, nx/2 and ny/2 modes are included, in a
c format compatible with the real to complex FFT used in traditional PIC
c sctx = scratch array for sines and cosines
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of field arrays, must be >= nxhm
c nyv = second dimension of field arrays, must be >= 2*nyhm
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxhm, nyhm, nxvh, nyv, ipbc
      real part, qbm, dt, dtc, ci, ek
      complex fxy, bxy, sctx
      dimension part(idimp,nop), fxy(3,nxvh,nyv), bxy(3,nxvh,nyv)
      dimension sctx(nxvh)
c local data
      integer i, j, k, k1, nxh, nyh, ny2
      real edgelx, edgely, edgerx, edgery, dnx, dny, dkx, dky, at1, at3
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9, anorm
      real dth, omti, epl, ci2, p2, gam, gami, qtmg, dtg, x, y
      double precision sum1, exl, eyl, ezl, ex, ey, ez
      double precision bxl, byl, bzl, bx, by, bz
      complex zt1, zt2, zt3, zt4
      dth = 0.5*dt
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = 2*nyhm + 2
      ci2 = ci*ci
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.0
         edgely = 0.0
         edgerx = real(nx)
         edgery = real(ny)
      else if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 0.0
         edgerx = real(nx-1)
         edgery = real(ny)
      endif
      do 50 i = 1, nop
      x = part(1,i)
      y = part(2,i)
c find electric and magnetic field
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*x
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
      bx = 0.0d0
      by = 0.0d0
      bz = 0.0d0
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k-1)*y
      zt1 = cmplx(cos(dky),sin(dky))
      zt2 = conjg(zt1)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      bxl = 0.0d0
      byl = 0.0d0
      bzl = 0.0d0
      do 20 j = 2, nxhm
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      exl = exl + (real(fxy(1,j,k)*zt3) + real(fxy(1,j,k1)*zt4))
      eyl = eyl + (real(fxy(2,j,k)*zt3) + real(fxy(2,j,k1)*zt4))
      ezl = ezl + (real(fxy(3,j,k)*zt3) + real(fxy(3,j,k1)*zt4))
      bxl = bxl + (real(bxy(1,j,k)*zt3) + real(bxy(1,j,k1)*zt4))
      byl = byl + (real(bxy(2,j,k)*zt3) + real(bxy(2,j,k1)*zt4))
      bzl = bzl + (real(bxy(3,j,k)*zt3) + real(bxy(3,j,k1)*zt4))
   20 continue
c mode number kx = 0
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         zt2 = zt2*real(sctx(nxhm))
         exl = exl + (real(fxy(1,1,k)*zt1) + real(fxy(1,1,k1)*zt2))
         eyl = eyl + (real(fxy(2,1,k)*zt1) + real(fxy(2,1,k1)*zt2))
         ezl = ezl + (real(fxy(3,1,k)*zt1) + real(fxy(3,1,k1)*zt2))
         bxl = bxl + (real(bxy(1,1,k)*zt1) + real(bxy(1,1,k1)*zt2))
         byl = byl + (real(bxy(2,1,k)*zt1) + real(bxy(2,1,k1)*zt2))
         bzl = bzl + (real(bxy(3,1,k)*zt1) + real(bxy(3,1,k1)*zt2))
      else
         exl = exl + real(fxy(1,1,k)*zt1)
         eyl = eyl + real(fxy(2,1,k)*zt1)
         ezl = ezl + real(fxy(3,1,k)*zt1)
         bxl = bxl + real(bxy(1,1,k)*zt1)
         byl = byl + real(bxy(2,1,k)*zt1)
         bzl = bzl + real(bxy(3,1,k)*zt1)
      endif
      ex = ex + exl
      ey = ey + eyl
      ez = ez + ezl
      bx = bx + bxl
      by = by + byl
      bz = bz + bzl
   30 continue
c mode number ky = 0
      k1 = nyhm + 1
      dky = dny*real(nyh)*y
      at1 = cos(dky)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      bxl = 0.0d0
      byl = 0.0d0
      bzl = 0.0d0
      do 40 j = 2, nxhm
      zt3 = sctx(j-1)
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         zt1 = at1*zt3
         exl = exl + (real(fxy(1,j,1)*zt3) + real(fxy(1,j,k1)*zt1))
         eyl = eyl + (real(fxy(2,j,1)*zt3) + real(fxy(2,j,k1)*zt1))
         ezl = ezl + (real(fxy(3,j,1)*zt3) + real(fxy(3,j,k1)*zt1))
         bxl = bxl + (real(bxy(1,j,1)*zt3) + real(bxy(1,j,k1)*zt1))
         byl = byl + (real(bxy(2,j,1)*zt3) + real(bxy(2,j,k1)*zt1))
         bzl = bzl + (real(bxy(3,j,1)*zt3) + real(bxy(3,j,k1)*zt1))
      else
         exl = exl + real(fxy(1,j,1)*zt3)
         eyl = eyl + real(fxy(2,j,1)*zt3)
         ezl = ezl + real(fxy(3,j,1)*zt3)
         bxl = bxl + real(bxy(1,j,1)*zt3)
         byl = byl + real(bxy(2,j,1)*zt3)
         bzl = bzl + real(bxy(3,j,1)*zt3)
      endif
   40 continue
      ex = 2.0d0*(ex + exl)
      ey = 2.0d0*(ey + eyl)
      ez = 2.0d0*(ez + ezl)
      bx = 2.0d0*(bx + bxl)
      by = 2.0d0*(by + byl)
      bz = 2.0d0*(bz + bzl)
      at3 = real(sctx(nxhm))
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         ex = ex + (real(fxy(1,1,1)) + aimag(fxy(1,1,1))*at3)
         ey = ey + (real(fxy(2,1,1)) + aimag(fxy(2,1,1))*at3)
         ez = ez + (real(fxy(3,1,1)) + aimag(fxy(3,1,1))*at3)
         bx = bx + (real(bxy(1,1,1)) + aimag(bxy(1,1,1))*at3)
         by = by + (real(bxy(2,1,1)) + aimag(bxy(2,1,1))*at3)
         bz = bz + (real(bxy(3,1,1)) + aimag(bxy(3,1,1))*at3)
      else
         ex = ex + real(fxy(1,1,1))
         ey = ey + real(fxy(2,1,1))
         ez = ez + real(fxy(3,1,1))
         bx = bx + real(bxy(1,1,1))
         by = by + real(bxy(2,1,1))
         bz = bz + real(bxy(3,1,1))
      endif
c special case to match conventional PIC code
      if (nxhm.eq.nyh) then
         if (nxhm.eq.nxh) then
            ex = ex + (real(fxy(1,1,k1)) + aimag(fxy(1,1,k1))*at3)*at1
            ey = ey + (real(fxy(2,1,k1)) + aimag(fxy(2,1,k1))*at3)*at1
            ez = ez + (real(fxy(3,1,k1)) + aimag(fxy(3,1,k1))*at3)*at1
            bx = bx + (real(bxy(1,1,k1)) + aimag(bxy(1,1,k1))*at3)*at1
            by = by + (real(bxy(2,1,k1)) + aimag(bxy(2,1,k1))*at3)*at1
            bz = bz + (real(bxy(3,1,k1)) + aimag(bxy(3,1,k1))*at3)*at1
         else
            ex = ex + real(fxy(1,1,k1))*at1
            ey = ey + real(fxy(2,1,k1))*at1
            ez = ez + real(fxy(3,1,k1))*at1
            bx = bx + real(bxy(1,1,k1))*at1
            by = by + real(bxy(2,1,k1))*at1
            bz = bz + real(bxy(3,1,k1))*at1
         endif
      endif
      dx = ex
      dy = ey
      dz = ez
      ox = bx
      oy = by
      oz = bz
c normalize electric field
      dx = qbm*dx
      dy = qbm*dy
      dz = qbm*dz
c half acceleration to calculate gamma
      acx = part(3,i)
      acy = part(4,i)
      acz = part(5,i)
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
      part(3,i) = dx
      part(4,i) = dy
      part(5,i) = dz
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = x + dx*dtg
      dy = y + dy*dtg
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i)
            part(3,i) = -part(3,i)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,i)
            part(4,i) = -part(4,i)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i)
            part(3,i) = -part(3,i)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,i) = dx
      part(2,i) = dy
   50 continue
c normalize kinetic energy
      ek = ek + sum1
      return
      end
