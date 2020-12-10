c Fortran Library for Skeleton 2D Electrostatic Gridless OpenMP PIC Code
c written by Viktor K. Decyk, UCLA
c-----------------------------------------------------------------------
      subroutine DISTR2(part,vtx,vty,vdx,vdy,npx,npy,idimp,nop,nx,ny,   
     1ipbc)
c for 2d code, this subroutine calculates initial particle co-ordinates
c and velocities with uniform density and maxwellian velocity with drift
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c vtx/vty = thermal velocity of electrons in x/y direction
c vdx/vdy = drift velocity of beam electrons in x/y direction
c npx/npy = initial number of particles distributed in x/y direction
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer npx, npy, idimp, nop, nx, ny, ipbc
      real vtx, vty, vdx, vdy
      real part
      dimension part(idimp,nop)
c local data
      integer j, k, k1, npxy
      real edgelx, edgely, at1, at2, at3, sum1, sum2
      double precision dsum1, dsum2
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
   30 continue
c add correct drift
      dsum1 = 0.0d0
      dsum2 = 0.0d0
      do 40 j = 1, npxy
      dsum1 = dsum1 + part(3,j)
      dsum2 = dsum2 + part(4,j)
   40 continue
      sum1 = dsum1
      sum2 = dsum2
      at1 = 1.0/real(npxy)
      sum1 = at1*sum1 - vdx
      sum2 = at1*sum2 - vdy
      do 50 j = 1, npxy
      part(3,j) = part(3,j) - sum1
      part(4,j) = part(4,j) - sum2
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPUSH2GL(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxhm,nyhm, 
     1nxvh,nyv,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time using gridless spectral
c version, periodic boundaries and various particle boundary conditions.
c OpenMP version, parallelization over particles
c 28*nxhm*nyhm + 14 flops/particle, NX*NY loads, 4 stores
c plus (nxhm + nyhm) sines and cosines/particle
c input: all, output: part, sctx, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are calculated from the expression
c fx(x,y) = sum(fxy(1,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c fy(x,y) = sum(fxy(2,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c fxy(1,n,m) = x component of force/charge at fourier mode (n-1,m-1)
c fxy(2,n,m) = y component of force/charge at fourier mode (n-1,m-1)
c if nxhm=nx/2 and nyhm=ny/2, nx/2 and ny/2 modes are included, in a
c format compatible with the real to complex FFT used in traditional PIC
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of field arrays, must be >= nxhm
c nyv = second dimension of field arrays, must be >= 2*nyhm
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nop, idimp, nx, ny, nxhm, nyhm, nxvh, nyv, ipbc
      real part, qbm, dt, ek
      complex fxy
      dimension part(idimp,nop), fxy(2,nxvh,nyv)
c local data
      integer i, j, k, k1, nxh, nyh, ny2
      real edgelx, edgely, edgerx, edgery, dnx, dny, dkx, dky, at1, at3
      real qtm, x, y, dx, dy
      real zr3, zi3, zr4, zi4, frp, fip, frm, fim
      double precision sum1, exl, eyl, ex, ey
      complex zt1, zt2, zt3, zt4
      complex ssctx
      dimension ssctx(nxhm)
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = 2*nyhm + 2
      qtm = qbm*dt
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
c loop over particles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,k1,x,y,dkx,dky,at1,at3,zr3,zi3,zr4,zi4,frp,fip,frm,
!$OMP& fim,exl,eyl,ex,ey,dx,dy,zt1,zt2,zt3,zt4,ssctx) REDUCTION(+:sum1)
      do 50 i = 1, nop
      x = part(1,i)
      y = part(2,i)
c find electric field
c store sine/cosines in x
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*x
      ssctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ex = 0.0d0
      ey = 0.0d0
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k-1)*y
      zt1 = cmplx(cos(dky),sin(dky))
      zt2 = conjg(zt1)
      exl = 0.0d0
      eyl = 0.0d0
      do 20 j = 2, nxhm
      zt3 = ssctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      zr3 = real(zt3)
      zi3 = aimag(zt3)
      zr4 = real(zt4)
      zi4 = aimag(zt4)
!     exl = exl + (real(fxy(1,j,k)*zt3) + real(fxy(1,j,k1)*zt4))
      zt3 = fxy(1,j,k)
      zt4 = fxy(1,j,k1)
      frp = real(zt3)
      fip = aimag(zt3)
      frm = real(zt4)
      fim = aimag(zt4)
      exl = exl + (frp*zr3 - fip*zi3 + frm*zr4 - fim*zi4)
!     eyl = eyl + (real(fxy(2,j,k)*zt3) + real(fxy(2,j,k1)*zt4))
      zt3 = fxy(2,j,k)
      zt4 = fxy(2,j,k1)
      frp = real(zt3)
      fip = aimag(zt3)
      frm = real(zt4)
      fim = aimag(zt4)
      eyl = eyl + (frp*zr3 - fip*zi3 + frm*zr4 - fim*zi4)
   20 continue
c mode number kx = 0
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         zt2 = zt2*real(ssctx(nxhm))
         exl = exl + (real(fxy(1,1,k)*zt1) + real(fxy(1,1,k1)*zt2))
         eyl = eyl + (real(fxy(2,1,k)*zt1) + real(fxy(2,1,k1)*zt2))
      else
         exl = exl + real(fxy(1,1,k)*zt1)
         eyl = eyl + real(fxy(2,1,k)*zt1)
      endif
      ex = ex + exl
      ey = ey + eyl
   30 continue
c mode number ky = 0
      k1 = nyhm + 1
      dky = dny*real(nyh)*y
      at1 = cos(dky)
      exl = 0.0d0
      eyl = 0.0d0
      do 40 j = 2, nxhm
      zt3 = ssctx(j-1)
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         zt1 = at1*zt3
         exl = exl + (real(fxy(1,j,1)*zt3) + real(fxy(1,j,k1)*zt1))
         eyl = eyl + (real(fxy(2,j,1)*zt3) + real(fxy(2,j,k1)*zt1))
      else
         exl = exl + real(fxy(1,j,1)*zt3)
         eyl = eyl + real(fxy(2,j,1)*zt3)
      endif
   40 continue
      ex = 2.0d0*(ex + exl)
      ey = 2.0d0*(ey + eyl)
      at3 = real(ssctx(nxhm))
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         ex = ex + (real(fxy(1,1,1)) + aimag(fxy(1,1,1))*at3)
         ey = ey + (real(fxy(2,1,1)) + aimag(fxy(2,1,1))*at3)
      else
         ex = ex + real(fxy(1,1,1))
         ey = ey + real(fxy(2,1,1))
      endif
c special case to match conventional PIC code
      if (nxhm.eq.nyh) then
         if (nxhm.eq.nxh) then
            ex = ex + (real(fxy(1,1,k1)) + aimag(fxy(1,1,k1))*at3)*at1
            ey = ey + (real(fxy(2,1,k1)) + aimag(fxy(2,1,k1))*at3)*at1
         else
            ex = ex + real(fxy(1,1,k1))*at1
            ey = ey + real(fxy(2,1,k1))*at1
         endif
      endif
      dx = ex
      dy = ey
c new velocity
      dx = part(3,i) + qtm*dx
      dy = part(4,i) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,i))**2 + (dy + part(4,i))**2
      part(3,i) = dx
      part(4,i) = dy
c new position
      dx = x + dx*dt
      dy = y + dy*dt
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
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine M1PUSH2GL(part,fxy,sctx,qbm,dt,ek,idimp,nop,nx,ny,nxhm,
     1nyhm,nxvh,nyv,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time using gridless spectral
c version, periodic boundaries and various particle boundary conditions.
c OpenMP version, parallelization over fourier modes
c 28*nxhm*nyhm + 14 flops/particle, NX*NY loads, 4 stores
c plus (nxhm + nyhm) sines and cosines/particle
c input: all, output: part, sctx, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are calculated from the expression
c fx(x,y) = sum(fxy(1,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c fy(x,y) = sum(fxy(2,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c fxy(1,n,m) = x component of force/charge at fourier mode (n-1,m-1)
c fxy(2,n,m) = y component of force/charge at fourier mode (n-1,m-1)
c if nxhm=nx/2 and nyhm=ny/2, nx/2 and ny/2 modes are included, in a
c format compatible with the real to complex FFT used in traditional PIC
c sctx = scratch array for sines and cosines
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxhm/nyhm = number of fourier modes kept in x/y
c nxvh = first dimension of field arrays, must be >= nxhm
c nyv = second dimension of field arrays, must be >= 2*nyhm
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nop, idimp, nx, ny, nxhm, nyhm, nxvh, nyv, ipbc
      real part, qbm, dt, ek
      complex fxy, sctx
      dimension part(idimp,nop), fxy(2,nxvh,nyv), sctx(nxvh)
c local data
      integer i, j, k, k1, nxh, nyh, ny2
      real edgelx, edgely, edgerx, edgery, dnx, dny, dkx, dky, at1, at3
      real qtm, x, y, dx, dy
      real zr3, zi3, zr4, zi4, frp, fip, frm, fim
      double precision sum1, exl, eyl, ex, ey
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = 2*nyhm + 2
      qtm = qbm*dt
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
c loop over particles
      do 50 i = 1, nop
      x = part(1,i)
      y = part(2,i)
c find electric field
c store sine/cosines in x
c !$OMP PARALLEL DO PRIVATE(j,dkx)
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*x
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
c !$OMP END PARALLEL DO
      ex = 0.0d0
      ey = 0.0d0
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,k1,dky,zr3,zi3,zr4,zi4,frp,fip,frm,fim,exl,eyl,zt1,  
!$OMP& zt2,zt3,zt4) REDUCTION(+:ex,ey)
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k-1)*y
      zt1 = cmplx(cos(dky),sin(dky))
      zt2 = conjg(zt1)
      exl = 0.0d0
      eyl = 0.0d0
      do 20 j = 2, nxhm
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      zr3 = real(zt3)
      zi3 = aimag(zt3)
      zr4 = real(zt4)
      zi4 = aimag(zt4)
!     exl = exl + (real(fxy(1,j,k)*zt3) + real(fxy(1,j,k1)*zt4))
      zt3 = fxy(1,j,k)
      zt4 = fxy(1,j,k1)
      frp = real(zt3)
      fip = aimag(zt3)
      frm = real(zt4)
      fim = aimag(zt4)
      exl = exl + (frp*zr3 - fip*zi3 + frm*zr4 - fim*zi4)
!     eyl = eyl + (real(fxy(2,j,k)*zt3) + real(fxy(2,j,k1)*zt4))
      zt3 = fxy(2,j,k)
      zt4 = fxy(2,j,k1)
      frp = real(zt3)
      fip = aimag(zt3)
      frm = real(zt4)
      fim = aimag(zt4)
      eyl = eyl + (frp*zr3 - fip*zi3 + frm*zr4 - fim*zi4)
   20 continue
c mode number kx = 0
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         zt2 = zt2*real(sctx(nxhm))
         exl = exl + (real(fxy(1,1,k)*zt1) + real(fxy(1,1,k1)*zt2))
         eyl = eyl + (real(fxy(2,1,k)*zt1) + real(fxy(2,1,k1)*zt2))
      else
         exl = exl + real(fxy(1,1,k)*zt1)
         eyl = eyl + real(fxy(2,1,k)*zt1)
      endif
      ex = ex + exl
      ey = ey + eyl
   30 continue
!$OMP END PARALLEL DO
c mode number ky = 0
      k1 = nyhm + 1
      dky = dny*real(nyh)*y
      at1 = cos(dky)
      exl = 0.0d0
      eyl = 0.0d0
c !$OMP PARALLEL DO PRIVATE(j,zt1,zt3) REDUCTION(+:exl,eyl)
      do 40 j = 2, nxhm
      zt3 = sctx(j-1)
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         zt1 = at1*zt3
         exl = exl + (real(fxy(1,j,1)*zt3) + real(fxy(1,j,k1)*zt1))
         eyl = eyl + (real(fxy(2,j,1)*zt3) + real(fxy(2,j,k1)*zt1))
      else
         exl = exl + real(fxy(1,j,1)*zt3)
         eyl = eyl + real(fxy(2,j,1)*zt3)
      endif
   40 continue
c !$OMP END PARALLEL DO
      ex = 2.0d0*(ex + exl)
      ey = 2.0d0*(ey + eyl)
      at3 = real(sctx(nxhm))
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         ex = ex + (real(fxy(1,1,1)) + aimag(fxy(1,1,1))*at3)
         ey = ey + (real(fxy(2,1,1)) + aimag(fxy(2,1,1))*at3)
      else
         ex = ex + real(fxy(1,1,1))
         ey = ey + real(fxy(2,1,1))
      endif
c special case to match conventional PIC code
      if (nxhm.eq.nyh) then
         if (nxhm.eq.nxh) then
            ex = ex + (real(fxy(1,1,k1)) + aimag(fxy(1,1,k1))*at3)*at1
            ey = ey + (real(fxy(2,1,k1)) + aimag(fxy(2,1,k1))*at3)*at1
         else
            ex = ex + real(fxy(1,1,k1))*at1
            ey = ey + real(fxy(2,1,k1))*at1
         endif
      endif
      dx = ex
      dy = ey
c new velocity
      dx = part(3,i) + qtm*dx
      dy = part(4,i) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,i))**2 + (dy + part(4,i))**2
      part(3,i) = dx
      part(4,i) = dy
c new position
      dx = x + dx*dt
      dy = y + dy*dt
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
      ek = ek + .125*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine MDPOST2GL(part,q,sctx,qm,nop,idimp,nx,ny,nxhm,nyhm,nxvh
     1,nyv)
c for 2d code, this subroutine calculates particle charge density
c using gridless spectral version, periodic boundaries
c OpenMP version, parallelization over fourier modes
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
c !$OMP PARALLEL DO PRIVATE(j,dkx)
      do 10 j = 1, nxhm
      dkx = dnx*real(j)*x
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c !$OMP END PARALLEL DO
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
!$OMP PARALLEL DO PRIVATE(j,k,k1,dky,zt1,zt2,zt3)
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
!$OMP END PARALLEL DO
c mode number ky = 0
      k1 = nyhm + 1
      dky = dny*real(nyh)*y
      at1 = qmn*cos(dky)
c !$OMP PARALLEL DO PRIVATE(j,zt3)
      do 40 j = 2, nxhm
      zt3 = sctx(j-1)
      q(j,1) = q(j,1) + qmn*zt3
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         q(j,k1) = q(j,k1) + at1*zt3
      endif
   40 continue
c !$OMP END PARALLEL DO
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
      subroutine POIS22GL(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxhm,nyhm,
     1nxvh,nyv)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions, for gridless spectral code.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxhm,nyhm,nxvh,nyv
c output: ffc
c for isign /= 0, input: q,ffc,isign,nx,ny,nxhm,nyhm,nxvh,nyv
c output: fxy,we
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2),
c q(j,k) = complex charge density for fourier mode (j-1,k-1)
c fxy(1,j,k) = x component of complex force/charge,
c fxy(2,j,k) = y component of complex force/charge,
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
      dimension q(nxvh,nyv), fxy(2,nxvh,nyv)
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
c ky < 0
      fxy(1,j,k1) = at2*zt2
      fxy(2,j,k1) = -at3*zt2
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
c ky < 0
      fxy(1,1,k1) = zero
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
c unused extra row
      fxy(1,j,k1) = zero
      fxy(2,j,k1) = zero
      wp = wp + at1*(q(j,1)*conjg(q(j,1)))
   70 continue
      swp = swp + wp
c mode number kx = 0, ky = 0
      fxy(1,1,1) = zero
      fxy(2,1,1) = zero
c unused extra element
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      we = real(nx*ny)*swp
      return
      end
c-----------------------------------------------------------------------
      subroutine MEVFIELD22GL(fxy,gxy,sctx,xp,yp,nx,ny,nxhm,nyhm,nxvh,  
     1nyv)
c for 2d code, this subroutine calculates real space 2d field gxy at
c location xp, yp, from complex fourier co-efficients fxy
c OpenMP version, parallelization over fourier modes
c input: all, output: gxy, sctx
c equations used are:
c gxy(1) = sum(fxy(1,n,m)*exp(sqrt(-1)*2*n*pi*xp/nx)*
c                             exp(sqrt(-1)*2*m*pi*yp/ny))
c gxy(2) = sum(fxy(2,n,m)*exp(sqrt(-1)*2*n*pi*xp/nx)*
c                             exp(sqrt(-1)*2*m*pi*yp/ny))
c fxy(1,n,m) = x component of field at fourier mode (n-1,m-1)
c fxy(2,n,m) = y component of field at fourier mode (n-1,m-1)
c if nxhm=nx/2 and nyhm=ny/2, nx/2 and ny/2 modes are included, in a
c format compatible with the real to complex FFT used in traditional PIC
c gxy(1) = x component of field at location (xp,yp)
c gxy(2) = y component of field at location (xp,yp)
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
      dimension gxy(2)
      dimension fxy(2,nxvh,nyv), sctx(nxvh)
c local data
      integer j, k, k1, nxh, nyh, ny2
      real dnx, dny, dkx, dky, at1, at3
      double precision exl, eyl, ex, ey
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
c mode numbers 0 < kx < nxhm and 0 < ky < nyhm
!$OMP PARALLEL DO PRIVATE(j,k,k1,dky,exl,eyl,zt1,zt2,zt3,zt4)
!$OMP& REDUCTION(+:ex,ey)
      do 30 k = 2, nyhm
      k1 = ny2 - k
      dky = dny*real(k-1)*yp
      zt1 = cmplx(cos(dky),sin(dky))
      zt2 = conjg(zt1)
      exl = 0.0d0
      eyl = 0.0d0
      do 20 j = 2, nxhm
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      exl = exl + (real(fxy(1,j,k)*zt3) + real(fxy(1,j,k1)*zt4))
      eyl = eyl + (real(fxy(2,j,k)*zt3) + real(fxy(2,j,k1)*zt4))
   20 continue
c mode numbers kx = 0
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         zt2 = zt2*real(sctx(nxhm))
         exl = exl + (real(fxy(1,1,k)*zt1) + real(fxy(1,1,k1)*zt2))
         eyl = eyl + (real(fxy(2,1,k)*zt1) + real(fxy(2,1,k1)*zt2))
      else
         exl = exl + real(fxy(1,1,k)*zt1)
         eyl = eyl + real(fxy(2,1,k)*zt1)
      endif
      ex = ex + exl
      ey = ey + eyl
   30 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0
      k1 = nyhm + 1
      dky = dny*real(nyh)*yp
      at1 = cos(dky)
      exl = 0.0d0
      eyl = 0.0d0
      do 40 j = 2, nxhm
      zt3 = sctx(j-1)
c special case to match conventional PIC code
      if (nyhm.eq.nyh) then
         zt1 = at1*zt3
         exl = exl + (real(fxy(1,j,1)*zt3) + real(fxy(1,j,k1)*zt1))
         eyl = eyl + (real(fxy(2,j,1)*zt3) + real(fxy(2,j,k1)*zt1))
      else
         exl = exl + (real(fxy(1,j,1)*zt3))
         eyl = eyl + (real(fxy(2,j,1)*zt3))
      endif
   40 continue
      ex = 2.0d0*(ex + exl)
      ey = 2.0d0*(ey + eyl)
      at3 = real(sctx(nxhm))
c special case to match conventional PIC code
      if (nxhm.eq.nxh) then
         ex = ex + (real(fxy(1,1,1)) + aimag(fxy(1,1,1))*at3)
         ey = ey + (real(fxy(2,1,1)) + aimag(fxy(2,1,1))*at3)
      else
         ex = ex + real(fxy(1,1,1))
         ey = ey + real(fxy(2,1,1))
      endif
c special case to match conventional PIC code
      if (nxhm.eq.nyh) then
         if (nxhm.eq.nxh) then
            ex = ex + (real(fxy(1,1,k1)) + aimag(fxy(1,1,k1))*at3)*at1
            ey = ey + (real(fxy(2,1,k1)) + aimag(fxy(2,1,k1))*at3)*at1
         else
            ex = ex + real(fxy(1,1,k1))*at1
            ey = ey + real(fxy(2,1,k1))*at1
         endif
      endif
      gxy(1) = real(ex)
      gxy(2) = real(ey)
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
