!-----------------------------------------------------------------------
! Skeleton 1D Electrostatic Gridless PIC code
! written by Viktor K. Decyk, UCLA
      program pic1gl
      use push1gl_h
      implicit none
! indx = exponent which determines grid points in x direction:
! nx = 2**indx.
      integer, parameter :: indx =   9
! npx = number of electrons distributed in x direction.
      integer, parameter :: npx = 18432
! tend = time at end of simulation, in units of plasma frequency.
! dt = time interval between successive calculations.
! qme = charge on electron, in units of e.
      real, parameter :: tend = 10.0, dt = 0.1, qme = -1.0
! vtx = thermal velocity of electrons in x direction
! vx0 = drift velocity of electrons in x direction.
      real, parameter :: vtx = 1.0, vx0 = 0.0
! ax = smoothed particle size in x direction
      real :: ax = 0.4
! idimp = number of particle coordinates = 2
! ipbc = particle boundary condition: 1 = periodic
      integer :: idimp = 2, ipbc = 1
! wke/we/wt = particle kinetic/electric field/total energy
      real :: wke = 0.0, we = 0.0, wt = 0.0
! declare scalars for standard code
      integer :: np, nx, nxh, nxe, nxeh, nxhm
      integer :: ntime, nloop, isign
      real :: qbme, affp
!
! declare arrays for standard code:
! part = particle array
      real, dimension(:,:), pointer :: part
! qe/qi = electron/ion charge density with guard cells
      complex, dimension(:), pointer :: qe, qi
! fxe = smoothed electric field with guard cells
      complex, dimension(:), pointer :: fxe
! ffc = form factor array for poisson solver
      complex, dimension(:), pointer :: ffc
! sctx = scratch array for sines and cosines
      complex, dimension(:), pointer :: sctx
!
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime
      real :: tdpost = 0.0, tfield = 0.0, tpush = 0.0
      double precision :: dtime
!
! initialize scalars for standard code
! np = total number of particles in simulation
! nx = number of grid points in x direction
      np = npx; nx = 2**indx; nxh = nx/2
      nxe = nx + 2
! nloop = number of time steps in simulation
! ntime = current time step
      nloop = tend/dt + .0001; ntime = 0
      qbme = qme
      affp = real(nx)/real(np)
! nxhm = number of fourier modes kept
      nxhm = nxh; nxeh = max(nxe/2,nxhm)
!
! allocate data for standard code
      allocate(part(idimp,np))
      allocate(qe(nxeh),qi(nxeh),fxe(nxeh))
      allocate(ffc(nxhm),sctx(nxeh))
!
! calculate form factors for gridless code: updates ffc
      call FFC1INITGL(ffc,ax,affp,nx,nxhm)
! calculate form factors for conventional PIC code: updates ffc
!     isign = 0
!     call POIS1GL(qe,fxe,isign,ffc,ax,affp,we,nx,nxhm,nxeh)
! initialize electrons
      call DISTR1(part,vtx,vx0,npx,idimp,np,nx,ipbc)
! initialize ion background
      qi = 0.0
      call DPOST1GL(part,qi,sctx,-qme,np,idimp,nx,nxhm,nxeh)
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= ntime) go to 2000
!     write (*,*) 'ntime = ', ntime
!
! deposit charge and add ion background with standard procedure:
! updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      call DPOST1GL(part,qe,sctx,qme,np,idimp,nx,nxhm,nxeh)
      qe = qe + qi
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdpost = tdpost + time
!
! calculate force/charge in fourier space with standard procedure:
! updates fxe, we
      call dtimer(dtime,itime,-1)
      isign = -1
      call POIS1GL(qe,fxe,isign,ffc,ax,affp,we,nx,nxhm,nxeh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! push particles with standard procedure: updates part, wke
      wke = 0.0
      call dtimer(dtime,itime,-1)
      call PUSH1GL(part,fxe,sctx,qbme,dt,wke,idimp,np,nx,nxhm,nxeh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tpush = tpush + time
!
      if (ntime==0) then
         write (*,*) 'Initial Field, Kinetic and Total Energies:'
         write (*,'(3e14.7)') we, wke, wke + we
      endif
      ntime = ntime + 1
      go to 500
 2000 continue
!
! * * * end main iteration loop * * *
!
      write (*,*) 'ntime = ', ntime
      write (*,*) 'nxhm = ', nxhm
      write (*,*) 'Final Field, Kinetic and Total Energies:'
      write (*,'(3e14.7)') we, wke, wke + we
!
      write (*,*)
      write (*,*) 'deposit time = ', tdpost
      write (*,*) 'solver time = ', tfield
      write (*,*) 'push time = ', tpush
      time = tdpost + tpush
      write (*,*) 'total particle time = ', time
      wt = time + tfield
      write (*,*) 'total time = ', wt
      write (*,*)
!
      wt = 1.0e+09/(real(nloop)*real(np))
      write (*,*) 'Push Time (nsec) = ', tpush*wt
      write (*,*) 'Deposit Time (nsec) = ', tdpost*wt
      write (*,*) 'Total Particle Time (nsec) = ', time*wt
!
      stop
      end program
