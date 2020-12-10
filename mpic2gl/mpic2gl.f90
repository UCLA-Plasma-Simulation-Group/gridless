!-----------------------------------------------------------------------
! Skeleton 2D Electrostatic Gridless OpenMP PIC code
! written by Viktor K. Decyk, UCLA
      program mpic2gl
      use mpush2gl_h
      use omplib_h
      implicit none
! indx/indy = exponent which determines grid points in x/y direction:
! nx = 2**indx, ny = 2**indy.
      integer, parameter :: indx =   7, indy =   7
! npx/npy = number of electrons distributed in x/y direction.
      integer, parameter :: npx =  768, npy =   768
! ndim = number of velocity coordinates = 2
      integer, parameter :: ndim = 2
! tend = time at end of simulation, in units of plasma frequency.
! dt = time interval between successive calculations.
! qme = charge on electron, in units of e.
      real, parameter :: tend = 10.0, dt = 0.1, qme = -1.0
! vtx/vty = thermal velocity of electrons in x/y direction
! vx0/vy0 = drift velocity of electrons in x/y direction.
      real, parameter :: vtx = 1.0, vty = 1.0, vx0 = 0.0, vy0 = 0.0
! ax/ay = smoothed particle size in x/y direction
      real :: ax = 0.040, ay = 0.040
! idimp = number of particle coordinates = 4
! ipbc = particle boundary condition: 1 = periodic
      integer :: idimp = 4, ipbc = 1
! wke/we/wt = particle kinetic/electric field/total energy
      real :: wke = 0.0, we = 0.0, wt = 0.0
! declare scalars for standard code
      integer :: np, nx, ny, nxh, nyh, nxe, nye, nxeh, nxhm, nyhm
      integer :: ntime, nloop, isign
      real :: qbme, affp
!
! declare scalars for OpenMP code
      integer :: irc, nvp
!
! declare arrays for standard code:
! part = particle array
      real, dimension(:,:), pointer :: part
! qe/qi = electron/ion charge density with guard cells
      complex, dimension(:,:), pointer :: qe, qi
! fxye = smoothed electric field with guard cells
      complex, dimension(:,:,:), pointer :: fxye
! ffc = form factor array for poisson solver
      complex, dimension(:,:), pointer :: ffc
! sctx = scratch array for sines and cosines
      complex, dimension(:), pointer :: sctx
!
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime
      real :: tdpost = 0.0, tfield = 0.0, tpush = 0.0
      double precision :: dtime
!
      irc = 0
! nvp = number of shared memory nodes (0=default)
      nvp = 0
!     write (*,*) 'enter number of nodes:'
!     read (5,*) nvp
! initialize for shared memory parallel processing
      call INIT_OMP(nvp)
!
! initialize scalars for standard code
! np = total number of particles in simulation
! nx/ny = number of grid points in x/y direction
      np = npx*npy; nx = 2**indx; ny = 2**indy
      nxh = nx/2; nyh = max(1,ny/2)
      nxe = nx + 2; nye = ny + 1
! nloop = number of time steps in simulation
! ntime = current time step
      nloop = tend/dt + .0001; ntime = 0
      qbme = qme
      affp = real(nx*ny)/real(np)
! nxhm/nyhm = number of fourier modes kept in x/y
      nxhm = nxh; nyhm = nyh
!     nxhm = nxh+1; nyhm = nyh+1
      nxeh = max(nxe/2,nxhm); nye = max(nye,2*nyhm)
!
! allocate data for standard code
      allocate(part(idimp,np))
      allocate(qe(nxeh,nye),qi(nxeh,nye),fxye(ndim,nxeh,nye))
      allocate(ffc(nxhm,nyhm),sctx(nxeh))
!
! calculate form factors for gridless code: updates ffc
      call FFC2INITGL(ffc,ax,ay,affp,nx,ny,nxhm,nyhm)
! calculate form factors for conventional PIC code: updates ffc
!     isign = 0
!     call POIS22GL(qe,fxye,isign,ffc,ax,ay,affp,we,nx,ny,nxhm,nyhm,nxeh&
!    &,nye)     
! initialize electrons
      call DISTR2(part,vtx,vty,vx0,vy0,npx,npy,idimp,np,nx,ny,ipbc)
! initialize ion background
      qi = 0.0
      call MDPOST2GL(part,qi,sctx,-qme,np,idimp,nx,ny,nxhm,nyhm,nxeh,nye&
     &)
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
      call MDPOST2GL(part,qe,sctx,qme,np,idimp,nx,ny,nxhm,nyhm,nxeh,nye)
      qe = qe + qi
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdpost = tdpost + time
!
! calculate force/charge in fourier space with standard procedure:
! updates fxye, we
      call dtimer(dtime,itime,-1)
      isign = -1
      call POIS22GL(qe,fxye,isign,ffc,ax,ay,affp,we,nx,ny,nxhm,nyhm,nxeh&
     &,nye)   
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! push particles with standard procedure: updates part, wke
      call dtimer(dtime,itime,-1)
      wke = 0.0
      call MPUSH2GL(part,fxye,qbme,dt,wke,idimp,np,nx,ny,nxhm,nyhm,nxeh,&
     &nye,ipbc)
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
      write (*,*) 'nxhm, nyhm = ', nxhm, nyhm
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
