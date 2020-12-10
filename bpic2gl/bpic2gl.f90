!-----------------------------------------------------------------------
! Skeleton 2-1/2D Electromagnetic Gridless PIC code
! written by Viktor K. Decyk, UCLA
      program bpic2gl
      use bpush2gl_h
      implicit none
! indx/indy = exponent which determines grid points in x/y direction:
! nx = 2**indx, ny = 2**indy.
      integer, parameter :: indx =   7, indy =   7
! npx/npy = number of electrons distributed in x/y direction.
      integer, parameter :: npx =  768, npy =   768
! ndim = number of velocity coordinates = 3
      integer, parameter :: ndim = 3
! tend = time at end of simulation, in units of plasma frequency.
! dt = time interval between successive calculations.
! qme = charge on electron, in units of e.
      real, parameter :: tend = 10.0, dt = 0.04, qme = -1.0
! vtx/vty = thermal velocity of electrons in x/y direction
! vx0/vy0 = drift velocity of electrons in x/y direction.
      real, parameter :: vtx = 1.0, vty = 1.0, vx0 = 0.0, vy0 = 0.0
! vtx/vz0 = thermal/drift velocity of electrons in z direction
      real, parameter :: vtz = 1.0, vz0 = 0.0
! ax/ay = smoothed particle size in x/y direction
! ci = reciprocal of velocity of light.
      real :: ax = 0.4, ay = 0.4, ci = 0.1
! idimp = number of particle coordinates = 5
! ipbc = particle boundary condition: 1 = periodic
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: idimp = 5, ipbc = 1, relativity = 1
! wke/we = particle kinetic/electrostatic field energy
! wf/wm/wt = magnetic field/transverse electric field/total energy
      real :: wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0
! declare scalars for standard code
      integer :: np, nx, ny, nxh, nyh, nxe, nye, nxeh, nxhm, nyhm
      integer :: ntime, nloop, isign
      real :: qbme, affp, at, dth
!
! declare arrays for standard code:
! part = particle array
      real, dimension(:,:), pointer :: part
! qe/qi = electron/ion charge density
      complex, dimension(:,:), pointer :: qe, qi
! cue = electron current density
! fxyze/bxyze = smoothed electric/magnetic field
      complex, dimension(:,:,:), pointer :: cue, fxyze, bxyze
! exyz/bxyz = transverse electric/magnetic field in fourier space
! dcu = time derivative of electron current density
      complex, dimension(:,:,:), pointer :: exyz, bxyz, dcu
! ffc = form factor array for poisson solver
      complex, dimension(:,:), pointer :: ffc
! sctx = scratch array for sines and cosines
      complex, dimension(:), pointer :: sctx
!
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime
      real :: tdpost = 0.0, tfield = 0.0, tdjpost = 0.0, tpush = 0.0
      double precision :: dtime
!
! initialize scalars for standard code
! np = total number of particles in simulation
! nx/ny = number of grid points in x/y direction
      np = npx*npy; nx = 2**indx; ny = 2**indy
      nxh = nx/2; nyh = max(1,ny/2)
      nxe = nx + 2; nye = ny + 1; nxeh = nxe/2
! nloop = number of time steps in simulation
! ntime = current time step
      nloop = tend/dt + .0001; ntime = 0
      qbme = qme
      affp = real(nx*ny)/real(np)
! nxhm/nyhm = number of fourier modes kept in x/y
      nxhm = nxh; nyhm = nyh
!     nxhm = nxh+1; nyhm = nyh+1
      nxeh = max(nxe/2,nxhm); nye = max(nye,2*nyhm)
      dth = 0.0
!
! allocate data for standard code
      allocate(part(idimp,np))
      allocate(qe(nxeh,nye),qi(nxeh,nye),fxyze(ndim,nxeh,nye))
      allocate(cue(ndim,nxeh,nye),bxyze(ndim,nxeh,nye))
      allocate(exyz(ndim,nxeh,nye),bxyz(ndim,nxeh,nye))
      allocate(ffc(nxhm,nyhm),sctx(nxeh))
!
! calculate form factors for gridless code: updates ffc
      call FFC2INITGL(ffc,ax,ay,affp,nx,ny,nxhm,nyhm)
! calculate form factors for conventional PIC code: updates ffc
!     isign = 0
!     call POIS23GL(qe,fxyze,isign,ffc,ax,ay,affp,we,nx,ny,nxhm,nyhm,   &
!    &nxeh,nye)
! initialize electrons
      call DISTR2H(part,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,idimp,np,nx,ny, &
     &ipbc)
! initialize ion background
      qi = 0.0
      call DPOST2GL(part,qi,sctx,-qme,np,idimp,nx,ny,nxhm,nyhm,nxeh,nye)
!
! initialize transverse electromagnetic fields
      exyz = cmplx(0.0,0.0)
      bxyz = cmplx(0.0,0.0)
!
      at = sqrt((real(nxhm)/real(nxh))**2+(real(nyhm)/real(nyh))**2)
      if (dt > (0.63*ci/at)) then
         write (*,*) 'Warning: Courant condition may be exceeded!'
      endif
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= ntime) go to 2000
!     write (*,*) 'ntime = ', ntime
!
! deposit current with standard procedure: updates part, cue
      call dtimer(dtime,itime,-1)
      cue = cmplx(0.0,0.0)
      if (relativity==1) then
         call RDJPOST2GL(part,cue,sctx,qme,dth,ci,np,idimp,nx,ny,nxhm,  &
     &nyhm,nxeh,nye,ipbc)
      else
         call DJPOST2GL(part,cue,sctx,qme,dth,np,idimp,nx,ny,nxhm,nyhm, &
     &nxeh,nye,ipbc)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdjpost = tdjpost + time
!
! deposit charge with standard procedure: updates qe
      call dtimer(dtime,itime,-1)
      qe = cmplx(0.0,0.0)
      call DPOST2GL(part,qe,sctx,qme,np,idimp,nx,ny,nxhm,nyhm,nxeh,nye)
      qe = qe + qi
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdpost = tdpost + time
!
! take transverse part of current with standard procedure: updates cue
      call dtimer(dtime,itime,-1)
      call CUPERP2GL(cue,nx,ny,nxhm,nyhm,nxeh,nye)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate electromagnetic fields in fourier space with standard
! procedure: updates exyz, bxyz, wf, wm
      call dtimer(dtime,itime,-1)
! initialization assumes initial acceleration of particles is zero
      if (ntime==0) then
         allocate(dcu(ndim,nxeh,nye))
         dcu = cmplx(0.0,0.0)
! initialize with darwin electric and magnetic fields
! deposit time derivative of current: updates dcu
         if (relativity==1) then
            call RD2JPOST2GL(part,dcu,sctx,qme,ci,np,idimp,nx,ny,nxhm,  &
     &nyhm,nxeh,nye)
         else
            call D2JPOST2GL(part,dcu,sctx,qme,np,idimp,nx,ny,nxhm,nyhm, &
     &nxeh,nye)
         endif
! initialize electromagnetic fields from free-streaming particles
! deposit scaled current and time derivative of current:
! updates cue, dcu
!        cue = cmplx(0.0,0.0)
!        if (relativity==1) then
!           call RDSJPOST2GL(part,cue,dcu,sctx,qme,ci,np,idimp,nx,ny,   &
!    &nxhm,nyhm,nxeh,nye)
!        else
!           call DSJPOST2GL(part,cue,dcu,sctx,qme,ci,np,idimp,nx,ny,nxhm&
!    &,nyhm,nxeh,nye)
!        endif
!        call CUPERP2GL(cue,nx,ny,nxhm,nyhm,nxeh,nye)
! initialize magnetic field: updates bxyz, wm
         call IBPOIS23GL(cue,bxyz,ffc,ci,wm,nx,ny,nxhm,nyhm,nxeh,nye)
         wf = 0.0
! calculates transverse part of the derivative of the current
! updates: dcu
         call CUPERP2GL(dcu,nx,ny,nxhm,nyhm,nxeh,nye)
! initialize transverse electric field: updates exyz
         isign = 1
         call EPOIS23GL(dcu,exyz,isign,ffc,ax,ay,affp,0.0,ci,wf,nx,ny,  &
     &nxhm,nyhm,nxeh,nye)
         deallocate(dcu)
         dth = 0.5*dt
! calculates transverse electromagnetic fields:
! updates exyz, bxyz, wf, wm
      else
         call AMAXWEL2GL(exyz,bxyz,cue,ffc,ci,dt,wf,wm,nx,ny,nxhm,nyhm, &
     &nxeh,nye)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate force/charge in fourier space with standard procedure:
! updates fxyze, we
      call dtimer(dtime,itime,-1)
      isign = -1
      call POIS23GL(qe,fxyze,isign,ffc,ax,ay,affp,we,nx,ny,nxhm,nyhm,   &
     &nxeh,nye)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! add longitudinal and transverse electric fields with standard
! procedure: updates fxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call EMFIELD2GL(fxyze,exyz,ffc,isign,nxhm,nyhm,nxeh,nye)
! copy magnetic field with standard procedure: updates bxyze
      isign = -1
      call EMFIELD2GL(bxyze,bxyz,ffc,isign,nxhm,nyhm,nxeh,nye)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! push particles with standard procedure: updates part, wke
      wke = 0.0
      call dtimer(dtime,itime,-1)
      if (relativity==1) then
! analytic Boris mover, gamma constant during time step
         call ARBPUSH23GL(part,fxyze,bxyze,sctx,qbme,dt,dth,ci,wke,idimp&
     &,np,nx,ny,nxhm,nyhm,nxeh,nye,ipbc)
! exact analytic Boris mover, gamma varies during time step
!        call EARBPUSH23GL(part,fxyze,bxyze,sctx,qbme,dt,dth,ci,wke,    &
!    &idimp,np,nx,ny,nxhm,nyhm,nxeh,nye,ipbc)
      else
         call ABPUSH23GL(part,fxyze,bxyze,sctx,qbme,dt,dth,wke,idimp,np,&
     &nx,ny,nxhm,nyhm,nxeh,nye,ipbc)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tpush = tpush + time
!
      if (ntime==0) then
         wt = we + wf + wm
         write (*,*) 'Initial Total Field, Kinetic and Total Energies:'
         write (*,'(3e14.7)') wt, wke, wke + wt
         write (*,*) 'Initial Electrostatic, Transverse Electric and Mag&
     &netic Field Energies:'
         write (*,'(3e14.7)') we, wf, wm
      endif
      ntime = ntime + 1
      go to 500
 2000 continue
!
! * * * end main iteration loop * * *
!
      write (*,*) 'ntime, relativity = ', ntime, relativity
      write (*,*) 'nxhm, nyhm = ', nxhm, nyhm
      wt = we + wf + wm
      write (*,*) 'Final Total Field, Kinetic and Total Energies:'
      write (*,'(3e14.7)') wt, wke, wke + wt
      write (*,*) 'Final Electrostatic, Transverse Electric and Magnetic&
     & Field Energies:'
      write (*,'(3e14.7)') we, wf, wm
!
      write (*,*)
      write (*,*) 'deposit time = ', tdpost
      write (*,*) 'current deposit time = ', tdjpost
      tdpost = tdpost + tdjpost
      write (*,*) 'total deposit time = ', tdpost
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
