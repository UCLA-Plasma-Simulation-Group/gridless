!-----------------------------------------------------------------------
! Skeleton 1-2/2D Electromagnetic Gridless OpenMP PIC code
! written by Viktor K. Decyk, UCLA
      program mbpic1gl
      use mbpush1gl_h
      use omplib_h
      implicit none
! indx = exponent which determines grid points in x direction:
! nx = 2**indx.
      integer, parameter :: indx =   9
! npx = number of electrons distributed in x direction.
      integer, parameter :: npx =  18432
! tend = time at end of simulation, in units of plasma frequency.
! dt = time interval between successive calculations.
! qme = charge on electron, in units of e.
      real, parameter :: tend = 10.0, dt = 0.05, qme = -1.0
! vtx/vty = thermal velocity of electrons in x/y direction
! vx0/vy0 = drift velocity of electrons in x/y direction.
      real, parameter :: vtx = 1.0, vty = 1.0, vx0 = 0.0, vy0 = 0.0
! vtx/vz0 = thermal/drift velocity of electrons in z direction
      real, parameter :: vtz = 1.0, vz0 = 0.0
! omx = magnetic field electron cyclotron frequency in x
      real :: omx = 0.0
! ax = smoothed particle size in x direction
! ci = reciprocal of velocity of light.
      real :: ax = 0.4, ci = 0.1
! idimp = number of particle coordinates = 4
! ipbc = particle boundary condition: 1 = periodic
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: idimp = 4, ipbc = 1, relativity = 1
! wke/we = particle kinetic/electrostatic field energy
! wf/wm/wt = magnetic field/transverse electric field/total energy
      real :: wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0
! declare scalars for standard code
      integer :: np, nx, nxh, nxe, nxeh, nxhm
      integer :: ntime, nloop, isign
      real :: qbme, affp, dth, cux0
!
! declare scalars for OpenMP code
      integer :: irc, nvp
!
! declare arrays for standard code:
! part = particle array
      real, dimension(:,:), pointer :: part
! qe/qi = electron charge density with guard cells
! fxe = smoothed longitudinal electric field with guard cells
      complex, dimension(:), pointer :: qe, qi, fxe
! cue = electron current density with guard cells
! fxyze/byze = smoothed electric/magnetic field with guard cells
      complex, dimension(:,:), pointer :: cue, fxyze, byze
! eyz/byz = transverse electric/magnetic field in fourier space
! dcu = time derivative of electron current density
      complex, dimension(:,:), pointer :: eyz, byz, dcu
! ffc = form factor array for poisson solver
      complex, dimension(:), pointer :: ffc
!
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime
      real :: tdpost = 0.0, tfield = 0.0, tdjpost = 0.0, tpush = 0.0
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
      dth = 0.0
!
! allocate data for standard code
      allocate(part(idimp,np))
      allocate(qe(nxeh),qi(nxeh),fxe(nxeh))
      allocate(fxyze(3,nxeh),cue(2,nxeh),byze(2,nxeh))
      allocate(eyz(2,nxeh),byz(2,nxeh))
      allocate(ffc(nxhm))
!
! calculate form factors for gridless code: updates ffc
      call FFC1INITGL(ffc,ax,affp,nx,nxhm)
! calculate form factors for conventional PIC code: updates ffc
!     isign = 0
!     call POIS1GL(qe,fxe,isign,ffc,ax,affp,we,nx,nxhm,nxeh)
! initialize electrons
      call DISTR1H(part,vtx,vty,vtz,vx0,vy0,vz0,npx,idimp,np,nx,ipbc)
! initialize ion background
      qi = cmplx(0.0,0.0)
      call MDPOST1GL(part,qi,-qme,np,idimp,nx,nxhm,nxeh)
!
! initialize transverse electromagnetic fields
      eyz = cmplx(0.0,0.0)
      byz = cmplx(0.0,0.0)
      cux0 = 0.0
!
      if (dt > 0.63*ci*real(nxh)/real(nxhm)) then
         write (*,*) 'Info: Courant condition may be exceeded!'
      endif
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= ntime) go to 2000
!     write (*,*) 'ntime = ', ntime
!
! deposit current with standard procedure: updates part, cue, cux0
      call dtimer(dtime,itime,-1)
      cue = cmplx(0.0,0.0)
      if (relativity==1) then
         call MRDJPOST1GL(part,cue,cux0,qme,dth,ci,np,idimp,nx,nxhm,    &
     &nxeh)
      else
         call MDJPOST1GL(part,cue,cux0,qme,dth,np,idimp,nx,nxhm,nxeh)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdjpost = tdjpost + time
!
! deposit charge with standard procedure: updates qe
      call dtimer(dtime,itime,-1)
      qe = cmplx(0.0,0.0)
      call MDPOST1GL(part,qe,qme,np,idimp,nx,nxhm,nxeh)
      qe = qe + qi
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdpost = tdpost + time
!
! calculate electromagnetic fields in fourier space with standard
! procedure: updates eyz, byz
      call dtimer(dtime,itime,-1)
! initialization assumes initial acceleration of particles is zero
      if (ntime==0) then
         allocate(dcu(2,nxeh))
         dcu = cmplx(0.0,0.0)
! initialize with darwin electric and magnetic fields
! deposit time derivative of current: updates dcu
         if (relativity==1) then
            call MRD2JPOST1GL(part,dcu,qme,ci,np,idimp,nx,nxhm,nxeh)
         else
            call MD2JPOST1GL(part,dcu,qme,np,idimp,nx,nxhm,nxeh)
         endif
! initialize electromagnetic fields from free-streaming particles
! deposit scaled current and time derivative of current:
! updates cue, dcu
!        cue = cmplx(0.0,0.0)
!        if (relativity==1) then
!           call MRDSJPOST1GL(part,cue,dcu,qme,ci,np,idimp,nx,nxhm,nxeh)
!        else
!           call MDSJPOST1GL(part,cue,dcu,qme,ci,np,idimp,nx,nxhm,nxeh)
!        endif
! initialize magnetic field: updates byz, wm
         call IBPOIS13GL(cue,byz,ffc,ci,wm,nx,nxhm,nxeh)
         wf = 0.0
! initialize transverse electric field: updates eyz
         isign = 1
         call EPOIS13GL(dcu,eyz,isign,ffc,ax,affp,0.0,ci,wf,nx,nxhm,nxeh&
     &)
         deallocate(dcu)
         dth = 0.5*dt
      else
         call AMAXWEL1GL(eyz,byz,cue,ffc,ci,dt,wf,wm,nx,nxhm,nxeh)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
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
! add longitudinal and transverse electric fields with standard
! procedure: updates fxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call EMFIELD1GL(fxyze,fxe,eyz,ffc,nxhm,nxeh)
! copy magnetic field with standard procedure: updates byze
      isign = -1
      call BMFIELD1GL(byze,byz,ffc,nxhm,nxeh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! push particles with standard procedure: updates part, wke
      wke = 0.0
      call dtimer(dtime,itime,-1)
      if (relativity==1) then
! analytic Boris mover, gamma constant during time step
         call MARBPUSH13GL(part,fxyze,byze,omx,qbme,dt,dth,ci,wke,idimp,&
     &np,nx,nxhm,nxeh)
! exact analytic Boris mover, gamma varies during time step
!        call MEARBPUSH13GL(part,fxyze,byze,omx,qbme,dt,dth,ci,wke,idimp&
!    &,np,nx,nxhm,nxeh)
      else
         call MABPUSH13GL(part,fxyze,byze,omx,qbme,dt,dth,wke,idimp,np, &
     &nx,nxhm,nxeh)
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
      write (*,*) 'nxhm = ', nxhm
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
