!-----------------------------------------------------------------------
! Interface file for bpush1gl.f
      module bpush1gl_h
      implicit none
!
      interface
         subroutine DISTR1H(part,vtx,vty,vtz,vdx,vdy,vdz,npx,idimp,nop, &
     &nx,ipbc)
         implicit none
         integer, intent(in) :: npx, idimp, nop, nx, ipbc
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,nop), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine BPUSH13GL(part,fxyz,byz,sctx,omx,qbm,dt,dtc,ek,idimp&
     &,nop,nx,nxhm,nxvh)
         implicit none
         integer, intent(in) :: idimp, nop, nx, nxhm, nxvh
         real, intent(in) :: omx, qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         complex, dimension(3,nxvh), intent(in) :: fxyz
         complex, dimension(2,nxvh), intent(in) :: byz
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine RBPUSH13GL(part,fxyz,byz,sctx,omx,qbm,dt,dtc,ci,ek, &
     &idimp,nop,nx,nxhm,nxvh)
         implicit none
         integer, intent(in) :: idimp, nop, nx, nxhm, nxvh
         real, intent(in) :: omx, qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         complex, dimension(3,nxvh), intent(in) :: fxyz
         complex, dimension(2,nxvh), intent(in) :: byz
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine ABPUSH13GL(part,fxyz,byz,sctx,omx,qbm,dt,dtc,ek,    &
     &idimp,nop,nx,nxhm,nxvh)
         implicit none
         integer, intent(in) :: idimp, nop, nx, nxhm, nxvh
         real, intent(in) :: omx, qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         complex, dimension(3,nxvh), intent(in) :: fxyz
         complex, dimension(2,nxvh), intent(in) :: byz
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine ARBPUSH13GL(part,fxyz,byz,sctx,omx,qbm,dt,dtc,ci,ek,&
     &idimp,nop,nx,nxhm,nxvh)
         implicit none
         integer, intent(in) :: idimp, nop, nx, nxhm, nxvh
         real, intent(in) :: omx, qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         complex, dimension(3,nxvh), intent(in) :: fxyz
         complex, dimension(2,nxvh), intent(in) :: byz
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine EARBPUSH13GL(part,fxyz,byz,sctx,omx,qbm,dt,dtc,ci,ek&
     &,idimp,nop,nx,nxhm,nxvh)
         implicit none
         integer, intent(in) :: idimp, nop, nx, nxhm, nxvh
         real, intent(in) :: omx, qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         complex, dimension(3,nxvh), intent(in) :: fxyz
         complex, dimension(2,nxvh), intent(in) :: byz
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine DPOST1GL(part,q,sctx,qm,nop,idimp,nx,nxhm,nxvh)
         implicit none
         integer, intent(in) :: nop, idimp, nx,nxhm,nxvh
         real, intent(in) :: qm
         real, dimension(idimp,nop), intent(in) :: part
         complex, dimension(nxvh), intent(inout) :: q, sctx
         end subroutine
      end interface
!
      interface
         subroutine DJPOST1GL(part,cu,sctx,cux0,qm,dt,nop,idimp,nx,nxhm,&
     &nxvh)
         implicit none
         integer, intent(in) :: nop, idimp, nx, nxhm, nxvh
         real, intent(in) :: qm, dt
         real, intent(inout) :: cux0
         real, dimension(idimp,nop), intent(inout) :: part
         complex, dimension(2,nxvh), intent(inout) :: cu
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine RDJPOST1GL(part,cu,sctx,cux0,qm,dt,ci,nop,idimp,nx, &
     &nxhm,nxvh)
         implicit none
         integer, intent(in) :: nop, idimp, nx, nxhm, nxvh
         real, intent(in) :: qm, dt, ci
         real, intent(inout) :: cux0
         real, dimension(idimp,nop), intent(inout) :: part
         complex, dimension(2,nxvh), intent(inout) :: cu
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine D2JPOST1GL(part,dcu,sctx,qm,nop,idimp,nx,nxhm,nxvh)
         implicit none
         integer , intent(in):: idimp, nop, nx, nxhm, nxvh
         real, intent(in) :: qm
         real, dimension(idimp,nop), intent(in) :: part
         complex, dimension(2,nxvh), intent(inout) :: dcu
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine RD2JPOST1GL(part,dcu,sctx,qm,ci,nop,idimp,nx,nxhm,  &
     &nxvh)
         implicit none
         integer , intent(in):: idimp, nop, nx, nxhm, nxvh
         real, intent(in) :: qm, ci
         real, dimension(idimp,nop), intent(in) :: part
         complex, dimension(2,nxvh), intent(inout) :: dcu
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine DSJPOST1GL(part,cu,dcu,sctx,qm,ci,nop,idimp,nx,nxhm,&
     &nxvh)
         implicit none
         integer , intent(in):: idimp, nop, nx, nxhm, nxvh
         real, intent(in) :: qm, ci
         real, dimension(idimp,nop), intent(in) :: part
         complex, dimension(2,nxvh), intent(inout) :: cu, dcu
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine RDSJPOST1GL(part,cu,dcu,sctx,qm,ci,nop,idimp,nx,nxhm&
     &,nxvh)
         implicit none
         integer , intent(in):: idimp, nop, nx, nxhm, nxvh
         real, intent(in) :: qm, ci
         real, dimension(idimp,nop), intent(in) :: part
         complex, dimension(2,nxvh), intent(inout) :: cu, dcu
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine FFC1INITGL(ffc,ax,affp,nx,nxhm)
         implicit none
         integer, intent(in) :: nx, nxhm
         real, intent(in) :: ax, affp
         complex, dimension(nxhm), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine POIS1GL(q,fx,isign,ffc,ax,affp,we,nx,nxhm,nxvh)
         implicit none
         integer, intent(in) :: isign, nx, nxhm, nxvh
         real, intent(in) :: ax, affp
         real, intent(inout) :: we
         complex, dimension(nxvh), intent(in) :: q
         complex, dimension(nxvh), intent(inout) :: fx
         complex, dimension(nxhm), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine IBPOIS13GL(cu,byz,ffc,ci,wm,nx,nxhm,nxvh)
         implicit none
         integer, intent(in) :: nx, nxhm, nxvh
         real, intent(in) :: ci
         real, intent(inout) :: wm
         complex, dimension(2,nxvh), intent(in) :: cu
         complex, dimension(2,nxvh), intent(inout) :: byz
         complex, dimension(nxhm), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine AMAXWEL1GL(eyz,byz,cu,ffc,ci,dt,wf,wm,nx,nxhm,nxvh)
         implicit none
         integer, intent(in) :: nx, nxhm, nxvh
         real, intent(in) :: ci, dt
         real, intent(inout) :: wf, wm
         complex, dimension(2,nxvh), intent(inout) :: eyz, byz
         complex, dimension(2,nxvh), intent(in) :: cu
         complex, dimension(nxhm), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine EMFIELD1GL(fxyz,fx,eyz,ffc,nxhm,nxvh)
         implicit none
         integer, intent(in) :: nxhm, nxvh
         complex, dimension(3,nxvh), intent(inout) :: fxyz
         complex, dimension(nxvh), intent(in) :: fx
         complex, dimension(3,nxvh), intent(in) :: eyz
         complex, dimension(nxhm), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine BMFIELD1GL(fyz,eyz,ffc,nxhm,nxvh)
         implicit none
         integer, intent(in) :: nxhm, nxvh
         complex, dimension(2,nxvh), intent(inout) :: fyz
         complex, dimension(2,nxvh), intent(in) :: eyz
         complex, dimension(nxhm), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine EPOIS13GL(dcu,eyz,isign,ffe,ax,affp,wp0,ci,wf,nx,   &
     &nxhm,nxvh)
         implicit none
         integer, intent(in) :: isign, nx, nxhm, nxvh
         real, intent(in) :: ax, affp, wp0, ci
         real, intent(inout) :: wf
         complex, dimension(2,nxvh), intent(in) :: dcu
         complex, dimension(2,nxvh), intent(inout) :: eyz
         complex, dimension(nxhm), intent(in) :: ffe
         end subroutine
      end interface
!
      interface
         subroutine EVFIELD13GL(fxy,gxy,sctx,xp,nx,nxhm,nxvh)
         implicit none
         integer, intent(in) :: nx, nxhm, nxvh
         real, intent(in) ::  xp
         complex, dimension(3,nxvh), intent(in) :: fxy
         real, dimension(3), intent(inout) ::  gxy
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         function ranorm()
         implicit none
         double precision :: ranorm
         end function
      end interface
!
      interface
         function randum()
         implicit none
         double precision :: randum
         end function
      end interface
!
      interface
         subroutine MJPOST1GL(part,amu,sctx,qm,nop,idimp,nx,nxhm,nxvh)
         implicit none
         integer, intent(in) :: nop, idimp, nx, nxhm, nxvh
         real, intent(in) :: qm
         real, dimension(idimp,nop), intent(in) :: part
         complex, dimension(2,nxvh), intent(inout) :: amu
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine RMJPOST1GL(part,amu,sctx,qm,ci,nop,idimp,nx,nxhm,   &
     &nxvh)
         implicit none
         integer , intent(in):: idimp, nop, nx, nxhm, nxvh
         real, intent(in) :: qm, ci
         real, dimension(idimp,nop), intent(in) :: part
         complex, dimension(2,nxvh), intent(inout) :: amu
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine DCUPERP13GL(dcu,amu,nx,nxhm,nxvh)
         implicit none
         integer, intent(in) :: nx, nxhm, nxvh
         complex, dimension(2,nxvh), intent(inout) :: dcu
         complex, dimension(2,nxvh), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine MAXWEL1GL(eyz,byz,cu,ffc,ci,dt,wf,wm,nx,nxhm,nxvh)
         implicit none
         integer, intent(in) :: nx, nxhm, nxvh
         real, intent(in) :: ci, dt
         real, intent(inout) :: wf, wm
         complex, dimension(2,nxvh), intent(inout) :: eyz, byz
         complex, dimension(2,nxvh), intent(in) :: cu
         complex, dimension(nxhm), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine ELFIELD1GL(q,fx,ffc,we,nx,nxhm,nxvh)
         implicit none
         integer, intent(in) :: nx, nxhm, nxvh
         real, intent(inout) :: we
         complex, dimension(nxvh), intent(in) :: q
         complex, dimension(nxvh), intent(inout) :: fx
         complex, dimension(nxhm), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
      subroutine PSMOOTH1GL(q,qs,ffc,nxhm,nxvh)
         implicit none
         integer, intent(in) :: nxhm, nxvh
         complex, dimension(nxvh), intent(in) :: q
         complex, dimension(nxvh), intent(inout) :: qs
         complex, dimension(nxhm), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine POYNT1GL(q,eyz,byz,ffc,ex0,sx,sy,sz,nx,nxhm,nxvh)
         implicit none
         integer, intent(in) :: nx, nxhm, nxvh
         real, intent(inout) :: ex0, sx, sy, sz
         complex, dimension(nxvh), intent(in) :: q
         complex, dimension(2,nxvh), intent(in) :: eyz, byz
         complex, dimension(nxhm), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine A0RBPUSH13GL(part,fxyz,byz,sctx,omx,qbm,dt,dtc,ci,ek&
     &,idimp,nop,nx,nxhm,nxvh)
         implicit none
         integer, intent(in) :: idimp, nop, nx, nxhm, nxvh
         real, intent(in) :: omx, qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         complex, dimension(3,nxvh), intent(in) :: fxyz
         complex, dimension(2,nxvh), intent(in) :: byz
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      end module
