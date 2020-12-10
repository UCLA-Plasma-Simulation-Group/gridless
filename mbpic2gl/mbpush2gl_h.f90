!-----------------------------------------------------------------------
! Interface file for mbpush2gl.f
      module mbpush2gl_h
      implicit none
!
      interface
         subroutine DISTR2H(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,idimp, &
     &nop,nx,ny,ipbc)
         implicit none
         integer, intent(in) :: npx, npy, idimp, nop, nx, ny, ipbc
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,nop), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine MBPUSH23GL(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx, &
     &ny,nxhm,nyhm,nxvh,nyv,ipbc)
         implicit none
         integer, intent(in) :: idimp, nop, nx, ny, nxhm, nyhm
         integer, intent(in) :: nxvh, nyv, ipbc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         complex, dimension(3,nxvh,nyv), intent(in) :: fxy, bxy
         end subroutine
      end interface
!
      interface
         subroutine MRBPUSH23GL(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,&
     &nx,ny,nxhm,nyhm,nxvh,nyv,ipbc)
         implicit none
         integer , intent(in):: idimp, nop, nx, ny, nxhm, nyhm
         integer, intent(in) :: nxvh, nyv, ipbc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         complex, dimension(3,nxvh,nyv), intent(in) :: fxy, bxy
         end subroutine
      end interface
!
      interface
         subroutine MABPUSH23GL(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,&
     &ny,nxhm,nyhm,nxvh,nyv,ipbc)
         implicit none
         integer, intent(in) :: idimp, nop, nx, ny, nxhm, nyhm
         integer, intent(in) :: nxvh, nyv, ipbc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         complex, dimension(3,nxvh,nyv), intent(in) :: fxy, bxy
         end subroutine
      end interface
!
      interface
         subroutine MARBPUSH23GL(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop&
     &,nx,ny,nxhm,nyhm,nxvh,nyv,ipbc)
         implicit none
         integer , intent(in):: idimp, nop, nx, ny, nxhm, nyhm
         integer, intent(in) :: nxvh, nyv, ipbc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         complex, dimension(3,nxvh,nyv), intent(in) :: fxy, bxy
         end subroutine
      end interface
!
      interface
         subroutine MEARBPUSH23GL(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,  &
     &nop,nx,ny,nxhm,nyhm,nxvh,nyv,ipbc)
         implicit none
         integer , intent(in):: idimp, nop, nx, ny, nxhm, nyhm
         integer, intent(in) :: nxvh, nyv, ipbc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         complex, dimension(3,nxvh,nyv), intent(in) :: fxy, bxy
         end subroutine
      end interface
!
      interface
         subroutine MDPOST2GL(part,q,sctx,qm,nop,idimp,nx,ny,nxhm,nyhm, &
     &nxvh,nyv)
         implicit none
         integer, intent(in) :: nop, idimp, nx, ny, nxhm, nyhm
         integer, intent(in) :: nxvh, nyv
         real, intent(in) :: qm
         real, dimension(idimp,nop), intent(in) :: part
         complex, dimension(nxvh,nyv), intent(inout) :: q
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine MDJPOST2GL(part,cu,sctx,qm,dt,nop,idimp,nx,ny,nxhm, &
     &nyhm,nxvh,nyv,ipbc)
         implicit none
         integer, intent(in) :: nop, idimp, nx, ny, nxhm, nyhm
         integer, intent(in) :: nxvh, nyv, ipbc
         real, intent(in) :: qm, dt
         real, dimension(idimp,nop), intent(inout) :: part
         complex, dimension(3,nxvh,nyv), intent(inout) :: cu
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine MRDJPOST2GL(part,cu,sctx,qm,dt,ci,nop,idimp,nx,ny,  &
     &nxhm,nyhm,nxvh,nyv,ipbc)
         implicit none
         integer , intent(in):: idimp, nop, nx, ny, nxhm, nyhm
         integer, intent(in) :: nxvh, nyv, ipbc
         real, intent(in) :: qm, dt, ci
         real, dimension(idimp,nop), intent(inout) :: part
         complex, dimension(3,nxvh,nyv), intent(inout) :: cu
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine MD2JPOST2GL(part,dcu,sctx,qm,nop,idimp,nx,ny,nxhm,  &
     &nyhm,nxvh,nyv)
         implicit none
         integer , intent(in):: idimp, nop, nx, ny, nxhm, nyhm
         integer, intent(in) :: nxvh, nyv
         real, intent(in) :: qm
         real, dimension(idimp,nop), intent(in) :: part
         complex, dimension(3,nxvh,nyv), intent(inout) :: dcu
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine MRD2JPOST2GL(part,dcu,sctx,qm,ci,nop,idimp,nx,ny,   &
     &nxhm,nyhm,nxvh,nyv)
         implicit none
         integer , intent(in):: idimp, nop, nx, ny, nxhm, nyhm
         integer, intent(in) :: nxvh, nyv
         real, intent(in) :: qm, ci
         real, dimension(idimp,nop), intent(in) :: part
         complex, dimension(3,nxvh,nyv), intent(inout) :: dcu
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine MDSJPOST2GL(part,cu,dcu,sctx,qm,ci,nop,idimp,nx,ny, &
     &nxhm,nyhm,nxvh,nyv)
         implicit none
         integer , intent(in):: idimp, nop, nx, ny, nxhm, nyhm
         integer, intent(in) :: nxvh, nyv
         real, intent(in) :: qm, ci
         real, dimension(idimp,nop), intent(in) :: part
         complex, dimension(3,nxvh,nyv), intent(inout) :: cu, dcu
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine MRDSJPOST2GL(part,cu,dcu,sctx,qm,ci,nop,idimp,nx,ny,&
     &nxhm,nyhm,nxvh,nyv)
         implicit none
         integer , intent(in):: idimp, nop, nx, ny, nxhm, nyhm
         integer, intent(in) :: nxvh, nyv
         real, intent(in) :: qm, ci
         real, dimension(idimp,nop), intent(in) :: part
         complex, dimension(3,nxvh,nyv), intent(inout) :: cu, dcu
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine FFC2INITGL(ffc,ax,ay,affp,nx,ny,nxhm,nyhm)
         implicit none
         integer, intent(in) :: nx, ny, nxhm, nyhm
         real, intent(in) :: ax, ay, affp
         complex, dimension(nxhm,nyhm), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine POIS23GL(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxhm,  &
     &nyhm,nxvh,nyv)
         implicit none
         integer, intent(in) :: isign, nx, ny, nxhm, nyhm, nxvh, nyv
         real, intent(in) :: ax, ay, affp
         real, intent(inout) :: we
         complex, dimension(nxvh,nyv), intent(in) :: q
         complex, dimension(3,nxvh,nyv), intent(inout) :: fxy
         complex, dimension(nxhm,nyhm), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine CUPERP2GL(cu,nx,ny,nxhm,nyhm,nxvh,nyv)
         implicit none
         integer, intent(in) :: nx, ny, nxhm, nyhm, nxvh, nyv
         complex, dimension(3,nxvh,nyv), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine IBPOIS23GL(cu,bxy,ffc,ci,wm,nx,ny,nxhm,nyhm,nxvh,nyv&
     &)
         implicit none
         integer, intent(in) :: nx, ny, nxhm, nyhm, nxvh, nyv
         real, intent(in) :: ci
         real, intent(inout) :: wm
         complex, dimension(3,nxvh,nyv), intent(in) :: cu
         complex, dimension(3,nxvh,nyv), intent(inout) :: bxy
         complex, dimension(nxhm,nyhm), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine AMAXWEL2GL(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxhm,   &
     &nyhm,nxvh,nyv)
         implicit none
         integer, intent(in) :: nx, ny, nxhm, nyhm, nxvh, nyv
         real, intent(in) :: ci, dt
         real, intent(inout) :: wf, wm
         complex, dimension(3,nxvh,nyv), intent(inout) :: exy, bxy
         complex, dimension(3,nxvh,nyv), intent(in) :: cu
         complex, dimension(nxhm,nyhm), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine EMFIELD2GL(fxy,exy,ffc,isign,nxhm,nyhm,nxvh,nyv)
         implicit none
         integer, intent(in) :: isign, nxhm, nyhm, nxvh, nyv
         complex, dimension(3,nxvh,nyv), intent(inout) :: fxy
         complex, dimension(3,nxvh,nyv), intent(in) :: exy
         complex, dimension(nxhm,nyhm), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine EPOIS23GL(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,&
     &ny,nxhm,nyhm,nxvh,nyv)
         implicit none
         integer, intent(in) :: isign, nx, ny, nxhm, nyhm, nxvh, nyv
         real, intent(in) :: ax, ay, affp, wp0, ci
         real, intent(inout) :: wf
         complex, dimension(3,nxvh,nyv), intent(in) :: dcu
         complex, dimension(3,nxvh,nyv), intent(inout) :: exy
         complex, dimension(nxhm,nyhm), intent(in) :: ffe
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
         subroutine M1BPUSH23GL(part,fxy,bxy,sctx,qbm,dt,dtc,ek,idimp,  &
     &nop,nx,ny,nxhm,nyhm,nxvh,nyv,ipbc)
         implicit none
         integer, intent(in) :: idimp, nop, nx, ny, nxhm, nyhm
         integer, intent(in) :: nxvh, nyv, ipbc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         complex, dimension(3,nxvh,nyv), intent(in) :: fxy, bxy
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine M1RBPUSH23GL(part,fxy,bxy,sctx,qbm,dt,dtc,ci,ek,    &
     &idimp,nop,nx,ny,nxhm,nyhm,nxvh,nyv,ipbc)
         implicit none
         integer , intent(in):: idimp, nop, nx, ny, nxhm, nyhm
         integer, intent(in) :: nxvh, nyv, ipbc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         complex, dimension(3,nxvh,nyv), intent(in) :: fxy, bxy
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine MD4JPOST2GL(part,dcu,dcu2,sctx,qm,nop,idimp,nx,ny,  &
     &nxhm,nyhm,nxvh,nyv)
         implicit none
         integer, intent(in) :: nop, idimp, nx, ny, nxhm, nyhm
         integer, intent(in) :: nxvh, nyv
         real, intent(in) :: qm
         real, dimension(idimp,nop), intent(in) :: part
         complex, dimension(3,nxvh,nyv), intent(inout) :: dcu, dcu2
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine MRD4JPOST2GL(part,dcu,dcu2,sctx,qm,ci,nop,idimp,nx, &
     &ny,nxhm,nyhm,nxvh,nyv)
         implicit none
         integer , intent(in):: idimp, nop, nx, ny, nxhm, nyhm
         integer, intent(in) :: nxvh, nyv
         real, intent(in) :: qm, ci
         real, dimension(idimp,nop), intent(in) :: part
         complex, dimension(3,nxvh,nyv), intent(inout) :: dcu, dcu2
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine A4MAXWEL2GL(exy,bxy,cu,dcu,dcu2,ffc,ci,dt,wf,wm,nx, &
     &ny,nxhm,nyhm,nxvh,nyv)
         implicit none
         integer, intent(in) :: nx, ny, nxhm, nyhm, nxvh, nyv
         real, intent(in) :: ci, dt
         real, intent(inout) :: wf, wm
         complex, dimension(3,nxvh,nyv), intent(inout) :: exy, bxy
         complex, dimension(3,nxvh,nyv), intent(in) :: cu, dcu, dcu2
         complex, dimension(nxhm,nyhm), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MAXWEL2GL(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxhm,nyhm&
     &,nxvh,nyv)
         implicit none
         integer, intent(in) :: nx, ny, nxhm, nyhm, nxvh, nyv
         real, intent(in) :: ci, dt
         real, intent(inout) :: wf, wm
         complex, dimension(3,nxvh,nyv), intent(inout) :: exy, bxy
         complex, dimension(3,nxvh,nyv), intent(in) :: cu
         complex, dimension(nxhm,nyhm), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine ELFIELD23GL(q,fxy,ffc,we,nx,ny,nxhm,nyhm,nxvh,nyv)
         implicit none
         integer, intent(in) :: nx, ny, nxhm, nyhm, nxvh, nyv
         real, intent(inout) :: we
         complex, dimension(nxvh,nyv), intent(in) :: q
         complex, dimension(3,nxvh,nyv), intent(inout) :: fxy
         complex, dimension(nxhm,nyhm), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MEVFIELD23GL(fxy,gxy,sctx,xp,yp,nx,ny,nxhm,nyhm,nxvh&
     &,nyv)
         implicit none
         integer, intent(in) :: nx, ny, nxhm, nyhm, nxvh, nyv
         real, intent(in) ::  xp, yp
         complex, dimension(3,nxvh,nyv), intent(in) :: fxy
         real, dimension(3), intent(inout) ::  gxy
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine MA0RBPUSH23GL(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,  &
     &nop,nx,ny,nxhm,nyhm,nxvh,nyv,ipbc)
         implicit none
         integer , intent(in):: idimp, nop, nx, ny, nxhm, nyhm
         integer, intent(in) :: nxvh, nyv, ipbc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         complex, dimension(3,nxvh,nyv), intent(in) :: fxy, bxy
         end subroutine
      end interface
!
      end module
