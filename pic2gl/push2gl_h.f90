!-----------------------------------------------------------------------
! Interface file for push2gl.f
      module push2gl_h
      implicit none
!
      interface
         subroutine DISTR2(part,vtx,vty,vdx,vdy,npx,npy,idimp,nop,nx,ny,&
     &ipbc)
         implicit none
         integer, intent(in) :: npx, npy, idimp, nop, nx, ny, ipbc
         real, intent(in) :: vtx, vty, vdx, vdy
         real, dimension(idimp,nop), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine PUSH2GL(part,fxy,sctx,qbm,dt,ek,idimp,nop,nx,ny,nxhm&
     &,nyhm,nxvh,nyv,ipbc)
         implicit none
         integer, intent(in) :: idimp, nop, nx, ny, nxhm, nyhm
         integer, intent(in) :: nxvh, nyv, ipbc
         real, intent(in) :: qbm, dt
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         complex, dimension(2,nxvh,nyv), intent(in) :: fxy
         complex, dimension(nxvh), intent(inout) :: sctx
         end subroutine
      end interface
!
      interface
         subroutine DPOST2GL(part,q,sctx,qm,nop,idimp,nx,ny,nxhm,nyhm,  &
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
         subroutine FFC2INITGL(ffc,ax,ay,affp,nx,ny,nxhm,nyhm)
         implicit none
         integer, intent(in) :: nx, ny, nxhm, nyhm
         real, intent(in) :: ax, ay, affp
         complex, dimension(nxhm,nyhm), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine POIS22GL(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxhm,  &
     &nyhm,nxvh,nyv)
         implicit none
         integer, intent(in) :: isign, nx, ny, nxhm, nyhm, nxvh, nyv
         real, intent(in) :: ax, ay, affp
         real, intent(inout) :: we
         complex, dimension(nxvh,nyv), intent(in) :: q
         complex, dimension(2,nxvh,nyv), intent(inout) :: fxy
         complex, dimension(nxhm,nyhm), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine EVFIELD22GL(fxy,gxy,sctx,xp,yp,nx,ny,nxhm,nyhm,nxvh,&
     &nyv)
         implicit none
         integer, intent(in) :: nx, ny, nxhm, nyhm, nxvh, nyv
         real, intent(in) ::  xp, yp
         complex, dimension(2,nxvh,nyv), intent(in) :: fxy
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
      end module
