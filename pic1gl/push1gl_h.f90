!-----------------------------------------------------------------------
! Interface file for push1gl.f
      module push1gl_h
      implicit none
!
      interface
         subroutine DISTR1(part,vtx,vdx,npx,idimp,nop,nx,ipbc)
         implicit none
         integer, intent(in) :: npx, idimp, nop, nx, ipbc
         real, intent(in) :: vtx, vdx
         real, dimension(idimp,nop), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine PUSH1GL(part,fx,sctx,qbm,dt,ek,idimp,nop,nx,nxhm,   &
     &nxvh)
         implicit none
         integer, intent(in) :: idimp, nop, nx, nxhm, nxvh
         real, intent(in) :: qbm, dt
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         complex, dimension(nxvh), intent(in) :: fx
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
