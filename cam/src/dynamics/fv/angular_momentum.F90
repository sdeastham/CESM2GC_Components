module angular_momentum

! Compute angular momentum using FV data structures.   

use shr_kind_mod,  only: r8 => shr_kind_r8
use dynamics_vars, only: T_FVDYCORE_GRID, T_FVDYCORE_CONSTANTS

#if defined( SPMD )
use mod_comm,      only: mp_send3d, mp_recv3d
#endif

implicit none
private
save

public :: &
   am_calc

!=========================================================================================
CONTAINS
!=========================================================================================

subroutine am_calc(constants, grid, dt, u, delp, &
                   pe, peln, phis, r_vir, cp,    &
                   rg, ame, am0)

   ! THT 16.11.04  Determines the column and globally integrated total AM.
   !               Filched and perverted from benergy.F90

   ! INPUT PARAMETERS:
   type (T_FVDYCORE_GRID)     , intent(in) :: grid         ! YZ decomposition
   type (T_FVDYCORE_CONSTANTS), intent(in) :: constants ! For convenience

   real(r8), intent(in) :: dt ! time-step
   real(r8), intent(in) :: &
      u   (grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km),   &! U-winds
      delp(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km),   &! Delta pressure
      pe  (grid%ifirstxy:grid%ilastxy, grid%km+1, grid%jfirstxy:grid%jlastxy), &! Edge pressure
      peln(grid%ifirstxy:grid%ilastxy, grid%km+1, grid%jfirstxy:grid%jlastxy), &! Edge pressure
      phis(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy)              ! Surface heights

   real(r8), intent(in) :: r_vir  ! Virtual effect constant ( rwv/rg-1 )
   real(r8), intent(in) :: cp     ! C_p ( = rg / cappa )
   real(r8), intent(in) :: rg     ! Gas constant for dry air

   ! OUTPUT PARAMETERS:
   real(r8), intent(out) :: ame(grid%jm) ! column integrated Total AM
   real(r8), intent(out) :: am0          ! globally integrated Total AM

   ! Local
   real (r8), parameter :: D0_0        = 0.0_r8
   real (r8), parameter :: D0_25       = 0.25_r8
   real (r8), parameter :: D0_5        = 0.5_r8
   real (r8), parameter :: D1_0        = 1.0_r8

   integer   :: im, jm, km, ifirstxy, ilastxy, jfirstxy, jlastxy
   integer   :: iam, myidxy_x, myidxy_y, nprxy_x, nprxy_y, dest, src    ! SPMD related
   integer   :: i, j, k, js1g0, js2g0, jn1g0, jn1g1, jn2g0, ktot, jtot, itot

   real (r8) :: u2(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy+1)

   real (r8) :: bte(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy)
   real (r8) :: te_sp(grid%ifirstxy:grid%ilastxy,grid%km)
   real (r8) :: te_np(grid%ifirstxy:grid%ilastxy,grid%km)
   real (r8) :: gztop(grid%ifirstxy:grid%ilastxy)
   real (r8) :: xsum(grid%jfirstxy:grid%jlastxy)
   real (r8) :: sp_sum(grid%km), np_sum(grid%km)
   real (r8) :: tmp
   real (r8) :: te(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)
   real(r8)  :: unorth(grid%ifirstxy:grid%ilastxy,grid%km)    ! North halo
   real(r8) :: pheast(grid%jfirstxy:grid%jlastxy), phwest(grid%jfirstxy:grid%jlastxy) ! more haloes
   real(r8) :: oma
   !----------------------------------------------------------------------------

   oma=constants%ae*constants%omega

   im     = grid%im
   jm     = grid%jm
   km     = grid%km

   ifirstxy = grid%ifirstxy
   ilastxy  = grid%ilastxy
   jfirstxy = grid%jfirstxy
   jlastxy  = grid%jlastxy

   iam      = grid%iam
   myidxy_x = grid%myidxy_x
   myidxy_y = grid%myidxy_y
   nprxy_x  = grid%nprxy_x
   nprxy_y  = grid%nprxy_y
      
   js1g0  = max(1,jfirstxy)
   js2g0  = max(2,jfirstxy)
   jn2g0  = min(jm-1,jlastxy)
   jn1g0  = min(jm,jlastxy)
   jn1g1  = min(jm,jlastxy+1)

   itot   = ilastxy - ifirstxy + 1
   jtot   = jlastxy - jfirstxy + 1

#if defined(SPMD)
   call mp_send3d( grid%commxy, iam-nprxy_x, iam+nprxy_x, im, jm, km,      &
                   ifirstxy, ilastxy, jfirstxy, jlastxy, 1, km,           &
                   ifirstxy, ilastxy, jfirstxy, jfirstxy, 1, km, u )
   call mp_recv3d( grid%commxy, iam+nprxy_x, im, jm, km,                   &
                   ifirstxy, ilastxy, jlastxy+1, jlastxy+1, 1, km,        &
                   ifirstxy, ilastxy, jlastxy+1, jlastxy+1, 1, km, unorth )
#endif

!$omp parallel do private(i, j, k, u2)
   do k = 1, km

      ! Check the poles for consistent values
      do j=js2g0,jlastxy
         do i=ifirstxy,ilastxy
            u2(i,j) = grid%cose(j)**2 *(u(i,j,k)+oma*grid%cose(j)) 
         enddo
      enddo

      if ( jlastxy /= jm ) then    ! Pull information out of northern halo
         do i=ifirstxy,ilastxy
            u2(i,jlastxy+1) = grid%cose(jlastxy+1)**2 *(unorth(i,k)+oma*grid%cose(jlastxy+1)) 
         enddo
      endif

      do j=js2g0,jn2g0
         do i=ifirstxy,ilastxy
            te(i,j,k) = D0_5*(u2(i,j)+u2(i,j+1))*grid%acosu(j) 
         enddo
      enddo

      do j=js2g0,jn2g0
         do i=ifirstxy, ilastxy
            te(i,j,k) = delp(i,j,k) * te(i,j,k) 
         enddo
      enddo

      if ( jfirstxy == 1 ) then
         do i=ifirstxy,ilastxy
            te_sp(i,k) =    u2(i,2)/grid%cose(2) 
         enddo
      endif

      if ( jlastxy == jm ) then
         do i=ifirstxy,ilastxy
            te_np(i,k)=     u2(i,jm)/grid%cose(jm) 
         enddo
      endif

   enddo


   if ( jfirstxy == 1 ) then
      call par_xsum( grid, te_sp, km, sp_sum )
!$omp parallel do private(i, k, tmp)
      do k=1,km
         tmp = delp(ifirstxy,1,k) * sp_sum(k)/real(im,r8)
         do i=ifirstxy,ilastxy
            te(i,1,k)  = tmp
         enddo
      enddo
   endif

   if ( jlastxy == jm ) then
      call par_xsum( grid, te_np, km, np_sum )
!$omp parallel do private(i, k, tmp)
      do k=1,km
         tmp = delp(ifirstxy,jm,k) * np_sum(k)/real(im,r8)
         do i=ifirstxy,ilastxy
            te(i,jm,k) = tmp
         enddo
      enddo
   endif

   ! source term (mountain torque)
#if defined(SPMD)
   dest = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
   src  = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
   call mp_send3d( grid%commxy, dest, src, im, jm, km,                  &
                   ifirstxy, ilastxy, jfirstxy, jlastxy, 1,  1,        &
                   ifirstxy, ifirstxy, jfirstxy, jlastxy, 1,  1, phis)
   call mp_recv3d( grid%commxy, src, im, jm, km,                        &
                   ilastxy+1, ilastxy+1, jfirstxy, jlastxy, 1,  1,     &
                   ilastxy+1, ilastxy+1, jfirstxy, jlastxy, 1,  1, pheast )
   dest = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
   src  = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
   call mp_send3d( grid%commxy, dest, src, im, jm, km,                  &
                   ifirstxy, ilastxy, jfirstxy, jlastxy, 1,  1,        &
                   ilastxy, ilastxy, jfirstxy, jlastxy, 1,  1, phis)
   call mp_recv3d( grid%commxy, src, im, jm, km,                        &
                   ifirstxy-1, ifirstxy-1, jfirstxy, jlastxy, 1,  1,     &
                   ifirstxy-1, ifirstxy-1, jfirstxy, jlastxy, 1,  1, phwest )
#endif
   bte = D0_0
!$omp parallel do private(i,j,k)
   do j=jfirstxy,jlastxy
      if (j == 1) then
         ame(1) = D0_0 
         do k=1,km
            ame(1) = ame(1) + te(ifirstxy,1,k)
         enddo
         ame(1)  = grid%acap * ame(1)
      elseif (j == jm) then
         ame(jm) = D0_0 
         do k=1,km
            ame(jm) = ame(jm) + te(ifirstxy,jm,k)
         enddo
         ame(jm) = grid%acap * ame(jm)
      else
         do i=ifirstxy+1,ilastxy-1
            !  bte(i,j) = D0_5*(pe(i+1,km+1,j)+pe(i,km+1,j))*(phis(i+1,j)-phis(i,j)) ! W'ward torque
            !  bte(i,j) = D0_5*(pe(i-1,km+1,j)-pe(i+1,km+1,j))*phis(i,j) ! equivalent form
            bte(i,j) = D0_5*pe(i,km+1,j)*(phis(i+1,j)-phis(i-1,j))    ! also equivalent 
         enddo
         if (itot.eq.im) then 
            pheast(j)=phis(1,j)
            phwest(j)=phis(im,j)
         endif
         i=ilastxy
         bte(i,j) = D0_5*pe(i,km+1,j)*(pheast(j)-phis(i-1,j))   ! equivalent
         i=ifirstxy
         bte(i,j) = D0_5*pe(i,km+1,j)*(phis(i+1,j)-phwest(j))   ! equivalent
         do i=ifirstxy,ilastxy
            bte(i,j) = dt/(constants%ae*grid%dl) * bte(i,j)   ! W'ward AM from form drag
         enddo
         do k=1,km
            do i=ifirstxy,ilastxy
               bte(i,j) = bte(i,j) + te(i,j,k)*grid%cosp(j) ! add back AM lost to form drag
            enddo
         enddo
      endif
   enddo

   call par_xsum(grid, bte, jtot, xsum)
  
!$omp parallel do private(j)
   do j=js2g0,jn2g0
      ame(j) = xsum(j)
   enddo

   call par_vecsum(jm, jfirstxy, jlastxy, ame, am0, grid%commxy_y, grid%nprxy_y)

end subroutine am_calc

!=========================================================================================

end module angular_momentum
