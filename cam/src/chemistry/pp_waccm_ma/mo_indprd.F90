      module mo_indprd
      use shr_kind_mod, only : r8 => shr_kind_r8
      private
      public :: indprd
      contains
      subroutine indprd( class, prod, nprod, y, extfrc, rxt, ncol )
      use chem_mods, only : gas_pcnst, extcnt, rxntot
      use ppgrid, only : pver
      implicit none
!--------------------------------------------------------------------
! ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: class
      integer, intent(in) :: ncol
      integer, intent(in) :: nprod
      real(r8), intent(in) :: y(ncol,pver,gas_pcnst)
      real(r8), intent(in) :: rxt(ncol,pver,rxntot)
      real(r8), intent(in) :: extfrc(ncol,pver,extcnt)
      real(r8), intent(inout) :: prod(ncol,pver,nprod)
!--------------------------------------------------------------------
! ... "independent" production for Explicit species
!--------------------------------------------------------------------
      if( class == 1 ) then
         prod(:,:,1) = 0._r8
         prod(:,:,2) =rxt(:,:,153)*y(:,:,10)*y(:,:,8)
         prod(:,:,3) = 0._r8
         prod(:,:,4) = 0._r8
         prod(:,:,5) = 0._r8
         prod(:,:,6) = 0._r8
         prod(:,:,7) = 0._r8
         prod(:,:,8) = 0._r8
         prod(:,:,9) = 0._r8
         prod(:,:,10) = 0._r8
         prod(:,:,11) = 0._r8
         prod(:,:,12) = 0._r8
         prod(:,:,13) = 0._r8
         prod(:,:,14) = 0._r8
         prod(:,:,15) = 0._r8
         prod(:,:,16) = 0._r8
         prod(:,:,17) = 0._r8
         prod(:,:,18) = 0._r8
         prod(:,:,19) = 0._r8
         prod(:,:,20) = 0._r8
         prod(:,:,21) = (rxt(:,:,233)*y(:,:,22) +rxt(:,:,234)*y(:,:,22))*y(:,:,19)
         prod(:,:,22) = 0._r8
         prod(:,:,23) = 0._r8
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,:,49) = 0._r8
         prod(:,:,44) = (rxt(:,:,58) +rxt(:,:,86))*y(:,:,63) +.180_r8*rxt(:,:,60) &
                 *y(:,:,15)
         prod(:,:,47) =rxt(:,:,5)*y(:,:,7)
         prod(:,:,32) = 0._r8
         prod(:,:,3) = 0._r8
         prod(:,:,2) = 0._r8
         prod(:,:,35) =1.440_r8*rxt(:,:,60)*y(:,:,15)
         prod(:,:,16) = (rxt(:,:,58) +rxt(:,:,86))*y(:,:,63) +.380_r8*rxt(:,:,60) &
                 *y(:,:,15) + extfrc(:,:,3)
         prod(:,:,28) = (rxt(:,:,70) +.800_r8*rxt(:,:,73) +rxt(:,:,82) + &
                 .800_r8*rxt(:,:,85)) + extfrc(:,:,9)
         prod(:,:,48) = + extfrc(:,:,1)
         prod(:,:,40) = + extfrc(:,:,2)
         prod(:,:,43) =.330_r8*rxt(:,:,60)*y(:,:,15) + extfrc(:,:,11)
         prod(:,:,45) = 0._r8
         prod(:,:,31) = 0._r8
         prod(:,:,12) = 0._r8
         prod(:,:,10) = 0._r8
         prod(:,:,50) =rxt(:,:,59)*y(:,:,15) +rxt(:,:,37)*y(:,:,41) +rxt(:,:,48) &
                 *y(:,:,42)
         prod(:,:,11) = 0._r8
         prod(:,:,34) =.180_r8*rxt(:,:,60)*y(:,:,15)
         prod(:,:,38) = (rxt(:,:,59) +.330_r8*rxt(:,:,60))*y(:,:,15)
         prod(:,:,41) = 0._r8
         prod(:,:,18) = 0._r8
         prod(:,:,33) =.050_r8*rxt(:,:,60)*y(:,:,15)
         prod(:,:,42) =rxt(:,:,37)*y(:,:,41) +2.000_r8*rxt(:,:,40)*y(:,:,43) &
                  +2.000_r8*rxt(:,:,41)*y(:,:,44) +2.000_r8*rxt(:,:,42)*y(:,:,45) &
                  +rxt(:,:,45)*y(:,:,46) +4.000_r8*rxt(:,:,38)*y(:,:,47) &
                  +3.000_r8*rxt(:,:,39)*y(:,:,48) +rxt(:,:,50)*y(:,:,50) +rxt(:,:,46) &
                 *y(:,:,51) +rxt(:,:,47)*y(:,:,52) +2.000_r8*rxt(:,:,43)*y(:,:,53) &
                  +rxt(:,:,44)*y(:,:,54)
         prod(:,:,6) = 0._r8
         prod(:,:,37) = 0._r8
         prod(:,:,4) = 0._r8
         prod(:,:,1) = 0._r8
         prod(:,:,46) = 0._r8
         prod(:,:,29) = 0._r8
         prod(:,:,30) = 0._r8
         prod(:,:,8) = 0._r8
         prod(:,:,39) =rxt(:,:,48)*y(:,:,42) +rxt(:,:,49)*y(:,:,49) +rxt(:,:,50) &
                 *y(:,:,50) +2.000_r8*rxt(:,:,53)*y(:,:,55) +2.000_r8*rxt(:,:,54) &
                 *y(:,:,56) +3.000_r8*rxt(:,:,51)*y(:,:,57) +2.000_r8*rxt(:,:,52) &
                 *y(:,:,58)
         prod(:,:,36) = 0._r8
         prod(:,:,27) = 0._r8
         prod(:,:,22) = 0._r8
         prod(:,:,17) = 0._r8
         prod(:,:,20) = (rxt(:,:,66) +rxt(:,:,78)) + extfrc(:,:,7)
         prod(:,:,23) = + extfrc(:,:,5)
         prod(:,:,19) = (rxt(:,:,70) +rxt(:,:,71) +rxt(:,:,82) +rxt(:,:,83)) &
                  + extfrc(:,:,6)
         prod(:,:,21) = + extfrc(:,:,4)
         prod(:,:,24) = 0._r8
         prod(:,:,15) = (rxt(:,:,71) +1.200_r8*rxt(:,:,73) +rxt(:,:,83) + &
                 1.200_r8*rxt(:,:,85)) + extfrc(:,:,8)
         prod(:,:,25) = (rxt(:,:,66) +rxt(:,:,70) +rxt(:,:,71) +rxt(:,:,78) + &
                 rxt(:,:,82) +rxt(:,:,83)) + extfrc(:,:,10)
         prod(:,:,13) = 0._r8
         prod(:,:,14) = 0._r8
         prod(:,:,5) =rxt(:,:,41)*y(:,:,44) +rxt(:,:,42)*y(:,:,45) +rxt(:,:,45) &
                 *y(:,:,46) +rxt(:,:,49)*y(:,:,49) +rxt(:,:,50)*y(:,:,50) +rxt(:,:,47) &
                 *y(:,:,52) +2.000_r8*rxt(:,:,43)*y(:,:,53) +2.000_r8*rxt(:,:,44) &
                 *y(:,:,54) +rxt(:,:,53)*y(:,:,55) +2.000_r8*rxt(:,:,54)*y(:,:,56)
         prod(:,:,7) =rxt(:,:,40)*y(:,:,43) +rxt(:,:,42)*y(:,:,45) +rxt(:,:,46) &
                 *y(:,:,51)
         prod(:,:,9) = 0._r8
         prod(:,:,26) =rxt(:,:,49)*y(:,:,49) +rxt(:,:,44)*y(:,:,54)
      end if
      end subroutine indprd
      end module mo_indprd
