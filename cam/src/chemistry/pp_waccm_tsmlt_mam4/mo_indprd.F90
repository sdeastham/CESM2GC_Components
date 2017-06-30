      module mo_indprd
      use shr_kind_mod, only : r8 => shr_kind_r8
      private
      public :: indprd
      contains
      subroutine indprd( class, prod, nprod, y, extfrc, rxt, chnkpnts )
      use chem_mods, only : gas_pcnst, extcnt, rxntot
      implicit none
!--------------------------------------------------------------------
! ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: class
      integer, intent(in) :: chnkpnts
      integer, intent(in) :: nprod
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(in) :: extfrc(chnkpnts,extcnt)
      real(r8), intent(inout) :: prod(chnkpnts,nprod)
!--------------------------------------------------------------------
! ... "independent" production for Explicit species
!--------------------------------------------------------------------
      if( class == 1 ) then
         prod(:,1) = 0._r8
         prod(:,2) = 0._r8
         prod(:,3) = 0._r8
         prod(:,4) = 0._r8
         prod(:,5) = 0._r8
         prod(:,6) = 0._r8
         prod(:,7) = 0._r8
         prod(:,8) = 0._r8
         prod(:,9) = 0._r8
         prod(:,10) = 0._r8
         prod(:,11) = 0._r8
         prod(:,12) = 0._r8
         prod(:,13) =.100_r8*rxt(:,349)*y(:,118)*y(:,25)
         prod(:,14) = 0._r8
         prod(:,15) = (rxt(:,309)*y(:,57) +rxt(:,311)*y(:,78) +rxt(:,316)*y(:,57) + &
                 rxt(:,336)*y(:,46) +.500_r8*rxt(:,337)*y(:,47) + &
                 .800_r8*rxt(:,342)*y(:,65) +rxt(:,343)*y(:,66) + &
                 .500_r8*rxt(:,392)*y(:,100) +1.800_r8*rxt(:,502)*y(:,150))*y(:,216) &
                  + (2.000_r8*rxt(:,332)*y(:,191) +.900_r8*rxt(:,333)*y(:,192) + &
                 rxt(:,335)*y(:,111) +2.000_r8*rxt(:,382)*y(:,204) + &
                 rxt(:,406)*y(:,200) +rxt(:,431)*y(:,224))*y(:,191) &
                  + (.200_r8*rxt(:,349)*y(:,25) +.100_r8*rxt(:,393)*y(:,102) + &
                 .270_r8*rxt(:,481)*y(:,3) +.270_r8*rxt(:,484)*y(:,101))*y(:,118) &
                  + (rxt(:,383)*y(:,192) +.450_r8*rxt(:,384)*y(:,198) + &
                 2.000_r8*rxt(:,385)*y(:,204))*y(:,204) &
                  + (.500_r8*rxt(:,491)*y(:,192) +.900_r8*rxt(:,493)*y(:,111)) &
                 *y(:,221) +rxt(:,38)*y(:,47) +.400_r8*rxt(:,61)*y(:,123) +rxt(:,66) &
                 *y(:,146) +.800_r8*rxt(:,70)*y(:,150)
         prod(:,16) = 0._r8
         prod(:,17) = 0._r8
         prod(:,18) = 0._r8
         prod(:,19) = 0._r8
         prod(:,20) =rxt(:,195)*y(:,112)*y(:,103)
         prod(:,21) = 0._r8
         prod(:,22) =rxt(:,519)*y(:,216)*y(:,108) +rxt(:,525)*y(:,109)
         prod(:,23) = (rxt(:,453)*y(:,193) +rxt(:,456)*y(:,203) +rxt(:,459)*y(:,205) + &
                 rxt(:,463)*y(:,125))*y(:,112) +.500_r8*rxt(:,392)*y(:,216)*y(:,100) &
                  +rxt(:,524)*y(:,110) +.200_r8*rxt(:,488)*y(:,211)*y(:,111) &
                  +.500_r8*rxt(:,500)*y(:,149)*y(:,113)
         prod(:,24) = 0._r8
         prod(:,25) = 0._r8
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,125) = 0._r8
         prod(:,124) = 0._r8
         prod(:,157) = 0._r8
         prod(:,49) = 0._r8
         prod(:,88) = 0._r8
         prod(:,50) = 0._r8
         prod(:,89) = 0._r8
         prod(:,96) = 0._r8
         prod(:,70) = 0._r8
         prod(:,118) = 0._r8
         prod(:,80) = 0._r8
         prod(:,62) = 0._r8
         prod(:,85) = 0._r8
         prod(:,185) =rxt(:,80)*y(:,30) +rxt(:,81)*y(:,31) +2.000_r8*rxt(:,87)*y(:,37) &
                  +rxt(:,88)*y(:,39) +3.000_r8*rxt(:,91)*y(:,51) +2.000_r8*rxt(:,99) &
                 *y(:,69)
         prod(:,64) = 0._r8
         prod(:,189) = 0._r8
         prod(:,112) = 0._r8
         prod(:,65) = 0._r8
         prod(:,81) = 0._r8
         prod(:,71) = 0._r8
         prod(:,115) = 0._r8
         prod(:,67) = 0._r8
         prod(:,82) = 0._r8
         prod(:,72) = 0._r8
         prod(:,161) = 0._r8
         prod(:,91) = 0._r8
         prod(:,40) = 0._r8
         prod(:,68) = 0._r8
         prod(:,201) =.180_r8*rxt(:,41)*y(:,50)
         prod(:,171) = 0._r8
         prod(:,39) = 0._r8
         prod(:,160) = 0._r8
         prod(:,176) = 0._r8
         prod(:,116) = 0._r8
         prod(:,110) = 0._r8
         prod(:,142) = 0._r8
         prod(:,92) = 0._r8
         prod(:,196) =4.000_r8*rxt(:,79)*y(:,29) +rxt(:,80)*y(:,30) &
                  +2.000_r8*rxt(:,82)*y(:,32) +2.000_r8*rxt(:,83)*y(:,33) &
                  +2.000_r8*rxt(:,84)*y(:,34) +rxt(:,85)*y(:,35) +2.000_r8*rxt(:,86) &
                 *y(:,36) +3.000_r8*rxt(:,89)*y(:,40) +rxt(:,90)*y(:,42) +rxt(:,101) &
                 *y(:,73) +rxt(:,102)*y(:,74) +rxt(:,103)*y(:,75)
         prod(:,52) = 0._r8
         prod(:,37) = 0._r8
         prod(:,192) = 0._r8
         prod(:,156) = 0._r8
         prod(:,168) = (rxt(:,42) +rxt(:,109))*y(:,58) +.380_r8*rxt(:,41)*y(:,50) &
                  + extfrc(:,9)
         prod(:,41) =rxt(:,80)*y(:,30) +rxt(:,81)*y(:,31) +rxt(:,83)*y(:,33) &
                  +2.000_r8*rxt(:,84)*y(:,34) +2.000_r8*rxt(:,85)*y(:,35) +rxt(:,86) &
                 *y(:,36) +2.000_r8*rxt(:,99)*y(:,69) +rxt(:,102)*y(:,74) +rxt(:,103) &
                 *y(:,75)
         prod(:,54) =rxt(:,82)*y(:,32) +rxt(:,83)*y(:,33) +rxt(:,101)*y(:,73)
         prod(:,56) = 0._r8
         prod(:,73) = 0._r8
         prod(:,42) = 0._r8
         prod(:,136) =rxt(:,81)*y(:,31) +rxt(:,85)*y(:,35)
         prod(:,165) = 0._r8
         prod(:,153) = 0._r8
         prod(:,195) = (rxt(:,40) +.330_r8*rxt(:,41))*y(:,50)
         prod(:,172) =1.440_r8*rxt(:,41)*y(:,50)
         prod(:,120) = 0._r8
         prod(:,43) = 0._r8
         prod(:,147) = 0._r8
         prod(:,184) = 0._r8
         prod(:,53) = 0._r8
         prod(:,143) = 0._r8
         prod(:,63) = 0._r8
         prod(:,194) = 0._r8
         prod(:,100) = 0._r8
         prod(:,135) = 0._r8
         prod(:,148) = 0._r8
         prod(:,167) = 0._r8
         prod(:,61) = 0._r8
         prod(:,169) = 0._r8
         prod(:,83) = 0._r8
         prod(:,44) = 0._r8
         prod(:,149) = 0._r8
         prod(:,114) = 0._r8
         prod(:,102) = 0._r8
         prod(:,174) = 0._r8
         prod(:,86) = 0._r8
         prod(:,132) = 0._r8
         prod(:,35) = 0._r8
         prod(:,175) = 0._r8
         prod(:,74) = 0._r8
         prod(:,106) = 0._r8
         prod(:,75) = 0._r8
         prod(:,113) = 0._r8
         prod(:,155) = 0._r8
         prod(:,180) = 0._r8
         prod(:,145) = (.800_r8*rxt(:,111) +rxt(:,114) +rxt(:,115) + &
                 .800_r8*rxt(:,117)) + extfrc(:,15)
         prod(:,69) = 0._r8
         prod(:,76) = 0._r8
         prod(:,146) = 0._r8
         prod(:,38) = 0._r8
         prod(:,1) = 0._r8
         prod(:,2) = 0._r8
         prod(:,188) = + extfrc(:,8)
         prod(:,200) = + extfrc(:,11)
         prod(:,198) = 0._r8
         prod(:,150) = 0._r8
         prod(:,77) = 0._r8
         prod(:,190) = (rxt(:,42) +rxt(:,109))*y(:,58) +.180_r8*rxt(:,41)*y(:,50)
         prod(:,187) = 0._r8
         prod(:,197) = 0._r8
         prod(:,78) = 0._r8
         prod(:,84) = 0._r8
         prod(:,45) = 0._r8
         prod(:,98) = 0._r8
         prod(:,46) = 0._r8
         prod(:,99) = 0._r8
         prod(:,51) = 0._r8
         prod(:,79) = 0._r8
         prod(:,111) = 0._r8
         prod(:,87) = 0._r8
         prod(:,131) = 0._r8
         prod(:,183) = 0._r8
         prod(:,154) = + extfrc(:,10)
         prod(:,66) = 0._r8
         prod(:,3) = 0._r8
         prod(:,4) = 0._r8
         prod(:,5) = 0._r8
         prod(:,6) = 0._r8
         prod(:,7) = 0._r8
         prod(:,8) = 0._r8
         prod(:,9) = 0._r8
         prod(:,10) = 0._r8
         prod(:,11) = 0._r8
         prod(:,12) = 0._r8
         prod(:,36) = + extfrc(:,12)
         prod(:,57) = 0._r8
         prod(:,119) = 0._r8
         prod(:,107) = 0._r8
         prod(:,101) = 0._r8
         prod(:,158) = 0._r8
         prod(:,163) = 0._r8
         prod(:,123) = 0._r8
         prod(:,55) = 0._r8
         prod(:,58) = 0._r8
         prod(:,59) = 0._r8
         prod(:,128) = 0._r8
         prod(:,60) = 0._r8
         prod(:,90) = 0._r8
         prod(:,13) = + extfrc(:,4)
         prod(:,14) = + extfrc(:,7)
         prod(:,15) = 0._r8
         prod(:,16) = 0._r8
         prod(:,17) = 0._r8
         prod(:,18) = 0._r8
         prod(:,19) = 0._r8
         prod(:,20) = 0._r8
         prod(:,21) = + extfrc(:,5)
         prod(:,22) = + extfrc(:,6)
         prod(:,23) = 0._r8
         prod(:,24) = + extfrc(:,14)
         prod(:,25) = + extfrc(:,3)
         prod(:,26) = + extfrc(:,13)
         prod(:,27) = + extfrc(:,1)
         prod(:,28) = + extfrc(:,2)
         prod(:,29) = 0._r8
         prod(:,30) = 0._r8
         prod(:,31) = 0._r8
         prod(:,32) = 0._r8
         prod(:,33) = 0._r8
         prod(:,34) = 0._r8
         prod(:,108) = 0._r8
         prod(:,159) = 0._r8
         prod(:,103) = 0._r8
         prod(:,93) = 0._r8
         prod(:,151) = 0._r8
         prod(:,144) = 0._r8
         prod(:,122) = 0._r8
         prod(:,182) = 0._r8
         prod(:,186) =rxt(:,88)*y(:,39) +rxt(:,90)*y(:,42) +rxt(:,40)*y(:,50)
         prod(:,134) = 0._r8
         prod(:,140) = (rxt(:,112) +rxt(:,113) +rxt(:,114) +rxt(:,115) +rxt(:,116) + &
                 rxt(:,118)) + extfrc(:,22)
         prod(:,117) = 0._r8
         prod(:,97) = 0._r8
         prod(:,137) = 0._r8
         prod(:,199) = 0._r8
         prod(:,94) = 0._r8
         prod(:,177) = 0._r8
         prod(:,178) = 0._r8
         prod(:,179) = 0._r8
         prod(:,129) = 0._r8
         prod(:,181) = 0._r8
         prod(:,164) = 0._r8
         prod(:,126) = 0._r8
         prod(:,109) = (1.200_r8*rxt(:,111) +rxt(:,112) +rxt(:,116) + &
                 1.200_r8*rxt(:,117)) + extfrc(:,21)
         prod(:,127) = (rxt(:,113) +rxt(:,118)) + extfrc(:,19)
         prod(:,139) = 0._r8
         prod(:,104) = (rxt(:,112) +rxt(:,114) +rxt(:,115) +rxt(:,116)) + extfrc(:,20)
         prod(:,162) = 0._r8
         prod(:,193) =rxt(:,12)*y(:,104)
         prod(:,47) = 0._r8
         prod(:,48) = 0._r8
         prod(:,138) = + extfrc(:,18)
         prod(:,191) =.330_r8*rxt(:,41)*y(:,50) + extfrc(:,16)
         prod(:,121) = + extfrc(:,17)
         prod(:,95) = 0._r8
         prod(:,141) = 0._r8
         prod(:,170) = 0._r8
         prod(:,166) = 0._r8
         prod(:,152) = 0._r8
         prod(:,130) = 0._r8
         prod(:,173) = 0._r8
         prod(:,133) = 0._r8
         prod(:,105) = 0._r8
         prod(:,202) =.050_r8*rxt(:,41)*y(:,50)
      end if
      end subroutine indprd
      end module mo_indprd
