
      module mo_setrxt

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setrxt
      public :: setrxt_hrates

      contains

      subroutine setrxt( rate, temp, m, ncol )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      real(r8)  ::  itemp(ncol,pver)
      real(r8)  ::  exp_fac(ncol,pver)
      real(r8)  :: ko(ncol,pver)
      real(r8)  :: kinf(ncol,pver)

      rate(:,:,90) = 8.00e-14_r8
      rate(:,:,91) = 3.90e-17_r8
      rate(:,:,94) = 4.20e-13_r8
      rate(:,:,95) = 8.50e-2_r8
      rate(:,:,96) = 1.30e-16_r8
      rate(:,:,98) = 1.00e-20_r8
      rate(:,:,99) = 2.58e-04_r8
      rate(:,:,106) = 1.20e-10_r8
      rate(:,:,107) = 2.02e-10_r8
      rate(:,:,108) = 1.204e-10_r8
      rate(:,:,109) = 1.50e-10_r8
      rate(:,:,110) = 9.75e-11_r8
      rate(:,:,111) = 1.50e-11_r8
      rate(:,:,112) = 7.20e-11_r8
      rate(:,:,113) = 1.794e-10_r8
      rate(:,:,114) = 1.628e-10_r8
      rate(:,:,115) = 2.84e-10_r8
      rate(:,:,116) = 1.674e-10_r8
      rate(:,:,117) = 9.60e-11_r8
      rate(:,:,118) = 4.10e-11_r8
      rate(:,:,119) = 1.012e-10_r8
      rate(:,:,120) = 1.20e-10_r8
      rate(:,:,121) = 4.49e-10_r8
      rate(:,:,122) = 2.57e-10_r8
      rate(:,:,123) = 2.14e-11_r8
      rate(:,:,124) = 1.90e-10_r8
      rate(:,:,125) = 1.31e-10_r8
      rate(:,:,126) = 3.50e-11_r8
      rate(:,:,127) = 9.00e-12_r8
      rate(:,:,128) = 1.20e-10_r8
      rate(:,:,129) = 1.50e-10_r8
      rate(:,:,130) = 1.20e-10_r8
      rate(:,:,133) = 7.20e-11_r8
      rate(:,:,134) = 6.90e-12_r8
      rate(:,:,135) = 1.60e-12_r8
      rate(:,:,139) = 1.80e-12_r8
      rate(:,:,142) = 1.80e-12_r8
      rate(:,:,148) = 5.00e-12_r8
      rate(:,:,149) = 7.00e-13_r8
      rate(:,:,150) = 5.00e-11_r8
      rate(:,:,167) = 1.00e-11_r8
      rate(:,:,168) = 2.20e-11_r8
      rate(:,:,169) = 3.50e-12_r8
      rate(:,:,194) = 1.70e-13_r8
      rate(:,:,266) = 9.0e-10_r8
      rate(:,:,267) = 1.0e-10_r8
      rate(:,:,268) = 4.4e-10_r8
      rate(:,:,269) = 4.0e-10_r8
      rate(:,:,270) = 2.0e-10_r8
      rate(:,:,271) = 1.0e-12_r8
      rate(:,:,272) = 6.0e-11_r8
      rate(:,:,273) = 5.0e-16_r8
      rate(:,:,277) = 4.8e-10_r8
      rate(:,:,278) = 1.0e-10_r8
      rate(:,:,279) = 4.0e-10_r8
      rate(:,:,282) = 5.0e-12_r8
      rate(:,:,283) = 7.0e-10_r8
      rate(:,:,284) = 8.0e-10_r8
      rate(:,:,286) = 4.7e-2_r8
      rate(:,:,287) = 1.71e-1_r8
      rate(:,:,288) = 7.7e-5_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,88) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,92) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:,93) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:,97) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:,100) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:,101) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:,102) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:,103) = 1.63e-10_r8 * exp( 60._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 20._r8 * itemp(:,:) )
      rate(:,:,104) = 7.25e-11_r8 * exp_fac(:,:)
      rate(:,:,105) = 4.63e-11_r8 * exp_fac(:,:)
      rate(:,:,132) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:,136) = 1.80e-11_r8 * exp( 180._r8 * itemp(:,:) )
      rate(:,:,137) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,138) = 4.80e-11_r8 * exp_fac(:,:)
      rate(:,:,204) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,141) = 2.80e-12_r8 * exp( -1800._r8 * itemp(:,:) )
      rate(:,:,143) = 1.60e-11_r8 * exp( -4570._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,144) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,212) = 5.50e-12_r8 * exp_fac(:,:)
      rate(:,:,240) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,145) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:,147) = 1.40e-12_r8 * exp( -2000._r8 * itemp(:,:) )
      rate(:,:,151) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:,152) = 2.10e-11_r8 * exp( 100._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 220._r8 * itemp(:,:) )
      rate(:,:,153) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,154) = 1.45e-12_r8 * exp_fac(:,:)
      rate(:,:,155) = 1.45e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,157) = 3.30e-12_r8 * exp_fac(:,:)
      rate(:,:,176) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,181) = 7.40e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1500._r8 * itemp(:,:) )
      rate(:,:,158) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,213) = 5.80e-12_r8 * exp_fac(:,:)
      rate(:,:,159) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2450._r8 * itemp(:,:) )
      rate(:,:,161) = 1.20e-13_r8 * exp_fac(:,:)
      rate(:,:,187) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,166) = 1.50e-11_r8 * exp( 170._r8 * itemp(:,:) )
      rate(:,:,171) = 1.30e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,173) = 2.30e-11_r8 * exp( -200._r8 * itemp(:,:) )
      rate(:,:,174) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:,:) )
      rate(:,:,175) = 1.10e-11_r8 * exp( -980._r8 * itemp(:,:) )
      rate(:,:,177) = 3.60e-11_r8 * exp( -375._r8 * itemp(:,:) )
      rate(:,:,178) = 8.10e-11_r8 * exp( -30._r8 * itemp(:,:) )
      rate(:,:,179) = 7.30e-12_r8 * exp( -1280._r8 * itemp(:,:) )
      rate(:,:,180) = 2.80e-11_r8 * exp( 85._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 230._r8 * itemp(:,:) )
      rate(:,:,182) = 6.00e-13_r8 * exp_fac(:,:)
      rate(:,:,203) = 1.90e-11_r8 * exp_fac(:,:)
      rate(:,:,211) = 1.50e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 290._r8 * itemp(:,:) )
      rate(:,:,183) = 2.60e-12_r8 * exp_fac(:,:)
      rate(:,:,185) = 6.40e-12_r8 * exp_fac(:,:)
      rate(:,:,210) = 4.10e-13_r8 * exp_fac(:,:)
      rate(:,:,184) = 3.3e-12_r8 * exp( -115._r8 * itemp(:,:) )
      rate(:,:,188) = 1.00e-12_r8 * exp( -1590._r8 * itemp(:,:) )
      rate(:,:,189) = 3.50e-13_r8 * exp( -1370._r8 * itemp(:,:) )
      rate(:,:,192) = 1.80e-12_r8 * exp( -250._r8 * itemp(:,:) )
      rate(:,:,193) = 1.00e-11_r8 * exp( -3300._r8 * itemp(:,:) )
      rate(:,:,195) = 3.40e-12_r8 * exp( -130._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -500._r8 * itemp(:,:) )
      rate(:,:,196) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,217) = 1.40e-10_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -840._r8 * itemp(:,:) )
      rate(:,:,197) = 3.60e-12_r8 * exp_fac(:,:)
      rate(:,:,228) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,198) = 1.20e-12_r8 * exp( -330._r8 * itemp(:,:) )
      rate(:,:,199) = 6.50e-12_r8 * exp( 135._r8 * itemp(:,:) )
      rate(:,:,200) = 1.60e-11_r8 * exp( -780._r8 * itemp(:,:) )
      rate(:,:,201) = 4.80e-12_r8 * exp( -310._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -800._r8 * itemp(:,:) )
      rate(:,:,202) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,230) = 6.30e-12_r8 * exp_fac(:,:)
      rate(:,:,205) = 4.50e-12_r8 * exp( 460._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 260._r8 * itemp(:,:) )
      rate(:,:,206) = 8.80e-12_r8 * exp_fac(:,:)
      rate(:,:,209) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,208) = 9.50e-13_r8 * exp( 550._r8 * itemp(:,:) )
      rate(:,:,214) = 1.20e-10_r8 * exp( -430._r8 * itemp(:,:) )
      rate(:,:,215) = 1.90e-11_r8 * exp( 215._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 0._r8 * itemp(:,:) )
      rate(:,:,216) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,266) = 9.0e-10_r8 * exp_fac(:,:)
      rate(:,:,267) = 1.0e-10_r8 * exp_fac(:,:)
      rate(:,:,268) = 4.4e-10_r8 * exp_fac(:,:)
      rate(:,:,269) = 4.0e-10_r8 * exp_fac(:,:)
      rate(:,:,270) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,271) = 1.0e-12_r8 * exp_fac(:,:)
      rate(:,:,272) = 6.0e-11_r8 * exp_fac(:,:)
      rate(:,:,273) = 5.0e-16_r8 * exp_fac(:,:)
      rate(:,:,277) = 4.8e-10_r8 * exp_fac(:,:)
      rate(:,:,278) = 1.0e-10_r8 * exp_fac(:,:)
      rate(:,:,279) = 4.0e-10_r8 * exp_fac(:,:)
      rate(:,:,282) = 5.0e-12_r8 * exp_fac(:,:)
      rate(:,:,283) = 7.0e-10_r8 * exp_fac(:,:)
      rate(:,:,284) = 8.0e-10_r8 * exp_fac(:,:)
      rate(:,:,286) = 4.7e-2_r8 * exp_fac(:,:)
      rate(:,:,287) = 1.71e-1_r8 * exp_fac(:,:)
      rate(:,:,288) = 7.7e-5_r8 * exp_fac(:,:)
      rate(:,:,218) = 1.60e-10_r8 * exp( -260._r8 * itemp(:,:) )
      rate(:,:,219) = 6.00e-12_r8 * exp( 400._r8 * itemp(:,:) )
      rate(:,:,220) = 2.17e-11_r8 * exp( -1130._r8 * itemp(:,:) )
      rate(:,:,221) = 2.40e-12_r8 * exp( -1250._r8 * itemp(:,:) )
      rate(:,:,222) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1600._r8 * itemp(:,:) )
      rate(:,:,223) = 1.05e-12_r8 * exp_fac(:,:)
      rate(:,:,226) = 1.25e-12_r8 * exp_fac(:,:)
      rate(:,:,237) = 3.40e-11_r8 * exp_fac(:,:)
      rate(:,:,224) = 2.35e-12_r8 * exp( -1300._r8 * itemp(:,:) )
      rate(:,:,225) = 1.40e-11_r8 * exp( -1030._r8 * itemp(:,:) )
      rate(:,:,227) = 1.30e-12_r8 * exp( -1770._r8 * itemp(:,:) )
      rate(:,:,229) = 1.35e-12_r8 * exp( -600._r8 * itemp(:,:) )
      rate(:,:,231) = 4.85e-12_r8 * exp( -850._r8 * itemp(:,:) )
      rate(:,:,232) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,235) = 6.00e-13_r8 * exp( -2058._r8 * itemp(:,:) )
      rate(:,:,236) = 5.50e-12_r8 * exp( 125._r8 * itemp(:,:) )
      rate(:,:,238) = 2.80e-12_r8 * exp( 300._r8 * itemp(:,:) )
      rate(:,:,239) = 4.10e-13_r8 * exp( 750._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( rate(1,1,131), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.90e-31_r8 * itemp(:,:)**1.0_r8
      kinf(:,:) = 2.60e-11_r8
      call jpl( rate(1,1,140), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 9.00e-32_r8 * itemp(:,:)**1.5_r8
      kinf(:,:) = 3.0e-11_r8
      call jpl( rate(1,1,156), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.50e-31_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 2.2e-11_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,160), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-30_r8 * itemp(:,:)**4.4_r8
      kinf(:,:) = 1.4e-12_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,162), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,164), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 2.9e-12_r8 * itemp(:,:)**1.1_r8
      call jpl( rate(1,1,170), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 1.5e-11_r8 * itemp(:,:)**1.9_r8
      call jpl( rate(1,1,186), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.60e-32_r8 * itemp(:,:)**4.5_r8
      kinf(:,:) = 3.0e-12_r8 * itemp(:,:)**2.0_r8
      call jpl( rate(1,1,190), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.20e-31_r8 * itemp(:,:)**3.2_r8
      kinf(:,:) = 6.9e-12_r8 * itemp(:,:)**2.9_r8
      call jpl( rate(1,1,207), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.90e-33_r8 * itemp(:,:)**1.4_r8
      kinf(:,:) = 1.10e-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,234), m, 0.6_r8, ko, kinf, n )

      end subroutine setrxt


      subroutine setrxt_hrates( rate, temp, m, ncol, kbot )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: kbot
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      real(r8)  ::  itemp(ncol,kbot)
      real(r8)  ::  exp_fac(ncol,kbot)
      real(r8)  :: ko(ncol,kbot)
      real(r8)  :: kinf(ncol,kbot)
      real(r8)  :: wrk(ncol,kbot)

      rate(:,:kbot,90) = 8.00e-14_r8
      rate(:,:kbot,91) = 3.90e-17_r8
      rate(:,:kbot,96) = 1.30e-16_r8
      rate(:,:kbot,98) = 1.00e-20_r8
      rate(:,:kbot,134) = 6.90e-12_r8
      rate(:,:kbot,148) = 5.00e-12_r8
      rate(:,:kbot,149) = 7.00e-13_r8
      rate(:,:kbot,267) = 1.0e-10_r8
      rate(:,:kbot,268) = 4.4e-10_r8
      rate(:,:kbot,269) = 4.0e-10_r8
      rate(:,:kbot,270) = 2.0e-10_r8
      rate(:,:kbot,271) = 1.0e-12_r8
      rate(:,:kbot,272) = 6.0e-11_r8
      rate(:,:kbot,277) = 4.8e-10_r8
      rate(:,:kbot,278) = 1.0e-10_r8
      rate(:,:kbot,279) = 4.0e-10_r8
      rate(:,:kbot,282) = 5.0e-12_r8
      rate(:,:kbot,283) = 7.0e-10_r8
      rate(:,:kbot,284) = 8.0e-10_r8
      rate(:,:kbot,286) = 4.7e-2_r8
      rate(:,:kbot,287) = 1.71e-1_r8
      rate(:,:kbot,288) = 7.7e-5_r8
      itemp(:ncol,:kbot) = 1._r8 / temp(:ncol,:kbot)
      n = ncol*kbot
      rate(:,:kbot,88) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:kbot,92) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:kbot,93) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:kbot,97) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:kbot,100) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:kbot,101) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:kbot,102) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:kbot,132) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:kbot,136) = 1.80e-11_r8 * exp( 180._r8 * itemp(:,:) )
      rate(:,:kbot,137) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:kbot,138) = 4.80e-11_r8 * exp( 250._r8 * itemp(:,:) )
      rate(:,:kbot,144) = 3.00e-11_r8 * exp( 200._r8 * itemp(:,:) )
      rate(:,:kbot,145) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:kbot,151) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:kbot,152) = 2.10e-11_r8 * exp( 100._r8 * itemp(:,:) )
      rate(:,:kbot,157) = 3.30e-12_r8 * exp( 270._r8 * itemp(:,:) )
      rate(:,:kbot,158) = 3.00e-12_r8 * exp( -1500._r8 * itemp(:,:) )
      rate(:,:kbot,159) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:,:kbot,131) = wrk(:,:)











      end subroutine setrxt_hrates

      end module mo_setrxt
