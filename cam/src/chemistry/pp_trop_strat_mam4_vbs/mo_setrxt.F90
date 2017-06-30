
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

      rate(:,:,143) = 1.20e-10_r8
      rate(:,:,144) = 2.02e-10_r8
      rate(:,:,145) = 1.204e-10_r8
      rate(:,:,146) = 1.50e-10_r8
      rate(:,:,147) = 9.75e-11_r8
      rate(:,:,148) = 1.50e-11_r8
      rate(:,:,149) = 7.20e-11_r8
      rate(:,:,150) = 1.794e-10_r8
      rate(:,:,151) = 1.628e-10_r8
      rate(:,:,152) = 2.84e-10_r8
      rate(:,:,153) = 1.674e-10_r8
      rate(:,:,154) = 9.60e-11_r8
      rate(:,:,155) = 4.10e-11_r8
      rate(:,:,156) = 1.012e-10_r8
      rate(:,:,157) = 1.20e-10_r8
      rate(:,:,158) = 4.49e-10_r8
      rate(:,:,159) = 2.57e-10_r8
      rate(:,:,160) = 1.31e-10_r8
      rate(:,:,161) = 3.50e-11_r8
      rate(:,:,162) = 9.00e-12_r8
      rate(:,:,163) = 1.20e-10_r8
      rate(:,:,164) = 1.50e-10_r8
      rate(:,:,165) = 1.20e-10_r8
      rate(:,:,169) = 7.20e-11_r8
      rate(:,:,170) = 6.90e-12_r8
      rate(:,:,171) = 1.60e-12_r8
      rate(:,:,175) = 1.80e-12_r8
      rate(:,:,178) = 1.80e-12_r8
      rate(:,:,202) = 1.00e-11_r8
      rate(:,:,203) = 2.20e-11_r8
      rate(:,:,204) = 3.50e-12_r8
      rate(:,:,229) = 1.70e-13_r8
      rate(:,:,276) = 4.50e-13_r8
      rate(:,:,288) = 1.00e-14_r8
      rate(:,:,291) = 7.00e-13_r8
      rate(:,:,294) = 2.00e-13_r8
      rate(:,:,295) = 6.80e-14_r8
      rate(:,:,304) = 1.00e-12_r8
      rate(:,:,305) = 1.00e-11_r8
      rate(:,:,306) = 1.15e-11_r8
      rate(:,:,309) = 4.00e-14_r8
      rate(:,:,326) = 3.00e-12_r8
      rate(:,:,329) = 6.7e-13_r8
      rate(:,:,330) = 5.40e-11_r8
      rate(:,:,333) = 3.5e-13_r8
      rate(:,:,344) = 2.40e-12_r8
      rate(:,:,347) = 1.40e-11_r8
      rate(:,:,350) = 5.00e-12_r8
      rate(:,:,358) = 2.0e-12_r8
      rate(:,:,364) = 4.e-11_r8
      rate(:,:,365) = 4.e-11_r8
      rate(:,:,366) = 2.40e-12_r8
      rate(:,:,367) = 2.40e-12_r8
      rate(:,:,371) = 1.3e-11_r8
      rate(:,:,374) = 1.40e-11_r8
      rate(:,:,375) = 1.40e-11_r8
      rate(:,:,379) = 2.40e-12_r8
      rate(:,:,381) = 1.4e-11_r8
      rate(:,:,383) = 4.0e-11_r8
      rate(:,:,384) = 7.0e-11_r8
      rate(:,:,385) = 1.0e-10_r8
      rate(:,:,388) = 2.40e-12_r8
      rate(:,:,394) = 3.50e-12_r8
      rate(:,:,395) = 6.7e-12_r8
      rate(:,:,397) = 1.6e-12_r8
      rate(:,:,405) = 2.1e-12_r8
      rate(:,:,406) = 2.8e-13_r8
      rate(:,:,417) = 4.7e-11_r8
      rate(:,:,435) = 1.7e-11_r8
      rate(:,:,436) = 8.4e-11_r8
      rate(:,:,444) = 2.0e-10_r8
      rate(:,:,446) = 1.2e-14_r8
      rate(:,:,448) = 1.9e-11_r8
      rate(:,:,452) = 3.3e-11_r8
      rate(:,:,453) = 2.3e-11_r8
      rate(:,:,454) = 5.7e-11_r8
      rate(:,:,455) = 1e-12_r8
      rate(:,:,459) = 3.4e-11_r8
      rate(:,:,463) = 2.4e-12_r8
      rate(:,:,464) = 2.0e-11_r8
      rate(:,:,465) = 2.0e-11_r8
      rate(:,:,468) = 2.0e-10_r8
      rate(:,:,469) = 1.34e-11_r8
      rate(:,:,470) = 1.34e-11_r8
      rate(:,:,471) = 1.34e-11_r8
      rate(:,:,472) = 1.34e-11_r8
      rate(:,:,475) = 1.7e-11_r8
      rate(:,:,478) = 1.2e-14_r8
      rate(:,:,480) = 2.5e-12_r8
      rate(:,:,481) = 2.5e-12_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,136) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,138) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      rate(:,:,139) = 3.30e-11_r8 * exp( 55._r8 * itemp(:,:) )
      rate(:,:,140) = 1.63e-10_r8 * exp( 60._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 20._r8 * itemp(:,:) )
      rate(:,:,141) = 7.25e-11_r8 * exp_fac(:,:)
      rate(:,:,142) = 4.63e-11_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 100._r8 * itemp(:,:) )
      rate(:,:,166) = 7.70e-11_r8 * exp_fac(:,:)
      rate(:,:,187) = 2.10e-11_r8 * exp_fac(:,:)
      rate(:,:,168) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 180._r8 * itemp(:,:) )
      rate(:,:,172) = 1.80e-11_r8 * exp_fac(:,:)
      rate(:,:,286) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,313) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,318) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,337) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,362) = 4.40e-12_r8 * exp_fac(:,:)
      rate(:,:,363) = 4.40e-12_r8 * exp_fac(:,:)
      rate(:,:,449) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,456) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,460) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,173) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,174) = 4.80e-11_r8 * exp_fac(:,:)
      rate(:,:,239) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,177) = 2.80e-12_r8 * exp( -1800._r8 * itemp(:,:) )
      rate(:,:,179) = 1.60e-11_r8 * exp( -4570._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,180) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,247) = 5.50e-12_r8 * exp_fac(:,:)
      rate(:,:,275) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,296) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,316) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,320) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,325) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,339) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,348) = 2.30e-11_r8 * exp_fac(:,:)
      rate(:,:,370) = 1.52e-11_r8 * exp_fac(:,:)
      rate(:,:,399) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,404) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,409) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,412) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,419) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,427) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,439) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,441) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,181) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2000._r8 * itemp(:,:) )
      rate(:,:,183) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,360) = 1.05e-14_r8 * exp_fac(:,:)
      rate(:,:,476) = 1.05e-14_r8 * exp_fac(:,:)
      rate(:,:,185) = 7.80e-13_r8 * exp( -1050._r8 * itemp(:,:) )
      rate(:,:,186) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 220._r8 * itemp(:,:) )
      rate(:,:,188) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,189) = 1.45e-12_r8 * exp_fac(:,:)
      rate(:,:,190) = 1.45e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,192) = 3.30e-12_r8 * exp_fac(:,:)
      rate(:,:,211) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,216) = 7.40e-12_r8 * exp_fac(:,:)
      rate(:,:,299) = 8.10e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1500._r8 * itemp(:,:) )
      rate(:,:,193) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,248) = 5.80e-12_r8 * exp_fac(:,:)
      rate(:,:,194) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2450._r8 * itemp(:,:) )
      rate(:,:,196) = 1.20e-13_r8 * exp_fac(:,:)
      rate(:,:,222) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,201) = 1.50e-11_r8 * exp( 170._r8 * itemp(:,:) )
      rate(:,:,206) = 1.30e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,208) = 2.30e-11_r8 * exp( -200._r8 * itemp(:,:) )
      rate(:,:,209) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:,:) )
      rate(:,:,210) = 1.10e-11_r8 * exp( -980._r8 * itemp(:,:) )
      rate(:,:,212) = 3.60e-11_r8 * exp( -375._r8 * itemp(:,:) )
      rate(:,:,213) = 8.10e-11_r8 * exp( -30._r8 * itemp(:,:) )
      rate(:,:,214) = 7.30e-12_r8 * exp( -1280._r8 * itemp(:,:) )
      rate(:,:,215) = 2.80e-11_r8 * exp( 85._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 230._r8 * itemp(:,:) )
      rate(:,:,217) = 6.00e-13_r8 * exp_fac(:,:)
      rate(:,:,238) = 1.90e-11_r8 * exp_fac(:,:)
      rate(:,:,246) = 1.50e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 290._r8 * itemp(:,:) )
      rate(:,:,218) = 2.60e-12_r8 * exp_fac(:,:)
      rate(:,:,220) = 6.40e-12_r8 * exp_fac(:,:)
      rate(:,:,245) = 4.10e-13_r8 * exp_fac(:,:)
      rate(:,:,414) = 7.5e-12_r8 * exp_fac(:,:)
      rate(:,:,424) = 7.5e-12_r8 * exp_fac(:,:)
      rate(:,:,430) = 7.5e-12_r8 * exp_fac(:,:)
      rate(:,:,432) = 7.5e-12_r8 * exp_fac(:,:)
      rate(:,:,219) = 3.3e-12_r8 * exp( -115._r8 * itemp(:,:) )
      rate(:,:,223) = 1.00e-12_r8 * exp( -1590._r8 * itemp(:,:) )
      rate(:,:,224) = 3.50e-13_r8 * exp( -1370._r8 * itemp(:,:) )
      rate(:,:,227) = 1.80e-12_r8 * exp( -250._r8 * itemp(:,:) )
      rate(:,:,228) = 1.00e-11_r8 * exp( -3300._r8 * itemp(:,:) )
      rate(:,:,230) = 3.40e-12_r8 * exp( -130._r8 * itemp(:,:) )
      rate(:,:,231) = 3.00e-12_r8 * exp( -500._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -840._r8 * itemp(:,:) )
      rate(:,:,232) = 3.60e-12_r8 * exp_fac(:,:)
      rate(:,:,259) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,233) = 1.20e-12_r8 * exp( -330._r8 * itemp(:,:) )
      rate(:,:,234) = 6.50e-12_r8 * exp( 135._r8 * itemp(:,:) )
      rate(:,:,235) = 1.60e-11_r8 * exp( -780._r8 * itemp(:,:) )
      rate(:,:,236) = 4.80e-12_r8 * exp( -310._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -800._r8 * itemp(:,:) )
      rate(:,:,237) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,261) = 6.30e-12_r8 * exp_fac(:,:)
      rate(:,:,240) = 4.50e-12_r8 * exp( 460._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 260._r8 * itemp(:,:) )
      rate(:,:,241) = 8.80e-12_r8 * exp_fac(:,:)
      rate(:,:,244) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,243) = 9.50e-13_r8 * exp( 550._r8 * itemp(:,:) )
      rate(:,:,249) = 1.20e-10_r8 * exp( -430._r8 * itemp(:,:) )
      rate(:,:,250) = 1.90e-11_r8 * exp( 215._r8 * itemp(:,:) )
      rate(:,:,251) = 2.17e-11_r8 * exp( -1130._r8 * itemp(:,:) )
      rate(:,:,252) = 2.40e-12_r8 * exp( -1250._r8 * itemp(:,:) )
      rate(:,:,253) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1600._r8 * itemp(:,:) )
      rate(:,:,254) = 1.05e-12_r8 * exp_fac(:,:)
      rate(:,:,257) = 1.25e-12_r8 * exp_fac(:,:)
      rate(:,:,268) = 3.40e-11_r8 * exp_fac(:,:)
      rate(:,:,255) = 2.35e-12_r8 * exp( -1300._r8 * itemp(:,:) )
      rate(:,:,256) = 1.40e-11_r8 * exp( -1030._r8 * itemp(:,:) )
      rate(:,:,258) = 1.30e-12_r8 * exp( -1770._r8 * itemp(:,:) )
      rate(:,:,260) = 1.35e-12_r8 * exp( -600._r8 * itemp(:,:) )
      rate(:,:,262) = 4.85e-12_r8 * exp( -850._r8 * itemp(:,:) )
      rate(:,:,263) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,266) = 6.00e-13_r8 * exp( -2058._r8 * itemp(:,:) )
      rate(:,:,267) = 5.50e-12_r8 * exp( 125._r8 * itemp(:,:) )
      rate(:,:,269) = 9.7e-15_r8 * exp( 625._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 300._r8 * itemp(:,:) )
      rate(:,:,270) = 2.80e-12_r8 * exp_fac(:,:)
      rate(:,:,322) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,271) = 4.10e-13_r8 * exp( 750._r8 * itemp(:,:) )
      rate(:,:,272) = 5.00e-13_r8 * exp( -424._r8 * itemp(:,:) )
      rate(:,:,273) = 1.90e-14_r8 * exp( 706._r8 * itemp(:,:) )
      rate(:,:,274) = 2.90e-12_r8 * exp( -345._r8 * itemp(:,:) )
      rate(:,:,277) = 2.40e12_r8 * exp( -7000._r8 * itemp(:,:) )
      rate(:,:,278) = 2.60e-12_r8 * exp( 265._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 700._r8 * itemp(:,:) )
      rate(:,:,279) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,287) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,293) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,314) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,319) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,323) = 8.60e-13_r8 * exp_fac(:,:)
      rate(:,:,338) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,345) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,368) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,369) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,380) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,389) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,398) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,403) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,408) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,411) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,418) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,426) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,438) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,440) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,450) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,457) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,461) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,282) = 7.20e-11_r8 * exp( -70._r8 * itemp(:,:) )
      rate(:,:,284) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:,:) )
      rate(:,:,289) = 1.60e11_r8 * exp( -4150._r8 * itemp(:,:) )
      rate(:,:,290) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 365._r8 * itemp(:,:) )
      rate(:,:,292) = 2.60e-12_r8 * exp_fac(:,:)
      rate(:,:,402) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,407) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,410) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,420) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,428) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,437) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,442) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,297) = 4.63e-12_r8 * exp( 350._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1900._r8 * itemp(:,:) )
      rate(:,:,298) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,311) = 6.50e-15_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 1040._r8 * itemp(:,:) )
      rate(:,:,301) = 4.30e-13_r8 * exp_fac(:,:)
      rate(:,:,351) = 4.30e-13_r8 * exp_fac(:,:)
      rate(:,:,415) = 4.3E-13_r8 * exp_fac(:,:)
      rate(:,:,425) = 4.3E-13_r8 * exp_fac(:,:)
      rate(:,:,429) = 4.3E-13_r8 * exp_fac(:,:)
      rate(:,:,431) = 4.3E-13_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 500._r8 * itemp(:,:) )
      rate(:,:,302) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,303) = 2.50e-12_r8 * exp_fac(:,:)
      rate(:,:,324) = 7.10e-13_r8 * exp_fac(:,:)
      rate(:,:,352) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,451) = 2e-12_r8 * exp_fac(:,:)
      rate(:,:,458) = 2e-12_r8 * exp_fac(:,:)
      rate(:,:,462) = 2e-12_r8 * exp_fac(:,:)
      rate(:,:,307) = 6.90e-12_r8 * exp( -230._r8 * itemp(:,:) )
      rate(:,:,312) = 4.60e-13_r8 * exp( -1156._r8 * itemp(:,:) )
      rate(:,:,315) = 3.75e-13_r8 * exp( -40._r8 * itemp(:,:) )
      rate(:,:,317) = 8.70e-12_r8 * exp( -615._r8 * itemp(:,:) )
      rate(:,:,327) = 8.40e-13_r8 * exp( 830._r8 * itemp(:,:) )
      rate(:,:,328) = 1.40e-12_r8 * exp( -1860._r8 * itemp(:,:) )
      rate(:,:,331) = 4.8e-12_r8 * exp( 120._r8 * itemp(:,:) )
      rate(:,:,332) = 5.1e-14_r8 * exp( 693._r8 * itemp(:,:) )
      rate(:,:,334) = 4.13e-12_r8 * exp( 452._r8 * itemp(:,:) )
      rate(:,:,335) = 7.52e-16_r8 * exp( -1521._r8 * itemp(:,:) )
      rate(:,:,336) = 2.30e-12_r8 * exp( -170._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 175._r8 * itemp(:,:) )
      rate(:,:,340) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,377) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,386) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,341) = 4.40e-15_r8 * exp( -2500._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 360._r8 * itemp(:,:) )
      rate(:,:,342) = 2.70e-12_r8 * exp_fac(:,:)
      rate(:,:,343) = 1.30e-13_r8 * exp_fac(:,:)
      rate(:,:,349) = 5.30e-12_r8 * exp_fac(:,:)
      rate(:,:,378) = 2.70e-12_r8 * exp_fac(:,:)
      rate(:,:,387) = 2.7e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 400._r8 * itemp(:,:) )
      rate(:,:,346) = 5.00e-13_r8 * exp_fac(:,:)
      rate(:,:,372) = 5.00e-13_r8 * exp_fac(:,:)
      rate(:,:,373) = 5.00e-13_r8 * exp_fac(:,:)
      rate(:,:,382) = 5.e-13_r8 * exp_fac(:,:)
      rate(:,:,390) = 5.e-13_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 530._r8 * itemp(:,:) )
      rate(:,:,353) = 4.60e-12_r8 * exp_fac(:,:)
      rate(:,:,354) = 2.30e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 410._r8 * itemp(:,:) )
      rate(:,:,359) = 2.54e-11_r8 * exp_fac(:,:)
      rate(:,:,466) = 2.54e-11_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -446._r8 * itemp(:,:) )
      rate(:,:,361) = 3.03e-12_r8 * exp_fac(:,:)
      rate(:,:,479) = 3.03e-12_r8 * exp_fac(:,:)
      rate(:,:,376) = 1.6e9_r8 * exp( -8300._r8 * itemp(:,:) )
      rate(:,:,391) = 1.3e-12_r8 * exp( 640._r8 * itemp(:,:) )
      rate(:,:,392) = 1.90e-12_r8 * exp( 190._r8 * itemp(:,:) )
      rate(:,:,396) = 5.4e-14_r8 * exp( 870._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -193._r8 * itemp(:,:) )
      rate(:,:,400) = 2.3e-12_r8 * exp_fac(:,:)
      rate(:,:,473) = 2.3e-12_r8 * exp_fac(:,:)
      rate(:,:,401) = 4.7e-13_r8 * exp( 1220._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 352._r8 * itemp(:,:) )
      rate(:,:,416) = 1.7e-12_r8 * exp_fac(:,:)
      rate(:,:,474) = 1.7e-12_r8 * exp_fac(:,:)
      rate(:,:,421) = 5.9e-12_r8 * exp( 225._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 440._r8 * itemp(:,:) )
      rate(:,:,443) = 1.2e-11_r8 * exp_fac(:,:)
      rate(:,:,467) = 1.2e-11_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -580._r8 * itemp(:,:) )
      rate(:,:,445) = 6.3e-16_r8 * exp_fac(:,:)
      rate(:,:,477) = 6.3e-16_r8 * exp_fac(:,:)
      rate(:,:,447) = 1.2e-12_r8 * exp( 490._r8 * itemp(:,:) )
      rate(:,:,488) = 9.60e-12_r8 * exp( -234._r8 * itemp(:,:) )
      rate(:,:,490) = 1.90e-13_r8 * exp( 520._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( rate(1,1,167), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.90e-31_r8 * itemp(:,:)**1.0_r8
      kinf(:,:) = 2.60e-11_r8
      call jpl( rate(1,1,176), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 4.28e-33_r8
      kinf(:,:) = 9.30e-15_r8 * itemp(:,:)**(-4.42_r8)
      call jpl( rate(1,1,184), m, 0.8_r8, ko, kinf, n )

      ko(:,:) = 9.00e-32_r8 * itemp(:,:)**1.5_r8
      kinf(:,:) = 3.0e-11_r8
      call jpl( rate(1,1,191), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.50e-31_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 2.2e-11_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,195), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-30_r8 * itemp(:,:)**4.4_r8
      kinf(:,:) = 1.4e-12_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,197), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,199), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 2.9e-12_r8 * itemp(:,:)**1.1_r8
      call jpl( rate(1,1,205), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 1.5e-11_r8 * itemp(:,:)**1.9_r8
      call jpl( rate(1,1,221), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.60e-32_r8 * itemp(:,:)**4.5_r8
      kinf(:,:) = 3.0e-12_r8 * itemp(:,:)**2.0_r8
      call jpl( rate(1,1,225), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.20e-31_r8 * itemp(:,:)**3.2_r8
      kinf(:,:) = 6.9e-12_r8 * itemp(:,:)**2.9_r8
      call jpl( rate(1,1,242), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.90e-33_r8 * itemp(:,:)**1.4_r8
      kinf(:,:) = 1.10e-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,265), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.20e-30_r8 * itemp(:,:)**2.4_r8
      kinf(:,:) = 2.2e-10_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,280), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.60e-29_r8 * itemp(:,:)**3.3_r8
      kinf(:,:) = 3.1e-10_r8 * itemp(:,:)
      call jpl( rate(1,1,281), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.50e-30_r8
      kinf(:,:) = 8.3e-13_r8 * itemp(:,:)**(-2.0_r8)
      call jpl( rate(1,1,283), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 8.60e-29_r8 * itemp(:,:)**3.1_r8
      kinf(:,:) = 9.00e-12_r8 * itemp(:,:)**0.85_r8
      call jpl( rate(1,1,285), m, 0.48_r8, ko, kinf, n )

      ko(:,:) = 9.70e-29_r8 * itemp(:,:)**5.6_r8
      kinf(:,:) = 9.30e-12_r8 * itemp(:,:)**1.5_r8
      call jpl( rate(1,1,300), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 8.00e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3.00e-11_r8
      call jpl( rate(1,1,310), m, 0.5_r8, ko, kinf, n )

      ko(:,:) = 8.00e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3.00e-11_r8
      call jpl( rate(1,1,357), m, 0.5_r8, ko, kinf, n )

      ko(:,:) = 9.70e-29_r8 * itemp(:,:)**5.6_r8
      kinf(:,:) = 9.30e-12_r8 * itemp(:,:)**1.5_r8
      call jpl( rate(1,1,413), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 9.70e-29_r8 * itemp(:,:)**5.6_r8
      kinf(:,:) = 9.30e-12_r8 * itemp(:,:)**1.5_r8
      call jpl( rate(1,1,422), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 9.70e-29_r8 * itemp(:,:)**5.6_r8
      kinf(:,:) = 9.30e-12_r8 * itemp(:,:)**1.5_r8
      call jpl( rate(1,1,433), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 9.70e-29_r8 * itemp(:,:)**5.6_r8
      kinf(:,:) = 9.30e-12_r8 * itemp(:,:)**1.5_r8
      call jpl( rate(1,1,434), m, 0.6_r8, ko, kinf, n )

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


      end subroutine setrxt_hrates

      end module mo_setrxt
