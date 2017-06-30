      module mo_nln_matrix
      use shr_kind_mod, only : r8 => shr_kind_r8
      private
      public :: nlnmat
      contains
      subroutine nlnmat01( ofl, ofu, chnkpnts, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = ofl,ofu
         mat(k,544) = -(rxt(k,397)*y(k,216))
         mat(k,1645) = -rxt(k,397)*y(k,1)
         mat(k,1445) = rxt(k,400)*y(k,185)
         mat(k,919) = rxt(k,400)*y(k,111)
         mat(k,533) = -(rxt(k,401)*y(k,216))
         mat(k,1644) = -rxt(k,401)*y(k,2)
         mat(k,918) = rxt(k,398)*y(k,198)
         mat(k,1999) = rxt(k,398)*y(k,185)
         mat(k,890) = -(rxt(k,480)*y(k,113) + rxt(k,481)*y(k,118) + rxt(k,482) &
                      *y(k,216))
         mat(k,1922) = -rxt(k,480)*y(k,3)
         mat(k,1868) = -rxt(k,481)*y(k,3)
         mat(k,1672) = -rxt(k,482)*y(k,3)
         mat(k,87) = -(rxt(k,439)*y(k,216))
         mat(k,1583) = -rxt(k,439)*y(k,4)
         mat(k,289) = -(rxt(k,442)*y(k,216))
         mat(k,1615) = -rxt(k,442)*y(k,5)
         mat(k,380) = rxt(k,440)*y(k,198)
         mat(k,1977) = rxt(k,440)*y(k,186)
         mat(k,88) = .120_r8*rxt(k,439)*y(k,216)
         mat(k,1584) = .120_r8*rxt(k,439)*y(k,4)
         mat(k,886) = .100_r8*rxt(k,481)*y(k,118)
         mat(k,845) = .100_r8*rxt(k,484)*y(k,118)
         mat(k,1856) = .100_r8*rxt(k,481)*y(k,3) + .100_r8*rxt(k,484)*y(k,101)
         mat(k,1432) = .500_r8*rxt(k,441)*y(k,186) + .200_r8*rxt(k,468)*y(k,223) &
                      + .060_r8*rxt(k,474)*y(k,225)
         mat(k,381) = .500_r8*rxt(k,441)*y(k,111)
         mat(k,596) = .200_r8*rxt(k,468)*y(k,111)
         mat(k,627) = .060_r8*rxt(k,474)*y(k,111)
         mat(k,1426) = .200_r8*rxt(k,468)*y(k,223) + .200_r8*rxt(k,474)*y(k,225)
         mat(k,595) = .200_r8*rxt(k,468)*y(k,111)
         mat(k,625) = .200_r8*rxt(k,474)*y(k,111)
         mat(k,1443) = .200_r8*rxt(k,468)*y(k,223) + .150_r8*rxt(k,474)*y(k,225)
         mat(k,597) = .200_r8*rxt(k,468)*y(k,111)
         mat(k,628) = .150_r8*rxt(k,474)*y(k,111)
         mat(k,1427) = .210_r8*rxt(k,474)*y(k,225)
         mat(k,626) = .210_r8*rxt(k,474)*y(k,111)
         mat(k,160) = -(rxt(k,402)*y(k,216))
         mat(k,1595) = -rxt(k,402)*y(k,12)
         mat(k,885) = .050_r8*rxt(k,481)*y(k,118)
         mat(k,844) = .050_r8*rxt(k,484)*y(k,118)
         mat(k,1855) = .050_r8*rxt(k,481)*y(k,3) + .050_r8*rxt(k,484)*y(k,101)
         mat(k,269) = -(rxt(k,368)*y(k,113) + rxt(k,369)*y(k,216))
         mat(k,1914) = -rxt(k,368)*y(k,13)
         mat(k,1612) = -rxt(k,369)*y(k,13)
         mat(k,1334) = -(rxt(k,254)*y(k,38) + rxt(k,255)*y(k,198) + rxt(k,256) &
                      *y(k,118))
         mat(k,2108) = -rxt(k,254)*y(k,14)
         mat(k,2042) = -rxt(k,255)*y(k,14)
         mat(k,1891) = -rxt(k,256)*y(k,14)
         mat(k,1512) = 4.000_r8*rxt(k,257)*y(k,16) + (rxt(k,258)+rxt(k,259))*y(k,55) &
                      + rxt(k,262)*y(k,111) + rxt(k,265)*y(k,116) + rxt(k,507) &
                      *y(k,131) + rxt(k,266)*y(k,216)
         mat(k,1727) = (rxt(k,258)+rxt(k,259))*y(k,16)
         mat(k,760) = rxt(k,267)*y(k,116) + rxt(k,273)*y(k,212) + rxt(k,268)*y(k,216)
         mat(k,1487) = rxt(k,262)*y(k,16)
         mat(k,1552) = rxt(k,265)*y(k,16) + rxt(k,267)*y(k,72)
         mat(k,1301) = rxt(k,507)*y(k,16)
         mat(k,1753) = rxt(k,273)*y(k,72)
         mat(k,1700) = rxt(k,266)*y(k,16) + rxt(k,268)*y(k,72)
         mat(k,1505) = rxt(k,260)*y(k,55)
         mat(k,1720) = rxt(k,260)*y(k,16)
         mat(k,1315) = (rxt(k,550)+rxt(k,555))*y(k,82)
         mat(k,650) = (rxt(k,550)+rxt(k,555))*y(k,76)
         mat(k,1515) = -(4._r8*rxt(k,257)*y(k,16) + (rxt(k,258) + rxt(k,259) + rxt(k,260) &
                      ) * y(k,55) + rxt(k,261)*y(k,198) + rxt(k,262)*y(k,111) &
                      + rxt(k,263)*y(k,112) + rxt(k,265)*y(k,116) + rxt(k,266) &
                      *y(k,216) + rxt(k,507)*y(k,131))
         mat(k,1731) = -(rxt(k,258) + rxt(k,259) + rxt(k,260)) * y(k,16)
         mat(k,2046) = -rxt(k,261)*y(k,16)
         mat(k,1491) = -rxt(k,262)*y(k,16)
         mat(k,2088) = -rxt(k,263)*y(k,16)
         mat(k,1556) = -rxt(k,265)*y(k,16)
         mat(k,1704) = -rxt(k,266)*y(k,16)
         mat(k,1304) = -rxt(k,507)*y(k,16)
         mat(k,1336) = rxt(k,256)*y(k,118)
         mat(k,443) = rxt(k,264)*y(k,116)
         mat(k,761) = rxt(k,274)*y(k,212)
         mat(k,654) = rxt(k,269)*y(k,116)
         mat(k,1556) = mat(k,1556) + rxt(k,264)*y(k,17) + rxt(k,269)*y(k,82)
         mat(k,1895) = rxt(k,256)*y(k,14)
         mat(k,1757) = rxt(k,274)*y(k,72)
         mat(k,440) = -(rxt(k,264)*y(k,116))
         mat(k,1534) = -rxt(k,264)*y(k,17)
         mat(k,1507) = rxt(k,263)*y(k,112)
         mat(k,2066) = rxt(k,263)*y(k,16)
         mat(k,169) = -(rxt(k,443)*y(k,216))
         mat(k,1596) = -rxt(k,443)*y(k,18)
         mat(k,1425) = rxt(k,446)*y(k,187)
         mat(k,319) = rxt(k,446)*y(k,111)
         mat(k,248) = -(rxt(k,445)*y(k,216))
         mat(k,1608) = -rxt(k,445)*y(k,19)
         mat(k,320) = rxt(k,444)*y(k,198)
         mat(k,1973) = rxt(k,444)*y(k,187)
         mat(k,194) = -(rxt(k,317)*y(k,52) + rxt(k,318)*y(k,216))
         mat(k,1814) = -rxt(k,317)*y(k,20)
         mat(k,1600) = -rxt(k,318)*y(k,20)
         mat(k,464) = -(rxt(k,319)*y(k,52) + rxt(k,320)*y(k,118) + rxt(k,345)*y(k,216))
         mat(k,1816) = -rxt(k,319)*y(k,21)
         mat(k,1858) = -rxt(k,320)*y(k,21)
         mat(k,1636) = -rxt(k,345)*y(k,21)
         mat(k,177) = -(rxt(k,325)*y(k,216))
         mat(k,1598) = -rxt(k,325)*y(k,22)
         mat(k,803) = .800_r8*rxt(k,321)*y(k,188) + .200_r8*rxt(k,322)*y(k,192)
         mat(k,1345) = .200_r8*rxt(k,322)*y(k,188)
         mat(k,253) = -(rxt(k,326)*y(k,216))
         mat(k,1609) = -rxt(k,326)*y(k,23)
         mat(k,804) = rxt(k,323)*y(k,198)
         mat(k,1974) = rxt(k,323)*y(k,188)
         mat(k,200) = -(rxt(k,327)*y(k,52) + rxt(k,328)*y(k,216))
         mat(k,1815) = -rxt(k,327)*y(k,24)
         mat(k,1601) = -rxt(k,328)*y(k,24)
         mat(k,943) = -(rxt(k,348)*y(k,113) + rxt(k,349)*y(k,118) + rxt(k,366) &
                      *y(k,216))
         mat(k,1926) = -rxt(k,348)*y(k,25)
         mat(k,1872) = -rxt(k,349)*y(k,25)
         mat(k,1676) = -rxt(k,366)*y(k,25)
         mat(k,782) = .130_r8*rxt(k,426)*y(k,118)
         mat(k,1872) = mat(k,1872) + .130_r8*rxt(k,426)*y(k,89)
         mat(k,307) = -(rxt(k,353)*y(k,216))
         mat(k,1617) = -rxt(k,353)*y(k,26)
         mat(k,724) = rxt(k,351)*y(k,198)
         mat(k,1979) = rxt(k,351)*y(k,189)
         mat(k,56) = -(rxt(k,354)*y(k,216))
         mat(k,1579) = -rxt(k,354)*y(k,27)
         mat(k,181) = -(rxt(k,449)*y(k,216))
         mat(k,1599) = -rxt(k,449)*y(k,28)
         mat(k,511) = rxt(k,447)*y(k,198)
         mat(k,1968) = rxt(k,447)*y(k,190)
         mat(k,2124) = -(rxt(k,218)*y(k,52) + rxt(k,254)*y(k,14) + rxt(k,298)*y(k,198) &
                      + rxt(k,299)*y(k,113) + rxt(k,300)*y(k,116) + rxt(k,301) &
                      *y(k,216))
         mat(k,1846) = -rxt(k,218)*y(k,38)
         mat(k,1343) = -rxt(k,254)*y(k,38)
         mat(k,2058) = -rxt(k,298)*y(k,38)
         mat(k,1964) = -rxt(k,299)*y(k,38)
         mat(k,1568) = -rxt(k,300)*y(k,38)
         mat(k,1716) = -rxt(k,301)*y(k,38)
         mat(k,553) = .400_r8*rxt(k,397)*y(k,216)
         mat(k,907) = .340_r8*rxt(k,481)*y(k,118)
         mat(k,276) = .500_r8*rxt(k,368)*y(k,113)
         mat(k,471) = rxt(k,320)*y(k,118)
         mat(k,957) = .500_r8*rxt(k,349)*y(k,118)
         mat(k,430) = .500_r8*rxt(k,337)*y(k,216)
         mat(k,719) = rxt(k,306)*y(k,216)
         mat(k,317) = .300_r8*rxt(k,307)*y(k,216)
         mat(k,1743) = rxt(k,225)*y(k,192)
         mat(k,998) = .800_r8*rxt(k,342)*y(k,216)
         mat(k,795) = .910_r8*rxt(k,426)*y(k,118)
         mat(k,463) = .300_r8*rxt(k,417)*y(k,216)
         mat(k,1117) = .800_r8*rxt(k,421)*y(k,192)
         mat(k,1129) = .120_r8*rxt(k,379)*y(k,118)
         mat(k,455) = .500_r8*rxt(k,392)*y(k,216)
         mat(k,866) = .340_r8*rxt(k,484)*y(k,118)
         mat(k,1244) = .600_r8*rxt(k,393)*y(k,118)
         mat(k,1503) = .100_r8*rxt(k,399)*y(k,185) + rxt(k,305)*y(k,192) &
                      + .500_r8*rxt(k,370)*y(k,195) + .500_r8*rxt(k,339)*y(k,197) &
                      + .920_r8*rxt(k,409)*y(k,200) + .250_r8*rxt(k,377)*y(k,202) &
                      + rxt(k,386)*y(k,204) + rxt(k,360)*y(k,219) + rxt(k,364) &
                      *y(k,220) + .340_r8*rxt(k,493)*y(k,221) + .320_r8*rxt(k,498) &
                      *y(k,222) + .250_r8*rxt(k,434)*y(k,224)
         mat(k,1964) = mat(k,1964) + .500_r8*rxt(k,368)*y(k,13) + rxt(k,410)*y(k,200) &
                      + .250_r8*rxt(k,376)*y(k,202) + rxt(k,387)*y(k,204)
         mat(k,1907) = .340_r8*rxt(k,481)*y(k,3) + rxt(k,320)*y(k,21) &
                      + .500_r8*rxt(k,349)*y(k,25) + .910_r8*rxt(k,426)*y(k,89) &
                      + .120_r8*rxt(k,379)*y(k,96) + .340_r8*rxt(k,484)*y(k,101) &
                      + .600_r8*rxt(k,393)*y(k,102)
         mat(k,354) = rxt(k,344)*y(k,216)
         mat(k,982) = .680_r8*rxt(k,502)*y(k,216)
         mat(k,932) = .100_r8*rxt(k,399)*y(k,111)
         mat(k,814) = .700_r8*rxt(k,322)*y(k,192)
         mat(k,734) = rxt(k,350)*y(k,192)
         mat(k,1295) = rxt(k,333)*y(k,192) + rxt(k,406)*y(k,200) + .250_r8*rxt(k,373) &
                      *y(k,202) + rxt(k,382)*y(k,204) + .250_r8*rxt(k,431)*y(k,224)
         mat(k,1391) = rxt(k,225)*y(k,55) + .800_r8*rxt(k,421)*y(k,92) + rxt(k,305) &
                      *y(k,111) + .700_r8*rxt(k,322)*y(k,188) + rxt(k,350)*y(k,189) &
                      + rxt(k,333)*y(k,191) + (4.000_r8*rxt(k,302)+2.000_r8*rxt(k,303)) &
                      *y(k,192) + 1.500_r8*rxt(k,407)*y(k,200) + .750_r8*rxt(k,412) &
                      *y(k,201) + .880_r8*rxt(k,374)*y(k,202) + 2.000_r8*rxt(k,383) &
                      *y(k,204) + .750_r8*rxt(k,486)*y(k,211) + .800_r8*rxt(k,362) &
                      *y(k,220) + .930_r8*rxt(k,491)*y(k,221) + .950_r8*rxt(k,496) &
                      *y(k,222) + .800_r8*rxt(k,432)*y(k,224)
         mat(k,483) = .500_r8*rxt(k,370)*y(k,111)
         mat(k,677) = .500_r8*rxt(k,339)*y(k,111)
         mat(k,2058) = mat(k,2058) + .450_r8*rxt(k,384)*y(k,204) + .150_r8*rxt(k,363) &
                      *y(k,220)
         mat(k,1166) = .920_r8*rxt(k,409)*y(k,111) + rxt(k,410)*y(k,113) + rxt(k,406) &
                      *y(k,191) + 1.500_r8*rxt(k,407)*y(k,192)
         mat(k,1200) = .750_r8*rxt(k,412)*y(k,192)
         mat(k,1222) = .250_r8*rxt(k,377)*y(k,111) + .250_r8*rxt(k,376)*y(k,113) &
                      + .250_r8*rxt(k,373)*y(k,191) + .880_r8*rxt(k,374)*y(k,192)
         mat(k,1262) = rxt(k,386)*y(k,111) + rxt(k,387)*y(k,113) + rxt(k,382)*y(k,191) &
                      + 2.000_r8*rxt(k,383)*y(k,192) + .450_r8*rxt(k,384)*y(k,198) &
                      + 4.000_r8*rxt(k,385)*y(k,204)
         mat(k,973) = .750_r8*rxt(k,486)*y(k,192)
         mat(k,1716) = mat(k,1716) + .400_r8*rxt(k,397)*y(k,1) + .500_r8*rxt(k,337) &
                      *y(k,47) + rxt(k,306)*y(k,48) + .300_r8*rxt(k,307)*y(k,49) &
                      + .800_r8*rxt(k,342)*y(k,65) + .300_r8*rxt(k,417)*y(k,90) &
                      + .500_r8*rxt(k,392)*y(k,100) + rxt(k,344)*y(k,123) &
                      + .680_r8*rxt(k,502)*y(k,150)
         mat(k,714) = rxt(k,360)*y(k,111)
         mat(k,1052) = rxt(k,364)*y(k,111) + .800_r8*rxt(k,362)*y(k,192) &
                      + .150_r8*rxt(k,363)*y(k,198)
         mat(k,1014) = .340_r8*rxt(k,493)*y(k,111) + .930_r8*rxt(k,491)*y(k,192)
         mat(k,826) = .320_r8*rxt(k,498)*y(k,111) + .950_r8*rxt(k,496)*y(k,192)
         mat(k,1093) = .250_r8*rxt(k,434)*y(k,111) + .250_r8*rxt(k,431)*y(k,191) &
                      + .800_r8*rxt(k,432)*y(k,192)
      end do
      end subroutine nlnmat01
      subroutine nlnmat02( ofl, ofu, chnkpnts, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = ofl,ofu
         mat(k,1055) = -(rxt(k,329)*y(k,113) + rxt(k,330)*y(k,216))
         mat(k,1936) = -rxt(k,329)*y(k,41)
         mat(k,1686) = -rxt(k,330)*y(k,41)
         mat(k,548) = .800_r8*rxt(k,397)*y(k,216)
         mat(k,272) = rxt(k,368)*y(k,113)
         mat(k,178) = rxt(k,325)*y(k,216)
         mat(k,255) = .500_r8*rxt(k,326)*y(k,216)
         mat(k,946) = .500_r8*rxt(k,349)*y(k,118)
         mat(k,1228) = .100_r8*rxt(k,393)*y(k,118)
         mat(k,1476) = .400_r8*rxt(k,399)*y(k,185) + rxt(k,324)*y(k,188) &
                      + .270_r8*rxt(k,352)*y(k,189) + rxt(k,370)*y(k,195) + rxt(k,389) &
                      *y(k,206) + rxt(k,360)*y(k,219)
         mat(k,1936) = mat(k,1936) + rxt(k,368)*y(k,13)
         mat(k,1880) = .500_r8*rxt(k,349)*y(k,25) + .100_r8*rxt(k,393)*y(k,102)
         mat(k,924) = .400_r8*rxt(k,399)*y(k,111)
         mat(k,807) = rxt(k,324)*y(k,111) + 3.200_r8*rxt(k,321)*y(k,188) &
                      + .800_r8*rxt(k,322)*y(k,192)
         mat(k,727) = .270_r8*rxt(k,352)*y(k,111)
         mat(k,1368) = .800_r8*rxt(k,322)*y(k,188)
         mat(k,479) = rxt(k,370)*y(k,111)
         mat(k,2029) = .200_r8*rxt(k,388)*y(k,206)
         mat(k,556) = rxt(k,389)*y(k,111) + .200_r8*rxt(k,388)*y(k,198)
         mat(k,1686) = mat(k,1686) + .800_r8*rxt(k,397)*y(k,1) + rxt(k,325)*y(k,22) &
                      + .500_r8*rxt(k,326)*y(k,23)
         mat(k,708) = rxt(k,360)*y(k,111)
         mat(k,53) = -(rxt(k,331)*y(k,216))
         mat(k,1578) = -rxt(k,331)*y(k,43)
         mat(k,933) = -(rxt(k,367)*y(k,216))
         mat(k,1675) = -rxt(k,367)*y(k,44)
         mat(k,547) = .800_r8*rxt(k,397)*y(k,216)
         mat(k,893) = .520_r8*rxt(k,481)*y(k,118)
         mat(k,271) = .500_r8*rxt(k,368)*y(k,113)
         mat(k,852) = .520_r8*rxt(k,484)*y(k,118)
         mat(k,1466) = .250_r8*rxt(k,399)*y(k,185) + .820_r8*rxt(k,352)*y(k,189) &
                      + .500_r8*rxt(k,370)*y(k,195) + .270_r8*rxt(k,493)*y(k,221) &
                      + .040_r8*rxt(k,498)*y(k,222)
         mat(k,1925) = .500_r8*rxt(k,368)*y(k,13)
         mat(k,1871) = .520_r8*rxt(k,481)*y(k,3) + .520_r8*rxt(k,484)*y(k,101)
         mat(k,974) = .500_r8*rxt(k,502)*y(k,216)
         mat(k,923) = .250_r8*rxt(k,399)*y(k,111)
         mat(k,726) = .820_r8*rxt(k,352)*y(k,111) + .820_r8*rxt(k,350)*y(k,192)
         mat(k,1359) = .820_r8*rxt(k,350)*y(k,189) + .150_r8*rxt(k,491)*y(k,221) &
                      + .025_r8*rxt(k,496)*y(k,222)
         mat(k,477) = .500_r8*rxt(k,370)*y(k,111)
         mat(k,1675) = mat(k,1675) + .800_r8*rxt(k,397)*y(k,1) + .500_r8*rxt(k,502) &
                      *y(k,150)
         mat(k,1002) = .270_r8*rxt(k,493)*y(k,111) + .150_r8*rxt(k,491)*y(k,192)
         mat(k,820) = .040_r8*rxt(k,498)*y(k,111) + .025_r8*rxt(k,496)*y(k,192)
         mat(k,1132) = -(rxt(k,355)*y(k,113) + rxt(k,356)*y(k,216))
         mat(k,1940) = -rxt(k,355)*y(k,45)
         mat(k,1691) = -rxt(k,356)*y(k,45)
         mat(k,1036) = rxt(k,357)*y(k,216)
         mat(k,1121) = .880_r8*rxt(k,379)*y(k,118)
         mat(k,1229) = .500_r8*rxt(k,393)*y(k,118)
         mat(k,1480) = .170_r8*rxt(k,452)*y(k,193) + .050_r8*rxt(k,415)*y(k,201) &
                      + .250_r8*rxt(k,377)*y(k,202) + .170_r8*rxt(k,458)*y(k,205) &
                      + .400_r8*rxt(k,468)*y(k,223) + .250_r8*rxt(k,434)*y(k,224) &
                      + .540_r8*rxt(k,474)*y(k,225) + .510_r8*rxt(k,477)*y(k,226)
         mat(k,1940) = mat(k,1940) + .050_r8*rxt(k,416)*y(k,201) + .250_r8*rxt(k,376) &
                      *y(k,202) + .250_r8*rxt(k,435)*y(k,224)
         mat(k,798) = rxt(k,358)*y(k,216)
         mat(k,1883) = .880_r8*rxt(k,379)*y(k,96) + .500_r8*rxt(k,393)*y(k,102)
         mat(k,1278) = .250_r8*rxt(k,373)*y(k,202) + .250_r8*rxt(k,431)*y(k,224)
         mat(k,1372) = .240_r8*rxt(k,374)*y(k,202) + .500_r8*rxt(k,362)*y(k,220) &
                      + .100_r8*rxt(k,432)*y(k,224)
         mat(k,644) = .170_r8*rxt(k,452)*y(k,111) + .070_r8*rxt(k,451)*y(k,198)
         mat(k,2034) = .070_r8*rxt(k,451)*y(k,193) + .070_r8*rxt(k,457)*y(k,205)
         mat(k,1185) = .050_r8*rxt(k,415)*y(k,111) + .050_r8*rxt(k,416)*y(k,113)
         mat(k,1210) = .250_r8*rxt(k,377)*y(k,111) + .250_r8*rxt(k,376)*y(k,113) &
                      + .250_r8*rxt(k,373)*y(k,191) + .240_r8*rxt(k,374)*y(k,192)
         mat(k,987) = .170_r8*rxt(k,458)*y(k,111) + .070_r8*rxt(k,457)*y(k,198)
         mat(k,1691) = mat(k,1691) + rxt(k,357)*y(k,86) + rxt(k,358)*y(k,114)
         mat(k,1045) = .500_r8*rxt(k,362)*y(k,192)
         mat(k,605) = .400_r8*rxt(k,468)*y(k,111)
         mat(k,1084) = .250_r8*rxt(k,434)*y(k,111) + .250_r8*rxt(k,435)*y(k,113) &
                      + .250_r8*rxt(k,431)*y(k,191) + .100_r8*rxt(k,432)*y(k,192)
         mat(k,636) = .540_r8*rxt(k,474)*y(k,111)
         mat(k,399) = .510_r8*rxt(k,477)*y(k,111)
         mat(k,472) = -(rxt(k,336)*y(k,216))
         mat(k,1637) = -rxt(k,336)*y(k,46)
         mat(k,939) = .120_r8*rxt(k,349)*y(k,118)
         mat(k,1859) = .120_r8*rxt(k,349)*y(k,25)
         mat(k,1267) = .100_r8*rxt(k,333)*y(k,192) + .150_r8*rxt(k,334)*y(k,198)
         mat(k,1350) = .100_r8*rxt(k,333)*y(k,191)
         mat(k,1994) = .150_r8*rxt(k,334)*y(k,191) + .150_r8*rxt(k,384)*y(k,204)
         mat(k,1248) = .150_r8*rxt(k,384)*y(k,198)
         mat(k,426) = -(rxt(k,337)*y(k,216))
         mat(k,1632) = -rxt(k,337)*y(k,47)
         mat(k,1266) = .400_r8*rxt(k,334)*y(k,198)
         mat(k,1992) = .400_r8*rxt(k,334)*y(k,191) + .400_r8*rxt(k,384)*y(k,204)
         mat(k,1246) = .400_r8*rxt(k,384)*y(k,198)
         mat(k,716) = -(rxt(k,306)*y(k,216))
         mat(k,1657) = -rxt(k,306)*y(k,48)
         mat(k,1097) = .200_r8*rxt(k,421)*y(k,192)
         mat(k,805) = .300_r8*rxt(k,322)*y(k,192)
         mat(k,1351) = .200_r8*rxt(k,421)*y(k,92) + .300_r8*rxt(k,322)*y(k,188) &
                      + 2.000_r8*rxt(k,303)*y(k,192) + .250_r8*rxt(k,407)*y(k,200) &
                      + .250_r8*rxt(k,412)*y(k,201) + .250_r8*rxt(k,374)*y(k,202) &
                      + .250_r8*rxt(k,486)*y(k,211) + .500_r8*rxt(k,362)*y(k,220) &
                      + .250_r8*rxt(k,491)*y(k,221) + .250_r8*rxt(k,496)*y(k,222) &
                      + .300_r8*rxt(k,432)*y(k,224)
         mat(k,1142) = .250_r8*rxt(k,407)*y(k,192)
         mat(k,1173) = .250_r8*rxt(k,412)*y(k,192)
         mat(k,1203) = .250_r8*rxt(k,374)*y(k,192)
         mat(k,961) = .250_r8*rxt(k,486)*y(k,192)
         mat(k,1042) = .500_r8*rxt(k,362)*y(k,192)
         mat(k,1001) = .250_r8*rxt(k,491)*y(k,192)
         mat(k,817) = .250_r8*rxt(k,496)*y(k,192)
         mat(k,1078) = .300_r8*rxt(k,432)*y(k,192)
         mat(k,313) = -(rxt(k,307)*y(k,216))
         mat(k,1618) = -rxt(k,307)*y(k,49)
         mat(k,1348) = rxt(k,304)*y(k,198)
         mat(k,1980) = rxt(k,304)*y(k,192)
         mat(k,1841) = -(rxt(k,218)*y(k,38) + rxt(k,220)*y(k,68) + rxt(k,221)*y(k,70) &
                      + (rxt(k,222) + rxt(k,223)) * y(k,198) + rxt(k,224)*y(k,118) &
                      + rxt(k,231)*y(k,56) + rxt(k,240)*y(k,83) + rxt(k,327)*y(k,24))
         mat(k,2119) = -rxt(k,218)*y(k,52)
         mat(k,1074) = -rxt(k,220)*y(k,52)
         mat(k,501) = -rxt(k,221)*y(k,52)
         mat(k,2053) = -(rxt(k,222) + rxt(k,223)) * y(k,52)
         mat(k,1902) = -rxt(k,224)*y(k,52)
         mat(k,876) = -rxt(k,231)*y(k,52)
         mat(k,773) = -rxt(k,240)*y(k,52)
         mat(k,204) = -rxt(k,327)*y(k,52)
         mat(k,1522) = rxt(k,259)*y(k,55)
         mat(k,1738) = rxt(k,259)*y(k,16) + (4.000_r8*rxt(k,226)+2.000_r8*rxt(k,228)) &
                      *y(k,55) + rxt(k,230)*y(k,111) + rxt(k,235)*y(k,116) &
                      + rxt(k,508)*y(k,131) + rxt(k,225)*y(k,192) + rxt(k,236) &
                      *y(k,216)
         mat(k,111) = rxt(k,280)*y(k,212)
         mat(k,1328) = rxt(k,238)*y(k,116) + rxt(k,250)*y(k,212) + rxt(k,239)*y(k,216)
         mat(k,1498) = rxt(k,230)*y(k,55)
         mat(k,1563) = rxt(k,235)*y(k,55) + rxt(k,238)*y(k,76)
         mat(k,1309) = rxt(k,508)*y(k,55)
         mat(k,1386) = rxt(k,225)*y(k,55)
         mat(k,1764) = rxt(k,280)*y(k,60) + rxt(k,250)*y(k,76)
         mat(k,1711) = rxt(k,236)*y(k,55) + rxt(k,239)*y(k,76)
         mat(k,1813) = rxt(k,231)*y(k,56)
         mat(k,1719) = 2.000_r8*rxt(k,227)*y(k,55)
         mat(k,868) = rxt(k,231)*y(k,52) + (rxt(k,548)+rxt(k,553)+rxt(k,558))*y(k,76)
         mat(k,1314) = (rxt(k,548)+rxt(k,553)+rxt(k,558))*y(k,56) + (rxt(k,543) &
                       +rxt(k,549)+rxt(k,554))*y(k,83)
         mat(k,767) = (rxt(k,543)+rxt(k,549)+rxt(k,554))*y(k,76)
         mat(k,1718) = 2.000_r8*rxt(k,252)*y(k,55)
         mat(k,1734) = -(rxt(k,225)*y(k,192) + (4._r8*rxt(k,226) + 4._r8*rxt(k,227) &
                      + 4._r8*rxt(k,228) + 4._r8*rxt(k,252)) * y(k,55) + rxt(k,229) &
                      *y(k,198) + rxt(k,230)*y(k,111) + rxt(k,232)*y(k,112) + rxt(k,235) &
                      *y(k,116) + (rxt(k,236) + rxt(k,237)) * y(k,216) + (rxt(k,258) &
                      + rxt(k,259) + rxt(k,260)) * y(k,16) + rxt(k,508)*y(k,131))
         mat(k,1383) = -rxt(k,225)*y(k,55)
         mat(k,2049) = -rxt(k,229)*y(k,55)
         mat(k,1494) = -rxt(k,230)*y(k,55)
         mat(k,2091) = -rxt(k,232)*y(k,55)
         mat(k,1559) = -rxt(k,235)*y(k,55)
         mat(k,1707) = -(rxt(k,236) + rxt(k,237)) * y(k,55)
         mat(k,1518) = -(rxt(k,258) + rxt(k,259) + rxt(k,260)) * y(k,55)
         mat(k,1307) = -rxt(k,508)*y(k,55)
         mat(k,1837) = rxt(k,240)*y(k,83) + rxt(k,224)*y(k,118) + rxt(k,223)*y(k,198)
         mat(k,874) = rxt(k,233)*y(k,116)
         mat(k,1324) = rxt(k,251)*y(k,212)
         mat(k,772) = rxt(k,240)*y(k,52) + rxt(k,241)*y(k,116) + rxt(k,242)*y(k,216)
         mat(k,1559) = mat(k,1559) + rxt(k,233)*y(k,56) + rxt(k,241)*y(k,83)
         mat(k,1898) = rxt(k,224)*y(k,52)
         mat(k,236) = rxt(k,513)*y(k,131)
         mat(k,1307) = mat(k,1307) + rxt(k,513)*y(k,120)
         mat(k,2049) = mat(k,2049) + rxt(k,223)*y(k,52)
         mat(k,1760) = rxt(k,251)*y(k,76)
         mat(k,1707) = mat(k,1707) + rxt(k,242)*y(k,83)
         mat(k,870) = -(rxt(k,231)*y(k,52) + rxt(k,233)*y(k,116) + rxt(k,234)*y(k,216) &
                      + (rxt(k,548) + rxt(k,553) + rxt(k,558)) * y(k,76))
         mat(k,1823) = -rxt(k,231)*y(k,56)
         mat(k,1547) = -rxt(k,233)*y(k,56)
         mat(k,1671) = -rxt(k,234)*y(k,56)
         mat(k,1318) = -(rxt(k,548) + rxt(k,553) + rxt(k,558)) * y(k,56)
         mat(k,1724) = rxt(k,232)*y(k,112)
         mat(k,2075) = rxt(k,232)*y(k,55)
         mat(k,1031) = -((rxt(k,309) + rxt(k,316)) * y(k,216))
         mat(k,1683) = -(rxt(k,309) + rxt(k,316)) * y(k,57)
         mat(k,896) = .230_r8*rxt(k,481)*y(k,118)
         mat(k,1333) = rxt(k,254)*y(k,38)
         mat(k,197) = .350_r8*rxt(k,318)*y(k,216)
         mat(k,467) = .630_r8*rxt(k,320)*y(k,118)
         mat(k,944) = .560_r8*rxt(k,349)*y(k,118)
         mat(k,2105) = rxt(k,254)*y(k,14) + rxt(k,218)*y(k,52) + rxt(k,299)*y(k,113) &
                      + rxt(k,300)*y(k,116) + rxt(k,301)*y(k,216)
         mat(k,1131) = rxt(k,355)*y(k,113) + rxt(k,356)*y(k,216)
         mat(k,1825) = rxt(k,218)*y(k,38)
         mat(k,829) = rxt(k,343)*y(k,216)
         mat(k,783) = .620_r8*rxt(k,426)*y(k,118)
         mat(k,1119) = .650_r8*rxt(k,379)*y(k,118)
         mat(k,855) = .230_r8*rxt(k,484)*y(k,118)
         mat(k,1226) = .560_r8*rxt(k,393)*y(k,118)
         mat(k,1473) = .170_r8*rxt(k,452)*y(k,193) + .220_r8*rxt(k,377)*y(k,202) &
                      + .400_r8*rxt(k,455)*y(k,203) + .350_r8*rxt(k,458)*y(k,205) &
                      + .225_r8*rxt(k,493)*y(k,221) + .250_r8*rxt(k,434)*y(k,224)
         mat(k,1933) = rxt(k,299)*y(k,38) + rxt(k,355)*y(k,45) + .220_r8*rxt(k,376) &
                      *y(k,202) + .500_r8*rxt(k,435)*y(k,224)
         mat(k,1548) = rxt(k,300)*y(k,38) + rxt(k,503)*y(k,121)
         mat(k,1877) = .230_r8*rxt(k,481)*y(k,3) + .630_r8*rxt(k,320)*y(k,21) &
                      + .560_r8*rxt(k,349)*y(k,25) + .620_r8*rxt(k,426)*y(k,89) &
                      + .650_r8*rxt(k,379)*y(k,96) + .230_r8*rxt(k,484)*y(k,101) &
                      + .560_r8*rxt(k,393)*y(k,102)
         mat(k,264) = rxt(k,503)*y(k,116) + rxt(k,504)*y(k,216)
         mat(k,976) = .700_r8*rxt(k,502)*y(k,216)
         mat(k,1273) = .220_r8*rxt(k,373)*y(k,202) + .250_r8*rxt(k,431)*y(k,224)
         mat(k,1365) = .110_r8*rxt(k,374)*y(k,202) + .125_r8*rxt(k,491)*y(k,221) &
                      + .200_r8*rxt(k,432)*y(k,224)
         mat(k,643) = .170_r8*rxt(k,452)*y(k,111) + .070_r8*rxt(k,451)*y(k,198)
         mat(k,2026) = .070_r8*rxt(k,451)*y(k,193) + .160_r8*rxt(k,454)*y(k,203) &
                      + .140_r8*rxt(k,457)*y(k,205)
         mat(k,1206) = .220_r8*rxt(k,377)*y(k,111) + .220_r8*rxt(k,376)*y(k,113) &
                      + .220_r8*rxt(k,373)*y(k,191) + .110_r8*rxt(k,374)*y(k,192)
         mat(k,591) = .400_r8*rxt(k,455)*y(k,111) + .160_r8*rxt(k,454)*y(k,198)
         mat(k,986) = .350_r8*rxt(k,458)*y(k,111) + .140_r8*rxt(k,457)*y(k,198)
         mat(k,1683) = mat(k,1683) + .350_r8*rxt(k,318)*y(k,20) + rxt(k,301)*y(k,38) &
                      + rxt(k,356)*y(k,45) + rxt(k,343)*y(k,66) + rxt(k,504)*y(k,121) &
                      + .700_r8*rxt(k,502)*y(k,150)
         mat(k,1006) = .225_r8*rxt(k,493)*y(k,111) + .125_r8*rxt(k,491)*y(k,192)
         mat(k,1081) = .250_r8*rxt(k,434)*y(k,111) + .500_r8*rxt(k,435)*y(k,113) &
                      + .250_r8*rxt(k,431)*y(k,191) + .200_r8*rxt(k,432)*y(k,192)
      end do
      end subroutine nlnmat02
      subroutine nlnmat03( ofl, ofu, chnkpnts, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = ofl,ofu
         mat(k,60) = -(rxt(k,279)*y(k,212))
         mat(k,1745) = -rxt(k,279)*y(k,59)
         mat(k,108) = -(rxt(k,280)*y(k,212))
         mat(k,1748) = -rxt(k,280)*y(k,60)
         mat(k,124) = -(rxt(k,450)*y(k,216))
         mat(k,1589) = -rxt(k,450)*y(k,61)
         mat(k,118) = .180_r8*rxt(k,470)*y(k,216)
         mat(k,1589) = mat(k,1589) + .180_r8*rxt(k,470)*y(k,152)
         mat(k,206) = -(rxt(k,517)*y(k,113) + (rxt(k,518) + rxt(k,520)) * y(k,216))
         mat(k,1912) = -rxt(k,517)*y(k,62)
         mat(k,1602) = -(rxt(k,518) + rxt(k,520)) * y(k,62)
         mat(k,668) = rxt(k,338)*y(k,198)
         mat(k,1966) = rxt(k,338)*y(k,197)
         mat(k,660) = -(rxt(k,276)*y(k,68) + rxt(k,277)*y(k,227) + rxt(k,278)*y(k,80))
         mat(k,1065) = -rxt(k,276)*y(k,64)
         mat(k,2129) = -rxt(k,277)*y(k,64)
         mat(k,1772) = -rxt(k,278)*y(k,64)
         mat(k,61) = 2.000_r8*rxt(k,279)*y(k,212)
         mat(k,109) = rxt(k,280)*y(k,212)
         mat(k,1749) = 2.000_r8*rxt(k,279)*y(k,59) + rxt(k,280)*y(k,60)
         mat(k,994) = -(rxt(k,342)*y(k,216))
         mat(k,1680) = -rxt(k,342)*y(k,65)
         mat(k,457) = .700_r8*rxt(k,417)*y(k,216)
         mat(k,375) = .500_r8*rxt(k,418)*y(k,216)
         mat(k,224) = rxt(k,429)*y(k,216)
         mat(k,1470) = .050_r8*rxt(k,415)*y(k,201) + .530_r8*rxt(k,377)*y(k,202) &
                      + .225_r8*rxt(k,493)*y(k,221) + .250_r8*rxt(k,434)*y(k,224)
         mat(k,1930) = .050_r8*rxt(k,416)*y(k,201) + .530_r8*rxt(k,376)*y(k,202) &
                      + .250_r8*rxt(k,435)*y(k,224)
         mat(k,1406) = rxt(k,341)*y(k,196)
         mat(k,1272) = .530_r8*rxt(k,373)*y(k,202) + .250_r8*rxt(k,431)*y(k,224)
         mat(k,1363) = .260_r8*rxt(k,374)*y(k,202) + .125_r8*rxt(k,491)*y(k,221) &
                      + .100_r8*rxt(k,432)*y(k,224)
         mat(k,344) = rxt(k,341)*y(k,117)
         mat(k,1177) = .050_r8*rxt(k,415)*y(k,111) + .050_r8*rxt(k,416)*y(k,113)
         mat(k,1204) = .530_r8*rxt(k,377)*y(k,111) + .530_r8*rxt(k,376)*y(k,113) &
                      + .530_r8*rxt(k,373)*y(k,191) + .260_r8*rxt(k,374)*y(k,192)
         mat(k,1680) = mat(k,1680) + .700_r8*rxt(k,417)*y(k,90) + .500_r8*rxt(k,418) &
                      *y(k,91) + rxt(k,429)*y(k,106)
         mat(k,1004) = .225_r8*rxt(k,493)*y(k,111) + .125_r8*rxt(k,491)*y(k,192)
         mat(k,1080) = .250_r8*rxt(k,434)*y(k,111) + .250_r8*rxt(k,435)*y(k,113) &
                      + .250_r8*rxt(k,431)*y(k,191) + .100_r8*rxt(k,432)*y(k,192)
         mat(k,828) = -(rxt(k,343)*y(k,216))
         mat(k,1668) = -rxt(k,343)*y(k,66)
         mat(k,196) = .650_r8*rxt(k,318)*y(k,216)
         mat(k,993) = .200_r8*rxt(k,342)*y(k,216)
         mat(k,753) = rxt(k,430)*y(k,216)
         mat(k,1463) = rxt(k,441)*y(k,186) + .050_r8*rxt(k,415)*y(k,201) &
                      + .400_r8*rxt(k,455)*y(k,203) + .170_r8*rxt(k,458)*y(k,205) &
                      + .700_r8*rxt(k,461)*y(k,218) + .600_r8*rxt(k,468)*y(k,223) &
                      + .250_r8*rxt(k,434)*y(k,224) + .340_r8*rxt(k,474)*y(k,225) &
                      + .170_r8*rxt(k,477)*y(k,226)
         mat(k,1919) = .050_r8*rxt(k,416)*y(k,201) + .250_r8*rxt(k,435)*y(k,224)
         mat(k,384) = rxt(k,441)*y(k,111)
         mat(k,1270) = .250_r8*rxt(k,431)*y(k,224)
         mat(k,1357) = .100_r8*rxt(k,432)*y(k,224)
         mat(k,2017) = .160_r8*rxt(k,454)*y(k,203) + .070_r8*rxt(k,457)*y(k,205)
         mat(k,1175) = .050_r8*rxt(k,415)*y(k,111) + .050_r8*rxt(k,416)*y(k,113)
         mat(k,590) = .400_r8*rxt(k,455)*y(k,111) + .160_r8*rxt(k,454)*y(k,198)
         mat(k,984) = .170_r8*rxt(k,458)*y(k,111) + .070_r8*rxt(k,457)*y(k,198)
         mat(k,1668) = mat(k,1668) + .650_r8*rxt(k,318)*y(k,20) + .200_r8*rxt(k,342) &
                      *y(k,65) + rxt(k,430)*y(k,107)
         mat(k,335) = .700_r8*rxt(k,461)*y(k,111)
         mat(k,602) = .600_r8*rxt(k,468)*y(k,111)
         mat(k,1079) = .250_r8*rxt(k,434)*y(k,111) + .250_r8*rxt(k,435)*y(k,113) &
                      + .250_r8*rxt(k,431)*y(k,191) + .100_r8*rxt(k,432)*y(k,192)
         mat(k,633) = .340_r8*rxt(k,474)*y(k,111)
         mat(k,398) = .170_r8*rxt(k,477)*y(k,111)
         mat(k,1805) = -((rxt(k,176) + rxt(k,177) + rxt(k,178)) * y(k,198) + rxt(k,179) &
                      *y(k,117) + rxt(k,182)*y(k,118))
         mat(k,2052) = -(rxt(k,176) + rxt(k,177) + rxt(k,178)) * y(k,67)
         mat(k,1417) = -rxt(k,179)*y(k,67)
         mat(k,1901) = -rxt(k,182)*y(k,67)
         mat(k,2118) = rxt(k,301)*y(k,216)
         mat(k,1840) = rxt(k,220)*y(k,68)
         mat(k,1033) = rxt(k,316)*y(k,216)
         mat(k,665) = rxt(k,276)*y(k,68)
         mat(k,1073) = rxt(k,220)*y(k,52) + rxt(k,276)*y(k,64) + rxt(k,174)*y(k,116) &
                      + rxt(k,157)*y(k,212) + rxt(k,183)*y(k,216)
         mat(k,765) = rxt(k,274)*y(k,212)
         mat(k,1327) = rxt(k,251)*y(k,212)
         mat(k,748) = rxt(k,206)*y(k,216)
         mat(k,1562) = rxt(k,174)*y(k,68) + rxt(k,186)*y(k,216)
         mat(k,268) = rxt(k,504)*y(k,216)
         mat(k,615) = rxt(k,509)*y(k,216)
         mat(k,1308) = rxt(k,514)*y(k,216)
         mat(k,1763) = rxt(k,157)*y(k,68) + rxt(k,274)*y(k,72) + rxt(k,251)*y(k,76)
         mat(k,1710) = rxt(k,301)*y(k,38) + rxt(k,316)*y(k,57) + rxt(k,183)*y(k,68) &
                      + rxt(k,206)*y(k,103) + rxt(k,186)*y(k,116) + rxt(k,504) &
                      *y(k,121) + rxt(k,509)*y(k,130) + rxt(k,514)*y(k,131)
         mat(k,1066) = -(rxt(k,157)*y(k,212) + rxt(k,174)*y(k,116) + rxt(k,183) &
                      *y(k,216) + rxt(k,220)*y(k,52) + rxt(k,276)*y(k,64))
         mat(k,1751) = -rxt(k,157)*y(k,68)
         mat(k,1549) = -rxt(k,174)*y(k,68)
         mat(k,1687) = -rxt(k,183)*y(k,68)
         mat(k,1827) = -rxt(k,220)*y(k,68)
         mat(k,661) = -rxt(k,276)*y(k,68)
         mat(k,1793) = rxt(k,176)*y(k,198)
         mat(k,2030) = rxt(k,176)*y(k,67)
         mat(k,497) = -(rxt(k,175)*y(k,116) + rxt(k,184)*y(k,216) + rxt(k,221)*y(k,52))
         mat(k,1535) = -rxt(k,175)*y(k,70)
         mat(k,1641) = -rxt(k,184)*y(k,70)
         mat(k,1817) = -rxt(k,221)*y(k,70)
         mat(k,1996) = 2.000_r8*rxt(k,190)*y(k,198)
         mat(k,1641) = mat(k,1641) + 2.000_r8*rxt(k,189)*y(k,216)
         mat(k,172) = rxt(k,516)*y(k,227)
         mat(k,2126) = rxt(k,516)*y(k,133)
         mat(k,759) = -(rxt(k,267)*y(k,116) + rxt(k,268)*y(k,216) + (rxt(k,273) &
                      + rxt(k,274)) * y(k,212))
         mat(k,1544) = -rxt(k,267)*y(k,72)
         mat(k,1662) = -rxt(k,268)*y(k,72)
         mat(k,1750) = -(rxt(k,273) + rxt(k,274)) * y(k,72)
         mat(k,1332) = rxt(k,254)*y(k,38) + rxt(k,255)*y(k,198)
         mat(k,2104) = rxt(k,254)*y(k,14)
         mat(k,2012) = rxt(k,255)*y(k,14)
         mat(k,1319) = -(rxt(k,238)*y(k,116) + rxt(k,239)*y(k,216) + (rxt(k,250) &
                      + rxt(k,251)) * y(k,212) + (rxt(k,543) + rxt(k,549) + rxt(k,554) &
                      ) * y(k,83) + (rxt(k,548) + rxt(k,553) + rxt(k,558)) * y(k,56) &
                      + (rxt(k,550) + rxt(k,555)) * y(k,82))
         mat(k,1551) = -rxt(k,238)*y(k,76)
         mat(k,1699) = -rxt(k,239)*y(k,76)
         mat(k,1752) = -(rxt(k,250) + rxt(k,251)) * y(k,76)
         mat(k,769) = -(rxt(k,543) + rxt(k,549) + rxt(k,554)) * y(k,76)
         mat(k,871) = -(rxt(k,548) + rxt(k,553) + rxt(k,558)) * y(k,76)
         mat(k,652) = -(rxt(k,550) + rxt(k,555)) * y(k,76)
         mat(k,202) = rxt(k,327)*y(k,52)
         mat(k,2107) = rxt(k,218)*y(k,52)
         mat(k,1829) = rxt(k,327)*y(k,24) + rxt(k,218)*y(k,38) + rxt(k,220)*y(k,68) &
                      + rxt(k,221)*y(k,70) + rxt(k,240)*y(k,83) + rxt(k,222)*y(k,198)
         mat(k,1726) = rxt(k,237)*y(k,216)
         mat(k,1067) = rxt(k,220)*y(k,52)
         mat(k,498) = rxt(k,221)*y(k,52)
         mat(k,769) = mat(k,769) + rxt(k,240)*y(k,52)
         mat(k,2041) = rxt(k,222)*y(k,52)
         mat(k,1699) = mat(k,1699) + rxt(k,237)*y(k,55)
         mat(k,104) = -(rxt(k,310)*y(k,216) + rxt(k,315)*y(k,212))
         mat(k,1586) = -rxt(k,310)*y(k,77)
         mat(k,1747) = -rxt(k,315)*y(k,77)
         mat(k,720) = -(rxt(k,311)*y(k,216))
         mat(k,1658) = -rxt(k,311)*y(k,78)
         mat(k,887) = .050_r8*rxt(k,481)*y(k,118)
         mat(k,195) = .350_r8*rxt(k,318)*y(k,216)
         mat(k,466) = .370_r8*rxt(k,320)*y(k,118)
         mat(k,941) = .120_r8*rxt(k,349)*y(k,118)
         mat(k,780) = .110_r8*rxt(k,426)*y(k,118)
         mat(k,1118) = .330_r8*rxt(k,379)*y(k,118)
         mat(k,846) = .050_r8*rxt(k,484)*y(k,118)
         mat(k,1224) = .120_r8*rxt(k,393)*y(k,118)
         mat(k,1456) = rxt(k,314)*y(k,199)
         mat(k,1863) = .050_r8*rxt(k,481)*y(k,3) + .370_r8*rxt(k,320)*y(k,21) &
                      + .120_r8*rxt(k,349)*y(k,25) + .110_r8*rxt(k,426)*y(k,89) &
                      + .330_r8*rxt(k,379)*y(k,96) + .050_r8*rxt(k,484)*y(k,101) &
                      + .120_r8*rxt(k,393)*y(k,102)
         mat(k,2010) = rxt(k,312)*y(k,199)
         mat(k,328) = rxt(k,314)*y(k,111) + rxt(k,312)*y(k,198)
         mat(k,1658) = mat(k,1658) + .350_r8*rxt(k,318)*y(k,20)
         mat(k,659) = rxt(k,276)*y(k,68) + rxt(k,278)*y(k,80) + rxt(k,277)*y(k,227)
         mat(k,1064) = rxt(k,276)*y(k,64)
         mat(k,1771) = rxt(k,278)*y(k,64)
         mat(k,2127) = rxt(k,277)*y(k,64)
         mat(k,1784) = -(rxt(k,215)*y(k,216) + rxt(k,278)*y(k,64))
         mat(k,1709) = -rxt(k,215)*y(k,80)
         mat(k,664) = -rxt(k,278)*y(k,80)
         mat(k,2117) = rxt(k,299)*y(k,113)
         mat(k,1059) = rxt(k,329)*y(k,113)
         mat(k,1135) = rxt(k,355)*y(k,113)
         mat(k,875) = (rxt(k,548)+rxt(k,553)+rxt(k,558))*y(k,76)
         mat(k,209) = rxt(k,517)*y(k,113)
         mat(k,1326) = (rxt(k,548)+rxt(k,553)+rxt(k,558))*y(k,56)
         mat(k,2093) = rxt(k,214)*y(k,216)
         mat(k,1957) = rxt(k,299)*y(k,38) + rxt(k,329)*y(k,41) + rxt(k,355)*y(k,45) &
                      + rxt(k,517)*y(k,62)
         mat(k,1709) = mat(k,1709) + rxt(k,214)*y(k,112)
         mat(k,359) = -(rxt(k,191)*y(k,216))
         mat(k,1624) = -rxt(k,191)*y(k,81)
         mat(k,2064) = rxt(k,212)*y(k,198)
         mat(k,1987) = rxt(k,212)*y(k,112)
         mat(k,651) = -(rxt(k,269)*y(k,116) + (rxt(k,550) + rxt(k,555)) * y(k,76))
         mat(k,1539) = -rxt(k,269)*y(k,82)
         mat(k,1316) = -(rxt(k,550) + rxt(k,555)) * y(k,82)
         mat(k,1508) = rxt(k,261)*y(k,198)
         mat(k,2007) = rxt(k,261)*y(k,16)
         mat(k,768) = -(rxt(k,240)*y(k,52) + rxt(k,241)*y(k,116) + rxt(k,242)*y(k,216) &
                      + (rxt(k,543) + rxt(k,549) + rxt(k,554)) * y(k,76))
         mat(k,1820) = -rxt(k,240)*y(k,83)
         mat(k,1545) = -rxt(k,241)*y(k,83)
         mat(k,1663) = -rxt(k,242)*y(k,83)
         mat(k,1317) = -(rxt(k,543) + rxt(k,549) + rxt(k,554)) * y(k,83)
         mat(k,1722) = rxt(k,229)*y(k,198)
         mat(k,869) = rxt(k,234)*y(k,216)
         mat(k,2013) = rxt(k,229)*y(k,55)
         mat(k,1663) = mat(k,1663) + rxt(k,234)*y(k,56)
         mat(k,1019) = -(rxt(k,372)*y(k,216))
         mat(k,1682) = -rxt(k,372)*y(k,84)
         mat(k,458) = .300_r8*rxt(k,417)*y(k,216)
         mat(k,376) = .500_r8*rxt(k,418)*y(k,216)
         mat(k,1472) = rxt(k,371)*y(k,195) + rxt(k,378)*y(k,202)
         mat(k,478) = rxt(k,371)*y(k,111)
         mat(k,1205) = rxt(k,378)*y(k,111)
         mat(k,1682) = mat(k,1682) + .300_r8*rxt(k,417)*y(k,90) + .500_r8*rxt(k,418) &
                      *y(k,91)
      end do
      end subroutine nlnmat03
      subroutine nlnmat04( ofl, ofu, chnkpnts, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = ofl,ofu
         mat(k,155) = -(rxt(k,403)*y(k,216))
         mat(k,1594) = -rxt(k,403)*y(k,85)
         mat(k,1035) = -(rxt(k,357)*y(k,216))
         mat(k,1684) = -rxt(k,357)*y(k,86)
         mat(k,459) = .700_r8*rxt(k,417)*y(k,216)
         mat(k,377) = .500_r8*rxt(k,418)*y(k,216)
         mat(k,449) = .500_r8*rxt(k,392)*y(k,216)
         mat(k,1474) = .050_r8*rxt(k,415)*y(k,201) + .220_r8*rxt(k,377)*y(k,202) &
                      + .250_r8*rxt(k,434)*y(k,224)
         mat(k,1934) = .050_r8*rxt(k,416)*y(k,201) + .220_r8*rxt(k,376)*y(k,202) &
                      + .250_r8*rxt(k,435)*y(k,224)
         mat(k,434) = .500_r8*rxt(k,361)*y(k,216)
         mat(k,1274) = .220_r8*rxt(k,373)*y(k,202) + .250_r8*rxt(k,431)*y(k,224)
         mat(k,1366) = .230_r8*rxt(k,374)*y(k,202) + .200_r8*rxt(k,362)*y(k,220) &
                      + .100_r8*rxt(k,432)*y(k,224)
         mat(k,1180) = .050_r8*rxt(k,415)*y(k,111) + .050_r8*rxt(k,416)*y(k,113)
         mat(k,1207) = .220_r8*rxt(k,377)*y(k,111) + .220_r8*rxt(k,376)*y(k,113) &
                      + .220_r8*rxt(k,373)*y(k,191) + .230_r8*rxt(k,374)*y(k,192)
         mat(k,1684) = mat(k,1684) + .700_r8*rxt(k,417)*y(k,90) + .500_r8*rxt(k,418) &
                      *y(k,91) + .500_r8*rxt(k,392)*y(k,100) + .500_r8*rxt(k,361) &
                      *y(k,128)
         mat(k,1043) = .200_r8*rxt(k,362)*y(k,192)
         mat(k,1082) = .250_r8*rxt(k,434)*y(k,111) + .250_r8*rxt(k,435)*y(k,113) &
                      + .250_r8*rxt(k,431)*y(k,191) + .100_r8*rxt(k,432)*y(k,192)
         mat(k,258) = -(rxt(k,404)*y(k,216))
         mat(k,1610) = -rxt(k,404)*y(k,87)
         mat(k,1428) = .870_r8*rxt(k,415)*y(k,201)
         mat(k,1913) = .950_r8*rxt(k,416)*y(k,201)
         mat(k,1264) = rxt(k,411)*y(k,201)
         mat(k,1347) = .750_r8*rxt(k,412)*y(k,201)
         mat(k,1169) = .870_r8*rxt(k,415)*y(k,111) + .950_r8*rxt(k,416)*y(k,113) &
                      + rxt(k,411)*y(k,191) + .750_r8*rxt(k,412)*y(k,192)
         mat(k,69) = -(rxt(k,405)*y(k,216))
         mat(k,1581) = -rxt(k,405)*y(k,88)
         mat(k,617) = .600_r8*rxt(k,428)*y(k,216)
         mat(k,1581) = mat(k,1581) + .600_r8*rxt(k,428)*y(k,94)
         mat(k,781) = -(rxt(k,419)*y(k,113) + rxt(k,426)*y(k,118) + rxt(k,427) &
                      *y(k,216))
         mat(k,1917) = -rxt(k,419)*y(k,89)
         mat(k,1864) = -rxt(k,426)*y(k,89)
         mat(k,1664) = -rxt(k,427)*y(k,89)
         mat(k,456) = -(rxt(k,417)*y(k,216))
         mat(k,1635) = -rxt(k,417)*y(k,90)
         mat(k,1441) = .080_r8*rxt(k,409)*y(k,200)
         mat(k,1140) = .080_r8*rxt(k,409)*y(k,111)
         mat(k,373) = -(rxt(k,418)*y(k,216))
         mat(k,1626) = -rxt(k,418)*y(k,91)
         mat(k,1435) = .080_r8*rxt(k,415)*y(k,201)
         mat(k,1170) = .080_r8*rxt(k,415)*y(k,111)
         mat(k,1105) = -(rxt(k,420)*y(k,191) + rxt(k,421)*y(k,192) + rxt(k,422) &
                      *y(k,198) + rxt(k,423)*y(k,111) + rxt(k,424)*y(k,113))
         mat(k,1276) = -rxt(k,420)*y(k,92)
         mat(k,1370) = -rxt(k,421)*y(k,92)
         mat(k,2032) = -rxt(k,422)*y(k,92)
         mat(k,1478) = -rxt(k,423)*y(k,92)
         mat(k,1938) = -rxt(k,424)*y(k,92)
         mat(k,784) = rxt(k,419)*y(k,113)
         mat(k,1938) = mat(k,1938) + rxt(k,419)*y(k,89)
         mat(k,277) = -(rxt(k,425)*y(k,216))
         mat(k,1613) = -rxt(k,425)*y(k,93)
         mat(k,1095) = rxt(k,422)*y(k,198)
         mat(k,1975) = rxt(k,422)*y(k,92)
         mat(k,618) = -(rxt(k,428)*y(k,216))
         mat(k,1651) = -rxt(k,428)*y(k,94)
         mat(k,2004) = rxt(k,408)*y(k,200) + rxt(k,413)*y(k,201)
         mat(k,1141) = rxt(k,408)*y(k,198)
         mat(k,1172) = rxt(k,413)*y(k,198)
         mat(k,40) = -(rxt(k,535)*y(k,216))
         mat(k,1575) = -rxt(k,535)*y(k,95)
         mat(k,1120) = -(rxt(k,379)*y(k,118) + rxt(k,380)*y(k,216))
         mat(k,1882) = -rxt(k,379)*y(k,96)
         mat(k,1690) = -rxt(k,380)*y(k,96)
         mat(k,785) = .300_r8*rxt(k,426)*y(k,118)
         mat(k,1479) = .360_r8*rxt(k,409)*y(k,200)
         mat(k,1939) = .400_r8*rxt(k,410)*y(k,200)
         mat(k,1882) = mat(k,1882) + .300_r8*rxt(k,426)*y(k,89)
         mat(k,1277) = .390_r8*rxt(k,406)*y(k,200)
         mat(k,1371) = .310_r8*rxt(k,407)*y(k,200)
         mat(k,1150) = .360_r8*rxt(k,409)*y(k,111) + .400_r8*rxt(k,410)*y(k,113) &
                      + .390_r8*rxt(k,406)*y(k,191) + .310_r8*rxt(k,407)*y(k,192)
         mat(k,212) = -(rxt(k,381)*y(k,216))
         mat(k,1603) = -rxt(k,381)*y(k,97)
         mat(k,1969) = rxt(k,375)*y(k,202)
         mat(k,1202) = rxt(k,375)*y(k,198)
         mat(k,404) = -(rxt(k,390)*y(k,216))
         mat(k,1629) = -rxt(k,390)*y(k,98)
         mat(k,1438) = .800_r8*rxt(k,399)*y(k,185)
         mat(k,917) = .800_r8*rxt(k,399)*y(k,111)
         mat(k,217) = -(rxt(k,391)*y(k,216))
         mat(k,1604) = -rxt(k,391)*y(k,99)
         mat(k,1970) = .800_r8*rxt(k,388)*y(k,206)
         mat(k,554) = .800_r8*rxt(k,388)*y(k,198)
         mat(k,448) = -(rxt(k,392)*y(k,216))
         mat(k,1634) = -rxt(k,392)*y(k,100)
         mat(k,2067) = rxt(k,395)*y(k,204)
         mat(k,1247) = rxt(k,395)*y(k,112)
         mat(k,849) = -(rxt(k,483)*y(k,113) + rxt(k,484)*y(k,118) + rxt(k,485) &
                      *y(k,216))
         mat(k,1921) = -rxt(k,483)*y(k,101)
         mat(k,1867) = -rxt(k,484)*y(k,101)
         mat(k,1670) = -rxt(k,485)*y(k,101)
         mat(k,1231) = -(rxt(k,393)*y(k,118) + rxt(k,394)*y(k,216))
         mat(k,1887) = -rxt(k,393)*y(k,102)
         mat(k,1695) = -rxt(k,394)*y(k,102)
         mat(k,788) = .200_r8*rxt(k,426)*y(k,118)
         mat(k,1484) = .560_r8*rxt(k,409)*y(k,200)
         mat(k,1944) = .600_r8*rxt(k,410)*y(k,200)
         mat(k,1887) = mat(k,1887) + .200_r8*rxt(k,426)*y(k,89)
         mat(k,1282) = .610_r8*rxt(k,406)*y(k,200)
         mat(k,1376) = .440_r8*rxt(k,407)*y(k,200)
         mat(k,1154) = .560_r8*rxt(k,409)*y(k,111) + .600_r8*rxt(k,410)*y(k,113) &
                      + .610_r8*rxt(k,406)*y(k,191) + .440_r8*rxt(k,407)*y(k,192)
         mat(k,741) = -(rxt(k,194)*y(k,111) + (rxt(k,195) + rxt(k,196) + rxt(k,197) &
                      ) * y(k,112) + rxt(k,198)*y(k,117) + rxt(k,206)*y(k,216) &
                      + rxt(k,568)*y(k,215))
         mat(k,1458) = -rxt(k,194)*y(k,103)
         mat(k,2072) = -(rxt(k,195) + rxt(k,196) + rxt(k,197)) * y(k,103)
         mat(k,1404) = -rxt(k,198)*y(k,103)
         mat(k,1660) = -rxt(k,206)*y(k,103)
         mat(k,681) = -rxt(k,568)*y(k,103)
         mat(k,1543) = rxt(k,192)*y(k,207) + rxt(k,565)*y(k,210)
         mat(k,1404) = mat(k,1404) + rxt(k,566)*y(k,210)
         mat(k,699) = 1.100_r8*rxt(k,561)*y(k,208) + .200_r8*rxt(k,559)*y(k,209)
         mat(k,422) = rxt(k,192)*y(k,116)
         mat(k,568) = 1.100_r8*rxt(k,561)*y(k,194)
         mat(k,689) = .200_r8*rxt(k,559)*y(k,194)
         mat(k,393) = rxt(k,565)*y(k,116) + rxt(k,566)*y(k,117)
         mat(k,2061) = rxt(k,213)*y(k,113)
         mat(k,1911) = rxt(k,213)*y(k,112)
         mat(k,222) = -(rxt(k,429)*y(k,216))
         mat(k,1605) = -rxt(k,429)*y(k,106)
         mat(k,1094) = .200_r8*rxt(k,421)*y(k,192)
         mat(k,1346) = .200_r8*rxt(k,421)*y(k,92)
         mat(k,751) = -(rxt(k,430)*y(k,216))
         mat(k,1661) = -rxt(k,430)*y(k,107)
         mat(k,1098) = rxt(k,423)*y(k,111) + rxt(k,424)*y(k,113) + rxt(k,420)*y(k,191) &
                      + .800_r8*rxt(k,421)*y(k,192)
         mat(k,1459) = rxt(k,423)*y(k,92)
         mat(k,1916) = rxt(k,424)*y(k,92)
         mat(k,1268) = rxt(k,420)*y(k,92)
         mat(k,1353) = .800_r8*rxt(k,421)*y(k,92)
         mat(k,50) = -(rxt(k,519)*y(k,216))
         mat(k,1577) = -rxt(k,519)*y(k,108)
         mat(k,1490) = -(rxt(k,194)*y(k,103) + rxt(k,203)*y(k,113) + rxt(k,207) &
                      *y(k,198) + rxt(k,208)*y(k,118) + rxt(k,209)*y(k,116) + rxt(k,230) &
                      *y(k,55) + rxt(k,262)*y(k,16) + rxt(k,305)*y(k,192) + rxt(k,314) &
                      *y(k,199) + rxt(k,324)*y(k,188) + rxt(k,335)*y(k,191) + rxt(k,339) &
                      *y(k,197) + rxt(k,352)*y(k,189) + rxt(k,360)*y(k,219) + rxt(k,364) &
                      *y(k,220) + (rxt(k,370) + rxt(k,371)) * y(k,195) + (rxt(k,377) &
                      + rxt(k,378)) * y(k,202) + rxt(k,386)*y(k,204) + rxt(k,389) &
                      *y(k,206) + (rxt(k,399) + rxt(k,400)) * y(k,185) + rxt(k,409) &
                      *y(k,200) + rxt(k,415)*y(k,201) + rxt(k,423)*y(k,92) + rxt(k,434) &
                      *y(k,224) + rxt(k,438)*y(k,184) + rxt(k,441)*y(k,186) + rxt(k,446) &
                      *y(k,187) + rxt(k,448)*y(k,190) + rxt(k,452)*y(k,193) + rxt(k,455) &
                      *y(k,203) + rxt(k,458)*y(k,205) + rxt(k,461)*y(k,218) + rxt(k,468) &
                      *y(k,223) + rxt(k,474)*y(k,225) + rxt(k,477)*y(k,226) + rxt(k,488) &
                      *y(k,211) + rxt(k,493)*y(k,221) + rxt(k,498)*y(k,222) + rxt(k,570) &
                      *y(k,215))
         mat(k,744) = -rxt(k,194)*y(k,111)
         mat(k,1951) = -rxt(k,203)*y(k,111)
         mat(k,2045) = -rxt(k,207)*y(k,111)
         mat(k,1894) = -rxt(k,208)*y(k,111)
         mat(k,1555) = -rxt(k,209)*y(k,111)
         mat(k,1730) = -rxt(k,230)*y(k,111)
         mat(k,1514) = -rxt(k,262)*y(k,111)
         mat(k,1381) = -rxt(k,305)*y(k,111)
         mat(k,329) = -rxt(k,314)*y(k,111)
         mat(k,810) = -rxt(k,324)*y(k,111)
         mat(k,1287) = -rxt(k,335)*y(k,111)
         mat(k,673) = -rxt(k,339)*y(k,111)
         mat(k,730) = -rxt(k,352)*y(k,111)
         mat(k,710) = -rxt(k,360)*y(k,111)
         mat(k,1048) = -rxt(k,364)*y(k,111)
         mat(k,480) = -(rxt(k,370) + rxt(k,371)) * y(k,111)
         mat(k,1215) = -(rxt(k,377) + rxt(k,378)) * y(k,111)
         mat(k,1254) = -rxt(k,386)*y(k,111)
         mat(k,558) = -rxt(k,389)*y(k,111)
         mat(k,928) = -(rxt(k,399) + rxt(k,400)) * y(k,111)
         mat(k,1158) = -rxt(k,409)*y(k,111)
         mat(k,1192) = -rxt(k,415)*y(k,111)
         mat(k,1111) = -rxt(k,423)*y(k,111)
         mat(k,1087) = -rxt(k,434)*y(k,111)
         mat(k,417) = -rxt(k,438)*y(k,111)
         mat(k,385) = -rxt(k,441)*y(k,111)
         mat(k,323) = -rxt(k,446)*y(k,111)
         mat(k,514) = -rxt(k,448)*y(k,111)
         mat(k,646) = -rxt(k,452)*y(k,111)
         mat(k,592) = -rxt(k,455)*y(k,111)
         mat(k,989) = -rxt(k,458)*y(k,111)
         mat(k,336) = -rxt(k,461)*y(k,111)
         mat(k,606) = -rxt(k,468)*y(k,111)
         mat(k,638) = -rxt(k,474)*y(k,111)
         mat(k,400) = -rxt(k,477)*y(k,111)
         mat(k,968) = -rxt(k,488)*y(k,111)
         mat(k,1010) = -rxt(k,493)*y(k,111)
         mat(k,822) = -rxt(k,498)*y(k,111)
         mat(k,683) = -rxt(k,570)*y(k,111)
         mat(k,744) = mat(k,744) + 2.000_r8*rxt(k,196)*y(k,112) + rxt(k,198)*y(k,117) &
                      + rxt(k,206)*y(k,216)
         mat(k,2087) = 2.000_r8*rxt(k,196)*y(k,103) + rxt(k,199)*y(k,116) + rxt(k,510) &
                      *y(k,131)
         mat(k,1555) = mat(k,1555) + rxt(k,199)*y(k,112)
         mat(k,1411) = rxt(k,198)*y(k,103) + rxt(k,193)*y(k,207)
         mat(k,1303) = rxt(k,510)*y(k,112)
         mat(k,424) = rxt(k,193)*y(k,117)
         mat(k,1703) = rxt(k,206)*y(k,103)
      end do
      end subroutine nlnmat04
      subroutine nlnmat05( ofl, ofu, chnkpnts, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = ofl,ofu
         mat(k,2099) = -((rxt(k,195) + rxt(k,196) + rxt(k,197)) * y(k,103) + (rxt(k,199) &
                      + rxt(k,201)) * y(k,116) + rxt(k,200)*y(k,118) + rxt(k,212) &
                      *y(k,198) + rxt(k,213)*y(k,113) + rxt(k,214)*y(k,216) + rxt(k,232) &
                      *y(k,55) + rxt(k,263)*y(k,16) + rxt(k,346)*y(k,191) + rxt(k,395) &
                      *y(k,204) + rxt(k,453)*y(k,193) + rxt(k,456)*y(k,203) + rxt(k,459) &
                      *y(k,205) + rxt(k,463)*y(k,125) + rxt(k,466)*y(k,184) + rxt(k,510) &
                      *y(k,131))
         mat(k,749) = -(rxt(k,195) + rxt(k,196) + rxt(k,197)) * y(k,112)
         mat(k,1567) = -(rxt(k,199) + rxt(k,201)) * y(k,112)
         mat(k,1906) = -rxt(k,200)*y(k,112)
         mat(k,2057) = -rxt(k,212)*y(k,112)
         mat(k,1963) = -rxt(k,213)*y(k,112)
         mat(k,1715) = -rxt(k,214)*y(k,112)
         mat(k,1742) = -rxt(k,232)*y(k,112)
         mat(k,1526) = -rxt(k,263)*y(k,112)
         mat(k,1294) = -rxt(k,346)*y(k,112)
         mat(k,1261) = -rxt(k,395)*y(k,112)
         mat(k,649) = -rxt(k,453)*y(k,112)
         mat(k,594) = -rxt(k,456)*y(k,112)
         mat(k,992) = -rxt(k,459)*y(k,112)
         mat(k,358) = -rxt(k,463)*y(k,112)
         mat(k,420) = -rxt(k,466)*y(k,112)
         mat(k,1312) = -rxt(k,510)*y(k,112)
         mat(k,552) = rxt(k,397)*y(k,216)
         mat(k,275) = rxt(k,368)*y(k,113)
         mat(k,1526) = mat(k,1526) + rxt(k,262)*y(k,111)
         mat(k,1742) = mat(k,1742) + rxt(k,230)*y(k,111)
         mat(k,364) = rxt(k,191)*y(k,216)
         mat(k,462) = .700_r8*rxt(k,417)*y(k,216)
         mat(k,1116) = rxt(k,423)*y(k,111) + rxt(k,424)*y(k,113)
         mat(k,1502) = rxt(k,262)*y(k,16) + rxt(k,230)*y(k,55) + rxt(k,423)*y(k,92) &
                      + 2.000_r8*rxt(k,203)*y(k,113) + rxt(k,209)*y(k,116) &
                      + rxt(k,208)*y(k,118) + rxt(k,438)*y(k,184) + rxt(k,399) &
                      *y(k,185) + rxt(k,441)*y(k,186) + rxt(k,446)*y(k,187) &
                      + rxt(k,324)*y(k,188) + rxt(k,352)*y(k,189) + rxt(k,448) &
                      *y(k,190) + rxt(k,335)*y(k,191) + rxt(k,305)*y(k,192) &
                      + rxt(k,452)*y(k,193) + rxt(k,370)*y(k,195) + rxt(k,339) &
                      *y(k,197) + rxt(k,207)*y(k,198) + rxt(k,314)*y(k,199) &
                      + .920_r8*rxt(k,409)*y(k,200) + .920_r8*rxt(k,415)*y(k,201) &
                      + rxt(k,377)*y(k,202) + rxt(k,455)*y(k,203) + rxt(k,386) &
                      *y(k,204) + rxt(k,458)*y(k,205) + rxt(k,389)*y(k,206) &
                      + 1.600_r8*rxt(k,488)*y(k,211) + rxt(k,461)*y(k,218) &
                      + rxt(k,360)*y(k,219) + rxt(k,364)*y(k,220) + .900_r8*rxt(k,493) &
                      *y(k,221) + .800_r8*rxt(k,498)*y(k,222) + rxt(k,468)*y(k,223) &
                      + rxt(k,434)*y(k,224) + rxt(k,474)*y(k,225) + rxt(k,477) &
                      *y(k,226)
         mat(k,1963) = mat(k,1963) + rxt(k,368)*y(k,13) + rxt(k,424)*y(k,92) &
                      + 2.000_r8*rxt(k,203)*y(k,111) + rxt(k,204)*y(k,116) &
                      + rxt(k,202)*y(k,198) + rxt(k,410)*y(k,200) + rxt(k,416) &
                      *y(k,201) + rxt(k,376)*y(k,202) + rxt(k,387)*y(k,204) &
                      + 2.000_r8*rxt(k,489)*y(k,211) + rxt(k,205)*y(k,216) &
                      + rxt(k,435)*y(k,224)
         mat(k,801) = rxt(k,358)*y(k,216)
         mat(k,1567) = mat(k,1567) + rxt(k,209)*y(k,111) + rxt(k,204)*y(k,113)
         mat(k,1906) = mat(k,1906) + rxt(k,208)*y(k,111)
         mat(k,413) = rxt(k,495)*y(k,216)
         mat(k,420) = mat(k,420) + rxt(k,438)*y(k,111)
         mat(k,931) = rxt(k,399)*y(k,111)
         mat(k,388) = rxt(k,441)*y(k,111)
         mat(k,326) = rxt(k,446)*y(k,111)
         mat(k,813) = rxt(k,324)*y(k,111)
         mat(k,733) = rxt(k,352)*y(k,111)
         mat(k,518) = rxt(k,448)*y(k,111)
         mat(k,1294) = mat(k,1294) + rxt(k,335)*y(k,111)
         mat(k,1390) = rxt(k,305)*y(k,111) + .500_r8*rxt(k,486)*y(k,211)
         mat(k,649) = mat(k,649) + rxt(k,452)*y(k,111)
         mat(k,482) = rxt(k,370)*y(k,111)
         mat(k,676) = rxt(k,339)*y(k,111)
         mat(k,2057) = mat(k,2057) + rxt(k,207)*y(k,111) + rxt(k,202)*y(k,113)
         mat(k,331) = rxt(k,314)*y(k,111)
         mat(k,1165) = .920_r8*rxt(k,409)*y(k,111) + rxt(k,410)*y(k,113)
         mat(k,1199) = .920_r8*rxt(k,415)*y(k,111) + rxt(k,416)*y(k,113)
         mat(k,1221) = rxt(k,377)*y(k,111) + rxt(k,376)*y(k,113)
         mat(k,594) = mat(k,594) + rxt(k,455)*y(k,111)
         mat(k,1261) = mat(k,1261) + rxt(k,386)*y(k,111) + rxt(k,387)*y(k,113)
         mat(k,992) = mat(k,992) + rxt(k,458)*y(k,111)
         mat(k,561) = rxt(k,389)*y(k,111)
         mat(k,972) = 1.600_r8*rxt(k,488)*y(k,111) + 2.000_r8*rxt(k,489)*y(k,113) &
                      + .500_r8*rxt(k,486)*y(k,192)
         mat(k,1715) = mat(k,1715) + rxt(k,397)*y(k,1) + rxt(k,191)*y(k,81) &
                      + .700_r8*rxt(k,417)*y(k,90) + rxt(k,205)*y(k,113) + rxt(k,358) &
                      *y(k,114) + rxt(k,495)*y(k,147)
         mat(k,339) = rxt(k,461)*y(k,111)
         mat(k,713) = rxt(k,360)*y(k,111)
         mat(k,1051) = rxt(k,364)*y(k,111)
         mat(k,1013) = .900_r8*rxt(k,493)*y(k,111)
         mat(k,825) = .800_r8*rxt(k,498)*y(k,111)
         mat(k,609) = rxt(k,468)*y(k,111)
         mat(k,1092) = rxt(k,434)*y(k,111) + rxt(k,435)*y(k,113)
         mat(k,641) = rxt(k,474)*y(k,111)
         mat(k,403) = rxt(k,477)*y(k,111)
         mat(k,1961) = -(rxt(k,202)*y(k,198) + rxt(k,203)*y(k,111) + rxt(k,204) &
                      *y(k,116) + rxt(k,205)*y(k,216) + rxt(k,213)*y(k,112) + rxt(k,299) &
                      *y(k,38) + rxt(k,329)*y(k,41) + rxt(k,348)*y(k,25) + rxt(k,355) &
                      *y(k,45) + rxt(k,368)*y(k,13) + rxt(k,376)*y(k,202) + rxt(k,387) &
                      *y(k,204) + rxt(k,410)*y(k,200) + rxt(k,416)*y(k,201) + rxt(k,419) &
                      *y(k,89) + rxt(k,424)*y(k,92) + rxt(k,435)*y(k,224) + rxt(k,480) &
                      *y(k,3) + rxt(k,483)*y(k,101) + rxt(k,489)*y(k,211) + rxt(k,500) &
                      *y(k,149) + rxt(k,517)*y(k,62))
         mat(k,2055) = -rxt(k,202)*y(k,113)
         mat(k,1500) = -rxt(k,203)*y(k,113)
         mat(k,1565) = -rxt(k,204)*y(k,113)
         mat(k,1713) = -rxt(k,205)*y(k,113)
         mat(k,2097) = -rxt(k,213)*y(k,113)
         mat(k,2121) = -rxt(k,299)*y(k,113)
         mat(k,1061) = -rxt(k,329)*y(k,113)
         mat(k,954) = -rxt(k,348)*y(k,113)
         mat(k,1137) = -rxt(k,355)*y(k,113)
         mat(k,274) = -rxt(k,368)*y(k,113)
         mat(k,1219) = -rxt(k,376)*y(k,113)
         mat(k,1259) = -rxt(k,387)*y(k,113)
         mat(k,1163) = -rxt(k,410)*y(k,113)
         mat(k,1197) = -rxt(k,416)*y(k,113)
         mat(k,793) = -rxt(k,419)*y(k,113)
         mat(k,1114) = -rxt(k,424)*y(k,113)
         mat(k,1090) = -rxt(k,435)*y(k,113)
         mat(k,904) = -rxt(k,480)*y(k,113)
         mat(k,863) = -rxt(k,483)*y(k,113)
         mat(k,970) = -rxt(k,489)*y(k,113)
         mat(k,915) = -rxt(k,500)*y(k,113)
         mat(k,210) = -rxt(k,517)*y(k,113)
         mat(k,446) = rxt(k,264)*y(k,116)
         mat(k,1843) = rxt(k,231)*y(k,56)
         mat(k,877) = rxt(k,231)*y(k,52) + rxt(k,233)*y(k,116) + rxt(k,234)*y(k,216)
         mat(k,666) = rxt(k,278)*y(k,80)
         mat(k,1788) = rxt(k,278)*y(k,64) + rxt(k,215)*y(k,216)
         mat(k,452) = .500_r8*rxt(k,392)*y(k,216)
         mat(k,2097) = mat(k,2097) + rxt(k,201)*y(k,116) + rxt(k,200)*y(k,118)
         mat(k,1565) = mat(k,1565) + rxt(k,264)*y(k,17) + rxt(k,233)*y(k,56) &
                      + rxt(k,201)*y(k,112)
         mat(k,1904) = rxt(k,200)*y(k,112)
         mat(k,352) = rxt(k,344)*y(k,216)
         mat(k,1713) = mat(k,1713) + rxt(k,234)*y(k,56) + rxt(k,215)*y(k,80) &
                      + .500_r8*rxt(k,392)*y(k,100) + rxt(k,344)*y(k,123)
         mat(k,797) = -(rxt(k,358)*y(k,216))
         mat(k,1665) = -rxt(k,358)*y(k,114)
         mat(k,942) = rxt(k,348)*y(k,113)
         mat(k,374) = .500_r8*rxt(k,418)*y(k,216)
         mat(k,279) = rxt(k,425)*y(k,216)
         mat(k,223) = rxt(k,429)*y(k,216)
         mat(k,752) = rxt(k,430)*y(k,216)
         mat(k,1918) = rxt(k,348)*y(k,25)
         mat(k,1665) = mat(k,1665) + .500_r8*rxt(k,418)*y(k,91) + rxt(k,425)*y(k,93) &
                      + rxt(k,429)*y(k,106) + rxt(k,430)*y(k,107)
         mat(k,227) = -(rxt(k,490)*y(k,216))
         mat(k,1606) = -rxt(k,490)*y(k,115)
         mat(k,1971) = rxt(k,487)*y(k,211)
         mat(k,959) = rxt(k,487)*y(k,198)
         mat(k,1557) = -(rxt(k,171)*y(k,118) + 4._r8*rxt(k,172)*y(k,116) + rxt(k,173) &
                      *y(k,117) + rxt(k,174)*y(k,68) + rxt(k,175)*y(k,70) + rxt(k,180) &
                      *y(k,198) + rxt(k,186)*y(k,216) + (rxt(k,199) + rxt(k,201) &
                      ) * y(k,112) + rxt(k,204)*y(k,113) + rxt(k,209)*y(k,111) &
                      + rxt(k,233)*y(k,56) + rxt(k,235)*y(k,55) + rxt(k,238)*y(k,76) &
                      + rxt(k,241)*y(k,83) + rxt(k,264)*y(k,17) + rxt(k,265)*y(k,16) &
                      + rxt(k,267)*y(k,72) + rxt(k,269)*y(k,82) + rxt(k,300)*y(k,38) &
                      + rxt(k,503)*y(k,121) + (rxt(k,563) + rxt(k,564)) * y(k,208) &
                      + rxt(k,565)*y(k,210))
         mat(k,1896) = -rxt(k,171)*y(k,116)
         mat(k,1413) = -rxt(k,173)*y(k,116)
         mat(k,1069) = -rxt(k,174)*y(k,116)
         mat(k,499) = -rxt(k,175)*y(k,116)
         mat(k,2047) = -rxt(k,180)*y(k,116)
         mat(k,1705) = -rxt(k,186)*y(k,116)
         mat(k,2089) = -(rxt(k,199) + rxt(k,201)) * y(k,116)
         mat(k,1953) = -rxt(k,204)*y(k,116)
         mat(k,1492) = -rxt(k,209)*y(k,116)
         mat(k,872) = -rxt(k,233)*y(k,116)
         mat(k,1732) = -rxt(k,235)*y(k,116)
         mat(k,1322) = -rxt(k,238)*y(k,116)
         mat(k,770) = -rxt(k,241)*y(k,116)
         mat(k,444) = -rxt(k,264)*y(k,116)
         mat(k,1516) = -rxt(k,265)*y(k,116)
         mat(k,762) = -rxt(k,267)*y(k,116)
         mat(k,655) = -rxt(k,269)*y(k,116)
         mat(k,2113) = -rxt(k,300)*y(k,116)
         mat(k,266) = -rxt(k,503)*y(k,116)
         mat(k,572) = -(rxt(k,563) + rxt(k,564)) * y(k,116)
         mat(k,395) = -rxt(k,565)*y(k,116)
         mat(k,1800) = rxt(k,178)*y(k,198)
         mat(k,745) = rxt(k,194)*y(k,111) + rxt(k,195)*y(k,112) + rxt(k,198)*y(k,117) &
                      + rxt(k,568)*y(k,215)
         mat(k,1492) = mat(k,1492) + rxt(k,194)*y(k,103)
         mat(k,2089) = mat(k,2089) + rxt(k,195)*y(k,103)
         mat(k,1413) = mat(k,1413) + rxt(k,198)*y(k,103) + rxt(k,505)*y(k,130) &
                      + rxt(k,511)*y(k,131) + rxt(k,567)*y(k,210) + (rxt(k,160) &
                       +rxt(k,161))*y(k,212) + rxt(k,573)*y(k,217)
         mat(k,613) = rxt(k,505)*y(k,117)
         mat(k,1305) = rxt(k,511)*y(k,117)
         mat(k,703) = rxt(k,559)*y(k,209) + 1.150_r8*rxt(k,560)*y(k,215)
         mat(k,2047) = mat(k,2047) + rxt(k,178)*y(k,67)
         mat(k,692) = rxt(k,559)*y(k,194)
         mat(k,395) = mat(k,395) + rxt(k,567)*y(k,117)
         mat(k,1758) = (rxt(k,160)+rxt(k,161))*y(k,117)
         mat(k,684) = rxt(k,568)*y(k,103) + 1.150_r8*rxt(k,560)*y(k,194)
         mat(k,1705) = mat(k,1705) + 2.000_r8*rxt(k,188)*y(k,216)
         mat(k,510) = rxt(k,573)*y(k,117)
         mat(k,1410) = -(rxt(k,160)*y(k,212) + rxt(k,165)*y(k,213) + rxt(k,173) &
                      *y(k,116) + rxt(k,179)*y(k,67) + rxt(k,193)*y(k,207) + rxt(k,198) &
                      *y(k,103) + rxt(k,341)*y(k,196) + rxt(k,505)*y(k,130) + rxt(k,511) &
                      *y(k,131) + rxt(k,562)*y(k,208) + (rxt(k,566) + rxt(k,567) &
                      ) * y(k,210) + rxt(k,573)*y(k,217))
         mat(k,1755) = -rxt(k,160)*y(k,117)
         mat(k,78) = -rxt(k,165)*y(k,117)
         mat(k,1554) = -rxt(k,173)*y(k,117)
         mat(k,1797) = -rxt(k,179)*y(k,117)
         mat(k,423) = -rxt(k,193)*y(k,117)
         mat(k,743) = -rxt(k,198)*y(k,117)
         mat(k,345) = -rxt(k,341)*y(k,117)
         mat(k,612) = -rxt(k,505)*y(k,117)
         mat(k,1302) = -rxt(k,511)*y(k,117)
         mat(k,570) = -rxt(k,562)*y(k,117)
         mat(k,394) = -(rxt(k,566) + rxt(k,567)) * y(k,117)
         mat(k,509) = -rxt(k,573)*y(k,117)
         mat(k,1335) = rxt(k,256)*y(k,118) + rxt(k,255)*y(k,198)
         mat(k,1513) = 2.000_r8*rxt(k,257)*y(k,16) + (rxt(k,259)+rxt(k,260))*y(k,55) &
                      + rxt(k,265)*y(k,116) + rxt(k,261)*y(k,198)
         mat(k,1832) = rxt(k,224)*y(k,118) + rxt(k,222)*y(k,198)
         mat(k,1729) = (rxt(k,259)+rxt(k,260))*y(k,16) + (2.000_r8*rxt(k,226) &
                       +2.000_r8*rxt(k,227))*y(k,55) + rxt(k,235)*y(k,116) &
                      + rxt(k,229)*y(k,198) + rxt(k,237)*y(k,216)
         mat(k,1797) = mat(k,1797) + rxt(k,182)*y(k,118) + rxt(k,176)*y(k,198)
         mat(k,360) = rxt(k,191)*y(k,216)
         mat(k,743) = mat(k,743) + rxt(k,197)*y(k,112)
         mat(k,1489) = rxt(k,208)*y(k,118) + rxt(k,570)*y(k,215)
         mat(k,2086) = rxt(k,197)*y(k,103) + rxt(k,199)*y(k,116) + rxt(k,200)*y(k,118)
         mat(k,1950) = rxt(k,204)*y(k,116) + rxt(k,202)*y(k,198)
         mat(k,1554) = mat(k,1554) + rxt(k,265)*y(k,16) + rxt(k,235)*y(k,55) &
                      + rxt(k,199)*y(k,112) + rxt(k,204)*y(k,113) &
                      + 2.000_r8*rxt(k,172)*y(k,116) + 2.000_r8*rxt(k,171)*y(k,118) &
                      + rxt(k,180)*y(k,198) + rxt(k,164)*y(k,213) + rxt(k,186) &
                      *y(k,216)
         mat(k,1410) = mat(k,1410) + 2.000_r8*rxt(k,165)*y(k,213)
         mat(k,1893) = rxt(k,256)*y(k,14) + rxt(k,224)*y(k,52) + rxt(k,182)*y(k,67) &
                      + rxt(k,208)*y(k,111) + rxt(k,200)*y(k,112) &
                      + 2.000_r8*rxt(k,171)*y(k,116) + rxt(k,506)*y(k,130) &
                      + rxt(k,512)*y(k,131) + 2.000_r8*rxt(k,181)*y(k,198) &
                      + 2.000_r8*rxt(k,162)*y(k,212) + rxt(k,187)*y(k,216)
         mat(k,612) = mat(k,612) + rxt(k,506)*y(k,118)
         mat(k,1302) = mat(k,1302) + rxt(k,512)*y(k,118)
         mat(k,809) = rxt(k,323)*y(k,198)
         mat(k,729) = rxt(k,351)*y(k,198)
         mat(k,1380) = rxt(k,304)*y(k,198)
         mat(k,2044) = rxt(k,255)*y(k,14) + rxt(k,261)*y(k,16) + rxt(k,222)*y(k,52) &
                      + rxt(k,229)*y(k,55) + rxt(k,176)*y(k,67) + rxt(k,202)*y(k,113) &
                      + rxt(k,180)*y(k,116) + 2.000_r8*rxt(k,181)*y(k,118) &
                      + rxt(k,323)*y(k,188) + rxt(k,351)*y(k,189) + rxt(k,304) &
                      *y(k,192) + 2.000_r8*rxt(k,190)*y(k,198) + rxt(k,185)*y(k,216) &
                      + rxt(k,359)*y(k,219)
         mat(k,1755) = mat(k,1755) + 2.000_r8*rxt(k,162)*y(k,118)
         mat(k,78) = mat(k,78) + rxt(k,164)*y(k,116) + 2.000_r8*rxt(k,165)*y(k,117)
         mat(k,682) = rxt(k,570)*y(k,111)
         mat(k,1702) = rxt(k,237)*y(k,55) + rxt(k,191)*y(k,81) + rxt(k,186)*y(k,116) &
                      + rxt(k,187)*y(k,118) + rxt(k,185)*y(k,198)
         mat(k,709) = rxt(k,359)*y(k,198)
      end do
      end subroutine nlnmat05
      subroutine nlnmat06( ofl, ofu, chnkpnts, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = ofl,ofu
         mat(k,1903) = -(rxt(k,162)*y(k,212) + rxt(k,171)*y(k,116) + rxt(k,181) &
                      *y(k,198) + rxt(k,182)*y(k,67) + rxt(k,187)*y(k,216) + rxt(k,200) &
                      *y(k,112) + rxt(k,208)*y(k,111) + rxt(k,224)*y(k,52) + rxt(k,256) &
                      *y(k,14) + rxt(k,320)*y(k,21) + rxt(k,349)*y(k,25) + rxt(k,379) &
                      *y(k,96) + rxt(k,393)*y(k,102) + rxt(k,426)*y(k,89) + rxt(k,464) &
                      *y(k,125) + rxt(k,481)*y(k,3) + rxt(k,484)*y(k,101) + rxt(k,506) &
                      *y(k,130) + rxt(k,512)*y(k,131))
         mat(k,1765) = -rxt(k,162)*y(k,118)
         mat(k,1564) = -rxt(k,171)*y(k,118)
         mat(k,2054) = -rxt(k,181)*y(k,118)
         mat(k,1807) = -rxt(k,182)*y(k,118)
         mat(k,1712) = -rxt(k,187)*y(k,118)
         mat(k,2096) = -rxt(k,200)*y(k,118)
         mat(k,1499) = -rxt(k,208)*y(k,118)
         mat(k,1842) = -rxt(k,224)*y(k,118)
         mat(k,1341) = -rxt(k,256)*y(k,118)
         mat(k,469) = -rxt(k,320)*y(k,118)
         mat(k,953) = -rxt(k,349)*y(k,118)
         mat(k,1127) = -rxt(k,379)*y(k,118)
         mat(k,1240) = -rxt(k,393)*y(k,118)
         mat(k,792) = -rxt(k,426)*y(k,118)
         mat(k,357) = -rxt(k,464)*y(k,118)
         mat(k,903) = -rxt(k,481)*y(k,118)
         mat(k,862) = -rxt(k,484)*y(k,118)
         mat(k,616) = -rxt(k,506)*y(k,118)
         mat(k,1310) = -rxt(k,512)*y(k,118)
         mat(k,1564) = mat(k,1564) + rxt(k,173)*y(k,117)
         mat(k,1419) = rxt(k,173)*y(k,116)
         mat(k,1291) = .150_r8*rxt(k,334)*y(k,198)
         mat(k,2054) = mat(k,2054) + .150_r8*rxt(k,334)*y(k,191) + .150_r8*rxt(k,384) &
                      *y(k,204)
         mat(k,1258) = .150_r8*rxt(k,384)*y(k,198)
         mat(k,232) = -(rxt(k,513)*y(k,131))
         mat(k,1297) = -rxt(k,513)*y(k,120)
         mat(k,1506) = rxt(k,258)*y(k,55)
         mat(k,1721) = rxt(k,258)*y(k,16) + 2.000_r8*rxt(k,228)*y(k,55)
         mat(k,261) = -(rxt(k,503)*y(k,116) + rxt(k,504)*y(k,216))
         mat(k,1531) = -rxt(k,503)*y(k,121)
         mat(k,1611) = -rxt(k,504)*y(k,121)
         mat(k,1016) = rxt(k,372)*y(k,216)
         mat(k,1424) = .100_r8*rxt(k,493)*y(k,221)
         mat(k,1582) = rxt(k,372)*y(k,84)
         mat(k,999) = .100_r8*rxt(k,493)*y(k,111)
         mat(k,348) = -(rxt(k,344)*y(k,216))
         mat(k,1622) = -rxt(k,344)*y(k,123)
         mat(k,2062) = rxt(k,346)*y(k,191)
         mat(k,1265) = rxt(k,346)*y(k,112)
         mat(k,2060) = rxt(k,466)*y(k,184)
         mat(k,414) = rxt(k,466)*y(k,112)
         mat(k,355) = -(rxt(k,463)*y(k,112) + rxt(k,464)*y(k,118))
         mat(k,2063) = -rxt(k,463)*y(k,125)
         mat(k,1857) = -rxt(k,464)*y(k,125)
         mat(k,126) = .070_r8*rxt(k,450)*y(k,216)
         mat(k,1434) = rxt(k,448)*y(k,190)
         mat(k,99) = .060_r8*rxt(k,462)*y(k,216)
         mat(k,151) = .070_r8*rxt(k,478)*y(k,216)
         mat(k,512) = rxt(k,448)*y(k,111)
         mat(k,1623) = .070_r8*rxt(k,450)*y(k,61) + .060_r8*rxt(k,462)*y(k,126) &
                      + .070_r8*rxt(k,478)*y(k,156)
         mat(k,97) = -(rxt(k,462)*y(k,216))
         mat(k,1585) = -rxt(k,462)*y(k,126)
         mat(k,89) = .530_r8*rxt(k,439)*y(k,216)
         mat(k,1585) = mat(k,1585) + .530_r8*rxt(k,439)*y(k,4)
         mat(k,237) = -(rxt(k,465)*y(k,216))
         mat(k,1607) = -rxt(k,465)*y(k,127)
         mat(k,1972) = rxt(k,460)*y(k,218)
         mat(k,333) = rxt(k,460)*y(k,198)
         mat(k,432) = -(rxt(k,361)*y(k,216))
         mat(k,1633) = -rxt(k,361)*y(k,128)
         mat(k,1993) = rxt(k,359)*y(k,219)
         mat(k,705) = rxt(k,359)*y(k,198)
         mat(k,283) = -(rxt(k,365)*y(k,216))
         mat(k,1614) = -rxt(k,365)*y(k,129)
         mat(k,1976) = .850_r8*rxt(k,363)*y(k,220)
         mat(k,1041) = .850_r8*rxt(k,363)*y(k,198)
         mat(k,610) = -(rxt(k,505)*y(k,117) + rxt(k,506)*y(k,118) + rxt(k,509) &
                      *y(k,216))
         mat(k,1400) = -rxt(k,505)*y(k,130)
         mat(k,1861) = -rxt(k,506)*y(k,130)
         mat(k,1650) = -rxt(k,509)*y(k,130)
         mat(k,1300) = -(rxt(k,507)*y(k,16) + rxt(k,508)*y(k,55) + rxt(k,510)*y(k,112) &
                      + rxt(k,511)*y(k,117) + rxt(k,512)*y(k,118) + rxt(k,513) &
                      *y(k,120) + rxt(k,514)*y(k,216))
         mat(k,1510) = -rxt(k,507)*y(k,131)
         mat(k,1725) = -rxt(k,508)*y(k,131)
         mat(k,2082) = -rxt(k,510)*y(k,131)
         mat(k,1408) = -rxt(k,511)*y(k,131)
         mat(k,1890) = -rxt(k,512)*y(k,131)
         mat(k,234) = -rxt(k,513)*y(k,131)
         mat(k,1698) = -rxt(k,514)*y(k,131)
         mat(k,1550) = rxt(k,503)*y(k,121)
         mat(k,1408) = mat(k,1408) + rxt(k,505)*y(k,130)
         mat(k,1890) = mat(k,1890) + rxt(k,506)*y(k,130)
         mat(k,265) = rxt(k,503)*y(k,116)
         mat(k,611) = rxt(k,505)*y(k,117) + rxt(k,506)*y(k,118) + rxt(k,509)*y(k,216)
         mat(k,1698) = mat(k,1698) + rxt(k,509)*y(k,130)
         mat(k,833) = -(rxt(k,515)*y(k,216))
         mat(k,1669) = -rxt(k,515)*y(k,132)
         mat(k,1509) = rxt(k,507)*y(k,131)
         mat(k,1723) = rxt(k,508)*y(k,131)
         mat(k,207) = rxt(k,517)*y(k,113) + (rxt(k,518)+.500_r8*rxt(k,520))*y(k,216)
         mat(k,2074) = rxt(k,510)*y(k,131)
         mat(k,1920) = rxt(k,517)*y(k,62)
         mat(k,1405) = rxt(k,511)*y(k,131)
         mat(k,1866) = rxt(k,512)*y(k,131)
         mat(k,233) = rxt(k,513)*y(k,131)
         mat(k,263) = rxt(k,504)*y(k,216)
         mat(k,1299) = rxt(k,507)*y(k,16) + rxt(k,508)*y(k,55) + rxt(k,510)*y(k,112) &
                      + rxt(k,511)*y(k,117) + rxt(k,512)*y(k,118) + rxt(k,513) &
                      *y(k,120) + rxt(k,514)*y(k,216)
         mat(k,1669) = mat(k,1669) + (rxt(k,518)+.500_r8*rxt(k,520))*y(k,62) &
                      + rxt(k,504)*y(k,121) + rxt(k,514)*y(k,131)
         mat(k,173) = -(rxt(k,516)*y(k,227))
         mat(k,2128) = -rxt(k,516)*y(k,133)
         mat(k,832) = rxt(k,515)*y(k,216)
         mat(k,1597) = rxt(k,515)*y(k,132)
         mat(k,46) = -(rxt(k,539)*y(k,216))
         mat(k,1576) = -rxt(k,539)*y(k,144)
         mat(k,119) = .100_r8*rxt(k,470)*y(k,216)
         mat(k,141) = .230_r8*rxt(k,472)*y(k,216)
         mat(k,1590) = .100_r8*rxt(k,470)*y(k,152) + .230_r8*rxt(k,472)*y(k,154)
         mat(k,488) = -(rxt(k,494)*y(k,216))
         mat(k,1640) = -rxt(k,494)*y(k,146)
         mat(k,1995) = rxt(k,492)*y(k,221)
         mat(k,1000) = rxt(k,492)*y(k,198)
         mat(k,409) = -(rxt(k,495)*y(k,216))
         mat(k,1630) = -rxt(k,495)*y(k,147)
         mat(k,1439) = .200_r8*rxt(k,488)*y(k,211) + .200_r8*rxt(k,498)*y(k,222)
         mat(k,1349) = .500_r8*rxt(k,486)*y(k,211)
         mat(k,960) = .200_r8*rxt(k,488)*y(k,111) + .500_r8*rxt(k,486)*y(k,192)
         mat(k,816) = .200_r8*rxt(k,498)*y(k,111)
         mat(k,366) = -(rxt(k,499)*y(k,216))
         mat(k,1625) = -rxt(k,499)*y(k,148)
         mat(k,1988) = rxt(k,497)*y(k,222)
         mat(k,815) = rxt(k,497)*y(k,198)
         mat(k,909) = -(rxt(k,500)*y(k,113) + rxt(k,501)*y(k,216))
         mat(k,1923) = -rxt(k,500)*y(k,149)
         mat(k,1673) = -rxt(k,501)*y(k,149)
         mat(k,891) = .330_r8*rxt(k,481)*y(k,118)
         mat(k,850) = .330_r8*rxt(k,484)*y(k,118)
         mat(k,1464) = .800_r8*rxt(k,488)*y(k,211) + .800_r8*rxt(k,498)*y(k,222)
         mat(k,1923) = mat(k,1923) + rxt(k,489)*y(k,211)
         mat(k,1869) = .330_r8*rxt(k,481)*y(k,3) + .330_r8*rxt(k,484)*y(k,101)
         mat(k,410) = rxt(k,495)*y(k,216)
         mat(k,1358) = .500_r8*rxt(k,486)*y(k,211) + rxt(k,496)*y(k,222)
         mat(k,962) = .800_r8*rxt(k,488)*y(k,111) + rxt(k,489)*y(k,113) &
                      + .500_r8*rxt(k,486)*y(k,192)
         mat(k,1673) = mat(k,1673) + rxt(k,495)*y(k,147)
         mat(k,819) = .800_r8*rxt(k,498)*y(k,111) + rxt(k,496)*y(k,192)
         mat(k,975) = -(rxt(k,502)*y(k,216))
         mat(k,1678) = -rxt(k,502)*y(k,150)
         mat(k,895) = .300_r8*rxt(k,481)*y(k,118)
         mat(k,854) = .300_r8*rxt(k,484)*y(k,118)
         mat(k,1468) = .900_r8*rxt(k,493)*y(k,221)
         mat(k,1874) = .300_r8*rxt(k,481)*y(k,3) + .300_r8*rxt(k,484)*y(k,101)
         mat(k,1361) = rxt(k,491)*y(k,221)
         mat(k,1003) = .900_r8*rxt(k,493)*y(k,111) + rxt(k,491)*y(k,192)
         mat(k,522) = -(rxt(k,469)*y(k,216))
         mat(k,1643) = -rxt(k,469)*y(k,151)
         mat(k,1998) = rxt(k,467)*y(k,223)
         mat(k,598) = rxt(k,467)*y(k,198)
         mat(k,117) = -(rxt(k,470)*y(k,216))
         mat(k,1588) = -rxt(k,470)*y(k,152)
         mat(k,133) = -(rxt(k,436)*y(k,216))
         mat(k,1591) = -rxt(k,436)*y(k,153)
         mat(k,1967) = rxt(k,433)*y(k,224)
         mat(k,1077) = rxt(k,433)*y(k,198)
         mat(k,142) = -(rxt(k,472)*y(k,216))
         mat(k,1592) = -rxt(k,472)*y(k,154)
         mat(k,578) = -(rxt(k,475)*y(k,216))
         mat(k,1647) = -rxt(k,475)*y(k,155)
         mat(k,2001) = rxt(k,473)*y(k,225)
         mat(k,629) = rxt(k,473)*y(k,198)
         mat(k,150) = -(rxt(k,478)*y(k,216))
         mat(k,1593) = -rxt(k,478)*y(k,156)
         mat(k,143) = .150_r8*rxt(k,472)*y(k,216)
         mat(k,1593) = mat(k,1593) + .150_r8*rxt(k,472)*y(k,154)
         mat(k,301) = -(rxt(k,479)*y(k,216))
         mat(k,1616) = -rxt(k,479)*y(k,157)
         mat(k,1978) = rxt(k,476)*y(k,226)
         mat(k,396) = rxt(k,476)*y(k,198)
         mat(k,880) = .2202005_r8*rxt(k,529)*y(k,118) + .2202005_r8*rxt(k,530) &
                      *y(k,216)
         mat(k,82) = .0023005_r8*rxt(k,531)*y(k,216)
         mat(k,775) = .0031005_r8*rxt(k,534)*y(k,216)
         mat(k,35) = .2381005_r8*rxt(k,535)*y(k,216)
         mat(k,839) = .0508005_r8*rxt(k,537)*y(k,118) + .0508005_r8*rxt(k,538) &
                      *y(k,216)
         mat(k,1848) = .2202005_r8*rxt(k,529)*y(k,3) + .0508005_r8*rxt(k,537)*y(k,101)
         mat(k,41) = .5931005_r8*rxt(k,539)*y(k,216)
         mat(k,112) = .1364005_r8*rxt(k,540)*y(k,216)
         mat(k,136) = .1677005_r8*rxt(k,541)*y(k,216)
         mat(k,1570) = .2202005_r8*rxt(k,530)*y(k,3) + .0023005_r8*rxt(k,531)*y(k,4) &
                      + .0031005_r8*rxt(k,534)*y(k,89) + .2381005_r8*rxt(k,535) &
                      *y(k,95) + .0508005_r8*rxt(k,538)*y(k,101) &
                      + .5931005_r8*rxt(k,539)*y(k,144) + .1364005_r8*rxt(k,540) &
                      *y(k,152) + .1677005_r8*rxt(k,541)*y(k,154)
         mat(k,881) = .2067005_r8*rxt(k,529)*y(k,118) + .2067005_r8*rxt(k,530) &
                      *y(k,216)
         mat(k,83) = .0008005_r8*rxt(k,531)*y(k,216)
         mat(k,776) = .0035005_r8*rxt(k,534)*y(k,216)
         mat(k,36) = .1308005_r8*rxt(k,535)*y(k,216)
         mat(k,840) = .1149005_r8*rxt(k,537)*y(k,118) + .1149005_r8*rxt(k,538) &
                      *y(k,216)
         mat(k,1849) = .2067005_r8*rxt(k,529)*y(k,3) + .1149005_r8*rxt(k,537)*y(k,101)
         mat(k,42) = .1534005_r8*rxt(k,539)*y(k,216)
         mat(k,113) = .0101005_r8*rxt(k,540)*y(k,216)
         mat(k,137) = .0174005_r8*rxt(k,541)*y(k,216)
         mat(k,1571) = .2067005_r8*rxt(k,530)*y(k,3) + .0008005_r8*rxt(k,531)*y(k,4) &
                      + .0035005_r8*rxt(k,534)*y(k,89) + .1308005_r8*rxt(k,535) &
                      *y(k,95) + .1149005_r8*rxt(k,538)*y(k,101) &
                      + .1534005_r8*rxt(k,539)*y(k,144) + .0101005_r8*rxt(k,540) &
                      *y(k,152) + .0174005_r8*rxt(k,541)*y(k,154)
      end do
      end subroutine nlnmat06
      subroutine nlnmat07( ofl, ofu, chnkpnts, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = ofl,ofu
         mat(k,882) = .0653005_r8*rxt(k,529)*y(k,118) + .0653005_r8*rxt(k,530) &
                      *y(k,216)
         mat(k,84) = .0843005_r8*rxt(k,531)*y(k,216)
         mat(k,777) = .0003005_r8*rxt(k,534)*y(k,216)
         mat(k,37) = .0348005_r8*rxt(k,535)*y(k,216)
         mat(k,841) = .0348005_r8*rxt(k,537)*y(k,118) + .0348005_r8*rxt(k,538) &
                      *y(k,216)
         mat(k,1850) = .0653005_r8*rxt(k,529)*y(k,3) + .0348005_r8*rxt(k,537)*y(k,101)
         mat(k,43) = .0459005_r8*rxt(k,539)*y(k,216)
         mat(k,114) = .0763005_r8*rxt(k,540)*y(k,216)
         mat(k,138) = .086_r8*rxt(k,541)*y(k,216)
         mat(k,1572) = .0653005_r8*rxt(k,530)*y(k,3) + .0843005_r8*rxt(k,531)*y(k,4) &
                      + .0003005_r8*rxt(k,534)*y(k,89) + .0348005_r8*rxt(k,535) &
                      *y(k,95) + .0348005_r8*rxt(k,538)*y(k,101) &
                      + .0459005_r8*rxt(k,539)*y(k,144) + .0763005_r8*rxt(k,540) &
                      *y(k,152) + .086_r8*rxt(k,541)*y(k,154)
         mat(k,883) = .1749305_r8*rxt(k,528)*y(k,113) + .1284005_r8*rxt(k,529) &
                      *y(k,118) + .1284005_r8*rxt(k,530)*y(k,216)
         mat(k,85) = .0443005_r8*rxt(k,531)*y(k,216)
         mat(k,778) = .0590245_r8*rxt(k,532)*y(k,113) + .0033005_r8*rxt(k,533) &
                      *y(k,118) + .0271005_r8*rxt(k,534)*y(k,216)
         mat(k,38) = .0076005_r8*rxt(k,535)*y(k,216)
         mat(k,842) = .1749305_r8*rxt(k,536)*y(k,113) + .0554005_r8*rxt(k,537) &
                      *y(k,118) + .0554005_r8*rxt(k,538)*y(k,216)
         mat(k,1909) = .1749305_r8*rxt(k,528)*y(k,3) + .0590245_r8*rxt(k,532)*y(k,89) &
                      + .1749305_r8*rxt(k,536)*y(k,101)
         mat(k,1851) = .1284005_r8*rxt(k,529)*y(k,3) + .0033005_r8*rxt(k,533)*y(k,89) &
                      + .0554005_r8*rxt(k,537)*y(k,101)
         mat(k,44) = .0085005_r8*rxt(k,539)*y(k,216)
         mat(k,115) = .2157005_r8*rxt(k,540)*y(k,216)
         mat(k,139) = .0512005_r8*rxt(k,541)*y(k,216)
         mat(k,1573) = .1284005_r8*rxt(k,530)*y(k,3) + .0443005_r8*rxt(k,531)*y(k,4) &
                      + .0271005_r8*rxt(k,534)*y(k,89) + .0076005_r8*rxt(k,535) &
                      *y(k,95) + .0554005_r8*rxt(k,538)*y(k,101) &
                      + .0085005_r8*rxt(k,539)*y(k,144) + .2157005_r8*rxt(k,540) &
                      *y(k,152) + .0512005_r8*rxt(k,541)*y(k,154)
         mat(k,884) = .5901905_r8*rxt(k,528)*y(k,113) + .114_r8*rxt(k,529)*y(k,118) &
                      + .114_r8*rxt(k,530)*y(k,216)
         mat(k,86) = .1621005_r8*rxt(k,531)*y(k,216)
         mat(k,779) = .0250245_r8*rxt(k,532)*y(k,113) + .0474005_r8*rxt(k,534) &
                      *y(k,216)
         mat(k,39) = .0113005_r8*rxt(k,535)*y(k,216)
         mat(k,843) = .5901905_r8*rxt(k,536)*y(k,113) + .1278005_r8*rxt(k,537) &
                      *y(k,118) + .1278005_r8*rxt(k,538)*y(k,216)
         mat(k,1910) = .5901905_r8*rxt(k,528)*y(k,3) + .0250245_r8*rxt(k,532)*y(k,89) &
                      + .5901905_r8*rxt(k,536)*y(k,101)
         mat(k,1852) = .114_r8*rxt(k,529)*y(k,3) + .1278005_r8*rxt(k,537)*y(k,101)
         mat(k,45) = .0128005_r8*rxt(k,539)*y(k,216)
         mat(k,116) = .0232005_r8*rxt(k,540)*y(k,216)
         mat(k,140) = .1598005_r8*rxt(k,541)*y(k,216)
         mat(k,1574) = .114_r8*rxt(k,530)*y(k,3) + .1621005_r8*rxt(k,531)*y(k,4) &
                      + .0474005_r8*rxt(k,534)*y(k,89) + .0113005_r8*rxt(k,535) &
                      *y(k,95) + .1278005_r8*rxt(k,538)*y(k,101) &
                      + .0128005_r8*rxt(k,539)*y(k,144) + .0232005_r8*rxt(k,540) &
                      *y(k,152) + .1598005_r8*rxt(k,541)*y(k,154)
         mat(k,415) = -(rxt(k,437)*y(k,198) + rxt(k,438)*y(k,111) + rxt(k,466) &
                      *y(k,112))
         mat(k,1991) = -rxt(k,437)*y(k,184)
         mat(k,1440) = -rxt(k,438)*y(k,184)
         mat(k,2065) = -rxt(k,466)*y(k,184)
         mat(k,170) = rxt(k,443)*y(k,216)
         mat(k,1631) = rxt(k,443)*y(k,18)
         mat(k,922) = -(rxt(k,398)*y(k,198) + (rxt(k,399) + rxt(k,400)) * y(k,111))
         mat(k,2019) = -rxt(k,398)*y(k,185)
         mat(k,1465) = -(rxt(k,399) + rxt(k,400)) * y(k,185)
         mat(k,536) = rxt(k,401)*y(k,216)
         mat(k,161) = rxt(k,402)*y(k,216)
         mat(k,1674) = rxt(k,401)*y(k,2) + rxt(k,402)*y(k,12)
         mat(k,382) = -(rxt(k,440)*y(k,198) + rxt(k,441)*y(k,111))
         mat(k,1989) = -rxt(k,440)*y(k,186)
         mat(k,1436) = -rxt(k,441)*y(k,186)
         mat(k,90) = .350_r8*rxt(k,439)*y(k,216)
         mat(k,291) = rxt(k,442)*y(k,216)
         mat(k,1627) = .350_r8*rxt(k,439)*y(k,4) + rxt(k,442)*y(k,5)
         mat(k,321) = -(rxt(k,444)*y(k,198) + rxt(k,446)*y(k,111))
         mat(k,1981) = -rxt(k,444)*y(k,187)
         mat(k,1429) = -rxt(k,446)*y(k,187)
         mat(k,249) = rxt(k,445)*y(k,216)
         mat(k,120) = .070_r8*rxt(k,470)*y(k,216)
         mat(k,144) = .060_r8*rxt(k,472)*y(k,216)
         mat(k,1619) = rxt(k,445)*y(k,19) + .070_r8*rxt(k,470)*y(k,152) &
                      + .060_r8*rxt(k,472)*y(k,154)
         mat(k,806) = -(4._r8*rxt(k,321)*y(k,188) + rxt(k,322)*y(k,192) + rxt(k,323) &
                      *y(k,198) + rxt(k,324)*y(k,111))
         mat(k,1355) = -rxt(k,322)*y(k,188)
         mat(k,2015) = -rxt(k,323)*y(k,188)
         mat(k,1461) = -rxt(k,324)*y(k,188)
         mat(k,254) = .500_r8*rxt(k,326)*y(k,216)
         mat(k,201) = rxt(k,327)*y(k,52) + rxt(k,328)*y(k,216)
         mat(k,1821) = rxt(k,327)*y(k,24)
         mat(k,1666) = .500_r8*rxt(k,326)*y(k,23) + rxt(k,328)*y(k,24)
         mat(k,725) = -(rxt(k,350)*y(k,192) + rxt(k,351)*y(k,198) + rxt(k,352) &
                      *y(k,111))
         mat(k,1352) = -rxt(k,350)*y(k,189)
         mat(k,2011) = -rxt(k,351)*y(k,189)
         mat(k,1457) = -rxt(k,352)*y(k,189)
         mat(k,308) = rxt(k,353)*y(k,216)
         mat(k,57) = rxt(k,354)*y(k,216)
         mat(k,1659) = rxt(k,353)*y(k,26) + rxt(k,354)*y(k,27)
         mat(k,513) = -(rxt(k,447)*y(k,198) + rxt(k,448)*y(k,111))
         mat(k,1997) = -rxt(k,447)*y(k,190)
         mat(k,1444) = -rxt(k,448)*y(k,190)
         mat(k,183) = rxt(k,449)*y(k,216)
         mat(k,1444) = mat(k,1444) + rxt(k,438)*y(k,184)
         mat(k,1860) = rxt(k,464)*y(k,125)
         mat(k,356) = rxt(k,464)*y(k,118)
         mat(k,416) = rxt(k,438)*y(k,111) + .400_r8*rxt(k,437)*y(k,198)
         mat(k,1997) = mat(k,1997) + .400_r8*rxt(k,437)*y(k,184)
         mat(k,1642) = rxt(k,449)*y(k,28)
         mat(k,1284) = -(4._r8*rxt(k,332)*y(k,191) + rxt(k,333)*y(k,192) + rxt(k,334) &
                      *y(k,198) + rxt(k,335)*y(k,111) + rxt(k,346)*y(k,112) + rxt(k,373) &
                      *y(k,202) + rxt(k,406)*y(k,200) + rxt(k,411)*y(k,201) + rxt(k,420) &
                      *y(k,92) + rxt(k,431)*y(k,224))
         mat(k,1378) = -rxt(k,333)*y(k,191)
         mat(k,2040) = -rxt(k,334)*y(k,191)
         mat(k,1486) = -rxt(k,335)*y(k,191)
         mat(k,2081) = -rxt(k,346)*y(k,191)
         mat(k,1213) = -rxt(k,373)*y(k,191)
         mat(k,1156) = -rxt(k,406)*y(k,191)
         mat(k,1190) = -rxt(k,411)*y(k,191)
         mat(k,1109) = -rxt(k,420)*y(k,191)
         mat(k,1085) = -rxt(k,431)*y(k,191)
         mat(k,899) = .060_r8*rxt(k,481)*y(k,118)
         mat(k,1056) = rxt(k,329)*y(k,113) + rxt(k,330)*y(k,216)
         mat(k,1133) = rxt(k,355)*y(k,113) + rxt(k,356)*y(k,216)
         mat(k,427) = .500_r8*rxt(k,337)*y(k,216)
         mat(k,789) = .080_r8*rxt(k,426)*y(k,118)
         mat(k,1124) = .100_r8*rxt(k,379)*y(k,118)
         mat(k,858) = .060_r8*rxt(k,484)*y(k,118)
         mat(k,1233) = .280_r8*rxt(k,393)*y(k,118)
         mat(k,1486) = mat(k,1486) + .530_r8*rxt(k,377)*y(k,202) + rxt(k,386)*y(k,204) &
                      + rxt(k,389)*y(k,206) + rxt(k,364)*y(k,220)
         mat(k,1946) = rxt(k,329)*y(k,41) + rxt(k,355)*y(k,45) + .530_r8*rxt(k,376) &
                      *y(k,202) + rxt(k,387)*y(k,204)
         mat(k,1889) = .060_r8*rxt(k,481)*y(k,3) + .080_r8*rxt(k,426)*y(k,89) &
                      + .100_r8*rxt(k,379)*y(k,96) + .060_r8*rxt(k,484)*y(k,101) &
                      + .280_r8*rxt(k,393)*y(k,102)
         mat(k,978) = .650_r8*rxt(k,502)*y(k,216)
         mat(k,1284) = mat(k,1284) + .530_r8*rxt(k,373)*y(k,202)
         mat(k,1378) = mat(k,1378) + .260_r8*rxt(k,374)*y(k,202) + rxt(k,383)*y(k,204) &
                      + .300_r8*rxt(k,362)*y(k,220)
         mat(k,2040) = mat(k,2040) + .450_r8*rxt(k,384)*y(k,204) + .200_r8*rxt(k,388) &
                      *y(k,206) + .150_r8*rxt(k,363)*y(k,220)
         mat(k,1213) = mat(k,1213) + .530_r8*rxt(k,377)*y(k,111) + .530_r8*rxt(k,376) &
                      *y(k,113) + .530_r8*rxt(k,373)*y(k,191) + .260_r8*rxt(k,374) &
                      *y(k,192)
         mat(k,1252) = rxt(k,386)*y(k,111) + rxt(k,387)*y(k,113) + rxt(k,383)*y(k,192) &
                      + .450_r8*rxt(k,384)*y(k,198) + 4.000_r8*rxt(k,385)*y(k,204)
         mat(k,557) = rxt(k,389)*y(k,111) + .200_r8*rxt(k,388)*y(k,198)
         mat(k,1697) = rxt(k,330)*y(k,41) + rxt(k,356)*y(k,45) + .500_r8*rxt(k,337) &
                      *y(k,47) + .650_r8*rxt(k,502)*y(k,150)
         mat(k,1046) = rxt(k,364)*y(k,111) + .300_r8*rxt(k,362)*y(k,192) &
                      + .150_r8*rxt(k,363)*y(k,198)
         mat(k,1379) = -(rxt(k,225)*y(k,55) + (4._r8*rxt(k,302) + 4._r8*rxt(k,303) &
                      ) * y(k,192) + rxt(k,304)*y(k,198) + rxt(k,305)*y(k,111) &
                      + rxt(k,322)*y(k,188) + rxt(k,333)*y(k,191) + rxt(k,350) &
                      *y(k,189) + rxt(k,362)*y(k,220) + rxt(k,374)*y(k,202) + rxt(k,383) &
                      *y(k,204) + rxt(k,407)*y(k,200) + rxt(k,412)*y(k,201) + rxt(k,421) &
                      *y(k,92) + rxt(k,432)*y(k,224) + rxt(k,486)*y(k,211) + rxt(k,491) &
                      *y(k,221) + rxt(k,496)*y(k,222))
         mat(k,1728) = -rxt(k,225)*y(k,192)
         mat(k,2043) = -rxt(k,304)*y(k,192)
         mat(k,1488) = -rxt(k,305)*y(k,192)
         mat(k,808) = -rxt(k,322)*y(k,192)
         mat(k,1285) = -rxt(k,333)*y(k,192)
         mat(k,728) = -rxt(k,350)*y(k,192)
         mat(k,1047) = -rxt(k,362)*y(k,192)
         mat(k,1214) = -rxt(k,374)*y(k,192)
         mat(k,1253) = -rxt(k,383)*y(k,192)
         mat(k,1157) = -rxt(k,407)*y(k,192)
         mat(k,1191) = -rxt(k,412)*y(k,192)
         mat(k,1110) = -rxt(k,421)*y(k,192)
         mat(k,1086) = -rxt(k,432)*y(k,192)
         mat(k,967) = -rxt(k,486)*y(k,192)
         mat(k,1009) = -rxt(k,491)*y(k,192)
         mat(k,821) = -rxt(k,496)*y(k,192)
         mat(k,949) = .280_r8*rxt(k,349)*y(k,118)
         mat(k,473) = rxt(k,336)*y(k,216)
         mat(k,314) = .700_r8*rxt(k,307)*y(k,216)
         mat(k,790) = .050_r8*rxt(k,426)*y(k,118)
         mat(k,1110) = mat(k,1110) + rxt(k,420)*y(k,191)
         mat(k,1488) = mat(k,1488) + rxt(k,335)*y(k,191) + .830_r8*rxt(k,452)*y(k,193) &
                      + .170_r8*rxt(k,458)*y(k,205)
         mat(k,1892) = .280_r8*rxt(k,349)*y(k,25) + .050_r8*rxt(k,426)*y(k,89)
         mat(k,1285) = mat(k,1285) + rxt(k,420)*y(k,92) + rxt(k,335)*y(k,111) &
                      + 4.000_r8*rxt(k,332)*y(k,191) + .900_r8*rxt(k,333)*y(k,192) &
                      + .450_r8*rxt(k,334)*y(k,198) + rxt(k,406)*y(k,200) + rxt(k,411) &
                      *y(k,201) + rxt(k,373)*y(k,202) + rxt(k,382)*y(k,204) &
                      + rxt(k,431)*y(k,224)
         mat(k,1379) = mat(k,1379) + .900_r8*rxt(k,333)*y(k,191)
         mat(k,645) = .830_r8*rxt(k,452)*y(k,111) + .330_r8*rxt(k,451)*y(k,198)
         mat(k,2043) = mat(k,2043) + .450_r8*rxt(k,334)*y(k,191) + .330_r8*rxt(k,451) &
                      *y(k,193) + .070_r8*rxt(k,457)*y(k,205)
         mat(k,1157) = mat(k,1157) + rxt(k,406)*y(k,191)
         mat(k,1191) = mat(k,1191) + rxt(k,411)*y(k,191)
         mat(k,1214) = mat(k,1214) + rxt(k,373)*y(k,191)
         mat(k,1253) = mat(k,1253) + rxt(k,382)*y(k,191)
         mat(k,988) = .170_r8*rxt(k,458)*y(k,111) + .070_r8*rxt(k,457)*y(k,198)
         mat(k,1701) = rxt(k,336)*y(k,46) + .700_r8*rxt(k,307)*y(k,49)
         mat(k,1086) = mat(k,1086) + rxt(k,431)*y(k,191)
         mat(k,642) = -(rxt(k,451)*y(k,198) + rxt(k,452)*y(k,111) + rxt(k,453) &
                      *y(k,112))
         mat(k,2006) = -rxt(k,451)*y(k,193)
         mat(k,1450) = -rxt(k,452)*y(k,193)
         mat(k,2070) = -rxt(k,453)*y(k,193)
         mat(k,698) = -(rxt(k,559)*y(k,209) + rxt(k,560)*y(k,215) + rxt(k,561) &
                      *y(k,208))
         mat(k,688) = -rxt(k,559)*y(k,194)
         mat(k,680) = -rxt(k,560)*y(k,194)
         mat(k,567) = -rxt(k,561)*y(k,194)
         mat(k,476) = -((rxt(k,370) + rxt(k,371)) * y(k,111))
         mat(k,1442) = -(rxt(k,370) + rxt(k,371)) * y(k,195)
         mat(k,270) = rxt(k,369)*y(k,216)
         mat(k,1638) = rxt(k,369)*y(k,13)
      end do
      end subroutine nlnmat07
      subroutine nlnmat08( ofl, ofu, chnkpnts, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = ofl,ofu
         mat(k,343) = -(rxt(k,341)*y(k,117))
         mat(k,1395) = -rxt(k,341)*y(k,196)
         mat(k,1433) = .750_r8*rxt(k,339)*y(k,197)
         mat(k,669) = .750_r8*rxt(k,339)*y(k,111)
         mat(k,670) = -(rxt(k,338)*y(k,198) + rxt(k,339)*y(k,111))
         mat(k,2008) = -rxt(k,338)*y(k,197)
         mat(k,1451) = -rxt(k,339)*y(k,197)
         mat(k,465) = rxt(k,345)*y(k,216)
         mat(k,1655) = rxt(k,345)*y(k,21)
         mat(k,2056) = -((rxt(k,176) + rxt(k,177) + rxt(k,178)) * y(k,67) + rxt(k,180) &
                      *y(k,116) + rxt(k,181)*y(k,118) + rxt(k,185)*y(k,216) &
                      + 4._r8*rxt(k,190)*y(k,198) + rxt(k,202)*y(k,113) + rxt(k,207) &
                      *y(k,111) + rxt(k,212)*y(k,112) + (rxt(k,222) + rxt(k,223) &
                      ) * y(k,52) + rxt(k,229)*y(k,55) + rxt(k,255)*y(k,14) + rxt(k,261) &
                      *y(k,16) + rxt(k,298)*y(k,38) + rxt(k,304)*y(k,192) + rxt(k,312) &
                      *y(k,199) + rxt(k,323)*y(k,188) + rxt(k,334)*y(k,191) + rxt(k,338) &
                      *y(k,197) + rxt(k,351)*y(k,189) + rxt(k,359)*y(k,219) + rxt(k,363) &
                      *y(k,220) + rxt(k,375)*y(k,202) + rxt(k,384)*y(k,204) + rxt(k,388) &
                      *y(k,206) + rxt(k,398)*y(k,185) + rxt(k,408)*y(k,200) + rxt(k,413) &
                      *y(k,201) + rxt(k,422)*y(k,92) + rxt(k,433)*y(k,224) + rxt(k,437) &
                      *y(k,184) + rxt(k,440)*y(k,186) + rxt(k,444)*y(k,187) + rxt(k,447) &
                      *y(k,190) + rxt(k,451)*y(k,193) + rxt(k,454)*y(k,203) + rxt(k,457) &
                      *y(k,205) + rxt(k,460)*y(k,218) + rxt(k,467)*y(k,223) + rxt(k,473) &
                      *y(k,225) + rxt(k,476)*y(k,226) + rxt(k,487)*y(k,211) + rxt(k,492) &
                      *y(k,221) + rxt(k,497)*y(k,222))
         mat(k,1809) = -(rxt(k,176) + rxt(k,177) + rxt(k,178)) * y(k,198)
         mat(k,1566) = -rxt(k,180)*y(k,198)
         mat(k,1905) = -rxt(k,181)*y(k,198)
         mat(k,1714) = -rxt(k,185)*y(k,198)
         mat(k,1962) = -rxt(k,202)*y(k,198)
         mat(k,1501) = -rxt(k,207)*y(k,198)
         mat(k,2098) = -rxt(k,212)*y(k,198)
         mat(k,1844) = -(rxt(k,222) + rxt(k,223)) * y(k,198)
         mat(k,1741) = -rxt(k,229)*y(k,198)
         mat(k,1342) = -rxt(k,255)*y(k,198)
         mat(k,1525) = -rxt(k,261)*y(k,198)
         mat(k,2122) = -rxt(k,298)*y(k,198)
         mat(k,1389) = -rxt(k,304)*y(k,198)
         mat(k,330) = -rxt(k,312)*y(k,198)
         mat(k,812) = -rxt(k,323)*y(k,198)
         mat(k,1293) = -rxt(k,334)*y(k,198)
         mat(k,675) = -rxt(k,338)*y(k,198)
         mat(k,732) = -rxt(k,351)*y(k,198)
         mat(k,712) = -rxt(k,359)*y(k,198)
         mat(k,1050) = -rxt(k,363)*y(k,198)
         mat(k,1220) = -rxt(k,375)*y(k,198)
         mat(k,1260) = -rxt(k,384)*y(k,198)
         mat(k,560) = -rxt(k,388)*y(k,198)
         mat(k,930) = -rxt(k,398)*y(k,198)
         mat(k,1164) = -rxt(k,408)*y(k,198)
         mat(k,1198) = -rxt(k,413)*y(k,198)
         mat(k,1115) = -rxt(k,422)*y(k,198)
         mat(k,1091) = -rxt(k,433)*y(k,198)
         mat(k,419) = -rxt(k,437)*y(k,198)
         mat(k,387) = -rxt(k,440)*y(k,198)
         mat(k,325) = -rxt(k,444)*y(k,198)
         mat(k,517) = -rxt(k,447)*y(k,198)
         mat(k,648) = -rxt(k,451)*y(k,198)
         mat(k,593) = -rxt(k,454)*y(k,198)
         mat(k,991) = -rxt(k,457)*y(k,198)
         mat(k,338) = -rxt(k,460)*y(k,198)
         mat(k,608) = -rxt(k,467)*y(k,198)
         mat(k,640) = -rxt(k,473)*y(k,198)
         mat(k,402) = -rxt(k,476)*y(k,198)
         mat(k,971) = -rxt(k,487)*y(k,198)
         mat(k,1012) = -rxt(k,492)*y(k,198)
         mat(k,824) = -rxt(k,497)*y(k,198)
         mat(k,905) = .570_r8*rxt(k,481)*y(k,118)
         mat(k,92) = .650_r8*rxt(k,439)*y(k,216)
         mat(k,1342) = mat(k,1342) + rxt(k,254)*y(k,38)
         mat(k,1525) = mat(k,1525) + rxt(k,266)*y(k,216)
         mat(k,199) = .350_r8*rxt(k,318)*y(k,216)
         mat(k,470) = .130_r8*rxt(k,320)*y(k,118)
         mat(k,180) = rxt(k,325)*y(k,216)
         mat(k,955) = .280_r8*rxt(k,349)*y(k,118)
         mat(k,2122) = mat(k,2122) + rxt(k,254)*y(k,14) + rxt(k,218)*y(k,52) &
                      + rxt(k,299)*y(k,113) + rxt(k,300)*y(k,116)
         mat(k,55) = rxt(k,331)*y(k,216)
         mat(k,718) = rxt(k,306)*y(k,216)
         mat(k,1844) = mat(k,1844) + rxt(k,218)*y(k,38) + rxt(k,221)*y(k,70)
         mat(k,1741) = mat(k,1741) + rxt(k,225)*y(k,192) + rxt(k,236)*y(k,216)
         mat(k,1034) = rxt(k,309)*y(k,216)
         mat(k,128) = .730_r8*rxt(k,450)*y(k,216)
         mat(k,211) = .500_r8*rxt(k,520)*y(k,216)
         mat(k,997) = rxt(k,342)*y(k,216)
         mat(k,831) = rxt(k,343)*y(k,216)
         mat(k,1809) = mat(k,1809) + rxt(k,179)*y(k,117)
         mat(k,502) = rxt(k,221)*y(k,52) + rxt(k,175)*y(k,116) + rxt(k,184)*y(k,216)
         mat(k,107) = rxt(k,310)*y(k,216)
         mat(k,722) = rxt(k,311)*y(k,216)
         mat(k,1027) = rxt(k,372)*y(k,216)
         mat(k,1039) = rxt(k,357)*y(k,216)
         mat(k,794) = .370_r8*rxt(k,426)*y(k,118)
         mat(k,461) = .300_r8*rxt(k,417)*y(k,216)
         mat(k,379) = rxt(k,418)*y(k,216)
         mat(k,1115) = mat(k,1115) + rxt(k,423)*y(k,111) + rxt(k,424)*y(k,113) &
                      + rxt(k,420)*y(k,191) + 1.200_r8*rxt(k,421)*y(k,192)
         mat(k,281) = rxt(k,425)*y(k,216)
         mat(k,1128) = .140_r8*rxt(k,379)*y(k,118)
         mat(k,216) = .200_r8*rxt(k,381)*y(k,216)
         mat(k,453) = .500_r8*rxt(k,392)*y(k,216)
         mat(k,864) = .570_r8*rxt(k,484)*y(k,118)
         mat(k,1242) = .280_r8*rxt(k,393)*y(k,118)
         mat(k,226) = rxt(k,429)*y(k,216)
         mat(k,757) = rxt(k,430)*y(k,216)
         mat(k,1501) = mat(k,1501) + rxt(k,423)*y(k,92) + rxt(k,399)*y(k,185) &
                      + rxt(k,441)*y(k,186) + rxt(k,446)*y(k,187) + rxt(k,324) &
                      *y(k,188) + rxt(k,352)*y(k,189) + rxt(k,305)*y(k,192) &
                      + .170_r8*rxt(k,452)*y(k,193) + rxt(k,370)*y(k,195) &
                      + .250_r8*rxt(k,339)*y(k,197) + rxt(k,314)*y(k,199) &
                      + .920_r8*rxt(k,409)*y(k,200) + .920_r8*rxt(k,415)*y(k,201) &
                      + .470_r8*rxt(k,377)*y(k,202) + .400_r8*rxt(k,455)*y(k,203) &
                      + .830_r8*rxt(k,458)*y(k,205) + rxt(k,461)*y(k,218) + rxt(k,360) &
                      *y(k,219) + .900_r8*rxt(k,493)*y(k,221) + .800_r8*rxt(k,498) &
                      *y(k,222) + rxt(k,468)*y(k,223) + rxt(k,434)*y(k,224) &
                      + rxt(k,474)*y(k,225) + rxt(k,477)*y(k,226)
         mat(k,1962) = mat(k,1962) + rxt(k,299)*y(k,38) + rxt(k,424)*y(k,92) &
                      + rxt(k,410)*y(k,200) + rxt(k,416)*y(k,201) + .470_r8*rxt(k,376) &
                      *y(k,202) + rxt(k,205)*y(k,216) + rxt(k,435)*y(k,224)
         mat(k,1566) = mat(k,1566) + rxt(k,300)*y(k,38) + rxt(k,175)*y(k,70)
         mat(k,1420) = rxt(k,179)*y(k,67) + rxt(k,341)*y(k,196)
         mat(k,1905) = mat(k,1905) + .570_r8*rxt(k,481)*y(k,3) + .130_r8*rxt(k,320) &
                      *y(k,21) + .280_r8*rxt(k,349)*y(k,25) + .370_r8*rxt(k,426) &
                      *y(k,89) + .140_r8*rxt(k,379)*y(k,96) + .570_r8*rxt(k,484) &
                      *y(k,101) + .280_r8*rxt(k,393)*y(k,102) + rxt(k,187)*y(k,216)
         mat(k,101) = .800_r8*rxt(k,462)*y(k,216)
         mat(k,837) = rxt(k,515)*y(k,216)
         mat(k,981) = .200_r8*rxt(k,502)*y(k,216)
         mat(k,123) = .280_r8*rxt(k,470)*y(k,216)
         mat(k,149) = .380_r8*rxt(k,472)*y(k,216)
         mat(k,154) = .630_r8*rxt(k,478)*y(k,216)
         mat(k,930) = mat(k,930) + rxt(k,399)*y(k,111)
         mat(k,387) = mat(k,387) + rxt(k,441)*y(k,111)
         mat(k,325) = mat(k,325) + rxt(k,446)*y(k,111)
         mat(k,812) = mat(k,812) + rxt(k,324)*y(k,111) + 2.400_r8*rxt(k,321)*y(k,188) &
                      + rxt(k,322)*y(k,192)
         mat(k,732) = mat(k,732) + rxt(k,352)*y(k,111) + rxt(k,350)*y(k,192)
         mat(k,1293) = mat(k,1293) + rxt(k,420)*y(k,92) + .900_r8*rxt(k,333)*y(k,192) &
                      + rxt(k,406)*y(k,200) + rxt(k,411)*y(k,201) + .470_r8*rxt(k,373) &
                      *y(k,202) + rxt(k,431)*y(k,224)
         mat(k,1389) = mat(k,1389) + rxt(k,225)*y(k,55) + 1.200_r8*rxt(k,421)*y(k,92) &
                      + rxt(k,305)*y(k,111) + rxt(k,322)*y(k,188) + rxt(k,350) &
                      *y(k,189) + .900_r8*rxt(k,333)*y(k,191) + 4.000_r8*rxt(k,302) &
                      *y(k,192) + rxt(k,407)*y(k,200) + rxt(k,412)*y(k,201) &
                      + .730_r8*rxt(k,374)*y(k,202) + rxt(k,383)*y(k,204) &
                      + .500_r8*rxt(k,486)*y(k,211) + .300_r8*rxt(k,362)*y(k,220) &
                      + rxt(k,491)*y(k,221) + rxt(k,496)*y(k,222) + .800_r8*rxt(k,432) &
                      *y(k,224)
         mat(k,648) = mat(k,648) + .170_r8*rxt(k,452)*y(k,111) + .070_r8*rxt(k,451) &
                      *y(k,198)
         mat(k,481) = rxt(k,370)*y(k,111)
         mat(k,346) = rxt(k,341)*y(k,117)
         mat(k,675) = mat(k,675) + .250_r8*rxt(k,339)*y(k,111)
         mat(k,2056) = mat(k,2056) + .070_r8*rxt(k,451)*y(k,193) + .160_r8*rxt(k,454) &
                      *y(k,203) + .330_r8*rxt(k,457)*y(k,205)
         mat(k,330) = mat(k,330) + rxt(k,314)*y(k,111)
         mat(k,1164) = mat(k,1164) + .920_r8*rxt(k,409)*y(k,111) + rxt(k,410)*y(k,113) &
                      + rxt(k,406)*y(k,191) + rxt(k,407)*y(k,192)
         mat(k,1198) = mat(k,1198) + .920_r8*rxt(k,415)*y(k,111) + rxt(k,416)*y(k,113) &
                      + rxt(k,411)*y(k,191) + rxt(k,412)*y(k,192)
         mat(k,1220) = mat(k,1220) + .470_r8*rxt(k,377)*y(k,111) + .470_r8*rxt(k,376) &
                      *y(k,113) + .470_r8*rxt(k,373)*y(k,191) + .730_r8*rxt(k,374) &
                      *y(k,192)
         mat(k,593) = mat(k,593) + .400_r8*rxt(k,455)*y(k,111) + .160_r8*rxt(k,454) &
                      *y(k,198)
         mat(k,1260) = mat(k,1260) + rxt(k,383)*y(k,192)
         mat(k,991) = mat(k,991) + .830_r8*rxt(k,458)*y(k,111) + .330_r8*rxt(k,457) &
                      *y(k,198)
         mat(k,971) = mat(k,971) + .500_r8*rxt(k,486)*y(k,192)
         mat(k,1714) = mat(k,1714) + .650_r8*rxt(k,439)*y(k,4) + rxt(k,266)*y(k,16) &
                      + .350_r8*rxt(k,318)*y(k,20) + rxt(k,325)*y(k,22) + rxt(k,331) &
                      *y(k,43) + rxt(k,306)*y(k,48) + rxt(k,236)*y(k,55) + rxt(k,309) &
                      *y(k,57) + .730_r8*rxt(k,450)*y(k,61) + .500_r8*rxt(k,520) &
                      *y(k,62) + rxt(k,342)*y(k,65) + rxt(k,343)*y(k,66) + rxt(k,184) &
                      *y(k,70) + rxt(k,310)*y(k,77) + rxt(k,311)*y(k,78) + rxt(k,372) &
                      *y(k,84) + rxt(k,357)*y(k,86) + .300_r8*rxt(k,417)*y(k,90) &
                      + rxt(k,418)*y(k,91) + rxt(k,425)*y(k,93) + .200_r8*rxt(k,381) &
                      *y(k,97) + .500_r8*rxt(k,392)*y(k,100) + rxt(k,429)*y(k,106) &
                      + rxt(k,430)*y(k,107) + rxt(k,205)*y(k,113) + rxt(k,187) &
                      *y(k,118) + .800_r8*rxt(k,462)*y(k,126) + rxt(k,515)*y(k,132) &
                      + .200_r8*rxt(k,502)*y(k,150) + .280_r8*rxt(k,470)*y(k,152) &
                      + .380_r8*rxt(k,472)*y(k,154) + .630_r8*rxt(k,478)*y(k,156)
         mat(k,338) = mat(k,338) + rxt(k,461)*y(k,111)
         mat(k,712) = mat(k,712) + rxt(k,360)*y(k,111)
         mat(k,1050) = mat(k,1050) + .300_r8*rxt(k,362)*y(k,192)
         mat(k,1012) = mat(k,1012) + .900_r8*rxt(k,493)*y(k,111) + rxt(k,491)*y(k,192)
         mat(k,824) = mat(k,824) + .800_r8*rxt(k,498)*y(k,111) + rxt(k,496)*y(k,192)
         mat(k,608) = mat(k,608) + rxt(k,468)*y(k,111)
         mat(k,1091) = mat(k,1091) + rxt(k,434)*y(k,111) + rxt(k,435)*y(k,113) &
                      + rxt(k,431)*y(k,191) + .800_r8*rxt(k,432)*y(k,192)
         mat(k,640) = mat(k,640) + rxt(k,474)*y(k,111)
         mat(k,402) = mat(k,402) + rxt(k,477)*y(k,111)
         mat(k,327) = -(rxt(k,312)*y(k,198) + rxt(k,314)*y(k,111))
         mat(k,1982) = -rxt(k,312)*y(k,199)
         mat(k,1430) = -rxt(k,314)*y(k,199)
         mat(k,2102) = rxt(k,298)*y(k,198)
         mat(k,1982) = mat(k,1982) + rxt(k,298)*y(k,38)
         mat(k,1152) = -(rxt(k,406)*y(k,191) + rxt(k,407)*y(k,192) + rxt(k,408) &
                      *y(k,198) + rxt(k,409)*y(k,111) + rxt(k,410)*y(k,113))
         mat(k,1279) = -rxt(k,406)*y(k,200)
         mat(k,1373) = -rxt(k,407)*y(k,200)
         mat(k,2035) = -rxt(k,408)*y(k,200)
         mat(k,1481) = -rxt(k,409)*y(k,200)
         mat(k,1941) = -rxt(k,410)*y(k,200)
         mat(k,786) = .600_r8*rxt(k,427)*y(k,216)
         mat(k,1692) = .600_r8*rxt(k,427)*y(k,89)
      end do
      end subroutine nlnmat08
      subroutine nlnmat09( ofl, ofu, chnkpnts, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = ofl,ofu
         mat(k,1186) = -(rxt(k,411)*y(k,191) + rxt(k,412)*y(k,192) + rxt(k,413) &
                      *y(k,198) + rxt(k,415)*y(k,111) + rxt(k,416)*y(k,113))
         mat(k,1280) = -rxt(k,411)*y(k,201)
         mat(k,1374) = -rxt(k,412)*y(k,201)
         mat(k,2036) = -rxt(k,413)*y(k,201)
         mat(k,1482) = -rxt(k,415)*y(k,201)
         mat(k,1942) = -rxt(k,416)*y(k,201)
         mat(k,787) = .400_r8*rxt(k,427)*y(k,216)
         mat(k,1693) = .400_r8*rxt(k,427)*y(k,89)
         mat(k,1211) = -(rxt(k,373)*y(k,191) + rxt(k,374)*y(k,192) + rxt(k,375) &
                      *y(k,198) + rxt(k,376)*y(k,113) + (rxt(k,377) + rxt(k,378) &
                      ) * y(k,111))
         mat(k,1281) = -rxt(k,373)*y(k,202)
         mat(k,1375) = -rxt(k,374)*y(k,202)
         mat(k,2037) = -rxt(k,375)*y(k,202)
         mat(k,1943) = -rxt(k,376)*y(k,202)
         mat(k,1483) = -(rxt(k,377) + rxt(k,378)) * y(k,202)
         mat(k,1122) = .500_r8*rxt(k,380)*y(k,216)
         mat(k,213) = .200_r8*rxt(k,381)*y(k,216)
         mat(k,1230) = rxt(k,394)*y(k,216)
         mat(k,1694) = .500_r8*rxt(k,380)*y(k,96) + .200_r8*rxt(k,381)*y(k,97) &
                      + rxt(k,394)*y(k,102)
         mat(k,589) = -(rxt(k,454)*y(k,198) + rxt(k,455)*y(k,111) + rxt(k,456) &
                      *y(k,112))
         mat(k,2002) = -rxt(k,454)*y(k,203)
         mat(k,1447) = -rxt(k,455)*y(k,203)
         mat(k,2069) = -rxt(k,456)*y(k,203)
         mat(k,1251) = -(rxt(k,382)*y(k,191) + rxt(k,383)*y(k,192) + rxt(k,384) &
                      *y(k,198) + 4._r8*rxt(k,385)*y(k,204) + rxt(k,386)*y(k,111) &
                      + rxt(k,387)*y(k,113) + rxt(k,395)*y(k,112))
         mat(k,1283) = -rxt(k,382)*y(k,204)
         mat(k,1377) = -rxt(k,383)*y(k,204)
         mat(k,2039) = -rxt(k,384)*y(k,204)
         mat(k,1485) = -rxt(k,386)*y(k,204)
         mat(k,1945) = -rxt(k,387)*y(k,204)
         mat(k,2080) = -rxt(k,395)*y(k,204)
         mat(k,1123) = .500_r8*rxt(k,380)*y(k,216)
         mat(k,214) = .500_r8*rxt(k,381)*y(k,216)
         mat(k,1696) = .500_r8*rxt(k,380)*y(k,96) + .500_r8*rxt(k,381)*y(k,97)
         mat(k,985) = -(rxt(k,457)*y(k,198) + rxt(k,458)*y(k,111) + rxt(k,459) &
                      *y(k,112))
         mat(k,2023) = -rxt(k,457)*y(k,205)
         mat(k,1469) = -rxt(k,458)*y(k,205)
         mat(k,2076) = -rxt(k,459)*y(k,205)
         mat(k,555) = -(rxt(k,388)*y(k,198) + rxt(k,389)*y(k,111))
         mat(k,2000) = -rxt(k,388)*y(k,206)
         mat(k,1446) = -rxt(k,389)*y(k,206)
         mat(k,405) = rxt(k,390)*y(k,216)
         mat(k,218) = rxt(k,391)*y(k,216)
         mat(k,1646) = rxt(k,390)*y(k,98) + rxt(k,391)*y(k,99)
         mat(k,421) = -(rxt(k,192)*y(k,116) + rxt(k,193)*y(k,117))
         mat(k,1533) = -rxt(k,192)*y(k,207)
         mat(k,1397) = -rxt(k,193)*y(k,207)
         mat(k,1533) = mat(k,1533) + rxt(k,563)*y(k,208)
         mat(k,694) = .900_r8*rxt(k,561)*y(k,208) + .800_r8*rxt(k,559)*y(k,209)
         mat(k,562) = rxt(k,563)*y(k,116) + .900_r8*rxt(k,561)*y(k,194)
         mat(k,686) = .800_r8*rxt(k,559)*y(k,194)
         mat(k,564) = -(rxt(k,561)*y(k,194) + rxt(k,562)*y(k,117) + (rxt(k,563) &
                      + rxt(k,564)) * y(k,116))
         mat(k,695) = -rxt(k,561)*y(k,208)
         mat(k,1399) = -rxt(k,562)*y(k,208)
         mat(k,1537) = -(rxt(k,563) + rxt(k,564)) * y(k,208)
         mat(k,687) = -(rxt(k,559)*y(k,194))
         mat(k,697) = -rxt(k,559)*y(k,209)
         mat(k,739) = rxt(k,568)*y(k,215)
         mat(k,1453) = rxt(k,570)*y(k,215)
         mat(k,1541) = rxt(k,563)*y(k,208)
         mat(k,1402) = rxt(k,567)*y(k,210)
         mat(k,566) = rxt(k,563)*y(k,116)
         mat(k,392) = rxt(k,567)*y(k,117)
         mat(k,679) = rxt(k,568)*y(k,103) + rxt(k,570)*y(k,111)
         mat(k,389) = -(rxt(k,565)*y(k,116) + (rxt(k,566) + rxt(k,567)) * y(k,117))
         mat(k,1532) = -rxt(k,565)*y(k,210)
         mat(k,1396) = -(rxt(k,566) + rxt(k,567)) * y(k,210)
         mat(k,963) = -(rxt(k,486)*y(k,192) + rxt(k,487)*y(k,198) + rxt(k,488) &
                      *y(k,111) + rxt(k,489)*y(k,113))
         mat(k,1360) = -rxt(k,486)*y(k,211)
         mat(k,2021) = -rxt(k,487)*y(k,211)
         mat(k,1467) = -rxt(k,488)*y(k,211)
         mat(k,1927) = -rxt(k,489)*y(k,211)
         mat(k,894) = rxt(k,480)*y(k,113)
         mat(k,853) = rxt(k,483)*y(k,113)
         mat(k,1927) = mat(k,1927) + rxt(k,480)*y(k,3) + rxt(k,483)*y(k,101) &
                      + .500_r8*rxt(k,500)*y(k,149)
         mat(k,229) = rxt(k,490)*y(k,216)
         mat(k,910) = .500_r8*rxt(k,500)*y(k,113)
         mat(k,1677) = rxt(k,490)*y(k,115)
         mat(k,1761) = -(rxt(k,157)*y(k,68) + rxt(k,158)*y(k,227) + (rxt(k,160) &
                      + rxt(k,161)) * y(k,117) + rxt(k,162)*y(k,118) + (rxt(k,250) &
                      + rxt(k,251)) * y(k,76) + (rxt(k,273) + rxt(k,274)) * y(k,72) &
                      + rxt(k,279)*y(k,59) + rxt(k,280)*y(k,60) + rxt(k,315)*y(k,77))
         mat(k,1071) = -rxt(k,157)*y(k,212)
         mat(k,2142) = -rxt(k,158)*y(k,212)
         mat(k,1416) = -(rxt(k,160) + rxt(k,161)) * y(k,212)
         mat(k,1899) = -rxt(k,162)*y(k,212)
         mat(k,1325) = -(rxt(k,250) + rxt(k,251)) * y(k,212)
         mat(k,764) = -(rxt(k,273) + rxt(k,274)) * y(k,212)
         mat(k,62) = -rxt(k,279)*y(k,212)
         mat(k,110) = -rxt(k,280)*y(k,212)
         mat(k,106) = -rxt(k,315)*y(k,212)
         mat(k,1416) = mat(k,1416) + rxt(k,193)*y(k,207)
         mat(k,704) = .850_r8*rxt(k,560)*y(k,215)
         mat(k,425) = rxt(k,193)*y(k,117)
         mat(k,685) = .850_r8*rxt(k,560)*y(k,194)
         mat(k,77) = -(rxt(k,164)*y(k,116) + rxt(k,165)*y(k,117))
         mat(k,1529) = -rxt(k,164)*y(k,213)
         mat(k,1393) = -rxt(k,165)*y(k,213)
         mat(k,1529) = mat(k,1529) + rxt(k,168)*y(k,214)
         mat(k,1393) = mat(k,1393) + rxt(k,169)*y(k,214)
         mat(k,1853) = rxt(k,170)*y(k,214)
         mat(k,79) = rxt(k,168)*y(k,116) + rxt(k,169)*y(k,117) + rxt(k,170)*y(k,118)
         mat(k,80) = -(rxt(k,168)*y(k,116) + rxt(k,169)*y(k,117) + rxt(k,170)*y(k,118))
         mat(k,1530) = -rxt(k,168)*y(k,214)
         mat(k,1394) = -rxt(k,169)*y(k,214)
         mat(k,1854) = -rxt(k,170)*y(k,214)
         mat(k,1394) = mat(k,1394) + rxt(k,160)*y(k,212)
         mat(k,1746) = rxt(k,160)*y(k,117)
         mat(k,678) = -(rxt(k,560)*y(k,194) + rxt(k,568)*y(k,103) + rxt(k,570) &
                      *y(k,111))
         mat(k,696) = -rxt(k,560)*y(k,215)
         mat(k,738) = -rxt(k,568)*y(k,215)
         mat(k,1452) = -rxt(k,570)*y(k,215)
         mat(k,1401) = rxt(k,562)*y(k,208) + rxt(k,566)*y(k,210) + rxt(k,573)*y(k,217)
         mat(k,565) = rxt(k,562)*y(k,117)
         mat(k,391) = rxt(k,566)*y(k,117)
         mat(k,505) = rxt(k,573)*y(k,117)
         mat(k,1706) = -(rxt(k,183)*y(k,68) + rxt(k,184)*y(k,70) + rxt(k,185)*y(k,198) &
                      + rxt(k,186)*y(k,116) + rxt(k,187)*y(k,118) + (4._r8*rxt(k,188) &
                      + 4._r8*rxt(k,189)) * y(k,216) + rxt(k,191)*y(k,81) + rxt(k,205) &
                      *y(k,113) + rxt(k,206)*y(k,103) + rxt(k,214)*y(k,112) + rxt(k,215) &
                      *y(k,80) + rxt(k,234)*y(k,56) + (rxt(k,236) + rxt(k,237) &
                      ) * y(k,55) + rxt(k,239)*y(k,76) + rxt(k,242)*y(k,83) + rxt(k,266) &
                      *y(k,16) + rxt(k,268)*y(k,72) + rxt(k,301)*y(k,38) + rxt(k,306) &
                      *y(k,48) + rxt(k,307)*y(k,49) + (rxt(k,309) + rxt(k,316) &
                      ) * y(k,57) + rxt(k,310)*y(k,77) + rxt(k,311)*y(k,78) + rxt(k,318) &
                      *y(k,20) + rxt(k,325)*y(k,22) + rxt(k,326)*y(k,23) + rxt(k,328) &
                      *y(k,24) + rxt(k,330)*y(k,41) + rxt(k,331)*y(k,43) + rxt(k,336) &
                      *y(k,46) + rxt(k,337)*y(k,47) + rxt(k,342)*y(k,65) + rxt(k,343) &
                      *y(k,66) + rxt(k,344)*y(k,123) + rxt(k,345)*y(k,21) + rxt(k,353) &
                      *y(k,26) + rxt(k,354)*y(k,27) + rxt(k,356)*y(k,45) + rxt(k,357) &
                      *y(k,86) + rxt(k,358)*y(k,114) + rxt(k,361)*y(k,128) + rxt(k,365) &
                      *y(k,129) + rxt(k,366)*y(k,25) + rxt(k,367)*y(k,44) + rxt(k,369) &
                      *y(k,13) + rxt(k,372)*y(k,84) + rxt(k,380)*y(k,96) + rxt(k,381) &
                      *y(k,97) + rxt(k,390)*y(k,98) + rxt(k,391)*y(k,99) + rxt(k,392) &
                      *y(k,100) + rxt(k,394)*y(k,102) + rxt(k,397)*y(k,1) + rxt(k,401) &
                      *y(k,2) + rxt(k,402)*y(k,12) + rxt(k,403)*y(k,85) + rxt(k,404) &
                      *y(k,87) + rxt(k,405)*y(k,88) + rxt(k,417)*y(k,90) + rxt(k,418) &
                      *y(k,91) + rxt(k,425)*y(k,93) + rxt(k,427)*y(k,89) + rxt(k,428) &
                      *y(k,94) + rxt(k,429)*y(k,106) + rxt(k,430)*y(k,107) + rxt(k,436) &
                      *y(k,153) + rxt(k,439)*y(k,4) + rxt(k,442)*y(k,5) + rxt(k,443) &
                      *y(k,18) + rxt(k,445)*y(k,19) + rxt(k,449)*y(k,28) + rxt(k,450) &
                      *y(k,61) + rxt(k,462)*y(k,126) + rxt(k,465)*y(k,127) + rxt(k,469) &
                      *y(k,151) + rxt(k,470)*y(k,152) + rxt(k,472)*y(k,154) + rxt(k,475) &
                      *y(k,155) + rxt(k,478)*y(k,156) + rxt(k,479)*y(k,157) + rxt(k,482) &
                      *y(k,3) + rxt(k,485)*y(k,101) + rxt(k,490)*y(k,115) + rxt(k,494) &
                      *y(k,146) + rxt(k,495)*y(k,147) + rxt(k,499)*y(k,148) + rxt(k,501) &
                      *y(k,149) + rxt(k,502)*y(k,150) + rxt(k,504)*y(k,121) + rxt(k,509) &
                      *y(k,130) + rxt(k,514)*y(k,131) + rxt(k,515)*y(k,132) + (rxt(k,518) &
                      + rxt(k,520)) * y(k,62) + rxt(k,519)*y(k,108))
         mat(k,1070) = -rxt(k,183)*y(k,216)
         mat(k,500) = -rxt(k,184)*y(k,216)
         mat(k,2048) = -rxt(k,185)*y(k,216)
         mat(k,1558) = -rxt(k,186)*y(k,216)
         mat(k,1897) = -rxt(k,187)*y(k,216)
         mat(k,361) = -rxt(k,191)*y(k,216)
         mat(k,1954) = -rxt(k,205)*y(k,216)
         mat(k,746) = -rxt(k,206)*y(k,216)
         mat(k,2090) = -rxt(k,214)*y(k,216)
         mat(k,1781) = -rxt(k,215)*y(k,216)
         mat(k,873) = -rxt(k,234)*y(k,216)
         mat(k,1733) = -(rxt(k,236) + rxt(k,237)) * y(k,216)
         mat(k,1323) = -rxt(k,239)*y(k,216)
         mat(k,771) = -rxt(k,242)*y(k,216)
         mat(k,1517) = -rxt(k,266)*y(k,216)
         mat(k,763) = -rxt(k,268)*y(k,216)
         mat(k,2114) = -rxt(k,301)*y(k,216)
         mat(k,717) = -rxt(k,306)*y(k,216)
         mat(k,315) = -rxt(k,307)*y(k,216)
         mat(k,1032) = -(rxt(k,309) + rxt(k,316)) * y(k,216)
         mat(k,105) = -rxt(k,310)*y(k,216)
         mat(k,721) = -rxt(k,311)*y(k,216)
         mat(k,198) = -rxt(k,318)*y(k,216)
         mat(k,179) = -rxt(k,325)*y(k,216)
         mat(k,256) = -rxt(k,326)*y(k,216)
         mat(k,203) = -rxt(k,328)*y(k,216)
         mat(k,1058) = -rxt(k,330)*y(k,216)
         mat(k,54) = -rxt(k,331)*y(k,216)
         mat(k,474) = -rxt(k,336)*y(k,216)
         mat(k,429) = -rxt(k,337)*y(k,216)
         mat(k,996) = -rxt(k,342)*y(k,216)
         mat(k,830) = -rxt(k,343)*y(k,216)
         mat(k,351) = -rxt(k,344)*y(k,216)
         mat(k,468) = -rxt(k,345)*y(k,216)
         mat(k,310) = -rxt(k,353)*y(k,216)
         mat(k,58) = -rxt(k,354)*y(k,216)
         mat(k,1134) = -rxt(k,356)*y(k,216)
         mat(k,1038) = -rxt(k,357)*y(k,216)
         mat(k,800) = -rxt(k,358)*y(k,216)
         mat(k,436) = -rxt(k,361)*y(k,216)
         mat(k,286) = -rxt(k,365)*y(k,216)
         mat(k,952) = -rxt(k,366)*y(k,216)
         mat(k,937) = -rxt(k,367)*y(k,216)
         mat(k,273) = -rxt(k,369)*y(k,216)
         mat(k,1026) = -rxt(k,372)*y(k,216)
         mat(k,1125) = -rxt(k,380)*y(k,216)
         mat(k,215) = -rxt(k,381)*y(k,216)
         mat(k,408) = -rxt(k,390)*y(k,216)
         mat(k,221) = -rxt(k,391)*y(k,216)
         mat(k,451) = -rxt(k,392)*y(k,216)
         mat(k,1237) = -rxt(k,394)*y(k,216)
         mat(k,550) = -rxt(k,397)*y(k,216)
         mat(k,540) = -rxt(k,401)*y(k,216)
         mat(k,162) = -rxt(k,402)*y(k,216)
         mat(k,158) = -rxt(k,403)*y(k,216)
         mat(k,260) = -rxt(k,404)*y(k,216)
         mat(k,71) = -rxt(k,405)*y(k,216)
         mat(k,460) = -rxt(k,417)*y(k,216)
         mat(k,378) = -rxt(k,418)*y(k,216)
         mat(k,280) = -rxt(k,425)*y(k,216)
         mat(k,791) = -rxt(k,427)*y(k,216)
         mat(k,622) = -rxt(k,428)*y(k,216)
         mat(k,225) = -rxt(k,429)*y(k,216)
         mat(k,756) = -rxt(k,430)*y(k,216)
         mat(k,135) = -rxt(k,436)*y(k,216)
         mat(k,91) = -rxt(k,439)*y(k,216)
         mat(k,293) = -rxt(k,442)*y(k,216)
         mat(k,171) = -rxt(k,443)*y(k,216)
         mat(k,251) = -rxt(k,445)*y(k,216)
         mat(k,184) = -rxt(k,449)*y(k,216)
         mat(k,127) = -rxt(k,450)*y(k,216)
         mat(k,100) = -rxt(k,462)*y(k,216)
         mat(k,240) = -rxt(k,465)*y(k,216)
         mat(k,530) = -rxt(k,469)*y(k,216)
         mat(k,122) = -rxt(k,470)*y(k,216)
         mat(k,148) = -rxt(k,472)*y(k,216)
         mat(k,587) = -rxt(k,475)*y(k,216)
         mat(k,153) = -rxt(k,478)*y(k,216)
         mat(k,305) = -rxt(k,479)*y(k,216)
         mat(k,902) = -rxt(k,482)*y(k,216)
         mat(k,861) = -rxt(k,485)*y(k,216)
         mat(k,230) = -rxt(k,490)*y(k,216)
         mat(k,494) = -rxt(k,494)*y(k,216)
         mat(k,411) = -rxt(k,495)*y(k,216)
         mat(k,370) = -rxt(k,499)*y(k,216)
         mat(k,914) = -rxt(k,501)*y(k,216)
         mat(k,980) = -rxt(k,502)*y(k,216)
         mat(k,267) = -rxt(k,504)*y(k,216)
         mat(k,614) = -rxt(k,509)*y(k,216)
         mat(k,1306) = -rxt(k,514)*y(k,216)
         mat(k,836) = -rxt(k,515)*y(k,216)
         mat(k,208) = -(rxt(k,518) + rxt(k,520)) * y(k,216)
         mat(k,51) = -rxt(k,519)*y(k,216)
         mat(k,902) = mat(k,902) + .630_r8*rxt(k,481)*y(k,118)
         mat(k,198) = mat(k,198) + .650_r8*rxt(k,318)*y(k,216)
         mat(k,468) = mat(k,468) + .130_r8*rxt(k,320)*y(k,118)
         mat(k,256) = mat(k,256) + .500_r8*rxt(k,326)*y(k,216)
         mat(k,952) = mat(k,952) + .360_r8*rxt(k,349)*y(k,118)
         mat(k,2114) = mat(k,2114) + rxt(k,300)*y(k,116)
         mat(k,315) = mat(k,315) + .300_r8*rxt(k,307)*y(k,216)
         mat(k,1836) = rxt(k,223)*y(k,198)
         mat(k,663) = rxt(k,277)*y(k,227)
         mat(k,1801) = rxt(k,182)*y(k,118) + 2.000_r8*rxt(k,177)*y(k,198)
         mat(k,1070) = mat(k,1070) + rxt(k,174)*y(k,116) + rxt(k,157)*y(k,212)
         mat(k,500) = mat(k,500) + rxt(k,175)*y(k,116)
         mat(k,763) = mat(k,763) + rxt(k,267)*y(k,116) + rxt(k,273)*y(k,212)
         mat(k,1323) = mat(k,1323) + rxt(k,238)*y(k,116) + rxt(k,250)*y(k,212)
         mat(k,105) = mat(k,105) + rxt(k,315)*y(k,212)
         mat(k,656) = rxt(k,269)*y(k,116)
         mat(k,771) = mat(k,771) + rxt(k,241)*y(k,116)
         mat(k,791) = mat(k,791) + .320_r8*rxt(k,426)*y(k,118)
         mat(k,622) = mat(k,622) + .600_r8*rxt(k,428)*y(k,216)
         mat(k,1125) = mat(k,1125) + .240_r8*rxt(k,379)*y(k,118)
         mat(k,215) = mat(k,215) + .100_r8*rxt(k,381)*y(k,216)
         mat(k,861) = mat(k,861) + .630_r8*rxt(k,484)*y(k,118)
         mat(k,1237) = mat(k,1237) + .360_r8*rxt(k,393)*y(k,118)
         mat(k,1493) = rxt(k,207)*y(k,198)
         mat(k,1954) = mat(k,1954) + rxt(k,202)*y(k,198)
         mat(k,1558) = mat(k,1558) + rxt(k,300)*y(k,38) + rxt(k,174)*y(k,68) &
                      + rxt(k,175)*y(k,70) + rxt(k,267)*y(k,72) + rxt(k,238)*y(k,76) &
                      + rxt(k,269)*y(k,82) + rxt(k,241)*y(k,83) + rxt(k,180)*y(k,198)
         mat(k,1897) = mat(k,1897) + .630_r8*rxt(k,481)*y(k,3) + .130_r8*rxt(k,320) &
                      *y(k,21) + .360_r8*rxt(k,349)*y(k,25) + rxt(k,182)*y(k,67) &
                      + .320_r8*rxt(k,426)*y(k,89) + .240_r8*rxt(k,379)*y(k,96) &
                      + .630_r8*rxt(k,484)*y(k,101) + .360_r8*rxt(k,393)*y(k,102) &
                      + rxt(k,181)*y(k,198)
         mat(k,436) = mat(k,436) + .500_r8*rxt(k,361)*y(k,216)
         mat(k,135) = mat(k,135) + .500_r8*rxt(k,436)*y(k,216)
         mat(k,418) = .400_r8*rxt(k,437)*y(k,198)
         mat(k,1288) = .450_r8*rxt(k,334)*y(k,198)
         mat(k,647) = .400_r8*rxt(k,451)*y(k,198)
         mat(k,2048) = mat(k,2048) + rxt(k,223)*y(k,52) + 2.000_r8*rxt(k,177)*y(k,67) &
                      + rxt(k,207)*y(k,111) + rxt(k,202)*y(k,113) + rxt(k,180) &
                      *y(k,116) + rxt(k,181)*y(k,118) + .400_r8*rxt(k,437)*y(k,184) &
                      + .450_r8*rxt(k,334)*y(k,191) + .400_r8*rxt(k,451)*y(k,193) &
                      + .450_r8*rxt(k,384)*y(k,204) + .400_r8*rxt(k,457)*y(k,205) &
                      + .200_r8*rxt(k,388)*y(k,206) + .150_r8*rxt(k,363)*y(k,220)
         mat(k,1255) = .450_r8*rxt(k,384)*y(k,198)
         mat(k,990) = .400_r8*rxt(k,457)*y(k,198)
         mat(k,559) = .200_r8*rxt(k,388)*y(k,198)
         mat(k,1759) = rxt(k,157)*y(k,68) + rxt(k,273)*y(k,72) + rxt(k,250)*y(k,76) &
                      + rxt(k,315)*y(k,77) + 2.000_r8*rxt(k,158)*y(k,227)
         mat(k,1706) = mat(k,1706) + .650_r8*rxt(k,318)*y(k,20) + .500_r8*rxt(k,326) &
                      *y(k,23) + .300_r8*rxt(k,307)*y(k,49) + .600_r8*rxt(k,428) &
                      *y(k,94) + .100_r8*rxt(k,381)*y(k,97) + .500_r8*rxt(k,361) &
                      *y(k,128) + .500_r8*rxt(k,436)*y(k,153)
         mat(k,1049) = .150_r8*rxt(k,363)*y(k,198)
         mat(k,2140) = rxt(k,277)*y(k,64) + 2.000_r8*rxt(k,158)*y(k,212)
      end do
      end subroutine nlnmat09
      subroutine nlnmat10( ofl, ofu, chnkpnts, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = ofl,ofu
         mat(k,504) = -(rxt(k,573)*y(k,117))
         mat(k,1398) = -rxt(k,573)*y(k,217)
         mat(k,1536) = rxt(k,564)*y(k,208) + rxt(k,565)*y(k,210)
         mat(k,563) = rxt(k,564)*y(k,116)
         mat(k,390) = rxt(k,565)*y(k,116)
         mat(k,334) = -(rxt(k,460)*y(k,198) + rxt(k,461)*y(k,111))
         mat(k,1983) = -rxt(k,460)*y(k,218)
         mat(k,1431) = -rxt(k,461)*y(k,218)
         mat(k,125) = .200_r8*rxt(k,450)*y(k,216)
         mat(k,98) = .140_r8*rxt(k,462)*y(k,216)
         mat(k,238) = rxt(k,465)*y(k,216)
         mat(k,1620) = .200_r8*rxt(k,450)*y(k,61) + .140_r8*rxt(k,462)*y(k,126) &
                      + rxt(k,465)*y(k,127)
         mat(k,706) = -(rxt(k,359)*y(k,198) + rxt(k,360)*y(k,111))
         mat(k,2009) = -rxt(k,359)*y(k,219)
         mat(k,1455) = -rxt(k,360)*y(k,219)
         mat(k,940) = rxt(k,366)*y(k,216)
         mat(k,433) = .500_r8*rxt(k,361)*y(k,216)
         mat(k,1656) = rxt(k,366)*y(k,25) + .500_r8*rxt(k,361)*y(k,128)
         mat(k,1044) = -(rxt(k,362)*y(k,192) + rxt(k,363)*y(k,198) + rxt(k,364) &
                      *y(k,111))
         mat(k,1367) = -rxt(k,362)*y(k,220)
         mat(k,2028) = -rxt(k,363)*y(k,220)
         mat(k,1475) = -rxt(k,364)*y(k,220)
         mat(k,897) = .060_r8*rxt(k,481)*y(k,118)
         mat(k,934) = rxt(k,367)*y(k,216)
         mat(k,856) = .060_r8*rxt(k,484)*y(k,118)
         mat(k,1879) = .060_r8*rxt(k,481)*y(k,3) + .060_r8*rxt(k,484)*y(k,101)
         mat(k,284) = rxt(k,365)*y(k,216)
         mat(k,977) = .150_r8*rxt(k,502)*y(k,216)
         mat(k,1685) = rxt(k,367)*y(k,44) + rxt(k,365)*y(k,129) + .150_r8*rxt(k,502) &
                      *y(k,150)
         mat(k,1005) = -(rxt(k,491)*y(k,192) + rxt(k,492)*y(k,198) + rxt(k,493) &
                      *y(k,111))
         mat(k,1364) = -rxt(k,491)*y(k,221)
         mat(k,2025) = -rxt(k,492)*y(k,221)
         mat(k,1471) = -rxt(k,493)*y(k,221)
         mat(k,1931) = .500_r8*rxt(k,500)*y(k,149)
         mat(k,492) = rxt(k,494)*y(k,216)
         mat(k,912) = .500_r8*rxt(k,500)*y(k,113) + rxt(k,501)*y(k,216)
         mat(k,1681) = rxt(k,494)*y(k,146) + rxt(k,501)*y(k,149)
         mat(k,818) = -(rxt(k,496)*y(k,192) + rxt(k,497)*y(k,198) + rxt(k,498) &
                      *y(k,111))
         mat(k,1356) = -rxt(k,496)*y(k,222)
         mat(k,2016) = -rxt(k,497)*y(k,222)
         mat(k,1462) = -rxt(k,498)*y(k,222)
         mat(k,888) = rxt(k,482)*y(k,216)
         mat(k,847) = rxt(k,485)*y(k,216)
         mat(k,367) = rxt(k,499)*y(k,216)
         mat(k,1667) = rxt(k,482)*y(k,3) + rxt(k,485)*y(k,101) + rxt(k,499)*y(k,148)
         mat(k,600) = -(rxt(k,467)*y(k,198) + rxt(k,468)*y(k,111))
         mat(k,2003) = -rxt(k,467)*y(k,223)
         mat(k,1448) = -rxt(k,468)*y(k,223)
         mat(k,524) = rxt(k,469)*y(k,216)
         mat(k,121) = .650_r8*rxt(k,470)*y(k,216)
         mat(k,1649) = rxt(k,469)*y(k,151) + .650_r8*rxt(k,470)*y(k,152)
         mat(k,1083) = -(rxt(k,431)*y(k,191) + rxt(k,432)*y(k,192) + rxt(k,433) &
                      *y(k,198) + rxt(k,434)*y(k,111) + rxt(k,435)*y(k,113))
         mat(k,1275) = -rxt(k,431)*y(k,224)
         mat(k,1369) = -rxt(k,432)*y(k,224)
         mat(k,2031) = -rxt(k,433)*y(k,224)
         mat(k,1477) = -rxt(k,434)*y(k,224)
         mat(k,1937) = -rxt(k,435)*y(k,224)
         mat(k,157) = rxt(k,403)*y(k,216)
         mat(k,259) = rxt(k,404)*y(k,216)
         mat(k,70) = rxt(k,405)*y(k,216)
         mat(k,619) = .400_r8*rxt(k,428)*y(k,216)
         mat(k,134) = .500_r8*rxt(k,436)*y(k,216)
         mat(k,1688) = rxt(k,403)*y(k,85) + rxt(k,404)*y(k,87) + rxt(k,405)*y(k,88) &
                      + .400_r8*rxt(k,428)*y(k,94) + .500_r8*rxt(k,436)*y(k,153)
         mat(k,631) = -(rxt(k,473)*y(k,198) + rxt(k,474)*y(k,111))
         mat(k,2005) = -rxt(k,473)*y(k,225)
         mat(k,1449) = -rxt(k,474)*y(k,225)
         mat(k,145) = .560_r8*rxt(k,472)*y(k,216)
         mat(k,580) = rxt(k,475)*y(k,216)
         mat(k,1652) = .560_r8*rxt(k,472)*y(k,154) + rxt(k,475)*y(k,155)
         mat(k,397) = -(rxt(k,476)*y(k,198) + rxt(k,477)*y(k,111))
         mat(k,1990) = -rxt(k,476)*y(k,226)
         mat(k,1437) = -rxt(k,477)*y(k,226)
         mat(k,152) = .300_r8*rxt(k,478)*y(k,216)
         mat(k,302) = rxt(k,479)*y(k,216)
         mat(k,1628) = .300_r8*rxt(k,478)*y(k,156) + rxt(k,479)*y(k,157)
         mat(k,2151) = -(rxt(k,158)*y(k,212) + rxt(k,277)*y(k,64) + rxt(k,516) &
                      *y(k,133))
         mat(k,1770) = -rxt(k,158)*y(k,227)
         mat(k,667) = -rxt(k,277)*y(k,227)
         mat(k,176) = -rxt(k,516)*y(k,227)
         mat(k,205) = rxt(k,328)*y(k,216)
         mat(k,312) = rxt(k,353)*y(k,216)
         mat(k,59) = rxt(k,354)*y(k,216)
         mat(k,2125) = rxt(k,301)*y(k,216)
         mat(k,1063) = rxt(k,330)*y(k,216)
         mat(k,938) = rxt(k,367)*y(k,216)
         mat(k,1139) = rxt(k,356)*y(k,216)
         mat(k,475) = rxt(k,336)*y(k,216)
         mat(k,431) = rxt(k,337)*y(k,216)
         mat(k,318) = rxt(k,307)*y(k,216)
         mat(k,1812) = rxt(k,178)*y(k,198)
         mat(k,1076) = rxt(k,183)*y(k,216)
         mat(k,503) = rxt(k,184)*y(k,216)
         mat(k,766) = rxt(k,268)*y(k,216)
         mat(k,1331) = (rxt(k,550)+rxt(k,555))*y(k,82) + (rxt(k,543)+rxt(k,549) &
                       +rxt(k,554))*y(k,83) + rxt(k,239)*y(k,216)
         mat(k,723) = rxt(k,311)*y(k,216)
         mat(k,1792) = rxt(k,215)*y(k,216)
         mat(k,365) = rxt(k,191)*y(k,216)
         mat(k,658) = (rxt(k,550)+rxt(k,555))*y(k,76)
         mat(k,774) = (rxt(k,543)+rxt(k,549)+rxt(k,554))*y(k,76) + rxt(k,242)*y(k,216)
         mat(k,1130) = .500_r8*rxt(k,380)*y(k,216)
         mat(k,52) = rxt(k,519)*y(k,216)
         mat(k,439) = rxt(k,361)*y(k,216)
         mat(k,288) = rxt(k,365)*y(k,216)
         mat(k,2059) = rxt(k,178)*y(k,67) + rxt(k,185)*y(k,216)
         mat(k,1717) = rxt(k,328)*y(k,24) + rxt(k,353)*y(k,26) + rxt(k,354)*y(k,27) &
                      + rxt(k,301)*y(k,38) + rxt(k,330)*y(k,41) + rxt(k,367)*y(k,44) &
                      + rxt(k,356)*y(k,45) + rxt(k,336)*y(k,46) + rxt(k,337)*y(k,47) &
                      + rxt(k,307)*y(k,49) + rxt(k,183)*y(k,68) + rxt(k,184)*y(k,70) &
                      + rxt(k,268)*y(k,72) + rxt(k,239)*y(k,76) + rxt(k,311)*y(k,78) &
                      + rxt(k,215)*y(k,80) + rxt(k,191)*y(k,81) + rxt(k,242)*y(k,83) &
                      + .500_r8*rxt(k,380)*y(k,96) + rxt(k,519)*y(k,108) + rxt(k,361) &
                      *y(k,128) + rxt(k,365)*y(k,129) + rxt(k,185)*y(k,198) &
                      + 2.000_r8*rxt(k,188)*y(k,216)
      end do
      end subroutine nlnmat10
      subroutine nlnmat_finit( ofl, ofu, chnkpnts, mat, lmat, dti )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: dti
      real(r8), intent(in) :: lmat(chnkpnts,nzcnt)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = ofl,ofu
         mat(k, 1) = lmat(k, 1)
         mat(k, 2) = lmat(k, 2)
         mat(k, 3) = lmat(k, 3)
         mat(k, 4) = lmat(k, 4)
         mat(k, 5) = lmat(k, 5)
         mat(k, 6) = lmat(k, 6)
         mat(k, 7) = lmat(k, 7)
         mat(k, 8) = lmat(k, 8)
         mat(k, 9) = lmat(k, 9)
         mat(k, 10) = lmat(k, 10)
         mat(k, 11) = lmat(k, 11)
         mat(k, 12) = lmat(k, 12)
         mat(k, 13) = lmat(k, 13)
         mat(k, 14) = lmat(k, 14)
         mat(k, 15) = lmat(k, 15)
         mat(k, 16) = lmat(k, 16)
         mat(k, 17) = lmat(k, 17)
         mat(k, 18) = lmat(k, 18)
         mat(k, 19) = lmat(k, 19)
         mat(k, 20) = lmat(k, 20)
         mat(k, 21) = lmat(k, 21)
         mat(k, 22) = lmat(k, 22)
         mat(k, 23) = lmat(k, 23)
         mat(k, 24) = lmat(k, 24)
         mat(k, 25) = lmat(k, 25)
         mat(k, 26) = lmat(k, 26)
         mat(k, 27) = lmat(k, 27)
         mat(k, 28) = lmat(k, 28)
         mat(k, 29) = lmat(k, 29)
         mat(k, 30) = lmat(k, 30)
         mat(k, 31) = lmat(k, 31)
         mat(k, 32) = lmat(k, 32)
         mat(k, 33) = lmat(k, 33)
         mat(k, 34) = lmat(k, 34)
         mat(k, 40) = mat(k, 40) + lmat(k, 40)
         mat(k, 46) = mat(k, 46) + lmat(k, 46)
         mat(k, 47) = lmat(k, 47)
         mat(k, 48) = lmat(k, 48)
         mat(k, 49) = lmat(k, 49)
         mat(k, 50) = mat(k, 50) + lmat(k, 50)
         mat(k, 53) = mat(k, 53) + lmat(k, 53)
         mat(k, 56) = mat(k, 56) + lmat(k, 56)
         mat(k, 60) = mat(k, 60) + lmat(k, 60)
         mat(k, 61) = mat(k, 61) + lmat(k, 61)
         mat(k, 63) = lmat(k, 63)
         mat(k, 64) = lmat(k, 64)
         mat(k, 65) = lmat(k, 65)
         mat(k, 66) = lmat(k, 66)
         mat(k, 67) = lmat(k, 67)
         mat(k, 68) = lmat(k, 68)
         mat(k, 69) = mat(k, 69) + lmat(k, 69)
         mat(k, 72) = lmat(k, 72)
         mat(k, 73) = lmat(k, 73)
         mat(k, 74) = lmat(k, 74)
         mat(k, 75) = lmat(k, 75)
         mat(k, 76) = lmat(k, 76)
         mat(k, 77) = mat(k, 77) + lmat(k, 77)
         mat(k, 78) = mat(k, 78) + lmat(k, 78)
         mat(k, 79) = mat(k, 79) + lmat(k, 79)
         mat(k, 80) = mat(k, 80) + lmat(k, 80)
         mat(k, 81) = lmat(k, 81)
         mat(k, 87) = mat(k, 87) + lmat(k, 87)
         mat(k, 93) = lmat(k, 93)
         mat(k, 94) = lmat(k, 94)
         mat(k, 95) = lmat(k, 95)
         mat(k, 96) = lmat(k, 96)
         mat(k, 97) = mat(k, 97) + lmat(k, 97)
         mat(k, 102) = lmat(k, 102)
         mat(k, 103) = lmat(k, 103)
         mat(k, 104) = mat(k, 104) + lmat(k, 104)
         mat(k, 108) = mat(k, 108) + lmat(k, 108)
         mat(k, 109) = mat(k, 109) + lmat(k, 109)
         mat(k, 111) = mat(k, 111) + lmat(k, 111)
         mat(k, 117) = mat(k, 117) + lmat(k, 117)
         mat(k, 124) = mat(k, 124) + lmat(k, 124)
         mat(k, 129) = lmat(k, 129)
         mat(k, 130) = lmat(k, 130)
         mat(k, 131) = lmat(k, 131)
         mat(k, 132) = lmat(k, 132)
         mat(k, 133) = mat(k, 133) + lmat(k, 133)
         mat(k, 135) = mat(k, 135) + lmat(k, 135)
         mat(k, 142) = mat(k, 142) + lmat(k, 142)
         mat(k, 150) = mat(k, 150) + lmat(k, 150)
         mat(k, 155) = mat(k, 155) + lmat(k, 155)
         mat(k, 156) = lmat(k, 156)
         mat(k, 158) = mat(k, 158) + lmat(k, 158)
         mat(k, 159) = lmat(k, 159)
         mat(k, 160) = mat(k, 160) + lmat(k, 160)
         mat(k, 163) = lmat(k, 163)
         mat(k, 164) = lmat(k, 164)
         mat(k, 165) = lmat(k, 165)
         mat(k, 166) = lmat(k, 166)
         mat(k, 167) = lmat(k, 167)
         mat(k, 168) = lmat(k, 168)
         mat(k, 169) = mat(k, 169) + lmat(k, 169)
         mat(k, 173) = mat(k, 173) + lmat(k, 173)
         mat(k, 174) = lmat(k, 174)
         mat(k, 175) = lmat(k, 175)
         mat(k, 177) = mat(k, 177) + lmat(k, 177)
         mat(k, 181) = mat(k, 181) + lmat(k, 181)
         mat(k, 182) = lmat(k, 182)
         mat(k, 184) = mat(k, 184) + lmat(k, 184)
         mat(k, 185) = lmat(k, 185)
         mat(k, 186) = lmat(k, 186)
         mat(k, 187) = lmat(k, 187)
         mat(k, 188) = lmat(k, 188)
         mat(k, 189) = lmat(k, 189)
         mat(k, 190) = lmat(k, 190)
         mat(k, 191) = lmat(k, 191)
         mat(k, 192) = lmat(k, 192)
         mat(k, 193) = lmat(k, 193)
         mat(k, 194) = mat(k, 194) + lmat(k, 194)
         mat(k, 200) = mat(k, 200) + lmat(k, 200)
         mat(k, 206) = mat(k, 206) + lmat(k, 206)
         mat(k, 212) = mat(k, 212) + lmat(k, 212)
         mat(k, 217) = mat(k, 217) + lmat(k, 217)
         mat(k, 219) = lmat(k, 219)
         mat(k, 220) = lmat(k, 220)
         mat(k, 221) = mat(k, 221) + lmat(k, 221)
         mat(k, 222) = mat(k, 222) + lmat(k, 222)
         mat(k, 227) = mat(k, 227) + lmat(k, 227)
         mat(k, 228) = lmat(k, 228)
         mat(k, 230) = mat(k, 230) + lmat(k, 230)
         mat(k, 231) = lmat(k, 231)
         mat(k, 232) = mat(k, 232) + lmat(k, 232)
         mat(k, 235) = lmat(k, 235)
         mat(k, 236) = mat(k, 236) + lmat(k, 236)
         mat(k, 237) = mat(k, 237) + lmat(k, 237)
         mat(k, 239) = lmat(k, 239)
         mat(k, 240) = mat(k, 240) + lmat(k, 240)
         mat(k, 241) = lmat(k, 241)
         mat(k, 242) = lmat(k, 242)
         mat(k, 243) = lmat(k, 243)
         mat(k, 244) = lmat(k, 244)
         mat(k, 245) = lmat(k, 245)
         mat(k, 246) = lmat(k, 246)
         mat(k, 247) = lmat(k, 247)
         mat(k, 248) = mat(k, 248) + lmat(k, 248)
         mat(k, 251) = mat(k, 251) + lmat(k, 251)
         mat(k, 252) = lmat(k, 252)
         mat(k, 253) = mat(k, 253) + lmat(k, 253)
         mat(k, 255) = mat(k, 255) + lmat(k, 255)
         mat(k, 256) = mat(k, 256) + lmat(k, 256)
         mat(k, 257) = lmat(k, 257)
         mat(k, 258) = mat(k, 258) + lmat(k, 258)
         mat(k, 261) = mat(k, 261) + lmat(k, 261)
         mat(k, 262) = lmat(k, 262)
         mat(k, 264) = mat(k, 264) + lmat(k, 264)
         mat(k, 269) = mat(k, 269) + lmat(k, 269)
         mat(k, 277) = mat(k, 277) + lmat(k, 277)
         mat(k, 278) = lmat(k, 278)
         mat(k, 281) = mat(k, 281) + lmat(k, 281)
         mat(k, 282) = lmat(k, 282)
         mat(k, 283) = mat(k, 283) + lmat(k, 283)
         mat(k, 285) = lmat(k, 285)
         mat(k, 286) = mat(k, 286) + lmat(k, 286)
         mat(k, 287) = lmat(k, 287)
         mat(k, 289) = mat(k, 289) + lmat(k, 289)
         mat(k, 290) = lmat(k, 290)
         mat(k, 292) = lmat(k, 292)
         mat(k, 293) = mat(k, 293) + lmat(k, 293)
         mat(k, 294) = lmat(k, 294)
         mat(k, 295) = lmat(k, 295)
         mat(k, 296) = lmat(k, 296)
         mat(k, 297) = lmat(k, 297)
         mat(k, 298) = lmat(k, 298)
         mat(k, 299) = lmat(k, 299)
         mat(k, 300) = lmat(k, 300)
         mat(k, 301) = mat(k, 301) + lmat(k, 301)
         mat(k, 303) = lmat(k, 303)
         mat(k, 304) = lmat(k, 304)
         mat(k, 305) = mat(k, 305) + lmat(k, 305)
         mat(k, 306) = lmat(k, 306)
         mat(k, 307) = mat(k, 307) + lmat(k, 307)
         mat(k, 309) = lmat(k, 309)
         mat(k, 310) = mat(k, 310) + lmat(k, 310)
         mat(k, 311) = lmat(k, 311)
         mat(k, 313) = mat(k, 313) + lmat(k, 313)
         mat(k, 315) = mat(k, 315) + lmat(k, 315)
         mat(k, 316) = lmat(k, 316)
         mat(k, 317) = mat(k, 317) + lmat(k, 317)
         mat(k, 321) = mat(k, 321) + lmat(k, 321)
         mat(k, 327) = mat(k, 327) + lmat(k, 327)
         mat(k, 330) = mat(k, 330) + lmat(k, 330)
         mat(k, 332) = lmat(k, 332)
         mat(k, 334) = mat(k, 334) + lmat(k, 334)
         mat(k, 340) = lmat(k, 340)
         mat(k, 341) = lmat(k, 341)
         mat(k, 342) = lmat(k, 342)
         mat(k, 343) = mat(k, 343) + lmat(k, 343)
         mat(k, 346) = mat(k, 346) + lmat(k, 346)
         mat(k, 347) = lmat(k, 347)
         mat(k, 348) = mat(k, 348) + lmat(k, 348)
         mat(k, 349) = lmat(k, 349)
         mat(k, 350) = lmat(k, 350)
         mat(k, 352) = mat(k, 352) + lmat(k, 352)
         mat(k, 353) = lmat(k, 353)
         mat(k, 355) = mat(k, 355) + lmat(k, 355)
         mat(k, 359) = mat(k, 359) + lmat(k, 359)
         mat(k, 361) = mat(k, 361) + lmat(k, 361)
         mat(k, 362) = lmat(k, 362)
         mat(k, 363) = lmat(k, 363)
         mat(k, 364) = mat(k, 364) + lmat(k, 364)
         mat(k, 366) = mat(k, 366) + lmat(k, 366)
         mat(k, 368) = lmat(k, 368)
         mat(k, 369) = lmat(k, 369)
         mat(k, 370) = mat(k, 370) + lmat(k, 370)
         mat(k, 371) = lmat(k, 371)
         mat(k, 372) = lmat(k, 372)
         mat(k, 373) = mat(k, 373) + lmat(k, 373)
         mat(k, 382) = mat(k, 382) + lmat(k, 382)
         mat(k, 389) = mat(k, 389) + lmat(k, 389)
         mat(k, 397) = mat(k, 397) + lmat(k, 397)
         mat(k, 404) = mat(k, 404) + lmat(k, 404)
         mat(k, 406) = lmat(k, 406)
         mat(k, 407) = lmat(k, 407)
         mat(k, 409) = mat(k, 409) + lmat(k, 409)
         mat(k, 410) = mat(k, 410) + lmat(k, 410)
         mat(k, 412) = lmat(k, 412)
         mat(k, 413) = mat(k, 413) + lmat(k, 413)
         mat(k, 415) = mat(k, 415) + lmat(k, 415)
         mat(k, 421) = mat(k, 421) + lmat(k, 421)
         mat(k, 426) = mat(k, 426) + lmat(k, 426)
         mat(k, 428) = lmat(k, 428)
         mat(k, 429) = mat(k, 429) + lmat(k, 429)
         mat(k, 432) = mat(k, 432) + lmat(k, 432)
         mat(k, 435) = lmat(k, 435)
         mat(k, 436) = mat(k, 436) + lmat(k, 436)
         mat(k, 437) = lmat(k, 437)
         mat(k, 438) = lmat(k, 438)
         mat(k, 440) = mat(k, 440) + lmat(k, 440)
         mat(k, 441) = lmat(k, 441)
         mat(k, 442) = lmat(k, 442)
         mat(k, 443) = mat(k, 443) + lmat(k, 443)
         mat(k, 445) = lmat(k, 445)
         mat(k, 446) = mat(k, 446) + lmat(k, 446)
         mat(k, 447) = lmat(k, 447)
         mat(k, 448) = mat(k, 448) + lmat(k, 448)
         mat(k, 450) = lmat(k, 450)
         mat(k, 454) = lmat(k, 454)
         mat(k, 456) = mat(k, 456) + lmat(k, 456)
         mat(k, 464) = mat(k, 464) + lmat(k, 464)
         mat(k, 472) = mat(k, 472) + lmat(k, 472)
         mat(k, 476) = mat(k, 476) + lmat(k, 476)
         mat(k, 484) = lmat(k, 484)
         mat(k, 485) = lmat(k, 485)
         mat(k, 486) = lmat(k, 486)
         mat(k, 487) = lmat(k, 487)
         mat(k, 488) = mat(k, 488) + lmat(k, 488)
         mat(k, 489) = lmat(k, 489)
         mat(k, 490) = lmat(k, 490)
         mat(k, 491) = lmat(k, 491)
         mat(k, 493) = lmat(k, 493)
         mat(k, 494) = mat(k, 494) + lmat(k, 494)
         mat(k, 495) = lmat(k, 495)
         mat(k, 496) = lmat(k, 496)
         mat(k, 497) = mat(k, 497) + lmat(k, 497)
         mat(k, 500) = mat(k, 500) + lmat(k, 500)
         mat(k, 504) = mat(k, 504) + lmat(k, 504)
         mat(k, 505) = mat(k, 505) + lmat(k, 505)
         mat(k, 506) = lmat(k, 506)
         mat(k, 507) = lmat(k, 507)
         mat(k, 508) = lmat(k, 508)
         mat(k, 513) = mat(k, 513) + lmat(k, 513)
         mat(k, 519) = lmat(k, 519)
         mat(k, 520) = lmat(k, 520)
         mat(k, 521) = lmat(k, 521)
         mat(k, 522) = mat(k, 522) + lmat(k, 522)
         mat(k, 526) = lmat(k, 526)
         mat(k, 529) = lmat(k, 529)
         mat(k, 530) = mat(k, 530) + lmat(k, 530)
         mat(k, 531) = lmat(k, 531)
         mat(k, 532) = lmat(k, 532)
         mat(k, 533) = mat(k, 533) + lmat(k, 533)
         mat(k, 537) = lmat(k, 537)
         mat(k, 538) = lmat(k, 538)
         mat(k, 540) = mat(k, 540) + lmat(k, 540)
         mat(k, 541) = lmat(k, 541)
         mat(k, 542) = lmat(k, 542)
         mat(k, 543) = lmat(k, 543)
         mat(k, 544) = mat(k, 544) + lmat(k, 544)
         mat(k, 547) = mat(k, 547) + lmat(k, 547)
         mat(k, 548) = mat(k, 548) + lmat(k, 548)
         mat(k, 551) = lmat(k, 551)
         mat(k, 552) = mat(k, 552) + lmat(k, 552)
         mat(k, 553) = mat(k, 553) + lmat(k, 553)
         mat(k, 555) = mat(k, 555) + lmat(k, 555)
         mat(k, 564) = mat(k, 564) + lmat(k, 564)
         mat(k, 574) = lmat(k, 574)
         mat(k, 575) = lmat(k, 575)
         mat(k, 576) = lmat(k, 576)
         mat(k, 577) = lmat(k, 577)
         mat(k, 578) = mat(k, 578) + lmat(k, 578)
         mat(k, 582) = lmat(k, 582)
         mat(k, 585) = lmat(k, 585)
         mat(k, 587) = mat(k, 587) + lmat(k, 587)
         mat(k, 588) = lmat(k, 588)
         mat(k, 589) = mat(k, 589) + lmat(k, 589)
         mat(k, 600) = mat(k, 600) + lmat(k, 600)
         mat(k, 610) = mat(k, 610) + lmat(k, 610)
         mat(k, 618) = mat(k, 618) + lmat(k, 618)
         mat(k, 620) = lmat(k, 620)
         mat(k, 621) = lmat(k, 621)
         mat(k, 623) = lmat(k, 623)
         mat(k, 624) = lmat(k, 624)
         mat(k, 631) = mat(k, 631) + lmat(k, 631)
         mat(k, 642) = mat(k, 642) + lmat(k, 642)
         mat(k, 651) = mat(k, 651) + lmat(k, 651)
         mat(k, 653) = lmat(k, 653)
         mat(k, 656) = mat(k, 656) + lmat(k, 656)
         mat(k, 659) = mat(k, 659) + lmat(k, 659)
         mat(k, 660) = mat(k, 660) + lmat(k, 660)
         mat(k, 662) = lmat(k, 662)
         mat(k, 670) = mat(k, 670) + lmat(k, 670)
         mat(k, 678) = mat(k, 678) + lmat(k, 678)
         mat(k, 679) = mat(k, 679) + lmat(k, 679)
         mat(k, 683) = mat(k, 683) + lmat(k, 683)
         mat(k, 687) = mat(k, 687) + lmat(k, 687)
         mat(k, 698) = mat(k, 698) + lmat(k, 698)
         mat(k, 706) = mat(k, 706) + lmat(k, 706)
         mat(k, 716) = mat(k, 716) + lmat(k, 716)
         mat(k, 720) = mat(k, 720) + lmat(k, 720)
         mat(k, 725) = mat(k, 725) + lmat(k, 725)
         mat(k, 736) = lmat(k, 736)
         mat(k, 740) = lmat(k, 740)
         mat(k, 741) = mat(k, 741) + lmat(k, 741)
         mat(k, 750) = lmat(k, 750)
         mat(k, 751) = mat(k, 751) + lmat(k, 751)
         mat(k, 757) = mat(k, 757) + lmat(k, 757)
         mat(k, 758) = lmat(k, 758)
         mat(k, 759) = mat(k, 759) + lmat(k, 759)
         mat(k, 760) = mat(k, 760) + lmat(k, 760)
         mat(k, 765) = mat(k, 765) + lmat(k, 765)
         mat(k, 768) = mat(k, 768) + lmat(k, 768)
         mat(k, 771) = mat(k, 771) + lmat(k, 771)
         mat(k, 773) = mat(k, 773) + lmat(k, 773)
         mat(k, 781) = mat(k, 781) + lmat(k, 781)
         mat(k, 797) = mat(k, 797) + lmat(k, 797)
         mat(k, 799) = lmat(k, 799)
         mat(k, 801) = mat(k, 801) + lmat(k, 801)
         mat(k, 802) = lmat(k, 802)
         mat(k, 806) = mat(k, 806) + lmat(k, 806)
         mat(k, 818) = mat(k, 818) + lmat(k, 818)
         mat(k, 827) = lmat(k, 827)
         mat(k, 828) = mat(k, 828) + lmat(k, 828)
         mat(k, 829) = mat(k, 829) + lmat(k, 829)
         mat(k, 831) = mat(k, 831) + lmat(k, 831)
         mat(k, 833) = mat(k, 833) + lmat(k, 833)
         mat(k, 834) = lmat(k, 834)
         mat(k, 835) = lmat(k, 835)
         mat(k, 849) = mat(k, 849) + lmat(k, 849)
         mat(k, 869) = mat(k, 869) + lmat(k, 869)
         mat(k, 870) = mat(k, 870) + lmat(k, 870)
         mat(k, 874) = mat(k, 874) + lmat(k, 874)
         mat(k, 875) = mat(k, 875) + lmat(k, 875)
         mat(k, 876) = mat(k, 876) + lmat(k, 876)
         mat(k, 877) = mat(k, 877) + lmat(k, 877)
         mat(k, 878) = lmat(k, 878)
         mat(k, 890) = mat(k, 890) + lmat(k, 890)
         mat(k, 909) = mat(k, 909) + lmat(k, 909)
         mat(k, 911) = lmat(k, 911)
         mat(k, 913) = lmat(k, 913)
         mat(k, 916) = lmat(k, 916)
         mat(k, 922) = mat(k, 922) + lmat(k, 922)
         mat(k, 933) = mat(k, 933) + lmat(k, 933)
         mat(k, 935) = lmat(k, 935)
         mat(k, 936) = lmat(k, 936)
         mat(k, 943) = mat(k, 943) + lmat(k, 943)
         mat(k, 963) = mat(k, 963) + lmat(k, 963)
         mat(k, 974) = mat(k, 974) + lmat(k, 974)
         mat(k, 975) = mat(k, 975) + lmat(k, 975)
         mat(k, 976) = mat(k, 976) + lmat(k, 976)
         mat(k, 977) = mat(k, 977) + lmat(k, 977)
         mat(k, 978) = mat(k, 978) + lmat(k, 978)
         mat(k, 981) = mat(k, 981) + lmat(k, 981)
         mat(k, 982) = mat(k, 982) + lmat(k, 982)
         mat(k, 985) = mat(k, 985) + lmat(k, 985)
         mat(k, 994) = mat(k, 994) + lmat(k, 994)
         mat(k, 995) = lmat(k, 995)
         mat(k, 997) = mat(k, 997) + lmat(k, 997)
         mat(k, 998) = mat(k, 998) + lmat(k, 998)
         mat(k,1005) = mat(k,1005) + lmat(k,1005)
         mat(k,1017) = lmat(k,1017)
         mat(k,1018) = lmat(k,1018)
         mat(k,1019) = mat(k,1019) + lmat(k,1019)
         mat(k,1020) = lmat(k,1020)
         mat(k,1021) = lmat(k,1021)
         mat(k,1023) = lmat(k,1023)
         mat(k,1024) = lmat(k,1024)
         mat(k,1027) = mat(k,1027) + lmat(k,1027)
         mat(k,1028) = lmat(k,1028)
         mat(k,1029) = lmat(k,1029)
         mat(k,1031) = mat(k,1031) + lmat(k,1031)
         mat(k,1035) = mat(k,1035) + lmat(k,1035)
         mat(k,1037) = lmat(k,1037)
         mat(k,1039) = mat(k,1039) + lmat(k,1039)
         mat(k,1040) = lmat(k,1040)
         mat(k,1044) = mat(k,1044) + lmat(k,1044)
         mat(k,1054) = lmat(k,1054)
         mat(k,1055) = mat(k,1055) + lmat(k,1055)
         mat(k,1057) = lmat(k,1057)
         mat(k,1062) = lmat(k,1062)
         mat(k,1066) = mat(k,1066) + lmat(k,1066)
         mat(k,1083) = mat(k,1083) + lmat(k,1083)
         mat(k,1105) = mat(k,1105) + lmat(k,1105)
         mat(k,1119) = mat(k,1119) + lmat(k,1119)
         mat(k,1120) = mat(k,1120) + lmat(k,1120)
         mat(k,1123) = mat(k,1123) + lmat(k,1123)
         mat(k,1124) = mat(k,1124) + lmat(k,1124)
         mat(k,1128) = mat(k,1128) + lmat(k,1128)
         mat(k,1129) = mat(k,1129) + lmat(k,1129)
         mat(k,1131) = mat(k,1131) + lmat(k,1131)
         mat(k,1132) = mat(k,1132) + lmat(k,1132)
         mat(k,1133) = mat(k,1133) + lmat(k,1133)
         mat(k,1138) = lmat(k,1138)
         mat(k,1152) = mat(k,1152) + lmat(k,1152)
         mat(k,1168) = lmat(k,1168)
         mat(k,1186) = mat(k,1186) + lmat(k,1186)
         mat(k,1198) = mat(k,1198) + lmat(k,1198)
         mat(k,1211) = mat(k,1211) + lmat(k,1211)
         mat(k,1225) = lmat(k,1225)
         mat(k,1226) = mat(k,1226) + lmat(k,1226)
         mat(k,1231) = mat(k,1231) + lmat(k,1231)
         mat(k,1233) = mat(k,1233) + lmat(k,1233)
         mat(k,1234) = lmat(k,1234)
         mat(k,1251) = mat(k,1251) + lmat(k,1251)
         mat(k,1284) = mat(k,1284) + lmat(k,1284)
         mat(k,1298) = lmat(k,1298)
         mat(k,1300) = mat(k,1300) + lmat(k,1300)
         mat(k,1305) = mat(k,1305) + lmat(k,1305)
         mat(k,1319) = mat(k,1319) + lmat(k,1319)
         mat(k,1327) = mat(k,1327) + lmat(k,1327)
         mat(k,1328) = mat(k,1328) + lmat(k,1328)
         mat(k,1334) = mat(k,1334) + lmat(k,1334)
         mat(k,1379) = mat(k,1379) + lmat(k,1379)
         mat(k,1398) = mat(k,1398) + lmat(k,1398)
         mat(k,1401) = mat(k,1401) + lmat(k,1401)
         mat(k,1403) = lmat(k,1403)
         mat(k,1410) = mat(k,1410) + lmat(k,1410)
         mat(k,1413) = mat(k,1413) + lmat(k,1413)
         mat(k,1416) = mat(k,1416) + lmat(k,1416)
         mat(k,1453) = mat(k,1453) + lmat(k,1453)
         mat(k,1454) = lmat(k,1454)
         mat(k,1458) = mat(k,1458) + lmat(k,1458)
         mat(k,1490) = mat(k,1490) + lmat(k,1490)
         mat(k,1492) = mat(k,1492) + lmat(k,1492)
         mat(k,1512) = mat(k,1512) + lmat(k,1512)
         mat(k,1515) = mat(k,1515) + lmat(k,1515)
         mat(k,1516) = mat(k,1516) + lmat(k,1516)
         mat(k,1536) = mat(k,1536) + lmat(k,1536)
         mat(k,1542) = lmat(k,1542)
         mat(k,1557) = mat(k,1557) + lmat(k,1557)
         mat(k,1580) = lmat(k,1580)
         mat(k,1587) = lmat(k,1587)
         mat(k,1700) = mat(k,1700) + lmat(k,1700)
         mat(k,1701) = mat(k,1701) + lmat(k,1701)
         mat(k,1706) = mat(k,1706) + lmat(k,1706)
         mat(k,1711) = mat(k,1711) + lmat(k,1711)
         mat(k,1714) = mat(k,1714) + lmat(k,1714)
         mat(k,1717) = mat(k,1717) + lmat(k,1717)
         mat(k,1732) = mat(k,1732) + lmat(k,1732)
         mat(k,1734) = mat(k,1734) + lmat(k,1734)
         mat(k,1738) = mat(k,1738) + lmat(k,1738)
         mat(k,1745) = mat(k,1745) + lmat(k,1745)
         mat(k,1748) = mat(k,1748) + lmat(k,1748)
         mat(k,1749) = mat(k,1749) + lmat(k,1749)
         mat(k,1751) = mat(k,1751) + lmat(k,1751)
         mat(k,1753) = mat(k,1753) + lmat(k,1753)
         mat(k,1754) = lmat(k,1754)
         mat(k,1755) = mat(k,1755) + lmat(k,1755)
         mat(k,1756) = lmat(k,1756)
         mat(k,1758) = mat(k,1758) + lmat(k,1758)
         mat(k,1759) = mat(k,1759) + lmat(k,1759)
         mat(k,1761) = mat(k,1761) + lmat(k,1761)
         mat(k,1763) = mat(k,1763) + lmat(k,1763)
         mat(k,1764) = mat(k,1764) + lmat(k,1764)
         mat(k,1767) = lmat(k,1767)
         mat(k,1769) = lmat(k,1769)
         mat(k,1781) = mat(k,1781) + lmat(k,1781)
         mat(k,1784) = mat(k,1784) + lmat(k,1784)
         mat(k,1790) = lmat(k,1790)
         mat(k,1805) = mat(k,1805) + lmat(k,1805)
         mat(k,1825) = mat(k,1825) + lmat(k,1825)
         mat(k,1829) = mat(k,1829) + lmat(k,1829)
         mat(k,1830) = lmat(k,1830)
         mat(k,1831) = lmat(k,1831)
         mat(k,1841) = mat(k,1841) + lmat(k,1841)
         mat(k,1844) = mat(k,1844) + lmat(k,1844)
         mat(k,1853) = mat(k,1853) + lmat(k,1853)
         mat(k,1893) = mat(k,1893) + lmat(k,1893)
         mat(k,1896) = mat(k,1896) + lmat(k,1896)
         mat(k,1899) = mat(k,1899) + lmat(k,1899)
         mat(k,1903) = mat(k,1903) + lmat(k,1903)
         mat(k,1950) = mat(k,1950) + lmat(k,1950)
         mat(k,1951) = mat(k,1951) + lmat(k,1951)
         mat(k,1953) = mat(k,1953) + lmat(k,1953)
         mat(k,1957) = mat(k,1957) + lmat(k,1957)
         mat(k,1961) = mat(k,1961) + lmat(k,1961)
         mat(k,1963) = mat(k,1963) + lmat(k,1963)
         mat(k,1996) = mat(k,1996) + lmat(k,1996)
         mat(k,2056) = mat(k,2056) + lmat(k,2056)
         mat(k,2087) = mat(k,2087) + lmat(k,2087)
         mat(k,2089) = mat(k,2089) + lmat(k,2089)
         mat(k,2090) = mat(k,2090) + lmat(k,2090)
         mat(k,2093) = mat(k,2093) + lmat(k,2093)
         mat(k,2099) = mat(k,2099) + lmat(k,2099)
         mat(k,2105) = mat(k,2105) + lmat(k,2105)
         mat(k,2106) = lmat(k,2106)
         mat(k,2118) = mat(k,2118) + lmat(k,2118)
         mat(k,2124) = mat(k,2124) + lmat(k,2124)
         mat(k,2131) = lmat(k,2131)
         mat(k,2139) = lmat(k,2139)
         mat(k,2140) = mat(k,2140) + lmat(k,2140)
         mat(k,2142) = mat(k,2142) + lmat(k,2142)
         mat(k,2144) = lmat(k,2144)
         mat(k,2151) = mat(k,2151) + lmat(k,2151)
         mat(k, 146) = 0._r8
         mat(k, 147) = 0._r8
         mat(k, 250) = 0._r8
         mat(k, 322) = 0._r8
         mat(k, 324) = 0._r8
         mat(k, 337) = 0._r8
         mat(k, 383) = 0._r8
         mat(k, 386) = 0._r8
         mat(k, 401) = 0._r8
         mat(k, 515) = 0._r8
         mat(k, 516) = 0._r8
         mat(k, 523) = 0._r8
         mat(k, 525) = 0._r8
         mat(k, 527) = 0._r8
         mat(k, 528) = 0._r8
         mat(k, 534) = 0._r8
         mat(k, 535) = 0._r8
         mat(k, 539) = 0._r8
         mat(k, 545) = 0._r8
         mat(k, 546) = 0._r8
         mat(k, 549) = 0._r8
         mat(k, 569) = 0._r8
         mat(k, 571) = 0._r8
         mat(k, 573) = 0._r8
         mat(k, 579) = 0._r8
         mat(k, 581) = 0._r8
         mat(k, 583) = 0._r8
         mat(k, 584) = 0._r8
         mat(k, 586) = 0._r8
         mat(k, 599) = 0._r8
         mat(k, 601) = 0._r8
         mat(k, 603) = 0._r8
         mat(k, 604) = 0._r8
         mat(k, 607) = 0._r8
         mat(k, 630) = 0._r8
         mat(k, 632) = 0._r8
         mat(k, 634) = 0._r8
         mat(k, 635) = 0._r8
         mat(k, 637) = 0._r8
         mat(k, 639) = 0._r8
         mat(k, 657) = 0._r8
         mat(k, 671) = 0._r8
         mat(k, 672) = 0._r8
         mat(k, 674) = 0._r8
         mat(k, 690) = 0._r8
         mat(k, 691) = 0._r8
         mat(k, 693) = 0._r8
         mat(k, 700) = 0._r8
         mat(k, 701) = 0._r8
         mat(k, 702) = 0._r8
         mat(k, 707) = 0._r8
         mat(k, 711) = 0._r8
         mat(k, 715) = 0._r8
         mat(k, 731) = 0._r8
         mat(k, 735) = 0._r8
         mat(k, 737) = 0._r8
         mat(k, 742) = 0._r8
         mat(k, 747) = 0._r8
         mat(k, 754) = 0._r8
         mat(k, 755) = 0._r8
         mat(k, 796) = 0._r8
         mat(k, 811) = 0._r8
         mat(k, 823) = 0._r8
         mat(k, 838) = 0._r8
         mat(k, 848) = 0._r8
         mat(k, 851) = 0._r8
         mat(k, 857) = 0._r8
         mat(k, 859) = 0._r8
         mat(k, 860) = 0._r8
         mat(k, 865) = 0._r8
         mat(k, 867) = 0._r8
         mat(k, 879) = 0._r8
         mat(k, 889) = 0._r8
         mat(k, 892) = 0._r8
         mat(k, 898) = 0._r8
         mat(k, 900) = 0._r8
         mat(k, 901) = 0._r8
         mat(k, 906) = 0._r8
         mat(k, 908) = 0._r8
         mat(k, 920) = 0._r8
         mat(k, 921) = 0._r8
         mat(k, 925) = 0._r8
         mat(k, 926) = 0._r8
         mat(k, 927) = 0._r8
         mat(k, 929) = 0._r8
         mat(k, 945) = 0._r8
         mat(k, 947) = 0._r8
         mat(k, 948) = 0._r8
         mat(k, 950) = 0._r8
         mat(k, 951) = 0._r8
         mat(k, 956) = 0._r8
         mat(k, 958) = 0._r8
         mat(k, 964) = 0._r8
         mat(k, 965) = 0._r8
         mat(k, 966) = 0._r8
         mat(k, 969) = 0._r8
         mat(k, 979) = 0._r8
         mat(k, 983) = 0._r8
         mat(k,1007) = 0._r8
         mat(k,1008) = 0._r8
         mat(k,1011) = 0._r8
         mat(k,1015) = 0._r8
         mat(k,1022) = 0._r8
         mat(k,1025) = 0._r8
         mat(k,1030) = 0._r8
         mat(k,1053) = 0._r8
         mat(k,1060) = 0._r8
         mat(k,1068) = 0._r8
         mat(k,1072) = 0._r8
         mat(k,1075) = 0._r8
         mat(k,1088) = 0._r8
         mat(k,1089) = 0._r8
         mat(k,1096) = 0._r8
         mat(k,1099) = 0._r8
         mat(k,1100) = 0._r8
         mat(k,1101) = 0._r8
         mat(k,1102) = 0._r8
         mat(k,1103) = 0._r8
         mat(k,1104) = 0._r8
         mat(k,1106) = 0._r8
         mat(k,1107) = 0._r8
         mat(k,1108) = 0._r8
         mat(k,1112) = 0._r8
         mat(k,1113) = 0._r8
         mat(k,1126) = 0._r8
         mat(k,1136) = 0._r8
         mat(k,1143) = 0._r8
         mat(k,1144) = 0._r8
         mat(k,1145) = 0._r8
         mat(k,1146) = 0._r8
         mat(k,1147) = 0._r8
         mat(k,1148) = 0._r8
         mat(k,1149) = 0._r8
         mat(k,1151) = 0._r8
         mat(k,1153) = 0._r8
         mat(k,1155) = 0._r8
         mat(k,1159) = 0._r8
         mat(k,1160) = 0._r8
         mat(k,1161) = 0._r8
         mat(k,1162) = 0._r8
         mat(k,1167) = 0._r8
         mat(k,1171) = 0._r8
         mat(k,1174) = 0._r8
         mat(k,1176) = 0._r8
         mat(k,1178) = 0._r8
         mat(k,1179) = 0._r8
         mat(k,1181) = 0._r8
         mat(k,1182) = 0._r8
         mat(k,1183) = 0._r8
         mat(k,1184) = 0._r8
         mat(k,1187) = 0._r8
         mat(k,1188) = 0._r8
         mat(k,1189) = 0._r8
         mat(k,1193) = 0._r8
         mat(k,1194) = 0._r8
         mat(k,1195) = 0._r8
         mat(k,1196) = 0._r8
         mat(k,1201) = 0._r8
         mat(k,1208) = 0._r8
         mat(k,1209) = 0._r8
         mat(k,1212) = 0._r8
         mat(k,1216) = 0._r8
         mat(k,1217) = 0._r8
         mat(k,1218) = 0._r8
         mat(k,1223) = 0._r8
         mat(k,1227) = 0._r8
         mat(k,1232) = 0._r8
         mat(k,1235) = 0._r8
         mat(k,1236) = 0._r8
         mat(k,1238) = 0._r8
         mat(k,1239) = 0._r8
         mat(k,1241) = 0._r8
         mat(k,1243) = 0._r8
         mat(k,1245) = 0._r8
         mat(k,1249) = 0._r8
         mat(k,1250) = 0._r8
         mat(k,1256) = 0._r8
         mat(k,1257) = 0._r8
         mat(k,1263) = 0._r8
         mat(k,1269) = 0._r8
         mat(k,1271) = 0._r8
         mat(k,1286) = 0._r8
         mat(k,1289) = 0._r8
         mat(k,1290) = 0._r8
         mat(k,1292) = 0._r8
         mat(k,1296) = 0._r8
         mat(k,1311) = 0._r8
         mat(k,1313) = 0._r8
         mat(k,1320) = 0._r8
         mat(k,1321) = 0._r8
         mat(k,1329) = 0._r8
         mat(k,1330) = 0._r8
         mat(k,1337) = 0._r8
         mat(k,1338) = 0._r8
         mat(k,1339) = 0._r8
         mat(k,1340) = 0._r8
         mat(k,1344) = 0._r8
         mat(k,1354) = 0._r8
         mat(k,1362) = 0._r8
         mat(k,1382) = 0._r8
         mat(k,1384) = 0._r8
         mat(k,1385) = 0._r8
         mat(k,1387) = 0._r8
         mat(k,1388) = 0._r8
         mat(k,1392) = 0._r8
         mat(k,1407) = 0._r8
         mat(k,1409) = 0._r8
         mat(k,1412) = 0._r8
         mat(k,1414) = 0._r8
         mat(k,1415) = 0._r8
         mat(k,1418) = 0._r8
         mat(k,1421) = 0._r8
         mat(k,1422) = 0._r8
         mat(k,1423) = 0._r8
         mat(k,1460) = 0._r8
         mat(k,1495) = 0._r8
         mat(k,1496) = 0._r8
         mat(k,1497) = 0._r8
         mat(k,1504) = 0._r8
         mat(k,1511) = 0._r8
         mat(k,1519) = 0._r8
         mat(k,1520) = 0._r8
         mat(k,1521) = 0._r8
         mat(k,1523) = 0._r8
         mat(k,1524) = 0._r8
         mat(k,1527) = 0._r8
         mat(k,1528) = 0._r8
         mat(k,1538) = 0._r8
         mat(k,1540) = 0._r8
         mat(k,1546) = 0._r8
         mat(k,1553) = 0._r8
         mat(k,1560) = 0._r8
         mat(k,1561) = 0._r8
         mat(k,1569) = 0._r8
         mat(k,1621) = 0._r8
         mat(k,1639) = 0._r8
         mat(k,1648) = 0._r8
         mat(k,1653) = 0._r8
         mat(k,1654) = 0._r8
         mat(k,1679) = 0._r8
         mat(k,1689) = 0._r8
         mat(k,1708) = 0._r8
         mat(k,1735) = 0._r8
         mat(k,1736) = 0._r8
         mat(k,1737) = 0._r8
         mat(k,1739) = 0._r8
         mat(k,1740) = 0._r8
         mat(k,1744) = 0._r8
         mat(k,1762) = 0._r8
         mat(k,1766) = 0._r8
         mat(k,1768) = 0._r8
         mat(k,1773) = 0._r8
         mat(k,1774) = 0._r8
         mat(k,1775) = 0._r8
         mat(k,1776) = 0._r8
         mat(k,1777) = 0._r8
         mat(k,1778) = 0._r8
         mat(k,1779) = 0._r8
         mat(k,1780) = 0._r8
         mat(k,1782) = 0._r8
         mat(k,1783) = 0._r8
         mat(k,1785) = 0._r8
         mat(k,1786) = 0._r8
         mat(k,1787) = 0._r8
         mat(k,1789) = 0._r8
         mat(k,1791) = 0._r8
         mat(k,1794) = 0._r8
         mat(k,1795) = 0._r8
         mat(k,1796) = 0._r8
         mat(k,1798) = 0._r8
         mat(k,1799) = 0._r8
         mat(k,1802) = 0._r8
         mat(k,1803) = 0._r8
         mat(k,1804) = 0._r8
         mat(k,1806) = 0._r8
         mat(k,1808) = 0._r8
         mat(k,1810) = 0._r8
         mat(k,1811) = 0._r8
         mat(k,1818) = 0._r8
         mat(k,1819) = 0._r8
         mat(k,1822) = 0._r8
         mat(k,1824) = 0._r8
         mat(k,1826) = 0._r8
         mat(k,1828) = 0._r8
         mat(k,1833) = 0._r8
         mat(k,1834) = 0._r8
         mat(k,1835) = 0._r8
         mat(k,1838) = 0._r8
         mat(k,1839) = 0._r8
         mat(k,1845) = 0._r8
         mat(k,1847) = 0._r8
         mat(k,1862) = 0._r8
         mat(k,1865) = 0._r8
         mat(k,1870) = 0._r8
         mat(k,1873) = 0._r8
         mat(k,1875) = 0._r8
         mat(k,1876) = 0._r8
         mat(k,1878) = 0._r8
         mat(k,1881) = 0._r8
         mat(k,1884) = 0._r8
         mat(k,1885) = 0._r8
         mat(k,1886) = 0._r8
         mat(k,1888) = 0._r8
         mat(k,1900) = 0._r8
         mat(k,1908) = 0._r8
         mat(k,1915) = 0._r8
         mat(k,1924) = 0._r8
         mat(k,1928) = 0._r8
         mat(k,1929) = 0._r8
         mat(k,1932) = 0._r8
         mat(k,1935) = 0._r8
         mat(k,1947) = 0._r8
         mat(k,1948) = 0._r8
         mat(k,1949) = 0._r8
         mat(k,1952) = 0._r8
         mat(k,1955) = 0._r8
         mat(k,1956) = 0._r8
         mat(k,1958) = 0._r8
         mat(k,1959) = 0._r8
         mat(k,1960) = 0._r8
         mat(k,1965) = 0._r8
         mat(k,1984) = 0._r8
         mat(k,1985) = 0._r8
         mat(k,1986) = 0._r8
         mat(k,2014) = 0._r8
         mat(k,2018) = 0._r8
         mat(k,2020) = 0._r8
         mat(k,2022) = 0._r8
         mat(k,2024) = 0._r8
         mat(k,2027) = 0._r8
         mat(k,2033) = 0._r8
         mat(k,2038) = 0._r8
         mat(k,2050) = 0._r8
         mat(k,2051) = 0._r8
         mat(k,2068) = 0._r8
         mat(k,2071) = 0._r8
         mat(k,2073) = 0._r8
         mat(k,2077) = 0._r8
         mat(k,2078) = 0._r8
         mat(k,2079) = 0._r8
         mat(k,2083) = 0._r8
         mat(k,2084) = 0._r8
         mat(k,2085) = 0._r8
         mat(k,2092) = 0._r8
         mat(k,2094) = 0._r8
         mat(k,2095) = 0._r8
         mat(k,2100) = 0._r8
         mat(k,2101) = 0._r8
         mat(k,2103) = 0._r8
         mat(k,2109) = 0._r8
         mat(k,2110) = 0._r8
         mat(k,2111) = 0._r8
         mat(k,2112) = 0._r8
         mat(k,2115) = 0._r8
         mat(k,2116) = 0._r8
         mat(k,2120) = 0._r8
         mat(k,2123) = 0._r8
         mat(k,2130) = 0._r8
         mat(k,2132) = 0._r8
         mat(k,2133) = 0._r8
         mat(k,2134) = 0._r8
         mat(k,2135) = 0._r8
         mat(k,2136) = 0._r8
         mat(k,2137) = 0._r8
         mat(k,2138) = 0._r8
         mat(k,2141) = 0._r8
         mat(k,2143) = 0._r8
         mat(k,2145) = 0._r8
         mat(k,2146) = 0._r8
         mat(k,2147) = 0._r8
         mat(k,2148) = 0._r8
         mat(k,2149) = 0._r8
         mat(k,2150) = 0._r8
         mat(k, 1) = mat(k, 1) - dti
         mat(k, 2) = mat(k, 2) - dti
         mat(k, 3) = mat(k, 3) - dti
         mat(k, 4) = mat(k, 4) - dti
         mat(k, 5) = mat(k, 5) - dti
         mat(k, 6) = mat(k, 6) - dti
         mat(k, 7) = mat(k, 7) - dti
         mat(k, 8) = mat(k, 8) - dti
         mat(k, 9) = mat(k, 9) - dti
         mat(k, 10) = mat(k, 10) - dti
         mat(k, 11) = mat(k, 11) - dti
         mat(k, 12) = mat(k, 12) - dti
         mat(k, 13) = mat(k, 13) - dti
         mat(k, 14) = mat(k, 14) - dti
         mat(k, 15) = mat(k, 15) - dti
         mat(k, 16) = mat(k, 16) - dti
         mat(k, 17) = mat(k, 17) - dti
         mat(k, 18) = mat(k, 18) - dti
         mat(k, 19) = mat(k, 19) - dti
         mat(k, 20) = mat(k, 20) - dti
         mat(k, 21) = mat(k, 21) - dti
         mat(k, 22) = mat(k, 22) - dti
         mat(k, 23) = mat(k, 23) - dti
         mat(k, 24) = mat(k, 24) - dti
         mat(k, 25) = mat(k, 25) - dti
         mat(k, 26) = mat(k, 26) - dti
         mat(k, 27) = mat(k, 27) - dti
         mat(k, 28) = mat(k, 28) - dti
         mat(k, 29) = mat(k, 29) - dti
         mat(k, 30) = mat(k, 30) - dti
         mat(k, 31) = mat(k, 31) - dti
         mat(k, 32) = mat(k, 32) - dti
         mat(k, 33) = mat(k, 33) - dti
         mat(k, 34) = mat(k, 34) - dti
         mat(k, 40) = mat(k, 40) - dti
         mat(k, 46) = mat(k, 46) - dti
         mat(k, 47) = mat(k, 47) - dti
         mat(k, 50) = mat(k, 50) - dti
         mat(k, 53) = mat(k, 53) - dti
         mat(k, 56) = mat(k, 56) - dti
         mat(k, 60) = mat(k, 60) - dti
         mat(k, 63) = mat(k, 63) - dti
         mat(k, 66) = mat(k, 66) - dti
         mat(k, 69) = mat(k, 69) - dti
         mat(k, 72) = mat(k, 72) - dti
         mat(k, 74) = mat(k, 74) - dti
         mat(k, 77) = mat(k, 77) - dti
         mat(k, 80) = mat(k, 80) - dti
         mat(k, 87) = mat(k, 87) - dti
         mat(k, 93) = mat(k, 93) - dti
         mat(k, 97) = mat(k, 97) - dti
         mat(k, 102) = mat(k, 102) - dti
         mat(k, 104) = mat(k, 104) - dti
         mat(k, 108) = mat(k, 108) - dti
         mat(k, 117) = mat(k, 117) - dti
         mat(k, 124) = mat(k, 124) - dti
         mat(k, 129) = mat(k, 129) - dti
         mat(k, 133) = mat(k, 133) - dti
         mat(k, 142) = mat(k, 142) - dti
         mat(k, 150) = mat(k, 150) - dti
         mat(k, 155) = mat(k, 155) - dti
         mat(k, 160) = mat(k, 160) - dti
         mat(k, 163) = mat(k, 163) - dti
         mat(k, 166) = mat(k, 166) - dti
         mat(k, 169) = mat(k, 169) - dti
         mat(k, 173) = mat(k, 173) - dti
         mat(k, 177) = mat(k, 177) - dti
         mat(k, 181) = mat(k, 181) - dti
         mat(k, 185) = mat(k, 185) - dti
         mat(k, 191) = mat(k, 191) - dti
         mat(k, 194) = mat(k, 194) - dti
         mat(k, 200) = mat(k, 200) - dti
         mat(k, 206) = mat(k, 206) - dti
         mat(k, 212) = mat(k, 212) - dti
         mat(k, 217) = mat(k, 217) - dti
         mat(k, 222) = mat(k, 222) - dti
         mat(k, 227) = mat(k, 227) - dti
         mat(k, 232) = mat(k, 232) - dti
         mat(k, 237) = mat(k, 237) - dti
         mat(k, 242) = mat(k, 242) - dti
         mat(k, 248) = mat(k, 248) - dti
         mat(k, 253) = mat(k, 253) - dti
         mat(k, 258) = mat(k, 258) - dti
         mat(k, 261) = mat(k, 261) - dti
         mat(k, 269) = mat(k, 269) - dti
         mat(k, 277) = mat(k, 277) - dti
         mat(k, 283) = mat(k, 283) - dti
         mat(k, 289) = mat(k, 289) - dti
         mat(k, 295) = mat(k, 295) - dti
         mat(k, 301) = mat(k, 301) - dti
         mat(k, 307) = mat(k, 307) - dti
         mat(k, 313) = mat(k, 313) - dti
         mat(k, 321) = mat(k, 321) - dti
         mat(k, 327) = mat(k, 327) - dti
         mat(k, 334) = mat(k, 334) - dti
         mat(k, 340) = mat(k, 340) - dti
         mat(k, 343) = mat(k, 343) - dti
         mat(k, 348) = mat(k, 348) - dti
         mat(k, 355) = mat(k, 355) - dti
         mat(k, 359) = mat(k, 359) - dti
         mat(k, 366) = mat(k, 366) - dti
         mat(k, 373) = mat(k, 373) - dti
         mat(k, 382) = mat(k, 382) - dti
         mat(k, 389) = mat(k, 389) - dti
         mat(k, 397) = mat(k, 397) - dti
         mat(k, 404) = mat(k, 404) - dti
         mat(k, 409) = mat(k, 409) - dti
         mat(k, 415) = mat(k, 415) - dti
         mat(k, 421) = mat(k, 421) - dti
         mat(k, 426) = mat(k, 426) - dti
         mat(k, 432) = mat(k, 432) - dti
         mat(k, 440) = mat(k, 440) - dti
         mat(k, 448) = mat(k, 448) - dti
         mat(k, 456) = mat(k, 456) - dti
         mat(k, 464) = mat(k, 464) - dti
         mat(k, 472) = mat(k, 472) - dti
         mat(k, 476) = mat(k, 476) - dti
         mat(k, 484) = mat(k, 484) - dti
         mat(k, 488) = mat(k, 488) - dti
         mat(k, 497) = mat(k, 497) - dti
         mat(k, 504) = mat(k, 504) - dti
         mat(k, 513) = mat(k, 513) - dti
         mat(k, 522) = mat(k, 522) - dti
         mat(k, 533) = mat(k, 533) - dti
         mat(k, 544) = mat(k, 544) - dti
         mat(k, 555) = mat(k, 555) - dti
         mat(k, 564) = mat(k, 564) - dti
         mat(k, 578) = mat(k, 578) - dti
         mat(k, 589) = mat(k, 589) - dti
         mat(k, 600) = mat(k, 600) - dti
         mat(k, 610) = mat(k, 610) - dti
         mat(k, 618) = mat(k, 618) - dti
         mat(k, 631) = mat(k, 631) - dti
         mat(k, 642) = mat(k, 642) - dti
         mat(k, 651) = mat(k, 651) - dti
         mat(k, 660) = mat(k, 660) - dti
         mat(k, 670) = mat(k, 670) - dti
         mat(k, 678) = mat(k, 678) - dti
         mat(k, 687) = mat(k, 687) - dti
         mat(k, 698) = mat(k, 698) - dti
         mat(k, 706) = mat(k, 706) - dti
         mat(k, 716) = mat(k, 716) - dti
         mat(k, 720) = mat(k, 720) - dti
         mat(k, 725) = mat(k, 725) - dti
         mat(k, 741) = mat(k, 741) - dti
         mat(k, 751) = mat(k, 751) - dti
         mat(k, 759) = mat(k, 759) - dti
         mat(k, 768) = mat(k, 768) - dti
         mat(k, 781) = mat(k, 781) - dti
         mat(k, 797) = mat(k, 797) - dti
         mat(k, 806) = mat(k, 806) - dti
         mat(k, 818) = mat(k, 818) - dti
         mat(k, 828) = mat(k, 828) - dti
         mat(k, 833) = mat(k, 833) - dti
         mat(k, 849) = mat(k, 849) - dti
         mat(k, 870) = mat(k, 870) - dti
         mat(k, 890) = mat(k, 890) - dti
         mat(k, 909) = mat(k, 909) - dti
         mat(k, 922) = mat(k, 922) - dti
         mat(k, 933) = mat(k, 933) - dti
         mat(k, 943) = mat(k, 943) - dti
         mat(k, 963) = mat(k, 963) - dti
         mat(k, 975) = mat(k, 975) - dti
         mat(k, 985) = mat(k, 985) - dti
         mat(k, 994) = mat(k, 994) - dti
         mat(k,1005) = mat(k,1005) - dti
         mat(k,1019) = mat(k,1019) - dti
         mat(k,1031) = mat(k,1031) - dti
         mat(k,1035) = mat(k,1035) - dti
         mat(k,1044) = mat(k,1044) - dti
         mat(k,1055) = mat(k,1055) - dti
         mat(k,1066) = mat(k,1066) - dti
         mat(k,1083) = mat(k,1083) - dti
         mat(k,1105) = mat(k,1105) - dti
         mat(k,1120) = mat(k,1120) - dti
         mat(k,1132) = mat(k,1132) - dti
         mat(k,1152) = mat(k,1152) - dti
         mat(k,1186) = mat(k,1186) - dti
         mat(k,1211) = mat(k,1211) - dti
         mat(k,1231) = mat(k,1231) - dti
         mat(k,1251) = mat(k,1251) - dti
         mat(k,1284) = mat(k,1284) - dti
         mat(k,1300) = mat(k,1300) - dti
         mat(k,1319) = mat(k,1319) - dti
         mat(k,1334) = mat(k,1334) - dti
         mat(k,1379) = mat(k,1379) - dti
         mat(k,1410) = mat(k,1410) - dti
         mat(k,1490) = mat(k,1490) - dti
         mat(k,1515) = mat(k,1515) - dti
         mat(k,1557) = mat(k,1557) - dti
         mat(k,1706) = mat(k,1706) - dti
         mat(k,1734) = mat(k,1734) - dti
         mat(k,1761) = mat(k,1761) - dti
         mat(k,1784) = mat(k,1784) - dti
         mat(k,1805) = mat(k,1805) - dti
         mat(k,1841) = mat(k,1841) - dti
         mat(k,1903) = mat(k,1903) - dti
         mat(k,1961) = mat(k,1961) - dti
         mat(k,2056) = mat(k,2056) - dti
         mat(k,2099) = mat(k,2099) - dti
         mat(k,2124) = mat(k,2124) - dti
         mat(k,2151) = mat(k,2151) - dti
      end do
      end subroutine nlnmat_finit
      subroutine nlnmat( ofl, ofu, chnkpnts, mat, y, rxt, lmat, dti )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: dti
      real(r8), intent(in) :: lmat(chnkpnts,nzcnt)
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
      call nlnmat01( ofl, ofu, chnkpnts, mat, y, rxt )
      call nlnmat02( ofl, ofu, chnkpnts, mat, y, rxt )
      call nlnmat03( ofl, ofu, chnkpnts, mat, y, rxt )
      call nlnmat04( ofl, ofu, chnkpnts, mat, y, rxt )
      call nlnmat05( ofl, ofu, chnkpnts, mat, y, rxt )
      call nlnmat06( ofl, ofu, chnkpnts, mat, y, rxt )
      call nlnmat07( ofl, ofu, chnkpnts, mat, y, rxt )
      call nlnmat08( ofl, ofu, chnkpnts, mat, y, rxt )
      call nlnmat09( ofl, ofu, chnkpnts, mat, y, rxt )
      call nlnmat10( ofl, ofu, chnkpnts, mat, y, rxt )
      call nlnmat_finit( ofl, ofu, chnkpnts, mat, lmat, dti )
      end subroutine nlnmat
      end module mo_nln_matrix
