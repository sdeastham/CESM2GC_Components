
      module efield
!------------------------------------------------------------------------------ 
! description: calculates the electric potential for a given year,
!      day of year,UT, F10.7, B_z(K_p)
!     - low/midlatitudes electric potential is from an empirical model from
!       L.Scherliess ludger@gaim.cass.usu.edu
!     - high latitude electric potential is from Weimer96 model
!     - the transition zone is smoothed
!     - output is the horizontal global electric field in magnetic coordinates direction
!      at every magnetic local time grid point expressed in degrees (0 deg-0MLT; 360deg 24 MLT)
!
! input 
!      integer :: iday,     ! day number of year
!                 iyear     ! year
!      real(r8):: ut,       ! universal time 
!                 F10.7,    ! solar flux       (see ionosphere module)
!                 bz        ! component of IMF (see ionosphere module)
! output
!      real(r8) ::               &
!       ed1(0:nmlon,0:nmlat),    &  ! zonal electric field Ed1  [V/m] 
!       ed2(0:nmlon,0:nmlat)        ! meridional electric field Ed2/sin I_m  [V/m]  
!
! notes:
!
! - !to be done (commented out): input S_a F10.7/ Kp from WACCM and calculate B_z 
!    from these inputs
! - assume regular geomagnetic grid 
! - uses average year 365.24 days/year 30.6001 day/mo s. Weimer
! - get_tilt works only for iyear >= 1900
! - Weimer model 1996, Dan Weimer (not with the updates from B.Emery)
! - fixed parameters: B_z, B_y units nT  CHANGE THIS
!                     F10.7
! - we assume that the reference height is 300km for the emperical potential model
! - as a first approximation the electric field is constant in height
!   WATCH what is the upper boundary condition in WACCM
! - for all the calculation done here we set the reference height to the same 
!   value as in tiegcm (hr=130km)
! - 12/15/03 input value iseasav : replaced by day -> month and day of month
! - 12/15/03 S_aM calculated according to Scherliess draft paper and added
!   S_aM(corrected) = 90*(S_aM+1) to get variation in fig 1 Scherliess draft
!
! Author: A. Maute Dec 2003  am 12/30/03 
!------------------------------------------------------------------------------ 

      use shr_kind_mod,      only: r8 => shr_kind_r8
      use physconst,         only: pi
      use cam_abortutils,    only: endrun
      use cam_logfile,       only: iulog
   
      implicit none

      public :: efield_init,  & ! interface routine
                get_efield      ! interface routine
      public :: ed1,          & ! zonal electric field Ed1  [V/m] 
                ed2,          & ! meridional electric field Ed2 [V/m] 
                potent,       & ! electric potential [V]
	        nmlon, nmlat, & ! dimension of mag. grid 
                dlatm, dlonm, & ! grid spacing of mag. grid 
	        ylonm, ylatm    ! magnetic longitudes/latitudes (deg)
      private

      integer ::  &
        iday,     &      ! day number of year
        iyear,    &      ! year
        iday_m,   &      ! day of month
        imo              ! month
      real(r8) ::  ut    ! universal time  

!------------------------------------------------------------------------------ 
! solar parameters
!------------------------------------------------------------------------------ 
      real(r8) ::   f107d           ! 10.7 cm solar flux
      real(r8) ::   by              ! By component of IMF [nT]
      real(r8) ::   bz              ! Bz component of IMF [nT]
!------------------------------------------------------------------------------ 
! mag. grid dimensions (assumed resolution of 2deg)
!------------------------------------------------------------------------------ 
      integer, parameter ::  &
      nmlon = 180,       &  ! mlon 
      nmlat = 90,        &  ! mlat
      nmlath= nmlat/2,   &  ! mlat/2
      nmlonh= nmlon/2,   &  ! mlon/2
      nmlonp1 = nmlon+1, &  ! mlon+1 
      nmlatp1 = nmlat+1     ! mlat+1

      real(r8) ::        &
        ylatm(0:nmlat),  &   ! magnetic latitudes (deg)
        ylonm(0:nmlon),  &   ! magnetic longitudes (deg)
        dlonm,	         &   ! delon lon grid spacing
        dlatm		     ! delat lat grid spacing

!------------------------------------------------------------------------------ 
! array on magnetic grid:    
!------------------------------------------------------------------------------ 
      real(r8) ::                &
        potent(0:nmlon,0:nmlat), &  ! electric potential   [V]  
        ed1(0:nmlon,0:nmlat),    &  ! zonal electric field Ed1  [V/m] 
        ed2(0:nmlon,0:nmlat)        ! meridional electric field Ed2/sin I_m  [V/m]  
       
      real(r8) :: &
       day      ! iday+ut

      logical, parameter :: iutav=.false.   ! .true.  means UT-averaging 
                                            ! .false. means no UT-averaging
      real(r8), parameter ::          &
        v_sw = 400._r8 	                   ! solar wind velocity [km/s]

!------------------------------------------------------------------------------ 
! boundary for Weimer
!------------------------------------------------------------------------------ 
      real(r8), parameter :: bnd_wei = 44._r8 ! colat. [deg]
      integer :: nmlat_wei
      
!------------------------------------------------------------------------------ 
! flag for choosing factors for empirical low latitude model      
!------------------------------------------------------------------------------ 
      integer, parameter ::  iseasav = 0  ! flag for season 

!------------------------------------------------------------------------------ 
! constants:
!------------------------------------------------------------------------------ 
      real(r8), parameter ::        & 
     	r_e  =  6.371e6_r8,         &  ! radius_earth [m] (same as for apex.F90)
        h_r  = 130.0e3_r8,          &  ! reference height [m] (same as for apex.F90)
     	dy2yr= 365.24_r8,           &  ! day per avg. year used in Weimer
     	dy2mo= 30.6001_r8  	       ! day per avg. month used in Weimer

      real(r8) :: &
     	rtd ,     &    ! radians -> deg
     	dtr,      &    ! deg -> radians
     	sqr2,     &      
     	dy2rd,    &    ! 2*pi/365.24  average year
     	deg2mlt,  &    ! for mlon to deg
        sinIm_mag(0:nmlat)    ! sinIm

      integer :: jmin, jmax   ! latitude index for interpolation of 
                              ! northward e-field ed2 at mag. equator

!------------------------------------------------------------------------------ 
!  for spherical harmonics
!------------------------------------------------------------------------------ 
      integer, parameter ::   &
     	nm   = 19,     &
     	mm   = 18,     &					
     	nmp  = nm + 1, &					       
     	mmp  = mm + 1	  

      real(r8) :: r(0:nm,0:mm)      ! R_n^m
      real(r8) :: pmopmmo(0:mm)     ! sqrt(1+1/2m)

!------------------------------------------------------------------------------ 
!  index for factors f_m(mlt),f_l(UT),f_-k(d)
!------------------------------------------------------------------------------ 
      integer, parameter :: ni = 1091  ! for n=12 m=-18:18
      integer :: imax                                         ! max number of index
      integer,dimension(0:ni) :: &
     	kf, &
     	lf, &
     	mf, &
     	nf, &
     	jf
      real(r8) :: ft(1:3,0:2)  ! used for f_-k(season,k)

      real(r8) ::  a_klnm(0:ni)        !  A_klm
      real(r8) ::  a_lf(0:ni)          ! A_klmn^lf for minimum &
      real(r8) ::  a_hf(0:ni)          ! A_klmn^hf for maximum

!------------------------------------------------------------------------------ 
! high_latitude boundary
!------------------------------------------------------------------------------ 
      real(r8), parameter ::  & 
     	ef_max  = 0.015_r8,   &  ! max e-field for high latitude boundary location [V/m]
     	lat_sft = 54._r8	 ! shift of highlat_bnd to 54 deg
      integer :: ilat_sft        ! index of shift for high latitude boundary
      integer, parameter :: nmax_sin = 2 ! max. wave number to be represented
      logical, parameter :: debug =.false.

      contains

      subroutine efield_init(efield_lflux_file, efield_hflux_file, efield_wei96_file)
      character(len=*), intent(in) :: efield_lflux_file
      character(len=*), intent(in) :: efield_hflux_file
      character(len=*), intent(in) :: efield_wei96_file
      end subroutine efield_init

      subroutine get_efield
      end subroutine get_efield

      subroutine GlobalElPotential
      end subroutine GlobalElPotential

      subroutine ff( ph, mt, f )                                                    
!-----------------------------------------------------------------------
! Purpose: calculate F for normalized associated Legendre polynomial P_n^m
!          Ref.: Richmond J.Atm.Ter.Phys. 1974
!
! Method:  f_m(phi) = sqrt(2) sin(m phi) m > 0
!                   = 1                  m = 0
!                   = sqrt(2) cos(m phi) m < 0
!
! Author: A. Maute Nov 2003  am 11/18/03
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! dummy arguments
!-----------------------------------------------------------------------
      integer,intent(in)   :: mt
      real(r8),intent(in)  :: ph	! geo. longitude of 0SLT (ut*15)
      real(r8),intent(out) :: f(-mt:mt)

!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
      integer  :: m, mmo
      real(r8) :: sp, cp    

      sp   = sin( ph/rtd )
      cp   = cos( ph/rtd )
      f(0) = 1.e0_r8
                                                                
      f(-1) = sqr2*cp
      f(1)  = sqr2*sp      								 
      do m = 2,mt
        mmo   = m - 1  
        f(m)  = f(-mmo)*sp + cp*f(mmo)
        f(-m) = f(-mmo)*cp - sp*f(mmo)
      end do      

      end subroutine ff                                                                      

      subroutine pnm( ct, p )
!----------------------------------------------------------------------------                                                                   
! Purpose: normalized associated Legendre polynomial P_n^m
!          Ref.: Richmond J.Atm.Ter.Phys. 1974
! Method:
!   P_m^m    = sqrt(1+1/2m)*si*P_m-1^m-1                  m>0
!   P_n^m    = [cos*P_n-1^m - R_n-1^m*P_n-2^m ]/R_n^m     n>m>=0
!   dP/d phi = n*cos*P_n^m/sin-(2*n+1)*R_n^m*P_n-1^m/sin  n>=m>=0
!   R_n^m    = sqrt[ (n^2-m^2)/(4n^2-1) ]
!
! Author: A. Maute Nov 2003  am 11/18/03
!----------------------------------------------------------------------------                                                                   

      implicit none

!-----------------------------------------------------------------------
! dummy arguments
!-----------------------------------------------------------------------
      real(r8), intent(inout) :: ct ! cos(colat)                 
      real(r8), intent(inout) :: p(0:nm,0:mm)

!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
      integer  :: mp, m, n
      real(r8) :: pm2, st

!      ct = min( ct,.99_r8 )		! cos(colat)
      st = sqrt( 1._r8 - ct*ct ) 	! sin(colat)

      p(0,0) = 1._r8  
      do mp = 1,mmp  ! m+1=1,mm+1
        m = mp - 1
	if( m >= 1 ) then
           p(m,m) = pmopmmo(m)*p(m-1,m-1)*st 			
	end if
	pm2 = 0._r8                                                                  
	do n = mp,nm                    ! n=m+1,N
	   p(n,m) = (ct*p(n-1,m) - r(n-1,m)*pm2)/r(n,m)
	   pm2    = p(n-1,m)
        end do
      end do

      end subroutine pnm                                                                         

      subroutine prep_pnm
!----------------------------------------------------------------------------                                                                   
! Purpose: constant factors for normalized associated Legendre polynomial P_n^m
!          Ref.: Richmond J.Atm.Ter.Phys. 1974
!
! Method:
!   PmoPmmo(m) = sqrt(1+1/2m)
!   R_n^m      = sqrt[ (n^2-m^2)/(4n^2-1) ]
!
! Author: A. Maute Nov 2003  am 11/18/03
!----------------------------------------------------------------------------                                                                   

      implicit none                

!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
      integer  :: mp, m, n
      real(r8) :: xms, xns, den

      do mp = 1, mmp            ! m+1 = 1,mm+1                                     
	m = mp - 1                                               
	xms = m*m                                                
	if( mp /= 1 ) then
           pmopmmo(m) = sqrt( 1._r8 + .5_r8/M )
	end if
	do n = m,nm      ! n = m,N                                     
	  xns    = n*n                                       
	  den    = max(4._r8*xns - 1._r8,1._r8)
	  r(n,m) = sqrt( (xns  - xms)/den )
	end do                 
      end do 

      end subroutine prep_pnm                                                                         

      subroutine index_quiet
!----------------------------------------------------------------------------                                                                   
! Purpose: set up index for factors f_m(mlt),f_l(UT),f_-k(d) to
!    describe the electric potential Phi for the empirical model   
!
! Method:
!    Phi = sum_k sum_l sum_m sum_n [ A_klmn * P_n^m *f_m(mlt)*f_l(UT)*f_-k(d)]
!    - since the electric potential is symmetric about the equator
!      n+m odd terms are set zero resp. not used
!    - in the summation for calculation Phi the index have the following
!      range n=1,12 and m=-n,n, k=0,2 l=-2,2
!
! Author: A. Maute Nov 2003  am 11/18/03
!----------------------------------------------------------------------------                                                                   

      implicit none

!----------------------------------------------------------------------------                                                                   
!	... local variables
!----------------------------------------------------------------------------                                                                   
      integer :: i, j, k, l, n, m

      i = 0 	! initialize
      j = 1 
      do k = 2,0,-1
        do l = -2,2
          if( k == 2 .and. abs(l) == 2 ) then
             cycle
          end if
          do n = 1,12
            do m = -18,18 
              if( abs(m) <= n ) then		    !  |m| < n
                if( (((n-m)/2)*2) == (n-m) ) then   ! only n+m even
             	  if( n-abs(m) <= 9 ) then	    ! n-|m| <= 9 why?
             	    kf(i) = 2-k
             	    lf(i) = l
             	    nf(i) = n
             	    mf(i) = m
             	    jf(i) = j
             	    i	  = i + 1	 ! counter
                  end if
                end if
              end if
            end do ! m
          end do ! n
        end do ! l
      end do ! k

      imax = i - 1  
      if(imax /= ni ) then    ! check if imax == ni 
        write(iulog,'(a19,i5,a18,i5)') 'index_quiet: imax= ',imax, &
          ' not equal to ni =',ni 
        call endrun('index_quiet ERROR')
      end if							
      if(debug) write(iulog,*) 'imax=',imax

      end subroutine index_quiet                                                           

      subroutine read_acoef (efield_lflux_file, efield_hflux_file)
      use ioFileMod,     only : getfil
      use units,         only : getunit, freeunit
      character(len=*), intent(in) :: efield_lflux_file
      character(len=*), intent(in) :: efield_hflux_file
      end subroutine read_acoef

      subroutine adj_S_a
      end subroutine adj_S_a

      subroutine constants
      end subroutine constants

      subroutine prep_fk
      end subroutine prep_fk

      subroutine set_fkflfs( fk, fl, fs )
      real(r8), intent(out) ::  &
     	fk(0:2),  &	                ! f_-k(day) 
     	fl(-2:2), &	                ! f_l(ut)  
     	fs(2)		                ! f_s(f10.7) 
      end subroutine set_fkflfs

      subroutine efield_mid( mlat, mlon, pot )
      real(r8), intent(in)  :: mlat, mlon
      real(r8), intent(out) :: pot               ! electric potential (V)
      end subroutine efield_mid                                              

      subroutine prep_weimer
      end subroutine prep_weimer

      subroutine pot_latsmo( pot, idlat )  ! pots == pot_highlats
      integer, intent(in)     :: idlat
      real(r8), intent(inout) :: pot(0:nmlon,0:nmlat)
      end subroutine pot_latsmo

      subroutine pot_latsmo2( pot, idlat ) 
      integer, intent(in)     :: idlat
      real(r8), intent(inout) :: pot(0:nmlon,0:nmlat)
      end subroutine pot_latsmo2

      subroutine pot_lonsmo( pot, idlon ) 
      integer, intent(in)     :: idlon
      real(r8), intent(inout) :: pot(0:nmlon,0:nmlat)
      end subroutine pot_lonsmo

      subroutine highlat_getbnd( ihlat_bnd ) 
      integer, intent(out) :: ihlat_bnd(0:nmlon)
      end subroutine highlat_getbnd

      subroutine bnd_sinus( ihlat_bnd, itrans_width )  
      integer, intent(inout) :: ihlat_bnd(0:nmlon)    ! loaction of boundary
      integer, intent(out)   :: itrans_width(0:nmlon) ! width of transition zone
      end subroutine bnd_sinus

      subroutine highlat_adjust( pot_highlats, pot_highlat, pot_midlat, ihlat_bnd )
      integer, intent(in)     :: ihlat_bnd(0:nmlon)	                  ! boundary mid to high latitude
      real(r8), intent(in)    :: pot_midlat(0:nmlon,0:nmlat)              ! low/mid latitude potentail
      real(r8), intent(inout) :: pot_highlat(0:nmlon,0:nmlat)             ! high_lat potential
      real(r8), intent(inout) :: pot_highlats(0:nmlon,0:nmlat)            ! high_lat potential! smoothed high_lat potential
      end subroutine highlat_adjust

      subroutine interp_poten( pot_highlats, pot_highlat, pot_midlat, &
                               ihlat_bnd, itrans_width ) 
      integer, intent(in)  :: ihlat_bnd(0:nmlon)
      integer, intent(in)  :: itrans_width(0:nmlon)
      real(r8), intent(in) :: pot_highlats(0:nmlon,0:nmlat)
      real(r8), intent(in) :: pot_highlat(0:nmlon,0:nmlat)
      real(r8), intent(in) :: pot_midlat(0:nmlon,0:nmlat)
      end subroutine interp_poten

      subroutine DerivPotential
      end subroutine DerivPotential

   end module efield
