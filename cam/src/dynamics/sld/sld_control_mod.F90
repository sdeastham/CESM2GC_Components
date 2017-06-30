module sld_control_mod

! SLD dycore shared data

use shr_kind_mod, only: r8=>shr_kind_r8
use pmgrid,       only: plat, plon, plev
use pspect,       only: pnmax

implicit none
public
save

real(r8) ::  tmass(plat)  ! Mass integral for each latitude pair
real(r8) ::  tmass0       ! Specified dry mass of atmosphere
real(r8) ::  tmassf       ! Global mass integral
real(r8) ::  qmassf       ! Global moisture integral
real(r8) ::  fixmas       ! Proportionality factor for ps in dry mass fixer
real(r8) ::  qmass1       ! Contribution to global moisture integral (mass
                          !  weighting is based upon the "A" part of the hybrid grid)
real(r8) ::  qmass2       ! Contribution to global moisture integral (mass
                          !  weighting is based upon the "B" part of the hybrid grid)
real(r8) ::  pdela(plon,plev) ! pressure difference between interfaces (pressure
                              !  defined using the "A" part of hybrid grid only)
real(r8) ::  zgsint       ! global integral of geopotential height

integer  :: pcray         ! length of vector register (words)
parameter (pcray=64)

real(r8) :: trig (3*plon/2+1,plat)  ! trigonometric funct values used by fft
integer  :: ifax(19,plat)           ! fft factorization of plon/2

real(r8) :: cnfac                ! Courant num factor(multiply by max |V|)
real(r8) :: cnlim                ! Maximum allowable courant number
real(r8) :: hdfsd2(pnmax)        ! Del^2 mult. for each wave (vort-div)
real(r8) :: hdfst2(pnmax)        ! Del^2 multiplier for each wave (t-q)
real(r8) :: hdfsd4(pnmax)        ! Del^4 mult. for each wave (vort-div)
real(r8) :: hdfst4(pnmax)        ! Del^4 multiplier for each wave (t-q)
real(r8) :: hdiftq(pnmax,plev)   ! Temp-tracer diffusion factors
real(r8) :: hdifzd(pnmax,plev)   ! Vort-div diffusion factors
integer  :: kmnhd4               ! Top level for del^4 diffusion
integer  :: kmxhd2               ! Bottom level for increased del^2 diffusion
integer  :: nindex(plev)         ! Starting index for spectral truncation
integer  :: nmaxhd               ! Maximum two dimensional wave number

! Spectral Namelist variables

real(r8) :: dif2               ! del2 horizontal diffusion coeff.
real(r8) :: dif4               ! del4 horizontal diffusion coeff.
real(r8) :: divdampn           ! Number of days to invoke divergence damper
real(r8) :: eps                ! time filter coefficient. Defaults to 0.06.
integer  :: kmxhdc             ! Number of levels to apply Courant limiter

end module sld_control_mod
