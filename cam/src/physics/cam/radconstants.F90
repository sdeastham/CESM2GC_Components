module radconstants

! provide stubs to allow building with no radiation scheme active

implicit none
private
save

integer, parameter, public :: nswbands = 1
integer, parameter, public :: nlwbands = 1
integer, parameter, public :: idx_sw_diag = 1
integer, parameter, public :: idx_lw_diag = 1
integer, parameter, public :: idx_nir_diag = 1
integer, parameter, public :: idx_uv_diag = 1
integer, parameter, public :: nrh = 1
integer, parameter, public :: ot_length = 1


public :: get_number_sw_bands
public :: get_sw_spectral_boundaries
public :: get_lw_spectral_boundaries
public :: get_ref_solar_band_irrad
public :: get_true_ref_solar_band_irrad
public :: get_ref_total_solar_irrad
public :: get_solar_band_fraction_irrad
public :: radconstants_init
public :: rad_gas_index

integer, public, parameter :: gasnamelength = 1
integer, public, parameter :: nradgas = 1
character(len=gasnamelength), public, parameter :: gaslist(nradgas) &
   = (/' '/)

!===============================================================================
contains
!===============================================================================

integer function rad_gas_index(gasname)

   ! return the index in the gaslist array of the specified gasname

   character(len=*),intent(in) :: gasname

   rad_gas_index = -1
end function rad_gas_index

end module radconstants
