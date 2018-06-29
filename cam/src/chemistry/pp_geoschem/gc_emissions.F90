!================================================================================================
! This is the "GEOS-Chem" chemistry emissions interface
!================================================================================================
module gc_emissions_mod

  use shr_kind_mod,        only: r8 => shr_kind_r8
  use spmd_utils,          only : masterproc, myCPU=>iam, nCPUs=>npes
  use cam_logfile,         only : iulog
  use cam_abortutils,      only : endrun

  use chem_mods,           only : ntracers
  use chem_mods,           only : tracernames
  use chem_mods,           only : map2gc

  use tracer_data,         only : trfld,trfile
 
  implicit none

  type :: emission
     integer           :: spc_ndx
     real(r8)          :: mw
     real(r8)          :: scalefactor
     character(len=256):: filename
     character(len=16) :: species
     character(len=8)  :: units
     integer                   :: nsectors
     character(len=32),pointer :: sectors(:)
     type(trfld), pointer      :: fields(:)
     type(trfile)              :: file
  end type emission

  private

  public :: gc_emissions_init
  public :: gc_emissions_calc
  public :: gc_emissions_final

  ! Stand-in: emissions
  type(emission), allocatable :: emissions(:)
  integer                     :: n_emis_files

!================================================================================================
contains
!================================================================================================
  
  subroutine gc_emissions_init

  integer :: ierr

  n_emis_files=1
  allocate(emissions(n_emis_files),stat=ierr)
  if (ierr.ne.0) call endrun('Could not allocate GC emissions')

  end subroutine gc_emissions_init

  subroutine gc_emissions_calc(eflx)

  ! Emissions in kg/m2/s
  ! Dimensions: [ N columns x K levels x C constituents ]
  real(r8), intent(out) :: eflx(:,:,:)
  integer :: i_trc, i_emis

  eflx(:,:,:) = 0.0e+0_r8
  do i_emis = 1,n_emis_files
     ! Read emissions file
     do i_trc = 1,ntracers
     end do
  end do

  end subroutine gc_emissions_calc

  subroutine gc_emissions_final
  if (allocated(emissions)) deallocate(emissions)
  end subroutine gc_emissions_final

end module gc_emissions_mod
