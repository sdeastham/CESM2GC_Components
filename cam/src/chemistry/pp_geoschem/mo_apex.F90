module mo_apex
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid,          only: pcols, begchunk, endchunk          ! physics grid
   use cam_abortutils,  only: endrun
   use cam_logfile,     only: iulog
   use spmd_utils,      only: masterproc
   implicit none

   private
   public :: mo_apex_readnl
   public :: mo_apex_init
   public :: alatm, alonm, bnorth, beast, bdown, bmag
   public :: d1vec, d2vec, colatp, elonp
   public :: maglon0 ! geographic longitude at the equator where geomagnetic longitude is zero (radians)

   ! year to initialize apex
   real(r8), public, protected :: geomag_year = -1._r8
   integer :: fixed_geomag_year = -1

!-------------------------------------------------------------------------------
! Magnetic field output arrays, chunked for physics:
! (these are allocated (pcols,begchunk:endchunk) by sub allocate_arrays)
!-------------------------------------------------------------------------------
   real(r8), protected, allocatable, dimension(:,:) :: & ! (pcols,begchunk:endchunk)
     alatm,  & ! apex mag latitude at each geographic grid point (radians)
     alonm,  & ! apex mag longitude at each geographic grid point (radians)
     bnorth, & ! northward component of magnetic field
     beast,  & ! eastward component of magnetic field
     bdown,  & ! downward component of magnetic field
     bmag      ! magnitude of magnetic field
   real(r8), protected, allocatable, dimension(:,:,:) :: & ! (3,pcols,begchunk:endchunk)
     d1vec,    & ! base vectors more-or-less magnetic eastward direction
     d2vec       ! base vectors more-or-less magnetic downward/equatorward direction
   real(r8), protected :: &
     colatp,   & ! geocentric colatitude of geomagnetic dipole north pole (deg)
     elonp	 ! East longitude of geomagnetic dipole north pole (deg)

   real(r8), protected :: maglon0

   character(len=256) :: igrf_geomag_coefs_file = 'igrf_geomag_coefs_file'

contains

!======================================================================
!======================================================================
subroutine mo_apex_readnl(nlfile)

  use namelist_utils, only : find_group_name
  use units,          only : getunit, freeunit
  use spmd_utils,     only : mpicom, masterprocid, mpi_integer, mpi_character

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

end subroutine mo_apex_readnl

!======================================================================
!======================================================================
subroutine mo_apex_init(phys_state)
   !use physconst,only : pi
   use physics_types, only: physics_state
   !use epp_ionization,only: epp_ionization_setmag
   !use time_manager,  only: get_curr_date
   !use dyn_grid,      only: get_horiz_grid_dim_d

   ! Input/output arguments
   type(physics_state), intent(in), dimension(begchunk:endchunk) :: phys_state
end subroutine mo_apex_init

subroutine allocate_arrays
end subroutine allocate_arrays

end module mo_apex
