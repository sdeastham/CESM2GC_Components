module mo_lightning
  !----------------------------------------------------------------------
  ! ... the lightning module
  !----------------------------------------------------------------------

  use shr_kind_mod,      only : r8 => shr_kind_r8
  use ppgrid,            only : begchunk, endchunk, pcols, pver
  use phys_grid,         only : ngcols_p
  use cam_abortutils,    only : endrun
  use cam_logfile,       only : iulog
  use spmd_utils,        only : masterproc, mpicom

  implicit none

  private
  public  :: lightning_inti
  public  :: lightning_no_prod
  public  :: prod_no

  save

  real(r8) :: csrf
  real(r8) :: factor = 0.1_r8              ! user-controlled scaling factor to achieve arbitrary no prod.
  real(r8) :: geo_factor                   ! grid cell area factor
  real(r8) :: vdist(16,3)                  ! vertical distribution of lightning
  real(r8), allocatable :: prod_no(:,:,:)
  real(r8), allocatable :: glob_prod_no_col(:,:)
  real(r8), allocatable :: flash_freq(:,:)
  integer :: no_ndx,xno_ndx
  logical :: has_no_lightning_prod = .false.

contains

  subroutine lightning_inti( lght_no_prd_factor )
    !----------------------------------------------------------------------
    !       ... initialize the lightning module
    !----------------------------------------------------------------------
    !use mo_constants,  only : pi
    use ioFileMod,     only : getfil
    !use mo_chem_utls,  only : get_spc_ndx

    use cam_history,   only : addfld, horiz_only
    use dyn_grid,      only : get_dyn_grid_parm

    implicit none

    !----------------------------------------------------------------------
    !	... dummy args
    !----------------------------------------------------------------------
    real(r8), intent(in) :: lght_no_prd_factor        ! lightning no production factor

    !call addfld( 'LNO_COL_PROD', horiz_only,  'I', 'TG N/YR', 'lighting column NO source' )
    !call addfld( 'LNO_PROD',     (/ 'lev' /), 'I', '/cm3/s',  'lighting insitu NO source' )
    !call addfld( 'FLASHFRQ',     horiz_only,  'I', '1/MIN',   'lighting flash rate' )        ! flash frequency in grid box per minute (PPP)
    !call addfld( 'FLASHENGY',    horiz_only,  'I', '   ',     'lighting flash rate' )          ! flash frequency in grid box per minute (PPP)
    !call addfld( 'CLDHGT',       horiz_only,  'I', 'KM',      'cloud top height' )              ! cloud top height
    !call addfld( 'DCHGZONE',     horiz_only,  'I', 'KM',      'depth of discharge zone' )       ! depth of discharge zone
    !call addfld( 'CGIC',         horiz_only,  'I', 'RATIO',   'ratio of cloud-ground/intracloud discharges' ) ! ratio of cloud-ground/intracloud discharges

  end subroutine lightning_inti

  subroutine lightning_no_prod( state, pbuf2d,  cam_in )
    !----------------------------------------------------------------------
    !	... set no production from lightning
    !----------------------------------------------------------------------
    use physics_types,    only : physics_state
    
    use physics_buffer,   only : pbuf_get_index, physics_buffer_desc, pbuf_get_field, pbuf_get_chunk
    use physconst,        only : rga
    use phys_grid,        only : get_rlat_all_p, get_lat_all_p, get_lon_all_p, get_wght_all_p
    use cam_history,      only : outfld
    use camsrfexch,       only : cam_in_t
    use shr_reprosum_mod, only : shr_reprosum_calc
    !use mo_constants,  only : rearth, d2r
    implicit none

    !----------------------------------------------------------------------
    !	... dummy args
    !----------------------------------------------------------------------
    type(physics_state), intent(in) :: state(begchunk:endchunk) ! physics state
    
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(cam_in_t), intent(in) :: cam_in(begchunk:endchunk) ! physics state

    !----------------------------------------------------------------------
    !	... local variables
    !----------------------------------------------------------------------

    !----------------------------------------------------------------------
    ! 	... parameters to determine cg/ic ratio [price and rind, 1993]
    !----------------------------------------------------------------------

    if (.not.has_no_lightning_prod) return

    ! < === INSERT CALCULATION HERE === >

    !!--------------------------------------------------------------------------------
    !!       ... output lightning no production to history file
    !!--------------------------------------------------------------------------------
    !do c = begchunk,endchunk
    !   lchnk = state(c)%lchnk
    !   call outfld( 'LNO_PROD',     prod_no(:,:,c),        pcols, lchnk )
    !   call outfld( 'LNO_COL_PROD', glob_prod_no_col(:,c), pcols, lchnk )
    !   call outfld( 'FLASHFRQ',     flash_freq(:,c),       pcols, lchnk )
    !   call outfld( 'FLASHENGY',    flash_energy(:,c),     pcols, lchnk )
    !   call outfld( 'CLDHGT',       cldhgt(:,c),           pcols, lchnk )
    !   call outfld( 'DCHGZONE',     dchgzone(:,c),         pcols, lchnk )
    !   call outfld( 'CGIC',         cgic(:,c),             pcols, lchnk )
    !enddo

  end subroutine lightning_no_prod

end module mo_lightning
