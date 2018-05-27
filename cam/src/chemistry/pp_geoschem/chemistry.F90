!================================================================================================
! This is the "GEOS-Chem" chemistry module.
!================================================================================================

module chemistry
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use physics_types,       only: physics_state, physics_ptend, physics_ptend_init
  use ppgrid,              only: begchunk, endchunk, pcols
  use ppgrid,              only: pver
  use constituents,        only: pcnst, cnst_add
  !use mo_gas_phase_chemdr, only: map2chm
  !use mo_constants,        only: pi
  use shr_const_mod,       only: molw_dryair=>SHR_CONST_MWDAIR
  !use mo_chem_utls,        only : get_spc_ndx
  !use chem_mods,           only : gas_pcnst, adv_mass
  !use mo_sim_dat, only: set_sim_dat
  use spmd_utils,          only : masterproc, myCPU=>iam, nCPUs=>npes
  use cam_logfile,         only : iulog

  use input_opt_mod,       only : optinput
  use state_met_mod,       only : metstate
  use state_chm_mod,       only : chmstate

  ! GEOS-Chem precision specifiers
  use precision_mod,       only : fp, f4
 
  ! GEOS-Chem error codes
  use errcode_mod

  ! Exit routine in CAM
  use cam_abortutils, only : endrun

  implicit none
  private
  save
  !
  ! Public interfaces
  !
  public :: chem_is                        ! identify which chemistry is being used
  public :: chem_register                  ! register consituents
  public :: chem_is_active                 ! returns true if this package is active (ghg_chem=.true.)
  public :: chem_implements_cnst           ! returns true if consituent is implemented by this package
  public :: chem_init_cnst                 ! initialize mixing ratios if not read from initial file
  public :: chem_init                      ! initialize (history) variables
  public :: chem_timestep_tend             ! interface to tendency computation
  public :: chem_final
  public :: chem_write_restart
  public :: chem_read_restart
  public :: chem_init_restart
  public :: chem_readnl                    ! read chem namelist 

  public :: chem_emissions
  public :: chem_timestep_init

  ! Private data
  !===== SDE DEBUG =====
  integer, parameter :: ntracersmax = 200    ! Must be equal to nadv_chem
  integer            :: ntracers
  character(len=255) :: tracernames(ntracersmax)
  integer            :: indices(ntracersmax)
  real(r8)           :: adv_mass(ntracersmax)
  real(r8)           :: ref_mmr(ntracersmax)

  ! Short-lived species (i.e. not advected)
  integer, parameter :: nslsmax = 500        ! UNadvected species only
  integer            :: nsls    
  character(len=255) :: slsnames(nslsmax)
  !===== SDE DEBUG =====

  ! Location of valid input.geos
  character(len=500) :: inputGeosPath

  ! Location of chemistry input (for now)
  character(len=500) :: chemInputsDir

  ! Mapping between constituents and GEOS-Chem tracers
  integer :: map2gc(pcnst)

  ! GEOS-Chem state variables
  Type(OptInput)                 :: Input_Opt
  Type(MetState),Allocatable     :: State_Met(:)
  Type(ChmState),Allocatable     :: State_Chm(:)

  ! Indices of critical species
  integer :: iH2O

  ! Indices in the physics buffer
  integer :: ndx_pblh

!================================================================================================
contains
!================================================================================================

  logical function chem_is (name)

    character(len=*), intent(in) :: name

    chem_is = .false.
    if (name == 'geoschem' ) then
       chem_is = .true.
    end if
    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_IS'

  end function chem_is

!================================================================================================

  subroutine chem_register

    use physics_buffer, only : pbuf_add_field, dtype_r8
    use physconst,      only : mwdry

    use state_chm_mod,  only : init_state_chm, cleanup_state_chm
    use state_chm_mod,  only : Ind_
    use input_opt_mod,  only : set_input_opt,  cleanup_input_opt
    USE Species_Mod,    only : species
    
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: register advected constituents for chemistry
    ! 
    !-----------------------------------------------------------------------
    ! Need to generate a temporary species database - therefore temp State_Chm
    Type(ChmState)         :: SC
    Type(OptInput)         :: IO
    TYPE(Species), POINTER :: ThisSpc

    integer            :: i, n
    real(r8)           :: cptmp
    real(r8)           :: qmin
    real(r8)           :: ref_vmr
    character(len=128) :: mixtype
    character(len=128) :: molectype
    character(len=128) :: lng_name
    logical            :: camout
    logical            :: ic_from_cam2
    logical            :: has_fixed_ubc
    logical            :: has_fixed_ubflx

    integer            :: rc
    ! SDE 2018-05-02: This seems to get called before anything else
    ! That includes CHEM_INIT
    ! At this point, mozart calls SET_SIM_DAT, which is specified by each
    ! mechanism separately (ie mozart/chemistry.F90 calls the subroutine
    ! set_sim_dat which is in pp_[mechanism]/mo_sim_dat.F90. That sets a lot of
    ! data in other places, notably in "chem_mods"

    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_REGISTER'

    ! Generate fake state_chm
    IO%Max_Diag            = 1000
    IO%Max_Trcs            = 500
    IO%Max_Memb            = 15
    IO%Max_Fams            = 250
    IO%Max_Dep             = 500

    IO%RootCPU             = .False.
    Call Set_Input_Opt(.False.,IO,RC)
    If (rc.ne.GC_SUCCESS) Call endrun('Could not generate reference input options object')

    ! Options needed by Init_State_Chm
    IO%ITS_A_FULLCHEM_SIM  = .True.
    IO%LUCX                = .True.
    IO%LPRT                = .False.
    IO%N_Advect            = nTracers
    Do I=1,nTracers
       IO%AdvectSpc_Name(I) = Trim(tracernames(I))
    End Do
    IO%SalA_rEdge_um(1)    = 0.01e+0_fp
    IO%SalA_rEdge_um(2)    = 0.50e+0_fp
    IO%SalC_rEdge_um(1)    = 0.50e+0_fp
    IO%SalC_rEdge_um(2)    = 8.00e+0_fp

    ! Prevent reporting
    IO%rootCPU             = .False.
    Input_Opt%myCPU        = myCPU
    Call Init_State_Chm( .False., 1, 1, 1, &
                         IO, SC, 10, RC )
    If (rc.ne.GC_SUCCESS) Call endrun('Could not generate reference species database')

    ! At the moment, we force nadv_chem=200 in the setup file
    ! Default
    map2gc = -1
    do i = 1, ntracersmax
       if (i.le.ntracers) then
          n = Ind_(tracernames(i))
          ThisSpc => SC%SpcData(N)%Info
          lng_name    = Trim(ThisSpc%FullName)
          adv_mass(i) = real(ThisSpc%MW_g,r8)
          ref_vmr     = real(ThisSpc%BackgroundVV,r8)
          ref_mmr(i)  = ref_vmr / (mwdry / adv_mass(i))
       else
          lng_name = trim(tracernames(i))
          adv_mass(i) = 1000.0e+0_r8 * (0.001e+0_r8)
          ref_mmr(i)  = 1.0e-38_r8
       endif
       ! dummy value for specific heat of constant pressure (Cp)
       cptmp = 666._r8
       ! minimum mixing ratio
       qmin = 1.e-38_r8
       ! mixing ratio type
       mixtype = 'dry'
       ! Used for ionospheric WACCM (WACCM-X)
       molectype = 'minor'
       ! Is an output field (?)
       camout = .false.
       ! Not true for O2(1-delta) or O2(1-sigma)
       ic_from_cam2  = .true.
       ! Use a fixed value at the upper boundary
       has_fixed_ubc = .false.
       ! Use a fixed flux condition at the upper boundary
       has_fixed_ubflx = .false.
       !write(tracernames(i),'(a,I0.4)') 'GCTRC_', i
       ! NOTE: In MOZART, this only gets called for tracers
       ! This is the call to add a "constituent"
       call cnst_add( trim(tracernames(i)), adv_mass(i), cptmp, qmin, n, &
                      readiv=ic_from_cam2, mixtype=mixtype, cam_outfld=camout, &
                      molectype=molectype, fixed_ubc=has_fixed_ubc, &
                      fixed_ubflx=has_fixed_ubflx, longname=trim(lng_name) )

       ! Add to GC mapping. When starting a timestep, we will want to update the
       ! concentration of State_Chm(x)%Species(1,iCol,iLev,i) with data from
       ! constituent n
       if (i.le.ntracers) map2gc(n) = i
       ! Nullify pointer
       ThisSpc => NULL()
    end do

    ! Now unadvected species
    ! MOZART uses this for short-lived species. Not certain exactly what it
    ! does, but note that the "ShortLivedSpecies" physics buffer already
    ! needs to have been initialized, which we haven't done. Physics buffers
    ! are fields which are available either across timesteps or for use to
    ! modules outside of chemistry
    ! More information:
    ! http://www.cesm.ucar.edu/models/atm-cam/docs/phys-interface/node5.html
    !call pbuf_add_field('ShortLivedSpecies','global',dtype_r8,(/pcols,pver,nslvd/),pbf_idx)
    ! returned values
    !  n : mapping in CAM
    ! map2chm is a mozart variable
    !map2chm(n) = i
    !indices(i) = 0
    ! ===== SDE DEBUG =====

    ! Clean up
    Call Cleanup_State_Chm( .False., SC, RC )
    Call Cleanup_Input_Opt( .False., IO, RC )

  end subroutine chem_register

  subroutine chem_readnl(nlfile)
    ! This is the FIRST routine to get called - so it should read in 
    ! GEOS-Chem options from input.geos without actually doing any 
    ! initialization

    use cam_abortutils, only : endrun
    use units,          only : getunit, freeunit
    use mpishorthand
    use gckpp_Model,   only : nspec, spc_names

    ! args
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: i, unitn, ierr
    character(len=500) :: line
    logical :: menuFound, validSLS

    ! Set paths
    !inputGeosPath='/n/regal/jacob_lab/seastham/CESM2/CESM2_GC2/ut_src/runs/4x5_standard/input.geos.template'
    !chemInputsDir='/n/holylfs/EXTERNAL_REPOS/GEOS-CHEM/gcgrid/gcdata/ExtData/CHEM_INPUTS/'
    inputGeosPath='/net/d05/data/SDEData/CESM2/ut_src/runs/4x5_standard/input.geos.template'
    chemInputsDir='/net/d06/data/GCdata/ExtData/CHEM_INPUTS/'

    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_READNL'

    ! TODO: Read in input.geos and get species names
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(inputGeosPath), status='old', iostat=ierr )
       if (ierr .ne. 0) then
          call endrun('chem_readnl: ERROR opening input.geos')
       end if
       ! Go to ADVECTED SPECIES MENU
       menuFound = .False.
       Do While (.not.menuFound)
          read( unitn, '(a)', iostat=ierr ) line
          if (ierr.ne.0) then
             call endrun('chem_readnl: ERROR finding advected species menu')
          else if (index(line,'ADVECTED SPECIES MENU') > 0) then
             menuFound=.True.
          end if
       end do
       ! Skip first line
       read(unitn,'(a)',iostat=ierr) line
       ! Read in tracer count
       read(unitn,'(26x,I)',iostat=ierr) ntracers
       ! Skip divider line
       read(unitn,'(a)',iostat=ierr) line
       ! Read in each tracer
       do i=1,ntracers
          read(unitn,'(26x,a)',iostat=ierr) line
          tracernames(i) = trim(line)
       end do
       close(unitn)
       call freeunit(unitn)
       ! Assign remaining tracers dummy names
       do i=(ntracers+1),ntracersmax
          write(tracernames(i),'(a,I0.4)') 'GCTRC_',i
       end do

       ! Now go through the KPP mechanism and add any species not implemented by
       ! the tracer list in input.geos
       if ( nspec > nslsmax ) then
          call endrun('chem_readnl: too many species - increase nslsmax')
       end If
       nsls = 0
       do i=1,nspec
          ! Get the name of the species from KPP
          line = adjustl(trim(spc_names(i)))
          ! Only add this 
          validSLS = ( (.not.any(trim(line).eq.tracernames)).and.&
                       (.not.(line(1:2) == 'RR')) )
          if (validSLS) then
             ! Genuine new short-lived species
             nsls = nsls + 1
             slsnames(nsls) = trim(line)
             write(iulog,'(a,I5,a,a)') ' --> GC species ',nsls, ': ',trim(line)
          end if
       end do
    end if

    ! Broadcast to all processors
#if defined( SPMD )
    call mpibcast(ntracers,    1,                               mpiint,  0, mpicom )
    call mpibcast(tracernames, len(tracernames(1))*ntracersmax, mpichar, 0, mpicom )
    call mpibcast(nsls,        1,                               mpiint,  0, mpicom )
    call mpibcast(slsnames,    len(slsnames(1))*nslsmax,        mpichar, 0, mpicom )
#endif

  end subroutine chem_readnl

!================================================================================================

  function chem_is_active()
    !-----------------------------------------------------------------------
    logical :: chem_is_active
    !-----------------------------------------------------------------------
    chem_is_active = .true.
  end function chem_is_active

!================================================================================================

  function chem_implements_cnst(name)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: return true if specified constituent is implemented by this package
    ! 
    ! Author: B. Eaton
    ! 
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------Arguments---------------------------------

    character(len=*), intent(in) :: name   ! constituent name
    logical :: chem_implements_cnst        ! return value

    integer :: i
    
    chem_implements_cnst = .false.

    do i = 1, ntracers
       if (trim(tracernames(i)) .eq. trim(name)) then
          chem_implements_cnst = .true.
          exit
       end if
    end do

    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_IMPLEMENTS_CNST'
  end function chem_implements_cnst

!===============================================================================
  
  subroutine chem_init(phys_state, pbuf2d)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: initialize GEOS-Chem parts (state objects, mainly)
    !          (and declare history variables)
    ! 
    !-----------------------------------------------------------------------
    use physics_buffer, only: physics_buffer_desc, pbuf_get_index
    use cam_history,    only: addfld, add_default, horiz_only

    use mpishorthand
    use cam_abortutils, only : endrun

    ! For getting grid area
    use physconst,    only : rearth
    use phys_grid,    only : get_area_all_p

    ! Use GEOS-Chem versions of physical constants
    use physconstants,  only : pi, pi_180

    use hycoef,         only : ps0, hyai, hybi

    use input_opt_mod
    use state_met_mod
    use state_chm_mod
    use gc_environment_mod

    use time_mod,      only : accept_external_date_time
    use time_mod,      only : Set_Begin_Time,   Set_End_Time
    use time_mod,      only : Set_Current_Time, Set_DiagB
    use time_mod,      only : Set_NDiagTime,    Get_Tau
    use gc_grid_mod,   only : Set_XOffset, Set_YOffset
    use gc_grid_mod,   only : compute_grid, init_grid
    use gc_grid_mod,   only : SetGridFromCtr
    !use transfer_mod,  only : init_transfer

    use cmn_o3_mod
    use cmn_size_mod

    use linoz_mod,     only : linoz_read
    use ucx_mod,       only : cfcyear
    use ucx_mod,       only : P_Ice_Supersat, T_NAT_Supercool

    use error_mod,     only : init_error
    use drydep_mod,    only : init_drydep
    use carbon_mod,    only : init_carbon
    use dust_mod,      only : init_dust
    use seasalt_mod,   only : init_seasalt
    use sulfate_mod,   only : init_sulfate
    use aerosol_mod,   only : init_aerosol
    use linoz_mod,     only : init_linoz
    use WetScav_Mod,   Only : init_wetscav
    use chemgrid_mod,  only : init_chemgrid
    use TOMS_mod,      only : init_TOMS
    use C2H6_mod,      only : init_C2H6
    use pressure_mod,  only : init_pressure
    use chemistry_mod, only : init_chemistry
    !use MODIS_LAI_Mod, Only : init_MODIS_LAI
    use ucx_mod,       only : init_ucx

    use pbl_mix_mod,   only : init_pbl_mix

    use state_chm_mod,    only: Ind_

    type(physics_state), intent(in):: phys_state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    ! Local variables
    integer               :: lchnk(begchunk:endchunk), ncol(begchunk:endchunk)
    integer               :: iwait, ierr

    integer               :: nX, nY, nZ
    integer               :: iX, iY, iZ
    integer               :: nlev, i, j, l, rc
    integer               :: nLinoz

    ! Grid setup
    !real(fp), allocatable :: DLonVec(:), DLatVec(:)
    real(fp)              :: lonVal, latVal, dLonFix, dLatFix
    real(fp), allocatable :: lonMidArr(:,:), latMidArr(:,:)
    real(r8), allocatable :: linozData(:,:,:,:)

    real(r8), allocatable :: col_area(:)
    real(fp), allocatable :: Ap_CAM_Flip(:), Bp_CAM_Flip(:)

    logical               :: rootCPU, rootChunk

    ! lchnk: which chunks we have on this process
    lchnk = phys_state%lchnk
    ! ncol: number of atmospheric columns for each chunk
    ncol  = phys_state%ncol
    ! nlev: number of vertical levels
    nlev  = pver

    write(iulog,'(2(a,x,I6,x))') 'chem_init called on PE ', myCPU, ' of ', nCPUs

    ! The GEOS-Chem grids on every "chunk" will all be the same size, otherwise
    ! the possibility of differently-sized chunks will become a huge headache
    nX = 1
    nY = MaxVal(ncol)
    nZ = nlev

    ! This ensures that each process allocates everything needed for its chunks
    Allocate(State_Met(begchunk:endchunk))
    Allocate(State_Chm(begchunk:endchunk))
   
    ! Set some basic flags
    Input_Opt%Max_Diag          = 1000
    Input_Opt%Max_Trcs          = 500
    Input_Opt%Max_Memb          = 15
    Input_Opt%Max_Fams          = 250
    Input_Opt%Max_Dep           = 500

    ! Note - this is called AFTER chem_readnl, after X, and after 
    ! every constituent has had its initial conditions read. Any
    ! constituent which is not found in the CAM restart file will
    ! then have already had a call to chem_implements_cnst, and will
    ! have then had a call to chem_init_cnst to set a default VMR
    call GC_Allocate_All ( am_I_Root      = masterproc,&
                           Input_Opt      = Input_Opt, &
                           value_I_Lo     = 1,         &
                           value_J_Lo     = 1,         &
                           value_I_Hi     = nX,        &
                           value_J_Hi     = nY,        &
                           value_IM       = nX,        &
                           value_JM       = nY,        &
                           value_LM       = nZ,        &
                           value_IM_WORLD = nX,        &
                           value_JM_WORLD = nY,        &
                           value_LM_WORLD = nZ,        &
                           RC             = RC )
 
    Input_Opt%myCPU             = myCPU
    Input_Opt%rootCPU           = masterproc

    ! TODO: Mimic GEOS-Chem's reading of input options
    !If (masterproc) then
    !   Call Read_Input_File( am_I_Root   = .True., &
    !                         Input_Opt   = Input_Opt, &
    !                         srcFile     = inputGeosPath,      &
    !                         RC          = RC )
    !End If
    !Call <broadcast data to other CPUs>

    ! For now just hard-code it
    ! First set up directories
    Input_Opt%Chem_Inputs_Dir      = Trim(chemInputsDir)

    ! Simulation menu
    ! Ignore the data directories for now
    Input_Opt%NYMDb                = 20000101
    Input_Opt%NHMSb                =   000000
    Input_Opt%NYMDe                = 20010101
    Input_Opt%NHMSe                =   000000
    Input_Opt%LUnzip               = .False. 
    Input_Opt%LWait                = .False.
    Input_Opt%LVarTrop             = .True.
    Input_Opt%Its_A_Nested_Grid    = .False.
    Input_Opt%Nested_I0            = 0
    Input_Opt%Nested_J0            = 0

    ! Now READ_ADVECTED_SPECIES_MENU
    Input_Opt%ITS_A_RnPbBe_SIM       = .False.
    Input_Opt%ITS_A_CH3I_SIM         = .False.
    Input_Opt%ITS_A_FULLCHEM_SIM     = .True.
    Input_Opt%ITS_A_HCN_SIM          = .False.
    Input_Opt%ITS_A_TAGO3_SIM        = .False.
    Input_Opt%ITS_A_TAGCO_SIM        = .False.
    Input_Opt%ITS_A_C2H6_SIM         = .False.
    Input_Opt%ITS_A_CH4_SIM          = .False.
    Input_Opt%ITS_A_CH4_SIM          = .False.
    Input_Opt%ITS_A_RnPbBe_SIM       = .False.
    Input_Opt%ITS_NOT_COPARAM_OR_CH4 = .True.
    Input_Opt%ITS_AN_AEROSOL_SIM     = .False.
    Input_Opt%ITS_A_MERCURY_SIM      = .False.
    Input_Opt%ITS_A_CO2_SIM          = .False.
    Input_Opt%ITS_A_H2HD_SIM         = .False.
    Input_Opt%ITS_A_POPS_SIM         = .False.
    Input_Opt%ITS_A_SPECIALTY_SIM    = .False.
    Input_Opt%SIM_NAME               = 'NOx-Ox-Hydrocarbon-Aerosol'
    Input_Opt%N_Advect               = NTracers
    If (Input_Opt%N_Advect.gt.Input_Opt%Max_Trcs) Then
       call endrun('Number of tracers exceeds max count')
    End If
    ! Assign tracer names
    Do j=1,Input_Opt%N_Advect
       Input_Opt%AdvectSpc_Name(j) = Trim(tracerNames(j))
    End Do
    ! No tagged species
    Input_Opt%LSplit = .False.

    ! when iam = -1, ntracers isn't received...
    write(iulog,'(a,3(x,I6))') ' NTRAC: ', Input_Opt%myCPU, Input_Opt%N_Advect, NTracers
    !if (ntracers<1) call endrun('Tracer count missing?')

    ! Now READ_AEROSOL_MENU
    Input_Opt%LSulf                  = .True.
    Input_Opt%LCryst                 = .False.
    Input_Opt%LCarb                  = .True.
    Input_Opt%LBrC                   = .False.
    Input_Opt%LSOA                   = .True.
    Input_Opt%LSVPOA                 = .False.
    Input_Opt%LDust                  = .True.
    Input_Opt%LDstUp                 = .False.
    Input_Opt%LSSalt                 = .True.
    Input_Opt%SalA_rEdge_um(1)       = 0.01e+0_fp
    Input_Opt%SalA_rEdge_um(2)       = 0.50e+0_fp
    Input_Opt%SalC_rEdge_um(1)       = 0.50e+0_fp
    Input_Opt%SalC_rEdge_um(2)       = 8.00e+0_fp
    Input_Opt%LMPOA                  = .False.
    Input_Opt%LDicarb                = .False.
    Input_Opt%LGravStrat             = .True.
    Input_Opt%LSolidPSC              = .True.
    Input_Opt%LHomNucNAT             = .False.
    Input_Opt%T_NAT_Supercool        = 3.0e+0_fp
    Input_Opt%P_Ice_Supersat         = 1.2e+0_fp
    Input_Opt%LPSCChem               = .True.
    Input_Opt%LStratOD               = .True.

    ! Update options in UCX_mod for this
    P_Ice_Supersat  = Input_Opt%P_Ice_Supersat
    T_NAT_Supercool = Input_Opt%T_NAT_Supercool

    ! Now READ_EMISSIONS_MENU
    Input_Opt%LEmis                  = .False.
    Input_Opt%HCOConfigFile          = 'HEMCO_Config.rc'
    Input_Opt%LFix_PBL_Bro           = .False. 

    ! Set surface VMRs - turn this off so that CAM can handle it
    Input_Opt%LCH4Emis               = .False. 
    Input_Opt%LCH4SBC                = .False.
    Input_Opt%LOCSEmis               = .False.
    Input_Opt%LCFCEmis               = .False.
    Input_Opt%LClEmis                = .False.
    Input_Opt%LBrEmis                = .False.
    Input_Opt%LN2OEmis               = .False.
    Input_Opt%LBasicEmis             = .False.

    ! Set initial conditions
    Input_Opt%LSetH2O                = .True.
    Input_Opt%LSetCH4                = .False.
    Input_Opt%LSetOCS                = .False.
    Input_Opt%LSetCFC                = .False.
    Input_Opt%LSetCl                 = .False.
    Input_Opt%LBrGCCM                = .False.
    Input_Opt%LSetBr                 = .False.
    Input_Opt%LSetBrStrat            = .False.
    Input_Opt%LSetNOyStrat           = .False.
    Input_Opt%LSetN2O                = .False.
    Input_Opt%LSetH2SO4              = .False.

    ! CFC control
    Input_Opt%CFCYear                = 0
    Input_Opt%LFutureCFC             = .False.

    ! Transmit to UCX_mod
    CFCYear = Input_Opt%CFCYear

    ! Now READ_CHEMISTRY_MENU
    Input_Opt%LChem                  = .True.
    Input_Opt%LSChem                 = .True.
    Input_Opt%LLinoz                 = .True.
    Input_Opt%LUCX                   = .True.
    Input_Opt%LCH4Chem               = .True.
    Input_Opt%LActiveH2O             = .True.
    Input_Opt%LO3FJX                 = .True.
    Input_Opt%Gamma_HO2              = 0.2e+0_fp

    ! Expect to get total overhead ozone, although it shouldn't make too much
    ! of a difference since we want to use "full-UCX"
    Input_Opt%Use_O3_From_Met        = .True.

    ! Now READ_TRANSPORT_MENU
    Input_Opt%LTran                  = .True.
    Input_Opt%LFill                  = .True.
    Input_Opt%TPCore_IOrd            = 3
    Input_Opt%TPCore_JOrd            = 3
    Input_Opt%TPCore_KOrd            = 3

    ! Now READ_CONVECTION_MENU
    ! For now,  no vertical motions in GC
    Input_Opt%LConv                  = .False.
    Input_Opt%LTurb                  = .False.
    Input_Opt%LNLPBL                 = .False.

    ! Now READ_DEPOSITION_MENU
    ! Disable dry/wet dep for now
    Input_Opt%LDryD                  = .False.
    Input_Opt%LWetD                  = .False.
    Input_Opt%Use_Olson_2001         = .True.

    ! Read in data for Linoz. All CPUs allocate one array to hold the data. Only
    ! the root CPU reads in the data; then we copy it out to a temporary array,
    ! broadcast to all other CPUs, and finally duplicate the data into every
    ! copy of Input_Opt
    If (Input_Opt%LLinoz) Then
       ! Allocate array for broadcast
       nLinoz = Input_Opt%Linoz_NLevels * &
                Input_Opt%Linoz_NLat    * &
                Input_Opt%Linoz_NMonths * &
                Input_Opt%Linoz_NFields 
       Allocate( linozData( Input_Opt%Linoz_NLevels,               &
                            Input_Opt%Linoz_NLat,                  &
                            Input_Opt%Linoz_NMonths,               &
                            Input_Opt%Linoz_NFields  ), Stat=ierr)
       If (ierr.ne.0) Call endrun('Failure while allocating linozData')
       linozData = 0.0e+0_r8
  
       If (masterproc) Then
          ! Read data in to Input_Opt%Linoz_TParm
          Call Linoz_Read( masterproc, Input_Opt, RC )
          ! Copy the data to a temporary array
          linozData = real(Input_Opt%LINOZ_TPARM,r8)
       End If
#if defined( SPMD )
       call mpibcast(linozData, nLinoz, mpir8, 0, mpicom )
#endif
       ! Now copy the data to all other Input_Opt copies
       if (.not.masterproc) then
          Input_Opt%LINOZ_TPARM = real(linozData,fp)
       end if
       Deallocate(linozData)
    End If

    ! Set the times held by "time_mod"
    Call Accept_External_Date_Time( am_I_Root   = masterproc,         &
                                    value_NYMDb = Input_Opt%NYMDb, &
                                    value_NHMSb = Input_Opt%NHMSb, &
                                    value_NYMDe = Input_Opt%NYMDe, &
                                    value_NHMSe = Input_Opt%NHMSe, &
                                    value_NYMD  = Input_Opt%NYMDb, &
                                    value_NHMS  = Input_Opt%NHMSb, &
                                    RC          = RC                  )

    !Call Set_NDiagTime( Input_Opt%NHMSe     )
    !Call Set_DiagB    ( Get_Tau()                     )
    !Call Set_XOffset  ( Input_Opt%Nested_I0 )
    !Call Set_YOffset  ( Input_Opt%Nested_J0 )
    !Call Initialize_GEOS_Grid( masterproc, Input_Opt, RC )
    Call Init_Grid   ( am_I_Root = masterproc, &
                       Input_Opt = Input_Opt,  &
                       IM        = nX,         &
                       JM        = nY,         &
                       LM        = nZ,         &
                       RC        = RC          )

    ! Use the GCHP method for setting up the grid - note that this does NOT set
    ! up the areas. In any case we will need to be constantly updating this grid
    ! to compensate for the "multiple chunks per processor" element
    Allocate(lonMidArr(nX,nY), Stat=ierr)
    If (ierr.ne.0) Call endrun('Failure while allocating lonMidArr')
    Allocate(latMidArr(nX,nY), Stat=ierr)
    If (ierr.ne.0) Call endrun('Failure while allocating latMidArr')

    ! We could try and get the data from CAM.. but the goal is to make this GC
    ! component completely grid independent. So for now, we set to arbitrary
    ! values
    dLonFix = 360.0e+0_fp / real(nX,fp)
    dLatFix = 180.0e+0_fp / real(nY,fp)
    if (masterproc) write(iulog,'(a,2(x,I4))') ' nX nY ', nX, nY
    do i=1,nX
       ! Center of box, assuming dateline edge
       lonVal = (-180.0e+0_fp + (0.5e+0_fp * dLonFix) + (real(i-1,fp)*dLonFix)) * pi_180
       do j=1,nY
          ! Center of box, assuming regular cells
          latVal = (-90.0e+0_fp + (0.5e+0_fp * dLatFix) + (real(j-1,fp)*dLatFix)) * pi_180
          lonMidArr(i,j) = lonVal
          latMidArr(i,j) = latVal
       end do
    end do
    Call SetGridFromCtr( masterproc, nX, nY, lonMidArr, latMidArr, RC )
    Deallocate(lonMidArr)
    Deallocate(latMidArr)

    ! Start by setting some dummy timesteps - go with 5 minutes for now
    call gc_update_timesteps(300.0e+0_r8)

    ! Initialize error mod
    call Init_Error( masterproc, Input_Opt, RC )

    ! Initialize the state objects for each chunk. This is like GC_Init_All in
    ! v11-01 of GC-Classic
    do i=begchunk, endchunk
       ! This would be gc_init_stateobj except that that is only in v11-02+ and
       ! also it requires history/diag components which aren't yet dealt with
       rootCPU = (masterproc .and. (i == begchunk))
       Call Init_State_Met( rootCPU, nX, nY, nZ, State_Met(i), RC )
       If (rc.ne.GC_SUCCESS) Call endrun('Could not initialize State_Met')
       Call Init_State_Chm( rootCPU, nX, nY, nZ, &
                            Input_Opt, State_Chm(i), &
                            nDust + nAer, RC )
       If (rc.ne.GC_SUCCESS) Call endrun('Could not initialize State_Chm')

       ! Start with v/v dry (CAM standard)
       State_Chm(i)%Spc_Units = 'v/v dry'
    end do

    ! === START GC_INIT_EXTRA ===
    ! NOTE: Could probably get away with just calling GC_Init_Extra, as long as
    ! the missing modules (eg MODIS_LAI) were sorted out
    If (Input_Opt%LDryD) Then
       Call Init_DryDep( masterproc, Input_Opt, State_Chm(begchunk), RC )
       If (rc.ne.GC_SUCCESS) Call endrun('Failed to initialize wet dry deposition')
    End If

    ! Note: just need *a* copy of State_Chm so that the species DB is accessible
    If (Input_Opt%LConv.or.Input_Opt%LWetD.or.Input_Opt%LChem) Then
       Call Init_WetScav( masterproc, Input_Opt, State_Chm(begchunk), RC)
       If (rc.ne.GC_SUCCESS) Call endrun('Failed to initialize wet scavenging')
    End If

    ! Init_MODIS_LAI..
    !Call Init_MODIS_LAI( masterproc, Input_Opt, RC )
    !If (rc.ne.GC_SUCCESS) Call endrun('Failed to initialize MODIS LAI')

    !Call Set_VDiff_Values( masterproc, Input_Opt, State_Chm(begchunk), RC )
    !If (rc.ne.GC_SUCCESS) Call endrun('Failed to set PBL mixing values')
    !Call Init_Get_NDep( masterproc, RC )
    !If (rc.ne.GC_SUCCESS) Call endrun('Failed to get N deposition data')

    If (Input_Opt%LCarb) Then
       Call Init_Carbon( masterproc, Input_Opt, RC )
       If (rc.ne.GC_SUCCESS) Call endrun('Failed to initialize carbon_mod')
    End If

    If (Input_Opt%LDust) Then
       Call Init_Dust( masterproc, Input_Opt, State_Chm(begchunk), RC )
       If (rc.ne.GC_SUCCESS) Call endrun('Failed to initialize dust_mod')
    End If

    If (Input_Opt%LSSalt) Then
       Call Init_SeaSalt( masterproc, Input_Opt, State_Chm(begchunk), RC )
       If (rc.ne.GC_SUCCESS) Call endrun('Failed to initialize seasalt_mod')
    End If

    If (Input_Opt%LSulf) Then
       Call Init_Sulfate( masterproc, Input_Opt, State_Chm(begchunk), RC )
       If (rc.ne.GC_SUCCESS) Call endrun('Failed to initialize sulfate_mod')
    End If

    If (Input_Opt%LSulf.or.Input_Opt%LCarb.or.Input_Opt%LDust.or.Input_Opt%LSSalt) Then
       Call Init_Aerosol( masterproc, Input_Opt, RC )
       If (rc.ne.GC_SUCCESS) Call endrun('Failed to initialize aerosol_mod')
    End If

    If (Input_Opt%LLinoz) Then
       Call Init_Linoz(masterproc, Input_Opt, RC)
       If (rc.ne.GC_SUCCESS) Call endrun('Failed to initialize linoz_mod')
    End If

    Call Init_ChemGrid( masterproc, Input_Opt, RC )
    If (rc.ne.GC_SUCCESS) Call endrun('Failed to initialize chemgrid')

    If (Input_Opt%LChem) Then
       Call Init_TOMS( masterproc, Input_Opt, RC )
       If (rc.ne.GC_SUCCESS) Call endrun('Failed to initialize chemgrid')
    End If

    Call Init_C2H6( masterproc, Input_Opt, RC )
    If (rc.ne.GC_SUCCESS) Call endrun('Failed to initialize chemgrid')

    ! This is a bare subroutine - no module
    Call NDXX_Setup( masterproc, Input_Opt, State_Chm(begchunk), RC )
    If (rc.ne.GC_SUCCESS) Call endrun('Failed to run NDXX setup')

    ! NOT DONE: TCPORE_BC, Diag_OH, Diag20, GAMAP, Mercury (all),
    !           Tagged simulations (all), Global_CH4

    ! === END GC_INIT_EXTRA ===

    ! Need this because some routines want to know where the PBL top is
    call init_pbl_mix

    ! Set area...
    do i=begchunk,endchunk
       allocate(col_area(ncol(i)),stat=ierr)
       If (ierr.ne.0) Call endrun('Failure while allocating area')
       call get_area_all_p(i, ncol(i), col_area)
       ! Set default value (in case of chunks with fewer columns)
       State_Met(i)%Area_M2 = 1.0e+10_fp
       Do iX=1,nX
       Do iY=1,ncol(i)
          State_Met(i)%Area_M2(iX,iY,:) = real(col_area(iY) * rearth**2,fp)
       End Do
       End Do
       deallocate(col_area)
    end do

    ! Initialize (mostly unused) diagnostic arrays
    ! WARNING: This routine likely calls on modules which are currently excluded
    ! from the GC-CESM build (eg diag03)
    ! Call Initialize( masterproc, Input_Opt, 2, RC )
    ! Call Initialize( masterproc, Input_Opt, 3, RC )
  
    ! Get Ap and Bp from CAM at pressure edges
    Allocate(Ap_CAM_Flip(nZ+1),stat=ierr)
       If (ierr.ne.0) Call endrun('Failure while allocating Ap_CAM_Flip')
    Allocate(Bp_CAM_Flip(nZ+1),stat=ierr)
       If (ierr.ne.0) Call endrun('Failure while allocating Bp_CAM_Flip')
    Ap_CAM_Flip = 0.0e+0_fp
    Bp_CAM_Flip = 0.0e+0_fp
    Do i=1,(nZ+1)
       Ap_CAM_Flip(i) = hyai(nZ+2-i) * ps0 * 100.0e+0_r8
       Bp_CAM_Flip(i) = hybi(nZ+2-i)
    End Do
    Call Init_Pressure( masterproc, Ap_CAM_Flip, Bp_CAM_Flip ) 
    Deallocate(Ap_CAM_Flip,Bp_CAM_Flip)

    If (Input_Opt%Its_A_FullChem_Sim.or.Input_Opt%Its_An_Aerosol_Sim) Then
       ! This also initializes Fast-JX
       Call Init_Chemistry( masterproc, Input_Opt, State_Chm(begchunk), RC )
       If (rc.ne.GC_SUCCESS) Call endrun('Failed to initialize chemistry')
    End If

    ! Initialize HEMCO?
!    CALL EMISSIONS_INIT ( am_I_Root, Input_Opt, State_Met, State_Chm, RC, &
!                          HcoConfig=HcoConfig )
!    ASSERT_(RC==GC_SUCCESS)

    If (Input_Opt%LChem.and.Input_Opt%LUCX) Then
       Call Init_UCX(masterproc, Input_Opt, State_Chm(begchunk)) 
    End If

    ! Init_Pressure...
    ! Init_PBL_Mix...
    ! Init_Chemistry...
    ! Init_TOMS...
    ! Emissions_Init...
    ! Init_UCX...
    ! Convert_Spc_Units...
    ! Can add history output here too with the "addfld" & "add_default" routines
    ! Note that constituents are already output by default

    ! Get the index of H2O
    iH2O = Ind_('H2O')

    ! Get indices for physical fields in physics buffer
    ndx_pblh    = pbuf_get_index('pblh')

    call addfld ( 'BCPI', (/'lev'/), 'A', 'mole/mole', trim('BCPI')//' mixing ratio' )
    call add_default ( 'BCPI',   1, ' ')
    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_INIT'

  end subroutine chem_init

!===============================================================================

  subroutine chem_timestep_init(phys_state, pbuf2d)
    use physics_buffer,   only: physics_buffer_desc

    type(physics_state), intent(in):: phys_state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_TIMESTEP_INIT'

    ! This is when we want to update State_Met and so on
    ! Note that here we have been passed MANY chunks

  end subroutine chem_timestep_init

!===============================================================================

  subroutine gc_update_timesteps(dt)

  use time_mod,       only : set_timesteps

  real(r8), intent(in) :: dt
  real(fp)             :: dt_min_real
  integer              :: dt_min
  integer, save        :: dt_min_last = -1
 
  ! At v11-01, GC still stored timesteps as integer minute counts 
  dt_min_real = real(dt/60.0e+0_r8,fp)

  dt_min = nint(dt_min_real)

  input_opt%ts_chem = dt_min
  input_opt%ts_emis = dt_min
  input_opt%ts_conv = dt_min
  input_opt%ts_dyn  = dt_min
  input_opt%ts_rad  = dt_min

  ! Only bother updating the module information if there's been a change
  if (dt_min.ne.dt_min_last) then
     if (masterproc) write(iulog,'(a,F7.1,a)') ' --> GC: updating dt to ', dt, ' seconds'
     call set_timesteps( masterproc,           &
                         chemistry   = dt_min, &
                         emission    = dt_min, &
                         dynamics    = dt_min, &
                         unit_conv   = dt_min, &
                         convection  = dt_min, &
                         diagnos     = dt_min, &
                         radiation   = dt_min  )
     dt_min_last = dt_min
  end if

  end subroutine

!===============================================================================

  subroutine chem_timestep_tend( state, ptend, cam_in, cam_out, dt, pbuf,  fh2o )

    use physics_buffer,   only: physics_buffer_desc, pbuf_get_field
    use cam_history,      only: outfld
    use camsrfexch,       only: cam_in_t, cam_out_t

    use phys_grid,        only: get_ncols_p, get_rlat_all_p, get_rlon_all_p

    use dao_mod,          only: set_dry_surface_pressure
    use dao_mod,          only: airqnt
    use pressure_mod,     only: set_floating_pressures
    use gc_grid_mod,      only: SetGridFromCtr
    use pressure_mod,     only: accept_external_pedge
    use time_mod,         only: accept_external_date_time
    use strat_chem_mod,   only: init_strat_chem
    use toms_mod,         only: compute_overhead_o3
    use chemistry_mod,    only: do_chemistry
    use wetscav_mod,      only: setup_wetscav

    use gc_grid_mod,      only: Area_M2
    use cmn_size_mod,     only: ptop

    ! For calculating SZA
    use orbit,            only: zenith
    use time_manager,     only: get_curr_calday

    ! Calculating relative humidity
    use wv_saturation,    only: qsat
    use physconst,        only: mwdry

    real(r8),            intent(in)    :: dt          ! time step
    type(physics_state), intent(in)    :: state       ! Physics state variables
    type(physics_ptend), intent(out)   :: ptend       ! indivdual parameterization tendencies
    type(cam_in_t),      intent(inout) :: cam_in
    type(cam_out_t),     intent(in)    :: cam_out
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), optional,  intent(out)   :: fh2o(pcols) ! h2o flux to balance source from chemistry

    ! Initial MMR for all species
    real(r8) :: mmr_beg(pcols,pver,nsls+ntracers)
    real(r8) :: mmr_end(pcols,pver,nsls+ntracers)
    real(r8) :: mmr_tend(pcols,pver,nsls+ntracers)

    ! Mapping (?)
    logical :: lq(pcnst)

    ! Indexing
    integer :: i,j,k,l,n,m
    integer :: nX, nY, nZ
    integer :: lchnk, ncol
    real(r8), dimension(state%ncol) :: &
         csza,      &                                      ! cosine of solar zenith angles
         zsurf, &                                          ! surface height (m)
         rlats, rlons                                      ! chunk latitudes and longitudes (radians)
    real(r8), pointer :: pblh(:)                           ! PBL height on each chunk
    real(r8) :: relhum(state%ncol,pver)                          ! Relative humidity (0-1)
    real(r8) :: satv(  state%ncol,pver)                            ! Work arrays
    real(r8) :: satq(  state%ncol,pver)                            ! Work arrays 
    real(r8) :: qh2o(  state%ncol,pver)                            ! Specific humidity (kg/kg)
    real(r8) :: h2ovmr(state%ncol,pver)                          ! H2O VMR 

    real(fp)     :: lonMidArr(1,pcols), latMidArr(1,pcols)
    integer      :: imaxloc(1)
 
    ! Calculating SZA
    real(r8)     :: calday

    ! Because of strat chem
    logical, save :: schem_ready=.false.

    logical      :: rootChunk
    integer      :: RC

    ! Here's where you'll call DO_CHEMISTRY
    ! NOTE: State_Met etc are in an ARRAY - so we will want to always pass
    ! State_Met%(lchnk) and so on

    ! lchnk: which chunk we have on this process
    lchnk = state%lchnk
    ! ncol: number of atmospheric columns on this chunk
    ncol  = state%ncol
  
    ! Am I the first chunk on the first CPU?
    rootChunk = ( masterproc.and.(lchnk==begchunk) )

    ! Need to update the timesteps throughout the code
    call gc_update_timesteps(dt)

    ! Doesn't actually matter, but for safety's sake...
    ptop = state%pint(1,1)*0.01e+0_fp

    ! Need to be super careful that the module arrays are updated and correctly
    ! set. NOTE: First thing - you'll need to flip all the data vertically

    nX = 1
    nY = ncol
    nZ = pver

    ! Update the grid lat/lons since they are module variables
    ! Assume (!) that area hasn't changed for now, as GEOS-Chem will
    ! retrieve this from State_Met which is chunked
    !call get_rlat_all_p( lchnk, ncol, rlats )
    !call get_rlon_all_p( lchnk, ncol, rlons )
    rlats(1:ncol) = state%lat(1:ncol)
    rlons(1:ncol) = state%lon(1:ncol)

    lonMidArr = 0.0e+0_fp
    latMidArr = 0.0e+0_fp
    do i=1,nX
    do j=1,nY
       lonMidArr(i,j) = rlons(j)
       latMidArr(i,j) = rlats(j)
    end do
    end do
   
    ! Update the grid 
    Call SetGridFromCtr( rootChunk, nX, nY, lonMidArr, latMidArr, RC )
    Area_M2 = State_Met(lchnk)%Area_M2

    ! 2. Copy tracers into State_Chm - again, remember to flip them
    ! Data was received in kg/kg dry
    State_Chm(lchnk)%Spc_Units = 'kg/kg dry'
    lq(:) = .false.
    mmr_beg = 0.0e+0_r8
    do n=1, pcnst
       m = map2gc(n)
       if (m > 0) then
          i=1
          do j=1,ncol
          do k=1,pver
          ! CURRENTLY KG/KG DRY
          mmr_beg(j,k,m) = state%q(j,pver+1-k,n)
          State_Chm(lchnk)%Species(1,j,k,m) = real(mmr_beg(j,k,m),fp)
          end do
          end do
          lq(n) = .true.
       end if
    end do
    call physics_ptend_init(ptend, state%psetcols, 'chemistry', lq=lq)

    ! Calculate cos(SZA)
    calday = get_curr_calday( )
    call zenith( calday, rlats, rlons, csza, ncol )
    !call outfld( 'SZA',   sza,    ncol, lchnk )

    ! Get PBL height (m)
    call pbuf_get_field(pbuf, ndx_pblh, pblh)

    ! Get VMR and MMR of H2O
    h2ovmr = 0.0e0_fp
    qh2o   = 0.0e0_fp
    ! Note MWDRY = 28.966 g/mol
    do j=1,nY
    do l=1,nZ
       qh2o(j,l) = State_Chm(lchnk)%Species(1,j,k,iH2O)
       h2ovmr(j,l) = qh2o(j,l) * mwdry / 18.016e+0_fp 
    end do
    end do

    ! Calculate RH (range 0-1, note still level 1 = TOA)
    relhum(:,:) = 0.0e+0_r8
    call qsat(state%t(:ncol,:), state%pmid(:ncol,:), satv, satq)
    do j=1,nY
    do l=1,nZ
       relhum(j,l) = 0.622e+0_r8 * h2ovmr(j,l) / satq(j,l)
       relhum(j,l) = max( 0.0e+0_r8, min( 1.0e+0_r8, relhum(j,l) ) )
    end do
    end do

    ! << === INCLUDES_BEFORE_RUN === >>
    State_Met(lchnk)%ALBD             (1,:) = cam_in%asdir(:)
    State_Met(lchnk)%CLDFRC           (1,:) = 0.0e+0_fp 
    State_Met(lchnk)%EFLUX            (1,:) = 0.0e+0_fp 
    State_Met(lchnk)%HFLUX            (1,:) = 0.0e+0_fp 
    State_Met(lchnk)%FRCLND           (1,:) = 0.0e+0_fp 
    State_Met(lchnk)%FRLAND           (1,:) = 1.0e+0_fp - cam_in%ocnfrac(:)
    State_Met(lchnk)%FROCEAN          (1,:) = cam_in%ocnfrac(:)
    State_Met(lchnk)%FRLAKE           (1,:) = 0.0e+0_fp 
    State_Met(lchnk)%FRLANDIC         (1,:) = 0.0e+0_fp 
    State_Met(lchnk)%PHIS             (1,:) = 0.0e+0_fp 
    State_Met(lchnk)%GWETROOT         (1,:) = 0.0e+0_fp 
    State_Met(lchnk)%GWETTOP          (1,:) = 0.0e+0_fp 
    State_Met(lchnk)%LAI              (1,:) = 0.0e+0_fp 
    State_Met(lchnk)%PARDR            (1,:) = 0.0e+0_fp 
    State_Met(lchnk)%PARDF            (1,:) = 0.0e+0_fp 
    State_Met(lchnk)%PBLH             (1,:) = pblh(:ncol)
    State_Met(lchnk)%PHIS             (1,:) = state%phis(:)
    State_Met(lchnk)%PRECANV          (1,:) = 0.0e+0_fp 
    State_Met(lchnk)%PRECCON          (1,:) = 0.0e+0_fp 
    State_Met(lchnk)%PRECLSC          (1,:) = 0.0e+0_fp  
    State_Met(lchnk)%PRECTOT          (1,:) = 0.0e+0_fp 
    State_Met(lchnk)%TROPP            (1,:) = 150.0e+0_fp
    State_Met(lchnk)%PS1_WET          (1,:) = state%ps(:)*0.01e+0_fp
    State_Met(lchnk)%PS2_WET          (1,:) = state%ps(:)*0.01e+0_fp
    State_Met(lchnk)%SLP              (1,:) = state%ps(:)*0.01e+0_fp
    State_Met(lchnk)%TS               (1,:) = cam_in%ts(:)
    State_Met(lchnk)%TSKIN            (1,:) = cam_in%ts(:)
    State_Met(lchnk)%SWGDN            (1,:) = 0.0e+0_fp 
    State_Met(lchnk)%TO3              (1,:) = 300.0e+0_fp ! Dummy value
    State_Met(lchnk)%SNODP            (1,:) = 0.0e+0_fp 
    State_Met(lchnk)%SNOMAS           (1,:) = 0.0e+0_fp 
    State_Met(lchnk)%SUNCOS           (1,:) = csza(:)
    State_Met(lchnk)%SUNCOSmid        (1,:) = csza(:)
    State_Met(lchnk)%U10M             (1,:) = state%u(:,nZ) 
    State_Met(lchnk)%USTAR            (1,:) = 0.0e+0_fp 
    State_Met(lchnk)%V10M             (1,:) = state%v(:,nZ)
    State_Met(lchnk)%Z0               (1,:) = 0.0e+0_fp 

    DO J=1,nY
    DO I=1,nX
       IMAXLOC = MAXLOC( (/ State_Met(lchnk)%FRLAND(I,J) +      &
                            State_Met(lchnk)%FRLANDIC(I,J) +    &
                            State_Met(lchnk)%FRLAKE(I,J),       &
                            State_Met(lchnk)%FRSEAICE(I,J),     &
                            State_Met(lchnk)%FROCEAN(I,J) -     &
                            State_Met(lchnk)%FRSEAICE(I,J) /) )
       IF ( IMAXLOC(1) == 3 ) IMAXLOC(1) = 0 ! reset ocean to 0
       State_Met(lchnk)%LWI(I,J) = FLOAT( IMAXLOC(1) )
    ENDDO
    ENDDO

    ! Three-dimensional fields on level edges
    do j=1,nY
    do l=1,nZ+1 
       State_Met(lchnk)%CMFMC   (1,j,l)        = 0.0e+0_fp 
       State_Met(lchnk)%PFICU   (1,j,l)        = 0.0e+0_fp 
       State_Met(lchnk)%PFILSAN (1,j,l)        = 0.0e+0_fp 
       State_Met(lchnk)%PFLCU   (1,j,l)        = 0.0e+0_fp 
       State_Met(lchnk)%PFLLSAN (1,j,l)        = 0.0e+0_fp 
       State_Met(lchnk)%PEDGE   (1,j,l)        = state%pint(j,nZ+2-l)*0.01e+0_fp
    end do
    end do

    ! These are set later
    State_Met(lchnk)%PS1_DRY (:,:) = 0.0e+0_fp
    State_Met(lchnk)%PS2_DRY (:,:) = 0.0e+0_fp

    ! Calculate CLDTOPS (highest location of CMFMC in the column)
    i = 1
    Do j=1,nY
       State_Met(lchnk)%CldTops(I,J) = 1
       Do l = nZ, 1, -1
          If ( State_Met(lchnk)%CMFMC(i,j,l) > 0.0e+0_fp ) Then
             State_Met(lchnk)%CldTops(i,j) = l + 1
             Exit
          End If
       End Do
    End Do

    ! Three-dimensional fields on level centers
    do j=1,nY
    do l=1,nZ 
       State_Met(lchnk)%U       (1,j,l)        = state%u(j,nZ+1-l)
       State_Met(lchnk)%V       (1,j,l)        = state%v(j,nZ+1-l)
       !State_Met(lchnk)%OMEGA   (1,j,l)        = state%omega(j,nZ+1-l)
       State_Met(lchnk)%CLDF    (1,j,l)        = 0.0e+0_fp 
       State_Met(lchnk)%DTRAIN  (1,j,l)        = 0.0e+0_fp
       State_Met(lchnk)%DQRCU   (1,j,l)        = 0.0e+0_fp 
       State_Met(lchnk)%DQRLSAN (1,j,l)        = 0.0e+0_fp 
       State_Met(lchnk)%QI      (1,j,l)        = 0.0e+0_fp 
       State_Met(lchnk)%QL      (1,j,l)        = 0.0e+0_fp 
       State_Met(lchnk)%RH      (1,j,l)        = relhum(j,nZ+1-l) * 100.0e+0_fp 
       State_Met(lchnk)%TAUCLI  (1,j,l)        = 0.0e+0_fp 
       State_Met(lchnk)%TAUCLW  (1,j,l)        = 0.0e+0_fp 
       State_Met(lchnk)%REEVAPCN(1,j,l)        = 0.0e+0_fp 
       State_Met(lchnk)%REEVAPLS(1,j,l)        = 0.0e+0_fp 
       State_Met(lchnk)%SPHU1   (1,j,l)        = qh2o(j,nZ+1-l) * 1.0e+3_fp ! g/kg
       State_Met(lchnk)%SPHU2   (1,j,l)        = qh2o(j,nZ+1-l) * 1.0e+3_fp ! g/kg
       State_Met(lchnk)%TMPU1   (1,j,l)        = state%t(j,nZ+1-l)
       State_Met(lchnk)%TMPU2   (1,j,l)        = state%t(j,nZ+1-l)
    end do
    end do

    ! Derived fields
    State_Met(lchnk)%T    = (State_Met(lchnk)%TMPU1  + State_Met(lchnk)%TMPU2)*0.5e+0_fp
    State_Met(lchnk)%SPHU = (State_Met(lchnk)%SPHU1  + State_Met(lchnk)%SPHU2)*0.5e+0_fp
  
    ! Calculate total OD as liquid cloud OD + ice cloud OD
    State_Met(lchnk)%OPTD =  State_Met(lchnk)%TAUCLI + State_Met(lchnk)%TAUCLW

    ! Nullify all pointers
    Nullify(pblh)
    ! << === INCLUDES_BEFORE_RUN === >> 

    ! Eventually initialize/reset wetdep
    IF ( Input_Opt%LConv .OR. Input_Opt%LChem .OR. Input_Opt%LWetD ) THEN
       CALL SETUP_WETSCAV( rootChunk, Input_Opt, State_Met(lchnk), State_Chm(lchnk), RC )
       If (rc.ne.GC_SUCCESS) Call endrun('Failed to set up wet scavenging')
    ENDIF

    ! Pass time values obtained from the ESMF environment to GEOS-Chem
    CALL Accept_External_Date_Time( am_I_Root      = rootChunk,  &
                                    value_NYMD     = 20000101,   &  
                                    value_NHMS     = 0000,       &  
                                    value_YEAR     = 2000,       &  
                                    value_MONTH    = 01,         &  
                                    value_DAY      = 01,         &  
                                    value_DAYOFYR  = 001,        &  
                                    value_HOUR     = 00,         &  
                                    value_MINUTE   = 00,         &  
                                    value_HELAPSED = 0.0e+0_f4,  & 
                                    value_UTC      = 0.0e+0_f4,  &
                                    RC             = RC         )
    If (rc.ne.GC_SUCCESS) Call endrun('Failed to update time in GEOS-Chem')

    Call Accept_External_PEdge( am_I_Root = rootChunk,        &
                                State_Met = State_Met(lchnk), &
                                RC        = RC                )
    If (rc.ne.GC_SUCCESS) Call endrun('Failed to update pressure edges')

    ! Calculate State_Met etc for this timestep
    ! Use the CAM psdry fields instead of using the GC calculation
    !Call Set_Dry_Surface_Pressure(State_Met(lchnk), 1)
    State_Met(lchnk)%PS1_DRY (1,:) = state%psdry(:)*0.01e+0_fp
    State_Met(lchnk)%PS2_DRY (1,:) = state%psdry(:)*0.01e+0_fp

    ! Set surface pressures to match those in input
    State_Met(lchnk)%PSC2_WET = State_Met(lchnk)%PS1_WET
    State_Met(lchnk)%PSC2_DRY = State_Met(lchnk)%PS1_DRY
    Call Set_Floating_Pressures( rootChunk, State_Met(lchnk), RC )
    If (rc.ne.GC_SUCCESS) Call endrun('Failed to set floating pressures')

    ! Set quantities of interest but do not change VMRs
    Call AirQnt( rootChunk, Input_Opt, State_Met(lchnk), &
                 State_Chm(lchnk), RC, update_mixing_ratio=.False. )
    If (rc.ne.GC_SUCCESS) Call endrun('Failed to calculate air properties')

    ! Initialize strat chem if not already done. This has to be done here because
    ! it needs to have non-zero values in State_Chm%AD, which only happens after
    ! the first call to AirQnt
    IF ( (.not.schem_ready).and. Input_Opt%LSCHEM ) THEN
       CALL INIT_STRAT_CHEM( rootChunk, Input_Opt, State_Chm(lchnk), State_Met(lchnk), RC )
       If (rc.ne.GC_SUCCESS) Call endrun('Could not initialize strat-chem')
       schem_ready = .True.
    ENDIF

    ! Run chemistry
    If (Input_Opt%LChem) Then
       Call Compute_Overhead_O3( rootChunk, 1, &
                Input_Opt%Use_O3_From_Met, State_Met(lchnk)%TO3 )
       Call Do_Chemistry( am_I_Root      = rootChunk,        &
                          Input_Opt      = Input_Opt,        &
                          State_Chm      = State_Chm(lchnk), &
                          State_Met      = State_Met(lchnk), &
                          RC             = RC                )
    End If

    !if (masterproc) write(iulog,*) ' --> TEND SIZE: ', size(state%ncol)
    !if (masterproc) write(iulog,'(a,2(x,I6))') ' --> TEND SIDE:  ', lbound(state%ncol),ubound(state%ncol)
    if (rootChunk) write(iulog,'(a)') 'GCCALL CHEM_TIMESTEP_TEND'

    ! Make sure State_Chm(lchnk) is back in kg/kg dry!

    ! Reset H2O MMR to the initial value (no chemistry tendency in H2O just yet)
    State_Chm(lchnk)%Species(1,:,:,iH2O) = mmr_beg(:,:,iH2O)

    ! NOTE: Re-flip all the arrays vertically or suffer the consequences
    ! ptend%q dimensions: [column, ?, species]
    ptend%q(:,:,:) = 0.0e+0_r8 
    mmr_end = 0.0e+0_r8
    do n=1, pcnst
       m = map2gc(n)
       if (m > 0) then
          i=1
          do j=1,ncol
          do k=1,pver
             ! CURRENTLY KG/KG
             mmr_end( j,k,m) = real(State_Chm(lchnk)%Species(1,j,k,m),r8)
             mmr_tend(j,k,m) = mmr_end(j,k,m) - mmr_beg(j,k,m)
             ptend%q(j,pver+1-k,n) = (mmr_end(j,k,m)-mmr_beg(j,k,m))/dt
          end do
          end do
       end if
    end do
    if (present(fh2o)) fh2o(:) = 0.0e+0_r8

    return
  end subroutine chem_timestep_tend

!===============================================================================
  subroutine chem_init_cnst(name, latvals, lonvals, mask, q)

    character(len=*), intent(in)  :: name       !  constituent name
    real(r8),         intent(in)  :: latvals(:) ! lat in degrees (ncol)
    real(r8),         intent(in)  :: lonvals(:) ! lon in degrees (ncol)
    logical,          intent(in)  :: mask(:)    ! Only initialize where .true.
    real(r8),         intent(out) :: q(:,:)     ! kg tracer/kg dry air (ncol, pver)
    ! Used to initialize tracer fields if desired.
    ! Will need a simple mapping structure as well as the CAM tracer registration
    ! routines.

    integer  :: ilev, nlev
    real(r8) :: qtemp

    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_INIT_CNST'

    nlev = size(q, 2)
    if ( any( tracernames .eq. name ) ) then
       ! Retrieve a "background value" for this from the database
       do ilev=1,nlev
          where(mask)
             ! Set to the minimum mixing ratio
             q(:,ilev) = 1.0e-38_r8
          end where
       end do
    end if

  end subroutine chem_init_cnst

!===============================================================================
  subroutine chem_final
   
    use input_opt_mod,   only : cleanup_input_opt
    use state_chm_mod,   only : cleanup_state_chm
    use state_met_mod,   only : cleanup_state_met
    use error_mod,       only : cleanup_error

    use ucx_mod,         only : cleanup_ucx
    use linoz_mod,       only : cleanup_linoz
    use drydep_mod,      only : cleanup_drydep
    use wetscav_mod,     only : cleanup_wetscav
    !use MODIS_LAI_Mod,   Only : cleanup_MODIS_LAI
    use carbon_mod,      only : cleanup_carbon
    use dust_mod,        only : cleanup_dust
    use seasalt_mod,     only : cleanup_seasalt
    use aerosol_mod,     only : cleanup_aerosol
    use chemgrid_mod,    only : cleanup_chemgrid
    use TOMS_mod,        only : cleanup_TOMS
    use C2H6_mod,        only : cleanup_C2H6
    use sulfate_mod,     only : cleanup_sulfate
    use pressure_mod,    only : cleanup_pressure
    use flexchem_mod,    only : cleanup_flexchem
    use strat_chem_mod,  only : cleanup_strat_chem

    use cmn_size_mod,    only : cleanup_cmn_size
    use cmn_o3_mod,      only : cleanup_cmn_o3
    use cmn_fjx_mod,     only : cleanup_cmn_fjx
    use ucx_mod,         only : cleanup_ucx
    use pbl_mix_mod,     only : cleanup_pbl_mix

    ! Special: cleans up after NDXX_Setup
    use Diag_mod,      only : cleanup_diag
 
    integer :: i, rc
    logical :: rootCPU

    ! Finalize GEOS-Chem
    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_FINAL'
    ! Clean up module variables
    Call Cleanup_Linoz
    Call Cleanup_Drydep
    Call Cleanup_Wetscav
    Call Cleanup_Carbon
    Call Cleanup_Dust
    !Call Cleanup_MODIS_LAI
    Call Cleanup_SeaSalt
    Call Cleanup_ChemGrid
    Call Cleanup_TOMS
    Call Cleanup_Sulfate
    Call Cleanup_Aerosol
    Call Cleanup_Diag
    Call Cleanup_C2H6
    Call Cleanup_FlexChem
    Call Cleanup_Strat_Chem
    Call Cleanup_UCX( masterproc )
    Call Cleanup_PBL_Mix
    ! Loop over each chunk and clean up the state variables
    Do i=begchunk,endchunk
       rootCPU = ((i.eq.begchunk) .and. MasterProc)
       Call Cleanup_State_Met( rootCPU, State_Met(i), RC )
       Call Cleanup_State_Chm( rootCPU, State_Chm(i), RC )
    End Do
    ! Clean up error modules and lower-level modules
    Call Cleanup_Error
    Call Cleanup_CMN_SIZE( masterproc, RC )
    Call Cleanup_CMN_FJX(  masterproc, RC )
    Call Cleanup_CMN_O3(   masterproc, RC )
    ! Clean up input_opt
    Call Cleanup_Input_Opt( masterproc, Input_Opt, RC )
    ! Finally deallocate the variables in full
    If (allocated(State_Met))     Deallocate(State_Met)
    If (allocated(State_Chm))     Deallocate(State_Chm)
    if (masterproc) write(iulog,'(a,2(x,L1))') ' --> DEALLOC CHECK : ', Allocated(state_met), Allocated(state_chm)

    return
  end subroutine chem_final
!===============================================================================
  subroutine chem_init_restart(File)
    use pio, only : file_desc_t
    type(file_desc_t) :: File
    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_INIT_RESTART'
    return
  end subroutine chem_init_restart
!===============================================================================
  subroutine chem_write_restart( File )
    !use tracer_cnst, only: write_tracer_cnst_restart
    !use tracer_srcs, only: write_tracer_srcs_restart
    !use linoz_data,  only: write_linoz_data_restart
    use pio, only : file_desc_t
    implicit none
    type(file_desc_t) :: File

    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_WRITE_RESTART'
    !
    ! data for offline tracers
    !
    !call write_tracer_cnst_restart(File)
    !call write_tracer_srcs_restart(File)
    !call write_linoz_data_restart(File)
  end subroutine chem_write_restart
!===============================================================================
  subroutine chem_read_restart( File )
    !use tracer_cnst, only: read_tracer_cnst_restart
    !use tracer_srcs, only: read_tracer_srcs_restart
    !use linoz_data,  only: read_linoz_data_restart

    use pio, only : file_desc_t
    implicit none
    type(file_desc_t) :: File

    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_READ_RESTART'
    !
    ! data for offline tracers
    !
    !call read_tracer_cnst_restart(File)
    !call read_tracer_srcs_restart(File)
    !call read_linoz_data_restart(File)
  end subroutine chem_read_restart
!================================================================================
  subroutine chem_emissions( state, cam_in )
    use camsrfexch,       only: cam_in_t     

    ! Arguments:

    type(physics_state),    intent(in)    :: state   ! Physics state variables
    type(cam_in_t),         intent(inout) :: cam_in  ! import state

    integer :: lchnk, ncol
    logical :: rootChunk

    ! lchnk: which chunk we have on this process
    lchnk = state%lchnk
    ! ncol: number of atmospheric columns on this chunk
    ncol  = state%ncol
    rootChunk = (masterproc.and.(lchnk.eq.begchunk))

    if (rootChunk) write(iulog,'(a)') 'GCCALL CHEM_EMISSIONS'

  end subroutine chem_emissions

end module chemistry
