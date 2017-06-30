module dyn_comp
!----------------------------------------------------------------------- 
! 
! SLD dycore interface module
!
!-----------------------------------------------------------------------

use shr_kind_mod,    only: r8=>shr_kind_r8

use spmd_utils,      only: masterproc, npes, mpicom, mpir8

use physconst,       only: rair
use pmgrid,          only: plon, plat, plev, plevp, plnlv, beglat, endlat
use dyn_grid,        only: ptimelevels

use prognostics,     only: ps, u3, v3, t3, tm, tl, q3, qm, ql, qm1, div, &
                           phis, phism, phisl, dpsm, dpsl, ed1

use cam_control_mod, only: initial_run, ideal_phys, aqua_planet, moist_physics, adiabatic
use constituents,    only: pcnst, cnst_name, cnst_longname, cnst_read_iv, qmin, &
                           tendnam, fixcnam, tottnam, hadvnam, vadvnam
use cam_initfiles,   only: initial_file_get_id, topo_file_get_id, pertlim
use cam_history,     only: addfld, add_default, horiz_only

use sld_control_mod, only: dif2, dif4, divdampn, eps, kmxhdc

use cam_pio_utils,   only: clean_iodesc_list
use pio,             only: file_desc_t, pio_noerr, pio_inq_attlen, pio_get_att, pio_inq_varid, &
                           pio_seterrorhandling, pio_bcast_error, pio_internal_error

#if (defined SPMD)
use spmd_dyn,        only: spmd_readnl
#endif

use cam_abortutils,  only: endrun
use cam_logfile,     only: iulog

implicit none
private
save

public :: &
   dyn_import_t, &
   dyn_export_t, &
   dyn_readnl,   &
   dyn_register, &
   dyn_init

! these structures are not used in this dycore, but are included
! for source code compatibility.  
type dyn_import_t
   integer :: placeholder
end type dyn_import_t

type dyn_export_t
   integer :: placeholder
end type dyn_export_t

real(r8), allocatable :: ps_tmp  (:,:  )
real(r8), allocatable :: phis_tmp(:,:  )
real(r8), allocatable :: q3_tmp  (:,:,:)
real(r8), allocatable :: t3_tmp  (:,:,:)
real(r8), allocatable :: arr3d_a (:,:,:)
real(r8), allocatable :: arr3d_b (:,:,:)

logical readvar            ! inquiry flag:  true => variable exists on netCDF file

!=============================================================================================
contains
!=============================================================================================

subroutine dyn_readnl(nlfile)

   ! Read SLD namelist group.

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_integer, mpi_real8

   ! args
   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
    
   ! local vars
   integer :: unitn, ierr

   real(r8) :: sld_dif2_coef          ! del2 horizontal diffusion coeff.
   real(r8) :: sld_dif4_coef          ! del4 horizontal diffusion coeff.
   real(r8) :: sld_divdampn      ! Number of days to invoke divergence damper
   real(r8) :: sld_tfilt_eps           ! time filter coefficient. Defaults to 0.06.
   integer  :: sld_kmxhdc        ! Number of levels to apply Courant limiter
    
   namelist /dyn_sld_inparm/ sld_dif2_coef, sld_dif4_coef, &
      sld_divdampn, sld_tfilt_eps, sld_kmxhdc

   character(len=*), parameter :: sub = 'dyn_readnl'
   !--------------------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'dyn_sld_inparm', status=ierr)
      if (ierr == 0) then
         read(unitn, dyn_sld_inparm, iostat=ierr)
         if (ierr /= 0) then
            call endrun(sub//': ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   call mpi_bcast(sld_dif2_coef, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: sld_dif2_coef")

   call mpi_bcast(sld_dif4_coef, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: sld_dif4_coef")

   call mpi_bcast(sld_divdampn, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: sld_divdampn")

   call mpi_bcast(sld_tfilt_eps, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: sld_tfilt_eps")

   call mpi_bcast(sld_kmxhdc, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: sld_kmxhdc")

   dif2     = sld_dif2_coef
   dif4     = sld_dif4_coef
   divdampn = sld_divdampn
   eps      = sld_tfilt_eps
   kmxhdc   = sld_kmxhdc

   ! Write namelist variables to logfile
   if (masterproc) then

      write(iulog,*) 'SLD Dycore Parameters:'

      if (divdampn > 0._r8) then
         write(iulog,*) '  Divergence damper for spectral dycore invoked for days 0. to ',divdampn,' of this case'
      elseif (divdampn < 0._r8) then
         call endrun ('READ_NAMELIST: divdampn must be a positive number')
      else
         write(iulog,*) '  Divergence damper for spectral dycore NOT invoked'
      endif

      if (kmxhdc >= plev .or. kmxhdc < 0) then
         call endrun ('DYN_SLD_READNL:  ERROR:  KMXHDC must be between 0 and plev-1')
      end if

      write(iulog,9108) eps, dif2, dif4, kmxhdc
   end if

#if (defined SPMD)
   call spmd_readnl(nlfile)
#endif 

9108 format('   Time filter coefficient (EPS)                 ',f10.3,/,&
            '   DEL2 Horizontal diffusion coefficient (DIF2)  ',e10.3/, &
            '   DEL4 Horizontal diffusion coefficient (DIF4)  ',e10.3/, &
            '   Number of levels Courant limiter applied      ',i10)

end subroutine dyn_readnl

!=========================================================================================

subroutine dyn_register()
end subroutine dyn_register

!=========================================================================================

subroutine dyn_init(dyn_in, dyn_out)

   use prognostics,         only: initialize_prognostics
   use scanslt,             only: slt_alloc
   use phys_control,        only: phys_getopts

#if (defined SPMD)
   use spmd_dyn,            only: spmdbuf
#endif

   ! Arguments are not used in this dycore, included for compatibility
   type(dyn_import_t), intent(out) :: dyn_in
   type(dyn_export_t), intent(out) :: dyn_out

   ! Local variables

   logical :: history_amwg       ! output for AMWG diagnostics
   integer m
   !----------------------------------------------------------------------------

   ! Initialize prognostics variables
   call initialize_prognostics
   call slt_alloc()

#if (defined SPMD)
   ! Allocate communication buffers for collective communications in realloc
   ! routines and in dp_coupling
   call spmdbuf ()
#endif

   if (initial_run) then
      call read_inidat()
      call clean_iodesc_list()
   end if

   call addfld ('ETADOT',(/ 'ilev' /),'A','1/s','Vertical (eta) velocity',gridname='gauss_grid')
   call addfld ('U&IC',  (/ 'lev' /), 'I','m/s','Zonal wind',             gridname='gauss_grid' )
   call addfld ('V&IC',  (/ 'lev' /), 'I','m/s','Meridional wind',        gridname='gauss_grid' )
   call add_default ('U&IC',0, 'I')
   call add_default ('V&IC',0, 'I')

   call addfld ('PS&IC',horiz_only, 'I','Pa','Surface pressure',          gridname='gauss_grid' )
   call addfld ('T&IC', (/ 'lev' /),'I','K', 'Temperature',               gridname='gauss_grid' )
   call add_default ('PS&IC',0, 'I')
   call add_default ('T&IC',0, 'I')

   do m = 1, pcnst
      call addfld (trim(cnst_name(m))//'&IC',(/ 'lev' /),'I', 'kg/kg',cnst_longname(m),                      gridname='gauss_grid')
      call add_default(trim(cnst_name(m))//'&IC',0, 'I')
      call addfld (hadvnam(m),(/ 'lev' /),'A','kg/kg/s',trim(cnst_name(m))//' horizontal advection tendency',gridname='gauss_grid')
      call addfld (vadvnam(m),(/ 'lev' /),'A','kg/kg/s',trim(cnst_name(m))//' vertical advection tendency',  gridname='gauss_grid')
      call addfld (tendnam(m),(/ 'lev' /),'A','kg/kg/s',trim(cnst_name(m))//' total tendency',               gridname='gauss_grid')
      call addfld (tottnam(m),(/ 'lev' /),'A','kg/kg/s',trim(cnst_name(m))//' horz + vert + fixer tendency', gridname='gauss_grid')
      call addfld (fixcnam(m),(/ 'lev' /),'A','kg/kg/s',trim(cnst_name(m))//' tendency due to slt fixer',    gridname='gauss_grid')
   end do

   call addfld ('DUH     ',(/ 'lev' /),'A','K/s     ','U horizontal diffusive heating',              gridname='gauss_grid')
   call addfld ('DVH     ',(/ 'lev' /),'A','K/s     ','V horizontal diffusive heating',              gridname='gauss_grid')
   call addfld ('DTH     ',(/ 'lev' /),'A','K/s     ','T horizontal diffusive heating',              gridname='gauss_grid')
   call addfld ('ENGYCORR',(/ 'lev' /),'A','W/m2    ','Energy correction for over-all conservation', gridname='gauss_grid')
   call addfld ('TFIX    ',horiz_only, 'A','K/s     ','T fixer (T equivalent of Energy correction)', gridname='gauss_grid')

   call phys_getopts(history_amwg_out=history_amwg)

   if (history_amwg) call add_default ('DTH     ', 1, ' ')

   call addfld ('FU      ',(/ 'lev' /),'A','m/s2    ','Zonal wind forcing term',            gridname='gauss_grid')
   call addfld ('FV      ',(/ 'lev' /),'A','m/s2    ','Meridional wind forcing term',       gridname='gauss_grid')
   call addfld ('UTEND   ',(/ 'lev' /),'A','m/s2    ','U tendency',                         gridname='gauss_grid')
   call addfld ('VTEND   ',(/ 'lev' /),'A','m/s2    ','V tendency',                         gridname='gauss_grid')
   call addfld ('TTEND   ',(/ 'lev' /),'A','K/s     ','T tendency',                         gridname='gauss_grid')
   call addfld ('LPSTEN  ',horiz_only, 'A','Pa/s    ','Surface pressure tendency',          gridname='gauss_grid')
   call addfld ('VAT     ',(/ 'lev' /),'A','K/s     ','Vertical advective tendency of T',   gridname='gauss_grid')
   call addfld ('KTOOP   ',(/ 'lev' /),'A','K/s     ','(Kappa*T)*(omega/P)',                gridname='gauss_grid')

end subroutine dyn_init

!=========================================================================================
! Private routines
!=========================================================================================

subroutine read_inidat()

   ! Read initial dataset and spectrally truncate as appropriate.
   ! Read and process the fields one at a time to minimize 
   ! memory usage.

   use cam_pio_utils,    only: cam_pio_get_var
   use inic_analytic,    only: analytic_ic_active, analytic_ic_set_ic
   use dyn_tests_utils,  only: vc_moist_pressure

   ! Local variables

   integer i,c,m,n                         ! indices
   integer ncol

   type(file_desc_t), pointer :: fh_ini, fh_topo

   real(r8), pointer, dimension(:,:,:)   :: convptr_2d
   real(r8), pointer, dimension(:,:,:,:) :: convptr_3d
   real(r8), pointer, dimension(:,:,:,:) :: cldptr
   real(r8), pointer, dimension(:,:    ) :: arr2d_tmp
   real(r8), pointer, dimension(:,:    ) :: arr2d
   character(len=3), parameter :: arraydims3(3) = (/ 'lon', 'lev', 'lat' /)
   character(len=3), parameter :: arraydims2(2) = (/ 'lon', 'lat' /)
   character(len=16) fieldname                  ! field name
   ! variables for analytic initial conditions
   integer,  allocatable       :: glob_ind(:)
   integer                     :: m_cnst(1)
   real(r8), allocatable       :: q4_tmp(:,:,:,:)

   character(len=*), parameter :: sub='read_inidat'
   !----------------------------------------------------------------------------

   fh_ini  => initial_file_get_id()
   fh_topo => topo_file_get_id()

   allocate ( ps_tmp  (plon,plat     ) )
   allocate ( phis_tmp(plon,plat     ) )
   allocate ( q3_tmp  (plon,plev,plat) )
   allocate ( t3_tmp  (plon,plev,plat) )
   allocate(  arr3d_a(plon,plev,plat)  )
   allocate(  arr3d_b(plon,plev,plat)  )

   if (analytic_ic_active()) then
      allocate(glob_ind(plon * plat))
      m = 1
      do c = 1, plat
         do i = 1, plon
            ! Create a global column index
            glob_ind(m) = i + (c-1)*plon
            m = m + 1
         end do
      end do
      call analytic_ic_set_ic(vc_moist_pressure, clat(:), clon(:,1),          &
           glob_ind(:), U=arr3d_a, V=arr3d_b, T=t3_tmp, PS=ps_tmp, PHIS=phis_tmp)
      readvar = .false.
      call process_inidat('PS')
      call process_inidat('UV')
      call process_inidat('T')
      ! Only use analytic value (zero) for PHIS if topography is not specified
      if (.not. associated(fh_topo)) then
        call process_inidat('PHIS')
      end if
      allocate(q4_tmp(plon,plev,plat,1))
      do m = 1, pcnst
         m_cnst(1) = m
         call analytic_ic_set_ic(vc_moist_pressure, clat(:), clon(:,1),        &
              glob_ind(:), Q=q4_tmp, m_cnst=m_cnst)
         arr3d_a(:,:,:) = q4_tmp(:,:,:,1)
         call process_inidat('CONSTS', m_cnst=m, fh=fh_ini)
      end do
      deallocate(q4_tmp)
      deallocate(glob_ind)
      deallocate(arr3d_a)
      deallocate(arr3d_b)
   else
      !---------------------
      ! Read required fields
      !---------------------

      !-----------
      ! 3-D fields
      !-----------

      fieldname = 'U'
      call cam_pio_get_var(fieldname, fh_ini, arraydims3, arr3d_a, found=readvar)
      if (.not. readvar) then
         call endrun(sub//': ERROR: reading '//trim(fieldname))
      end if

      fieldname = 'V'
      call cam_pio_get_var(fieldname, fh_ini, arraydims3, arr3d_b, found=readvar)
      if (.not. readvar) then
         call endrun(sub//': ERROR: reading '//trim(fieldname))
      end if

      call process_inidat('UV')


      fieldname = 'T'
      call cam_pio_get_var(fieldname, fh_ini, arraydims3, t3_tmp, found=readvar)
      if (.not. readvar) then
         call endrun(sub//': ERROR: reading '//trim(fieldname))
      end if
      call process_inidat('T')

      ! Constituents (read and process one at a time)

      do m = 1, pcnst

         readvar   = .false.
         fieldname = cnst_name(m)
         if (cnst_read_iv(m)) &
              call cam_pio_get_var(fieldname, fh_ini, arraydims3, arr3d_a, found=readvar)
         call process_inidat('CONSTS', m_cnst=m, fh=fh_ini)

      end do

      deallocate ( arr3d_a  )
      deallocate ( arr3d_b  )

      !-----------
      ! 2-D fields
      !-----------

      fieldname = 'PS'
      call cam_pio_get_var(fieldname, fh_ini, arraydims2, ps_tmp, found=readvar)
      if (.not. readvar) then
         call endrun(sub//': ERROR: reading '//trim(fieldname))
      end if
      call process_inidat('PS')
   end if

   ! PHIS processing for cases without analytic initial conditions
   fieldname = 'PHIS'
   readvar   = .false.
   if (associated(fh_topo)) then
      call cam_pio_get_var(fieldname, fh_topo, arraydims2, phis_tmp, found=readvar)
      if (.not. readvar) then
         call endrun(sub//': ERROR: reading '//trim(fieldname))
      end if
      call process_inidat('PHIS', fh=fh_topo)
   else if (.not. analytic_ic_active()) then
      phis_tmp(:,:) = 0._r8
      call process_inidat('PHIS')
   end if

   ! Integrals of mass, moisture and geopotential height
   ! (fix mass of moisture as well)
   call global_int

   deallocate ( ps_tmp   )
   deallocate ( phis_tmp )

   ! Initialization of other misc. required fields (this could be moved elsewhere)
   ed1(:,:,:) = 0._r8

   deallocate ( q3_tmp  )
   deallocate ( t3_tmp  )

   call copytimelevels()

end subroutine read_inidat

!=========================================================================================

subroutine process_inidat(fieldname, m_cnst, fh)

   ! Post-process input fields

   use pspect,       only: psp
   use commap
   use comspe
   use spetru,       only: spetru_phis, spetru_ps, spetru_3d_scalar, spetru_uv
   use const_init,          only: cnst_init_default

#if ( defined SPMD )
   use spmd_dyn, only: compute_gsfactors
#endif

   ! arguments
   character(len=*),  intent(in)              :: fieldname ! fields to be processed
   integer,           intent(in),    optional :: m_cnst    ! constituent index
   type(file_desc_t), intent(inout), optional :: fh

   !---------------------------Local workspace-----------------------------

   integer i,ii,j,k,n,lat,irow            ! grid and constituent indices
   real(r8) pertval                       ! perturbation value
   real(r8) tmp1                          ! tmp space
   real(r8) phi(2,psp/2)                  ! used in spectral truncation of phis
                                          ! using "B" part of hybrid grid only
   integer  varid                         ! netCDF variable id
   integer  ret, attlen                   ! netcdf return values
   logical  phis_hires                    ! true => PHIS came from hi res topo
   character*256 text
   character*256 trunits                  ! tracer untis

   real(r8), pointer, dimension(:,:,:) :: q_tmp
   real(r8), pointer, dimension(:,:,:) :: tmp3d_a, tmp3d_b, tmp3d_c, tmp3d_extend
   real(r8), pointer, dimension(:,:  ) :: tmp2d_a, tmp2d_b

#if ( defined SPMD )
   integer :: numperlat                   ! number of values per latitude band
   integer :: numsend(0:npes-1)           ! number of items to be sent
   integer :: numrecv                     ! number of items to be received
   integer :: displs(0:npes-1)            ! displacement array
#endif

   character(len=*), parameter :: sub='process_inidat'
   !---------------------------------------------------------------------------

   select case (fieldname)

   !------------
   ! Process U/V
   !------------

   case ('UV')

      allocate(tmp3d_a(plon,plev,plat))

      ! Spectral truncation

      call spetru_uv(arr3d_a, arr3d_b, div=tmp3d_a)

#if ( defined SPMD )
      numperlat = plnlv
      call compute_gsfactors(numperlat, numrecv, numsend, displs)
      call mpiscatterv(arr3d_a, numsend, displs, mpir8, u3 (:,:,beglat:endlat,1), numrecv, mpir8,0,mpicom)
      call mpiscatterv(arr3d_b, numsend, displs, mpir8, v3 (:,:,beglat:endlat,1), numrecv, mpir8,0,mpicom)
      call mpiscatterv(tmp3d_a, numsend, displs, mpir8, div(:,:,beglat:endlat,1), numrecv, mpir8,0,mpicom)
#else
      u3    (:,:,:,1) = arr3d_a(:plon,:plev,:plat)
      v3    (:,:,:,1) = arr3d_b(:plon,:plev,:plat)
      div   (:,:,:,1) = tmp3d_a(:,:,:)
#endif

      deallocate(tmp3d_a)

   !----------
   ! Process T
   !----------

   case ('T')

      allocate ( tmp3d_a(plon,plev,plat) )
      allocate ( tmp3d_b(plon,plev,plat) )

      ! Add random perturbation to temperature if required

      if (pertlim .ne. 0.0_r8) then
         if (masterproc) write(iulog,*) sub//': INFO:  Adding random perturbation bounded by +/-', &
            pertlim,' to initial temperature field'
         do lat = 1, plat
            do k = 1, plev
               do i = 1, plon
                  call random_number (pertval)
                  pertval = 2._r8*pertlim*(0.5_r8 - pertval)
                  t3_tmp(i,k,lat) = t3_tmp(i,k,lat)*(1._r8 + pertval)
               end do
            end do
         end do
      end if

      ! Spectral truncation

      call spetru_3d_scalar(t3_tmp, dl=tmp3d_a, dm=tmp3d_b)

#if ( defined SPMD )
      numperlat = plnlv
      call compute_gsfactors(numperlat, numrecv, numsend, displs)
      call mpiscatterv(t3_tmp,  numsend, displs, mpir8, t3(:,:,beglat:endlat,1) ,numrecv, mpir8,0,mpicom)
      call mpiscatterv(tmp3d_a, numsend, displs, mpir8, tl(:,:,beglat:endlat) ,numrecv, mpir8,0,mpicom)
      call mpiscatterv(tmp3d_b, numsend, displs, mpir8, tm(:,:,beglat:endlat) ,numrecv, mpir8,0,mpicom)
#else
      t3(:,:,:,1) = t3_tmp(:,:,:)
      tl(:,:,:)   = tmp3d_a(:,:,:)
      tm(:,:,:)   = tmp3d_b(:,:,:)
#endif

      deallocate ( tmp3d_a )
      deallocate ( tmp3d_b )

   !---------------------
   ! Process Constituents
   !---------------------

   case ('CONSTS')

      if (.not. present(m_cnst)) then
         call endrun(sub//': ERROR:  m_cnst needs to be present in the'// &
                      ' argument list')
      end if

      allocate(tmp3d_extend(plon,plev,beglat:endlat))

      ! If "Q", then allocate extra space for spectral truncation
      if (m_cnst == 1) then
         allocate(tmp3d_a(plon,plev,plat))
         allocate(tmp3d_b(plon,plev,plat))
         allocate(tmp3d_c(plon,plev,plat))
      end if

      if (readvar) then

         ! Check that all tracer units are in mass mixing ratios
         ret = pio_inq_varid(fh, cnst_name(m_cnst), varid)
         ret = pio_get_att(fh, varid, 'units', trunits)
         if (trunits(1:5) .ne. 'KG/KG' .and. trunits(1:5) .ne. 'kg/kg') then
            call endrun(sub//': ERROR:  Units for tracer ' &
                  //trim(cnst_name(m_cnst))//' must be in KG/KG')
         end if

      else
          ! Constituents not read from initial file are initialized by the package that implements them.

         if (m_cnst == 1 .and. moist_physics) then
            call endrun(sub//': ERROR:  Q must be on Initial File')
         end if

         call cnst_init_default(m_cnst, clat, clon(:,1), arr3d_a)
      end if

!$omp parallel do private(lat)
      do lat = 1,plat
         call qneg3(sub, lat   , plon, plon   ,plev    , &
            m_cnst, m_cnst, qmin(m_cnst) ,arr3d_a(1,1,lat))
      end do

      ! if "Q", "CLDLIQ", or "CLDICE", save off for later use
      if (m_cnst == 1) q3_tmp(:plon,:,:) = arr3d_a(:plon,:,:)

      ! Spectral truncation of "Q" (only to get spectral derivatives)

      if (m_cnst == 1) then
         tmp3d_a(:plon,:,:) = arr3d_a(:plon,:,:)
         call spetru_3d_scalar(tmp3d_a, dl=tmp3d_b, dm=tmp3d_c)
      end if


#if ( defined SPMD )
      numperlat = plnlv
      call compute_gsfactors(numperlat, numrecv, numsend, displs)
      call mpiscatterv (arr3d_a, numsend, displs, mpir8, tmp3d_extend, numrecv, mpir8,0,mpicom)
      q3(:,:,m_cnst,beglat:endlat,1) = tmp3d_extend(:,:,beglat:endlat)
      if (m_cnst == 1) then
         call mpiscatterv(tmp3d_b, numsend, displs, mpir8, ql, numrecv, mpir8,0,mpicom)
         call mpiscatterv(tmp3d_c, numsend, displs, mpir8, qm, numrecv, mpir8,0,mpicom)
      end if
#else
      q3(:plon,:plev,m_cnst,:plat,1) = arr3d_a(:plon,:plev,:plat)
      if (m_cnst == 1) then
         ql(:,:,:) = tmp3d_b(:,:,:)
         qm(:,:,:) = tmp3d_c(:,:,:)
      end if
#endif
      deallocate(tmp3d_extend)
      if (m_cnst == 1) then
         deallocate(tmp3d_a)
         deallocate(tmp3d_b)
         deallocate(tmp3d_c)
      end if

   !-----------
   ! Process PS
   !-----------

   case ('PS')

      allocate ( tmp2d_a(plon,plat) )
      allocate ( tmp2d_b(plon,plat) )

      ! Spectral truncation
      call spetru_ps(ps_tmp, tmp2d_a, tmp2d_b)

#if ( defined SPMD )
      numperlat = plon
      call compute_gsfactors(numperlat, numrecv, numsend, displs)
      call mpiscatterv(tmp2d_a, numsend, displs, mpir8, dpsl, numrecv, mpir8,0,mpicom)
      call mpiscatterv(tmp2d_b, numsend, displs, mpir8, dpsm, numrecv, mpir8,0,mpicom)
#else
      dpsl(:,:) = tmp2d_a(:,:)
      dpsm(:,:) = tmp2d_b(:,:)
#endif

      deallocate ( tmp2d_a )
      deallocate ( tmp2d_b )

   !-------------
   ! Process PHIS
   !-------------

   case ('PHIS')

      allocate ( tmp2d_a(plon,plat) )
      allocate ( tmp2d_b(plon,plat) )

      ! Check for presence of 'from_hires' attribute to decide whether to filter

      if (readvar) then
         ret = pio_inq_varid (fh, 'PHIS', varid)
         call pio_seterrorhandling(fh, PIO_BCAST_ERROR)

         ret = pio_inq_attlen (fh, varid, 'from_hires', attlen)
         if (ret.eq.PIO_NOERR .and. attlen.gt.256) then
            call endrun(sub//': ERROR:  from_hires attribute length is too long')
         end if
         ret = pio_get_att (fh, varid, 'from_hires', text)
         if (ret.eq.PIO_NOERR .and. text(1:4).eq.'true') then
            phis_hires = .true.
            if (masterproc) write(iulog,*) sub//': INFO: Will filter input PHIS: attribute from_hires is true'
         else
            phis_hires = .false.
            if (masterproc) write(iulog,*) sub//': INFO: Will not filter input PHIS: attribute ', &
                  'from_hires is either false or not present'
         end if
         call pio_seterrorhandling(fh, PIO_INTERNAL_ERROR)

      else
         phis_hires = .false.
      end if

      ! Spectral truncation
      call spetru_phis(phis_tmp, phis_hires, phisl=tmp2d_a, phism=tmp2d_b, phi_out=phi)

      ! Compute ln(Ps*) (Ritchie & Tanguay, 1995) in spectral space
      tmp1 = 1._r8/(rair*t0(plev))
      do ii = 1,psp/2
         i = 2*ii - 1
         lnpstar(i  ) = -phi(1,ii)*tmp1
         lnpstar(i+1) = -phi(2,ii)*tmp1
      end do
          
#if ( defined SPMD )
      call mpibcast(lnpstar, psp, mpir8, 0, mpicom)
      numperlat = plon
      call compute_gsfactors(numperlat, numrecv, numsend, displs)
      call mpiscatterv(phis_tmp, numsend, displs, mpir8, phis,  numrecv, mpir8,0,mpicom)
      call mpiscatterv(tmp2d_a,  numsend, displs, mpir8, phisl, numrecv, mpir8,0,mpicom)
      call mpiscatterv(tmp2d_b,  numsend, displs, mpir8, phism, numrecv, mpir8,0,mpicom)
#else
      phis  = phis_tmp
      phisl = tmp2d_a
      phism = tmp2d_b
#endif

      deallocate ( tmp2d_a )
      deallocate ( tmp2d_b )

   end select

end subroutine process_inidat

!=========================================================================================

subroutine global_int()

! Compute global integrals of mass, moisture and geopotential height
! and fix mass of atmosphere

   use commap
   use physconst,       only: gravit
   use hycoef,          only: hyai, ps0
#if ( defined SPMD )
   use mpishorthand
   use spmd_dyn,        only: compute_gsfactors
#endif
   use sld_control_mod, only: pdela, tmassf, fixmas, qmass1, qmass2, &
                              qmassf, zgsint, tmass0

   !---------------------------Local workspace-----------------------------

   integer i,k,lat,ihem,irow  ! grid indices
   real(r8) pdelb(plon,plev)  ! pressure diff between interfaces
                               ! using "B" part of hybrid grid only
   real(r8) pssum             ! surface pressure sum
   real(r8) dotproda          ! dot product
   real(r8) dotprodb          ! dot product
   real(r8) zgssum            ! partial sums of phis
   real(r8) hyad (plev)       ! del (A)
   real(r8) tmassf_tmp        ! Global mass integral
   real(r8) qmass1_tmp        ! Partial Global moisture mass integral
   real(r8) qmass2_tmp        ! Partial Global moisture mass integral
   real(r8) qmassf_tmp        ! Global moisture mass integral
   real(r8) zgsint_tmp        ! Geopotential integral

#if ( defined SPMD )
   integer :: numperlat         ! number of values per latitude band
   integer :: numsend(0:npes-1) ! number of items to be sent
   integer :: numrecv           ! number of items to be received
   integer :: displs(0:npes-1)  ! displacement array
#endif
   character(len=*), parameter :: sub='global_int'
   !-----------------------------------------------------------------------

   if (masterproc) then

      ! Initialize mass and moisture integrals for summation
      ! in a third calculation loop (assures bit-for-bit compare
      ! with non-random history tape).

      tmassf_tmp = 0._r8
      qmass1_tmp = 0._r8
      qmass2_tmp = 0._r8
      zgsint_tmp = 0._r8

      ! Compute pdel from "A" portion of hybrid vertical grid for later use in global integrals
      do k = 1, plev
         hyad(k) = hyai(k+1) - hyai(k)
      end do
      do k = 1, plev
         do i = 1, plon
            pdela(i,k) = hyad(k)*ps0
         end do
      end do

      ! Compute integrals of mass, moisture, and geopotential height
      do irow = 1, plat/2
         do ihem = 1, 2
            if (ihem == 1) then
               lat = irow
            else
               lat = plat - irow + 1
            end if

            ! Accumulate average mass of atmosphere
            call pdelb0(ps_tmp(1,lat), pdelb, plon)
            pssum  = 0._r8
            zgssum = 0._r8
            do i = 1, plon
               pssum  = pssum  + ps_tmp  (i,lat)
               zgssum = zgssum + phis_tmp(i,lat)
            end do
            tmassf_tmp = tmassf_tmp + w(irow)*pssum/plon
            zgsint_tmp = zgsint_tmp + w(irow)*zgssum/plon

            ! Calculate global integrals needed for water vapor adjustment
            do k = 1, plev
               dotproda = 0._r8
               dotprodb = 0._r8
               do i = 1, plon
                  dotproda = dotproda + q3_tmp(i,k,lat)*pdela(i,k)
                  dotprodb = dotprodb + q3_tmp(i,k,lat)*pdelb(i,k)
               end do
               qmass1_tmp = qmass1_tmp + w(irow)*dotproda/plon
               qmass2_tmp = qmass2_tmp + w(irow)*dotprodb/plon
            end do
         end do
      end do

      ! Normalize average mass, height
      tmassf_tmp = tmassf_tmp*.5_r8/gravit
      qmass1_tmp = qmass1_tmp*.5_r8/gravit
      qmass2_tmp = qmass2_tmp*.5_r8/gravit
      zgsint_tmp = zgsint_tmp*.5_r8/gravit
      qmassf_tmp = qmass1_tmp + qmass2_tmp

      ! Globally avgd sfc. partial pressure of dry air (i.e. global dry mass):
      tmass0 = 98222._r8/gravit
      if (adiabatic)   tmass0 =  tmassf_tmp
      if (ideal_phys ) tmass0 =  100000._r8/gravit
      if (aqua_planet) tmass0 = (101325._r8-245._r8)/gravit

      if (masterproc) then
         write(iulog,*) sub//': INFO:'
         write(iulog,*) '  Mass of initial data before correction      = ', tmassf_tmp
         write(iulog,*) '  Dry mass will be held at                    = ', tmass0
         write(iulog,*) '  Mass of moisture after removal of negatives = ', qmassf_tmp
         write(iulog,*) '  Globally averaged geopotential height (m)   = ', zgsint_tmp
      end if

      ! Compute and apply an initial mass fix factor which preserves horizontal
      ! gradients of ln(ps).
      if (.not. moist_physics) then
         fixmas = tmass0/tmassf_tmp
      else
         fixmas = (tmass0 + qmass1_tmp)/(tmassf_tmp - qmass2_tmp)
      end if
      do lat = 1,plat
         do i = 1, plon
            ps_tmp(i,lat) = ps_tmp(i,lat)*fixmas
         end do
      end do

      ! Global integerals
      tmassf = tmassf_tmp
      qmass1 = qmass1_tmp
      qmass2 = qmass2_tmp
      qmassf = qmassf_tmp
      zgsint = zgsint_tmp

   end if   ! end of if-masterproc

#if ( defined SPMD )
   call mpibcast (tmass0,1,mpir8,0,mpicom)
   call mpibcast (tmassf,1,mpir8,0,mpicom)
   call mpibcast (qmass1,1,mpir8,0,mpicom)
   call mpibcast (qmass2,1,mpir8,0,mpicom)
   call mpibcast (qmassf,1,mpir8,0,mpicom)
   call mpibcast (zgsint,1,mpir8,0,mpicom)

   numperlat = plon
   call compute_gsfactors(numperlat, numrecv, numsend, displs)
   call mpiscatterv(ps_tmp, numsend, displs, mpir8, ps(:,beglat:endlat,1), numrecv, &
                    mpir8, 0, mpicom)
#else
   ps(:,:,1) = ps_tmp(:,:)
#endif

end subroutine global_int

!=========================================================================================

subroutine copytimelevels()


   ! Local variables

   integer :: n
   !----------------------------------------------------------------------------

   ! Make all time levels of prognostics contain identical data.
   ! Fields to be convectively adjusted only *require* n3 time
   ! level since copy gets done in linems.
   do n = 2, ptimelevels
      ps(:,:,n)     = ps(:,:,1)
      u3(:,:,:,n)   = u3(:,:,:,1)
      v3(:,:,:,n)   = v3(:,:,:,1)
      t3(:,:,:,n)   = t3(:,:,:,1)
      q3(:,:,:,:,n) = q3(:,:,:,:,1)
      div(:,:,:,n)  = div(:,:,:,1)
   end do
   qm1(:,:,:,:) = q3(:,:,:,:,1)

end subroutine copytimelevels

!=========================================================================================

end module dyn_comp
