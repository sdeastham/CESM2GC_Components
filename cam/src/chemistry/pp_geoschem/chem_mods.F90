      module chem_mods
!--------------------------------------------------------------
! ... Basic chemistry parameters and arrays
!--------------------------------------------------------------
      use shr_kind_mod, only : r8 => shr_kind_r8
      use constituents, only : pcnst
      implicit none
      save

      integer, parameter :: ntracersmax = 200    ! Must be equal to nadv_chem
      integer            :: ntracers
      character(len=255) :: tracernames(ntracersmax)
      character(len=255) :: tracerlongnames(ntracersmax)
      real(r8)           :: adv_mass(ntracersmax)
      real(r8)           :: mwratio(ntracersmax)
      real(r8)           :: ref_mmr(ntracersmax)

      ! Short-lived species (i.e. not advected)
      integer, parameter :: nslsmax = 500        ! UNadvected species only
      integer            :: nsls    
      character(len=255) :: slsnames(nslsmax)
      character(len=255) :: slslongnames(nslsmax)
      real(r8)           :: sls_ref_mmr(nslsmax)
      real(r8)           :: slsmwratio(nslsmax)

      ! Mapping between constituents and GEOS-Chem tracers
      integer :: map2gc(pcnst)
      integer :: map2gc_sls(nslsmax)

      ! Mapping from constituents to raw index
      integer :: map2idx(pcnst)

      integer, parameter :: phtcnt = 40, & ! number of photolysis reactions
                            rxntot = 212, & ! number of total reactions
                            gascnt = 172, & ! number of gas phase reactions
                            nabscol = 2, & ! number of absorbing column densities
                            gas_pcnst = 103, & ! number of "gas phase" species
                            nfs = 4, & ! number of "fixed" species
                            relcnt = 0, & ! number of relationship species
                            grpcnt = 0, & ! number of group members
                            nzcnt = 824, & ! number of non-zero matrix entries
                            extcnt = 4, & ! number of species with external forcing
                            clscnt1 = 8, & ! number of species in explicit class
                            clscnt2 = 0, & ! number of species in hov class
                            clscnt3 = 0, & ! number of species in ebi class
                            clscnt4 = 95, & ! number of species in implicit class
                            clscnt5 = 0, & ! number of species in rodas class
                            indexm = 1, & ! index of total atm density in invariant array
                            indexh2o = 4, & ! index of water vapor density
                            clsze = 1, & ! loop length for implicit chemistry
                            rxt_tag_cnt = 95, &
                            enthalpy_cnt = 0
                            !nslvd = 0
     
      integer :: clscnt(5) = 0
      integer :: cls_rxt_cnt(4,5) = 0
      integer :: clsmap(gas_pcnst,5) = 0
      integer :: permute(gas_pcnst,5) = 0
      integer :: diag_map(clscnt4) = 0
!      real(r8) :: adv_mass(gas_pcnst) = 0._r8
      real(r8) :: crb_mass(gas_pcnst) = 0._r8
      real(r8) :: fix_mass(max(1,nfs))
      real(r8), allocatable :: cph_enthalpy(:)
      integer, allocatable :: cph_rid(:)
      integer, allocatable :: rxt_tag_map(:)
      real(r8), allocatable :: pht_alias_mult(:,:)
      character(len=16), allocatable :: rxt_tag_lst(:)
      character(len=16), allocatable :: pht_alias_lst(:,:)
      character(len=16) :: inv_lst(max(1,nfs))
      character(len=16) :: extfrc_lst(max(1,extcnt))
      logical :: frc_from_dataset(max(1,extcnt))
      logical :: is_vector
      logical :: is_scalar

      ! New SLS
      integer :: nslvd
      character(len=255), allocatable :: slvd_lst(:)
      real(r8), allocatable :: slvd_ref_mmr(:)
      end module chem_mods
