
      module mo_sim_dat

      private
      public :: set_sim_dat

      contains

      subroutine set_sim_dat

      use chem_mods,     only : clscnt, cls_rxt_cnt, clsmap, permute, adv_mass, fix_mass, crb_mass
      use chem_mods,     only : diag_map
      use chem_mods,     only : phtcnt, rxt_tag_cnt, rxt_tag_lst, rxt_tag_map
      use chem_mods,     only : pht_alias_lst, pht_alias_mult
      use chem_mods,     only : extfrc_lst, inv_lst, slvd_lst
      use chem_mods,     only : enthalpy_cnt, cph_enthalpy, cph_rid
      use cam_abortutils,only : endrun
      use mo_tracname,   only : solsym
      use chem_mods,     only : frc_from_dataset
      use chem_mods,     only : is_scalar, is_vector
      use shr_kind_mod,  only : r8 => shr_kind_r8
      use cam_logfile,   only : iulog

      implicit none

!--------------------------------------------------------------
!      ... local variables
!--------------------------------------------------------------
      integer :: ios

      is_scalar = .true.
      is_vector = .false.

      clscnt(:) = (/     24,     0,     0,   213,     0 /)

      cls_rxt_cnt(:,1) = (/     29,    57,     0,    24 /)
      cls_rxt_cnt(:,4) = (/     23,   175,   309,   213 /)

      solsym(:237) = (/ 'O3              ','O               ','O1D             ','N2O             ','N               ', &
                        'NO              ','NO2             ','NO3             ','HNO3            ','HO2NO2          ', &
                        'N2O5            ','CH4             ','CH3O2           ','CH3OOH          ','CH3OH           ', &
                        'CH2O            ','CO              ','HCN             ','CH3CN           ','H2              ', &
                        'H               ','OH              ','HO2             ','H2O2            ','CLY             ', &
                        'BRY             ','CL              ','CL2             ','CLO             ','OCLO            ', &
                        'CL2O2           ','HCL             ','HOCL            ','CLONO2          ','BRCL            ', &
                        'BR              ','BRO             ','HBR             ','HOBR            ','BRONO2          ', &
                        'C2H2            ','HCOOH           ','HOCH2OO         ','C2H4            ','C2H6            ', &
                        'C2H5O2          ','C2H5OOH         ','CH3CO3          ','CH3COOH         ','CH3CHO          ', &
                        'C2H5OH          ','GLYALD          ','GLYOXAL         ','CH3COOOH        ','EO2             ', &
                        'EO              ','EOOH            ','PAN             ','C3H6            ','C3H8            ', &
                        'C3H7O2          ','C3H7OOH         ','CH3COCH3        ','PO2             ','POOH            ', &
                        'HYAC            ','RO2             ','CH3COCHO        ','ROOH            ','BIGENE          ', &
                        'BIGALK          ','MEK             ','ENEO2           ','MEKO2           ','MEKOOH          ', &
                        'MCO3            ','MVK             ','MACR            ','MACRO2          ','MACROOH         ', &
                        'MPAN            ','ALKNIT          ','NOA             ','HONITR          ','ISOPNITA        ', &
                        'ISOPNITB        ','ISOPNOOH        ','NC4CHO          ','NC4CH2OH        ','TERPNIT         ', &
                        'NTERPOOH        ','ISOP            ','ALKO2           ','ALKOOH          ','BIGALD          ', &
                        'HYDRALD         ','ISOPAO2         ','ISOPBO2         ','HPALD           ','IEPOX           ', &
                        'ISOPNO3         ','ONITR           ','XO2             ','XOOH            ','ISOPOOH         ', &
                        'TOLUENE         ','CRESOL          ','TOLO2           ','TOLOOH          ','BENZENE         ', &
                        'PHENOL          ','BEPOMUC         ','BENZO2          ','PHENO2          ','PHENO           ', &
                        'PHENOOH         ','C6H5O2          ','C6H5OOH         ','BENZOOH         ','BIGALD1         ', &
                        'BIGALD2         ','BIGALD3         ','BIGALD4         ','MALO2           ','TEPOMUC         ', &
                        'BZOO            ','BZOOH           ','BZALD           ','ACBZO2          ','DICARBO2        ', &
                        'MDIALO2         ','PBZNIT          ','XYLENES         ','XYLOL           ','XYLOLO2         ', &
                        'XYLOLOOH        ','XYLENO2         ','XYLENOOH        ','IVOCbb          ','IVOCff          ', &
                        'MTERP           ','BCARY           ','TERPO2          ','TERPOOH         ','TERPROD1        ', &
                        'TERP2O2         ','SVOCbb          ','SVOCff          ','TERPROD2        ','TERP2OOH        ', &
                        'NTERPO2         ','CH3CL           ','CH3BR           ','CFC11           ','CFC12           ', &
                        'CFC113          ','HCFC22          ','CCL4            ','CH3CCL3         ','CF3BR           ', &
                        'CF2CLBR         ','HCFC141B        ','HCFC142B        ','CFC114          ','CFC115          ', &
                        'H1202           ','H2402           ','CHBR3           ','CH2BR2          ','SO2             ', &
                        'DMS             ','H2SO4           ','SOAGff0         ','SOAGff1         ','SOAGff2         ', &
                        'SOAGff3         ','SOAGff4         ','SOAGbb0         ','SOAGbb1         ','SOAGbb2         ', &
                        'SOAGbb3         ','SOAGbb4         ','SOAGbg0         ','SOAGbg1         ','SOAGbg2         ', &
                        'SOAGbg3         ','SOAGbg4         ','soaff1_a1       ','soaff2_a1       ','soaff3_a1       ', &
                        'soaff4_a1       ','soaff5_a1       ','soabb1_a1       ','soabb2_a1       ','soabb3_a1       ', &
                        'soabb4_a1       ','soabb5_a1       ','soabg1_a1       ','soabg2_a1       ','soabg3_a1       ', &
                        'soabg4_a1       ','soabg5_a1       ','soaff1_a2       ','soaff2_a2       ','soaff3_a2       ', &
                        'soaff4_a2       ','soaff5_a2       ','soabb1_a2       ','soabb2_a2       ','soabb3_a2       ', &
                        'soabb4_a2       ','soabb5_a2       ','soabg1_a2       ','soabg2_a2       ','soabg3_a2       ', &
                        'soabg4_a2       ','soabg5_a2       ','so4_a1          ','pom_a1          ','pom_a4          ', &
                        'bc_a1           ','dst_a1          ','ncl_a1          ','num_a1          ','so4_a2          ', &
                        'dst_a2          ','ncl_a2          ','num_a2          ','so4_a3          ','dst_a3          ', &
                        'ncl_a3          ','num_a3          ','bc_a4           ','num_a4          ','O3S             ', &
                        'CO2             ','H2O             ' /)

      adv_mass(:237) = (/    47.998200_r8,    15.999400_r8,    15.999400_r8,    44.012880_r8,    14.006740_r8, &
                             30.006140_r8,    46.005540_r8,    62.004940_r8,    63.012340_r8,    79.011740_r8, &
                            108.010480_r8,    16.040600_r8,    47.032000_r8,    48.039400_r8,    32.040000_r8, &
                             30.025200_r8,    28.010400_r8,    27.025140_r8,    41.050940_r8,     2.014800_r8, &
                              1.007400_r8,    17.006800_r8,    33.006200_r8,    34.013600_r8,   100.916850_r8, &
                             99.716850_r8,    35.452700_r8,    70.905400_r8,    51.452100_r8,    67.451500_r8, &
                            102.904200_r8,    36.460100_r8,    52.459500_r8,    97.457640_r8,   115.356700_r8, &
                             79.904000_r8,    95.903400_r8,    80.911400_r8,    96.910800_r8,   141.908940_r8, &
                             26.036800_r8,    46.024600_r8,    63.031400_r8,    28.051600_r8,    30.066400_r8, &
                             61.057800_r8,    62.065200_r8,    75.042400_r8,    60.050400_r8,    44.051000_r8, &
                             46.065800_r8,    60.050400_r8,    58.035600_r8,    76.049800_r8,    77.057200_r8, &
                             61.057800_r8,    78.064600_r8,   121.047940_r8,    42.077400_r8,    44.092200_r8, &
                             75.083600_r8,    76.091000_r8,    58.076800_r8,    91.083000_r8,    92.090400_r8, &
                             74.076200_r8,    89.068200_r8,    72.061400_r8,    90.075600_r8,    56.103200_r8, &
                             72.143800_r8,    72.102600_r8,   105.108800_r8,   103.094000_r8,   104.101400_r8, &
                            101.079200_r8,    70.087800_r8,    70.087800_r8,   119.093400_r8,   120.100800_r8, &
                            147.084740_r8,   133.141340_r8,   119.074340_r8,   133.100140_r8,   147.125940_r8, &
                            147.125940_r8,   163.125340_r8,   145.111140_r8,   147.125940_r8,   215.240140_r8, &
                            231.239540_r8,    68.114200_r8,   103.135200_r8,   104.142600_r8,    98.098200_r8, &
                            100.113000_r8,   117.119800_r8,   117.119800_r8,   116.112400_r8,   118.127200_r8, &
                            162.117940_r8,   147.125940_r8,   149.118600_r8,   150.126000_r8,   118.127200_r8, &
                             92.136200_r8,   108.135600_r8,   173.140600_r8,   174.148000_r8,    78.110400_r8, &
                             94.109800_r8,   126.108600_r8,   159.114800_r8,   175.114200_r8,   159.114800_r8, &
                            176.121600_r8,   109.101800_r8,   110.109200_r8,   160.122200_r8,    84.072400_r8, &
                             98.098200_r8,    98.098200_r8,   112.124000_r8,   115.063800_r8,   140.134400_r8, &
                            123.127600_r8,   124.135000_r8,   106.120800_r8,   137.112200_r8,   129.089600_r8, &
                            117.078600_r8,   183.117740_r8,   106.162000_r8,   122.161400_r8,   203.165800_r8, &
                            204.173200_r8,   187.166400_r8,   188.173800_r8,   184.350200_r8,   184.350200_r8, &
                            136.228400_r8,   204.342600_r8,   185.234000_r8,   186.241400_r8,   168.227200_r8, &
                            199.218600_r8,   310.582400_r8,   310.582400_r8,   154.201400_r8,   200.226000_r8, &
                            230.232140_r8,    50.485900_r8,    94.937200_r8,   137.367503_r8,   120.913206_r8, &
                            187.375310_r8,    86.467906_r8,   153.821800_r8,   133.402300_r8,   148.910210_r8, &
                            165.364506_r8,   116.948003_r8,   100.493706_r8,   170.921013_r8,   154.466716_r8, &
                            209.815806_r8,   259.823613_r8,   252.730400_r8,   173.833800_r8,    64.064800_r8, &
                             62.132400_r8,    98.078400_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8, &
                            250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8, &
                            250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8, &
                            250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8, &
                            250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8, &
                            250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8, &
                            250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8, &
                            250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8, &
                            250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8,   250.445000_r8, &
                            250.445000_r8,   250.445000_r8,   115.107340_r8,    12.011000_r8,    12.011000_r8, &
                             12.011000_r8,   135.064039_r8,    58.442468_r8,     1.007400_r8,   115.107340_r8, &
                            135.064039_r8,    58.442468_r8,     1.007400_r8,   115.107340_r8,   135.064039_r8, &
                             58.442468_r8,     1.007400_r8,    12.011000_r8,     1.007400_r8,    47.998200_r8, &
                             44.009800_r8,    18.014200_r8 /)

      crb_mass(:237) = (/     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8, &
                             12.011000_r8,    12.011000_r8,    12.011000_r8,    24.022000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,    12.011000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                             24.022000_r8,    12.011000_r8,    12.011000_r8,    24.022000_r8,    24.022000_r8, &
                             24.022000_r8,    24.022000_r8,    24.022000_r8,    24.022000_r8,    24.022000_r8, &
                             24.022000_r8,    24.022000_r8,    24.022000_r8,    24.022000_r8,    24.022000_r8, &
                             24.022000_r8,    24.022000_r8,    24.022000_r8,    36.033000_r8,    36.033000_r8, &
                             36.033000_r8,    36.033000_r8,    36.033000_r8,    36.033000_r8,    36.033000_r8, &
                             36.033000_r8,    36.033000_r8,    36.033000_r8,    36.033000_r8,    48.044000_r8, &
                             60.055000_r8,    48.044000_r8,    48.044000_r8,    48.044000_r8,    48.044000_r8, &
                             48.044000_r8,    48.044000_r8,    48.044000_r8,    48.044000_r8,    48.044000_r8, &
                             48.044000_r8,    60.055000_r8,    36.033000_r8,    48.044000_r8,    60.055000_r8, &
                             60.055000_r8,    60.055000_r8,    60.055000_r8,    60.055000_r8,   120.110000_r8, &
                            120.110000_r8,    60.055000_r8,    60.055000_r8,    60.055000_r8,    60.055000_r8, &
                             60.055000_r8,    60.055000_r8,    60.055000_r8,    60.055000_r8,    60.055000_r8, &
                             60.055000_r8,    60.055000_r8,    60.055000_r8,    60.055000_r8,    60.055000_r8, &
                             84.077000_r8,    84.077000_r8,    84.077000_r8,    84.077000_r8,    72.066000_r8, &
                             72.066000_r8,    72.066000_r8,    72.066000_r8,    72.066000_r8,    72.066000_r8, &
                             72.066000_r8,    72.066000_r8,    72.066000_r8,    72.066000_r8,    48.044000_r8, &
                             60.055000_r8,    60.055000_r8,    72.066000_r8,    48.044000_r8,    84.077000_r8, &
                             84.077000_r8,    84.077000_r8,    84.077000_r8,    84.077000_r8,    60.055000_r8, &
                             48.044000_r8,    84.077000_r8,    96.088000_r8,    96.088000_r8,    96.088000_r8, &
                             96.088000_r8,    96.088000_r8,    96.088000_r8,   156.143000_r8,   156.143000_r8, &
                            120.110000_r8,   180.165000_r8,   120.110000_r8,   120.110000_r8,   120.110000_r8, &
                            120.110000_r8,   264.242000_r8,   264.242000_r8,   108.099000_r8,   120.110000_r8, &
                            120.110000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8, &
                             24.022000_r8,    12.011000_r8,    12.011000_r8,    24.022000_r8,    12.011000_r8, &
                             12.011000_r8,    24.022000_r8,    24.022000_r8,    24.022000_r8,    24.022000_r8, &
                             12.011000_r8,    24.022000_r8,    12.011000_r8,    12.011000_r8,     0.000000_r8, &
                             24.022000_r8,     0.000000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8, &
                            180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8, &
                            180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8, &
                            180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8, &
                            180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8, &
                            180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8, &
                            180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8, &
                            180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8, &
                            180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8,   180.165000_r8, &
                            180.165000_r8,   180.165000_r8,     0.000000_r8,    12.011000_r8,    12.011000_r8, &
                             12.011000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,    12.011000_r8,     0.000000_r8,     0.000000_r8, &
                             12.011000_r8,     0.000000_r8 /)

      fix_mass(:  3) = (/ 0.00000000_r8, 28.0134800_r8, 31.9988000_r8 /)

      clsmap(: 24,1) = (/  235,  12,   4, 152, 153, 154, 155, 156, 164, 165, &
                           157, 162, 163, 158, 159, 160, 161, 166, 167, 168, &
                           169, 236,  25,  26 /)
      clsmap(:213,4) = (/    1,   2,   3,  20,  17,   5,   6,   7,  22,   8, &
                             9,  10,  11,  13,  14,  18,  19,  16,  21,  23, &
                            24, 237,  27,  28,  29,  30,  31,  32,  33,  34, &
                            35,  36,  37,  38,  39,  40,  59,  92,  64,  50, &
                            49,  65,  48,  54,  58,  45,  44,  71,  81,  82, &
                            83,  84,  85,  86,  87,  88,  89,  90,  91,  70, &
                            73,  93,  94,  72,  74,  75, 106, 107, 108, 109, &
                           143, 144,  95,  53, 110, 133, 111, 112, 113, 114, &
                           115, 116, 117, 118, 119, 120, 121, 122, 123, 124, &
                           125, 126, 127, 128, 129, 130, 131, 132, 134, 135, &
                           136, 137, 138, 141, 142, 140, 139, 148, 147, 145, &
                           146, 149, 150, 151,  97,  98,  99, 100,  77,  78, &
                            79,  80,  76,  46,  47,  60,  61,  62,  63,  69, &
                            15,  51,  52,  66,  55,  56,  57,  96,  67,  68, &
                           101, 102, 103, 104, 105,  41,  42,  43, 170, 171, &
                           172, 173, 174, 175, 176, 177, 178, 179, 180, 181, &
                           182, 183, 184, 185, 186, 187, 218, 219, 220, 221, &
                           222, 223, 224, 225, 188, 189, 190, 191, 192, 193, &
                           194, 195, 196, 197, 198, 199, 200, 201, 202, 203, &
                           204, 205, 206, 207, 208, 209, 210, 211, 212, 213, &
                           214, 215, 216, 217, 226, 227, 228, 229, 230, 231, &
                           232, 233, 234 /)

      permute(:213,4) = (/  202, 212, 213, 176, 182,  81, 209, 205, 204, 211, &
                            161, 116,  98, 199, 121,  82,  68, 201, 198, 208, &
                            140, 197, 207,  80, 203,  76,  69, 206, 162, 168, &
                             91, 210, 200, 157, 156, 133, 175, 165, 155, 186, &
                            134, 135, 196, 130, 126, 100, 136,  92, 137, 145, &
                            166, 183, 138, 128, 120, 164, 107, 131, 108, 110, &
                            139, 172, 146, 132, 147, 109,  86,  87, 150, 141, &
                            167, 125, 112, 169,  77,  83,  79,  78, 123, 117, &
                            124, 105, 144,  93, 118, 122,  99, 143, 101, 149, &
                             85, 119, 106,  88, 129, 154, 181,  72,  84, 127, &
                            111, 153, 148, 171, 170,   7,  14,   8,  15, 174, &
                            180, 179, 142, 178, 191, 194,  89,  73, 192, 188, &
                            193, 102, 195, 163, 103,  70, 158, 113, 173, 114, &
                            159,  94, 177, 184, 151,  95,  74, 104, 185, 189, &
                            187,  75, 190,  96, 152,  97, 160, 115,  71,  90, &
                              1,   2,   3,   4,   5,   6,   9,  10,  11,  12, &
                             13,  16,  17,  18,  19,  20,  21,  22,  23,  24, &
                             25,  26,  27,  28,  29,  30,  31,  32,  33,  34, &
                             35,  36,  37,  38,  39,  40,  41,  42,  43,  44, &
                             45,  46,  47,  48,  49,  50,  51,  52,  53,  54, &
                             55,  56,  57,  58,  59,  60,  61,  62,  63,  64, &
                             65,  66,  67 /)

      diag_map(:213) = (/    1,   2,   3,   4,   5,   6,  12,  18,  19,  20, &
                            21,  22,  23,  29,  35,  36,  37,  38,  39,  40, &
                            41,  42,  43,  44,  45,  46,  47,  48,  49,  50, &
                            51,  52,  53,  54,  55,  56,  57,  58,  59,  60, &
                            61,  62,  63,  64,  65,  66,  67,  68,  69,  70, &
                            71,  72,  73,  74,  75,  76,  77,  78,  79,  80, &
                            81,  82,  83,  84,  85,  86,  87,  88,  91,  94, &
                            99, 101, 104, 107, 110, 112, 120, 126, 130, 135, &
                           137, 141, 150, 157, 162, 172, 180, 185, 188, 194, &
                           199, 202, 205, 209, 213, 217, 221, 227, 233, 236, &
                           242, 247, 252, 257, 260, 266, 271, 276, 281, 286, &
                           294, 300, 306, 312, 318, 324, 331, 337, 345, 351, &
                           357, 363, 368, 375, 379, 386, 394, 401, 409, 415, &
                           421, 426, 431, 439, 443, 451, 459, 467, 475, 483, &
                           492, 501, 510, 516, 523, 534, 545, 556, 567, 578, &
                           590, 598, 611, 622, 631, 641, 649, 657, 667, 671, &
                           675, 681, 691, 700, 714, 730, 739, 751, 761, 775, &
                           804, 828, 838, 844, 856, 871, 880, 889, 901, 916, &
                           928, 936, 944, 956, 965, 977, 996,1012,1025,1040, &
                          1063,1084,1106,1138,1156,1188,1200,1209,1253,1273, &
                          1294,1348,1372,1524,1564,1586,1619,1713,1791,1810, &
                          1866,1892,1912 /)

      extfrc_lst(: 16) = (/ 'NO              ','NO2             ','CO              ','SO2             ','SVOCbb          ', &
                            'so4_a1          ','so4_a2          ','pom_a1          ','pom_a4          ','bc_a1           ', &
                            'bc_a4           ','num_a1          ','num_a2          ','num_a4          ','N               ', &
                            'OH              ' /)

      frc_from_dataset(: 16) = (/ .true., .true., .true., .true., .true., &
                                  .true., .true., .true., .true., .true., &
                                  .true., .true., .true., .true., .false., &
                                  .false. /)

      inv_lst(:  3) = (/ 'M               ', 'N2              ', 'O2              ' /)

      if( allocated( rxt_tag_lst ) ) then
         deallocate( rxt_tag_lst )
      end if
      allocate( rxt_tag_lst(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_lst; error = ',ios
         call endrun
      end if
      if( allocated( rxt_tag_map ) ) then
         deallocate( rxt_tag_map )
      end if
      allocate( rxt_tag_map(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_map; error = ',ios
         call endrun
      end if
      rxt_tag_lst(:rxt_tag_cnt) = (/ 'jo2_b           ', 'jo3_a           ', 'jo3_b           ', 'jn2o            ', &
                                     'jno             ', 'jno2            ', 'jn2o5_a         ', 'jn2o5_b         ', &
                                     'jhno3           ', 'jno3_a          ', 'jno3_b          ', 'jho2no2_a       ', &
                                     'jho2no2_b       ', 'jch3ooh         ', 'jch2o_a         ', 'jch2o_b         ', &
                                     'jh2o_a          ', 'jh2o_b          ', 'jh2o_c          ', 'jh2o2           ', &
                                     'jcl2            ', 'jclo            ', 'joclo           ', 'jcl2o2          ', &
                                     'jhocl           ', 'jhcl            ', 'jclono2_a       ', 'jclono2_b       ', &
                                     'jbrcl           ', 'jbro            ', 'jhobr           ', 'jhbr            ', &
                                     'jbrono2_a       ', 'jbrono2_b       ', 'jch3cl          ', 'jccl4           ', &
                                     'jch3ccl3        ', 'jcfcl3          ', 'jcf2cl2         ', 'jcfc113         ', &
                                     'jhcfc22         ', 'jcfc114         ', 'jcfc115         ', 'jhcfc141b       ', &
                                     'jhcfc142b       ', 'jch3br          ', 'jcf3br          ', 'jh1202          ', &
                                     'jh2402          ', 'jcf2clbr        ', 'jchbr3          ', 'jch2br2         ', &
                                     'jco2            ', 'jch4_a          ', 'jch4_b          ', 'jch3cho         ', &
                                     'jpooh           ', 'jch3co3h        ', 'jpan            ', 'jmpan           ', &
                                     'jmacr_a         ', 'jmacr_b         ', 'jmvk            ', 'jc2h5ooh        ', &
                                     'jeooh           ', 'jc3h7ooh        ', 'jrooh           ', 'jacet           ', &
                                     'jmgly           ', 'jxooh           ', 'jonitr          ', 'jalknit         ', &
                                     'jisopnooh       ', 'jnc4cho         ', 'jterpnit        ', 'jnterpooh       ', &
                                     'jnoa            ', 'jhonitr         ', 'jisopooh        ', 'jhyac           ', &
                                     'jglyald         ', 'jmek            ', 'jbigald         ', 'jglyoxal        ', &
                                     'jalkooh         ', 'jmekooh         ', 'jtolooh         ', 'jterpooh        ', &
                                     'jterp2ooh       ', 'jterprd1        ', 'jterprd2        ', 'jbigald1        ', &
                                     'jbepomuc        ', 'jtepomuc        ', 'jbigald2        ', 'jbigald3        ', &
                                     'jbigald4        ', 'jhpald          ', 'jphenooh        ', 'jc6h5ooh        ', &
                                     'jbenzooh        ', 'jbzooh          ', 'jxylolooh       ', 'jxylenooh       ', &
                                     'jsoabb1_a1      ', 'jsoabb2_a1      ', 'jsoabb3_a1      ', 'jsoabb4_a1      ', &
                                     'jsoabb5_a1      ', 'jsoabb1_a2      ', 'jsoabb2_a2      ', 'jsoabb3_a2      ', &
                                     'jsoabb4_a2      ', 'jsoabb5_a2      ', 'jsoaff1_a1      ', 'jsoaff2_a1      ', &
                                     'jsoaff3_a1      ', 'jsoaff4_a1      ', 'jsoaff5_a1      ', 'jsoaff1_a2      ', &
                                     'jsoaff2_a2      ', 'jsoaff3_a2      ', 'jsoaff4_a2      ', 'jsoaff5_a2      ', &
                                     'jsoabg1_a1      ', 'jsoabg2_a1      ', 'jsoabg3_a1      ', 'jsoabg4_a1      ', &
                                     'jsoabg5_a1      ', 'jsoabg1_a2      ', 'jsoabg2_a2      ', 'jsoabg3_a2      ', &
                                     'jsoabg4_a2      ', 'jsoabg5_a2      ', 'usr_O_O2        ', 'O_O3            ', &
                                     'usr_O_O         ', 'O1D_N2          ', 'O1D_O2b         ', 'O1D_H2O         ', &
                                     'O1D_N2Oa        ', 'O1D_N2Ob        ', 'O1D_O3          ', 'O1D_CFC11       ', &
                                     'O1D_CFC12       ', 'O1D_CFC113      ', 'O1D_CFC114      ', 'O1D_CFC115      ', &
                                     'O1D_HCFC22      ', 'O1D_HCFC141B    ', 'O1D_HCFC142B    ', 'O1D_CCL4        ', &
                                     'O1D_CH3BR       ', 'O1D_CF2CLBR     ', 'O1D_CF3BR       ', 'O1D_H1202       ', &
                                     'O1D_H2402       ', 'O1D_CHBR3       ', 'O1D_CH2BR2      ', 'O1D_CH4a        ', &
                                     'O1D_CH4b        ', 'O1D_CH4c        ', 'O1D_H2          ', 'O1D_HCL         ', &
                                     'O1D_HBR         ', 'O1D_HCN         ', 'H_O2            ', 'H_O3            ', &
                                     'H_HO2a          ', 'H_HO2b          ', 'H_HO2c          ', 'OH_O            ', &
                                     'OH_O3           ', 'OH_HO2          ', 'OH_OH           ', 'OH_OH_M         ', &
                                     'OH_H2           ', 'OH_H2O2         ', 'H2_O            ', 'HO2_O           ', &
                                     'HO2_O3          ', 'usr_HO2_HO2     ', 'H2O2_O          ', 'HCN_OH          ', &
                                     'CH3CN_OH        ', 'N_O2            ', 'N_NO            ', 'N_NO2a          ', &
                                     'N_NO2b          ', 'N_NO2c          ', 'NO_O_M          ', 'NO_HO2          ', &
                                     'NO_O3           ', 'NO2_O           ', 'NO2_O_M         ', 'NO2_O3          ', &
                                     'tag_NO2_NO3     ', 'usr_N2O5_M      ', 'tag_NO2_OH      ', 'usr_HNO3_OH     ', &
                                     'NO3_NO          ', 'NO3_O           ', 'NO3_OH          ', 'NO3_HO2         ', &
                                     'tag_NO2_HO2     ', 'HO2NO2_OH       ', 'usr_HO2NO2_M    ', 'CL_O3           ', &
                                     'CL_H2           ', 'CL_H2O2         ', 'CL_HO2a         ', 'CL_HO2b         ', &
                                     'CL_CH2O         ', 'CL_CH4          ', 'CLO_O           ', 'CLO_OHa         ', &
                                     'CLO_OHb         ', 'CLO_HO2         ', 'CLO_CH3O2       ', 'CLO_NO          ', &
                                     'CLO_NO2_M       ', 'CLO_CLOa        ', 'CLO_CLOb        ', 'CLO_CLOc        ', &
                                     'tag_CLO_CLO_M   ', 'usr_CL2O2_M     ', 'HCL_OH          ', 'HCL_O           ', &
                                     'HOCL_O          ', 'HOCL_CL         ', 'HOCL_OH         ', 'CLONO2_O        ', &
                                     'CLONO2_OH       ', 'CLONO2_CL       ', 'BR_O3           ', 'BR_HO2          ', &
                                     'BR_CH2O         ', 'BRO_O           ', 'BRO_OH          ', 'BRO_HO2         ', &
                                     'BRO_NO          ', 'BRO_NO2_M       ', 'BRO_CLOa        ', 'BRO_CLOb        ', &
                                     'BRO_CLOc        ', 'BRO_BRO         ', 'HBR_OH          ', 'HBR_O           ', &
                                     'HOBR_O          ', 'BRONO2_O        ', 'CH3CL_CL        ', 'CH3CL_OH        ', &
                                     'CH3CCL3_OH      ', 'HCFC22_OH       ', 'CH3BR_OH        ', 'CH3BR_CL        ', &
                                     'HCFC141B_OH     ', 'HCFC142B_OH     ', 'CH2BR2_OH       ', 'CHBR3_OH        ', &
                                     'CH2BR2_CL       ', 'CHBR3_CL        ', 'CH4_OH          ', 'usr_CO_OH_b     ', &
                                     'CO_OH_M         ', 'CH2O_NO3        ', 'CH2O_OH         ', 'CH2O_O          ', &
                                     'CH2O_HO2        ', 'CH3O2_NO        ', 'CH3O2_HO2       ', 'CH3O2_CH3O2a    ', &
                                     'CH3O2_CH3O2b    ', 'CH3OH_OH        ', 'CH3OOH_OH       ', 'HCOOH_OH        ', &
                                     'HOCH2OO_M       ', 'HOCH2OO_NO      ', 'HOCH2OO_HO2     ', 'C2H2_CL_M       ', &
                                     'C2H4_CL_M       ', 'C2H6_CL         ', 'C2H2_OH_M       ', 'C2H6_OH         ', &
                                     'tag_C2H4_OH     ', 'EO2_NO          ', 'EO2_HO2         ', 'EO_O2           ', &
                                     'EO_M            ', 'C2H4_O3         ', 'CH3COOH_OH      ', 'C2H5O2_NO       ', &
                                     'C2H5O2_HO2      ', 'C2H5O2_CH3O2    ', 'C2H5O2_C2H5O2   ', 'C2H5OOH_OH      ', &
                                     'CH3CHO_OH       ', 'CH3CHO_NO3      ', 'CH3CO3_NO       ', 'tag_CH3CO3_NO2  ', &
                                     'CH3CO3_HO2      ', 'CH3CO3_CH3O2    ', 'CH3CO3_CH3CO3   ', 'CH3COOOH_OH     ', &
                                     'GLYALD_OH       ', 'GLYOXAL_OH      ', 'C2H5OH_OH       ', 'usr_PAN_M       ', &
                                     'PAN_OH          ', 'tag_C3H6_OH     ', 'C3H6_O3         ', 'C3H6_NO3        ', &
                                     'C3H7O2_NO       ', 'C3H7O2_HO2      ', 'CH3H7O2_CH3O2   ', 'CH3H7OOH_OH     ', &
                                     'C3H8_OH         ', 'PO2_NO          ', 'PO2_HO2         ', 'POOH_OH         ', &
                                     'usr_CH3COCH3_OH ', 'RO2_NO          ', 'RO2_HO2         ', 'RO2_CH3O2       ', &
                                     'ROOH_OH         ', 'HYAC_OH         ', 'CH3COCHO_OH     ', 'CH3COCHO_NO3    ', &
                                     'NOA_OH          ', 'BIGENE_OH       ', 'ENEO2_NO        ', 'ENEO2_NOb       ', &
                                     'BIGENE_NO3      ', 'MVK_OH          ', 'MVK_O3          ', 'MEK_OH          ', &
                                     'MEKO2_NO        ', 'MEKO2_HO2       ', 'MEKOOH_OH       ', 'MACR_OH         ', &
                                     'MACR_O3         ', 'MACRO2_NOa      ', 'MACRO2_NOb      ', 'MACRO2_NO3      ', &
                                     'MACRO2_HO2      ', 'MACRO2_CH3O2    ', 'MACRO2_CH3CO3   ', 'MACROOH_OH      ', &
                                     'MCO3_NO         ', 'MCO3_NO3        ', 'MCO3_HO2        ', 'MCO3_CH3O2      ', &
                                     'MCO3_CH3CO3     ', 'MCO3_MCO3       ', 'usr_MCO3_NO2    ', 'usr_MPAN_M      ', &
                                     'MPAN_OH_M       ', 'HONITR_OH       ', 'ISOP_OH         ', 'ISOP_O3         ', &
                                     'ISOP_NO3        ', 'ISOPAO2_NO      ', 'ISOPBO2_NO      ', 'ISOPNITA_OH     ', &
                                     'ISOPNITB_OH     ', 'ISOPAO2_NO3     ', 'ISOPBO2_NO3     ', 'ISOPAO2_HO2     ', &
                                     'ISOPBO2_HO2     ', 'ISOPOOH_OH      ', 'IEPOX_OH        ', 'ISOPAO2_CH3O2   ', &
                                     'ISOPBO2_CH3O2   ', 'ISOPAO2_CH3CO3  ', 'ISOPBO2_CH3CO3  ', 'ISOPBO2_M       ', &
                                     'HPALD_OH        ', 'ISOPNO3_NO      ', 'ISOPNO3_NO3     ', 'ISOPNO3_HO2     ', &
                                     'ISOPNO3_CH3CO3  ', 'ISOPNO3_CH3O2   ', 'ISOPNOOH_OH     ', 'NC4CH2OH_OH     ', &
                                     'NC4CHO_OH       ', 'HYDRALD_OH      ', 'XO2_NO          ', 'XO2_NO3         ', &
                                     'XO2_HO2         ', 'XO2_CH3O2       ', 'XO2_CH3CO3      ', 'XOOH_OHa        ', &
                                     'usr_XOOH_OH     ', 'BIGALK_OH       ', 'ALKO2_NO        ', 'ALKO2_NOb       ', &
                                     'ALKNIT_OH       ', 'ALKO2_HO2       ', 'ALKOOH_OH       ', 'BENZENE_OH      ', &
                                     'PHENOL_OH       ', 'PHENO2_NO       ', 'PHENO2_HO2      ', 'PHENOOH_OH      ', &
                                     'PHENO_NO2       ', 'PHENO_O3        ', 'C6H5O2_NO       ', 'C6H5O2_HO2      ', &
                                     'C6H5OOH_OH      ', 'BENZO2_NO       ', 'BENZO2_HO2      ', 'BENZOOH_OH      ', &
                                     'MALO2_NO2       ', 'MALO2_NO        ', 'MALO2_HO2       ', 'TOLUENE_OH      ', &
                                     'CRESOL_OH       ', 'BZOO_HO2        ', 'BZOOH_OH        ', 'BZOO_NO         ', &
                                     'BZALD_OH        ', 'tag_ACBZO2_NO2  ', 'usr_PBZNIT_M    ', 'ACBZO2_NO       ', &
                                     'ACBZO2_HO2      ', 'TOLO2_HO2       ', 'TOLO2_OH        ', 'TOLO2_NO        ', &
                                     'DICARBO2_HO2    ', 'DICARBO2_NO     ', 'MDIALO2_HO2     ', 'MDIALO2_NO      ', &
                                     'DICARBO2_NO2    ', 'MDIALO2_NO2     ', 'XYLENES_OH      ', 'XYLOL_OH        ', &
                                     'XYLOLO2_NO      ', 'XYLOLO2_HO2     ', 'XYLOLOOH_OH     ', 'XYLENO2_HO2     ', &
                                     'XYLENOOH_OH     ', 'XYLENO2_NO      ', 'MTERP_OH        ', 'BCARY_OH        ', &
                                     'MTERP_O3        ', 'BCARY_O3        ', 'MTERP_NO3       ', 'BCARY_NO3       ', &
                                     'TERPO2_NO       ', 'TERPO2_HO2      ', 'TERPO2_CH3O2    ', 'TERPOOH_OH      ', &
                                     'TERP2OOH_OH     ', 'TERPROD1_OH     ', 'TERPROD1_NO3    ', 'TERP2O2_NO      ', &
                                     'TERP2O2_HO2     ', 'TERP2O2_CH3O2   ', 'TERPROD2_OH     ', 'NTERPO2_NO      ', &
                                     'NTERPO2_HO2     ', 'NTERPO2_CH3O2   ', 'NTERPO2_NO3     ', 'TERPNIT_OH      ', &
                                     'NTERPOOH_OH     ', 'usr_N2O5_aer    ', 'usr_NO3_aer     ', 'usr_NO2_aer     ', &
                                     'usr_HO2_aer     ', 'usr_GLYOXAL_aer ', 'usr_SO2_OH      ', 'DMS_OHa         ', &
                                     'usr_DMS_OH      ', 'DMS_NO3         ', 'het1            ', 'het2            ', &
                                     'het3            ', 'het4            ', 'het5            ', 'het6            ', &
                                     'het7            ', 'het8            ', 'het9            ', 'het10           ', &
                                     'het11           ', 'het12           ', 'het13           ', 'het14           ', &
                                     'het15           ', 'het16           ', 'het17           ' /)
      rxt_tag_map(:rxt_tag_cnt) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                                       11,  12,  13,  14,  15,  16,  17,  18,  19,  20, &
                                       21,  22,  23,  24,  25,  26,  27,  28,  29,  30, &
                                       31,  32,  33,  34,  35,  36,  37,  38,  39,  40, &
                                       41,  42,  43,  44,  45,  46,  47,  48,  49,  50, &
                                       51,  52,  53,  54,  55,  56,  57,  58,  59,  60, &
                                       61,  62,  63,  64,  65,  66,  67,  68,  69,  70, &
                                       71,  72,  73,  74,  75,  76,  77,  78,  79,  80, &
                                       81,  82,  83,  84,  85,  86,  87,  88,  89,  90, &
                                       91,  92,  93,  94,  95,  96,  97,  98,  99, 100, &
                                      101, 102, 103, 104, 105, 106, 107, 108, 109, 110, &
                                      111, 112, 113, 114, 115, 116, 117, 118, 119, 120, &
                                      121, 122, 123, 124, 125, 126, 127, 128, 129, 130, &
                                      131, 132, 133, 134, 135, 136, 137, 138, 139, 140, &
                                      141, 142, 143, 144, 145, 146, 147, 148, 149, 150, &
                                      151, 152, 153, 154, 155, 156, 157, 158, 159, 160, &
                                      161, 162, 163, 164, 165, 166, 167, 168, 169, 170, &
                                      171, 172, 173, 174, 175, 176, 177, 178, 179, 180, &
                                      181, 182, 183, 184, 185, 186, 187, 188, 189, 190, &
                                      191, 192, 193, 194, 195, 196, 197, 198, 199, 200, &
                                      201, 202, 203, 204, 205, 206, 207, 208, 209, 210, &
                                      211, 212, 213, 214, 215, 216, 217, 218, 219, 220, &
                                      221, 222, 223, 224, 225, 226, 227, 228, 229, 230, &
                                      231, 232, 233, 234, 235, 236, 237, 238, 239, 240, &
                                      241, 242, 243, 244, 245, 246, 247, 248, 249, 250, &
                                      251, 252, 253, 254, 255, 256, 257, 258, 259, 260, &
                                      261, 262, 263, 264, 265, 266, 267, 268, 269, 270, &
                                      271, 272, 273, 274, 275, 276, 277, 278, 279, 280, &
                                      281, 282, 283, 284, 285, 286, 287, 288, 289, 290, &
                                      291, 292, 293, 294, 295, 296, 297, 298, 299, 300, &
                                      301, 302, 303, 304, 305, 306, 307, 308, 309, 310, &
                                      311, 312, 313, 314, 315, 316, 317, 318, 319, 320, &
                                      321, 322, 323, 324, 325, 326, 327, 328, 329, 330, &
                                      331, 332, 333, 334, 335, 336, 337, 338, 339, 340, &
                                      341, 342, 343, 344, 345, 346, 347, 348, 349, 350, &
                                      351, 352, 353, 354, 355, 356, 357, 358, 359, 360, &
                                      361, 362, 363, 364, 365, 366, 367, 368, 369, 370, &
                                      371, 372, 373, 374, 375, 376, 377, 378, 379, 380, &
                                      381, 382, 383, 384, 385, 386, 387, 388, 389, 390, &
                                      391, 392, 393, 394, 395, 396, 397, 398, 399, 400, &
                                      401, 402, 403, 404, 405, 406, 407, 408, 409, 410, &
                                      411, 412, 413, 414, 415, 416, 417, 418, 419, 420, &
                                      421, 422, 423, 424, 425, 426, 427, 428, 429, 430, &
                                      431, 432, 433, 434, 435, 436, 437, 438, 439, 440, &
                                      441, 442, 443, 444, 445, 446, 447, 448, 449, 450, &
                                      451, 452, 453, 454, 455, 456, 457, 458, 459, 460, &
                                      461, 462, 463, 464, 465, 482, 483, 484, 485, 486, &
                                      487, 488, 489, 490, 491, 492, 493, 494, 495, 496, &
                                      497, 498, 499, 500, 501, 502, 503, 504, 505, 506, &
                                      507 /)
      if( allocated( pht_alias_lst ) ) then
         deallocate( pht_alias_lst )
      end if
      allocate( pht_alias_lst(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate pht_alias_lst; error = ',ios
         call endrun
      end if
      if( allocated( pht_alias_mult ) ) then
         deallocate( pht_alias_mult )
      end if
      allocate( pht_alias_mult(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate pht_alias_mult; error = ',ios
         call endrun
      end if
      pht_alias_lst(:,1) = (/ 'userdefined     ', '                ', '                ', '                ', &
                              'userdefined     ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ' /)
      pht_alias_lst(:,2) = (/ '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              'jch3ooh         ', 'jh2o2           ', '                ', 'jpan            ', &
                              '                ', '                ', '                ', 'jch3ooh         ', &
                              'jch3ooh         ', 'jch3ooh         ', 'jch3ooh         ', '                ', &
                              '                ', 'jch3ooh         ', 'jch3cho         ', 'jch3ooh         ', &
                              'jch3ooh         ', 'jch2o_a         ', 'jch3ooh         ', 'jch3ooh         ', &
                              'jch2o_a         ', 'jch2o_a         ', 'jch3ooh         ', '                ', &
                              '                ', 'jacet           ', 'jno2            ', 'jmgly           ', &
                              'jch3ooh         ', 'jch3ooh         ', 'jch3ooh         ', 'jch3ooh         ', &
                              'jch3ooh         ', 'jch3cho         ', 'jch3cho         ', 'jno2            ', &
                              'jno2            ', 'jno2            ', 'jno2            ', 'jno2            ', &
                              'jno2            ', 'jno2            ', 'jch3ooh         ', 'jch3ooh         ', &
                              'jch3ooh         ', 'jch3ooh         ', 'jch3ooh         ', 'jch3ooh         ', &
                              'jno2            ', 'jno2            ', 'jno2            ', 'jno2            ', &
                              'jno2            ', 'jno2            ', 'jno2            ', 'jno2            ', &
                              'jno2            ', 'jno2            ', 'jno2            ', 'jno2            ', &
                              'jno2            ', 'jno2            ', 'jno2            ', 'jno2            ', &
                              'jno2            ', 'jno2            ', 'jno2            ', 'jno2            ', &
                              'jno2            ', 'jno2            ', 'jno2            ', 'jno2            ', &
                              'jno2            ', 'jno2            ', 'jno2            ', 'jno2            ', &
                              'jno2            ', 'jno2            ' /)
      pht_alias_mult(:,1) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8 /)
      pht_alias_mult(:,2) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 0.28_r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 0.2_r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, .14_r8, .10_r8, .10_r8, .20_r8, &
                          .20_r8, .006_r8, .006_r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, .0004_r8, &
                          .0004_r8, .0004_r8, .0004_r8, .0004_r8, .0004_r8, &
                          .0004_r8, .0004_r8, .0004_r8, .0004_r8, .0004_r8, &
                          .0004_r8, .0004_r8, .0004_r8, .0004_r8, .0004_r8, &
                          .0004_r8, .0004_r8, .0004_r8, .0004_r8, .00004_r8, &
                          .00004_r8, .00004_r8, .00004_r8, .00004_r8, .00004_r8, &
                          .00004_r8, .00004_r8, .00004_r8, .00004_r8 /)

      end subroutine set_sim_dat

      end module mo_sim_dat
