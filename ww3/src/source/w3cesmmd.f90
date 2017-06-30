!/ ------------------------------------------------------------------- /
      MODULE W3CESMMD

!/ ------------------------------------------------------------------- /
!/
      private

      ! runtype is used by W3SRCE (values are startup, branch, continue)
      character(len=16),public :: runtype

      ! if a run is a startup or branch run, then initfile is used
      ! to construct the initial file and used in W3IORSMD
      character(len=256), public :: initfile

      ! if a run is a continue run, then casename is used to construct
      ! the restart filename in W3IORSMD
      character(len=256), public :: casename

      logical, public :: rstwr   ! true => write restart at end of day
      logical, public :: histwr  ! true => write history file (snapshot)

      integer, public :: stdout  ! output log file
!/
!/ End of module W3CESMMD -------------------------------------------- /
!/
      END MODULE W3CESMMD
