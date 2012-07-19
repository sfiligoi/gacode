!--------------------------------------------------------------
! prgen_read_ufilef90
!
! PURPOSE:
!  Extract data from UFILE using external library (read
!--------------------------------------------------------------

subroutine prgen_read_ufile

  use prgen_read_globals
  !use data_interface

  implicit none

  character (len=100) :: t
  
  ! UFILE parameters
  character (len=10) :: tok
  character (len=40) :: shot
  character (len=6) :: phase
  character (len=50) :: cudir
  real :: xp_time
  real :: endtime
  integer :: time_series
  integer :: itorque
  integer :: iptot
  integer :: mxgrid
  integer :: nsmooth
  integer :: idatzero
  integer :: iproc_d

  open(unit=1,file=raw_data_file,status='old')
  read(1,*) tok
  read(1,*) shot
  close(1)

  tok   = trim(tok)
  shot  = trim(shot)
  phase = ''
  cudir = './'
  xp_time = 2.1
  endtime = 3.0
  time_series = 0
  itorque = 0
  iptot = 0
  mxgrid = 20
  nsmooth = 2
  idatzero = 1 
  iproc_d = 0

  !call readufiles(&
  !     tok,&
  !     shot,&
  !     phase,&
  !     cudir,&
  !     xp_time,&
  !     endtime,&
  !     time_series,&
  !     itorque,&
  !     iptot,&
  !     mxgrid,&
  !     nsmooth,&
  !     idatzero,&
  !     iproc_d)

end subroutine prgen_read_ufile

