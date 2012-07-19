!--------------------------------------------------------------
! prgen_read_ufilef90
!
! PURPOSE:
!  Extract data from UFILE using external library (read
!--------------------------------------------------------------

program prgen_read_ufile

  !use prgen_read_globals
  use data_interface

  implicit none

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

  !open(unit=1,file=raw_data_file,status='old')
  !read(1,*) tok
  !read(1,*) shot
  !close(1)

  tok   = 'd3d'
  shot  = '90108'
  phase = ''
  cudir = '/home/staebler/TESTCASES/90108'
  xp_time = 3.09
  endtime = 1000.0
  time_series = 0
  itorque = 1
  iptot = 1
  mxgrid = 50
  nsmooth = 0
  idatzero = 0
  iproc_d = 0

  call readufiles(&
       tok,&
       shot,&
       phase,&
       cudir,&
       xp_time,&
       endtime,&
       time_series,&
       itorque,&
       iptot,&
       mxgrid,&
       nsmooth,&
       idatzero,&
       iproc_d)

end program prgen_read_ufile

