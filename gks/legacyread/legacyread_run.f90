subroutine data_run(idata,shot,tok,cudir,xp_time,endtime,&
     time_series,itorque,iptot,ncl_flag,&
     mxgrid,ismooth,idatzero,iproc_d)

  use legacyread_interface

  implicit none

  integer :: idata,itorque,iptot,mxgrid,ismooth
  integer :: iproc_d, time_series, ncl_flag
  integer :: idatzero
  character(len=50) :: cudir
  character(len=40) :: shot
  character(len=6) :: phase
  character(len=10) :: tok
  real :: xp_time, endtime  

  if (idata == 0) then
     call readufiles(tok,shot,phase,cudir,&
          xp_time,endtime,time_series, itorque,iptot,&
          mxgrid,ismooth,idatzero,iproc_d)
  else if (idata == 1) then
     call readiterdb(tok,shot,cudir,iptot,itorque,ncl_flag)
     call gridsetup(mxgrid)
  else
     write(*,*) 'Error: idata out of range',idata
     stop
  endif

  if (ismooth /= 0) call datavg(ismooth,ncl_flag)

end subroutine data_run
