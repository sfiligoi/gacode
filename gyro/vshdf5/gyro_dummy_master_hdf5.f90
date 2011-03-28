!------------------------------------------------
! write_hdf5_master_hdf5.f90 
!
! PURPOSE:
!  Write a bunch of data to hdf5 file
!------------------------------------------------
subroutine write_hdf5_data(datafile,action)
  !------------------------------------------
  !  Data that does not change with time.  
  !  It is equivalent to:
  !    profile_vugyro.out
  !------------------------------------------
  use gyro_globals
  !------------------------------------------
  implicit none
  integer, intent(in) :: action
  character (len=*), intent(in) :: datafile
  !
 return
 end subroutine write_hdf5_data

!------------------------------------------------------
! write_hdf5_timedata
! PURPOSE:
!  This is an hdf5 version of gyro_write_master.f90
!-----------------------------------------------------

subroutine write_hdf5_timedata(action)
  use gyro_globals

  !---------------------------------------------------
  implicit none
  include 'mpif.h'
  !
  integer :: mode
  integer, intent(in) :: action
  return
 
end subroutine write_hdf5_timedata

!------------------------------------------------------
! write_hdf5_fine_timedata
! PURPOSE:
!  This is like the above, only it is for just the 
!  fine data
!-----------------------------------------------------

subroutine write_hdf5_fine_timedata(action)
  use gyro_globals

  !---------------------------------------------------
  implicit none
  include 'mpif.h'
  !
  integer :: mode
  integer, intent(in) :: action
  return

end subroutine write_hdf5_fine_timedata

  !------------------------------------------------
  ! write_restart
  ! PURPOSE:
  !  File that can be used 
  !------------------------------------------------
subroutine write_hdf5_restart
  use gyro_globals
  !---------------------------------------------------
  implicit none
 return
 end subroutine write_hdf5_restart
