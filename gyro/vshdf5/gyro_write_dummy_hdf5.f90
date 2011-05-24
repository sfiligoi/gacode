!------------------------------------------------
! Dummy routines when HDF5 not available
!------------------------------------------------

subroutine gyro_write_initdata_hdf5(datafile,action)

  implicit none
  integer, intent(in) :: action
  character (len=*), intent(in) :: datafile

end subroutine gyro_write_initdata_hdf5

subroutine gyro_write_timedata_hdf5(action)

  implicit none
  integer, intent(in) :: action

end subroutine gyro_write_timedata_hdf5

subroutine gyro_write_timedata_wedge_hdf5(action)

  implicit none
  integer, intent(in) :: action

end subroutine gyro_write_timedata_wedge_hdf5

subroutine write_hdf5_restart

  implicit none

end subroutine write_hdf5_restart
