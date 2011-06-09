!-----------------------------------------------------
! gyro_write_precision.f90
!
! PURPOSE:
!  This routine prints high-precision "checksum" data.
!-----------------------------------------------------

subroutine gyro_write_precision(io,checksum)

  use mpi
  use gyro_globals
  use gyro_pointers

  !---------------------------------------------------
  implicit none
  !
  integer, intent(in) :: io
  real, intent(in) :: checksum
  !---------------------------------------------------

  if (output_flag == 0) return

  if (i_proc == 0) then

     if (step == 0) then
        open(unit=io,file=trim(precfile),status='replace')
     else
        open(unit=io,file=trim(precfile),status='old',position='append')
        write(io,10) checksum
        close(io)
     endif

  endif

  if (i_proc == 0 .and. debug_flag == 1) print *,'[gyro_write_precision called]'

10 format(t2,4(1pe22.15,1x))

end subroutine gyro_write_precision
