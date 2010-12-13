!-----------------------------------------------------
! gyro_write_precision.f90
!
! PURPOSE:
!  This routine prints high-precision "checksum" data.
!-----------------------------------------------------

subroutine gyro_write_precision(io,checksum)

  use gyro_globals
  use gyro_pointers

  !---------------------------------------------------
  implicit none
  !
  integer, intent(in) :: io
  integer :: io_mode
  real, intent(in) :: checksum
  !---------------------------------------------------

  include 'mpif.h'

  select case (output_flag)

  case (0)

     io_mode = 0

  case (1)

     io_mode = 1

  end select

  if (step == 0) then

     if (i_proc == 0 .and. io_mode == 1) then
        open(unit=io,file=trim(precfile),status='replace')
     endif

  else

     if (i_proc == 0 .and. io_mode == 1) then

        open(unit=io,file=trim(precfile),status='old',position='append')
        write(io,10) checksum
        close(io)

     endif

  endif

  if (i_proc == 0 .and. debug_flag == 1) print *,'[gyro_write_precision called]'

10 format(t2,4(1pe22.15,1x))

end subroutine gyro_write_precision
