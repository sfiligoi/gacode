!------------------------------------------------------------
! gyro_write_timers.f90
!
! PURPOSE:
!  Control final calculation and output of code timing data.
!------------------------------------------------------------
 
subroutine gyro_write_timers(datafile,io)

  use mpi
  use gyro_globals

  !-----------------------------------------------
  implicit none
  !
  integer, intent(in) :: io
  character (len=*), intent(in) :: datafile
  integer :: loop
  character (len=11) :: sep='--------------'
  real :: frac
  !-------------------------------------------------

  select case (io_control)

  case(0)

     return

  case(1)

     ! Initial open
     if (i_proc == 0) then
        open(unit=io,file=datafile,status='replace')
        write(io,'(64(a))') (cpu_tag(loop)//' ',loop=1,cpu_maxindx)
        write(io,'(64(a))') (sep//' ',loop=1,cpu_maxindx)
        close(io)
     endif

  case(2)

     if (step == 0) then
        cpu(:) = 0.0
        frac = 0.0
     else
        frac = sum(cpu(2:cpu_maxindx))/cpu(1)
     endif
     if (i_proc == 0) then
        open(unit=io,file=datafile,status='old',position='append')
        write(io,10) cpu(1:cpu_maxindx),frac
        close(io)
     endif

     ! Reset all timers
     cpu(1:cpu_maxindx) = 0.0

  case(3)

     ! Reposition after restart

     if (i_proc == 0) then

        open(unit=io,file=datafile,status='old')
        read(io,'(a)') cpu_tag(1)
        read(io,'(a)') cpu_tag(1)
        do loop=0,data_step
           read(io,10) frac
        enddo
        endfile(io)
        close(io)

     endif

  end select

  if (i_proc == 0 .and. debug_flag == 1) print *,'[gyro_write_timers called]'

10 format(64(1pe11.5,1x))

end subroutine gyro_write_timers
