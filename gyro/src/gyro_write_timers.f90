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
  character (len=9), dimension(64) :: a,b
  integer :: loop
  real :: frac
  !-------------------------------------------------

  select case (io_control)

  case(0)

     return

  case(1)

     ! Initial open
     if (i_proc == 0) then
        open(unit=io,file=datafile,status='replace')
        do loop=1,cpu_maxindx
           a(loop) = cpu_tag(loop)(1:scan(cpu_tag(loop),'-',.false.)-1)
           b(loop) = cpu_tag(loop)(scan(cpu_tag(loop),'-',.false.)+1:len(trim(cpu_tag(loop))))
        enddo
        write(io,'(64(a))') (a(loop)//' ',loop=1,cpu_maxindx)
        write(io,'(64(a))') (b(loop)//' ',loop=1,cpu_maxindx)
        write(io,'(64(a))') ('--------- ',loop=1,cpu_maxindx)
        close(io)
     endif

  case(2)

     if (step == 0) then
        cpu(:) = 0.0
        frac = 0.0
     else
        frac = sum(cpu(1:cpu_maxindx-1))/cpu(cpu_maxindx)
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
        read(io,'(a)') cpu_tag(1)
        do loop=0,data_step
           read(io,10) frac
        enddo
        endfile(io)
        close(io)

     endif

  end select

  if (i_proc == 0 .and. debug_flag == 1) print *,'[gyro_write_timers called]'

10 format(64(1pe9.3,1x))

end subroutine gyro_write_timers
