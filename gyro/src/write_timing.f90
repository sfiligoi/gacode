!------------------------------------------------------------
! write_timing.f90
!
! PURPOSE:
!  Control final calculation and output of code timing data.
!------------------------------------------------------------
 
subroutine write_timing(datafile,io)

  use mpi
  use gyro_globals

  !-----------------------------------------------
  implicit none
  !
  integer, intent(in) :: io
  !
  integer, parameter :: n_col=11
  !
  character (len=*), intent(in) :: datafile
  character (len=10), dimension(n_col) :: tag1
  character (len=10), dimension(n_col) :: tag2
  character (len=10), dimension(n_col) :: sep
  !
  character (len=1) :: dummy
  integer :: data_loop
  !
  real :: z_out(n_col)
  !-------------------------------------------------

  tag1(1) = 'Nonlinear'
  tag2(1) = 'Eval.'

  tag1(2) = 'Nonlinear'
  tag2(2) = 'Commun.'

  tag1(3) = 'Collision'
  tag2(3) = 'Eval.'

  tag1(4) = 'Collision'
  tag2(4) = 'Commun.'

  tag1(5) = 'RHS'
  tag2(5) = 'Linear'

  tag1(6) = 'Field'
  tag2(6) = 'Solve'

  tag1(7) = 'Field'
  tag2(7) = 'Interp'

  tag1(8) = 'Field'
  tag2(8) = 'Vel. Sum'

  tag1(9) = 'Diag.'
  tag2(9) = 'Extras'

  tag1(10) = 'I/O'
  tag2(10) = 'Calls'

  tag1(11) = 'INTERVAL'
  tag2(11) = 'TOTAL'

  sep(1) = '--------- '
  sep(:) = sep(1)

  select case (io_control)

  case(0)

     return

  case(1)

     ! Initial open

     if (i_proc == 0) then

        open(unit=io,file=datafile,status='replace')
        write(io,10) tag1(:)
        write(io,10) tag2(:)
        write(io,10) sep(:)
        close(io)

     endif

  case(2)

     ! Nonlinear  (NL)
     z_out(1) = CPU_NLt
     CPU_NLt = 0.0

     ! Nonlinear transpose (NL_tr)
     z_out(2) = CPU_NL-z_out(1)
     CPU_NL = 0.0

     ! Collision (Coll)
     z_out(3) = CPU_Ct
     CPU_Ct = 0.0

     ! Collision transpose (Coll_tr)
     z_out(4) = CPU_C-z_out(3)
     CPU_C = 0.0

     ! Linear RHS
     z_out(5) = CPU_RHS
     CPU_RHS = 0.0

     ! Field_core
     z_out(6) = CPU_field2
     CPU_field2 = 0.0

     ! Field_interp
     z_out(7) = CPU_interp
     CPU_interp = 0.0

     ! Field_vel_sum
     z_out(8) = CPU_field-z_out(6)-z_out(7)
     CPU_field = 0.0

     ! Extras
     z_out(9) = CPU_diag_a+ &
          (CPU_ts-(z_out(1)+z_out(2)+sum(z_out(5:8))))
     CPU_diag_a = 0.0
     CPU_ts = 0.0

     ! I/O
     z_out(10) = CPU_diag_b
     CPU_diag_b = 0.0

     call system_clock(clock_count,clock_rate,clock_max)

     elapsed_time = clock_count*1.0/clock_rate-elapsed_time

     z_out(11) = elapsed_time

     !----------------------------------------
     ! Exceptions
     !
     if (step == 0) z_out(:) = 0.0
     if (nonlinear_flag == 0) z_out(1:2) = 0.0
     if (collision_flag == 0) z_out(3:4) = 0.0
     !----------------------------------------

     if (i_proc == 0) then
        open(unit=io,file=datafile,status='old',position='append')
        write(io,15) z_out(:)
        close(io)
     endif

     elapsed_time = clock_count*1.0/clock_rate

  case(3)

     ! Reposition after restart

     if (i_proc == 0) then

        open(unit=io,file=datafile,status='old')

        read(io,20) dummy
        read(io,20) dummy
        read(io,20) dummy
        do data_loop=0,data_step
           read(io,20) dummy
        enddo

        endfile(io)
        close(io)

     endif

  end select

  if (i_proc == 0 .and. debug_flag == 1) print *,'[write_timing called]'

10 format(t2,12(a))
15 format(t2,12(es9.3,1x))
20 format(a)

end subroutine write_timing
