subroutine alloc_add(io,n_size,bytes,name)

  use gyro_globals

  implicit none
  !
  integer, intent(in) :: io
  integer, intent(in) :: n_size
  integer, intent(in) :: bytes
  character (len=*), intent(in) :: name
  !
  real :: this_memory

  select case (output_flag)

  case (1)

     if (i_proc == 0) then

        this_memory  = 1.0*n_size*bytes
        total_memory = total_memory+this_memory 

        write(io,10) this_memory/1048576.0,' MB',name

     endif

  end select

10 format(f7.3,a,3x,a)

end subroutine alloc_add
