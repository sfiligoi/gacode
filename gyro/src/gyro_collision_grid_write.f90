subroutine gyro_collision_grid_write(datafile,io)

  use gyro_globals
  use gyro_collision_private
  use math_constants

  !----------------------------------------
  implicit none
  !
  integer, intent(in) :: io
  character (len=*), intent(in) :: datafile
  !
  integer, parameter :: n_r_write = 1
  !
  integer :: i_r_write(n_r_write)
  integer :: i_stencil
  !----------------------------------------

  select case (output_flag)

  case (1)

     !--------------------------------------
     ! Specify at which radii to print the 
     ! stencil data:
     !
     i_r_write(1) = ir_norm
     !---------------------------------------

     open(unit=io,file=datafile,status='replace')

     write(io,*) n_lambda
     write(io,*) n_pass
     write(io,*) n_stack
     write(io,*) n_theta_section
     write(io,*) n_r_write
     write(io,10) i_r_write(:)

     do i_stencil=1,n_r_write
        i = i_r_write(i_stencil)
        do k=1,n_lambda
           do m=1,n_stack
              write(io,*) xi(i,k,m)
           enddo
        enddo
     enddo
     do i_stencil=1,n_r_write
        i = i_r_write(i_stencil)
        do k=1,n_lambda
           do m=1,n_stack
              write(io,*) theta_t(i,k,m)/pi
           enddo
        enddo
     enddo

     close(io)

  end select

  if (debug_flag == 1) print *,'[gyro_collision_grid_write called]'

10 format(4(i4,1x))

end subroutine gyro_collision_grid_write
