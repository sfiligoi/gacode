!--------------------------------------------------------------
! get_field_r0_plot.f90
!
! PURPOSE:
!  Generate interpolation of fields at r=r0 suitable for 
!  plotting.
!--------------------------------------------------------------

subroutine get_field_r0_plot

  use gyro_globals
  use gyro_pointers
  use math_constants

  !---------------------------------------------------
  implicit none
  !
  integer :: j_plot
  complex, dimension(field_r0_grid) :: field_tmp
  !---------------------------------------------------

  if (alltime_index == 0) then
     field_r0_plot(:,:) = (0.0,0.0)
  endif

  i = ir_norm

  !---------------------------------------------------------
  do ix=1,n_field

     do j_plot=1,field_r0_grid

        field_tmp(j_plot) = sum(field_blend(:,i,ix)*blend_r0_plot(:,j_plot))

     enddo ! j_plot

     field_r0_plot(j_plot,ix) = field_r0_plot(j_plot,ix)+&
          w_time(alltime_index+1)*field_tmp(j_plot)

  enddo ! ix
  !---------------------------------------------------------

  if (i_proc == 0 .and. debug_flag == 1) then
     print *,'[get_field_r0_plot done]' 
  endif

end subroutine get_field_r0_plot
