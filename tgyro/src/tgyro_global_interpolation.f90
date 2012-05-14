!-----------------------------------------------------------------
! tgyro_global_inerpolation.f90
!
! PURPOSE:
!  Linearly interpolates TGYRO "coarse" z-profile onto "fine"
!  EXPRO grid, then calls LOGINT to calculate corresponding
!  profiles on fine grid.  Used to generate updated input.profiles
!  in global GYRO-TGYRO runs.
!
! NOTES:
!  i -> coarse
!  j -> fine
!  Values of z on coarse grid are zero at i=1.
!-----------------------------------------------------------------

subroutine tgyro_global_interpolation(x_coarse,z_coarse,t_coarse,n_coarse,n_fine,x_fine,t_fine)

  implicit none

  integer, intent(in) :: n_coarse
  real, intent(in), dimension(n_coarse) :: x_coarse
  real, intent(in), dimension(n_coarse)  :: z_coarse
  real, intent(in), dimension(n_coarse)  :: t_coarse

  integer, intent(in) :: n_fine
  real, intent(inout), dimension(n_fine) :: x_fine
  real, intent(inout), dimension(n_fine) :: t_fine

  real, dimension(n_fine) :: z_fine
  integer :: i,j,i_bc

  !set z_fine(x > MAX(x_coarse)) = z_coarse(MAX(x_coarse)
  !set BC for integration at x = MAX(x_coarse)
  do j=n_fine,1,-1
     if (x_fine(j) >= x_coarse(n_coarse)) then
        z_fine(j) = z_coarse(n_coarse)
        i_bc = j
     else
        do i=1,n_coarse-1
           if (x_fine(j) >= x_coarse(i) .and. x_fine(j) < x_coarse(i+1)) then
              z_fine(j) = z_coarse(i)*(x_coarse(i+1)-x_fine(j))+z_coarse(i+1)*(x_fine(j)-x_coarse(i))
              z_fine(j) = z_fine(j)/(x_coarse(i+1)-x_coarse(i))
           endif
        enddo ! i
     endif
  enddo ! j

! CH: use these calls to put BC at outermost TGYRO radial instance
!  t_fine = t_coarse(n_coarse)
!  call logint(t_fine,z_fine,x_fine,n_fine,i_bc)

! CH: use the calls to put logint BC at r=0
! seems to reproduce TGYRO profiles slightly better
  t_fine(1) = t_coarse(1)
  call logint(t_fine,z_fine,x_fine,n_fine,1)

end subroutine tgyro_global_interpolation
