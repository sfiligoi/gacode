!-----------------------------------------------------------------
! tgyro_comm_sync.f90
!
! PURPOSE:
!  Synchronization (gather, broadcast) of variables returned by 
!  calls to flux routines.
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

  t_fine = t_coarse(n_coarse)

  call logint(t_fine,z_fine,x_fine,n_fine,i_bc)

end subroutine tgyro_global_interpolation
