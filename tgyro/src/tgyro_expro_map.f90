! r (IN)
! z (IN)
! n_r (IN)
!
! r_exp (IN)
! p_exp (OUT)
! n_exp

subroutine tgyro_expro_map(r,z,n_r,r_exp,p_exp,n_exp)

  integer, intent(in) :: n_r,n_exp
  real, intent(in), dimension(n_r) :: r,z
  real, intent(in), dimension(n_exp) :: r_exp
  real, intent(inout), dimension(n_exp) :: p_exp
  real, dimension(n_exp) :: z_exp
  integer :: i,i_exp


  ! Compute z's on exp grid
  z_exp(:) = 0.0
  i = 1
  do i_exp=2,n_exp
     if (r(i+1) < r_exp(i_exp)) i = i+1
     dr = r(i+1)-r(i)
     if (r_exp(i_exp) < r(n_r)) then
        z_exp(i_exp) = z(i+1)*(r_exp(i_exp)-r(i))/dr + z(i)*(r(i+1)-r_exp(i_exp))/dr
     else
        exit
     endif
  enddo

  do i_exp=n_exp,2,-1
     if (r_exp(i_exp) < r(n_r)) then
        p_exp(i_exp-1) = p_exp(i_exp)*exp(0.5*(z_exp(i_exp)+z_exp(i_exp-1))* &
             (r_exp(i_exp)-r_exp(i_exp-1)))
     endif
  enddo

end subroutine tgyro_expro_map
