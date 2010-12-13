subroutine func(y,yp,nv)

  use input_data

  implicit none

  integer, intent(in) :: nv
  real, intent(in), dimension(nv) :: y
  real, intent(inout), dimension(nv) :: yp

  integer :: i_n
  integer :: p
  integer :: j
  integer :: jj
  integer :: n_xd
  integer :: n_nd
  complex, dimension(:,:), allocatable :: anp_in
  real, dimension(:), allocatable :: g_theta_in ! WG
  complex :: kern
  complex :: kernp
  complex :: kernm


  !-----------------------------------------------------

  if (n_n == 1) then
     allocate(anp_in(-n_x/2:n_x/2-1,0:1))
  else
     allocate(anp_in(-n_x/2:n_x/2-1,-n_n+1:n_n-1))
  endif

  allocate(g_theta_in(-n_x/2:n_x/2-1)) ! WG

  j = 1+(y(3)+pi)/dtheta

  if (j == n_theta_plot+1) j=j-1

  jj = j+1 

  ! Perform linear interpolation in theta.

  anp_in(:,:) = ((y(3)-theta(j))*anp(jj,:,:)+&
       (theta(jj)-y(3))*anp(j,:,:))/dtheta

  ! Linearly interpolate G_theta(theta) ! WG
  g_theta_in(:) = ((y(3)-theta(j))*g_theta(jj,:)+& ! WG
       (theta(jj)-y(3))*g_theta(j,:))/dtheta       ! WG


  yp(:) = 0.0

  n_xd = n_x-n_x_decimate
  n_nd = n_n-n_n_decimate

  if (n_n == 1) then

     do p=-n_xd/2,n_xd/2-1

        kernp = i_c*lambda*g_theta_in(p)*anp_in(p,0)*exp(-i_c*(y(2)-p*y(1))) ! WG
        kernm = i_c*lambda*g_theta_in(p)*anp_in(p,1)*exp(-i_c*(-y(2)-p*y(1))) ! WG

        yp(1) = yp(1)+real((kernp-kernm))
        yp(2) = yp(2)+real((kernp+kernm)*p)

     enddo ! p

  else

     do i_n=-(n_nd-1),n_nd-1

        do p=-n_xd/2,n_xd/2-1

           kern = i_c*lambda*g_theta_in(p)*anp_in(p,i_n)*exp(-i_c*(i_n*y(2)-p*y(1))) ! WG

           yp(1) = yp(1)+real(kern*i_n)
           yp(2) = yp(2)+real(kern*p)

        enddo ! p

     enddo ! i_n

  endif

  yp(3) = 1.0

  deallocate(anp_in)
  deallocate(g_theta_in) ! WG

end subroutine func
