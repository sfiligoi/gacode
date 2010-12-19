!---------------------------------------------------------
! make_lambda_grid.f90
!
! PURPOSE:
!  Lambda grid (and tau/theta spacing) setup.
!---------------------------------------------------------

subroutine make_lambda_grid

  use gyro_globals
  use math_constants

  !------------------------------------------
  implicit none
  !
  real :: s_tp
  !
  real, dimension(:), allocatable :: s_temp
  real, dimension(:), allocatable :: w_temp
  !------------------------------------------

  if (n_lambda == 1) then
     lambda(:,:) = 0.0
     w_lambda(:,:) = 0.5
     return
  endif

  do i=1,n_x

     call gyro_to_geo(i)

     call gyro_banana_init(nint_ORB_s)

     ! Passing

     call gyro_banana_getlambda(lambda_tp(i),lambda_max(i))

     ! two regions

     call ORB_s(lambda_tp(i),s_tp)

     ! s_tp => s at trapped-passing boundary

     allocate(s_temp(n_pass))
     allocate(w_temp(n_pass))
     call gauss_legendre(0.0,s_tp,s_temp,w_temp,n_pass)
     s_lambda(1:n_pass) = s_temp
     w_lambda(i,1:n_pass) = w_temp
     deallocate(s_temp)
     deallocate(w_temp)

     allocate(s_temp(n_trap))
     allocate(w_temp(n_trap))
     call gauss_legendre(s_tp,1.0,s_temp,w_temp,n_trap)
     s_lambda(n_pass+1:n_lambda) = s_temp
     w_lambda(i,n_pass+1:n_lambda) = w_temp
     deallocate(s_temp)
     deallocate(w_temp)

     do k=1,n_lambda
        call ORB_s2lambda(s_lambda(k),lambda(i,k))
     enddo

     ! Two signs of velocity

     do k=1,n_lambda
        w_lambda(i,k) = 0.5*w_lambda(i,k)
     enddo

  enddo ! i


  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[make_lambda_grid done]'
  endif

end subroutine make_lambda_grid
