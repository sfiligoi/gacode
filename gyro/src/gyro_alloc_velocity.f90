!----------------------------------------------------------
! gyro_alloc_velocity
!
! PURPOSE:
!  Create and destroy velocity-space arrays.
!
! NOTES:
!  flag=0: deallocate
!  flag=1: allocate
!----------------------------------------------------------

subroutine gyro_alloc_velocity(flag)

  use gyro_globals

  integer, intent(in) :: flag

  if (flag == 1 .and. allocated(lambda)) then
     if (i_proc == 0) then
        print *,'WARNING: already allocated arrays in gyro_alloc_velocity'
     endif
     return
  endif

  if (flag == 0 .and. .not.allocated(lambda)) then
     if (i_proc == 0) then
        print *,'WARNING: cannot deallocate arrays in gyro_alloc_velocity'
     endif
     return
  endif

  if (flag == 1) then

     ! Lambda grid 
     allocate(lambda(n_x,n_lambda))
     allocate(s_lambda(n_lambda))
     allocate(w_lambda(n_x,n_lambda))
     allocate(lambda_tp(n_x))
     allocate(lambda_max(n_x))

     ! Weights
     allocate(w_p(n_energy,n_x,n_lambda))

     allocate(class(n_lambda))

  else 

     deallocate(lambda)
     deallocate(s_lambda)
     deallocate(w_lambda)
     deallocate(lambda_tp)
     deallocate(lambda_max)

     deallocate(w_p)

     deallocate(class)

  endif

end subroutine gyro_alloc_velocity
