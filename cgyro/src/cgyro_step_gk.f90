subroutine cgyro_step_gk

  use timer_lib
  use cgyro_globals
  use cgyro_globals_math

  implicit none

  ! RK4 time-advance for the distribution
  !
  !           z e             vpar            z e  vperp^2
  !  h = H - ----- G0 ( phi - ----- Apar ) + ----- ---------- Gperp Bpar
  !            T               c               T   omega_a c
  !
  ! After time advance, we will have 
  !
  ! h    -> h_x
  ! H    -> cap_h_c
  ! phi  -> field(1)
  ! Apar -> field(2)
  ! Bpar -> field(3)

  call timer_lib_in('str_mem')
  call cgyro_vel_copy(h0_x, h_x)
  call timer_lib_out('str_mem')

  
  ! Stage 1
  call cgyro_rhs(1)
  call timer_lib_in('str')
  call cgyro_vel_fma2(h_x, h0_x, 0.5 * delta_t, rhs(:,:,1))
  call timer_lib_out('str')
  call cgyro_field_c

  ! Stage 2
  call cgyro_rhs(2)
  call timer_lib_in('str')
  call cgyro_vel_fma2(h_x, h0_x, 0.5 * delta_t, rhs(:,:,2))
  call timer_lib_out('str')
  call cgyro_field_c

  ! Stage 3
  call cgyro_rhs(3)
  call timer_lib_in('str')
  call cgyro_vel_fma2(h_x, h0_x, delta_t, rhs(:,:,3))
  call timer_lib_out('str')
  call cgyro_field_c

  ! Stage 4
  call cgyro_rhs(4)
  call timer_lib_in('str')
  call cgyro_vel_fmaN(4, h_x, &
          h0_x, &
          (/ delta_t/6, 2*delta_t/6, 2*delta_t/6, delta_t/6 /), &
          rhs(:,:,1:4))
  call timer_lib_out('str')
  call cgyro_field_c

  ! rhs(1) = 3rd-order error estimate
  call timer_lib_in('str')
  ! cannot use cgyro_vel_fmaN, as the 3 arrays are not contiguous
  call cgyro_vel_fma4(rhs(:,:,1), &
          h0_x, &
          delta_t/3, rhs(:,:,2), &
          2*delta_t/3, rhs(:,:,3), &
          -1.0, h_x)
  call timer_lib_out('str')

end subroutine cgyro_step_gk

