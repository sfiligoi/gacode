subroutine make_theta_grid

  use gyro_globals
  use math_constants

  implicit none

  !--------------------------------------------------------
  ! NOTE: 
  !
  !  (1) tau's go all the way around a CLOSED orbit.
  !  (2) tau is NORMALIZED according to: 
  !
  !              max[tau(1)] = 1 (passing)
  !              max[tau(2)] = 2 (trapped)
  !
  !                        / 
  !      This ensures that | dtau = 1 for a given sign 
  !                        / 
  !      of velocity.
  !
  !  (3) n_theta_section is the number of thetas in an 
  !      "orbit section".  An orbit section is that length 
  !      required to do an orbit integral -- 1/2 the orbit
  !      for passing, and 1/4 the orbit for trapped.
  !   
  !
  !  (4) n_theta(1) and n_theta(2) are the number of 
  !      different values of theta touched by an orbit
  !      of given type.
  !
  ! Passing [ n_tau = n_theta ]:
  !
  n_theta(1) = 2*n_theta_section-2
  n_tau(1)   = n_theta(1)
  !
  ! Trapped [ n_tau = 2(n_theta-1) ]:
  !
  n_theta(2) = 2*n_theta_section-1
  n_tau(2)   = 2*n_theta(2)-2
  !
  ! ** Note that n_stack is the same for both 
  ! types of particles:
  !
  !  n_stack_pass = 2*n_theta_pass = 4*n_theta_section-4
  !  n_stack_trap = n_tau_trap = 4*n_theta_section-4
  !
  n_stack = n_tau(2)
  !
  ! These should be equal, despite separate calculations!
  !
  d_tau(1) = 1.0/n_tau(1)
  d_tau(2) = 2.0/n_tau(2)
  !-------------------------------------------------------

  d_theta = pi_2/n_blend

  ! Real theta grid (arbitrary size)
  n_theta_int = 2*n_blend

  ! Dimension of RBF collision matrices
  !
  n_rbf = n_stack*n_lambda

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[make_theta_grid done]'
  endif

end subroutine make_theta_grid
