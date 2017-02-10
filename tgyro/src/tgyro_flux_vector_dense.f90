subroutine tgyro_flux_vector_dense(x_vec,f_vec)

  use tgyro_globals

  implicit none

  real, intent(in), dimension(p_max) :: x_vec
  real, intent(inout), dimension(p_max) :: f_vec

  call tgyro_profile_set(x_vec,0.0,0)

  ! Make profile consistent with new gradients
  call tgyro_profile_functions
  call tgyro_flux
  call tgyro_comm_sync
  call tgyro_flux_set(f_vec)

end subroutine tgyro_flux_vector_dense
