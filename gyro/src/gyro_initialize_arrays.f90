!----------------------------------------------------
! gyro_initialize_arrays.f90
!
! PURPOSE:
!  Post-allocation initialization of arrays.
!------------------------------------------------

subroutine gyro_initialize_arrays

  use gyro_globals

  !---------------------------------------------------
  implicit none
  !---------------------------------------------------

  !---------------------------------------------------
  ! Do not reset distribution if we are
  ! making a "warm restart".
  !
  if (lskipinit == 0) h = 0.0
  !
  h_old = 0.0
  h_0   = 0.0
  h_err = 0.0
  !
  RHS    = 0.0
  RHS_dr = 0.0
  RHS_dt = 0.0
  !---------------------------------------------------

  !---------------------------------------------------
  ! Source-related variables
  !
  h0_mod = 0.0
  h0_n   = 0.0
  h0_e   = 0.0
  source_n = 0.0
  source_e = 0.0
  !---------------------------------------------------

  !---------------------------------------------------
  freq_n       = 0.0
  time_error   = 0.0
  diff         = 0.0
  elapsed_time = 0.0

  !  field_blend     = 0.0
  field_blend_dot  = 0.0
  field_blend_old  = 0.0
  field_blend_old2 = 0.0

  !  cap_h     = 0.0
  h_cap_dot  = 0.0
  h_cap_old  = 0.0
  h_cap_old2 = 0.0

  a_fluxave = 0.0

  gyro_uv_old  = 0.0
  gyro_uv_old2 = 0.0

  if (transport_method == 2) then
     diff_vec = 0.0
     gbflux_vec = 0.0
  endif

  entropy(:,:) = 0.0
  !---------------------------------------------------

  !---------------------------------------------------
  ! Initialization to prevent MPI-grid errors
  field_tau = 0.0
  if (collision_flag == 1) then
     f_coll  = 0.0
     fb_coll = 0.0
  endif
  !---------------------------------------------------

  if (i_proc == 0 .and. debug_flag == 1) then
     print *,'[gyro_initialize_arrays done]' 
  endif

end subroutine gyro_initialize_arrays
