subroutine gyro_rhs_nonlinear

  use gyro_globals
  use gyro_pointers
  use math_constants

  !---------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------


  if (nonlinear_flag == 0) then
     return
  endif

  ! Is zeroing necessary when transpose has "unused" elements?
  h_tran = 0.0
  gyro_u_tran = 0.0

  call gyro_timer_in('Nonlinear-comm')
  do is=1,n_kinetic
     call fSSUB(h(:,:,:,is),h_tran(:,:,:,is))
     call fSSUB(gyro_u(:,:,:,is),gyro_u_tran(:,:,:,is))
  enddo
  call gyro_timer_out('Nonlinear-comm')

  !------------------
  ! state:  p_ekj,n
  !------------------

  !--------------------------
  ! Stage 2: nonlinear step
  ! 
  call gyro_timer_in('Nonlinear-step')
  if (nl_method == 1) then
     call gyro_nl_direct
  else 
     call gyro_nl_fft
  endif
  call gyro_timer_out('Nonlinear-step')
  !--------------------------

  call gyro_timer_in('Nonlinear-comm')
  do is=1,n_kinetic
     call rSSUB(h_tran(:,:,:,is),rhs(:,:,:,is))
  enddo
  call gyro_timer_out('Nonlinear-comm')

  !-----------------
  ! state: p_nek,j
  !-----------------

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'*[gyro_rhs_nonlinear done]'
  endif

end subroutine gyro_rhs_nonlinear
