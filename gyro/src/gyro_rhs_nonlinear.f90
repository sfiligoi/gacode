subroutine gyro_rhs_nonlinear

  use gyro_globals
  use gyro_pointers
  use math_constants

  !---------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------

  if (nonlinear_flag == 0) then
     CPU_NL      = 0.0
     CPU_NLt     = 0.0
     CPU_NL_in   = 0.0
     CPU_NL_out  = 0.0
     CPU_NLt_in  = 0.0
     CPU_NLt_out = 0.0
     return
  endif

  call proc_time(CPU_NL_in)

  ! Is zeroing necessary when transpose has "unused" elements?
  h_tran = 0.0
  gyro_u_tran = 0.0

  do is=1,n_kinetic
     call fSSUB(h(:,:,:,is),h_tran(:,:,:,is))
     call fSSUB(gyro_u(:,:,:,is),gyro_u_tran(:,:,:,is))
  enddo

  !------------------
  ! state:  p_ekj,n
  !------------------

  call proc_time(CPU_NLt_in)

  !--------------------------
  ! Stage 2: nonlinear step
  ! 
  if (nl_method == 1) then
     call gyro_nl_direct
  else 
     call gyro_nl_fft
  endif
  !--------------------------

  call proc_time(CPU_NLt_out)
  CPU_NLt = CPU_NLt + (CPU_NLt_out - CPU_NLt_in)

  do is=1,n_kinetic
     call rSSUB(h_tran(:,:,:,is),rhs(:,:,:,is))
  enddo

  !-----------------
  ! state: p_nek,j
  !-----------------

  call proc_time(CPU_NL_out)
  CPU_NL = CPU_NL + (CPU_NL_out - CPU_NL_in)

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'*[gyro_rhs_nonlinear done]'
  endif

end subroutine gyro_rhs_nonlinear
