subroutine get_nonlinear_advance

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

  do is=1,n_kinetic
     call fSSUB(h(:,:,:,is),h_M(:,:,:,is))
     call fSSUB(gyro_u(:,:,:,is),gyro_u_M(:,:,:,is))
  enddo

  !------------------
  ! state:  p_ekj,n
  !------------------

  call proc_time(CPU_NLt_in)

  !--------------------------
  ! Stage 2: nonlinear step
  ! 
  if (boundary_method == 1) then

     if (nl_method == 1) then
        call do_nlfast_p
     else 
        call do_nlfft_p
     endif

  else

     if (nl_method == 1) then
        call do_nlfast_match
     else 
        call do_nlfft_match
     endif

  endif
  !--------------------------

  call proc_time(CPU_NLt_out)
  CPU_NLt = CPU_NLt + (CPU_NLt_out - CPU_NLt_in)

  do is=1,n_kinetic
     call rSSUB(h_M(:,:,:,is),rhs(:,:,:,is))
  enddo

  !-----------------
  ! state: p_nek,j
  !-----------------

  call proc_time(CPU_NL_out)
  CPU_NL = CPU_NL + (CPU_NL_out - CPU_NL_in)

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'*[get_nonlinear_advance done]'
  endif

end subroutine get_nonlinear_advance
