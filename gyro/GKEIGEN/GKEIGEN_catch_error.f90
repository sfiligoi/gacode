!---------------------------------------
! GKEIGEN_catch_error.f90
!
! PURPOSE:
!  Check for self consistency in all GKEIGEN
!  input parameters. All inconsistencies are 
!  checked before terminating.
!---------------------------------------

subroutine GKEIGEN_catch_error

  use gyro_globals
  use gyro_pointers
  use GKEIGEN_globals

  logical :: l_terminate = .false.

  character (len=80) :: error_message

10 format(A, I7, '.')
20 format(A, I7, A, I7, A)
30 format(A, I3, A)

!-- Temporarily disable the seldom used matrix reading/writing.
!   When it comes, the new implementation will be based on
!   pre-made PETSc matrix reading/writing.
!
  If ((gkeigen_matrixonly==1) .or. (gkeigen_mwrite_flag==1)) &
    call send_line('Matrix writing and restart temporarily disabled.')

  If (gkeigen_matrixonly == 1) Then
    call send_line('GKEIGEN_MATRIXONLY=1 means simulation will abort.')
    l_terminate = .true.
  EndIf
  If (gkeigen_mwrite_flag == 1) Then
   gkeigen_mwrite_flag = 0
   call send_line('GKEIGEN_MWRITE_FLAG reset to 0.')
  EndIf
      
  eigensolve_restart_flag = 0
!
!---------------------------------------

!-------------------------------------------------
! Check for consistent dimensions of the eigenvalue problem.
!
  If (gkeigen_kspace_dim > h_length) Then
    write(error_message,10) 'INCONSISTENT: GKEIGEN_KSPACE_DIM must be <= matrix dimension ', h_length
    call send_line(error_message)
    write(error_message,10) 'GKEIGEN_KSPACE_DIM reset to ', h_length
    call send_line(error_message)
    gkeigen_kspace_dim = h_length
  EndIf
  If (gkeigen_n_values > gkeigen_kspace_dim) Then
    call send_line('INCONSISTENT: GKEIGEN_N_VALUES must be <= GKEIGEN_KSPACE_DIM.') 
    write(error_message,10) 'GKEIGEN_N_VALUES reset to ', gkeigen_kspace_dim
    call send_line(error_message)
    gkeigen_n_values = gkeigen_kspace_dim
  EndIf

!--------------------------------------------------

!--------------------------------------------------
! Check for consistent processor counts. GKEIGEN must
! have a uniform distribution of matrix elements across
! processors. This implies several restrictions.
!
  If (mod(h_length,n_proc_tot)/=0) Then
    write(error_message,20) 'Matrix dimension (', h_length, &
      ') must be an integer multiple of total processors (', n_proc_tot, ')'
    call send_line(error_message)
    l_terminate = .true.
  EndIf

  If ((mod(n_proc_tot,gkeigen_proc_mult)/=0) .or. (h_length_loc*n_proc/=h_length)) Then
    write(error_message,20) 'Error(s) in n_proc = Total processors (', n_proc_tot, &
      ') / GKEIGEN_PROC_MULT (', gkeigen_proc_mult, '):'
    call send_line(error_message)
    l_terminate = .true.
  EndIf
  If (mod(n_proc_tot,gkeigen_proc_mult)/=0) &
    call send_line('   * n_proc must be an integer.')
  If (h_length_loc*n_proc/=h_length) Then
    write(error_message,30) '   * n_proc must divide evenly into ', n_energy*n_lambda/2, &
         ' (from velocity grid).'
    call send_line(error_message)
  EndIf

  If (l_terminate) call catch_error('ABORTING JOB (see above).')

end subroutine
