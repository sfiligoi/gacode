!-----------------------------------------------------
! GKEIGEN_matrix_write.f90
!
! PURPOSE:
!  Write up to 500 columns of the linear time-derivative
!  matrix for use in restarting eigensolve runs.
!  Applicable only in eigensolve mode.
!-----------------------------------------------------

subroutine GKEIGEN_matrix_write(restart_cols,n_last_col)

  use gyro_globals
  use gyro_pointers
  use GKEIGEN_globals
  use mpi

  !---------------------------------
  implicit none
  !
  integer :: n_last_col
  integer :: istate_max
  integer :: irestart
  !
  complex, dimension(500,h_length_loc) :: restart_cols
  complex, dimension(h_length_loc) :: one_col_w
  !---------------------------------

  10 format(1x,'(',g14.6,',',g14.6,')')
  20 format('----------------- Column: ',I6,'  ----------------------')

  If (i_proc == 0) Then
    open(unit=1,file=trim(path)//file_eigen_restart,status='old',position='append')
  EndIf

  istate_max = Min(n_last_col+499,h_length-1)
  Do istate_1= n_last_col, istate_max
    irestart = istate_1 - n_last_col + 1
    If (i_proc == 0) Then
      write (1,20) istate_1 + 1
      one_col_w(:) = restart_cols(irestart,:)
      write (1,10) one_col_w
      write (1,*)
    EndIf
    Do i_proc_w = 1, n_proc_1-1
      If (i_proc == 0) Then
        call MPI_RECV(one_col_w,&
             size(one_col_w),&
             MPI_DOUBLE_COMPLEX,&
             i_proc_w,&
             i_proc_w,&
             GYRO_COMM_WORLD,&
             recv_status,&
             i_err)

        write (1,10) one_col_w
        write (1,*)

      Else If (i_proc == i_proc_w) Then

        one_col_w(:) = restart_cols(irestart,:)
        call MPI_SEND(one_col_w,&
             size(one_col_w),&
             MPI_DOUBLE_COMPLEX,&
             0,&
             i_proc_w,&
             GYRO_COMM_WORLD,&
             i_err)
      EndIf

    EndDo ! i_proc_w

   EndDo ! istate_1

   If (i_proc == 0) Then
     close (unit=1)
     Print *, 'Called GKEIGEN_matrix_write successfully.'
   EndIf

end subroutine GKEIGEN_matrix_write

