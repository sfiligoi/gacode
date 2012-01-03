!-----------------------------------------------------------
! rTRANSP_DO.f90
!
! PURPOSE:
!  This routine takes as input a distributed matrix g and 
!  copies it into a "transposed" matrix gT also distributed 
!  across processors:
!
!             g({ij},k) -> gT({ki},j)
!
!  The curly braces denote a subset of indices.  The 
!  distributed index can be combination of many indices 
!  stacked together, as indicated by the notation.
!
! NOTES:
!  Communication routine originally used in old GYRO.
!-----------------------------------------------------------

subroutine rTRANSP_DO(g,gT)

  use mpi
  use rTRANSP_GLOBALS

  implicit none

  complex, intent(in), dimension(n_k,n_ij_loc) :: g
  complex, intent(inout), dimension(n_j,n_ki_loc) :: gT


  ! Sort g into packages to be sent to each
  ! processor: q_send(:,i_recv)

  s = 0
  p_ij_loc = 0
  do p_ij=1+i_proc,n_ij,n_proc

     i = i_ij(p_ij)
     p_ij_loc = p_ij_loc+1

     do k=1,n_k

        i_recv    = i_rc(k,i)
        s(i_recv) = s(i_recv)+1

        q_send(s(i_recv),i_recv) = g(k,p_ij_loc)

     enddo ! k

  enddo ! p_ij

  ! Do all-to-all exchange of q_send into q_recv

  call MPI_ALLTOALL(q_send, &
       s_dim, &
       MPI_DOUBLE_COMPLEX, &
       q_recv, &
       s_dim, &
       MPI_DOUBLE_COMPLEX, &
       TRANSP_COMM, &
       i_err)

  ! Build "transposed" g (gT) from packages sent by 
  ! all other processors.

  do i_from=0,n_proc-1
     do s0=1,s_dim
        j = j_map(s0,i_from)
        p_ki_loc = p_ki_loc_map(s0,i_from)
        if (p_ki_loc > 0) gT(j,p_ki_loc) = q_recv(s0,i_from)
     enddo
  enddo

  return

end subroutine rTRANSP_DO
