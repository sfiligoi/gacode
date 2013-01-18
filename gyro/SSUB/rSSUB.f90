!---------------------------------------------
! rSSUB.f90
!
! PURPOSE:
!  Reverse the process in fSSUB
!
! NOTES:
!  This scheme is most efficient when n_j 
!  is a multiple of n_n.  
!
! REVISIONS
! 08 Aug 02: jc
!  Created.
! 04 Nov 05: jc
!  New blocking scheme.
!---------------------------------------------

subroutine rSSUB(xt,x_OUT)

  use mpi
  use SSUB_private

  !-------------------------------------------------------
  implicit none
  !
  integer :: j,p,s
  complex, intent(in), dimension(nv1,jsplit,nn) :: xt
  complex, intent(inout), dimension(nj,nv1,nv2) :: x_OUT
  complex, dimension(nv1,jsplit*nn) :: x
  !-------------------------------------------------------  


  call MPI_ALLTOALL(xt, &
       nv1*jsplit, &
       MPI_DOUBLE_COMPLEX, &
       x, &
       nv1*jsplit, &
       MPI_DOUBLE_COMPLEX, &
       SSUB_COMM, &
       i_err)

!$omp parallel do private(s) schedule(static)
  do j=1,nj
     do p=1,nv2
        s = p + (j - 1)*nv2
        x_OUT(j,:,p) = x(:,s)
     enddo
  enddo

end subroutine rSSUB
