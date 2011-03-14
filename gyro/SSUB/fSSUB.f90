!---------------------------------------------
! fSSUB.f90
!
! PURPOSE:
!  Within a processor subgroup (defined in 
!  SSUB_init), with one value of n per 
!  processor in the subgroup, collect all
!  n's on one processor and distribute the 
!  j's.  Since n_j is typically larger than
!  n_n, we need jsplit rows is j's.
!            
!       x(:,:,n_j) -> xT(:,:,jsplit,n_n)
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

subroutine fSSUB(x_IN,xt)

  use mpi
  use SSUB_private

  !-------------------------------------------------------
  implicit none
  !
  integer :: j,p,s
  complex, intent(in), dimension(nj,nv1,nv2) :: x_IN
  complex, dimension(nv1,jsplit*nn) :: x
  complex, intent(inout), dimension(nv1,jsplit,nn) :: xt
  !-------------------------------------------------------


  s = 0
  do j=1,nj
     do p=1,nv2
        s = s+1 
        x(:,s) = x_IN(j,:,p)
     enddo
  enddo

  do s=nj*nv2+1,jsplit*nn
     x(:,s) = (0.0,0.0)
  enddo

  call MPI_ALLTOALL(x(:,:), &
       nv1*jsplit, &
       MPI_DOUBLE_COMPLEX, &
       xt(:,:,:), &
       nv1*jsplit, &
       MPI_DOUBLE_COMPLEX, &
       SSUB_COMM, &
       i_err)

end subroutine fSSUB
