!-----------------------------------------------------
! BLEND_do.f90
!
! REVISIONS
! 22 Dec 00: jc
!  Created.
! 14 Nov 03: jc
!  This routine is no longer used by GYRO.
!---------------------------------------------

subroutine BLEND_do(c_IN)

  use BLEND_private

  implicit none

  complex, intent(in), dimension(n_fit) :: c_IN

  ! Block-cyclic LU decomposition and solve:

  c0 = c_IN

  call ZGETRS('N',n_fit,1,cs,n_fit,i_piv,c0,n_fit,info)

end subroutine BLEND_do



