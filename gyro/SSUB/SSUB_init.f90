!---------------------------------------------
! SSUB_init.f90
!
! PURPOSE:
!  Initialize SSUB library.
!
! NOTES:
!  This scheme is most efficient when n_j 
!  is a multiple of n_n.  
!
! REVISIONS
! 08 Aug 02: jc
!  Created.
!---------------------------------------------

subroutine SSUB_init(nn_i,nj_i,nv1_i,nv2_i,COMM)

  use SSUB_private

  !-------------------------------------------
  implicit none
  !
  integer, intent(in) :: nn_i
  integer, intent(in) :: nj_i
  integer, intent(in) :: nv1_i
  integer, intent(in) :: nv2_i
  integer, intent(in) :: COMM
  !-------------------------------------------

  include 'mpif.h'

  SSUB_COMM = COMM

  !----------------------------------------------
  ! Get rank and number of processors from 
  ! SSUB_COMM (which is assumed to already exist).
  !
  call MPI_COMM_RANK(SSUB_COMM,i_proc,i_err)
  call MPI_COMM_SIZE(SSUB_COMM,n_proc,i_err)
  !
  !-----------------------------------------------

  nn  = nn_i
  nj  = nj_i
  nv1 = nv1_i
  nv2 = nv2_i

  jsplit = 1+(nj*nv2-1)/nn

end subroutine SSUB_init
