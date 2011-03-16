!-----------------------------------------------------
! fTRANSP_INIT.f90
!
! PURPOSE:
!  Initialization stage for calls fTRANSP.
!
! NOTES:
!  Communication routine originally used in old GYRO.
!-----------------------------------------------------

subroutine fTRANSP_INIT(n0_i,n0_j,n0_k,COMM)

  use mpi
  use fTRANSP_GLOBALS

  implicit none

  integer, intent(in) :: n0_i
  integer, intent(in) :: n0_j
  integer, intent(in) :: n0_k
  integer, intent(in) :: COMM

  integer, external :: parallel_dim


  TRANSP_COMM = COMM

  !----------------------------------------------
  ! Get rank and number of processors from 
  ! MPI_COMM_WORLD (which is assumed to already 
  ! exist).
  !
  call MPI_COMM_RANK(TRANSP_COMM,i_proc,i_err)
  call MPI_COMM_SIZE(TRANSP_COMM,n_proc,i_err)
  !
  !-----------------------------------------------

  !----------------------------
  ! Allocate send-pointer array 
  !
  allocate(s(0:n_proc-1))
  !
  !----------------------------

  n_i = n0_i
  n_j = n0_j
  n_k = n0_k

  n_ij = n_i*n_j
  n_jk = n_j*n_k

  n_ij_loc = parallel_dim(n_ij,n_proc)
  n_jk_loc = parallel_dim(n_jk,n_proc)

  !----------------------------------
  ! Pointer alocation and definitions
  !
  allocate(ij(n_i,n_j))
  allocate(i_ij(n_ij))
  allocate(j_ij(n_ij))
  !
  p_ij = 0
  do i=1,n_i
     do j=1,n_j
        p_ij = p_ij+1
        ij(i,j) = p_ij
        i_ij(p_ij) = i
        j_ij(p_ij) = j
     enddo
  enddo
  !
  allocate(jk(n_j,n_k))
  allocate(i_rc(n_j,n_k))
  allocate(j_jk(n_jk))
  allocate(k_jk(n_jk))
  !
  p_jk = 0
  do j=1,n_j
     do k=1,n_k
        p_jk = p_jk+1
        jk(j,k) = p_jk
        j_jk(p_jk) = j
        k_jk(p_jk) = k
     enddo
  enddo
  !
  !-----------------------------------


  !-----------------------------------------
  ! Compute maximum length of send-pointer 
  ! array, s.
  ! 
  s = 0 
  do p_ij=1,n_ij,n_proc
     j = j_ij(p_ij)
     do k=1,n_k
        p_jk = jk(j,k)
        i_recv = modulo(p_jk-1,n_proc)
        s(i_recv) = s(i_recv)+1
     enddo ! k
  enddo ! p_ij
  !
  s_dim = maxval(s)
  !
  !------------------------------------------

  !-----------------------------------------
  ! Allocate send and receive package arrays
  ! (deallocated in TRANSP_CLEANUP)
  !
  allocate(q_send(s_dim,0:n_proc-1))
  allocate(q_recv(s_dim,0:n_proc-1))
  !
  !-----------------------------------------

  !--------------------------------------------------
  ! Build mapping arrays i_map, p_jk_loc_map
  !
  allocate(i_map(s_dim,0:n_proc-1))
  allocate(p_jk_loc_map(s_dim,0:n_proc-1))
  allocate(i_map2(s_dim,0:n_proc-1))
  allocate(p_jk_loc_map2(s_dim,0:n_proc-1))

  i_map = 0
  p_jk_loc_map = 0

  s = 0
  p_ij_loc = 0
  do p_ij=1+i_proc,n_ij,n_proc

     j = j_ij(p_ij)
     i = i_ij(p_ij)
     p_ij_loc = p_ij_loc+1

     do k=1,n_k

        p_jk      = jk(j,k)
        p_jk_loc  = (p_jk-1)/n_proc+1 
        i_recv    = modulo(p_jk-1,n_proc)
        s(i_recv) = s(i_recv)+1

        i_map(s(i_recv),i_recv) = i
        p_jk_loc_map(s(i_recv),i_recv) = p_jk_loc

     enddo ! k

  enddo ! p_ij
  !
  !--------------------------------------------------

  call MPI_ALLTOALL(i_map, &
       s_dim, &
       MPI_INTEGER, &
       i_map2, &
       s_dim, &
       MPI_INTEGER, &
       TRANSP_COMM, &
       i_err)

  call MPI_ALLTOALL(p_jk_loc_map, &
       s_dim, &
       MPI_INTEGER, &
       p_jk_loc_map2, &
       s_dim, &
       MPI_INTEGER, &
       TRANSP_COMM, &
       i_err)


  i_map = i_map2
  p_jk_loc_map = p_jk_loc_map2

  deallocate(i_map2)
  deallocate(p_jk_loc_map2)

  do p_ij=1+i_proc,n_ij,n_proc

     j = j_ij(p_ij)

     do k=1,n_k
        i_rc(j,k) = modulo(jk(j,k)-1,n_proc)
     enddo ! k

  enddo ! p_ij

end subroutine fTRANSP_INIT

