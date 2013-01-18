!------------------------------------------------------
! rTRANSP_INIT.f90
!
! PURPOSE:
!  Initialization stage for calls rTRANSP.
!
! NOTES:
!  Communication routine originally used in old GYRO.
!------------------------------------------------------

subroutine rTRANSP_INIT(n0_i,n0_j,n0_k,n0_m,COMM)

  use mpi
  use rTRANSP_GLOBALS

  implicit none

  integer, intent(in) :: n0_i
  integer, intent(in) :: n0_j
  integer, intent(in) :: n0_k
  integer, intent(in) :: n0_m
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
  n_m = n0_m

  n_ij = n_i*n_j
  n_ki = n_k*n_i

  n_ij_loc = parallel_dim(n_ij,n_proc)
  n_ki_loc = parallel_dim(n_ki,n_proc)

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
  allocate(ki(n_k,n_i))
  allocate(i_rc(n_k,n_i))
  allocate(k_ki(n_ki))
  allocate(i_ki(n_ki))
  !
  p_ki = 0
  do k=1,n_k
     do i=1,n_i
        p_ki = p_ki+1
        ki(k,i) = p_ki
        k_ki(p_ki) = k
        i_ki(p_ki) = i
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
     i = i_ij(p_ij)
     do k=1,n_k
        p_ki = ki(k,i)
        i_recv = modulo(p_ki-1,n_proc)
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
  allocate(q_send(n_m,s_dim,0:n_proc-1))
  allocate(q_recv(n_m,s_dim,0:n_proc-1))
  !
  !-----------------------------------------


  !--------------------------------------------------
  ! Build mapping arrays i_map, p_jk_loc_map
  !
  allocate(j_map(s_dim,0:n_proc-1))
  allocate(p_ki_loc_map(s_dim,0:n_proc-1))
  allocate(j_map2(s_dim,0:n_proc-1))
  allocate(p_ki_loc_map2(s_dim,0:n_proc-1))

  j_map = 0
  p_ki_loc_map = 0

  s = 0
  p_ij_loc = 0
  do p_ij=1+i_proc,n_ij,n_proc

     j = j_ij(p_ij)
     i = i_ij(p_ij)
     p_ij_loc = p_ij_loc+1

     do k=1,n_k

        p_ki      = ki(k,i)
        p_ki_loc  = (p_ki-1)/n_proc+1 
        i_recv    = modulo(p_ki-1,n_proc)
        s(i_recv) = s(i_recv)+1

        j_map(s(i_recv),i_recv) = j
        p_ki_loc_map(s(i_recv),i_recv) = p_ki_loc

     enddo ! k

  enddo ! p_ij
  !
  !--------------------------------------------------

  call MPI_ALLTOALL(j_map, &
       s_dim, &
       MPI_INTEGER, &
       j_map2, &
       s_dim, &
       MPI_INTEGER, &
       TRANSP_COMM, &
       i_err)

  call MPI_ALLTOALL(p_ki_loc_map, &
       s_dim, &
       MPI_INTEGER, &
       p_ki_loc_map2, &
       s_dim, &
       MPI_INTEGER, &
       TRANSP_COMM, &
       i_err)


  j_map = j_map2
  p_ki_loc_map = p_ki_loc_map2

  deallocate(j_map2)
  deallocate(p_ki_loc_map2)

  do p_ij=1+i_proc,n_ij,n_proc

     i = i_ij(p_ij)

     do k=1,n_k
        i_rc(k,i) = modulo(ki(k,i)-1,n_proc)
     enddo ! k

  enddo ! p_ij


end subroutine rTRANSP_INIT

