!-----------------------------------------------------------------
! cgyro_mpi_grid.f90
!
! PURPOSE:
!  Subroutinized main cgyro program.  
!-----------------------------------------------------------------

subroutine cgyro_mpi_grid

  use timer_lib
  use mpi

  use cgyro_globals

  implicit none

  integer :: lc,lv,iv,ic,nc,nv
  integer :: ie,ix,is,ir,it
  complex, dimension(:,:,:), allocatable :: h
  complex, dimension(:,:,:), allocatable :: hc
  real :: error

  integer :: n_n
  integer :: splitkey
  integer, external :: parallel_dim

  n_n = 1

  nv = n_energy*n_xi*n_species
  nc = n_radial*n_theta

  allocate(ie_v(nv))
  allocate(ix_v(nv))
  allocate(is_v(nv))

  allocate(ir_c(nc))
  allocate(it_c(nc))

  !------------------------------------------------
  ! Check for grid validity
  !
  if (modulo(n_proc,n_n) /= 0) then
     call catch_error('ERROR: (CGYRO) bad processor count.')
  endif
  !------------------------------------------------

  !-------------------------------
  ! Assign subgroup dimensions:
  !
  n_proc_2 = n_n
  n_proc_1 = n_proc/n_n
  !-------------------------------

  !------------------------------------------------
  ! Local group indices:
  !
  i_group_1 = i_proc/n_proc_1
  i_group_2 = modulo(i_proc,n_proc_1)
  !------------------------------------------------

  !-----------------------------------------------------------
  ! Split up GYRO_COMM_WORLD into groups and adjoint:
  !
  !             NEW_COMM_1  and  NEW_COMM_2
  !
  splitkey = i_proc
  call MPI_COMM_SPLIT(CGYRO_COMM_WORLD,&
       i_group_1,& 
       splitkey,&
       NEW_COMM_1, &
       i_err)
  if (i_err /= 0) then
     call CATCH_ERROR('ERROR: (CGYRO) NEW_COMM_1 not created')
  endif

  ! Local adjoint Group number

  call MPI_COMM_SPLIT(CGYRO_COMM_WORLD,&
       i_group_2,& 
       splitkey,&
       NEW_COMM_2, &
       i_err)
  if (i_err /= 0) then
     call CATCH_ERROR('ERROR: (CGYRO) NEW_COMM_2 not created')
  endif
  !
  call MPI_COMM_RANK(NEW_COMM_1,i_proc_1,i_err)
  call MPI_COMM_RANK(NEW_COMM_2,i_proc_2,i_err)
  !
  !-----------------------------------------------------------

  ! Velocity pointers
  lv = 0
  do ie=1,n_energy
     do ix=1,n_xi
        do is=1,n_species
           lv = lv+1
           ie_v(lv) = ie
           ix_v(lv) = ix
           is_v(lv) = is
        enddo
     enddo
  enddo

  ! Configuration pointers
  lc = 0
  do ir=1,n_radial
     do it=1,n_theta
        lc = lc+1
        ir_c(lc) = ir
        it_c(lc) = it
     enddo
  enddo

  ! Parallelization dimensions
  !
  ! h(1,ic,iv_loc) -> hc(1,iv,ic_loc)
  nv_loc = parallel_dim(nv,n_proc_1)
  allocate(h(1,nc,nv_loc))
  nc_loc = parallel_dim(nc,n_proc_1)
  allocate(hc(1,nv,nc_loc))

  ! Distributed-in-v loop
  !
  !iv_loc = 0
  !do iv=1+i_proc,nv,n_proc_1
  !   iv_loc = iv_loc+1
  !   do ic=1,nc
  !      h(1,ic,iv_loc) = &
  !   enddo
  !enddo

  !call rTRANSP_INIT(1,nv,nc,1,NEW_COMM_1)
  !call rTRANSP_DO(h,hc)
  !call rTRANSP_CLEANUP

  ! Distributed-in-c loop
  !
  !ic_loc = 0
  !do ic=1+i_proc,nc,n_proc_1
  !   ic_loc = ic_loc+1
  !   do iv=1,nv
  !       hc(1,iv,ic_loc) = &
  !   enddo
  !enddo

end subroutine cgyro_mpi_grid
