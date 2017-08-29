program pneo

  use mpi
  use pneo_globals
  use neo_interface
  use EXPRO_interface

  implicit none

  integer :: i,j,p

  !---------------------------------------------------------------------
  ! Initialize MPI_COMM_WORLD communicator.
  !
  call MPI_INIT(i_err)
  call MPI_COMM_RANK(MPI_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc,i_err)
  !---------------------------------------------------------------------

  ! Path is cwd:
  path= './'

  ! Read vgen control parameters

  call neo_init_serial(path)
  neo_silent_flag_in = 1

  ! pointers
  ni=8
  nj=8

  ntot = ni*nj
  allocate(ic(ntot))
  allocate(jc(ntot))
  allocate(data_vec(3,ntot))
  allocate(data_tot(3,ntot))

  p = 0
  do i=1,ni
     do j=1,nj
        p = p+1
        ic(p) = i
        jc(p) = j
     enddo
  enddo

  do p=1+i_proc,ntot,n_proc

     i = ic(p) 
     j = jc(p)

     neo_delta_in = i/10.0
     neo_q_in     = 1.0+j*0.5

     call neo_run()

     data_vec(1,p) = neo_delta_in
     data_vec(2,p) = neo_q_in
     data_vec(3,p) = neo_jpar_dke_out

  enddo

  ! Collect all data 
  call MPI_ALLREDUCE(data_vec,data_tot,size(data_tot), &
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,i_err)

  if (i_proc == 0) then
     print *,data_tot(3,:)
  endif

  call MPI_finalize(i_err)

end program pneo
