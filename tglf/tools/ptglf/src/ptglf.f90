program ptglf

  use mpi
  use ptglf_globals
  use tglf_interface

  implicit none

  integer :: p

  !---------------------------------------------------------------------
  ! Initialize MPI_COMM_WORLD communicator.
  !
  call MPI_INIT(i_err)
  call MPI_COMM_RANK(MPI_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc,i_err)
  !---------------------------------------------------------------------

  ntot = n_proc
  splitkey = i_proc
  color = 0
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,&
       color,&
       splitkey,&
       loc_comm,&
       i_err)

  ! Path is cwd:
  path= './'

  call tglf_init(path,loc_comm)
  call ptglf_init()
  
  if (i_proc == 0) print '(a,i5)','NTOT = ',ntot

  do p=1+i_proc,ntot,n_proc

     if (i_proc == 0) print '(i5,a,i5)',p,' - ',p+n_proc-1

     call tglf_run()
     
  enddo

  ! Collect all data 
  !call MPI_ALLREDUCE(indata_loc,indata,size(indata), &
  !     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,i_err)
  !call MPI_ALLREDUCE(outdata_j_loc,outdata_j,size(outdata_j), &
  !     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,i_err)
  
  !if (i_proc == 0) then

   !  open(unit=1,file='out.pneo.indata',status='replace')
   !  write(1,10) indata(:,:)
   !  close(1)

   !  open(unit=1,file='out.pneo.c_j',status='replace')
   !  write(1,10) outdata_j(:,:)
   !  close(1)
     
  !endif

  call MPI_finalize(i_err)

10 format(1pe17.10)

end program ptglf
