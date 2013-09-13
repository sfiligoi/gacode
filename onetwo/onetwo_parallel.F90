!-------------------------------------------------------------------------------
!
PROGRAM onetwo_parallel
!
!  jmp may 2006
!

   USE mpi12
   USE param, ONLY: kj
   
   IMPLICIT NONE
   INTEGER ierr
   CHARACTER*256 processor_name
   INTEGER namelen
   INTEGER iflag_mpistop(1)
   REAL time_e,time_s
   INTEGER setvbuf3f
   
   CALL MPI_Init(mpierr)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,mpierr)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,mpierr )
   CALL MPI_GET_PROCESSOR_NAME(processor_name,namelen,mpierr)
   
   time_s = mpi_wtime()
   
#ifdef OPT_JAGUAR
   ierr = setvbuf3f(6,0,1024)
#endif   

   PRINT *, 'Process ', myid, ierr,' of ', numprocs, ' is starting....'
   
   CALL init_mpi12(kj)
   
   IF (myid==0) then

      CALL onetwo

      print *,'final call at root start'
      iflag_mpistop(1) = 1
      CALL mpi12_send_int(iflag_mpistop,1,99)
      print *,'final call at root end'

   ELSE

      CALL onetwo_slave

   END IF
   
   time_e = mpi_wtime()
   
   PRINT *, 'Process ', myid, time_e-time_s
   PRINT *, 'MPI',mpi_etimes(1),mpi_etimes(2),mpi_etimes(3),mpi_etimes(4)
   PRINT *, 'GLF ', myid,mpi_etimes(4)
   PRINT *, 'COUNTS', myid, mpi_counts(1)

   CALL MPI_FINALIZE(ierr)

END PROGRAM

