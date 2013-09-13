  MODULE MPI_data
    USE nrtype,                ONLY : I4B,I2B,DP


#if defined  (USEMPI)
    USE MPI
    IMPLICIT NONE
    INTEGER(I4B) myid,numprocs,mpiierr,master,jac_tree,jac_comm,  &
                 comm_groups,jac_grp_no,mpipass,parallel_model,   &
                 loc_nump,loc_rank,loc_comm,local_id
    INTEGER(I4B) mpiversion, mpisubversion,mstr_window,sizeof_dp
    INTEGER(I4b) mpi_status(MPI_STATUS_SIZE)
    INTEGER(I4B),ALLOCATABLE,DIMENSION(:)   :: p_jcg,proc_map
    REAL(DP),ALLOCATABLE,DIMENSION(:)       :: proc_time,comm_time
    REAL(DP) mpi_start_time,mpi_end_time
    CHARACTER*(MPI_MAX_PROCESSOR_NAME),ALLOCATABLE,DIMENSION(:)   &
                                            :: proc_name
    CHARACTER*(MPI_MAX_PROCESSOR_NAME)      :: local_proc_name,lpc_w
    LOGICAL,SAVE                            :: PROC_DAT_LOADED,mstr_active,&
                                               initialized_mpi
#ifdef NFREYA
    CHARACTER*1,ALLOCATABLE,DIMENSION(:) ::      mpi_buffer ! used for data passing in P_Nfreya
    INTEGER(I4B) buff_size
#endif

    DATA PROC_DAT_LOADED /.FALSE./
    DATA initialized_mpi /.FALSE./



    CONTAINS
      SUBROUTINE MPI_profile(iounit,include_tglf)
        USE common_constants,             ONLY :  zeroc
        USE MPI
        USE solcon_gcnmp,                 ONLY : t_slave_compt_alltot

        IMPLICIT NONE
        REAL(DP) dumy
        INTEGER(I2B)j,iounit,lstr,mlstr,include_tglf
        INTEGER status(MPI_STATUS_SIZE)

        IF(myid == master .AND. .NOT. ALLOCATED(proc_name))THEN                 
           ALLOCATE(proc_name(0:numprocs-1))
           proc_name(0) = ADJUSTL(local_proc_name)
        ENDIF


        loc_nump = numprocs 
        IF( .NOT. PROC_DAT_LOADED) THEN

           PROC_DAT_LOADED = .TRUE.
           DO j=0,loc_nump-1
              dumy = proc_time(j)
              proc_time(j) = zeroc
              CALL MPI_REDUCE(dumy,proc_time(j),1,MPI_DOUBLE_PRECISION,   &
                   MPI_SUM,master,MPI_COMM_WORLD,mpiierr) 
              dumy = comm_time(j)
              comm_time(j) = zeroc
              CALL MPI_REDUCE(dumy,comm_time(j),1,MPI_DOUBLE_PRECISION,   &
                   MPI_SUM,master,MPI_COMM_WORLD,mpiierr) 
           ENDDO

           IF(myid .NE. master)                                           &
                CALL MPI_SEND(local_proc_name,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,master,    &
                myid,MPI_COMM_WORLD,mpiierr)

           IF(myid == master)THEN
              DO j =1,loc_nump-1
                 CALL MPI_RECV(local_proc_name,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,MPI_ANY_SOURCE, &
                      MPI_ANY_TAG,MPI_COMM_WORLD,status,mpiierr)
                 proc_name(status(MPI_SOURCE)) = ADJUSTL(local_proc_name)
              ENDDO
           ENDIF
              
        ENDIF



        IF(myid == master)THEN
           WRITE(iounit,3)
           mlstr = 0
           DO j=0,loc_nump-1
              lstr = LEN_TRIM(proc_name(j))
              mlstr = MAX(lstr,mlstr)
           ENDDO

           DO j=0,loc_nump-1
              IF(include_tglf == 1)proc_time(j) = t_slave_compt_alltot(j)
              WRITE(iounit,4)j,proc_name(j)(1:mlstr),proc_time(j),comm_time(j)
           ENDDO
        ENDIF

3       FORMAT(2x,'process #:   process name:  sec in diffuse:   sec in collective communication:')
4       FORMAT(2x,i3,12x,a,9x,1pe12.4,14x,1pe12.4)

      END SUBROUTINE MPI_profile



#else
    IMPLICIT NONE
    INTEGER(I4B) myid,numprocs,mpiierr,master,jac_tree,jac_comm,  &
                 comm_groups,jac_grp_no,mpipass,parallel_model,&
                 loc_nump,loc_rank,loc_comm,local_id
    INTEGER(I4B),ALLOCATABLE,DIMENSION(:)    :: p_jcg,proc_map
    INTEGER(I4B) mpiversion, mpisubversion, mstr_window, mpi_status(1)
    REAL(DP),ALLOCATABLE,DIMENSION(:)        :: proc_time,comm_time
    REAL(DP) mpi_start_time,mpi_end_time
!    CHARACTER*(127),ALLOCTATBLE,DIMENSION(:) :: proc_name
    CHARACTER*(127)                          :: local_proc_name
    LOGICAL initialized_mpi, mstr_active
    DATA initialized_mpi /.FALSE./
#endif


  END MODULE MPI_DATA
 
