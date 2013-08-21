!-------------------------------------------------------------------------------
!
MODULE mpi12
!
!  jmp may 2006
!
   USE mpi
   IMPLICIT NONE
   INTEGER myid,numprocs,mpierr
   INTEGER, PARAMETER :: MAXPROC = 16
   INTEGER numj
   INTEGER js_mpi12(0:MAXPROC)
   INTEGER recvcounts(MAXPROC),displs(MAXPROC)
   REAL*8 mpi_etimes(10)
   INTEGER mpi_counts(10)
   INTEGER iloadbal
   REAL*8 glf_load(MAXPROC)

   CONTAINS

!-------------------------------------------------------------------------------
!  
   SUBROUTINE init_mpi12(numj_in)
  
   IMPLICIT NONE
   INTEGER numj_in
   INTEGER p
  
   numj = numj_in
  
   js_mpi12(0) = 1
   DO p=1,numprocs-1
  	  js_mpi12(p) = 1 + p*(numj/numprocs)
   END DO
   js_mpi12(numprocs) = numj	

   IF((iloadbal.eq.-1) .and. (numprocs.eq.4)) THEN
      !print *,'STATIC DECOMPOSITION'
      js_mpi12(1) = js_mpi12(1)+1+1!+1
      js_mpi12(2) = js_mpi12(2)+1+2!+1
      js_mpi12(3) = js_mpi12(3)+1+3!+1
   ENDIF

   mpi_etimes(1:10) = 0.0
  
   END SUBROUTINE

!-------------------------------------------------------------------------------
!  
   SUBROUTINE mpi12_g_send(vec,tag)
  
   IMPLICIT NONE
   REAL*8 vec(*) 
   INTEGER tag
   INTEGER :: request,istatus(MPI_STATUS_SIZE)
   INTEGER n,jin,dest,ierr
   
   jin = js_mpi12(myid)
   n = js_mpi12(myid+1)-js_mpi12(myid)
   if(myid.eq.numprocs-1) n = n + 1
   
   CALL MPI_ISEND(vec(jin), n, MPI_REAL8, 0, tag, MPI_COMM_WORLD, request, ierr)
   CALL MPI_WAIT(request, istatus, ierr)
   
   END SUBROUTINE  

!-------------------------------------------------------------------------------
!  
   SUBROUTINE mpi12_g_recv(vec,tag)
   
   IMPLICIT NONE
   REAL*8 vec(*) 
   INTEGER tag
   INTEGER, DIMENSION(MAXPROC) :: request
   INTEGER, DIMENSION(MPI_STATUS_SIZE,MAXPROC) :: istatus
   INTEGER n,jin,dest,ierr
   
   DO dest=1,numprocs-1
      jin = js_mpi12(dest)
   	  n = js_mpi12(dest+1)-js_mpi12(dest)
      if(dest.eq.numprocs-1) n = n + 1
      CALL MPI_IRECV(vec(jin), n, MPI_REAL8, dest, tag, MPI_COMM_WORLD, request(dest), ierr)
   END DO
   
   DO dest=1,numprocs-1
      CALL MPI_WAIT(request(dest), istatus(1,dest), ierr)
   END DO
   
   END SUBROUTINE  

!-------------------------------------------------------------------------------
!  
   SUBROUTINE mpi12_d_send(vec,tag)
   
   IMPLICIT NONE
   REAL*8 vec(*) 
   INTEGER tag
   INTEGER, DIMENSION(MAXPROC) :: request
   INTEGER, DIMENSION(MPI_STATUS_SIZE,MAXPROC) :: istatus
   INTEGER n,jin,jout,dest,ierr
   
   DO dest=1,numprocs-1
   	 jin = js_mpi12(dest)-1
   	 n = js_mpi12(dest+1)-js_mpi12(dest)+2
   	 CALL MPI_ISEND(vec(jin), n, MPI_REAL8, dest, tag, MPI_COMM_WORLD, request(dest), ierr)
   END DO
   
   DO dest=1,numprocs-1
      CALL MPI_WAIT(request(dest), istatus(1,dest), ierr)
   END DO
   
   END SUBROUTINE  

!-------------------------------------------------------------------------------
!  
   SUBROUTINE mpi12_d_recv(vec,tag)
   
   IMPLICIT NONE
   REAL*8 vec(*) 
   INTEGER tag
   INTEGER :: request,istatus(MPI_STATUS_SIZE)
   INTEGER n,jin,jout,dest,ierr
   
   jin = js_mpi12(myid)-1
   n = js_mpi12(myid+1)-js_mpi12(myid)+2
   
   CALL MPI_IRECV(vec(jin), n, MPI_REAL8, 0, tag, MPI_COMM_WORLD, request, ierr)
   CALL MPI_WAIT(request, istatus, ierr)
   
   END SUBROUTINE  

!-------------------------------------------------------------------------------
!  
   SUBROUTINE mpi12_send_int(vec,n,tag)
   
   IMPLICIT NONE
   INTEGER n
   INTEGER vec(n)
   INTEGER tag
   INTEGER, DIMENSION(MAXPROC) :: request
   INTEGER, DIMENSION(MPI_STATUS_SIZE,MAXPROC) :: istatus
   INTEGER dest,ierr
   
   DO dest=1,numprocs-1
   	 CALL MPI_ISEND(vec, n, MPI_INTEGER, dest, tag, MPI_COMM_WORLD, request(dest), ierr)
   END DO
   
   DO dest=1,numprocs-1
      CALL MPI_WAIT(request(dest), istatus(1,dest), ierr)
   END DO
   
   END SUBROUTINE
  
!-------------------------------------------------------------------------------
!  
   SUBROUTINE mpi12_recv_int(vec,n,tag)
   
   IMPLICIT NONE
   INTEGER n
   INTEGER vec(n)
   INTEGER tag
   INTEGER :: request,istatus(MPI_STATUS_SIZE)
   INTEGER dest,ierr
   
   CALL MPI_IRECV(vec, n, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, request, ierr)
   CALL MPI_WAIT(request, istatus, ierr)
    
   END SUBROUTINE
  
!-------------------------------------------------------------------------------
!  
   SUBROUTINE mpi12_send_real(vec,n,tag)
  
   IMPLICIT NONE
   INTEGER n
   REAL*8 vec(n)
   INTEGER tag
   INTEGER, DIMENSION(MAXPROC) :: request
   INTEGER, DIMENSION(MPI_STATUS_SIZE,MAXPROC) :: istatus
   INTEGER dest,ierr
   
   DO dest=1,numprocs-1
   	 CALL MPI_ISEND(vec, n, MPI_REAL8, dest, tag, MPI_COMM_WORLD, request(dest), ierr)
   END DO
   
   DO dest=1,numprocs-1
      CALL MPI_WAIT(request(dest), istatus(1,dest), ierr)
   END DO
  
   END SUBROUTINE
  
!-------------------------------------------------------------------------------
!  
   SUBROUTINE mpi12_recv_real(vec,n,tag)

   IMPLICIT NONE
   INTEGER n
   REAL*8 vec(n)
   INTEGER tag
   INTEGER :: request,istatus(MPI_STATUS_SIZE)
   INTEGER dest,ierr
   
   CALL MPI_IRECV(vec, n, MPI_REAL8, 0, tag, MPI_COMM_WORLD, request, ierr)
   CALL MPI_WAIT(request, istatus, ierr)
  
   END SUBROUTINE

!-------------------------------------------------------------------------------
!  
   SUBROUTINE dynamic_balance(glf_load_buf)
  
   IMPLICIT NONE
   REAL*8 glf_load_buf,glf_load_min,glf_load_max
   INTEGER ierr,p,p_min,p_max
   INTEGER, SAVE :: ncount = 0
   REAL*8 , SAVE :: glf_load_avg(10)
   INTEGER :: i
  
  !IF (ncount .eq. 0) glf_load_avg(1:10) = 0.0
  !
  !ncount = ncount + 1
  !IF(ncount .eq. 11) THEN
  !    ncount = 1
  !ENDIF
  !glf_load_avg(ncount) = glf_load_buf
  !
  !glf_load_buf = 0
  !DO i = 1, 10
  !  glf_load_buf = glf_load_buf + glf_load_avg(i)
  !END DO 
  !
  !IF(ncount .ne. 1) RETURN
  !
  !ncount = ncount + 1
  !IF(ncount .eq. 11) ncount = 1
  !IF(ncount .ne. 1) RETURN

   CALL MPI_GATHER(glf_load_buf,1,MPI_REAL8,glf_load,1,MPI_REAL8,0,MPI_COMM_WORLD, ierr)

   p_min = 1
   p_max = 1
   glf_load_min = glf_load(1)
   glf_load_max = glf_load(1)
      
   IF (myid .eq. 0) THEN 

     DO p = 2, numprocs

        IF(glf_load(p) > glf_load_max) THEN
      	   p_max = p
      	   glf_load_max = glf_load(p)
        END IF
      
        IF(glf_load(p) < glf_load_max) THEN
      	   p_min = p
      	   glf_load_min = glf_load(p)
        END IF

     END DO
  
     IF(p_min < p_max) THEN
    	
        DO p = p_min, p_max-1
      	   js_mpi12(p) = js_mpi12(p)+1
        END DO 	  
    
     ELSEIF(p_min > p_max) THEN
    	
        DO p = p_min-1, p_max, -1
      	   js_mpi12(p) = js_mpi12(p)-1
        END DO
     
     ENDIF  	 
  
   ENDIF  

   CALL MPI_BCAST(js_mpi12(0),MAXPROC+1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)

   END SUBROUTINE

END MODULE mpi12
