  MODULE plasma_shape
!---------------------------------------------------------
! -- put subroutines for shape related calcs here
!---------------------------------------------------------
 
    USE nrtype,              ONLY : I4B,DP

    USE vector_class,        ONLY : Vector_loc_min,Vector_loc_max

    USE plasma_properties,   ONLY : mhd_dat,dischg

    USE common_constants,    ONLY : pisq,zeroc,izero

    USE error_handler,       ONLY : lerrno,terminate,iomaxerr

    USE MPI_data,            ONLY : myid,master

#ifdef GCNMP
    USE io_gcnmp,            ONLY : ncrt,nlog 
#elif defined  NFREYA
    USE io_gcnmp,            ONLY : ncrt,nlog 
#else
    USE  io,                 ONLY : ncrt,nlog => nlog_gcnmp
#endif

    IMPLICIT NONE


  CONTAINS
    SUBROUTINE process_shape
     REAL(DP) rcross(5)
     REAL(DP) dzm,dzp,rmax,rmin
     INTEGER(I4B) k,j


    k = 0
    rcross(:) = zeroc
    DO j=2,dischg%nplasbdry-1
       ! search for change in sign
      dzm = dischg%zplasbdry%data(j)   - dischg%zma
      dzp = dischg%zplasbdry%data(j+1) - dischg%zma
      IF(dzm*dzp .LE. zeroc)THEN
         k = k+1
         rcross(k) = (dischg%rplasbdry%data(j+1)  +  dischg%rplasbdry%data(j))/2._DP
      ENDIF
    ENDDO
     rmin = MINVAL(rcross,MASK = rcross .GT. zeroc)
     rmax = MAXVAL(rcross)
     dischg%rminor = (rmax-rmin)/2._DP
     IF(k .GT. 2)THEN
        IF(myid == master)THEN
           PRINT *,' ERROR: detected more than 2 crossings of zma '
           PRINT *,'        during contour search for rmin,rmax '
           PRINT *,'        R Contour values are =',rcross(:)
           PRINT *,'        values of rplasbdry,zplasbdry in state file not consistent'
        ENDIF
        lerrno  = iomaxerr + 156
        CALL terminate( lerrno,nlog)
     ENDIF

    RETURN
    END     SUBROUTINE process_shape


  END MODULE plasma_shape
