MODULE echdat_module
 
!-------------------------------------------------------------------------
!   This module is created by the cdf2fortran code 
!   It should be regenerated each time the netcdf file changes 
!   The purpose of this module is to enable Onetwo to use
!   time dependent ech dat for gyrotrons given in the form
!   of a netcdf file .
!-----------------------------------------------------------------HSJ-----
      USE nrtype,                            ONLY : DP,SP,I4B,I2B

      USE typeSizes                  !part of netcdf 
      USE netcdf                     !ditto
      ! define dimensions as integer. Values to be read from netcdf file later.

       PRIVATE

       PUBLIC set_irf,nsys,echdat_read,get_ech_params
       PUBLIC gauszone,nharm,netcdfdat,modelc,model_globl

      INTEGER*4     nvars, ndims          ! number of dims and number of vars
 
      INTEGER      nsys
      INTEGER      ntimes
      INTEGER      namelength
      INTEGER      stringlength
      INTEGER ncid,nctype,nvdim,nvatts,model_globl 
      INTEGER lenstr,ndsize
      INTEGER*4  ::  shot
 
      INTEGER*4      rcode                   ! error code
      INTEGER*4     recdim                    ! record dimension
      INTEGER       dimid
 
      ! define variables as allocatable
      INTEGER*4,ALLOCATABLE,DIMENSION(:) :: dimsiz,start,count,vdims
      CHARACTER*128,ALLOCATABLE,DIMENSION(:) :: dimnam
 

 
     REAL(DP)  ::  dt_fast
     REAL(DP)  ::  t_avg
     CHARACTER,    ALLOCATABLE,DIMENSION(:)    ::  date
     CHARACTER,    ALLOCATABLE,DIMENSION(:)    ::  user
     REAL(DP),     ALLOCATABLE,DIMENSION(:)    ::  FREQ
     REAL(DP),     ALLOCATABLE,DIMENSION(:)    ::  HLWEC
     REAL(DP),     ALLOCATABLE,DIMENSION(:)    ::  time
     REAL(DP),     ALLOCATABLE,DIMENSION(:)    ::  stime
     CHARACTER,    ALLOCATABLE,DIMENSION(:,:)  ::  gyrotron
     REAL(DP),     ALLOCATABLE,DIMENSION(:,:)  ::  RFPOW
     REAL(DP),     ALLOCATABLE,DIMENSION(:,:)  ::  WRFO
     REAL(DP),     ALLOCATABLE,DIMENSION(:,:)  ::  XEC
     REAL(DP),     ALLOCATABLE,DIMENSION(:,:)  ::  ZEC
     REAL(DP),     ALLOCATABLE,DIMENSION(:,:)  ::  PHAIEC
     REAL(DP),     ALLOCATABLE,DIMENSION(:,:)  ::  THETEC
     REAL(DP),     ALLOCATABLE,DIMENSION(:)    ::  RATWEC
     REAL(DP),     ALLOCATABLE,DIMENSION(:)    ::  IRFCUR


     INTEGER(nf90_int), ALLOCATABLE, DIMENSION(:)   ::  igafit
     INTEGER(I4B), ALLOCATABLE, DIMENSION(:)   ::  gauszone
     INTEGER(I4B), ALLOCATABLE, DIMENSION(:)   ::  nray
     INTEGER(I4B), ALLOCATABLE, DIMENSION(:)   ::  idamp
     INTEGER(I4B), ALLOCATABLE, DIMENSION(:)   ::  nharm
     INTEGER(I4B), ALLOCATABLE, DIMENSION(:)   ::  modelc
     INTEGER(I4B), ALLOCATABLE, DIMENSION(:)   ::  netcdfdat
     
     
     
     

   CONTAINS 

   SUBROUTINE echdat_read(input_file)  
!-----------------------------------------------------------------------
!
!  This file is a fortran template file designed to read the netCDF file 
!  into memory.
!-----------------------------------------------------------HSJ---------

  USE param,                                  ONLY : krf

  USE io,                                     ONLY : nout,ncrt,nitre,    &
                                                     io_toray_hist

  USE rf,                                     ONLY : nmodel,&
                                                     rfmode,&
                                                     nray_12   => nray,  &
                                                     idamp_12  => idamp, &
                                                     ratwec_12 => ratwec,&
                                                     no_rf

                                                 
     IMPLICIT NONE

      INTEGER i,j,newd


      INTEGER(I4B)                 &
           id_date,                &
           id_user,                &
           id_shot,                &
           id_nsys,                &
           id_freq,                &
           id_hlwec,               &
           id_ratwec,              &
           id_time,                &
           id_gyrotron,            &
           id_rfpow,               &
           id_wrfo,                &
           id_xec,                 &
           id_zec,                 &
           id_phaiec,              &
           id_thetec,              &
           id_igafit,              &
           id_idamp,               &
           id_gauszone,            &
           id_nray,                &
           id_nharm,               &
           id_modelc,              &
           id_irfcur,              &
           id_dt_fast,             &
           id_t_avg,               &
           id_netcdfdat


      CHARACTER*(*) input_file
      CHARACTER*1024 strbuf, strbuf2            ! string buffer for var
                                                !  and attr names

!     Open netCDF file:
      rcode =  NF90_OPEN (input_file,NF90_NOWRITE, ncid)
      IF(rcode .NE. NF90_NOERR)THEN
         WRITE(ncrt,1000 )input_file
         WRITE(nitre,1000)input_file
1000      FORMAT("  ERROR, netcdf file not found: ",a,/, &
              "  Onetwo must stop ")
         CALL STOP('ech input netcdf file not found',1)
      ENDIF
!     Inquire about the number of dims, vars, and atts.
      CALL ncinq(ncid,ndims,nvars,nvatts,recdim,rcode)

      ALLOCATE(dimsiz(ndims))
      ALLOCATE(vdims(ndims))
      ALLOCATE(start(ndims))
      ALLOCATE(COUNT(ndims))
      ALLOCATE(dimnam(ndims))


      !     Store the dimension names and sizes in arrays.
      DO j=1,ndims !number of dimensions returned from ncinq
        dimid=j
        CALL ncdinq(ncid,dimid,dimnam(j),dimsiz(j),rcode)
        IF (dimnam(j) == "nsys") nsys = dimsiz(j)
        IF (dimnam(j) == "ntimes") ntimes = dimsiz(j)
        IF (dimnam(j) == "namelength") namelength = dimsiz(j)
        IF (dimnam(j) == "stringlength") stringlength = dimsiz(j)
      END DO

     ALLOCATE(date(stringlength))
     ALLOCATE(user(stringlength))
     ALLOCATE(FREQ(nsys))
     ALLOCATE(HLWEC(nsys))
     ALLOCATE(RATWEC(nsys))
     ALLOCATE(time(ntimes))
     ALLOCATE(gyrotron(nsys,namelength)) 

     ALLOCATE(RFPOW(ntimes,nsys))
     ALLOCATE(WRFO(ntimes,nsys))
     ALLOCATE(XEC(ntimes,nsys))
     ALLOCATE(ZEC(ntimes,nsys))
     ALLOCATE(PHAIEC(ntimes,nsys))
     ALLOCATE(THETEC(ntimes,nsys))
     ALLOCATE(igafit(nsys))
     ALLOCATE(idamp(nsys))
     ALLOCATE(gauszone(nsys))
     ALLOCATE(nray(nsys))
     ALLOCATE(nharm(nsys))
     ALLOCATE(modelc(nsys))
     ALLOCATE(irfcur(nsys))
     ALLOCATE(netcdfdat(nsys))


      igafit(1:nsys) = -1 ; nray(1:nsys) = -1  ; idamp (1:nsys) = -1 

      CALL ech_err( nf90_inq_varid(ncid,"date", id_date),-id_date)
      CALL ech_err( nf90_get_var(ncid,id_date,date),-id_date)

      CALL ech_err( nf90_inq_varid(ncid,"user", id_user),-id_user)
      CALL ech_err( nf90_get_var(ncid,id_user,user),-id_user)

      CALL ech_err( nf90_inq_varid(ncid,"shot", id_shot),-id_shot)
      CALL ech_err( nf90_get_var(ncid,id_shot,shot),-id_shot)

      CALL ech_err( nf90_inq_varid(ncid,"nsys", id_nsys),-id_nsys)
      CALL ech_err( nf90_get_var(ncid,id_nsys,nsys),-id_nsys)

      CALL ech_err( nf90_inq_varid(ncid,"dt_fast", id_dt_fast),-id_dt_fast)
      CALL ech_err( nf90_get_var(ncid,id_dt_fast,dt_fast),-id_dt_fast)

      CALL ech_err( nf90_inq_varid(ncid,"t_avg", id_t_avg),-id_t_avg)
      CALL ech_err( nf90_get_var(ncid,id_t_avg,t_avg),-id_t_avg)

      CALL ech_err( nf90_inq_varid(ncid,"gyrotron", id_gyrotron),-id_gyrotron)
      CALL ech_err( nf90_get_var(ncid,id_gyrotron,gyrotron),-id_gyrotron)

      CALL ech_err( nf90_inq_varid(ncid,"FREQ", id_freq),-id_freq)
      CALL ech_err( nf90_get_var(ncid,id_freq,FREQ),-id_freq)

      CALL ech_err( nf90_inq_varid(ncid,"HLWEC", id_hlwec),-id_hlwec)
      CALL ech_err( nf90_get_var(ncid,id_hlwec,HLWEC),-id_hlwec)

      CALL ech_err( nf90_inq_varid(ncid,"time", id_time),-id_time)
      CALL ech_err( nf90_get_var(ncid,id_time,time),-id_time)

      CALL ech_err( nf90_inq_varid(ncid,"RFPOW", id_rfpow),-id_rfpow)
      CALL ech_err( nf90_get_var(ncid,id_rfpow,RFPOW),-id_rfpow)

      CALL ech_err( nf90_inq_varid(ncid,"WRFO", id_wrfo),-id_wrfo)
      CALL ech_err( nf90_get_var(ncid,id_wrfo,WRFO),-id_wrfo)


      CALL ech_err( nf90_inq_varid(ncid,"XEC", id_xec),-id_xec)
      CALL ech_err( nf90_get_var(ncid,id_xec,XEC),-id_xec)

      CALL ech_err( nf90_inq_varid(ncid,"ZEC", id_zec),-id_zec)
      CALL ech_err( nf90_get_var(ncid,id_zec,ZEC),-id_zec)

      CALL ech_err( nf90_inq_varid(ncid,"PHAIEC", id_phaiec),-id_phaiec)
      CALL ech_err( nf90_get_var(ncid,id_phaiec,PHAIEC),-id_phaiec)

      CALL ech_err( nf90_inq_varid(ncid,"THETEC", id_thetec),-id_thetec)
      CALL ech_err( nf90_get_var(ncid,id_thetec,THETEC),-id_thetec)

      CALL ech_err( nf90_inq_varid(ncid,"IGAFIT", id_igafit),-id_igafit)
      CALL ech_err(nf90_get_var(ncid,id_igafit,IGAFIT),-id_igafit)

      CALL ech_err( nf90_inq_varid(ncid,"GAUSZONE", id_gauszone),-id_gauszone)
      CALL ech_err(nf90_get_var(ncid,id_gauszone,GAUSZONE),-id_gauszone)

      CALL ech_err( nf90_inq_varid(ncid,"NRAY", id_nray),-id_nray)
      CALL ech_err(nf90_get_var(ncid,id_nray,NRAY),-id_nray)

      CALL ech_err( nf90_inq_varid(ncid,"RATWEC", id_ratwec),-id_ratwec)
      CALL ech_err(nf90_get_var(ncid,id_ratwec,RATWEC),-id_ratwec)

      CALL ech_err( nf90_inq_varid(ncid,"IDAMP", id_idamp),-id_idamp)
      CALL ech_err(nf90_get_var(ncid,id_idamp,IDAMP),-id_idamp)

      CALL ech_err( nf90_inq_varid(ncid,"NHARM", id_nharm),-id_nharm)
      CALL ech_err(nf90_get_var(ncid,id_nharm,NHARM),-id_nharm)

      CALL ech_err( nf90_inq_varid(ncid,"MODELC", id_modelc),-id_modelc)
      CALL ech_err(nf90_get_var(ncid,id_modelc,MODELC),-id_modelc)

      CALL ech_err( nf90_inq_varid(ncid,"IRFCUR", id_irfcur),-id_irfcur)
      CALL ech_err(nf90_get_var(ncid,id_irfcur,IRFCUR),-id_irfcur)

      CALL ech_err( nf90_inq_varid(ncid,"NETCDFDAT", id_netcdfdat),-id_netcdfdat)
      CALL ech_err(nf90_get_var(ncid,id_netcdfdat,NETCDFDAT),-id_netcdfdat)
 

     zec(:,:) = 100.0*zec(:,:)  ! need these in cm
     xec(:,:) = 100.0*xec(:,:)



    !disable any ech input present in inone file.
    !The hard and fast rule is that either ech is input 
    ! in inone or it is given in the ech input netcdf file.

     DO j=1,krf
        IF(rfmode(j) == 'ech')THEN
           WRITE(ncrt,10) input_file   ! screen
           WRITE(nout,10) input_file   ! outone
           WRITE(nitre,10) input_file  ! runlog
           CALL STOP ('subroutine echdat_read', 1)
        ENDIF
     ENDDO
 10  FORMAT(2x,'ERROR: ech input can not be present in both inone',/,  &
            2x,'and file ',a)

     nmodel =0
     DO j =1,krf
       IF(rfmode(j) ==  no_rf)EXIT
       nmodel = nmodel + 1   ! this gives the model index for use in sub get_ech_params
     ENDDO

     i = nmodel
     DO j=1,nsys
        i  = i  +1
        IF(i .GT. krf)THEN
           WRITE(ncrt,15)   ! screen
           WRITE(nout,15)   ! outone
           WRITE(nitre,15)  ! runlog
           CALL STOP ('subroutine echdat_read', 2)
        ENDIF
        rfmode(i)    = 'ech'          ! not present in netcdf file
        nray_12(i)   =  nray(i)       ! current set default
        idamp_12(i)  =  idamp(i)      ! current set default
        ratwec_12(i) =  ratwec(i)
     ENDDO
  15           FORMAT(2x,' error ,krf too small to store ech netcdf input')


      ! create toray.in history  file
      ! reults will be written to this file in sub get_ech_params
      open  (unit = io_toray_hist, file = "toray.in_hist", status = 'UNKNOWN')

     RETURN

     END SUBROUTINE echdat_read

   
    SUBROUTINE get_ech_params(time_w,model)
!--------------------------------------------------------------------------
! -- this subbroutine works with the data read in 
! -- echdat.nc type file ( read in sub echdat_read).
! -- for each gyrotron we scan the data and construct input for toray
! -- NOTE:  wrt_toray_in2 requires the model index to select the 
! -- appropriate gyrotron. This inde is not known until get_ech_params
! -- is called from sub source.
! -- Here we capture this value and make it available by storing it
! -- in this module.
!--------------------------------------------------------------HSJ---------

     USE  nrtype,                            ONLY : DP,I4B,SP
 

     USE rf,                                 ONLY :                         &
                                                 rfpow_12  => rfpow,        &
                                                 freq_12   => freq,         &
                                                 xec_12    => xec,          &
                                                 zec_12    => zec,          &
                                                 phaiec_12 => phaiec,       &
                                                 thetec_12 => thetec,       &
                                                 wrfo_12   => wrfo,         &
                                                 hlwec_12  => hlwec,        &
                                                 rfmode_12 => rfmode,       &
                                                 nray_12   => nray,         &
                                                 nmodel

     USE io,                                 ONLY : nout,ncrt,nitre,        &
                                                    io_toray_hist

     IMPLICIT NONE
     REAL(DP) time_w
     REAL(DP) delp,delt
     INTEGER ndtp1,ndtlo,model,j
     REAL(DP) time_wms
          time_wms = SNGL(time_w)*1000._SP
          model_globl = model
 
         IF(rfmode_12(model) .NE.  'ech')THEN
           WRITE(ncrt,15)   ! screen
           WRITE(nout,15)   ! outone
           WRITE(nitre,15)  ! runlog
           CALL STOP ('subroutine get_ech_params', 2)
         ENDIF 
15       FORMAT(2X,'ERROR,mismatch in rf models') 

         ! find bracketing terms for time_wms in array time
         CALL tableintrp_sp (time, ntimes, time_wms, ndtlo)

         ndtp1              = ndtlo + 1
         delt               = time(ndtp1)-time(ndtlo)
         j                  = model - nmodel

         !delp               = rfpow(j,ndtp1)-rfpow(j,ndtlo)
         !rfpow_12(model)    = rfpow(j,ndtlo)+ (delp/delt)*(time_wms -time(ndtlo)) 
         delp               = rfpow(ndtp1,j)-rfpow(ndtlo,j)
         rfpow_12(model)    = rfpow(ndtlo,j)+ (delp/delt)*(time_wms -time(ndtlo)) 


         CALL tableintrp_sp (time, ntimes, time_wms, ndtlo)
         ndtp1              = ndtlo + 1
         delt               = time(ndtp1)-time(ndtlo)


         delp               = thetec(ndtp1,j)-thetec(ndtlo,j)
         thetec_12(model)   = thetec(ndtlo,j)+ (delp/delt)*(time_wms -time(ndtlo)) 

         delp               = phaiec(ndtp1,j)-phaiec(ndtlo,j)
         phaiec_12(model)   = phaiec(ndtlo,j)+ (delp/delt)*(time_wms -time(ndtlo))

         delp               = xec(ndtp1,j)-xec(ndtlo,j)
         xec_12(model)      = xec(ndtlo,j)+ (delp/delt)*(time_wms -time(ndtlo)) 

         delp               = zec(ndtp1,j)-zec(ndtlo,j)
         zec_12(model)      = zec(ndtlo,j)+ (delp/delt)*(time_wms -time(ndtlo))

         delp               = wrfo(ndtp1,j)-wrfo(ndtlo,j)
         wrfo_12(model)     = wrfo(ndtlo,j)+ (delp/delt)*(time_wms -time(ndtlo))



         freq_12(model)     = freq(j)
         hlwec_12(model)    = hlwec(j)
         nray_12(model)     = nray(j)

   
         write(io_toray_hist,FMT='("model,time = ",i5,x,1pe12.4)'),model, time_wms
         write(io_toray_hist,FMT='("rfpow_12(model)   = ",1pe12.4)')rfpow_12(model)
         write(io_toray_hist,FMT='("thetec_12(model)  = ",1pe12.4)')thetec_12(model)
         write(io_toray_hist,FMT='("phaiec_12(model ) = ",1pe12.4)')phaiec_12(model)
         write(io_toray_hist,FMT='("xec_12(model)     = ",1pe12.4)')xec_12(model)
         write(io_toray_hist,FMT='("zec_12(model)     = ",1pe12.4)')zec_12(model)
         write(io_toray_hist,FMT='("wrfo_12(model)    = ",1pe12.4)')wrfo_12(model)
         write(io_toray_hist,FMT='("nray_12(model)    = ",i5)')nray_12(model)

   RETURN


   END SUBROUTINE get_ech_params

  SUBROUTINE set_irf
!----------------------------------------------------------------------
! --   adds nsys ech cases. 
! --   Note that irf(k) is set to 3, indicating that gyrotron is always on
! --   it may be off , indicated by power level rfpow
!----------------------------------------------------------------------
     USE  nrtype,                            ONLY : DP,I4B

     USE rf,                                 ONLY : krf,irf,irfmodel,rfmode, &
                                                    irfcur_12 =>irfcur
     IMPLICIT NONE
     INTEGER(I4B) j,k
     k = 0
     DO j=1,krf
        IF(irf(j) == 0)THEN
           k = k+1
           IF(k .LE. nsys)THEN 
              irf(j) = 3 
              irfmodel(j) = 'ech'
              rfmode(j)   = 'ech'
              irfcur_12(j) = irfcur(k) 
           ENDIF
        ENDIF
     ENDDO

  END SUBROUTINE set_irf


      SUBROUTINE tableintrp_SP (xx, n, x, jlo)
! ----------------------------------------------------------------------
! --- correlated table search routine. USE jlo from previous CALL to
! --- get jlo for current value of x. IF jlo from previous CALL is
! --- no good(jlo=0 or jlo=n) THEN DO a binary search.
! --- this routine returns jlo so that
!   xx(jlo) .LE. x .LE. xx(jlo+1) IF jlo = 1
!   xx(jlo) .LT. x .LE. xx(jlo+1) IF jlo =2,3...n-1
!   it is assumed that xx(j),j = 1,2..n is monotonically
!   increasing; this is NOT checked for.
!   this is a modified version of the Numerical Recipes SUBROUTINE HUNT
! ------------------------------------------------------------------ HSJ
!
!
      USE nrtype,               ONLY : DP,SP,I4B
      IMPLICIT NONE
!
      INTEGER(I4B)  n, jlo, jhi, jmid, inc
      REAL(DP)   x, xx(n)
!
      IF      (x .LT. xx(1)) THEN
        jlo = 0                     ! indicates out of range below xx(1)
      ELSE IF (x .LE. xx(2)) THEN
        jlo = 1                     ! xx(1) .le. x .le. xx(2)
      ELSE IF (x .LE. xx(n)) THEN   ! xx(2) .lt. x .le. xx(n)
!
!       check IF jlo from previous CALL is usable:
!
        IF (jlo .LE. 0 .OR. jlo .GE. n) THEN
          jlo = 2
          jhi = n
          go to 15           ! no correlation go directly to bisection
        ELSE                 ! 1 .le. jlo .lt. n
!
!         bracket x, THEN USE bisection:
!         start WITH  jlo  from previous CALL:
!
          inc = 1  
          IF (x .GT. xx(jlo)) THEN    ! search up
    4       jhi = jlo + inc
            IF (jhi .GE. n) THEN
              jhi = n
            ELSE IF (x .GT. xx(jhi)) THEN
              inc = inc + inc
              jlo = jhi
              go to 4
            END IF
          ELSE
    5       jhi = jlo
            jlo = jlo-inc
            IF (jlo .LE. 1) THEN
              jlo = 1
            ELSE IF (x .LE. xx(jlo)) THEN
              inc = inc + inc
              go to 5
            END IF
          END IF
        END IF
!
!         bisection:
!
   10     IF (jhi-jlo .EQ. 1)  RETURN
   15     jmid = (jhi + jlo) / 2
          IF (x .GT. xx(jmid)) THEN
            jlo = jmid
          ELSE           ! x .le. xx(jmid)
            jhi = jmid
          END IF
          go to 10
!
      ELSE               ! x .gt. xx(n)
        jlo = n
      END IF
!
      RETURN
!
      END SUBROUTINE tableintrp_SP



SUBROUTINE ech_err(status,flag)

  USE nrtype,                            ONLY : I4B,DP

 
  USE error_handler,                     ONLY : lerrno,terminate

#ifndef ONETWO
      USE io_gcnmp,                      ONLY : ncrt,nlog
#else
      USE io,                            ONLY : ncrt,nlog => ncrt
#endif

  IMPLICIT NONE
  INTEGER, INTENT ( in) :: status
  INTEGER(I4B)  action,vartyp,nvdims,nvatts,rcode
  INTEGER,PARAMETER                       :: MAXNCNAM = 128 !these are supposed to be
  INTEGER,PARAMETER                       :: MAXVDIMS = 32  
  INTEGER(I4B) ,DIMENSION(MAXVDIMS)       :: vdims
  INTEGER(I4B), OPTIONAL                  :: flag
  CHARACTER(LEN=MAXNCNAM)                 :: varname
  action = 10000
  IF(PRESENT(flag))THEN
        action = flag
  ENDIF
  IF(status /= nf90_noerr .AND. action .NE.  0 ) THEN
     WRITE(ncrt,FMT='("  netdcf erro = ",i5)')status
     WRITE(nlog,FMT='("  netdcf erro = ",i5)')status
     WRITE(ncrt,1)TRIM(nf90_strerror(status))
     WRITE(ncrt,FMT='(a)')trim(nf90_strerror(status))
     WRITE(nlog,1)TRIM(nf90_strerror(status)) 
     WRITE(nlog,FMT='(a)')trim(nf90_strerror(status))
     IF(action > 0)THEN
        CALL ncvinq(ncid,action,varname,vartyp,nvdims,vdims,nvatts,rcode)
        WRITE(ncrt,2)action,varname(1:LEN_TRIM(varname))
        WRITE(nlog,2)action,varname(1:LEN_TRIM(varname))
     ELSE IF(action < 0)THEN
        action = -action
        CALL ncvinq(ncid,action,varname,vartyp,nvdims,vdims,nvatts,rcode)
        WRITE(ncrt,3)action,varname(1:LEN_TRIM(varname))
        WRITE(ncrt,FMT='(a)')trim(nf90_strerror(status))
        WRITE(nlog,3)action,varname(1:LEN_TRIM(varname))
        WRITE(nlog,FMT='(a)')trim(nf90_strerror(status))
     ENDIF
!     lerrno = 32
!     CALL  terminate(lerrno,nlog)
  ELSE IF(status /= nf90_noerr .AND. action ==  0 ) THEN
      flag =  -1
  END IF

  RETURN

3    FORMAT(2x,'READ  error caused by netcdf variable no :',i5,' name: ',a)
2    FORMAT(2x,'WRITE error caused by netcdf variable no :',i5,' name: ',a)
1    FORMAT(2x,a)

END SUBROUTINE ech_err

 
  END MODULE echdat_module
