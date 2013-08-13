!
  MODULE fast_ion_diffusion
    USE nrtype,                      ONLY : DP,I4B,SP
    USE common_constants,            ONLY : zeroc,izero
    USE nf_param,                    ONLY : ke
    USE dnubeam_mod,                 ONLY : ntimemax
    IMPLICIT NONE
    TYPE,PUBLIC :: fidif
       REAL(DP),DIMENSION(ntimemax,ke) :: adiff_0,adiff_a,adiff_xpin,adiff_xpout
       REAL(DP),DIMENSION(ntimemax)    :: adiff_time
       !  these values are used in NUBEAM to generate a profile:
       !  difb =  adiff_a +(adiff_0-adiff_a)*(1.-rho^adiff_xpin)^adiff_xpout
       !  IF adiff_xpin le 0.0 then difb = adiff _a
       !  IF adiff_xpout .le. 0 then difb = adiff_0
       !  IF addif_a .le. 0 and adiff_0 .le. 0 then dont use fast ion diffusion

       INTEGER fidif_on,    &                !=1 if model is used
               nkdifb,      &                !=3,2,1 determines what fast ion
                                             ! species are diffused
               ndifbe,      &                ! no of bins in fdifbe,edfibe
                                             ! 20 is assumed and hard wired 
                                             ! in may places
               adiff_ntime                   ! number of time values for adiff_a,etc
       REAL (DP),DIMENSION(:),POINTER :: fdifbe,  & ! multiplier for xchi fi  in
                                         edifbe     ! each energy bin edifbe IN EV
       REAL (DP), DIMENSION(:,:,:),POINTER :: fi_d,fi_den ! space,energy,beam no
       REAL (DP) fi_na
       REAL(DP) adiff_al,adiff_0l,adiff_xpinl,adiff_xpoutl,anal_time
    END TYPE fidif


    !tridiagonal system storage:
    REAL(DP), ALLOCATABLE,PRIVATE,DIMENSION(:) :: ld,        & ! lower diagona
                                                   d,        & ! diagonal
                                                  ud,        & ! upper diagonal
                                               alpha,        & ! grid factors
                                               beta
    LOGICAL, PUBLIC ::  fi_refactor_c



   CONTAINS 
     
  

   SUBROUTINE check_d_fast_ion_inpt(fidat,fdifbe,edifbe,ndifbep)

     USE nf_param,                             ONLY : ke,kb

     USE grid_class,                           ONLY : nj

     USE common_constants,                     ONLY : izero,ZEROC
          
     USE error_handler,                        ONLY : lerrno,terminate,iomaxerr
        
     USE io_gcnmp,                             ONLY : ncrt,nlog



     INTEGER(I4B) jt,ie,ndifbep
     REAL(DP) fdifbe(ndifbep),edifbe(ndifbep)
     TYPE(fidif),INTENT(inout) :: fidat


     fidat%fidif_on   = izero
     IF(fidat%ndifbe .GT. 0)THEN
         IF(ndifbep .LT. fidat%ndifbe)THEN
            PRINT *,'ERROR, ndifbep = ',ndifbep
            PRINT *,'But requested array is size',fidat%ndifbe
            lerrno = 257 + iomaxerr
            CALL terminate(lerrno,nlog)
         ENDIF
         IF(ASSOCIATED(fidat%fdifbe))DEALLOCATE(fidat%fdifbe)
         IF(ASSOCIATED(fidat%edifbe))DEALLOCATE(fidat%edifbe)
         ALLOCATE(fidat%edifbe(fidat%ndifbe),fidat%fdifbe(fidat%ndifbe))
         fidat%fdifbe(1:fidat%ndifbe) = fdifbe(1:fidat%ndifbe)
         fidat%edifbe(1:fidat%ndifbe) = edifbe(1:fidat%ndifbe)
     ENDIF
!     if(fidifin%ndifbe .gt. 0)fidat%fidif_on = 1 
!     P_Nfreya doesn't currently use fdifbe,edifbe



!    P_Nfreya may use diffusion coefficient:
! define tese arrays only when neede on output hsj 9/22/11
!     IF(ASSOCIATED(fidat%fi_d))DEALLOCATE(fidat%fi_d)
     DO jt = 1,ntimemax
        DO ie =1,ke
           IF(fidat%adiff_0(jt,ie) .GT. zeroc .OR. fidat%adiff_a(jt,ie) .GT. zeroc &
                .AND. .NOT. ASSOCIATED(fidat%fi_d)) THEN
              fidat%fidif_on = 1
              !ALLOCATE(fidat%fi_d(nj,ke,kb))
              !ALLOCATE(fidat%fi_den(nj,ke,kb))
           ENDIF
        ENDDO
     ENDDO
!     IF(ASSOCIATED(fidat%fi_d)) fidat%fidif_on = 1

     RETURN

  END SUBROUTINE check_d_fast_ion_inpt





   SUBROUTINE fi_steady_state(fi_s,fi_d,fi_na,fi_n)
!-----------------------------------------------------------------------------------------
! -- subroutine solves the steady state fast ion diffusion equation in space (not energy)
! -- using the source density fi_s   (#/(M^3 sec)) and
! -- a diffusion coefficient fi_d
! -- the equation is 
! -- (1/hr)d/dr(-hr fi_d  dn/dr) = fi_source
! --
! INPUT 
!    fidif_on         =1 if diffusion is to be applied , = 0 otherwise
!    nkdifb           = 1 is assume if fidif+on = 1
!    fi_d(:)          M^2/sec diffusion coefficient
!    fi_s(:)             source of fast ions  
!    fi_na            fast ion density at rho =a ,boundary condition
!    fi_refactor_c      true or false. if true setup matrix for finite diff diffusion equation
! from grid_class:
!    nj                   size of radial grid
!    hcap(1..nj)          geometric factor (dv = 4pisq ro hcap rho drho)
!    r(1..nj)             radial grid in meters
!    dr(1..nj-1)          r(i+1)-r(i)
!
! OUTPUT
!    fi_n(1..nj)          steady state fast ion density 
!--------------------------------------------------------------------HSJ-------------------

   USE grid_class,                        ONLY : r,dr,hcap,nj  !  note must be consistent with grid

   USE tension_spline,                    ONLY : progon3       !  tridiagonal solver

   INTEGER(I4B) i

   REAL(DP) rhs(nj)                                            ! temp arrays

   REAL(DP) fi_s(:),fi_d(:),fi_na,fi_n(:)

   IF(fi_refactor_c)THEN                 ! grid dependent quantites, set only when grid changes

      IF( ALLOCATED(d))THEN  ! assume grid change occured or first entry
           DEALLOCATE(alpha,beta,ld,d,ud)
      ENDIF
      ALLOCATE(alpha(nj),beta(nj),ld(nj),d(nj),ud(nj))



      fi_refactor_c  = .FALSE.        ! turn off for subsequent calls until turned on again externally
      alpha(:) = hcap(:)*r(:)
      beta(1) = (hcap(2)*r(2)*fi_d(2)+hcap(1)*r(1)*fi_d(1))/(2._DP*dr(1)) !( r(1) ==0)
      ld(1)    = zeroc ! lower diagonal not used
      d(1)     = 1._DP ; ud(1)    = -d(1) ! forces dn/drho =0 at rho =0
      DO i=2,nj
         ld(i) = beta(i-1)
         IF( i .lt. nj) beta(i) =0.5_DP*(hcap(i+1)*r(i+1)*fi_d(i+1)+hcap(i)*r(i)*fi_d(i))/dr(i)
         d(i)  = -(beta(i)+beta(i-1))
         ud(i) = beta(i)
      ENDDO
      ud(nj) = zeroc ! upper diagonal not used
      ld(nj) = zeroc 
      d(nj)  = 1._DP ! forces fi_n(at rho = a) to value fi_na
   ENDIF

   rhs(1) = zeroc
   rhs(nj) = fi_na
   DO i=2,nj-1
      rhs(i) = -fi_s(i)*alpha(i)*(r(i+1)-r(i-1))/2._DP
   ENDDO

   CALL Progon3 (d,ud,ld,rhs,fi_n,nj)

   RETURN 

   END SUBROUTINE  fi_steady_state


   SUBROUTINE get_d_fast_ion(fidifin,ie,jb)
!-----------------------------------------------------------------------
! -- set up diffusion coefficeint array fi_d
! -- using nubeam parameterization
! -- INPUT:
!    arg list:
!         ie                               energy 1,2,3 (full half third)
!    from grid_class:
!         nj                               size of radial grid
!         r(1..nj)                         radial grid in meters
!    from fast_ion_diffusion               (local module)
!         adiff_0(1:ke)                     [m^2/sec] on axis d
!         adiff_a(1:ke)                     [m^2/sec] edge d 
!         adiff_xpin(1:ke)                 (pure number)
!         adiff_xpout(1:ke)                (pure number)
!
! -- OUTPUT
!     fidifin%fi_d(1:nj,ie,jb)  =  adiff_a +(adiff_0-adiff_a)*(1.-rho^adiff_xpin)^adiff_xpout
!                                  it is assumed taht adiff_a and adiff_0 are in m^2/sec 
!------------------------------------------------------------------------------HSJ--------
     USE grid_class,                              ONLY : roa,nj 
     USE solcon_gcnmp,                            ONLY : time        ! this time values from statefile
     USE neutral_beams,                           ONLY : beam_sim_time_start,beam_sim_time_end


     IMPLICIT NONE

 
     REAL(DP) dt,dti,df

     TYPE(fidif),INTENT(inout) :: fidifin


     INTEGER(I4B)  j,ie,jb,jlo

  ! NOTE currently no dependance on beamlet nmbr jb:
  ! But introduced time dependence 2/15/2012  HSJ

     fidifin%adiff_al       = fidifin%adiff_a(1,ie)
     fidifin%adiff_0l       = fidifin%adiff_0(1,ie)
     fidifin%adiff_xpinl    = fidifin%adiff_xpin(1,ie)
     fidifin%adiff_xpoutl   = fidifin%adiff_xpout(1,ie)

     fidifin%anal_time = (beam_sim_time_start + beam_sim_time_end)*0.5_DP
     IF(fidifin%adiff_ntime .GT. 1)THEN
        IF(fidifin%adiff_time(1) .LE. fidifin%anal_time   .AND.    &
           fidifin%anal_time .LE.  fidifin%adiff_time(fidifin%adiff_ntime))THEN
           jlo = izero
           CALL tableintrp (fidifin%adiff_time,fidifin%adiff_ntime,fidifin%anal_time, jlo)
           dt     = fidifin%anal_time - fidifin%adiff_time(jlo)
           dti    = fidifin%adiff_time(jlo+1)-fidifin%adiff_time(jlo)
           df     = (fidifin%adiff_a(jlo+1,ie)-fidifin%adiff_a(jlo,ie))*dt/dti
           fidifin%adiff_al = fidifin%adiff_a(jlo,ie) + df
           df     = (fidifin%adiff_0(jlo+1,ie)-fidifin%adiff_0(jlo,ie))*dt/dti
           fidifin%adiff_0l = fidifin%adiff_0(jlo,ie) + df
           df     = (fidifin%adiff_xpin(jlo+1,ie)-fidifin%adiff_xpin(jlo,ie))*dt/dti
           fidifin%adiff_xpinl = fidifin%adiff_xpin(jlo,ie) + df
           df     = (fidifin%adiff_xpout(jlo+1,ie)-fidifin%adiff_xpout(jlo,ie))*dt/dti
           fidifin%adiff_xpoutl = fidifin%adiff_xpout(jlo,ie) + df
        ENDIF
     ENDIF

     IF(ASSOCIATED(fidifin%fi_d))THEN
        ! fidifin%fi_d will not be associated yet  when this routine is called from
        ! Check_P_Nfreya_beamlets (nfreya_load.f90) (data isnt needed until 
        ! beam deposition is done so this is OK)
        DO j=1,nj
           fidifin%fi_d(j,ie,jb) = fidifin%adiff_al+(fidifin%adiff_0l-fidifin%adiff_al)* &
                          (1._DP-roa(j)**fidifin%adiff_xpinl)**fidifin%adiff_xpoutl
        ENDDO
     ENDIF

     RETURN
   END SUBROUTINE get_d_fast_ion


  END MODULE fast_ion_diffusion
