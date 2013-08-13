   
    SUBROUTINE check_en_bc
! -------------------------------------------------------------------
! --- check for setting of en_bc,ren_bc
! --- reallocate these arrays as appropriate
! -------------------------------------------------------------------
       USE nrtype,                      ONLY : I4B,DP
       USE param,                       ONLY : ksplin,kion,kbctim,kj
       USE numbrs,                      ONLY : nj,nprim,nion,nbctim
       USE bd_condtn,                   ONLY : en_bc_inpt,ren_bc_inpt,  &
                                               ken_bc_inpt,renpin,      &
                                               reniin,enp,eni,          &
                                               fix_edge_ni_bc_inpt
       USE mesh,                        ONLY : fix_edge_ni
       USE io,                          ONLY : ncrt,nout
       USE common_constants,            ONLY : zeroc

       IMPLICIT NONE
       INTEGER(I4B) j,i,k,ierr,knots,itry,l
       REAL(DP) pos_def,one
       REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: ent
       INTEGER(I4B),DIMENSION(:,:),ALLOCATABLE :: knt

       pos_def = 5._DP*EPSILON(0.0_DP)
       one = 1.0_DP
       !en_bc,ren_bc were loaded with -1. on startup.
       !check if user has changed values:

       IERR = 0 
       !first determine the number of knots,ken_bc:

       DO i = 1,nion
            DO j=nbctim,1,-1 
               itry =0
10                CAll get_knot_num(ren_bc_inpt(1,i,j),ken_bc_inpt(i,j),ksplin)
                IF(ken_bc_inpt(i,j) .LT. 3)THEN  ! assume profiles are not set
                   IF(i .LE. nprim)THEN
                      ren_bc_inpt(1:kj,i,j) = renpin(1:kj,i)
                      en_bc_inpt(1:kj,i,j)  = enp(1:kj,i)
                   ELSEIF(i .gt. nprim .and. i .LE. nion)THEN
                      l=i-nprim
                      ren_bc_inpt(1:kj,i,j) = reniin(1:kj,l)
                      en_bc_inpt(1:kj,i,j)  = eni(1:kj,l)
                   ENDIF
                   itry =itry +1
                   IF(itry .lt. 2)THEN 
                      go to 10
                   ELSE
                      IERR=1
                      EXIT
                   ENDIF
                 ENDIF
              !check for edge boundary value settings. Assume these values are not
              !set if knots were not set (itry >0)
              !In this case the values may be time dependent but are the
              !same for all ions at each time:
              IF(itry .GT. 0) fix_edge_ni_bc_inpt(i,j) = fix_edge_ni(j)      
            ENDDO
       ENDDO
        
        !now check for correct values:
     IF(IERR == 0)THEN
        DO i=1,nion
           DO j=1,nbctim
              IF( ken_bc_inpt(i,j) .GT. 3)THEN
                 IF(ABS(ren_bc_inpt(1,i,j)) .GT. pos_def)THEN
                    WRITE(ncrt,2)
                    WRITE(nout,2)
                    IERR =2
                 ELSE
                    ren_bc_inpt(1,i,j) = zeroc
                 ENDIF
                 IF(ABS(ren_bc_inpt(ken_bc_inpt(i,j),i,j)- one) .GT. pos_def )THEN
                    IERR =3
                    WRITE(ncrt,3)
                    WRITE(nout,3)
                 ELSE
                    ren_bc_inpt(ken_bc_inpt(i,j),i,j)  = 1.0_DP
                 ENDIF
              ELSE
                 IERR =1
              ENDIF
              DO k=1,ken_bc_inpt(i,j)
                 IF(en_bc_inpt(k,i,j) .LT. pos_def)THEN 
                    print *,'knots ,ion,nbctim =',k,i,j
                    print *,'en_bc_inpt(k,i,j) =',en_bc_inpt(k,i,j)
                    WRITE(ncrt,4)
                    WRITE(nout,4)
                    IERR = 4
                 ENDIF
              ENDDO

           ENDDO
        ENDDO
     ENDIF







        !now redimension en_bc_inpt,ren_bc_inpt,ken_bc_inpt to the sizes actually used:
        !We dont go to disk because the arrays arent that large
        
!!$        IF(IERR  == 0)THEN 
!!$           
!!$           ALLOCATE(ent(ksplin,nion,nbctim))
!!$
!!$           ent(1:nj,1:nion,1:nbctim) =en_bc_inpt((1:nj,1:nion,1:nbctim)
!!$           DEALLOCATE(en_bc_inpt)
!!$           ALLOCATE(en_bc_inpt(nj,nion,nbctim))
!!$           en_bc_inpt(:,:,:) = ent(:,:,:)
!!$
!!$           ent(1:nj,1:nion,1:nbctim) =ren_bc_inpt((1:nj,1:nion,1:nbctim)
!!$           DEALLOCATE(ren_bc_inpt)
!!$           ALLOCATE(ren_bc_inpt(nj,nion,nbctim))
!!$           ren_bc_inpt(:,:,:) = ent(:,:,:)
!!$
!!$           DEALLOCATE(ent)
!!$           ALLOCATE(knt(nion,nbctim))
!!$
!!$           knt(1:nion,1:nbctim) =ken_bc_inpt(1:nion,1:nbctim)
!!$           DEALLOCATE(ken_bc_inpt)
!!$           ALLOCATE(ken_bc_inpt(nion,nbctim))
!!$           ken_bc_inpt(:,:) = knt(:,:)
!!$
!!$           DEALLOCATE(knt)
!!$
!!$        ENDIF
        

       IF(IERR == 1)CALL STOP('error in check_en_bc,ren_bc not set',1)
       IF(IERR == 2)CALL STOP('error in check_en_bc,ren_bc 0.0  not set',1)
       IF(IERR == 3)CALL STOP('error in check_en_bc,ren_bc 1.0  not set',1)
       IF(IERR == 4)CALL STOP('error in check_en_bc,en_bc_inpt not positive',1)
1      FORMAT(2x,'Error in spline input of ren_bc',/, &
              2x,'Must have at least 3 knots in spline')
2      FORMAT(2x,'Error in spline input of ren_bc',/, &
              2x,'First knot must be at = 0.0')
3      FORMAT(2x,'Error in spline input of ren_bc',/, &
              2x,'First last knot  must be at = 1.0')
4      FORMAT(2x,'Error in spline input of en_bc',/, &
              2x,'All values must be >= 0.0')

       RETURN
    END SUBROUTINE check_en_bc


    SUBROUTINE get_knot_num(r,knots,ksplin)
      !   check array r for # knots
      USE nrtype,                     ONLY : I4B,DP
      IMPLICIT NONE
      REAL(DP) r(ksplin),pos_def
      INTEGER(I4B) knots,j,ksplin
      pos_def = 5._DP*epsilon(0.0_DP)
      knots = -1
      DO j=2,ksplin
        IF(ABS(r(j)-1.0_DP)  .LT. pos_def)THEN
           knots = j
           EXIT
        ENDIF
      ENDDO
      IF(ABS(r(1)) .GT. pos_def)knots = -2

      RETURN

    END SUBROUTINE get_knot_num


 
