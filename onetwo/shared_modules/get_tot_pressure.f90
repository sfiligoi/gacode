   SUBROUTINE get_tot_pressure(fi_inc,p_type)
! ------------------------------------------------------------------------------------------
! -- get presure in joules/m**3
! -- fi_inc:  include or exclude fast ions(fusion alphas, beam - poor approximation for stored energy density)
! -- p_type determines which set of vriables are used to do the calcualtions.
! -- eq tot_press,ion_press or profile%(tot_press,ion_press,elect_press)
! -- these may not be exchangeable because they exist on different sized grids
! ----------------------------------------------------------------------------HSJ-----------
      USE nrtype,                  ONLY : DP,I4B

      USE dep_var,                 ONLY : te,ti,en,ene,tot_press,ion_press
      
      USE grid_class,              ONLY : nj,hcap,r,dr

      USE ions_gcnmp,              ONLY : nion

      USE fast_ion_data_gcnmp,     ONLY : w_alpha

      USE neutral_beams,           ONLY : nbion,wbeam

      USE common_constants,        ONLY : joupkev,zeroc,PISQ

      USE Plasma_properties,       ONLY : profile,mhd_dat,dischg

      USE error_handler,           ONLY : iomaxerr, lerrno,terminate

      USE io_gcnmp,                ONLY : nlog

      USE vector_class,            ONLY : crossprod_Vector,zero_Vector,length_Vector

      IMPLICIT NONE

      INTEGER(I4B) j,k,fi_inc,p_type,p_nj


      IF(p_type ==0)THEN
        DO j =1,nj
           tot_press(j) = ene(j)*te(j)                                          ! Kev/m**3
           DO k=1,nion
              tot_press(j) = tot_press(j)+en(j,k)*ti(j)
           ENDDO
           IF(fi_inc == 1)                                        &
           tot_press(j) = tot_press(j)+0.67*(wbeam(j) + w_alpha(j))

           tot_press(j) = tot_press(j)*joupkev                                  ! joules/m**3
           ion_press(j) = tot_press(j) -ene(j)*te(j)*joupkev
        ENDDO

     ELSE

        p_nj = length_vector(profile%te)
        IF(p_nj .NE. nj)THEN
           PRINT *,'p_nj .NE. nj in get_tot_pressure'
           lerrno = iomaxerr + 154
           CALL terminate (lerrno,nlog)
        ENDIF

        profile%elect_press  = crossprod_Vector(profile%ene,profile%te,joupkev)  ! joules/m**3
        profile%tot_press  = zero_vector(nj)
        profile%ion_press  = zero_vector(nj)


        DO j =1,nj
           DO k=1,nion
              profile%ion_press%data(j) = profile%ion_press%data(j)+  &
                                          profile%en(k)%data(j)*profile%ti%data(j)*joupkev 
           ENDDO
           IF(fi_inc == 1)                                            &
           profile%ion_press%data(j)    = profile%ion_press%data(j)   +  &
                                          0.67*(profile%wbeam%data(j) +  &  ! note crude approx
                                             profile%walp%data(j))*joupkev  ! joules/m**3

           profile%tot_press%data(j)    = profile%tot_press%data(j)         +      &
                                                profile%elect_press%data(j) +      &
                                                profile%ion_press%data(j)
        ENDDO


    ! integrate over volume

        profile%p_intg = zeroc
        Do j=1,nj-1

           profile%p_intg = profile%p_intg + (profile%tot_press%data(j)*hcap(j)*r(j) +  &
                                  profile%tot_press%data(j+1)*hcap(j+1)*r(j+1)) *dr(j)

        ENDDO

        profile%p_intg = profile%p_intg*2.*PISQ*dischg%rmajor  ! note 2pi not 4pi because of sum above
     ENDIF


       RETURN
    END SUBROUTINE get_tot_pressure
