!----- this is for Onetwo only-------------------------------------------
! --   (shared modules also has such a module for other codes)
!-------------------------------------------------------------------------
  MODULE fast_ion_diffusion
    IMPLICIT NONE
    TYPE,PUBLIC :: fidif
       REAL*8 adiff_0,adiff_a,adiff_xpin,adiff_xpout
       INTEGER fidif_on,    &                !=1 if model is used
               nkdifb,      &                !=3,2,1 determines what fast ion
                                             ! species are diffused
               ndifbe                        ! no of bins in fdifbe,edfibe
                                             ! 20 is assumed and hard wired 
                                             ! in may places
       REAL *8,DIMENSION(:),POINTER :: fdifbe,  & ! multiplier for xchi fi  in
                                       edifbe    !each energy bin edifbe IN EV
    END TYPE fidif
 
!  thesed values are used in NUBEAM to generate a profile:
!  difb =  adiff_a +(adiff_0-adiff_a)*(1.-rho^adiff_xpin)^adiff_xpout
!  IF adiff_xpin le 0.0 then difb = adiff _a
!  IF adiff_xpout .le. 0 then difb = adiff_0

   CONTAINS 
     
   SUBROUTINE fi_init(fidifin,size)
     TYPE(fidif),INTENT(out) :: fidifin
     INTEGER,INTENT(IN) :: size
     fidifin%ndifbe = size
     IF(fidifin%ndifbe .GT.0)THEN
          ALLOCATE(fidifin%fdifbe(fidifin%ndifbe))
          ALLOCATE(fidifin%edifbe(fidifin%ndifbe))
          fidifin%fdifbe(:) = 0.0
          fidifin%edifbe(:) = 0.0
     ENDIF
     if(fidifin%ndifbe .gt. 0)fidifin%fidif_on = 1

     RETURN
   END SUBROUTINE fi_init



  END MODULE fast_ion_diffusion
