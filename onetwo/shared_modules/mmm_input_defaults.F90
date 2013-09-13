
    SUBROUTINE mmm_input_defaults
    !---------------------------------------------------------------------------
    ! -- set input defaults for mmm model
    !----------------------------------------------------------HSJ-2//29//2012--
        USE nrtype,                           ONLY : DP,I4B

        USE mmm_in,                           ONLY : mmm_d_mult,mmm_v_mult,  &
                                                     mmm_cmodel,mmm_cw20,    &
                                                     mmm_cDBM,mmm_cW20,      &
                                                     mmm_cETG,mmm_lw20,      &
                                                     mmm_lDBM,mmm_lETG,      &
                                                     mmm_lswitch,mmm_cswitch,&
                                                     BADREAL,BADINT,use_mmm_zct

        USE modmmm7_1,                        ONLY : set_mmm7_1_switches

        USE common_constants,                 ONLY : izero,zeroc             

        IMPLICIT NONE

            mmm_d_mult  =  1.0_DP   ;   mmm_v_mult = 1.0_DP

            mmm_cmodel(:)  = (/1._DP,1._DP,1._DP /)   
            ! mmm_cmodel weights:      
            !        mmm_cmodel(1)  = 1.D0  ==>    Weiland
            !        mmm_cmodel(2)  = 1.D0  ==>    DRIBM
            !        mmm_cmodel(3)  = 1.D0  ==>     ETG

            ! cW20, cDBM, cETG ! To be passed to mmm_cswitch
            mmm_cW20(:)  = BADREAL
            ! 1D0 ExB shear coefficient
            ! 1D0 Momentum pinch scaling factor
            ! 1D-4 Lower bound of electron thermal diffusivity
            ! 1D2 Upper bound of electron thermal diffusivity
            ! 1D-4 Lower bound of ion thermal diffusivity
            ! 1D2 Upper bound of ion thermal diffusivity

            mmm_cDBM(:) = BADREAL
            mmm_cETG(:) = BADREAL

            ! lW20, lDBM, lETG ! To be passed to mmm_lswitch
            mmm_lW20(:)     = BADINT
            mmm_lDBM(:)     = BADINT
            mmm_lETG(:)     = BADINT

!write(940,FMT='("mmm_input_defaults strt")') ! 8888889999
!write(940,FMT='((5(1pe12.2,x)))')mmm_cswitch
!write(940,FMT='((5(i3,x)))')mmm_lswitch
            !Fill mmm parameter arrays with default  values
            CALL set_mmm7_1_switches( cmmm = mmm_cswitch, lmmm = mmm_lswitch )
            use_mmm_zct = .TRUE.
!write(940 ,FMT='("mmm_input_defaults end")') ! 8888889999
!write(940,FMT='((5(1pe12.2,x)))')mmm_cswitch
!write(940,FMT='((5(i3,x)))')mmm_lswitch
        RETURN

    END SUBROUTINE mmm_input_defaults

    SUBROUTINE mmm_check_input( vval, vname )
!------------------------------------------------------------------------
! This internal subroutine checks for variable with BADREAL
! If so, print an error message and increase the error count
!-------------------------------------------------------------------------
    USE nrtype,                           ONLY : DP,I4B

    USE mmm_in,                           ONLY :  BADREAL
    
    USE io_gcnmp,                         ONLY :  nlog,ncrt

    USE error_handler,                    ONLY : lerrno,iomaxerr,terminate

      IMPLICIT NONE

      REAL(DP), INTENT(In) :: vval          ! the value to be checked
      CHARACTER(len=*), INTENT(In) :: vname ! the corresponding variable name
      CHARACTER(2) :: vname_len             ! length of the variable name

      IF ( ABS( vval - DBLE(BADREAL) ) < 1.e-6 ) THEN
         WRITE(vname_len,'(I2)') LEN_TRIM(vname)
         WRITE(ncrt,'("** ERROR ** """A'//TRIM(ADJUSTL(vname_len))//&
              '""" is not set correctly ")') vname
         lerrno = iomaxerr + 185_I4B
         CALL terminate(lerrno,nlog)
      END IF

    END SUBROUTINE mmm_check_input 

