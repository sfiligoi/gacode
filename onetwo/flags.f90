   MODULE flags
! --- INCLUDE file flags.i
      USE param,only : kk
!
      INTEGER                                                          &
                    testing_NTCC, steps_per_plot,do_eqplot,            &
                    continuation_method,analysis_check,                &
                    first_step,no_w_convection,                        &
                    itran(kk), inenez, inrad, iten,                    &
                    no_te_convection, no_ti_convection,                &
                    set_te_to_ti,set_ti_to_te
!
      CHARACTER (len=5)  itranflag(kk)
      CHARACTER (len=16)  interp_prof
!
!
      REAL*8         zfrac
!
   END MODULE flags
