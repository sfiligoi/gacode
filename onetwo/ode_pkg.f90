    MODULE ode_pkg

     IMPLICIT NONE

!     module  for ODE solver package
!     flux => g1e, g1i, energy density => g2
!
      REAL*8  g1e_ode, g1i_ode, g2c_ode, ene_scale, eni_scale,        &
              dummy_ode(48), sdummy_ode(4)
      INTEGER  nlc_ode, nlu_ode, idummy_ode(38), iodeflg1
!
      DATA g1e_ode, g1i_ode, g2c_ode ,ene_scale, eni_scale   &   ! flux => 5/2, energy density => 3/2
        /    2.5,     2.5,     1.5,   1.0,   1.0/
    END MODULE ode_pkg
