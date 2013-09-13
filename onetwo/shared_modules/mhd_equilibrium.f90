MODULE mhd_equilibrium

  USE nrtype,                          ONLY : DP,I4B

  USE Plasma_properties,               ONLY : mhd_dat
  
  USE common_constants,                ONLY : izero,zeroc



  CONTAINS

    SUBROUTINE equilibrium_driver
        
        IMPLICIT NONE
        REAL(DP) delta_P
        IF(ABS(mhd_dat%run_mhd) .EQ. izero )THEN ! equilibirum calculations 
                                                 ! are not active
           RETURN

        ELSE                ! equilibirum calculations are active.
                            ! check if new equilibrium should be calculated

            mhd_dat%equil_avail(1) = 'fixedbd_var_grid'
            mhd_dat%equil_avail(2) = 'fixedbd_unf_grid'


            If(mhd_dat%run_mhd .GT. izero)THEN
               ! mhd solver will be called if more than mhd_dat%run_mhd 
               ! time steps have been taken
               IF(mhd_dat%run_mhd .LE. mhd_dat%mhd_ctr)THEN 
                   mhd_dat%mhd_ctr =0               ! reset ctr for next pass
                   CALL gs_solver
               ENDIF

            ELSEIF(mhd_dat%run_mhd .LT. izero)THEN
               ! mhd solver will be called if max change in pressure
               ! since last equilibrium is exceeded
                 CALL pressure_change(delta_P)

                 IF(FLOAT(ABS(mhd_dat%run_mhd)) .GE. delta_P)  &
                      CALL gs_solver

            ENDIF

        ENDIF

    END  SUBROUTINE equilibrium_driver


    SUBROUTINE  gs_solver

        USE error_handler,              ONLY : mhd_offset, lerrno,terminate

        USE io_gcnmp,                   ONLY : ncrt,nlog

        USE MPI_data,                   ONLY : myid,master

        IMPLICIT NONE

         IF(mhd_dat%equil_type == mhd_dat%equil_avail(1))THEN
            ! 'fixedbd_var_grid'
            CALL gs_solver_var_grid
         ELSEIF(mhd_dat%equil_type == mhd_dat%equil_avail(2))THEN
            ! 'fixedbd_unf_grid'
            CALL gs_solver_unf_grid
         ELSE
            IF(myid == master)THEN
               WRITE(ncrt,FMT='("Error in input value of mhd_dat%equil_type",/,&
                    " value not recognized : ",a)')mhd_dat%equil_type
            ENDIF

            lerrno = mhd_offset
            CALL terminate(lerrno,nlog)

         ENDIF

       RETURN

    END SUBROUTINE  gs_solver


    SUBROUTINE  gs_solver_var_grid

       IMPLICIT NONE

       RETURN

    END SUBROUTINE  gs_solver_var_grid



    SUBROUTINE  gs_solver_unf_grid

       IMPLICIT NONE

       RETURN

    END SUBROUTINE  gs_solver_unf_grid



END MODULE mhd_equilibrium
