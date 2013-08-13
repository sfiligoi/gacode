    SUBROUTINE pressure_change(delta_P)
!----------------------------------------------------------------------------
! -- return the largest relative change in pressure encountered over the grid
!----------------------------------------------------------------------------
      USE nrtype,                                  ONLY : DP,I4B
             
      USE dep_var,                                 ONLY : tot_press


      IMPLICIT NONE
      REAL(DP) delta_P
      REAL(DP),ALLOCATABLE,DIMENSION(:) :: tot_press_local
      INTEGER(I4B) fi_inc,p_type

      ALLOCATE(tot_press_local(SIZE(tot_press)))


      tot_press_local(:) = tot_press(:)               ! save previous pressure
 
      
      CALL get_tot_pressure(fi_inc,p_type)            ! get current pressure

      
      
      DEALLOCATE(tot_press_local)

      RETURN

    END SUBROUTINE pressure_change
