  MODULE fast_ion_data_gcnmp
    USE nrtype ,   ONLY : DP, I4b
    IMPLICIT NONE
    REAL(DP),SAVE,ALLOCATABLE,DIMENSION(:)   :: enalp,walp,w_alpha
!   walp is initial condition for w_alpha, see fusion routines
  END MODULE fast_ion_data_gcnmp
