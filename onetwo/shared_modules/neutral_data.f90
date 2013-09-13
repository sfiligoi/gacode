  MODULE neutral_data
    USE nrtype ,   ONLY : DP, I4b
    IMPLICIT NONE
    
    REAL(DP),SAVE,ALLOCATABLE,DIMENSION(:,:) :: enn,ennw,ennv,    &
              volsn


  END MODULE neutral_data
