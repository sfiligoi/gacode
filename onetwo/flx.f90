 
  MODULE flx
      USE param, ONLY: kk,kjm1
      IMPLICIT NONE 
      REAL *8 ,DIMENSION(:,:)  :: flux(kk,kjm1)
      REAL *8 ,DIMENSION(:)    :: fluxe(kjm1), fluxi(kjm1), qieneo(kjm1) 
      REAL *8 ,DIMENSION(:)    :: anal_eng_flux_e(kjm1),anal_eng_flux_i(kjm1)
 
  END MODULE flx

