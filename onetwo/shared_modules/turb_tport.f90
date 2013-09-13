   MODULE turb_tport  

     USE nrtype,                                     ONLY : Dp,I4B

     USE solcon_gcnmp,                               ONLY : itran_max

     IMPLICIT NONE
     REAL(DP)  wturbd(itran_max,itran_max), wturbv(itran_max,itran_max)

     REAL(DP),ALLOCATABLE,DIMENSION(:,:,:)  :: turb_stab

     REAL(DP) stab_mult,gmax,umax,gmax_crit

     INTEGER(I4B) jgmax,kgmax, jumax,kumax

   END  MODULE turb_tport  
