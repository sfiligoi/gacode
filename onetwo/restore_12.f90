  MODULE restore_12
    USE param,  ONLY :   kj,kion
    USE nrtype, ONLY :   I4B ,DP
    IMPLICIT NONE
    
    REAL(DP), DIMENSION(kj)  ::                                         &
                        storque_nub,storqueb_nubr,angrot_nub,           &
                        te_nub,ti_nub,rbp_nub,ene_nub,                  &
                        curden_nub,etor_nub,curdri_nub,curohm_nub,      &
                        currf_nub,curboot_nub,fcap_nub,hcap_nub,        &
                        gcap_nub

    REAL(DP), DIMENSION(kj,kion)  ::  en_nub  


   END MODULE restore_12
