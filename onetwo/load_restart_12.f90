 SUBROUTINE load_restart_12

       USE tordlrot,    ONLY : storque,storqueb,angrot
       USE soln,        ONLY : te,ti,rbp,ene,en,curden,etor,curden
       USE sourc,       ONLY : curdri,curohm,currf,curboot
       USE geom,        ONLY : fcap,gcap,hcap
       USE restore_12 
       USE numbrs,  ONLY : nj,nion
       IMPLICIT NONE

      
                te_nub(1:nj)                  = te(1:nj)
                ti_nub(1:nj)                  = ti(1:nj)
                ene_nub(1:nj)                 = ene(1:nj)
                en_nub(1:nj,1:nion)           = en(1:nj,1:nion)
                etor_nub(1:nj)                = etor(1:nj)
                rbp_nub(1:nj)                 = rbp(1:nj)
                curohm_nub(1:nj)              = curohm(1:nj)
                curden_nub(1:nj)              = curden(1:nj)
                currf_nub(1:nj)               = currf(1:nj)
                curboot_nub(1:nj)             = curboot(1:nj)
                curdri_nub(1:nj)              = curdri(1:nj)
                angrot_nub(1:nj)              = angrot(1:nj)
                storque_nub(1:nj)             = storque(1:nj)
                storqueb_nubr(1:nj)           = storqueb(1:nj)
                hcap_nub(1:nj)                = hcap(1:nj)
                gcap_nub(1:nj)                = gcap(1:nj)
                fcap_nub(1:nj)                = fcap(1:nj)
        RETURN
 END SUBROUTINE load_restart_12
