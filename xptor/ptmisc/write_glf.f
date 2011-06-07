      subroutine write_glf(zdate,j,j1,j2)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c     This routine writes out results in namelist for reading into
c     the stand-alone code for the GLF23 model
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      implicit none
c
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/ptor.m'
      include '../inc/glf.m'
c
      character*12 zdate
      integer j, j1, j2, jmx
c
        open (26,file='glf.out',status='unknown')
c
        jmx=j2-j1
        write(26,*) '&nlglf'
        write(26,55) zdate, tok, shot
        write(26,*) 'leigen = 1,'
        write(26,395) nroot_gf
        write(26,390) bt_flag
        write(26,*) 'lprint = 1,'
        write(26,*) 'jmm = 0,'
        write(26,400) jmx
        write(26,*) 'jshoot = 0,'
        write(26,*) 'igrad = 0,'
        write(26,*) 'idengrad = 2,'
        write(26,*) 'itport_pt(1) = 0,'
        write(26,*) 'itport_pt(2) = 1,'
        write(26,*) 'itport_pt(3) = 1,'
        write(26,*) 'itport_pt(4) = 1,'
        write(26,*) 'itport_pt(5) = 0,'
        write(26,401) (rho(j),j=j1,j2)
        write(26,402) (ne_m(j),j=j1,j2)
        write(26,403) (ni_m(j),j=j1,j2)
        write(26,404) (ns_m(j),j=j1,j2)
        write(26,405) (te_m(j),j=j1,j2)
        write(26,406) (ti_m(j),j=j1,j2)
        write(26,407) (zeff_exp(j),j=j1,j2)
        write(26,408) (angrotp_exp(j),j=j1,j2)
        write(26,409) (egamma_m(j),j=j1,j2)
        write(26,410) (gamma_p_exp(j),j=j1,j2)
        write(26,411) (vphi_m(j),j=j1,j2)
        write(26,412) (vpar_m(j),j=j1,j2)
        write(26,413) (vper_m(j),j=j1,j2)
        write(26,450) (gamma_p_m(j),j=j1,j2)
        write(26,414) bt_exp
        write(26,415) arho_exp
        write(26,416) (gradrho_exp(j),j=j1,j2)
        write(26,417) (gradrhosq_exp(j),j=j1,j2)
        write(26,418) (rmin_exp(j),j=j1,j2)
        write(26,419) (rmaj_exp(j),j=j1,j2)
        write(26,420) rmajor_exp
        write(26,421) (q_exp(j),j=j1,j2)
        write(26,422) (shat_exp(j),j=j1,j2)
        write(26,423) (alpha_exp(j),j=j1,j2)
        write(26,424) (elong_exp(j),j=j1,j2)
        write(26,425) amassgas_exp
        write(26,426) zimp_exp
        write(26,427) amassimp_exp
        write(26,428) alpha_e_gf
        write(26,429) xalpha
        write(26,*) '/'
        write(26,*) 'Results from XPTOR:'
        write(26,430) (rho(j),j=j1,j2)
        write(26,431) (chiigb_glf(j),j=j1,j2)
        write(26,432) (chiegb_glf(j),j=j1,j2)
        write(26,433) (etagb_phi_glf(j),j=j1,j2)
        write(26,434) (anrate_m(j),j=j1,j2)
        write(26,455) (anfreq_m(j),j=j1,j2)
        write(26,435) (egamma_m(j),j=j1,j2)
        write(26,436) (zpmti_glf(j),j=j1,j2)
        write(26,437) (zpmte_glf(j),j=j1,j2)
        write(26,438) (zpmne_glf(j),j=j1,j2)
        write(26,439) (zpmni_glf(j),j=j1,j2)
        write(26,440) (bteff_exp(j),j=j1,j2)
        close(26)
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
 55   format(' !',2x,a10,2x,a4,a6)
 390  format(' bt_flag = ',i2,',')
 395  format(' nroot = ',i2,',')
 400  format(' jmaxm = ',i2,',')
 401  format(' rho  = ',0p301f10.5)
 402  format(' ne_m = ',0p301f10.6)
 403  format(' ni_m = ',0p301f10.6)
 404  format(' ns_m = ',0p301f10.6)
 405  format(' te_m = ',0p301f10.6)
 406  format(' ti_m = ',0p301f10.6)
 407  format(' zeff_exp    = ',0p301f10.5)
 408  format(' angrotp_exp = ',1p51e13.5)
 409  format(' egamma_exp  = ',1p51e13.5)
 410  format(' gamma_p_exp = ',1p51e13.5)
 411  format(' vphi_m = ',1p51e13.5)
 412  format(' vpar_m = ',1p51e13.5)
 413  format(' vper_m = ',1p51e13.5)
 414  format(' bt_exp   = ',1p51e13.5)
 415  format(' arho_exp = ',1p51e13.5)
 416  format(' gradrho_exp   = ',1p51e13.5)
 417  format(' gradrhosq_exp = ',1p51e13.5)
 418  format(' rmin_exp   = ',1p51e13.5)
 419  format(' rmaj_exp   = ',1p51e13.5)
 420  format(' rmajor_exp = ',1p51e13.5)
 421  format(' q_exp      = ',1p51e13.5)
 422  format(' shat_exp   = ',1p51e13.5)
 423  format(' alpha_exp  = ',1p51e13.5)
 424  format(' elong_exp  = ',1p51e13.5)
 425  format(' amassgas_exp  = ',1p51e13.5)
 426  format(' zimp_exp      = ',1p51e13.5)
 427  format(' amassimp_exp  = ',1p51e13.5)
 428  format(' alpha_e  = ',f10.5)
 429  format(' x_alpha  = ',f10.5)
 430  format('  rho     = ',0p301e13.5)
 431  format('  chii    = ',1p51e13.5,' (m^2/s)')
 432  format('  chie    = ',1p51e13.5,' (m^2/s)')
 433  format('  eta-phi = ',1p51e13.5,' (m^2/s)')
 434  format('  anrate  = ',1p51e13.5)
 435  format('  egamma  = ',1p51e13.5)
 436  format('  zpmti  = ',1p51e13.5)
 437  format('  zpmte  = ',1p51e13.5)
 438  format('  zpmne  = ',1p51e13.5)
 439  format('  zpmni  = ',1p51e13.5)
 440  format('  bteff  = ',1p51e13.5)
 450  format(' gamma_p_m = ',1p51e13.5)
 455  format(' anfreq = ',1p51e13.5)
c
      return
      end
