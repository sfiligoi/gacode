       subroutine derivedmodprofiles       
************************************************************************
c 1/22/01 jek added +1.D-10 to zpne_m, zpni_m
cmnt   profiles of quantities derived from model profiles
cmnt   exp values substituted for model zeros      
************************************************************************
       implicit none
c
       include '../inc/input.m'
       include '../inc/tport.m'
       include '../inc/model.m'
c
       integer j, jin_p
       real*8 tem_p, tem_m, tem_b, tim_p, tim_m, tim_b
       real*8 nem_p, nem_m, nim_p, nim_m, nzm_p, nzm_m
       real*8 ve(0:jmaxm),vpar(0:jmaxm)
       real*8 eps, wstre, wstri, destr, distr, voltot
       real*8 aveniti, aveniti2, vnewk3, aveti, avete, dr
       real*8 dvoldr_p, dvoldr_m, destr_st, distr_st
       real*8 wstre_core_m,  wstri_core_m
c
       eps = 1.D-6
c   
c model confinement times and physics quantities for whole profile
c
        taue_m(0)=0.D0
        taui_m(0)=0.D0
        taupe_m(0)=0.D0     !taux(0)=infinity...set scale here
        taupi_m(0)=0.D0
c
        wstre=0.D0
        wstri=0.D0
        wstr_m=0.D0
        wstr_inc_m=0.D0
        wstre_core_m=0.D0
        wstri_core_m=0.D0
        wstr_core_m=0.D0
        destr=0.D0
        distr=0.D0
        voltot=0.D0
        aveti=0.D0
        avete=0.D0
        aveniti=0.D0
        aveniti2=0.D0
c           
        do j=1,jmaxm
         tem_p=te_exp(j)
         tem_m=te_exp(j-1)
         tim_p=ti_exp(j)
         tim_m=ti_exp(j-1)
         nem_p=ne_exp(j)
         nem_m=ne_exp(j-1)
         nim_p=ni_exp(j)
         nim_m=ni_exp(j-1)
         nzm_m=nz_exp(j-1)
c         
         if (te_m(j).gt.0.) tem_p=te_m(j)
         if (te_m(j-1).gt.0.) tem_m=te_m(j-1)
         if (ti_m(j).gt.0.) tim_p=ti_m(j)
         if (ti_m(j-1).gt.0.) tim_m=ti_m(j-1)
         if (ne_m(j).gt.0.) nem_p=ne_m(j)
         if (ne_m(j-1).gt.0.) nem_m=ne_m(j-1)
         if (ni_m(j).gt.0.) nim_p=ni_m(j)
         if (ni_m(j-1).gt.0.) nim_m=ni_m(j-1)
         if (nz_m(j).gt.0.) nzm_p=nz_m(j)
         if (ni_m(j-1).gt.0.) nzm_m=nz_m(j-1)
c         
c         zpte_m(j)=-(dlog(tem_m)-dlog(tem_p))/(rho(j-1)-rho(j))
c         zpti_m(j)=-(dlog(tim_m)-dlog(tim_p))/(rho(j-1)-rho(j))
c         zpne_m(j)=-(dlog(nem_m)-dlog(nem_p))/(rho(j-1)-rho(j))
c         zpni_m(j)=-(dlog(nim_m)-dlog(nim_p))/(rho(j-1)-rho(j))
c         zpnz_m(j)=-(dlog(nzm_m)-dlog(nzm_p))/(rho(j-1)-rho(j))
c         zpnitot_m(j)=-(dlog(nim_m+nzm_m)-dlog(nim_p+nzm_p))/
c     >                (rho(j-1)-rho(j))
c         if(i_dengrad.eq.0) then
c           zpni_m(j)=zpne_m(j)
c           zpnitot_m(j)=zpne_m(j)
c         endif
c         if(jbk_m.ne.0.and.j.ne.jmaxm) then
c          tem_b=te_exp(j+1)
c          tim_b=ti_exp(j+1)
c          if(te_m(j).gt.0) tem_b=te_m(j+1)
c          if(ti_m(j).gt.0) tim_b=ti_m(j+1)
c          zpte_m(j)=-(dlog(tem_m)-dlog(tem_b))/(rho(j-1)-rho(j+1))
c          zpti_m(j)=-(dlog(tim_m)-dlog(tim_b))/(rho(j-1)-rho(j+1))
c         endif
c  
         tem=tem_p
         tim=tim_p
         nem=nem_p
         nim=nim_p
         nzm=nzm_p
         if (DABS(tim).lt.eps) tim=eps
         if (DABS(tem).lt.eps) tem=eps
c         
cmnt recompute for whole profile
c
c      cgyrobohm_m(j)=1.D-4*979.D3*dsqrt(tem*1.D3)/
c    >  (arho_exp*100.D0)*(102.D0*dsqrt(tem*1.D3)/bt_exp/
c    >  1.D4)**2.D0*dsqrt(amassgas_exp)
c
       csda_m(j)=979.D3*dsqrt(tem*1.D3)/(arho_exp*100.D0)
     >    /dsqrt(amassgas_exp)
cjek
       if (bt_flag .gt. 0) then
         rhosda_m(j)=((102.D0*dsqrt(tem*1.D3))/bteff_exp(j)/1.D4)
     >    *dsqrt(amassgas_exp)/(arho_exp*100.D0)
       else
         rhosda_m(j)=((102.D0*dsqrt(tem*1.D3))/bt_exp/1.D4)
     >    *dsqrt(amassgas_exp)/(arho_exp*100.D0)
       endif
c            
       betae_m(j) = 400.D0*nem*tem/(1.D5*bt_exp**2.D0)
       betai_m(j) = 400.D0*nim*tim/(1.D5*bt_exp**2.D0)
c         
crew    gks collisionality (xnu/w_star_i)*(ky*rho_i)
c       vnewk3=  
c     >   0.117D0*nem*tem**(-1.5D0)/(tim**0.5D0)*(arho_exp)*
c     >   (amassgas_exp/2.D0)**0.5D0
c       xnu_m(j) =vnewk3/(2.D0*tem/tim)**0.5D0
c       xnu_m(j)=xnu_m(j)*(zeff_exp(j)+zeff_e)
c       
       vnewstare_m(j)=zeff_exp(j)*2.91D-6*nem*1.D13*15.D0/
     >  (tem*1.D3)**2.D0*rmaj_exp(j)*100.D0*q_exp(j)
     >     /(rmin_exp(j)/rmaj_exp(j)+1.D-10)**1.5D0/419.D5
       vnewstari_m(j)=4.78D-8*nem*1.D13*15.D0/
     >  (tim*1.D3)**2.D0*rmaj_exp(j)*100.D0*q_exp(j)
     >     /(rmin_exp(j)/rmaj_exp(j)+1.D-10)**1.5D0/979.D3
c     
         dr=(rho(j)-rho(j-1))*arho_exp
         dvoldr_p=sfactor(j)
         dvoldr_m=sfactor(j-1)
c
c... Total stored energies
c
         wstre=wstre+kevdsecpmw*1.D19*
     >       0.5D0*(dvoldr_p*1.5D0*nem_p*tem_p+
     >       dvoldr_m*1.5D0*nem_m*tem_m)*dr
c     
         taue_m(j)=wstre/(powe_exp(j)-pow_ei_cor_m(j)
     >   +xfus_m*powe_fus_cor_m(j)-xbr_m*pow_br_cor_m(j))
c        
         wstri=wstri+kevdsecpmw*1.D19*
     >       0.5D0*(dvoldr_p*1.5D0*nim_p*tim_p+
     >            dvoldr_m*1.5D0*nim_m*tim_m)*dr
c
c...  Core stored energies inside jout_m (rho_bc)
c
         if (j.le.jout_m) then
           wstre_core_m = wstre_core_m + kevdsecpmw*1.D19*
     >       0.5D0*(dvoldr_p*1.5D0*nem_p*tem_p+
     >       dvoldr_m*1.5D0*nem_m*tem_m)*dr
           wstri_core_m = wstri_core_m + kevdsecpmw*1.D19*
     >       0.5D0*(dvoldr_p*1.5D0*nim_p*tim_p+
     >       dvoldr_m*1.5D0*nim_m*tim_m)*dr
         endif
c
         taui_m(j)=wstri/(powi_exp(j)+pow_ei_cor_m(j)
     >    +xfus_m*powi_fus_cor_m(j))
c       
         taut_m(j)=(wstre+wstri)/(powe_exp(j)+powi_exp(j)
     >    +xfus_m*powi_fus_cor_m(j)+xfus_m*powe_fus_cor_m(j)-
     >    xbr_m*pow_br_cor_m(j))
c       
         destr=destr+kevdsecpmw*1.D19*
     >       0.5D0*(dvoldr_p*nem_p+dvoldr_m*nem_m)*dr
c    
         taupe_m(j)=destr/(flow_exp(j)+1.D-10)
c
         distr=distr+kevdsecpmw*1.D19*
     >       0.5D0*(dvoldr_p*nim_p+
     >            dvoldr_m*nim_m)*dr
c
         taupi_m(j)=distr/(flow_exp(j)+1.D-10)
c 
         voltot=voltot+0.5D0*(dvoldr_p+dvoldr_m)*dr
         aveniti=aveniti+
     >        0.5D0*(dvoldr_p*nim_p*tim_p+
     >             dvoldr_m*nim_m*tim_m)*dr
         aveniti2=aveniti2+
     >        0.5D0*(dvoldr_p*nim_p**2.D0*tim_p**2.D0+
     >             dvoldr_m*nim_m**2.D0*tim_m**2.D0)*dr
         aveti=aveti+0.5D0*(dvoldr_p*tim_p+
     >         dvoldr_m*tim_m)*dr
         avete=avete+0.5D0*(dvoldr_p*tem_p+
     >         dvoldr_m*tem_m)*dr
c         
        enddo
c 
        wstre_m=wstre
        wstri_m=wstri
        wstr_m=wstre+wstri
        wstr_inc_m=wstr_m-wstr_ped_exp
        wstr_core_m=wstre_core_m+wstri_core_m
c        volaveniti_m=aveniti2/aveniti
        gtaut_m=taut_m(jmaxm)
        gtaui_m=taui_m(jmaxm)
        n_t_tau_m=1.D19*volaveniti_m*gtaut_m
        n0_t0_tau_m=
     >       1.D19*ni_m(3*jmaxm/10)*ti_m(3*jmaxm/10)*gtaut_m
        jin_p=jmaxm/5
        te20_m=te_m(jin_p)
        ti20_m=ti_m(jin_p)
        volaveniti_m=aveniti/voltot
        volaveniti2_m=aveniti2/voltot
        volaveti_m=aveti/voltot
        volavete_m=avete/voltot
c
c      do j=1,jmaxm
c         zpti2_m(j)=(ti_m(j+1)+ti_m(j-1)-2.D0*ti_m(j))/
c     >              (rho(j)-rho(j-1))**2.D0
c         write(*,55) j, rho(j), ti_exp(j), zpti2_m(j)
c      enddo
c
 55   format(i3,2x,0p1f6.4,0p6e13.5)
c
       return
       end
