       subroutine derivedmodlocal       
************************************************************************
cmnt   profiles of quantities derived from model profiles of te,ti,ne      
************************************************************************
       implicit none
c
       include '../inc/input.m'
       include '../inc/tport.m'
       include '../inc/model.m'
c
       integer j
       real*8 ve(0:jmaxm),vpar(0:jmaxm),vmode(0:jmaxm)
       real*8 vnewk3x, egeo_local, pgeo_local, rdrho_local
       real*8 rdrho_local_p1, vstar_sign, alpha_neo_hold
       real*8 fc, akappa1
c
      do j=jm, jm+2
       csda_m(j)=979.D3*(te_m(j)*1.D3)**.5D0/
     >    (arho_exp*100.D0)/(amassgas_exp)**0.5D0
cjek
       if (bt_flag .gt. 0) then
         rhosda_m(j)=((102.D0*(te_m(j)*1.D3)**.5D0)/bteff_exp(j)
     >    /1.D4)*(amassgas_exp)**.5D0/(arho_exp*100.D0)
       else
         rhosda_m(j)=((102.D0*(te_m(j)*1.D3)**.5D0)/bt_exp/1.D4)
     >    *(amassgas_exp)**.5D0/(arho_exp*100.D0)
       endif
      enddo
c
c   local gyrobohm unit of diffusion
c
       cgyrobohm_m(jm)=1.D-4*979.D3*(tem*1.D3)**.5D0/
     >  (arho_exp*100.D0)*(1.02e2*(tem*1.D3)**.5D0/bt_exp/
     >  1.D4)**2.D0*(amassgas_exp)**.5D0
c          
       betae_m(jm) = 400.D0*nem*tem/(1.D5*bt_exp**2.D0)
       betai_m(jm) = 400.D0*nim*tim/(1.D5*bt_exp**2.D0)
c        
crew    gks collisionality (xnu/w_star_i)*(ky*rho_i)
       vnewk3x=  
     >   0.117D0*nem*tem**(-1.5D0)/(tim**0.5D0)*(arho_exp)*
     >   (amassgas_exp/2.D0)**0.5D0
crew   as used in gks multiply by 1/2 and takout any 1/2 factor in solfp
crew          vnewk3x=vnewk3x/2.
       xnu_m(jm) =vnewk3x/(2.D0*tem/tim)**0.5D0
crew  10/25/95 fixed zeff+1 factor: zeff col with ions;1 col with elecs.
       xnu_m(jm) = xnu_m(jm)*(zeff_exp(jm)+zeff_e)
c
       vnewstare_m(jm)=zeff_exp(jm) *2.91D-6*nem*1.D13*15.D0/
     >  (tem*1.D3)**2.D0*rmaj_exp(jm)*100.D0*q_exp(jm)
     >     /(rmin_exp(jm)/rmaj_exp(jm)+1.D-10)**1.5D0/419.D5
       vnewstari_m(jm)=4.78e-8*nem*1.D13*15.D0/
     >  (tim*1.D3)**2.D0*rmaj_exp(jm)*100.D0*q_exp(jm)
     >     /(rmin_exp(jm)/rmaj_exp(jm)+1.D-10)**1.5D0/979.D3
c
       if (jm.le.jout_m-4.and.ishoot.ge.0) then
c
         j=jm+1
         zpte_m(j)=-(dlog(te_m(j-1))-dlog(te_m(j)))/(rho(j-1)-rho(j))
         zpti_m(j)=-(dlog(ti_m(j-1))-dlog(ti_m(j)))/(rho(j-1)-rho(j))
         zpne_m(j)=-(dlog(ne_m(j-1))-dlog(ne_m(j)))/(rho(j-1)-rho(j))
         zpni_m(j)=-(dlog(ni_m(j-1))-dlog(ni_m(j)))/(rho(j-1)-rho(j))
         if(i_dengrad.eq.0) zpni_m(j)=zpne_m(j)
c
         j=jm+2
         zpte_m(j)=-(dlog(te_m(j-1))-dlog(te_m(j)))/(rho(j-1)-rho(j))
         zpti_m(j)=-(dlog(ti_m(j-1))-dlog(ti_m(j)))/(rho(j-1)-rho(j))
         zpne_m(j)=-(dlog(ne_m(j-1))-dlog(ne_m(j)))/(rho(j-1)-rho(j))
         zpni_m(j)=-(dlog(ni_m(j-1))-dlog(ni_m(j)))/(rho(j-1)-rho(j))
         if(i_dengrad.eq.0) zpni_m(j)=zpne_m(j)    
c
c shift values at j=jm to jm+1
c note elongation dependence of alpha
c note zpne used as better represenation of total pressure than zpni
c 
        j=jm+1
       if(igeo_m.ge.0) then
c2/27/98        alpha_m(jm)=relx*alpha_m(jm)+(1.D0-relx)*(
c2/27/98     > sqrt(elong_exp(j))/((1.D0+elong_exp(j)**2.D0)/2.D0)*
        alpha_m(jm)=relx*alpha_m(jm)+(1.D0-relx)*(
     > drhodr(j)*
     > q_exp(j)**2.D0*rmaj_exp(j)/arho_exp*
     >    betae_m(j)*((ti_m(j)*ni_m(j)/te_m(j)/ne_m(j))*
     >    (zpni_m(j)+zpti_m(j))
     >    +zpne_m(j)+zpte_m(j))   )
c2/27/98     >    betae_m(j)*((ti_m(j)/te_m(j))*(zpne_m(j)+zpti_m(j))
c2/27/98     >    +zpne_m(j)+zpte_m(j)) )
       endif
       if(igeo_m.ge.4) then
        alpha_m(jm)=relx*alpha_m(jm)+(1.D0-relx)*(
     > geoalpha(j)*
     > q_exp(j)**2.D0*rmaj_exp(j)/arho_exp*
     >    betae_m(j)*((ti_m(j)*ni_m(j)/te_m(j)/ne_m(j))*
     >    (zpni_m(j)+zpti_m(j))
     >    +zpne_m(j)+zpte_m(j))   )
c2/27/98     >    betae_m(j)*((ti_m(j)/te_m(j))*(zpne_m(j)+zpti_m(j))
c2/27/98     >    +zpne_m(j)+zpte_m(j))  )
       endif
       if(igeo_m.eq.-1) then
        alpha_m(jm)=relx*alpha_m(jm)+(1.D0-relx)*(
     > q_exp(j)**2.D0*rmajor_exp/arho_exp*
     > sqrt(elonga_exp)/((1.D0+elonga_exp**2.D0)/2.D0)*
     >    betae_m(j)*((ti_m(j)*ni_m(j)/te_m(j)/ne_m(j))*
     >    (zpni_m(j)+zpti_m(j))
     >    +zpne_m(j)+zpte_m(j))   )
c2/27/98     >    betae_m(j)*((ti_m(j)/te_m(j))*(zpne_m(j)+zpti_m(j))
c2/27/98     >    +zpne_m(j)+zpte_m(j)) )
c
       endif
c
        egeo_local=1.D0
        pgeo_local=1.D0
        rdrho_local=1.D0
        if(igeo_m.ge.5) egeo_local=georotrate(jm)
        if(igeo_m.ge.3) then
          pgeo_local=drhodr(jm)
          rdrho_local=rmin_exp(jm)/arho_exp/rho(jm)
          rdrho_local_p1=rmin_exp(jm+1)/arho_exp/rho(jm+1)
        endif
c
c vstarp_m is diamagnetic part of egamma_m (doppler shear rate)
c vstar_sign is negative for negative vstar_i. Thus for co-injection 
c or positive angrot toroidal rotation cancels the diamgnetic rotation
c
       vstar_sign=-1.*alpha_dia
c
        j=jm+1
        egeo_local=1.D0
        pgeo_local=1.D0
        rdrho_local=1.D0
        rdrho_local_p1=1.D0
        if(igeo_m.ge.5) egeo_local=georotrate(j)
        if(igeo_m.ge.3) then
          pgeo_local=drhodr(j)
          rdrho_local=rmin_exp(j)/arho_exp/rho(j)
          rdrho_local_p1=rmin_exp(j+1)/arho_exp/rho(j+1)
        endif   
c
        vstarp_m(jm)=relx*vstarp_m(jm)+(1.D0-relx)*(
     > +egeo_local*vstar_sign*
     > (rho(j+1)*rdrho_local_p1+rho(j)*rdrho_local)/2.D0*(
     >(ti_m(j+1)/te_m(j+1))*csda_m(j+1)*(zpti_m(j+1)+zpni_m(j+1))
     >  *pgeo_local*rhosda_m(j+1)/rho(j+1)/rdrho_local_p1
     >-(ti_m(j)/te_m(j))*csda_m(j)*(zpti_m(j)+zpni_m(j))
     >   *pgeo_local*rhosda_m(j)/rho(j)/rdrho_local
     >  )/(rho(j+1)-rho(j))/csda_m(j)
     >  )
c
c banana regime ie collisionless limit formulas
        alpha_neo_hold=alpha_neo
        j=jm+1
        fc=1-1.46D0*(rmin_exp(j)/rmaj_exp(j))**0.5D0+
     >      0.46D0*(rmin_exp(j)/rmaj_exp(j))**1.5D0
        akappa1=0.8839D0*fc/(0.3477+0.4058*fc)
        if(irot1.eq.1) alpha_neo=-akappa1+1.D0
        j=jm+1
        egeo_local=1.D0
        pgeo_local=1.D0
        rdrho_local=1.D0
        rdrho_local_p1=1.D0
        if(igeo_m.ge.5) egeo_local=georotrate(j)
        if(igeo_m.ge.3) then
          pgeo_local=drhodr(j)
          rdrho_local=rmin_exp(j)/arho_exp/rho(j)
          rdrho_local_p1=rmin_exp(j+1)/arho_exp/rho(j+1)
        endif   
c
        ve(j)=-(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
     > (zpni_m(j)+alpha_neo*zpti_m(j))*vstar_sign*pgeo_local
     > -rho(j)*rdrho_local*
     >  arho_exp/rmajor_exp/q_exp(j)*rmajor_exp*angrotp_exp(j)
        vpar(j)=rmajor_exp*angrotp_exp(j)-vstar_sign*
     > (ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*pgeo_local*
     >((alpha_neo-1.D0)*zpti_m(j))*rho(j)*rdrho_local*
     > arho_exp/rmajor_exp/q_exp(j)
        vmode(j)=anfreq_m(j)/(ky_j(j)+1.D-10)*
     >    csda_m(j)*arho_exp*rhosda_m(j)
c
        vexb_m(j)=ve(j)
        vmode_m(j)=vmode(j)
        vstar_m(j)=-(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
     > (zpni_m(j)+zpti_m(j))*vstar_sign*pgeo_local
c
        j=jm+2
        egeo_local=1.D0
        pgeo_local=1.D0
        rdrho_local=1.D0
        rdrho_local_p1=1.D0
        if(igeo_m.ge.5) egeo_local=georotrate(j)
        if(igeo_m.ge.3) then
          pgeo_local=drhodr(j)
          rdrho_local=rmin_exp(j)/arho_exp/rho(j)
          rdrho_local_p1=rmin_exp(j+1)/arho_exp/rho(j+1)
        endif   
c
        ve(j)=-(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
     > (zpni_m(j)+alpha_neo*zpti_m(j))*vstar_sign*pgeo_local
     > -rho(j)*rdrho_local*
     >  arho_exp/rmajor_exp/q_exp(j)*rmajor_exp*angrotp_exp(j)
        vpar(j)=rmajor_exp*angrotp_exp(j)-vstar_sign*
     > (ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*pgeo_local*
     >((alpha_neo-1.D0)*zpti_m(j))*rho(j)*rdrho_local*
     > arho_exp/rmajor_exp/q_exp(j)
        vmode(j)=anfreq_m(j)/(ky_j(j)+1.D-10)*
     >    csda_m(j)*arho_exp*rhosda_m(j)
c
        j=jm+1
        egeo_local=1.D0
        pgeo_local=1.D0
        rdrho_local=1.D0
        rdrho_local_p1=1.D0
        if(igeo_m.ge.5) egeo_local=georotrate(j)
        if(igeo_m.ge.3) then
          pgeo_local=drhodr(j)
          rdrho_local=rmin_exp(j)/arho_exp/rho(j)
          rdrho_local_p1=rmin_exp(j+1)/arho_exp/rho(j+1)
        endif 
      if(igeo_m.le.2) then
        egamma_m(jm)=relx*egamma_m(jm)+(1.D0-relx)*(
     >  (rho(j)+rho(j+1))/2.D0*(ve(j+1)/rho(j+1)-
     >  ve(j)/rho(j))/(rho(j+1)-rho(j))/arho_exp/csda_m(j)
     >   -(vpar(j+1)/(rmajor_exp*q_exp(j+1)/shat_exp(j+1))
     >   +vpar(j)/(rmajor_exp*q_exp(j)/shat_exp(j)))/2.D0/csda_m(j)
     >  )
c
        gamma_mode_m(jm)=relx*gamma_mode_m(jm)+(1.D0-relx)*(
     >  (rho(j)+rho(j+1))/2.D0*(vmode(j+1)/rho(j+1)-
     >  vmode(j)/rho(j))/(rho(j+1)-rho(j))/arho_exp/csda_m(j)
     >  )
c
        gamma_p_m(jm)=relx*gamma_p_m(jm)+(1.D0-relx)*(
     >   -(vpar(j+1)-vpar(j))/(rho(j+1)-rho(j))/arho_exp
     >   /csda_m(j)
     >  )
      endif
      if(igeo_m.ge.3) then
        egamma_m(jm)=relx*egamma_m(jm)+(1.D0-relx)*(
     >   egeo_local*drhodrrrho(j)*
     >  (rho(j)+rho(j+1))/(q_exp(j)+q_exp(j+1))*
     >  (ve(j+1)*q_exp(j+1)/rho(j+1)/rdrho_local_p1-
     >  ve(j)*q_exp(j)/rho(j)/rdrho_local)/
     > (rho(j+1)-rho(j))/arho_exp/csda_m(j)
     >  )
c
        gamma_mode_m(jm)=relx*gamma_mode_m(jm)+(1.D0-relx)*(
     >  egeo_local*drhodrrrho(j)*
     >  (rho(j)+rho(j+1))/2.D0*(vmode(j+1)/rho(j+1)/rdrho_local_p1-
     >  vmode(j)/rho(j)/rdrho_local)/
     >  (rho(j+1)-rho(j))/arho_exp/csda_m(j)
     >  )
c
        gamma_p_m(jm)=relx*gamma_p_m(jm)+(1.D0-relx)*(
     >   -drhodr(j)*
     >   (vpar(j+1)-vpar(j))/(rho(j+1)-rho(j))/arho_exp
     >   /csda_m(j)
     >  )
      endif
      if(igeo_m.ge.6) then
        egamma_m(jm)=relx*egamma_m(jm)+(1.D0-relx)*(
     >   -georotrate(j)*
     >  (rmin_exp(j)+rmin_exp(j+1))/
     >  (q_exp(j)+q_exp(j+1))*
     >  (angrot_exp(j+1)-angrot_exp(j))/
     >  (rmin_exp(j+1)-rmin_exp(j))/csda_m(j)
     >  )
c
        vstarp_m(jm)=relx*vstarp_m(jm)+(1.D0-relx)*(
     >  -corot*
     >  (zpti_m(j)+zpne_m(j))*
     >  drhodr(j)*zpti_m(j)*
     >  (ti_m(j)/te_m(j))*rhosda_m(j)
     >  )
      endif 
c
        alpha_neo=alpha_neo_hold
c
       endif
c
       return
       end
