       subroutine derivedmodlocal_ptor
************************************************************************
c 1/18/00 jek added egamma_vphi, egamma_vpol, egamma_vstar
c 1/17/00 jek added corot to calculation in ve(j)
cmnt   profiles of quantities derived from model profiles of te,ti,ne
cmnt   revised derivedmodlocal to center between jm and jm+1 for ptor
************************************************************************
       implicit none
c
       include '../inc/input.m'
       include '../inc/tport.m'
       include '../inc/model.m'
       include '../inc/ptor.m'
c
       integer j, jptf, jpt, jptt, jj
       real*8 ve(0:jmaxm),vpar(0:jmaxm),vmode(0:jmaxm)
       real*8 vnewk3x, zpte_m_hold, zpti_m_hold, zpne_m_hold
       real*8 zpni_m_hold, vstar_sign, egeo_local, pgeo_local
       real*8 rdrho_local, rdrho_local_p1, alpha_neo_hold
       real*8 fc, akappa1
c
      if(jm.eq.0.or.jm.eq.jmaxm) then
       write(6,*) 'can not call derivedmodlocal_ptor for this jm'
      endif
c
      tem=(te_m(jm+1)+te_m(jm))/2.D0
      tim=(ti_m(jm+1)+ti_m(jm))/2.D0
      nem=(ne_m(jm+1)+ne_m(jm))/2.D0
      nim=(ni_m(jm+1)+ni_m(jm))/2.D0
c
      jptf=iparam_pt(8)
       if(jm.eq.1) jptf=0
      jpt=1+iparam_pt(7)
      jptt=jpt
       if(jptt.lt.1) jptt=1
       if(jm.eq.jmaxm-1) jptt=1
      jpt=jptt
c
c  note: if iparam_pt(8)=0 and iparam_pt(7)=0; then the
c  dependence of shear's  on zpxx_m(jm-1) and zpxx_m(jm+1)
c  and alpha(jm) depends on zpxx_m(jm+1)
c
      do j=jm-jptf,jm+jpt
c   local rate unit
       csda_m(j)=979.D3*(te_m(j)*1.D3)**.5D0/
     >    (arho_exp*100.D0)/(amassgas_exp)**0.5D0
      enddo
      do j=jm-jptf,jm+jpt
c   local rho_star
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
       cgyrobohm_m(jm)=1.D-4*979.D3*(tem*1.D3)**.5D0/
     >  (arho_exp*100.D0)*(102.D0*(tem*1.D3)**.5D0/
     >  bt_exp/1.D4)**2.D0*(amassgas_exp)**.5D0
c
       betae_m(jm) = 400.D0*nem*tem/(1.D5*bt_exp**2.D0)
       betai_m(jm) = 400.D0*nim*tim/(1.D5*bt_exp**2.D0)
c
crew    gks collisionality (xnu/w_star_i)*(ky*rho_i)
       vnewk3x=
     >   0.117D0*nem*tem**(-1.5D0)/(tim**0.5D0)*(arho_exp)*
     >   (amassgas_exp/2.D0)**0.5D0
crew   as used in gks multiply by 1/2 and takout any 1/2 factor in solfp
crew          vnewk3x=vnewk3x/2.D0
       xnu_m(jm) =vnewk3x/(2.D0*tem/tim)**0.5D0
crew  10/25/95 fixed zeff+1 factor: zeff col with ions;1 col with elecs.
       xnu_m(jm) = xnu_m(jm)*(zeff_exp(jm)+zeff_e)
c
       vnewstare_m(jm)=zeff_exp(jm) *2.91D-6*nem*1.D13*15.D0/
     >  (tem*1.D3)**2.D0*rmaj_exp(jm)*100.D0*q_exp(jm)
     >     /(rmin_exp(jm)/rmaj_exp(jm)+1.D-10)**1.5D0/419.D3
       vnewstari_m(jm)=4.78D-8*nem*1.D13*15.D0/
     >  (tim*1.D3)**2.D0*rmaj_exp(jm)*100.D0*q_exp(jm)
     >     /(rmin_exp(jm)/rmaj_exp(jm)+1.D-10)**1.5D0/979.D3
c
      do j=jm-jptf, jm+jpt  
         zpte_m(j)=-(dlog(te_m(j-1))-dlog(te_m(j)))/(rho(j-1)-rho(j))
         zpti_m(j)=-(dlog(ti_m(j-1))-dlog(ti_m(j)))/(rho(j-1)-rho(j))
         zpne_m(j)=-(dlog(ne_m(j-1))-dlog(ne_m(j)))/(rho(j-1)-rho(j))
         zpni_m(j)=-(dlog(ni_m(j-1))-dlog(ni_m(j)))/(rho(j-1)-rho(j))
         if(i_dengrad.eq.0) zpni_m(j)=zpne_m(j) 
      enddo
c
        zpte_m_hold=zpte_m(jm+1)
        zpti_m_hold=zpti_m(jm+1)
        zpne_m_hold=zpne_m(jm+1)
        zpni_m_hold=zpni_m(jm+1)
c
        if(jptf.eq.0.and.iparam_pt(6).eq.2.D0) then
         zpte_m(jm+1)=zpmte
         zpti_m(jm+1)=zpmti
         zpne_m(jm+1)=zpmne
         zpni_m(jm+1)=zpmni
        endif
c
       if(igeo_m.ge.0) then
        alpha_m(jm)=relx*alpha_m(jm)+(1.D0-relx)*(
     > drhodr(jm)*
     > q_exp(jm)**2.D0*rmaj_exp(jm)/arho_exp*
     >    betae_m(jm)*((tim*nim/tem/nem)*
     >    (zpni_m(jm+1)+zpti_m(jm+1))
     >    +zpne_m(jm+1)+zpte_m(jm+1))   )
       endif
       if(igeo_m.ge.4) then
        alpha_m(jm)=relx*alpha_m(jm)+(1.D0-relx)*(
     > geoalpha(j)*
     > q_exp(j)**2.D0*rmaj_exp(j)/arho_exp*
     >    betae_m(jm)*((tim*nim/tem/nem)*
     >    (zpni_m(jm+1)+zpti_m(jm+1))
     >    +zpne_m(jm+1)+zpte_m(jm+1))   )
       endif
       if(igeo_m.eq.-1) then
        alpha_m(jm)=relx*alpha_m(jm)+(1.D0-relx)*(
     > q_exp(j)**2.D0*rmajor_exp/arho_exp*
     > sqrt(elonga_exp)/((1.D0+elonga_exp**2.D0)/2.D0)*
     >    betae_m(jm)*((tim*nim/tem/nem)*
     >    (zpni_m(jm+1)+zpti_m(jm+1))
     >    +zpne_m(jm+1)+zpte_m(jm+1))   )
       endif
c       write(*,*) 'alpha_m(jm) = ',alpha_m(jm)
c
c vstarp_m is diamagnetic part of egamma_m (doppler shear rate)
c vstar_sign is negative for negative vstar_i. Thus for co-injection or positive
c angrot toroidal rotation cancels the diamgnetic rotation
c
       vstar_sign=-1.D0*alpha_dia
c
        j=jm
        egeo_local=1.D0
        pgeo_local=1.D0
        rdrho_local=1.D0
        rdrho_local_p1=1.D0
        if(igeo_m.ge.5) egeo_local=georotrate(j)
        if(igeo_m.ge.3) then
          pgeo_local=drhodr(j)
          rdrho_local=rmin_exp(j-jptf)/arho_exp/rho(j-jptf)
          rdrho_local_p1=rmin_exp(j+jpt)/arho_exp/rho(j+jpt)
        endif
c
        vstarp_m(jm)=relx*vstarp_m(jm)+(1.D0-relx)*(
     > +egeo_local*vstar_sign*
     > (rho(j+jpt)*rdrho_local_p1+rho(j-jptf)*rdrho_local)/2.D0*(
     >(ti_m(j+jpt)/te_m(j+jpt))*csda_m(j+jpt)*
     > (zpti_m(j+jpt)+zpni_m(j+jpt))
     >  *pgeo_local*rhosda_m(j+jpt)/rho(j+jpt)/rdrho_local_p1
     >-(ti_m(j-jptf)/te_m(j-jptf))*csda_m(j-jptf)*
     > (zpti_m(j-jptf)+zpni_m(j-jptf))
     >  *pgeo_local*rhosda_m(j-jptf)/rho(j-jptf)/rdrho_local
     >  )/(rho(j+jpt)-rho(j-jptf))/csda_m(j)
     >  )    
c
      do jj=1,2
       if(jj.eq.1) j=jm-jptf
       if(jj.eq.2.D0) j=jm+jpt
c inside j=jm-jptf outside j=jm+jpt
c banana regime ie collisionless limit formulas
        alpha_neo_hold=alpha_neo
        fc=1-1.46D0*(rmin_exp(j)/rmaj_exp(j))**0.5D0+
     >      0.46D0*(rmin_exp(j)/rmaj_exp(j))**1.5D0
        akappa1=0.8839D0*fc/(0.3477D0+0.4058*fc)
        if(irot1.eq.1) alpha_neo=-akappa1+1.D0
c
        egeo_local=1.D0
        pgeo_local=1.D0
        rdrho_local=1.D0
        rdrho_local_p1=1.D0
        if(igeo_m.ge.5) egeo_local=georotrate(j)
        if(igeo_m.ge.3) then
          pgeo_local=drhodr(j)
          rdrho_local=rmin_exp(j)/arho_exp/rho(j)
        endif
c
        ve(j)=-(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
     > (zpni_m(j)+alpha_neo*zpti_m(j))*vstar_sign*pgeo_local
     > -rho(j)*rdrho_local*
     >  arho_exp/rmajor_exp/q_exp(j)*rmajor_exp*angrotp_exp(j)
        vpar(j)=rmajor_exp*angrotp_exp(j)-vstar_sign*
     > (ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*pgeo_local*
     >((alpha_neo-1.)*zpti_m(j))*rho(j)*rdrho_local*
     > arho_exp/rmajor_exp/q_exp(j)
        vmode(j)=anfreq_m(j)/(ky_j(j)+1.D-10)*
     >    csda_m(j)*arho_exp*rhosda_m(j)
c
       if(abs(itport_pt(4)).eq.0.and.itport_pt(5).eq.0) then
       vphi_m(j)=rmajor_exp*angrotp_exp(j)
       vpar_m(j)=vpar(j)
c       vper_m(j)=ve(j)
c     >+(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
c     > (zpni_m(j)+zpti_m(j))*vstar_sign*pgeo_local 
       endif
c
       if(abs(itport_pt(4)).eq.1.and.itport_pt(5).eq.1) then
         vpar(j)=vpar_m(j)
       endif             
c
       if(abs(itport_pt(4)).eq.1.and.itport_pt(5).eq.0) then
c this option vpar is vphi vexb from neo+vphi
        ve(j)=-(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
     > (zpni_m(j)+alpha_neo*zpti_m(j))*vstar_sign*pgeo_local
     > -rho(j)*rdrho_local*
     >  arho_exp/rmajor_exp/q_exp(j)*corot*vphi_m(j)
        vpar(j)=corot*vphi_m(j)-vstar_sign*
     > (ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*pgeo_local*
     >((alpha_neo-1.D0)*zpti_m(j))*rho(j)*rdrho_local*
     > arho_exp/rmajor_exp/q_exp(j) 
        vpar_m(j)=vpar(j)
c        vper_m(j)=ve(j)
c     >+(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
c     > (zpni_m(j)+zpti_m(j))*vstar_sign*pgeo_local
       endif
c
       if(itport_pt(5).eq.1) then
c this option vexb from vper and vpar with neo dampng built into
c vpar and vper tranport equations
c         ve(j)=vper_m(j)
c     >        -(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
c     >         (zpni_m(j)+zpti_m(j))*vstar_sign*pgeo_local 
       endif
c
        vexb_m(j)=ve(j)
        vmode_m(j)=vmode(j)
        vstar_m(j)=(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
     > (zpni_m(j)+zpti_m(j))*pgeo_local
        vetor_m(j)=-rho(j)*rdrho_local*
     >  arho_exp/rmajor_exp/q_exp(j)*corot*vphi_m(j)
        vepol_m(j)=-(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
     >  akappa1*zpti_m(j)*pgeo_local
c        write(*,*) j, vepol_m(j), -ve(j)*bt_exp, ' Er'
c       write(*,*) j, ve(j), -(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*
c    >    rhosda_m(j)*(zpni_m(j)+alpha_neo*zpti_m(j))*
c    >    vstar_sign*pgeo_local, vetor_m(j)
c
       enddo
c
c get shears outside minus inside
        j=jm
        egeo_local=1.D0
        pgeo_local=1.D0
        rdrho_local=1.D0
        rdrho_local_p1=1.D0
        if(igeo_m.ge.5) egeo_local=georotrate(j)
        if(igeo_m.ge.3) then
          pgeo_local=drhodr(j)
          rdrho_local=rmin_exp(j-jptf)/arho_exp/rho(j-jptf)
          rdrho_local_p1=rmin_exp(j+jpt)/arho_exp/rho(j+jpt)
        endif
      if(igeo_m.le.2.D0) then
        egamma_m(jm)=relx*egamma_m(jm)+(1.D0-relx)*(
     >  (rho(j-jptf)+rho(j+jpt))/2.D0*(ve(j+jpt)/rho(j+jpt)-
     >  ve(j-jptf)/rho(j-jptf))/
     >    (rho(j+jpt)-rho(j-jptf))/arho_exp/csda_m(j)
     >   -(vpar(j+jpt)/(rmajor_exp*q_exp(j+jpt)/shat_exp(j+jpt))
     >   +vpar(j-jptf)/(rmajor_exp*q_exp(j-jptf)/shat_exp(j-jptf)))/
     >   2.D0/csda_m(j)
     >  )
c
        gamma_mode_m(jm)=relx*gamma_mode_m(jm)+(1.D0-relx)*(
     >  (rho(j-jptf)+rho(j+jpt))/2.D0*(vmode(j+jpt)/rho(j+jpt)-
     >  vmode(j-jptf)/rho(j-jptf))/(rho(j+jpt)-rho(j-jptf))/
     >  arho_exp/csda_m(j)
     >  )
c
        gamma_p_m(jm)=relx*gamma_p_m(jm)+(1.D0-relx)*(
     >   -(vpar(j+jpt)-vpar(j-jptf))/(rho(j+jpt)-rho(j-jptf))/arho_exp
     >   /csda_m(j)
     >   )
      endif
      if(igeo_m.ge.3) then
        egamma_m(jm)=relx*egamma_m(jm)+(1.D0-relx)*(
     >   egeo_local*drhodrrrho(j)*
     >  (rho(j-jptf)+rho(j+jpt))/(q_exp(j-jptf)+q_exp(j+jpt))*
     >  (ve(j+jpt)*q_exp(j+jpt)/rho(j+jpt)/rdrho_local_p1-
     >  ve(j-jptf)*q_exp(j-jptf)/rho(j-jptf)/rdrho_local)/
     > (rho(j+jpt)-rho(j-jptf))/arho_exp/csda_m(j)
     >  )
c
        gamma_mode_m(jm)=relx*gamma_mode_m(jm)+(1.D0-relx)*(
     >  egeo_local*drhodrrrho(j)*
     >  (rho(j-jptf)+rho(j+jpt))/2.D0*
     >   (vmode(j+jpt)/rho(j+jpt)/rdrho_local_p1-
     >  vmode(j-jptf)/rho(j-jptf)/rdrho_local)/
     >  (rho(j+jpt)-rho(j-jptf))/arho_exp/csda_m(j)
     >  )
c
        gamma_p_m(jm)=relx*gamma_p_m(jm)+(1.D0-relx)*(
     >   -drhodr(j)*
     >   (vpar(j+jpt)-vpar(j-jptf))/(rho(j+jpt)-rho(j-jptf))/arho_exp
     >   /csda_m(j)
     >  )
c
c.. components of egamma_m
c
        egamma_vphi(jm)=relx*egamma_vphi(jm)+(1.D0-relx)*(
     >   egeo_local*drhodrrrho(j)*
     >  (rho(j-jptf)+rho(j+jpt))/(q_exp(j-jptf)+q_exp(j+jpt))*
     >  (vetor_m(j+jpt)*q_exp(j+jpt)/rho(j+jpt)/rdrho_local_p1-
     >  vetor_m(j-jptf)*q_exp(j-jptf)/rho(j-jptf)/rdrho_local)/
     > (rho(j+jpt)-rho(j-jptf))/arho_exp/csda_m(j)
     >  )
c
        egamma_vpol(jm)=relx*egamma_vpol(jm)+(1.D0-relx)*(
     >   egeo_local*drhodrrrho(j)*
     >  (rho(j-jptf)+rho(j+jpt))/(q_exp(j-jptf)+q_exp(j+jpt))*
     >  (vepol_m(j+jpt)*q_exp(j+jpt)/rho(j+jpt)/rdrho_local_p1-
     >  vepol_m(j-jptf)*q_exp(j-jptf)/rho(j-jptf)/rdrho_local)/
     > (rho(j+jpt)-rho(j-jptf))/arho_exp/csda_m(j)
     >  )
c
        egamma_vstar(jm)=relx*egamma_vstar(jm)+(1.D0-relx)*(
     >   egeo_local*drhodrrrho(j)*
     >  (rho(j-jptf)+rho(j+jpt))/(q_exp(j-jptf)+q_exp(j+jpt))*
     >  (vstar_m(j+jpt)*q_exp(j+jpt)/rho(j+jpt)/rdrho_local_p1-
     >  vstar_m(j-jptf)*q_exp(j-jptf)/rho(j-jptf)/rdrho_local)/
     > (rho(j+jpt)-rho(j-jptf))/arho_exp/csda_m(j)
     >  )
c      write(*,*) jm, egamma_m(jm), egamma_vphi(jm), 
c     >   egamma_vpol(jm), egamma_vstar(jm)
c
      endif
      if(igeo_m.ge.6) then
        egamma_m(jm)=relx*egamma_m(jm)+(1.D0-relx)*(
     >   -georotrate(j)*
     >  (rmin_exp(j-jptf)+rmin_exp(j+jpt))/
     >  (q_exp(j-jptf)+q_exp(j+jpt))*
     >  (angrot_exp(j+jpt)-angrot_exp(j-jptf))/
     >  (rmin_exp(j+jpt)-rmin_exp(j-jptf))/csda_m(j)
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
        zpte_m(jm+1)=zpte_m_hold
        zpti_m(jm+1)=zpti_m_hold
        zpne_m(jm+1)=zpne_m_hold
        zpni_m(jm+1)=zpni_m_hold
c
       return
       end 
