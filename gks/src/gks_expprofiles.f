      subroutine expprofiles(iz,njj,amin)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c     06-nov-01 jek added amin for minor radius - used in taue scalings
c     01-nov-00 jek added nzm, zpmnz, zpnz_exp
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      implicit none
      include 'input.m'
      include 'data_exp.m'
      include 'glf.m'
c
       integer j, jin_p, iz, njj
       real*8 ve(0:jmaxm),vpar(0:jmaxm),vper(0:jmaxm)
       real*8 vnewk3 
       real*8 vstar_sign, egeo_local, pgeo_local
       real*8 rdrho_local, rdrho_local_p1, alpha_neo_hold, fc
       real*8 akappa1, akappa2
       real*8 voltot, dr, drho
       real*8 dvoldr_p, dvoldr_m
       real*8 amin, eps_exp
       real*8 alpha_neo,alpha_dia
       real*8 tem,tim,nem,nim,nzm,nitotm
       real*8 zpmte,zpmti,zpmne,zpmni,zpmnitot,zpmnz
       real*8 qm,rmajm,zeff_e
c
c...setup
c
       alpha_neo=1.D0
       alpha_dia=1.D0
       zeff_e=0.D0
c
       eps_exp = amin / rmajor_exp
c
c       if(lprint_gf.eq.98) open(2)
c       if(lprint_gf.eq.-98) close(2)
c
c experimental effective diffusion in gyrobohm units
c
        do j=1,jmaxm-1
c         tem=dabs(te_exp(j))
c         tim=dabs(ti_exp(j))
c         nem=dabs(ne_exp(j))
c         nim=dabs(ni_exp(j))
c         nzm=dabs(nz_exp(j))
         tem = 0.5*(te_exp(j+1)+te_exp(j))
         tim = 0.5*(ti_exp(j+1)+ti_exp(j))
         nem = 0.5*(ne_exp(j+1)+ne_exp(j))
         nim = 0.5*(ni_exp(j+1)+ni_exp(j))
         nzm = 0.5*(nz_exp(j+1)+nz_exp(j))
         nitotm = 0.5*(nitot_exp(j+1)+nitot_exp(j))   
c
c         zpmte=-(dlog(te_exp(j-1))-dlog(te_exp(j)))/(rho(j-1)-rho(j))
c     >    + 1.D-10
c         zpmti=-(dlog(ti_exp(j-1))-dlog(ti_exp(j)))/(rho(j-1)-rho(j))
c     >    + 1.D-10
c         zpmne=-(dlog(ne_exp(j-1))-dlog(ne_exp(j)))/(rho(j-1)-rho(j))
c     >    + 1.D-10
c         write(*,*) j, rho(j), zpmne, ' zpmne-derivedexpprofiles'
c         zpmni=-(dlog(ni_exp(j-1))-dlog(ni_exp(j)))/(rho(j-1)-rho(j))
c     >    + 1.D-10
c         zpmnitot=-(dlog(nitot_exp(j-1))-dlog(nitot_exp(j)))/
c     >    (rho(j-1)-rho(j))+ 1.D-10
c         zpmnz=-(dlog(nz_exp(j-1))-dlog(nz_exp(j)))/
c     >    (rho(j-1)-rho(j)) + 1.D-10
         drho=rho(j+1)-rho(j)
         zpmte = -((te_exp(j+1)-te_exp(j))/drho)/tem
         zpmti = -((ti_exp(j+1)-ti_exp(j))/drho)/tim
         zpmne = -((ne_exp(j+1)-ne_exp(j))/drho)/nem
         zpmni = -((ni_exp(j+1)-ni_exp(j))/drho)/nim
         zpmnz = -((nz_exp(j+1)-nz_exp(j))/drho)/nzm
         zpmnitot =-((nitot_exp(j+1)-nitot_exp(j))/drho)/nitotm
         if(ABS(zpmte).lt.1.D-10)zpmte=1.D-10
         if(ABS(zpmti).lt.1.D-10)zpmti=1.D-10
         if(ABS(zpmne).lt.1.D-10)zpmne=1.D-10
         if(ABS(zpmni).lt.1.D-10)zpmni=1.D-10
         if(ABS(zpmnz).lt.1.D-10)zpmnz=1.D-10
         if(ABS(zpmnitot).lt.1.D-10)zpmnitot=1.D-10
c
         if(i_dengrad.eq.0) zpmni=zpmne
c
         zpte_exp(j)=zpmte
         zpti_exp(j)=zpmti
         zpne_exp(j)=zpmne
         zpni_exp(j)=zpmni
         zpnitot_exp(j)=zpmnitot
         zpnz_exp(j)=zpmnz
c
        rlne_exp(j)=rmajor_exp/arho_exp*zpmne*drhodr(j)
        rlni_exp(j)=rmajor_exp/arho_exp*zpmni*drhodr(j)
        rlte_exp(j)=rmajor_exp/arho_exp*zpmte*drhodr(j)
        rlti_exp(j)=rmajor_exp/arho_exp*zpmti*drhodr(j)
c
        cgyrobohm_exp(j)=1.D-4*979.D3*dsqrt(tem*1.D3)/
     >  (arho_exp*100.D0)*(1.02D2*dsqrt(tem*1.D3)/bt_exp/1.D4)**2.D0*
     >  dsqrt(amassgas_exp)
c 
        csda_exp(j)=979.D3*dsqrt(tem*1.D3)/(arho_exp*100.D0)/
     >  dsqrt(amassgas_exp)
cjek
        if (bt_flag .gt. 0) then
          rhosda_exp(j)=((102.D0*dsqrt(tem*1.D3))/bteff_exp(j)/1.D4)
     >    *dsqrt(amassgas_exp)/(arho_exp*100.D0)
        else
          rhosda_exp(j)=((102.D0*dsqrt(tem*1.D3))/bt_exp/1.D4)
     >    *dsqrt(amassgas_exp)/(arho_exp*100.D0)
        endif
c
        betae_exp(j) = 403.D0*nem*tem/(1.D5*bt_exp**2D0)
        betai_exp(j) = 403.D0*nim*tim/(1.D5*bt_exp**2D0)
c
crew    gks collisionality (xnu/w_star_i)*(ky*rho_i)
        vnewk3=
     >   0.117D0*nem*tem**(-1.5D0)/dsqrt(tim)*(arho_exp)*
     >   dsqrt(amassgas_exp/2.D0)
crew    as used in gks multiply by 1/2 and takout any 1/2 factor in solfp
crew          vnewk3=vnewk3/2.D0
        xnu_exp(j) =vnewk3/dsqrt(2.D0*tem/tim)
        xnu_exp(j)=xnu_exp(j)*(zeff_exp(j)+zeff_e)
c
       vnewstare_exp(j)=zeff_exp(j)*2.91D-6*nem*1.D13*15.D0/
     >  (tem*1.D3)**2*rmaj_exp(jm)*100.D0*q_exp(j)
     >     /(rmin_exp(j)/rmaj_exp(jm)+1.D-10)**1.5D0/419.D5
       vnewstari_exp(j)=4.78D-8*nem*1.D13*15.D0/
     >  (tim*1.D3)**2.D0*rmaj_exp(jm)*100.D0*q_exp(j)
     >     /(rmin_exp(j)/rmaj_exp(jm)+1.D-10)**1.5D0/979.D3
c
        enddo
        rhosda_exp(0)=rhosda_exp(1)
        csda_exp(0)=csda_exp(1)
        csda_exp(jmaxm)=979.D3*dsqrt(te_exp(jmaxm)*1.D3)
     >  /(arho_exp*100.D0*dsqrt(amassgas_exp))
        zpte_exp(0)=0.0
        zpti_exp(0)=0.0
        zpne_exp(0)=0.0
        zpni_exp(0)=0.0
        zpnitot_exp(0)=0.0
        zpnz_exp(0)=0.0
        xnu_exp(0)=xnu_exp(1)
        betae_exp(0)=betae_exp(1)
        betai_exp(0)=betai_exp(1)
c
c note elongation dependence of alpha
c note zpne used as better represenation of total pressure than zpni
       do j=1,jmaxm-1
        if(igeo_m.ge.0) then
          qm=0.5*(q_exp(j+1)+q_exp(j))
          rmajm=0.5*(rmaj_exp(j+1)+rmaj_exp(j))
          alpha_exp(j)=drhodr(j)*qm**2*rmajm/
     >    arho_exp*(-400.D0)*(ptot_exp(j+1)-ptot_exp(j))/
     >    (rho(j+1)-rho(j))/(1.D5*bt_exp**2)
        elseif(igeo_m.ge.4) then
          alpha_exp(j)=geoalpha(j)*
     >    qm**2*rmajm/arho_exp*
     >    400.D0*(ptot_exp(j+1)-ptot_exp(j))/
     >    (rho(j+1)-rho(j))/(1.D5*bt_exp**2)
        elseif(igeo_m.eq.-1) then 
          alpha_exp(j)=qm**2*rmajor_exp/arho_exp*
     >    dsqrt(elonga_exp)/((1.D0+elonga_exp**2)/2.D0)*
     >    400.D0*(ptot_exp(j+1)-ptot_exp(j))/
     >    (rho(j+1)-rho(j))/(1.D5*bt_exp**2)
        endif
       enddo
       alpha_exp(0)=0.0
       q_exp(0)=q_exp(1)
c
c errors in pgeo_local and rdrho_local fixed 11/15/96
c
c rotational shear calculation is circular with radius=rho
c vstarp_exp is diamagnetic part of egamma_exp (doppler shear rate)
c vstar_sign is negative for negative vtar_i. Thus for co-injection 
c or positive angrot toroidal rotation cancels the diamgnetic rotation
c
       vstar_sign=-1.D0*alpha_dia
       do j=1,jmaxm-1
        egeo_local=1.D0
        pgeo_local=1.D0
        rdrho_local=1.D0
        if(igeo_m.ge.5) egeo_local=georotrate(j)
        if(igeo_m.ge.3) pgeo_local=drhodr(j)  
        if(igeo_m.ge.3) rdrho_local=rmin_exp(j)/arho_exp/rho(j)
        vstarp_exp(j)=vstar_sign*(rho(j+1)+rho(j))/2.D0*(
     >  (ti_exp(j+1)/te_exp(j+1))*csda_exp(j+1)*
     >  (zpti_exp(j+1)+zpni_exp(j+1))
     >  *pgeo_local*rhosda_exp(j+1)/rho(j+1)
     >  -(ti_exp(j)/te_exp(j))*csda_exp(j)*(zpti_exp(j)+zpni_exp(j))
     >  *pgeo_local*rhosda_exp(j)/rho(j)
     >  )/(rho(j+1)-rho(j))/csda_exp(j)
c
        vstar_exp(j)=pgeo_local*
     >  (ti_exp(j)/te_exp(j))*csda_exp(j)*arho_exp*rhosda_exp(j)*
     >  (zpni_exp(j)+zpti_exp(j))*vstar_sign*pgeo_local
       enddo

       alpha_neo_hold=alpha_neo
       do j=1,jmaxm
        egeo_local=1.D0
        pgeo_local=1.D0
        rdrho_local=1.D0
        if(igeo_m.ge.5) egeo_local=georotrate(j)
        if(igeo_m.ge.3) then
         pgeo_local=drhodr(j)
         rdrho_local=rmin_exp(j)/arho_exp/rho(j) 
        endif

c bananna regime ie collisionless limit formulas 
        fc=1.D0-1.46D0*dsqrt(rmin_exp(j)/rmaj_exp(j)) +
     >      0.46D0*(rmin_exp(j)/rmaj_exp(j))**1.5D0
        akappa2=0.D0
c error fixed 5/28/96  if(irot2.eq.1) akappa2=1.-fc/(1.+1.1671*fc)
        if(irot2.eq.1) akappa2=(1.D0-fc)/(1.D0+1.1671*fc)
        akappa1=0.8839D0*fc/(0.3477D0+0.4058D0*fc)
        if(irot1.eq.1) alpha_neo=-akappa1+1.D0
        if(irot1.eq.0) akappa1=1.D0-alpha_neo
c trace impurity limit to go from impurity rotation angrot_exp to
c plasma rotation angrotp_exp
        angrotp_exp(j)=angrot_exp(j)+
     >     akappa2*1.5D0*csda_exp(j)*zpti_exp(j)*(ti_exp(j)/te_exp(j))
     >     /rho(j)*q_exp(j)*rhosda_exp(j)*pgeo_local/rdrho_local
        if(angrot_exp(j).eq.0.) angrotp_exp(j)=0.D0
        angrotp_exp(j)=corot*angrotp_exp(j)
c check on exb rotation
crew 5/11/98
        vexbc_exp(j)=vstar_sign*(rdrho_local*
     > rho(j)*arho_exp/rmajor_exp/q_exp(j)*rmajor_exp*angrotp_exp(j)+
     > csda_exp(j)*rhosda_exp(j)*arho_exp*(ti_exp(j)/te_exp(j))*
     > ((akappa1)*zpti_exp(j)*pgeo_local
     > -(zpti_exp(j)+zpni_exp(j))*pgeo_local))     
c
        ve(j)=-(ti_exp(j)/te_exp(j))*csda_exp(j)*arho_exp*
     >  rhosda_exp(j)*(zpni_exp(j)+zpti_exp(j)+
     >  (1.D0-rdrho_local*rho(j)*arho_exp/rmajor_exp/q_exp(j))*
     >  (alpha_neo-1.D0)*zpti_exp(j))*vstar_sign*pgeo_local
c
        ve(j)=ve(j)
     > -rdrho_local*
     >  rho(j)*arho_exp/rmajor_exp/q_exp(j)*rmajor_exp*angrotp_exp(j)
c
        vexb_exp(j)=ve(j)
c
c note vexb=vstar_sign*(-d phi/dr)*c/bt=vstar_sign*Er*c/bt
c
        vpar(j)=rmajor_exp*angrotp_exp(j)+
     >  vstar_sign*((1.D0-alpha_neo)*zpti_exp(j))*
     >  (ti_exp(j)/te_exp(j))*csda_exp(j)*arho_exp*rhosda_exp(j)*
     >  pgeo_local*rho(j)*arho_exp/rmajor_exp/q_exp(j)*rdrho_local
c      
        vper(j)=-(ti_exp(j)/te_exp(j))*csda_exp(j)*arho_exp*
     >  rhosda_exp(j)*((1.D0-rdrho_local*rho(j)*arho_exp/rmajor_exp/
     >  q_exp(j))*(alpha_neo-1.D0)*zpti_exp(j))*vstar_sign*pgeo_local
c
        vper(j)=vper(j)    
     > -rdrho_local*
     >  rho(j)*arho_exp/rmajor_exp/q_exp(j)*rmajor_exp*angrotp_exp(j)  
c
       vexbcc_exp(j)=vper(j)
     >-(ti_exp(j)/te_exp(j))*csda_exp(j)*arho_exp*rhosda_exp(j)*
     > (zpni_exp(j)+zpti_exp(j))*vstar_sign*pgeo_local 
c      
        vpar_exp(j)=vpar(j)
        vper_exp(j)=vper(j)
        vphi_exp(j)=rmajor_exp*angrotp_exp(j)
        vecur_exp(j) = -curden_exp(j)/(1.6022*ne_exp(j))  !m/sec
        mach_i_exp(j) = alpha_mach*vpar(j)/(csda_exp(j)*arho_exp)
        mach_e_exp(j) = alpha_cur*vecur_exp(j)/(csda_exp(j)*arho_exp)
     >  +mach_i_exp(j)
       enddo
c
       angrotp_exp(0)=angrotp_exp(1)
       vphi_exp(0)=vphi_exp(1)
       vpar_exp(0)=vpar_exp(1)
       vper_exp(0)=vper_exp(1)
       vecur_exp(0)=vecur_exp(1)
       mach_i_exp(0)=mach_i_exp(1)
       mach_e_exp(0)=mach_e_exp(1)
c
       do j=1,jmaxm-1
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
        egamma_exp(j)=
     >  (rho(j)+rho(j+1))/2.D0*
     >  (ve(j+1)/rho(j+1)-
     >  ve(j)/rho(j))
     >    /(rho(j+1)-rho(j))/arho_exp/csda_exp(j)
     >   -(vpar(j+1)/(rmajor_exp*q_exp(j+1)/shat_exp(j+1))
     >   +vpar(j)/(rmajor_exp*q_exp(j)/shat_exp(j)))/2.D0/csda_exp(j)
        gamma_p_i_exp(j)=
     >   -alpha_p*(vpar(j+1)-vpar(j))/(rho(j+1)-rho(j))/arho_exp
     >   /csda_exp(j)
        gamma_p_e_exp(j)=
     >   -alpha_p_cur*(vecur_exp(j+1)-vecur_exp(j))/(rho(j+1)-rho(j))
     >    /(arho_exp*csda_exp(j)) + gamma_p_i_exp(j)
      endif
      if(igeo_m.ge.3) then
        egamma_exp(j)=egeo_local*drhodrrrho(j)*
     >  (rho(j)+rho(j+1))/
     >  (q_exp(j)+q_exp(j+1))*
     >  (ve(j+1)*q_exp(j+1)/rho(j+1)/rdrho_local_p1-
     >  ve(j)*q_exp(j)/rho(j)/rdrho_local)/
     >  (rho(j+1)-rho(j))/arho_exp/csda_exp(j)
        gamma_p_i_exp(j)=
     >   -alpha_p*drhodr(j)*(vpar(j+1)-vpar(j))/(rho(j+1)-rho(j))
     >   /(arho_exp*csda_exp(j))
        gamma_p_e_exp(j)=
     >   -alpha_p_cur*(vecur_exp(j+1)-vecur_exp(j))/(rho(j+1)-rho(j))
     >   /(arho_exp*csda_exp(j)) + gamma_p_i_exp(j)
      endif
      if(igeo_m.ge.6) then
        egamma_exp(j)=-georotrate(j)*
     >  (rmin_exp(j)+rmin_exp(j+1))/
     >  (q_exp(j)+q_exp(j+1))*
     >  (angrot_exp(j+1)-angrot_exp(j))/
     >  (rmin_exp(j+1)-rmin_exp(j))/csda_exp(j)

        vstarp_exp(j)=-corot*
     >  (zpti_exp(j)+zpne_exp(j))*
     >  drhodr(j)*zpti_exp(j)*
     >  (ti_exp(j)/te_exp(j))*rhosda_exp(j)
      endif   
c
c      if(iexb.eq.1) egamma_exp(j)=megamma_exp(j)
c      if(use_xneo_m.ge.2 .or. iexb.eq.2) then
c        egamma_ncl_exp(j)=omexb_ncl_exp(j)/csda_exp(j)
c      endif
c     write(*,50) j, rho(j), egamma_exp(j), egamma_ncl_exp(j)
c      if(iexb.eq.2) egamma_exp(j)=egamma_ncl_exp(j)
c
       enddo
c
       alpha_neo=alpha_neo_hold
        voltot=0.D0
c
        do j=1,jmaxm
         dr=(rho(j)-rho(j-1))*arho_exp
         dvoldr_p=sfactor(j)
         dvoldr_m=sfactor(j-1)
c   
c
         voltot=voltot+0.5D0*(dvoldr_p+dvoldr_m)*dr
c
        enddo
c
c
       tem=0.D0
       tim=0.D0
       nem=0.D0
       nim=0.D0
       zpmte=0.D0
       zpmti=0.D0
       zpmne=0.D0
       zpmni=0.D0
c
      do j=1,jmaxm-1
         zpti2_exp(j)=(ti_exp(j+1)+ti_exp(j-1)-2.D0*ti_exp(j))/
     >                (rho(j)-rho(j-1))**2
c         write(*,55) j, rho(j), ti_exp(j), zpti2_exp(j)
      enddo
c
 50   format(i2,2x,0p1f4.2,0p6e13.5)
 55   format(i3,2x,0p1f6.4,0p6e13.5)
      return
      end
