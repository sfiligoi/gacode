      subroutine derivedexpprofiles(iz,njj,nx,amin,totcur,
     >           rho_d,q_d,pow_ped,w_ped,ismooth_vrot)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c     06-nov-01 jek added amin for minor radius - used in taue scalings
c     01-nov-00 jek added nzm, zpmnz, zpnz_exp
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      implicit none
      include 'mpif.h'
      include '../inc/input.m'
c      include 'data.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/glf.m'
c
       character*120 msg
       integer j, jin_p, iz, njj, nx, iflag
       integer ismooth_vrot
       real*8 ve(0:jmaxm),vpar(0:jmaxm),vper(0:jmaxm)
       real*8 rho_d(nx), q_d(nx)
       real*8 r_tar(1), q_tar(1)
       real*8 vnewk3, thresitg, chigb_simpair, chib_simpair
       real*8 chig_simpair, bpol_h, r_h, grad_ti_h, eps_h
       real*8 chi_h, vstar_sign, egeo_local, pgeo_local
       real*8 rdrho_local, rdrho_local_p1, alpha_neo_hold, fc
       real*8 akappa1, akappa2, wstre, wstri, wstre_ped, wstri_ped
       real*8 wstri_tot
       real*8 wstre_core_exp, wstri_core_exp
       real*8 destr, distr, aveniti, aveniti2, voltot, avene, dr
       real*8 dvoldr_p, dvoldr_m, gradte
       real*8 amin, tocur, totcur, eps_exp, area_exp
       real*8 zpi, q95, q_cyl, q_sh, volume, vol_ped, convert,
     >        n_ped, powe_ped, powi_ped, pow_ped, w_ped(8)
       real*8 zeffline_exp,drho,qm,rmajm
c
c...setup
c
       zpi = atan2(0.D0,-1.D0)
       kevdsecpmw=1.6022D-19*1.D3*1.D-6
       convert = 1.6022D-22
       zpticrit=1.D0
       epsthres=0.5D0
       thresitg0=1.D0
c
       tocur=dabs(totcur)*1.D-6  ! Ip (MA)
       eps_exp = amin / rmajor_exp
       area_exp = zpi*elonga_exp*amin**2.D0
c
       if(lprint_gf.eq.98) open(2)
       if(lprint_gf.eq.-98) close(2)
c
c experimental effective diffusion in gyrobohm units
c
        if(igyro.eq.1) then
         do j=2,jmaxm
           drhodr(j) = (rho(j)-rho(j-1))*arho_exp/
     >                 (rmin_exp(j)-rmin_exp(j-1))
           zpti_exp(j) = -(log(ti_exp(j-1))-log(ti_exp(j)))/
     >                   (rho(j-1)-rho(j))
         enddo
         drhodr(1) =  drhodr(2)
         zpti_exp(1) = zpti_exp(2)
        endif
c
        do j=1,jmaxm
c
         egamma_ncl_exp(j)=0.D0
c
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
c         if(i_dengrad.eq.0) zpmni=zpmne
         drho = rho(j+1)-rho(j)
         zpmte=-((te_exp(j+1)-te_exp(j))/drho)/tem
         zpmti=-((ti_exp(j+1)-ti_exp(j))/drho)/tim
         zpmne=-((ne_exp(j+1)-ne_exp(j))/drho)/nem
         zpmni=-((ni_exp(j+1)-ni_exp(j))/drho)/nim
         zpmnz=-((nz_exp(j+1)-nz_exp(j))/drho)/nzm
c
         if(igyro.ne.1) then
           zpte_exp(j)=zpmte
           zpti_exp(j)=zpmti
           zpne_exp(j)=zpmne
           zpni_exp(j)=zpmni
c          zpnitot_exp(j)=zpmnitot
           zpnz_exp(j)=zpmnz
         endif
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
        if(igyro.eq.1) then
          csda_exp(j)=979.D3*dsqrt(te_exp(j)*1.D3)/(arho_exp*100.D0)/
     >    dsqrt(amassgas_exp)
          rhosda_exp(j)=((102.D0*dsqrt(te_exp(j)*1.D3))/bt_exp/1.D4)
     >    *dsqrt(amassgas_exp)/(arho_exp*100.D0)
        endif
c
        betae_exp(j) = 400.D0*nem*tem/(1.D5*bt_exp**2D0)
        betai_exp(j) = 400.D0*nim*tim/(1.D5*bt_exp**2D0)
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
         chiegb_exp(j)=(powe_exp(j)-1.5D0*xconv*tem*flow_exp(j))
     >            /(kevdsecpmw*tem*nem*1.D19/arho_exp*
     >     gradrhosq_exp(j)*sfactor(j)*zpmte)/cgyrobohm_exp(j)
         chiigb_exp(j)=(powi_exp(j)-1.5D0*xconv*tim*flow_exp(j))
     >           /(kevdsecpmw*tim*nim*1.D19/arho_exp*
     >     gradrhosq_exp(j)*sfactor(j)*zpmti)/cgyrobohm_exp(j)
         diffgb_exp(j)=flow_exp(j)/(kevdsecpmw*nem*1.D19/arho_exp*
     >     gradrhosq_exp(j)*sfactor(j)*zpmne)/cgyrobohm_exp(j)

         chieneogb_exp(j)=chieneo_exp(j)/cgyrobohm_exp(j)
         chiineogb_exp(j)=chiineo_exp(j)/cgyrobohm_exp(j)
         powineo_exp(j)=
     >    kevdsecpmw*ti_exp(j)*ni_exp(j)*1.D19/arho_exp*
     >    gradrhosq_exp(j)*sfactor(j)*(chiineo_exp(j)*zpti_exp(j))
c
c construct mathematical gyrobohm and bohm fitting models
c 
      if (epsthres.ge.0.) then
c  
        thresitg= (1.D0-zpticrit/zpmti/dsqrt(elong_exp(j)))
      if(igeo_m.ge.1)
     >   thresitg= (1.D0-zpticrit/zpmti/drhodr(j))
        if (thresitg.gt.0.)  thresitg=thresitg**epsthres
        if (thresitg.le.0.)  thresitg=thresitg0
        if (epsthres.eq.0.)  thresitg=1.D0
c
         chigb_simpair=q_exp(j)**2.D0*
     >  (arho_exp/rmaj_exp(j))**(1.D0-epsthres)
     > *dabs(zpmti*dsqrt(elong_exp(j)))**epsthres*thresitg
        if(igeo_m.ge.1) 
     >   chigb_simpair=q_exp(j)**2.D0*
     >  (arho_exp/rmaj_exp(j))**(1.D0-epsthres)
     > *dabs(zpmti*drhodr(j))**epsthres*thresitg
        chib_simpair=chigb_simpair/
     >  (q_exp(j)*arho_exp/rmaj_exp(j)*rhosda_exp(j))
        chig_simpair=chigb_simpair/
     >  (q_exp(j)*arho_exp/rmaj_exp(j)*rhosda_exp(j))**1.5D0
c
      endif
c
      if(epsthres.lt.0.) then
c
        chigb_simpair=q_exp(j)**2.D0*(arho_exp/rmaj_exp(j))*
     >   (rmaj_exp(j)/arho_exp*dsqrt(elong_exp(j))*zpmti)**2.D0
        if(igeo_m.ge.1) 
     >  chigb_simpair=q_exp(j)**2.D0*(arho_exp/rmaj_exp(j))*
     >   (rmaj_exp(j)/arho_exp*drhodr(j)*zpmti)**2.D0
        chib_simpair=chigb_simpair/
     >  (q_exp(j)*arho_exp/rmaj_exp(j)*rhosda_exp(j))
        chig_simpair=chigb_simpair/
     >  (q_exp(j)*arho_exp/rmaj_exp(j)*rhosda_exp(j))**1.5D0
c  
c hseih model as example of goldston scaling
c
       bpol_h=bt_exp*arho_exp*rho(j)/rmaj_exp(j)/q_exp(j)+1.D-10
       r_h=rmin_exp(j)
       grad_ti_h=zpmti/arho_exp*dsqrt(elong_exp(j))
       eps_h=rmin_exp(j)/rmaj_exp(j)
       chi_h=9.D-2*(nem*tem/bpol_h**2.D0)*(eps_h)*
     > (r_h*grad_ti_h)**2.D0
c
      if(epsthres.eq.-10.) chig_simpair=chi_h/cgyrobohm_exp(j)
c
      endif            
 
         fgbe_m(j)=dabs(chiegb_exp(j)/chigb_simpair)
         fgbi_m(j)=dabs(chiigb_exp(j)/chigb_simpair) 
         fbe_m(j)=dabs(chiegb_exp(j)/chib_simpair)
         fbi_m(j)=dabs(chiigb_exp(j)/chib_simpair)
         fge_m(j)=dabs(chiegb_exp(j)/chig_simpair)
         fgi_m(j)=dabs(chiigb_exp(j)/chig_simpair)
c
c  chie_eff and chii_eff in m**2/sec as normally defined by conduction
c
         chiecond_exp(j)=(powe_exp(j)-2.5D0*te_exp(j)*flow_exp(j))
     >    /(kevdsecpmw*tem*nem*1.D19/arho_exp*
     >     gradrhosq_exp(j)*sfactor(j)*zpmte)
c     
         chiicond_exp(j)=(powi_exp(j)-2.5D0*ti_exp(j)*flow_exp(j))
     >    /(kevdsecpmw*tim*nim*1.D19/arho_exp*
     >     gradrhosq_exp(j)*sfactor(j)*zpmti)  
c
        enddo
c note elongation dependence of alpha
c note zpne used as better represenation of total pressure than zpni
       do j=1,jmaxm-1
        qm=0.5*(q_exp(j+1)+q_exp(j))
        rmajm=0.5*(rmaj_exp(j+1)+rmaj_exp(j))
        if(igeo_m.ge.0) then
          alpha_exp(j)=drhodr(j)*qm**2*rmajm/
     >    arho_exp*(-400.D0)*(ptot_exp(j+1)-ptot_exp(j))/
     >    (rho(j+1)-rho(j))/(1.D5*bt_exp**2)
        elseif(igeo_m.ge.4) then
          alpha_exp(j)=geoalpha(j)*qm**2*rmajm/arho_exp*
     >    400.D0*(ptot_exp(j+1)-ptot_exp(j))/
     >    (rho(j+1)-rho(j))/(1.D5*bt_exp**2)
        elseif(igeo_m.eq.-1) then 
          alpha_exp(j)=qm**2*rmajor_exp/arho_exp*
     >    dsqrt(elonga_exp)/((1.D0+elonga_exp**2.D0)/2.D0)*
     >    400.D0*(ptot_exp(j+1)-ptot_exp(j))/
     >    (rho(j+1)-rho(j))/(1.D5*bt_exp**2)
        endif
       enddo
       alpha_exp(0)=0.0
       alpha_exp(jmaxm)=alpha_exp(jmaxm-1)
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
c
c        write(*,111) j, rho(j), rmin_exp(j), angrotp_exp(j)
c
c check on exb rotation
crew 5/11/98
c        vexbc_exp(j)=vstar_sign*(rdrho_local*
c     >rho(j)*arho_exp/rmajor_exp/q_exp(j)*rmajor_exp*angrot_exp(j)+
c     > csda_exp(j)*rhosda_exp(j)*arho_exp*(ti_exp(j)/te_exp(j))*
c     > ((akappa1+1.5D0*akappa2*corot)*zpti_exp(j)*pgeo_local
c     > -(zpti_exp(j)+zpni_exp(j))*pgeo_local))
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
c        vper_exp(j)=vper(j)
        if(use_xneo_m .ge. 2 .and.irot2.gt.0) then
          vphi_exp(j)=vphi_ncl_exp(j)
        else
          vphi_exp(j)=rmajor_exp*angrotp_exp(j)
        endif
c        write(*,111) j, rho(j), angrot_exp(j)
       enddo
c
       if(ismooth_vrot.gt.0) then
         if(i_proc.eq.0.and.iz.eq.1) write(*,*) 'Smoothing vphi_exp'
         call average7_1d(vphi_exp,jmaxm)
       endif
c
       angrotp_exp(0)=angrotp_exp(1)
       vphi_exp(0)=vphi_exp(1)
       vpar_exp(0)=vpar_exp(1)
c       vper_exp(0)=vper_exp(1)
c
c       do j=0,jmaxm
c         write(*,111) j, rho(j), te_exp(j),ti_exp(j)
c       enddo
 111   format(2x,i2,2x,0p1f5.2,1p6e12.4)
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
        gamma_p_exp(j)=
     >   -(vpar(j+1)-vpar(j))/(rho(j+1)-rho(j))/arho_exp
     >   /csda_exp(j)
      endif
      if(igeo_m.ge.3) then
        egamma_exp(j)=egeo_local*drhodrrrho(j)*
     >  (rho(j)+rho(j+1))/
     >  (q_exp(j)+q_exp(j+1))*
     >  (ve(j+1)*q_exp(j+1)/rho(j+1)/rdrho_local_p1-
     >  ve(j)*q_exp(j)/rho(j)/rdrho_local)/
     >  (rho(j+1)-rho(j))/arho_exp/csda_exp(j)
        gamma_p_exp(j)=
     >   -drhodr(j)*(vpar(j+1)-vpar(j))/(rho(j+1)-rho(j))/arho_exp
     >   /csda_exp(j)
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
      if(iexb.eq.1) egamma_exp(j)=megamma_exp(j)
      if(use_xneo_m.ge.2 .or. iexb.eq.3) then
        egamma_ncl_exp(j)=omexb_ncl_exp(j)/csda_exp(j)
      endif
c      write(*,50) j, rho(j), egamma_exp(j), abs(egamma_ncl_exp(j))
      if(iexb.eq.3) egamma_exp(j)=abs(egamma_ncl_exp(j))
c
       enddo
c
       alpha_neo=alpha_neo_hold
c
c experimental confinement times
c
        taue_exp(0)=0.D0
        taui_exp(0)=0.D0
        taupe_exp(0)=0.D0     !taux(0)=infinity...set scale here
        taupi_exp(0)=0.D0
c
        neline_exp=0.D0
        niline_exp=0.D0
        zeffline_exp=0.D0
c
        wstre=0.D0
        wstri=0.D0
        wstri_tot=0.D0
        wstr_exp=0.D0
        wstr_inc_exp=0.D0
        wstr_ped_exp=0.D0
        wstr_core_exp=0.D0
        wstre_core_exp=0.D0
        wstri_core_exp=0.D0
        wstre_ped=0.D0
        wstri_ped=0.D0
        destr=0.D0
        distr=0.D0
        aveniti=0.D0
        aveniti2=0.D0
        avene=0.D0
        voltot=0.D0
c
        pheat_exp=powe_beam_exp(jmaxm)+powi_beam_exp(jmaxm)+
     >            powe_rf_exp(jmaxm)+powi_rf_exp(jmaxm)+
     >            powe_oh_exp(jmaxm)
c
        do j=1,jmaxm
         dr=(rho(j)-rho(j-1))*arho_exp
         dvoldr_p=sfactor(j)
         dvoldr_m=sfactor(j-1)
c   
         wstre=wstre+kevdsecpmw*1.D19*
     >       0.5D0*(dvoldr_p*1.5D0*ne_exp(j)*te_exp(j)+
     >       dvoldr_m*1.5D0*ne_exp(j-1)*te_exp(j-1))*dr
c
         taue_exp(j)=wstre/(powe_exp(j)+1.D-10)
c    
         wstri=wstri+kevdsecpmw*1.D19*
     >       0.5D0*(dvoldr_p*1.5D0*nitot_exp(j)*ti_exp(j)+
     >       dvoldr_m*1.5D0*nitot_exp(j-1)*ti_exp(j-1))*dr
         wstri_tot=wstri_tot+kevdsecpmw*1.D19*
     >       0.5D0*(dvoldr_p*1.5D0*(nitot_exp(j)+
     >       nfast_exp(j))*ti_exp(j)+dvoldr_m*1.5D0*
     >       (nitot_exp(j-1)+nfast_exp(j-1))*ti_exp(j-1))*dr
c
c...  Core and pedestal We and Wi
c
         if (j.le.jout_m) then
           wstre_core_exp = wstre_core_exp + kevdsecpmw*1.D19*
     >       0.5D0*(dvoldr_p*1.5D0*ne_exp(j)*te_exp(j)+
     >       dvoldr_m*1.5D0*ne_exp(j-1)*te_exp(j-1))*dr
           wstri_core_exp = wstri_core_exp + kevdsecpmw*1.D19*
     >       0.5D0*(dvoldr_p*1.5D0*ni_exp(j)*ti_exp(j)+
     >       dvoldr_m*1.5D0*ni_exp(j-1)*ti_exp(j-1))*dr
           wstre_ped = wstre_ped + kevdsecpmw*1.D19*
     >       0.5D0*(dvoldr_p*1.5D0*ne_exp(j)*te_exp(jout_m)+
     >       dvoldr_m*1.5D0*ne_exp(j-1)*te_exp(jout_m))*dr
           wstri_ped = wstri_ped + kevdsecpmw*1.D19*
     >       0.5D0*(dvoldr_p*1.5D0*ni_exp(j)*ti_exp(jout_m)+
     >       dvoldr_m*1.5D0*ni_exp(j-1)*ti_exp(jout_m))*dr
         else
           wstre_ped = wstre_ped + kevdsecpmw*1.D19*
     >       0.5D0*(dvoldr_p*1.5D0*ne_exp(j)*te_exp(j)+
     >       dvoldr_m*1.5D0*ne_exp(j-1)*te_exp(j-1))*dr
           wstri_ped = wstri_ped + kevdsecpmw*1.D19*
     >       0.5D0*(dvoldr_p*1.5D0*ni_exp(j)*ti_exp(j)+
     >       dvoldr_m*1.5D0*ni_exp(j-1)*ti_exp(j-1))*dr
         endif
c
c... taue (don't correct for radiation in total)
c
         taui_exp(j)=wstri/(powi_exp(j)+1.D-10)
         taut_exp(j)=(wstre+wstri)/(powe_exp(j)-
     >               xrad_exp*powe_rad_exp(j)+powi_exp(j))
         tautot_exp(j)=(wstre+wstri_tot)/(powe_exp(j)-
     >               xrad_exp*powe_rad_exp(j)+powi_exp(j))
c
         destr=destr+kevdsecpmw*1.D19*
     >       0.5D0*(dvoldr_p*ne_exp(j)+
     >       dvoldr_m*ne_exp(j-1))*dr
c   
         taupe_exp(j)=destr/(flow_exp(j)+1.D-10)
c
         distr=distr+kevdsecpmw*1.D19*
     >       0.5D0*(dvoldr_p*ni_exp(j)+
     >       dvoldr_m*ni_exp(j-1))*dr

         taupi_exp(j)=distr/(flow_exp(j)+1.D-10)
c
         voltot=voltot+0.5D0*(dvoldr_p+dvoldr_m)*dr
c
         aveniti=aveniti+
     >        0.5D0*(dvoldr_p*ni_exp(j)*ti_exp(j)+
     >        dvoldr_m*ni_exp(j-1)*ti_exp(j-1))*dr
         aveniti2=aveniti2+
     >        0.5D0*(dvoldr_p*ni_exp(j)**2.D0*ti_exp(j)**2.D0+
     >        dvoldr_m*ni_exp(j-1)**2.D0*ti_exp(j-1)**2.D0)*dr
c
         avene=avene+0.5D0*(dvoldr_p*ne_exp(j)+
     >         dvoldr_m*ne_exp(j-1))*dr
c
c... line averaged electron, ion density, zeff
c
         neline_exp=neline_exp+ne_exp(j)
         niline_exp=niline_exp+ni_exp(j)
         zeffline_exp=zeffline_exp+zeff_exp(j)
        enddo
c
        wstre_exp=wstre
        wstri_exp=wstri
        wstre_ped_exp=wstre_ped
        wstri_ped_exp=wstri_ped
        wstr_exp=wstre_exp+wstri_exp
        wstr_ped_exp=wstre_ped_exp+wstri_ped_exp
        wstr_core_exp=wstre_core_exp+wstri_core_exp
        wstr_inc_exp=wstr_exp-wstr_ped_exp
c        write(*,*) 'wtot = ',wstre+wstri_tot
c
        neline_exp=neline_exp/jmaxm
        niline_exp=niline_exp/jmaxm
        zeffline_exp=zeffline_exp/jmaxm
        volaveniti_exp=aveniti2/aveniti
c       gtaut_exp=taut_exp(jmaxm)
        gtaut_exp=wstr_exp/pheat_exp
        gtaui_exp=taui_exp(jmaxm)
c        gtautot_exp=tautot_exp(jmaxm)
        n_t_tau_exp=1.D19*volaveniti_exp*gtaut_exp
        n0_t0_tau_exp=
     >       1.D19*ni_exp(3*jmaxm/10)*ti_exp(3*jmaxm/10)*gtaut_exp
        jin_p=jmaxm/5
        te20_exp=te_exp(jin_p)
        ti20_exp=ti_exp(jin_p)
        rhosda20_exp=rhosda_exp(jin_p)
c     
        volaveniti2_exp=aveniti2/voltot
        volavene_exp=avene/voltot
c
        if (iz.eq.1 .and. i_proc.eq.0) then
        write(6,'(a26,2F8.4)') ' cross. area     (m^2)  = ',area_exp
        write(6,'(a26,2F8.4)') ' amin            (m)    = ',amin
        write(6,'(a26,2F8.4)') ' nebar    (10^19 m^-3)  = ',neline_exp
        write(6,'(a26,2F8.4)') ' volavene (10^19 m^-3)  = ',volavene_exp
        write(6,'(a26,2F8.4)') ' zeff (line-avg)        = ',zeffline_exp
        endif
c        write(*,*) 'voltot = ',voltot
c        write(*,*) 'nebar = ',neline_exp
c        write(*,*) 'volavene = ',volavene_exp
c
c...if no 0d data available for total input power
c...then do the best we can with the 1d data
c
        if (p_glob_exp.le.1.e-3) then
c         p_glob_exp =
c    >      powe_beam_exp(jmaxm)+powe_rf_exp(jmaxm)+powe_oh_exp(jmaxm)
c    >     -(1.D0-xwdot)*powe_wdot_exp(jmaxm)
c    >     +powi_beam_exp(jmaxm)+powi_rf_exp(jmaxm)
c    >     -(1.D0-xwdot)*powi_wdot_exp(jmaxm)
          p_glob_exp =
     >      powe_beam_exp(njj-1)+powe_rf_exp(njj-1)
     >     +powe_oh_exp(njj-1)-xwdot*powe_wdot_exp(njj-1)
     >     +powi_beam_exp(njj-1)+powi_rf_exp(njj-1)
     >     -xwdot*powi_wdot_exp(njj-1)
c          write(6,*) 'p_glob_exp calculated from 1D data',
        endif
c
c... global energy confinements scalings
c    ITER-89P, 97L, 98(y,2), ELMFree DB3
c
       if (iz.eq.1) then
        w_glob_exp = 0.D0
        tauescale_exp(1) = 0.D0
        tauescale_exp(2) = 0.D0
        tauescale_exp(3) = 0.D0
        tauescale_exp(4) = 0.D0
        if (tauscale(1).eq.1) then
          tauescale_exp(1) = 0.048*(tocur**0.85)*             ![MA] ITER-89P
     >                (dabs(rmajor_exp)**1.2D0)*              ![m]
     >                (dabs(amin)**0.3D0)*                    ![m]
     >                (dabs(elonga_exp)**0.5D0)*              !
     >                (dabs(neline_exp/10.D0)**0.1D0)*        ![10^20 m-3]
     >                (dabs(bt_exp)**0.2D0)*                  ![T]
     >                (dabs(amassgas_exp)**0.5D0)*            ![amu]
     >                (dabs(p_glob_exp)**(-0.5D0))            ![MW]
          w_glob_exp = tauescale_exp(1) * p_glob_exp
        endif
c
        if (tauscale(2).eq.1) then
          tauescale_exp(2) = 0.023*(tocur**0.96)*             ![MA] ITER-97L
     >                (dabs(rmajor_exp)**1.89D0)*             ![m]
     >                (dabs(amin)**(-0.06D0))*                ![m]
     >                (dabs(elonga_exp)**0.64D0)*             !
     >                (dabs(neline_exp)**0.4D0)*              ![10^19 m-3]
     >                (dabs(bt_exp)**0.03D0)*                 ![T]
     >                (dabs(amassgas_exp)**0.2D0)*            ![amu]
     >                (dabs(p_glob_exp)**(-0.73D0))           ![MW]
          w_glob_exp = tauescale_exp(2) * p_glob_exp
        endif
        if (tauscale(3).eq.1) then
          tauescale_exp(3) = 0.0562*(tocur**0.93)*            ![MA]  ITER-98(y,2)
     >                (dabs(rmajor_exp)**1.97)*               ![m]
     >                (dabs(eps_exp)**(0.58))*                !
     >                (dabs(elonga_exp)**0.78)*               !
     >                (dabs(neline_exp)**0.41)*               ![10^19 m-3]
     >                (dabs(bt_exp)**0.15)*                   ![T]
     >                (dabs(amassgas_exp)**0.19)*             ![amu]
     >                (dabs(p_glob_exp)**(-0.69))             ![MW]
          w_glob_exp = tauescale_exp(3) * p_glob_exp
        endif
c
        if (tauscale(4).eq.1) then
          tauescale_exp(4) = 0.0198*(tocur**0.85)*            ![MA]  IAEA04 DB3v13, no beta
     >                (dabs(rmajor_exp)**1.21D0)*             ![m]
     >                (dabs(amin)**(-1.25D0))*                !
     >                (dabs(area_exp**0.82))*                 ![m^2]
     >                (dabs(neline_exp)**0.26D0)*             ![10^19 m-3]
     >                (dabs(bt_exp)**0.17D0)*                 ![T]
     >                (dabs(amassgas_exp)**0.11D0)*           ![amu]
     >                (dabs(p_glob_exp)**(-0.45D0))           ![MW]
          w_glob_exp = tauescale_exp(4) * p_glob_exp
        endif
       endif
c
c      if (i_proc.eq.0) then
c         write(*,*) 'R = ',rmajor_exp
c         write(*,*) 'a = ',amin
c         write(*,*) 'kappa = ',elonga_exp
c         write(*,*) 'Bt = ',bt_exp
c         write(*,*) 'Ip = ',tocur
c         write(*,*) 'M = ',amassgas_exp
c         write(*,*) 'nebar = ',neline_exp
c         write(*,*) 'P = ',p_glob_exp
c         write(*,*) 'Pbe = ',powe_beam_exp(njj-1)
c         write(*,*) 'Pbi = ',powi_beam_exp(njj-1)
c         write(*,*) 'Poh = ',powe_oh_exp(njj-1)
c         write(*,*) 'gtaut = ',gtaut_exp
c         write(*,*) 'taueth = ',(wstre+wstri)/(
c    >      powe_beam_exp(njj-1)+powe_rf_exp(njj-1)
c    >     +powe_oh_exp(njj-1)-(1.D0-xwdot)*powe_wdot_exp(njj-1)
c    >     +powi_beam_exp(njj-1)+powi_rf_exp(njj-1)
c    >     -(1.D0-xwdot)*powi_wdot_exp(njj-1))
c         write(*,*) 'Wth = ',wstre+wstri
c      endif
c
       tem=0.D0
       tim=0.D0
       nem=0.D0
       nim=0.D0
       zpmte=0.D0
       zpmti=0.D0
       zpmne=0.D0
       zpmni=0.D0
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... H-mode pedestal scalings 
c
c get q95 first
c
       r_tar(1)=0.95D0
       call w_lin_interp_r8(njj,rho_d,q_d,1,r_tar,q_tar,iflag,msg)
c      do j=1,njj
c         write(*,*) j, rho_d(j), q_d(j)
c      enddo
       q95=q_tar(1)
c
c Power indep MHD scaling - Eqn 12 from 2-term Thomsen paper
c
       q_cyl = 5.D0*elonga_exp*amin**2.D0*bt_exp/rmajor_exp/tocur
       q_sh = q95/q_cyl
       volume = 2.D0*elonga_exp*
     >          zpi**2.D0*rmajor_exp*amin**2.D0
       volume = vol_exp(njj-1) ! above formula not accurate
       vol_ped = rho(jout_m)*volume
       n_ped = ne_exp(jout_m)
       w_ped(1) = dexp(-4.61D0)*tocur**2.D0*rmajor_exp*
     >            (amassgas_exp/n_ped/rmajor_exp**2.D0)**0.13D0
       w_ped(1) = w_ped(1)*q_sh**1.28D0/(amin/rmajor_exp)**1.68D0
       t_ped(1) = w_ped(1)/(convert*3.D0*n_ped*1.D19*vol_ped)
c
c Power dependent scaling - Eqn 1 from 2-term Thomsen paper
c Auxiliary heating + Ohmic only
c
       powe_ped = powe_beam_exp(njj-1)+powe_rf_exp(njj-1)+
     >            powe_lh_exp(njj-1)+powe_oh_exp(njj-1)
       powi_ped = powi_beam_exp(njj-1)+powi_rf_exp(njj-1)
       pow_ped = powe_ped + powi_ped
       w_ped(2) = dexp(-3.74D0)*tocur**1.71D0*rmajor_exp**1.16D0
       w_ped(2) = w_ped(2)*pow_ped**0.31D0*amassgas_exp**0.30D0*
     >            q_sh**1.20D0
       t_ped(2) = w_ped(2)/(convert*3.D0*n_ped*1.D19*vol_ped)
c
c Power dependent scaling - Wped1 from 2-term Cordey NF paper
c Eqn 2 w/ type I, III ELMy data, RMS=23.5%
c Auxiliary heating + Ohmic only, use line-averaged ne
c
       w_ped(3) = 6.43D-4*tocur**1.58D0*rmajor_exp**1.08D0*
     >            pow_ped**0.42D0*neline_exp**(-0.08D0)*
     >            bt_exp**0.06D0*elonga_exp**1.81D0*
     >            (amin/rmajor_exp)**(-2.13D0)*
     >            amassgas_exp**0.2D0*q_sh**2.09D0
       t_ped(3) = w_ped(3)/(convert*3.D0*n_ped*1.D19*vol_ped)
c
c Power dependent scaling - Wped2 from 2-term Cordey NF paper
c Type I ELMy only (no MAST data), no aspect ratio
c Auxiliary heating + Ohmic only, use line-averaged ne
c
       w_ped(4) = 8.07D-3*tocur**1.41D0*rmajor_exp**1.37D0*
     >            pow_ped**0.50D0*neline_exp**(-0.15D0)*
     >            bt_exp**0.32D0*amassgas_exp**0.2D0*q_sh**1.61D0
       t_ped(4) = w_ped(4)/(convert*3.D0*n_ped*1.D19*vol_ped)
c
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c      write(*,*) 'R = ',rmajor_exp
c      write(*,*) 'a =',amin
c      write(*,*) 'kappa = ',elonga_exp
c      write(*,*) 'Bt = ',bt_exp
c      write(*,*) 'Ip = ',tocur
c      write(*,*) 'amass = ',amassgas_exp
c      write(*,*) 'q95 = ',q95
c      write(*,*) 'qcyl = ',q_cyl
c      if(i_proc.eq.0) write(*,*) 'q95 = ',q95
c      if(i_proc.eq.0) write(*,*) 'qsh = ',q_sh
c       write(*,*) 'vol = ',volo_d
c      write(*,*) 'vol_exp = ',vol_exp(njj-1)
c      write(*,*) 'vol-compute = ',volume
c      write(*,*) 'ped = ',rho(jout_m)
c      write(*,*) 'vol_ped  = ',vol_ped
c      write(*,*) 'nped = ',n_ped
c      write(*,*) 'Wped_eq1 = ',w_ped(2)
c      write(*,*) 'Wped_eq12 = ',w_ped(1)
c      write(*,*) 'Wped-eq2(IAEA02) = ',w_ped(3)
c      write(*,*) 'Wped-eq3(IAEA02) = ',w_ped(4)
c      write(*,*) 'pow_ped = ',pow_ped
c      write(*,*) 'Tiped = ',ti_exp(jout_m)
c      write(*,*) 'Teped = ',te_exp(jout_m)
c      write(*,*) 'Tped_eq1 = ',t_ped(2)
c      write(*,*) 'Tped_eq12 = ',t_ped(1)
c      write(*,*) 'Tped-eq2(IAEA02) = ',t_ped(3)
c      write(*,*) 'Tped-eq3(IAEA02) = ',t_ped(4)
c
      if (iz.eq.1 .and. i_proc.eq.0) then
       write(6,'(a26,2F8.4)') ' q95                    = ',q95
      endif
      do j=1,jmaxm
         zpti2_exp(j)=(ti_exp(j+1)+ti_exp(j-1)-2.D0*ti_exp(j))/
     >                (rho(j)-rho(j-1))**2.D0
         zpte2_exp(j)=(te_exp(j+1)+te_exp(j-1)-2.D0*te_exp(j))/
     >                (rho(j)-rho(j-1))**2.D0
         zpne2_exp(j)=(ne_exp(j+1)+ne_exp(j-1)-2.D0*ne_exp(j))/
     >                (rho(j)-rho(j-1))**2.D0
c         write(*,55) j, rho(j), ti_exp(j), zpti2_exp(j)
      enddo
c
 50   format(i2,2x,0p1f4.2,0p6e13.5)
 55   format(i3,2x,0p1f6.4,0p6e13.5)
      return
      end












