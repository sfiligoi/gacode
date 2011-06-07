       subroutine readtestprofiles       
************************************************************************
c      Creates test profiles
c      converts to proper units for code       
************************************************************************
      implicit none
c
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/glf.m'
c
      integer j, jj, iread
      integer pgasa, pgasz, bgasa, bgasz,
     &        nnexp_d, ntexp_d, ntixp_d
      real*8 sum, pohsum, pradsum, peisum
      real*8 prfesum, prfisum, dr, ds, dvoldr_m, dvoldr_p
      real*8 zhalf, zthird, z2thrd, z4thrd, z3qtr
      real*8 xxq, xnq, zeff_t, alpha, cur0, cursum, zbrac, alamda,
     &       tau_e, rfexp, rfixp, rhostar, shift_shaf, aj
c
c     tok = 'test'
c     shot= '999'
c     phase='X'
c     xp_time=99.D0
      nstk=1
c
c     Assume the working and beam ions are same and the impurity
c     is carbon for the moment.
c
      pgasa=amassgas_exp
      pgasz=1
      bgasa=pgasa
      bgasz=1
c     pimpa=12
c     pimpz=6
c
      apgasa=pgasa
      apgasz=pgasz
      abgasa=bgasa
      abgasz=bgasz
      apimpa=pimpa
      apimpz=pimpz
c
      pi_m=3.1415926
      zhalf = 1.D0/2.D0
      zthird = 1.D0/3.D0
      z2thrd = 2.D0/3.D0
      z4thrd = 4.D0/3.D0
      z3qtr = 3.D0/4.D0
c
c... setup grids and constants
c
      do j=0,jmaxm
        aj=j
        rho(j)=aj/jmaxm
      enddo
      rho(0)=1.D-6
      rho(jmaxm)=rho(jmaxm)-1.D-6
c
c nj : the size of the vectors printed in this file
c      set by nj in data.m
c      nj_d=nj
       nj_d=jmaxm+1
c
c time :  time at which data is printed
       time_d=99.D0
c
c R0 : major radius of vaccuum btor ref
       rmajor_d=rmajor_exp     !meters
c
c kappa : plasma elongation
       kappa_d=elonga_exp
c
c Btor : vaccuum toroidal field at rmajor, tesla
       btor_d=bt_exp            !tesla
c
c---rho grid, meters
c---hcap, (ie (dvolume/drho)/(4*pi*pi*R0*rho))
       do j=1,nj_d
         r_d(j)=arho_exp*rho(j-1)            !meters
         hcap_d(j)=1.D0
       enddo
c
c---q (ie safety factor) profile
       do j=1,nj_d
         xxq=rho(j-1)**2
         q_d(j)=1.D0/((1.D0/q0_exp-1.D0/qa_exp)/xxq/alfj_exp*
     >     ((1.D0-xxq)-(1.D0-xxq)**(alfj_exp+1))+1.D0/qa_exp)
       enddo
c
c---zeff: assume flat profile
       do j=1,nj_d
         zeff_d(j)=zeff_sc
       enddo
       zeff_t=zeff_sc
c
c---electron density,#/m**3 and temperatures,kev
c
      do j=0,jmaxm
       xxq=rho(j)**2
       xnq=rho(j)**xpnq_exp
       te_exp(j)=te0_exp*(1.D0-ftea_exp)*(1.D0-xxq)**alfte_exp
     >    +te0_exp*ftea_exp
       ti_exp(j)=ti0_exp*(1.D0-ftia_exp)*(1.D0-xxq)**alfti_exp
     >    +ti0_exp*ftia_exp
       ne_exp(j)=ne0_exp*(1.D0-fnea_exp)*(1.D0-xnq)**alfne_exp
     >    +ne0_exp*fnea_exp
      enddo
c
      do j=1,nj_d
        ene_d(j)=ne_exp(j-1)*1.D19    !#/meter**3
        te_d(j)=te_exp(j-1)           !kev
        ti_d(j)=ti_exp(j-1)           !kev
      enddo
c
c---primary ion density,#/m**3,species
c---impurity ion density,#/m**3,species
      nprim_d=1
      do j=1,nj_d
       en_d(j,1)=(apimpz-zeff_t)/(apimpz-apgasz)/apgasz
     >          *ene_d(j)
       en_d(j,2)=(zeff_t-apgasz)/(apimpz-apgasz)/apimpz
     >          *ene_d(j)
      enddo
c
c...   JAK 9/06/95 Reset the <real> experimental data arrays
c...   as used fior chisquare between te_exp and te_m
c
       ntexp_d=0
       ntixp_d=0
       nnexp_d=0
c
c---sion   : source due to ionization,        #/(m**3*sec),species
c---srecom : source due to recombination,     #/(m**3*sec),species
c---scx    : source due to cx thermal neut.,  #/(m**3*sec),species
c---sbcx   : sink due to cx with beam neut.   #/(m**3*sec),species
c---s      : total source rate,               #/(m**3*sec),species
c---dudt   : s dot,                           #/(m**3*sec),species
      do j=1,nj_d
       do jj=1,nprim_d
         sion_d(j,jj)=0.D0    !#/meter**3/sec
         srecom_d(j,jj)=0.D0  !#/meter**3/sec
         scx_d(j,jj)=0.D0     !#/meter**3/sec
         sbcx_d(j,jj)=0.D0    !#/meter**3/sec
         s_d(j,jj)=0.D0       !#/meter**3/sec
         dudtsv_d(j,jj)=0.D0  !#/meter**3/sec
       enddo
c
c---fast ion density, #/m**3, species
       enbeam_d(j)=0.D0     !#/meter**3
c      
c---neutral density,  #/m**3,species
       enn_d(j,1)=0.D0      !#/meter**3
c
c---neutral density from  wall source, #/m**3, species
       ennw_d(j,1)=0.D0     !#/meter**3
c
c---neutral density from volume source, #/m**3,species (recomb and beam cx)
       ennv_d(j,1)=0.D0     !#/meter**3
c
c---sbeam : beam thermal ion source,  #/(m**3*sec) 
       sbeam_d(j)=0.D0      !#/(m**3*sec)
c
c---angular rotation speed profile, rad/sec
       angrot_d(j)=0.D0       !rad/sec
c
c---wdot,electrons, watts/m**3 d(electron energy)/dt profile:
       dpedtc_d(j)=0.D0       !watts/meter**3
c    
c---wdot,ions,watts/m**3 d(ion energy)/dt profile:
       dpidtc_d(j)=0.D0       !watts/meter**3
c
c---electron conduction, watts/m**3
       qconde_d(j)=0.D0        !watts/meter**3
c      
c---ion conduction, watts/m**3
       qcondi_d(j)=0.D0        !watts/meter**3
c
c---electron convection,watts/m**3
       qconve_d(j)=0.D0        !watts/meter**3
c
c---ion convection, watts/m**3
       qconvi_d(j)=0.D0        !watts/meter**3
c
c---power to elec.from beam, watts/m**3
       qbeame_d(j)=0.D0     !watts/meter**3
c
c---power to ions from beam, watts/m**3:
       qbeami_d(j)=0.D0        !watts/meter**3
c    
c---qrfe, rf electron heating, watts/m**3
       qrfe_d(j)=0.D0          !watts/meter**3
c     
c---qrfi, rf ion heating, watts/m**3
       qrfi_d(j)=0.D0          !watts/meter**3
c      
c---qione, recombination and impact ionization, watts/m**3
       qione_d(j)=0.D0         !watts/meter**3
c      
c---qioni, recombination and impact ionization, watts/m**3
       qioni_d(j)=0.D0         !watts/meter**3
c      
c---qcx,  neutral-ion charge exchange, watts/m**3
       qcx_d(j)=0.D0           !watts/meter**3
c     
c---fusion electron heating, watts/m**3
       qfuse_d(j)=0.D0         !watts/meter**3
c               
c---fusion ion heating, watts/m**3
       qfusi_d(j)=0.D0         !watts/meter**3
c
c---qdelt : electron-ion equilibration, watts/m**3
       qdelt_d(j)=0.D0
c
      enddo
c    
c---radiated power density, watts/m**3
       if(i_proc.eq.0) write(6,*) 'Calculating radiated power profile'
       pradsum=0.D0
       if (frad_exp.ne.0) then
       do j=2,nj_d
         dr=r_d(j)-r_d(j-1)
         dvoldr_p=2.D0*pi_m*r_d(j)*2.D0*pi_m*rmajor_d*hcap_d(j)
         dvoldr_m=2.D0*pi_m*r_d(j-1)*2.D0*pi_m*rmajor_d*hcap_d(j-1)
         ds=2.D0*pi_m*(r_d(j)*hcap_d(j)+r_d(j-1)*hcap_d(j-1))/2.D0*dr
         qrad_d(j)=r_d(j)*(frad_exp*r_d(j)**alfrad_exp+frad_exp)
         pradsum=pradsum+
     >     0.5D0*(dvoldr_p*qrad_d(j)+dvoldr_m*qrad_d(j-1))*dr
        write(*,100) j, r_d(j)/r_d(nj_d), qrad_d(j)
       enddo
       qrad_d(1)=r_d(1)*(frad_exp*r_d(1)**alfrad_exp+frad_exp)
       else
       do j=2,nj_d
         qrad_d(j)=0.D0
       enddo
       endif
       if(i_proc.eq.0)
     >  write(6,'(a27,F10.3)') 'Total integrated Prad [MW]:',
     >  pradsum*1.D-6
c
c---total current density, amps/m**2
c...   Assume j(r)=j(0)*(1-rhohat**2)**(alpha)
c...   Note that by approximation (circular geometry)
c...   alpha = qa/qa0-1
c
       if(i_proc.eq.0) write(6,*) 'Calculating j(r) from q-profile'
c
       alpha = qa_exp/q0_exp-1.D0
       sum=0.D0
       do j=2,nj_d
         dr=r_d(j)-r_d(j-1)
         ds=2.D0*pi_m*(r_d(j)*hcap_d(j)+r_d(j-1)*hcap_d(j-1))/2.D0*dr
         curden_d(j) = (1.D0-(r_d(j)/r_d(nj_d))**2)**alpha
         sum = sum + curden_d(j)*ds
       enddo
       cur0 = abs(tocur_d)/sum
c...normalize current density
       do j=1,nj_d
         curden_d(j)=cur0*(1.D0-(r_d(j)/r_d(nj_d))**2)**alpha
       enddo
c
       sum=0.D0
       do j=2,nj_d
         dr=r_d(j)-r_d(j-1)
         ds=2.D0*pi_m*(r_d(j)*hcap_d(j)+r_d(j-1)*hcap_d(j-1))/2.D0*dr
         sum = sum + curden_d(j)*ds
       enddo
       if(i_proc.eq.0)
     >  write(6,'(a27,F10.3)') 'Total integrated Ip [MA]:',
     >  sum*1.D-6
c
c---qdelt : electron-ion equilibration, watts/m**3
c NRL formulary
       peisum=0.D0
       if(qdelt_d(nj_d).eq.0.or.iexp_exch.ge.1) then
         if(i_proc.eq.0) 
     >   write(6,*) 'Calculating electron-ion equilibration'
         do j=1,nj_d
         zbrac=
     > en_d(j,1)/ene_d(j)+apgasa/apimpa*en_d(j,2)/ene_d(j)*apimpz**2
         alamda=24.D0-dlog((ene_d(j)*1.D-6)**0.5D0/(te_d(j)*1.D3))
         tau_e=3.44D5*(te_d(j)*1.D3)**1.5/(ene_d(j)*1.D-6)/alamda
         qdelt_d(j)=-10.D0**6*1.6022D-19*zbrac*3.D0/
     >  (1836.D0*apgasa)/tau_e*ene_d(j)*1.D-6*(te_d(j)-ti_d(j))*1.D3
         enddo
         do j=2,nj_d
         dr=r_d(j)-r_d(j-1)
         dvoldr_p=2.D0*pi_m*r_d(j)*2.D0*pi_m*rmajor_d*hcap_d(j)
         dvoldr_m=2.D0*pi_m*r_d(j-1)*2.D0*pi_m*rmajor_d*hcap_d(j-1)
         peisum=peisum+
     >     0.5D0*(dvoldr_p*qdelt_d(j)+dvoldr_m*qdelt_d(j-1))*dr
         enddo
       endif     !watts/meter**3     
c
       if(iexp_exch.eq.-1) then
        do j=1,nj_d
         qdelt_d(j)=0.D0
        enddo
        peisum=0.D0
       endif
c
       if(i_proc.eq.0) then
        write(6,'(a27,F10.3)') 'Equilibration term zbrac',
     >  zbrac
        write(6,'(a27,F10.3)') 'Total integrated Pei [MW]:',
     >  peisum*1.D-6
       endif
c
c---(electron) ohmic power density,watts/m**3
c...   calculate qohm_d from j(r) and loop voltage
c...   by: qohm = E.j
c...   Also calculate volume-integrated power 
c
       if(i_proc.eq.0) write(6,*) 'Calculating ohmic power profile'
c
       cursum=0.D0
       pohsum=0.D0
       do j=2,nj_d
         dr=r_d(j)-r_d(j-1)
         dvoldr_p=2.D0*pi_m*r_d(j)*2.D0*pi_m*rmajor_d*hcap_d(j)
         dvoldr_m=2.D0*pi_m*r_d(j-1)*2.D0*pi_m*rmajor_d*hcap_d(j-1)
         ds=2.D0*pi_m*(r_d(j)*hcap_d(j)+r_d(j-1)*hcap_d(j-1))/2.D0*dr
         cursum = cursum + curden_d(j)*ds
         qohm_d(j)=(vsurf_d*curden_d(j))/(2.D0*pi_m*rmajor_d)
         pohsum=pohsum+
     >     0.5D0*(dvoldr_p*qohm_d(j)+dvoldr_m*qohm_d(j-1))*dr
       enddo
       qohm_d(1)=(vsurf_d*curden_d(1))/(2.D0*pi_m*rmajor_exp) 
       pohm_d=pohsum
       if(i_proc.eq.0)
     >  write(6,'(a27,F10.3)') 'Total integrated Poh [MW]:',
     >  pohsum*1.D-6
c
c ECH power density, watts/m**3
c
       if(arfe_exp.gt.0) then
       if(i_proc.eq.0) write(6,*) 'Calculating ECH power deposition'
c
       prfesum=0.D0
       do j=2,nj_d
         dr=r_d(j)-r_d(j-1)
         dvoldr_p=2.D0*pi_m*r_d(j)*2.D0*pi_m*rmajor_d*hcap_d(j)
         dvoldr_m=2.D0*pi_m*r_d(j-1)*2.D0*pi_m*rmajor_d*hcap_d(j-1)
         rfexp=rho(j-1)-wrfe_exp
         qrfe_d(j)=arfe_exp*exp(-brfe_exp*(rfexp*rfexp))
         prfesum=prfesum+
     >     0.5D0*(dvoldr_p*qrfe_d(j)+dvoldr_m*qrfe_d(j-1))*dr
c        write(*,100) j, r_d(j)/r_d(nj_d), qrfe_d(j)
       enddo
c      qrfe_d(1)=arfe_exp*exp(-brfe_exp*(r_d(1)
c    >  -wrfe_exp)**2.D0)
	   else
       do j=2,nj_d
         qrfe_d(j)=0.D0
       enddo
       endif
       if(i_proc.eq.0)
     >  write(6,'(a27,F10.3)') 'Total integrated Pech [MW]:',
     >  prfesum*1.D-6
c
c ICRH power density, watts/m**3
c
       if(arfi_exp.gt.0) then
       if(i_proc.eq.0) write(6,*) 'Calculating ICRH power deposition'
c
       prfisum=0.D0
       do j=2,nj_d
         dr=r_d(j)-r_d(j-1)
         dvoldr_p=2.D0*pi_m*r_d(j)*2.D0*pi_m*rmajor_d*hcap_d(j)
         dvoldr_m=2.D0*pi_m*r_d(j-1)*2.D0*pi_m*rmajor_d*hcap_d(j-1)
         rfixp=rho(j)-wrfi_exp
         qrfi_d(j)=arfi_exp*exp(-brfi_exp*(rfixp*rfixp))
         prfisum=prfisum+
     >     0.5D0*(dvoldr_p*qrfi_d(j)+dvoldr_m*qrfi_d(j-1))*dr
       enddo
c      qrfi_d(1)=arfi_exp*exp(-brfi_exp*(r_d(1)
c    >  -wrfi_exp)**2.D0)
	   else
       do j=2,nj_d
         qrfi_d(j)=0.D0
       enddo
       prfisum=0.D0
       endif
       if(i_proc.eq.0)
     >  write(6,'(a27,F10.3)') 'Total integrated Picrh [MW]:',
     >  prfisum*1.D-6
c	     
c---avg major radius of each flux surface, meters, at elevation of mag. axis
c         rmajavnpsi_d(j)   !meters
c
c---avg minor radius of each flux surface,meters, at elevation of mag. axis
c         rminavnpsi_d(j)   !meters
c         
c---volume of each flux surface,  meters**3
c         psivolp_d(j)      !meters**3
c               
c--- elongation of each flux surface
c         elongx_d(j)
c           
c--- flux surface area,  meters**2 4.*pi*pi*R0*hcap*rho*<abs(grad rho)>
c         sfareanpsi_d(j)  !meters**2
c                        
c---flux surface average absolute grad rho  
c        grho1npsi_d(j)
c
c---flux surface average ( grad rho)**2  
c        grho2npsi_d(j)       
c
c ------------------------------------------------------------------
c
c...Neoclassical chii computed outside of subroutine
c       
c ------------------------------------------------
c
c......................................................................       
cmnt   It is assumed that 0:jmaxm greater than or equal 1:nj_dd
cmnt   The volume source data is integrated to power and flows
cmnt   then all the data may be splined to a more dense jmaxm grid

c chk rho(a)
c      do j=2,nj_d
c       rho_a_chk(j)=(psivolp_d(j)-psivolp_d(j-1))*r_d(nj_d)
c    >   /(r_d(j)-r_d(j-1))/(sfareanpsi_d(j)+sfareanpsi_d(j-1))
c    >   *(grho1npsi_d(j)+grho1npsi_d(j-1))
c      enddo
c      
c dimensionally similar discharges
c
      if(ascale.eq.0.) ascale=1.D0
      if(bscale.eq.0.) bscale=1.D0
c    
c      if(iscale.eq.1) pscale=bscale**(1.D0)*ascale**(zhalf)
c      if(iscale.eq.2) pscale=bscale**(5.D0/3.D0)*ascale**(z4thrd)
c
c      iscale=1 gyroBohm
c            =2 Bohm
c            =3 Goldston
c            =4 stochastic
c
       rhostar=bscale**(z2thrd)*ascale**(-5.D0/6.D0)
c
       if(iscale.eq.1) pscale=ascale**(-z3qtr)*rhostar**(-3.D0/2.D0)
       if(iscale.eq.2) pscale=ascale**(-z3qtr)*rhostar**(-5.D0/2.D0)
       if(iscale.eq.3) pscale=ascale**(-z3qtr)*rhostar**(-3.D0)
       if(iscale.eq.4) pscale=ascale**(-z3qtr)*rhostar**(-7.D0/2.D0)
c   
      if(pscale.eq.0) pscale=1.D0
     
c    global transfers
c    read in as inputs
c
c     arho_exp = r_d(nj_d)*ascale 
c     rmajor_exp = rmajor_d*ascale                             
c     bt_exp = abs(btor_d)*bscale                                 
c     elonga_exp=elongx_d(nj_d)
c
c    direct profile transfers
c
      do j=1,nj_d
       te_exp(j-1)=te_d(j)*bscale**(z2thrd)*ascale**(zthird)
       ti_exp(j-1)=ti_d(j)*bscale**(z2thrd)*ascale**(zthird)
       ne_exp(j-1)=1.D-19*ene_d(j)*bscale**(z4thrd)/ascale**(zthird)
       ni_exp(j-1)=1.D-19*en_d(j,1)*bscale**(z4thrd)/ascale**(zthird)
       nz_exp(j-1)=1.D-19*en_d(j,nprim_d+1)*
     >      bscale**(z4thrd)/ascale**(zthird)
       nfst_exp(j-1)=1.D-19*enbeam_d(j)*bscale**(z4thrd)/
     >      ascale**(zthird)
c   we are assuming dln(ni)/dr=dln(ne)/dr 
c   and will use ion plasma sources as electron 
c   plasma sources. zeff enters only into collisionality.
c   ni_m is optained from ne_m assuming constant ni_exp/ne_exp
       zeff_exp(j-1)=zeff_d(j)
       q_exp(j-1)=q_d(j)
       rmin_exp(j-1)=arho_exp*rho(j-1)/dsqrt(elonga_exp)
       rmaj_exp(j-1)=rmajor_exp
       gradrhosq_exp(j-1)=(1.D0+elonga_exp**2)/2.D0/elonga_exp
	   gradrho_exp(j-1)=dsqrt(gradrhosq_exp(j-1))
       angrot_exp(j-1)=angrot_d(j)*
     >     (bscale**(z2thrd)*ascale**(zthird))**(zhalf)/ascale
cc
cc...JAK 6/21/95 To be consistent with "across the street" database, compute
cc...            hcap from surfacearea
cc
c       IF (j.gt.1)
c     >   hcap_d(j)=sfareanpsi_d(j)/
c     >   (grho1npsi_d(j)*r_d(j)*4.D0*(3.1416D0*3.1416D0)*rmajor_d) 
c
       h_exp(j-1)=hcap_d(j)
       elong_exp(j-1)=elonga_exp
c      chiineo_exp(j-1)=xkineo_d(j)/ni_exp(j-1)/1.e19/
c    >    (bscale*ascale**(1./2.))
c    > /gradrhosq_exp(j-1)
c      zptineomax(j-1)=arho_exp*100./zfluxlim_nc(j)
c
      enddo
c 
        do j=1,jmaxm
         shat_exp(j)=(rho(j)+rho(j-1))/(q_exp(j)+q_exp(j-1))*
     >       (q_exp(j)-q_exp(j-1))/(rho(j)-rho(j-1))
c
c... JAK 7/06/95 zero check added
c
         if (ABS(shat_exp(j)).LT.1.D-6) shat_exp(j)=1.D-6
        enddo
         shat_exp(0)=0.D0
c
        jstinv=1
c
        do j=2,jmaxm
         if (q_exp(j).gt.1..and.q_exp(j-1).le.1.) jstinv=j
        enddo
c
        jstmix=dsqrt(2.D0)*jstinv
c
c  surface factor
c
       do j=0,jmaxm
        sfactor(j)=
     >      2.D0*pi_m*arho_exp*rho(j)*h_exp(j)*2.D0*pi_m*rmajor_exp
       enddo
c
c    source integration
c
cmnt the neutral source is broaken into wall and volume
cmnt volume is from beams and recombination where as wall is from wall
cmnt recycled neutrals. The latter are estimated. To chamge the estimate
cmnt the parameter wallneut can be changed from 1.0. The ion and cx 
cmnt energy sources proportional to the neutral density will be corrected
c 
       powe_beam_exp(0)=0.D0
       powe_rf_exp(0)=0.D0
       powe_oh_exp(0)=0.D0
       powe_rad_exp(0)=0.D0
       powe_ion_exp(0)=0.D0
       powe_wdot_exp(0)=0.D0
       powe_fus_exp(0)=0.D0
c           
       powi_beam_exp(0)=0.D0
       powi_rf_exp(0)=0.D0
       powi_ion_exp(0)=0.D0
       powi_cx_exp(0)=0.D0
       powi_wdot_exp(0)=0.D0
       powi_fus_exp(0)=0.D0
c   
       pow_ei_exp(0)=0.D0
c      
       powe_exp(0)=0.D0
       powi_exp(0)=0.D0
c
       flow_wall_exp(0)=0.D0
       flow_beam_exp(0)=0.D0
       flow_sdot_exp(0)=0.D0
       flow_exp(0)=0.D0
c      
       vol_exp(0)=0.D0
c  	   
       do j=2,nj_d       
        dr=r_d(j)-r_d(j-1)
        dvoldr_p=2.D0*pi_m*r_d(j)*2.D0*pi_m*rmajor_d*hcap_d(j)
        dvoldr_m=2.D0*pi_m*r_d(j-1)*2.D0*pi_m*rmajor_d*hcap_d(j-1)      
        vol_exp(j-1)=vol_exp(j-2)+
     >       0.5D0*(dvoldr_p+dvoldr_m)*dr
c             
        powe_beam_exp(j-1)=0.D0
        powe_rf_exp(j-1)=powe_rf_exp(j-2)+ 
     >   1.D-6*0.5D0*(dvoldr_p*qrfe_d(j)+dvoldr_m*qrfe_d(j-1))*dr
        powe_oh_exp(j-1)=powe_oh_exp(j-2)+
     >   1.D-6*0.5D0*(dvoldr_p*qohm_d(j)+dvoldr_m*qohm_d(j-1))*dr
        powe_rad_exp(j-1)=powe_rad_exp(j-2)+ 
     >   1.D-6*0.5D0*(dvoldr_p*qrad_d(j)+dvoldr_m*qrad_d(j-1))*dr
        powe_ion_exp(j-1)=0.D0
        powe_wdot_exp(j-1)=0.D0
        powe_fus_exp(j-1)=0.D0   
c
        powi_beam_exp(j-1)=0.D0
        powi_rf_exp(j-1)=powi_rf_exp(j-2)+ 
     >   1.D-6*0.5D0*(dvoldr_p*qrfi_d(j)+dvoldr_m*qrfi_d(j-1))*dr
        powi_ion_exp(j-1)=0.D0
        powi_cx_exp(j-1)=0.D0
        powi_wdot_exp(j-1)=0.D0
        powi_fus_exp(j-1)=0.D0
c  
        pow_ei_exp(j-1)=pow_ei_exp(j-2)-
     >   1.D-6*0.5D0*(dvoldr_p*qdelt_d(j)+dvoldr_m*qdelt_d(j-1))*dr
c
        powe_exp(j-1)=powe_oh_exp(j-1)-pow_ei_exp(j-1)
     >   -xrad_exp*powe_rad_exp(j-1)+powe_rf_exp(j-1)
        powi_exp(j-1)=pow_ei_exp(j-1)+powi_rf_exp(j-1)
        if(powi_exp(j-1).eq.0) powi_exp(j-1)=1.D-6
        flow_exp(j-1)=0.D0
c
       enddo
c	   
c pow_e,i_exp (total power) in MW= MJ/sec
c flow_exp in MW/kev=kamps
c
c      do j=0,jmaxm
c        xxq=rho(j)**2
c        powe_exp(j)=rho(j)*(powe0_exp-powe0_exp*(1.D0-xxq)**alfpowe_exp)       
c        powi_exp(j)=rho(j)*(powi0_exp-powi0_exp*(1.D0-xxq)**alfpowi_exp)       
c        flow_exp(j)=rho(j)*(flow0_exp-flow0_exp*(1.D0-xxq)**alfflow_exp)
c      enddo
c
       do j=2,nj_d
        vol_exp(j-1)=vol_exp(j-1)*ascale**3
        powi_exp(j-1)=powi_exp(j-1)*pscale
        powe_exp(j-1)=powe_exp(j-1)*pscale
        pow_ei_exp(j-1)=pow_ei_exp(j-1)*pscale
        powi_beam_exp(j-1)=powi_beam_exp(j-1)*pscale
        flow_beam_exp(j-1)=flow_beam_exp(j-1)*pscale
c note if iexp_exch=-1 pow_ei_exp=0 and pow_ei_cor is the actual
c corrected exchange
c
        flow_exp(j-1)=flow_exp(j-1)*
     >     pscale/(bscale**(z2thrd)*ascale**(zthird))
       enddo        
c
c flow_exch_exp(j) is the expansion cooling in the ion channel integrated
c to a power flow (MW); negative ie cooling for positive plasma flow.
c It estimates the possible size of "anomalous" energy exchange.  It should
c be compared with pow_ei_exp(j).
c For electron directed waves we expect electron channel cooling and ion 
c channel heating by this amount       
c For ion directed waves we expect electron channel heating  and ion 
c channel cooling by this amount
c
       flow_exch_exp(0)=0.D0
       do j=1,nj_d-1
        flow_exch_exp(j)=flow_exch_exp(j-1)+
     >      dr*flow_exp(j)*ti_exp(j)*(dlog(ti_exp(j)*ni_exp(j))
     >      -dlog(ti_exp(j-1)*ni_exp(j-1)))/dr
       enddo
c
c noncircular geometric factor
c
       do j=1,jmaxm-1
        geofac(j)=gradrho_exp(j)*(rho(j+1)-rho(j))*arho_exp
     >   /(rmin_exp(j+1)-rmin_exp(j))/gradrhosq_exp(j)
        drhodr(j)=(rho(j+1)-rho(j))*arho_exp/
     >   (rmin_exp(j+1)-rmin_exp(j))
        drhodrrrho(j)=drhodr(j)*rmin_exp(j)/arho_exp/rho(j)
        shift_shaf=(rmaj_exp(j+1)-rmaj_exp(j))/
     >   (rmin_exp(j+1)-rmin_exp(j))
        georotrate(j)=elong_exp(j)/(1.D0+shift_shaf)**2
        geoalpha(j)=1.D0/rmajor_exp*(vol_exp(j+1)-vol_exp(j))
     >   /(rho(j+1)-rho(j))/arho_exp/(2.D0*pi_m*rho(j)*arho_exp)**2
     >  *(vol_exp(j)/(2.D0*pi_m**2*rmajor_exp))**0.5D0
       enddo
       georotrate(jmaxm)=georotrate(jmaxm-1)
       georotrate(0)=georotrate(1)
       geoalpha(jmaxm)=geoalpha(jmaxm-1)
       geoalpha(0)=geoalpha(1)
       geofac(jmaxm)=geofac(jmaxm-1)
       drhodr(jmaxm)=drhodr(jmaxm-1)
       drhodrrrho(jmaxm)=drhodrrrho(jmaxm-1)
       geofac(0)=geofac(1)
       drhodr(0)=drhodr(1)
       drhodrrrho(0)=drhodrrrho(1)    
c07/16/98
c noncircular geometric factor
c
       do j=1,jmaxm
        geofac(j)=gradrho_exp(j)*(rho(j)-rho(j-1))*arho_exp
     >   /(rmin_exp(j)-rmin_exp(j-1))/gradrhosq_exp(j)
        drhodr(j)=(rho(j)-rho(j-1))*arho_exp/
     >   (rmin_exp(j)-rmin_exp(j-1))
        drhodrrrho(j)=drhodr(j)*rmin_exp(j)/arho_exp/rho(j)
        shift_shaf=(rmaj_exp(j)-rmaj_exp(j-1))/
     >   (rmin_exp(j)-rmin_exp(j-1))
        georotrate(j)=elong_exp(j)/(1.D0+shift_shaf)**2
        geoalpha(j)=1./rmajor_exp*(vol_exp(j)-vol_exp(j-1))
     >   /(rho(j)-rho(j-1))/arho_exp/(2.D0*pi_m*rho(j)*arho_exp)**2
     >  *(vol_exp(j)/(2.D0*pi_m**2*rmajor_exp))**0.5D0
c        write(*,100) j, rho(j), rmaj_exp(j), vphi_exp(j)
       enddo
       georotrate(0)=georotrate(1)
       geoalpha(0)=geoalpha(1)
       geofac(0)=geofac(1)
       drhodr(0)=drhodr(1)
       drhodrrrho(0)=drhodrrrho(1)
c
c... jek added 999 return and error message
c
c      return
c      end
c
 999   return
c
 89    continue
       if(i_proc.eq.0) write(*,*) 'Error reading iter namelist'
       close(50)
       goto 999
c
 100  format(i2,2x,0p1f4.2,1p8e13.5)
       end 
