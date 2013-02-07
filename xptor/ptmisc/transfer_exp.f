      subroutine transfer_exp(utime,istep)
****************************************************************
* 08-Aug-00 jek freeze profiles if simulation time > xp_time
* Transfers data stored in u*d arrays to *_exp arrays
* derived from rep_iter.f G.M. Staebler 6/28/2000
* istep used so printout only on first time-step in ptor
*
****************************************************************
      implicit none
      integer nrmax, ntmax, n1d, n2d, n3d, nur, istep
      parameter(nrmax=51) !max number of radial points in experimental grid
      parameter(ntmax=801) !max number of time points in experimental grid
      parameter(n1d=13)   !number of 1d fields read;
      parameter(n2d=44)   !number of 2d fields read;
      parameter(n3d=6)   !number of fields read for experimental data 
c
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/glf.m'
      include '../inc/ptor.m'
      include '../inc/ut.m'
c
      integer i, ig, is, i2d, ierr, time_channel
     &  , j, i1d, id, ifchar, ichar, ifail, k
     &  , nnexp_d, ntexp_d, ntixp_d, jj, jn, niterdb, lnblnk, ier
c
      real*8 utime
     &  , rfix, sfix, q01_exp, qa1_exp, zeff_t, zbrac, alamda
     &  , tau_e, ds, z2, alfj1_exp, xxq, alpha, sum, cur0, cursum
     &  , rhostar, shift_shaf, dummy1, aomega
     &  , zhalf, zthird, z2thrd, z4thrd, z3qtr
c
c      real*8 xnexp_d(mxgrid+1), ynexp_d(mxgrid+1)
c     &  , xnexpeb_d(mxgrid+1), ynexpeb_d(mxgrid+1)
c     &  , xtexp_d(mxgrid+1), ytexp_d(mxgrid+1)
c     &  , xtexpeb_d(mxgrid+1), ytexpeb_d(mxgrid+1)
c     &  , xtixp_d(mxgrid+1), ytixp_d(mxgrid+1)
c     &  , xtixpeb_d(mxgrid+1), ytixpeb_d(mxgrid+1)
       real*8 eps, sigmas
       real*8 wt(nj), ys(nj),ysp(nj)
       real*8 td(nj),tsd1(nj),hd(nj),hsd1(nj),hsd2(nj),
     &        rd(nj),rsd1(nj),rsd2(nj),v(nj)
c
c Need to input a list of global quantities that characterize the shot. 
c I will assume the following typical set might be available:
c     pgasa      plasma working gas atomic number
c     pgasz      plasma working gas charge
c     bgasa      beam species atomic number
c     bgasz      beam species charge
c     pimpa      impurity species atomic number
c     pimpz      impurity species charge
c 
c Need to input a list of global quantities vs. time
c I will assume the following typical set might be available:
c
c     amin       horizontal minor radius (m)
c     bt         vacuum toroidal field (T) at RGEO
c     delta      triangularity 
c     indent     indentation 
c     ip         plasma current (A)
c     kappa      elongation
c     li         induction 
c     q95        plasma safety factor evaluated at the flux surface that 
c                encloses 95% of the total poloidal flux
c     rgeo       plasma geometrical major radius (m)
c     wtot       total stored energy
c     zeff       ave zeff
c
c Need to input a list of profiles vs. time that are needed.
c I will assume the following typical set might be available:
c     chie         experimentally inferred chi_e
c     chii         experimentally inferred chi_i
c     curnbi       NBI currents??
c     curtot       Total current??
c     grho1        < abs(grad rho)>
c     grho2        <(grad rho)**2>
c     ne           electron density
c     nfast        fast ion density
c     nimp         low-Z impurity density??
c     nm1          main ion density??
c     qnbie        power from NBI to electrons
c     qnbii        power from NBI to ions
c     qohm         ohmic heating power
c     qrad         radiated power
c     rmajor       avg major radius of flux surface
c     rminor       avg minor radius of flux surface
c     snbie        electron particle source from NBI
c     snbii        ion particle source from NBI
c     swall        thermal ion particle source due to ionisation
c     surf         flux surface surface area
c     te           T_e
c     ti           T_i
c     volume       flux surface volume
c     zeffr        Z_eff profile
cgms  ptot         total pressure profile
c
c
      integer iread, iflag
      real*8 rho_norm_loc, aj
c     real*8 xp_time
      real*8 lnlam,eta,pohsum,pradsum,drm,dvoldr_p,dvoldr_m
      real*8 u1d(n1d),u2d(mxgrid+1,n2d)
c      real*8 u3dx(mxgrid+1,n3d),u3dy(mxgrid+1,n3d)
      real*8 prfesum, prfisum
c      integer nx3(n3d)
      integer time_flag, itt, itime_shift
c      character*7 fields1d(n1d)
c      character*7 fields2d(n2d)
c      character*7 fields3d(n3d)
c      character*50 ufile
c      character*100 udfile
c     character*50 cudir
c      character*67 u0names(200)
c      character*11 u0values(200)
c
c      character*6 u0phase
      integer pgasa,pgasz,bgasa,bgasz
c
      real, allocatable, dimension(:,:) :: rhob_nc
      real, allocatable, dimension(:) :: zfluxlim_nc, xnstari_nc,
     &      xnstare_nc, xkapi_nc
      save iflag
      allocate (rhob_nc(nj,3),zfluxlim_nc(nj),xnstari_nc(nj),
     &          xnstare_nc(nj), xkapi_nc(nj))
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c...jek misc switches
c
      rho_norm_loc=1.D0
c     jmaxm=50
      jmaxm=mxgrid
      nr = mxgrid+1
c     iexp_exch=-1
      imyneoclass=0   ! Use neoclassical transport in ufile
c     nkimod_nc=2
c     aimp_nco=12
c     xzimp_nco=6
c      amassgas_exp=2.D0
      pi_m=3.1415926D0
      kevdsecpmw=1.6022D-19*1.D3*1.D-6
      zhalf = 1.D0/2.D0
      zthird = 1.D0/3.D0
      z2thrd = 2.D0/3.D0
      z4thrd = 4.D0/3.D0
      z3qtr = 3.D0/4.D0
c
c...  Note that here p_glob_exp can still be .le.0
c...  In that case it is assigned a value (according
c...  to 1D data) just before ITER89P is calculated.
c
c     Assume the working and beam ions are same and the impurity
c     is carbon for the moment. ** no longer used **
c     The masses and charges are set in rep_iter using either
c     the 0D file values, user input values, or default values
c     given in mltin.
c
c     pgasa=amassgas_exp
c     pgasz=1
c     bgasa=pgasa
c     bgasz=1
c     pimpa=12
c     pimpz=6
c
c     apgasa=pgasa
c     apgasz=pgasz
c     abgasa=bgasa
c     abgasz=bgasz
c     apimpa=pimpa
c     apimpz=pimpz
c
c..interpolated 1D and 2D variables
c
      if(utime.le.ut(nut)) then
      do k=1,nr
        do i=1,n2d
          u2d(k,i) = 0.D0
        enddo
      enddo
      do i=1,n1d
        u1d(i) = 0.D0
      enddo
      endif
c
c loop over time of data
c
      if(utime.le.ut(nut)) then
      do j=2,nut
        if (ut(j).ge.utime .and. ut(j-1).lt.utime) then
c         write(*,*)'utime = ',utime,'ut = ',ut(j),ut(j-1),j
c loop over 2d data
           do i=1,n2d
c interpolate to the time of interest
            do k=1,nr
c              u2d(k,i) = u2d_t(k,j,i)
              u2d(k,i) = u2d_t(k,j-1,i) + 
     >                   (u2d_t(k,j,i)-u2d_t(k,j-1,i))*
     >                   (utime-ut(j-1))/(ut(j)-ut(j-1))
            enddo
           enddo
c loop over 1d data
           do i=1,n1d  
c interpolate to the time of interest
c             u1d(i) = u1d_t(j,i)
             u1d(i) = u1d_t(j-1,i) + (u1d_t(j,i)-u1d_t(j-1,i))*
     >                (utime-ut(j-1))/(ut(j)-ut(j-1))
           enddo
           iflag=0
         endif   
       enddo
       endif
c
c..check to see if ptor time > data time
c  if so, then give warning and freeze profiles at last diagnostic time
c
      if (utime.gt.ut(nut)) then
        iflag=iflag+1
        if (iflag.eq.1 .and. i_proc.eq.0) write(*,500)
           do i=1,n2d
            do k=1,nr
              u2d(k,i)=u2d_t(k,nut,i)
            enddo
           enddo
      endif
c
c...jek check on Ti
c     write(*,*) 'nr,nut = ',nr,nut
c     do k=1,nr
c       write(*,53) k, rho_d(k), u2d_t(k,nut,22), u2d_t(k,nut-1,22)
c     enddo
c     write(*,*) (u2d(j,22)/1.D3,j=1,nr)
c          
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c
c ishot : shot number          
cjak      ishot_d=10000
c
c nj : the size of the vectors printed in this file
          nj_d=nr
c
c nion : the number of ion species          
          nion_d=2
c
c nprim : the number of primary ion species
          nprim_d=1
c
c nimp : the number of impurity ion species          
          nimp_d=1
c
c nneu : the number of neutral ion species
          nneu_d=1
c
c ibion : index of beam species          
          ibion_d=1
c
c namep : name(s) of primary ion species          
c         (namep_d(i),i=1,nprim_d)
c
c namei : name(s) of impurity ion species          
c         (namei_d(i),i=1,nimp_d)
c
c namen : name(s) of neutral ion species         
c         (namen_d(i),i=1,nneu_d)
c              
c time :  time at which data is printed          
          time_d=utime
c
c Rgeom : major radius of geometric          
         rgeom_d=u1d(9)                  !meters
c 
c Rmag :  major radius of mag axis, meters        
         rmag_d=u2d(1,15)                    !meters
         if (rmag_d.eq.0.) rmag_d=u1d(9)
c
c R0 : major radius of vaccuum btor ref            
         rmajor_d=u1d(9)                  !meters
         amin_d=u1d(1)                  !meters
c
c JAK 960325 warning output added
         if (rmajor_d.eq.0.0 .and. istep.eq.1)
     >     write(6,*) 'Warning transfer_exp: rmajor zero'
c
c kappa : plasma elongation         
         if(u1d(6).gt.0.) then
           kappa_d=u1d(6)
         else
           u1d(6)=1.D0
           kappa_d=1.D0
         endif
c delta : plasma triangularity         
         deltao_d=u1d(3)
c  
c pindent : plasma indentation         
         pindento_d=u1d(4)
c
c volo : plasma volume,meters**3         
         volo_d=u2d(nr,21)                    !meters**3
c
c Btor : vaccuum toroidal field at rmajor, tesla         
         btor_d=u1d(2)                    !tesla
c
c total,ohmic,bootstrap,beam,and rf currents, amps       
c        tocur_d
c        totohm_d
c        totboot_d
c        totbeam_d
c        totrf_d
c
c totohm_d: total ohmic current               !A        (JAK)
c JAK CHECK THIS: IP is deterined from an external Rogowski loop
c does this include bootstrap and non-inductive components?
c        totohm_d=u1d(4)
         tocur_d=u1d(5)
         curtot=tocur_d/1.D6
c
c betap : poloidal beta         
c         betap_d
c
c beta : toroidal beta         
c         beta_d
c
c ali : plasma inductance         
         ali_d=u1d(7)
c
c te0 : central electron temperature         
         te0_d=u2d(1,23)/1.D3                        !kev
c
c ti0 : central ion temperature         
         ti0_d=u2d(1,22)/1.D3                        !kev
c
c pohm: total ohmic input power                      !Wa  (JAK)
         pohm_d=u1d(12)
c
c vsurf: loop voltage at the plasma boundary         !V   (JAK)
         vsurf_d=u1d(13)
c
c---psi on rho grid,volt*sec/rad 
c              
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(psir_d(j),j=1,nj_d)      !volt*se/rad        
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c   
c---rho grid, meters
c
c get rho_a from edge kappa
c r_d = rho*a*sqrt(kappa)
          do j=1,nj_d
             r_d(j)=rho_d(j)*u1d(1)*sqrt(u1d(6))    !meters
             r_d(j)=r_d(j)*rho_norm_loc
	  enddo
c JAK 960325 warning output added
          if (r_d(nj_d).eq.0.0) then
            if(i_proc.eq.0 .and. istep.eq.1) 
     >         write(6,*) 'Warning rep_iter: r_d zero'
          endif
c
c---fcap, (ie f(psilim)/f(psi) )
c
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(fcap_d(j),j=1,nj_d)
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c
c---gcap, (ie <(grad rho)**2*(R0/R)**2> )
c
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(gcap_d(j),j=1,nj_d)
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c       
c---hcap, (ie (surface)/(4*pi*pi*R0*gradrho*rho)) ---JAK
c
           rfix=1.D0
           if(istk.eq.88) rfix=4.46
           sfix=1.D0
           if(istk.eq.88) sfix=22.825
           if(istk.eq.88) write(6,*) 'sfix=22.825'
c
          do j=2,nj_d
            hcap_d(j)=sfix*u2d(j,19)/(u2d(j,5)*
     >      r_d(nj_d)/rfix)/(4.D0*pi_m*pi_m*rmajor_d*r_d(j))
          enddo
          hcap_d(1)=hcap_d(2)
c
c---read original experimental data for ne, te and ti ---JAK
c
c ne and errorbar
c
c          nnexp_d = nx3(1)
c          do j=1,nnexp_d
c           xnexp_d(j)=u3dx(j,1)
c           ynexp_d(j)=u3dy(j,1)*1.D-19
c           xnexpeb_d(j)=u3dx(j,2)
c           ynexpeb_d(j)=u3dy(j,2)*1.D-19
c          enddo
c
c te and errorbar
c
c          ntexp_d = nx3(3)
c          do j=1,ntexp_d
c           xtexp_d(j)=u3dx(j,3)
c           ytexp_d(j)=u3dy(j,3)*1.D-3        !keV
c           xtexpeb_d(j)=u3dx(j,4)
c           ytexpeb_d(j)=u3dy(j,4)*1.D-3      !keV
c          enddo
c
c ti and errorbar
c
c          ntixp_d = nx3(5)
c          do j=1,ntixp_d
c           xtixp_d(j)=u3dx(j,5)
c           ytixp_d(j)=u3dy(j,5)*1.D-3        !keV
c           xtixpeb_d(j)=u3dx(j,6)
c           ytixpeb_d(j)=u3dy(j,6)*1.D-3      !keV
c          enddo
c
c---electron temperature, kev
c
          do j=1,nj_d
           te_d(j)=u2d(j,23)/1.D3            !kev
          enddo
c
c---ion temperature, kev
c
          do j=1,nj_d
           ti_d(j)=u2d(j,22)/1.D3          !kev
          enddo
c
       if( iexp_exch.eq.2) then
        do j=1,nj_d
         ti_d(j)=te_d(j)-.001
        enddo
       endif    
c
c---q (ie safety factor) profile
c
c... JAK 7/20/95 section to calculate q-profile if not available
c... moved to section q_ohm (since elongation is needed here 
c... that is read later) 
          do j=1,nj_d
           q_d(j)=u2d(j,27)
          enddo
          q01_exp=q_d(1)
          qa1_exp=q_d(nj_d)
c
c---electron density,#/m**3
c
          do j=1,nj_d
           ene_d(j)=u2d(j,7)           !#/meter**3
          enddo                         
c      
c---primary ion density,#/m**3,species (NM1,NM2,NM3)
c---impurity ion density,#/m**3,species (NIMP)
c   assuming abgasz=apgasz
            do j=1,nj_d
               en_d(j,1)=u2d(j,10) + u2d(j,42) + u2d(j,43) !#/meter**3
c              if (imodel.lt.0) then
c                en_d(j,1)=u2d(j,10) + u2d(j,42) + 
c    >                     u2d(j,43) + u2d(j,9) ! chi=1 test
c              endif
               en_d(j,2)=u2d(j,9)    !#/meter**3
               zeff_t=u2d(j,20)
               if (zeff_t.eq.0.) zeff_t=u1d(11)
               if (zeff_t.eq.0.) then
                 zeff_t=1.
                 if(i_proc.eq.0 .and. istep.eq.1) 
     >              write(6,*) 'zeff_t=1.'
               endif
c
               if(en_d(j,1).eq.0.) 
     >          en_d(j,1)=(apimpz-zeff_t)/(apimpz-apgasz)/apgasz
     >                    *u2d(j,7)-u2d(j,8)
               if(en_d(j,2).eq.0.) 
     >          en_d(j,2)=(zeff_t-apgasz)/(apimpz-apgasz)/apimpz
     >		              *u2d(j,7)
               if(en_d(j,1).lt.0..or.en_d(j,2).lt.0.) then
                 if(i_proc.eq.0 .and. istep.eq.1) 
     >              write(6,*) 'negative density fixed'
                 en_d(j,1)=dabs(en_d(j,1))
                 en_d(j,2)=dabs(en_d(j,2))
               endif
             enddo
c
c---sion   : source due to ionization,        #/(m**3*sec),species
c---srecom : source due to recombination,     #/(m**3*sec),species
c---scx    : source due to cx thermal neut.,  #/(m**3*sec),species
c---sbcx   : sink due to cx with beam neut.   #/(m**3*sec),species
c---s      : total source rate,               #/(m**3*sec),species
c---dudt   : s dot,                           #/(m**3*sec),species
c
       do jj=1,nprim_d       
            do j=1,nj_d
              sion_d(j,jj)=0.
              srecom_d(j,jj)=0.
              scx_d(j,jj)=0.
              sbcx_d(j,jj)=0.
              s_d(j,jj)=0.
              dudtsv_d(j,jj)=0.
             enddo
          enddo
c
c---fast ion density, #/m**3, species
c
          do j=1,nj_d
             enbeam_d(j)=u2d(j,8)         !#/meter**3
          enddo
c      
c---neutral density,  #/m**3,species
c
       do jn=1,nneu_d
           do j=1,nj_d
              enn_d(j,jn)=0.     !#/meter**3
          enddo
       enddo
c
c---neutral density from  wall source, #/m**3, species
c
       do jn=1,nneu_d        
c               read(niterdb,'(a)')stflg
c               read(niterdb,10)(ennw_d(j,jn),j=1,nj_d)     !#/meter**3     
         do j=1,nj_d
          ennw_d(j,jn)=0.
         enddo
       enddo
c
c---neutral density from volume source, #/m**3,species (recomb and beam cx)
c
       do jn=1,nneu_d
c               read(niterdb,'(a)')stflg
c               read(niterdb,10)(ennv_d(j,jn),j=1,nj_d)     !#/meter**3
        do j=1,nj_d
         ennv_d(j,jn)=0.
        enddo
       enddo
c
c---volume source of neutrals, #/(m**3*sec),species
c
c       do jn=1,nneu_d   
c               read(niterdb,'(a)')stflg
cx               read(niterdb,10)(volsn_d(j,jn),j=1,nj_d)!#/(m**3*sec)
c               read(niterdb,10)(blank_d(j),j=1,nj_d) 
c       enddo
c
c---sbion : beam electron source,  #/(m**3*sec
c        
c               read(niterdb,'(a)')stflg
cx               read(niterdb,10)(sbion_d(j),j=1,nj_d)    !#/(m**3*sec)
c               read(niterdb,10)(blank_d(j),j=1,nj_d)          
c
c---sbeam : beam thermal ion source,  #/(m**3*sec)
c
c               read(niterdb,'(a)')stflg   
c               read(niterdb,10)(sbeam_d(j),j=1,nj_d)      !#/(m**3*sec)
c
c...JAK 6/19/95 changed to read from datafile: sbeam_d and sion_d
c...Note that the total source rate (s) and s dot (dudt) are not 
c...present in the database.
c
c...Note that SNBII as read contains all sources and sinks due to
c...beam thermal ion source and charge exchange processes; 
c...so: field SNBII = sbeam + scx + sbcx + srecom    (note sbcx<0 sink)
c...Read the source in sbeam_d and let all other source terms be 0.          
c
             do j=1,nj_d
              sbeam_d(j)=u2d(j,18)                       !#/(m**3*sec)
             enddo
c
c...thermal ion particle source due to ionisation; field SWALL
c
             do j=1,nj_d
              sion_d(j,1)=u2d(j,30)                       !#/(m**3*sec)
             enddo
c dn/dt
             do j=1,nj_d
               dudtsv_d(j,1)=u2d(j,41)
             enddo
c
c---total current density, amps/m**2
c     
c... JAK 6/20/95 read from field CURTOT
c... JAK 8/21/95 avoid negative current density
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(curden_d(j),j=1,nj_d)        !amps/meter**2
             do j=1,nj_d
              curden_d(j)=DABS(u2d(j,4))
             enddo
c    
c---ohmic current density, amps/m**2
c     
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(curohm_d(j),j=1,nj_d)        !amps/meter**2
c          read(niterdb,10)(blank_d(j),j=1,nj_d) 
c          
c---bootstrap current density, amps/m**2
c      
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(curboot_d(j),j=1,nj_d)       !amps/meter**2
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c      
c--- beam driven current density, amps/m**2
c     
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(curdbeam_d(j),j=1,nj_d)      !amps/meter**2
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c               
c---rf current density, amps/m**2
c    
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(currf_d(j),j=1,nj_d)         !amps/meter**2
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c          
c---rho*bp0*fcap*gcap*hcap, tesla*meters
c     
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(rbp_d(j),j=1,nj_d)
c          read(niterdb,10)(blank_d(j),j=1,nj_d)         !tesla*meters      
c
c---zeff profile:
c
c FIX THIS
          do j=1,nj_d
           zeff_d(j)=u2d(j,20)
           if (zeff_d(j).eq.0.) zeff_d(j)=u1d(11)
          enddo
c
c---angular rotation speed profile, rad/sec
c
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(angrot_d(j),j=1,nj_d)       !rad/sec
          do j=1,nj_d
            angrot_d(j)=u2d(j,29)
          enddo
c
c---electron thermal diffusivity, meters*2*/sec on half grid
c
          do j=1,nj_d
           chieinv_d(j)=u2d(j,1)       !meters**2/sec
          enddo
c                 
c---ion thermal diffusivity,meters*2*/sec  on half grid
c     
          do j=1,nj_d
           chiinv_d(j)=u2d(j,2)       !meters**2/sec
          enddo
c         
c---ion neocl.  thermal conductivity, 1/(m*sec) on half grid 
c        
          do j=1,nj_d
           xkineo_d(j)=0.       !meters**2/sec
          enddo
c
c...JAK 7/20/95 calculate neoclassical aafter section where q-profile
c...is recalculated
c
c---wdot,electrons, watts/m**3 d(electron energy)/dt profile:
c       
c... JAK 6/20/95 read from field DWER
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(dpedtc_d(j),j=1,nj_d)        !watts/meter**3      
          do j=1,nj_d
           dpedtc_d(j)=u2d(j,31)
          enddo
c    
c---wdot,ions,watts/m**3 d(ion energy)/dt profile:
c       
c... JAK 6/20/95 read from field DWIR
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(dpidtc_d(j),j=1,nj_d)        !watts/meter**3
          do j=1,nj_d
           dpidtc_d(j)=u2d(j,32)
          enddo
c
c---electron conduction, watts/m**3
cc      
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(qconde_d(j),j=1,nj_d)        !watts/meter**3 
          do j=1,nj_d
           qconde_d(j)=0.
          enddo
c   
c---ion conduction, watts/m**3
c
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(qcondi_d(j),j=1,nj_d)        !watts/meter**3
          do j=1,nj_d
           qcondi_d(j)=0.
          enddo
c
c---electron convection,watts/m**3 
c     
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(qconve_d(j),j=1,nj_d)        !watts/meter**3
          do j=1,nj_d
           qconve_d(j)=0.
          enddo
c
c---ion convection, watts/m**3 
c     
c         read(niterdb,'(a)')stflg
c         read(niterdb,10)(qconvi_d(j),j=1,nj_d)        !watts/meter**3
          do j=1,nj_d
           qconvi_d(j)=0.
          enddo
c
c---power to elec.from beam, watts/m**3 
c     
           do j=1,nj_d
             qbeame_d(j)=u2d(j,11)        !watts/meter**3
           enddo
c          
c---qdelt : electron-ion equilibration, watts/m**3
c
         do j=1,nj_d
          qdelt_d(j)=u2d(j,28)
         enddo
c      
       if(qdelt_d(nj_d).eq.0.or.iexp_exch.ge.1) then       
         do j=1,nj_d
crew          zeff_t=u2d(j,20)
crew          if (zeff_t.eq.0.) zeff_t=u1d(11)
crew         zbrac=zeff_t/2.
crewc        zbrac=sum_i n_i Z_i**2/(A_i*ne)
crew            qdelt_d(j) = alog(te_d(j)*3.95e16/sqrt(ene_d(j)))*
crew     >      1.503e-19*ene_d(j)**2/te_d(j)**1.5*
crew     >      zbrac*(te_d(j)-ti_d(j))*1.6021e-16
c
c NRL formulary 
         zbrac=en_d(j,1)/ene_d(j)+apgasa/
     >         apimpa*en_d(j,2)/ene_d(j)*apimpz**2
         alamda=24.D0-dlog((dabs(ene_d(j))*1.D-6)**zhalf/
     >         (dabs(te_d(j))*1.D3))
         tau_e=3.44D5*(dabs(te_d(j))*1.D3)**1.5D0/
     >         (ene_d(j)*1.D-6)/alamda
         qdelt_d(j)=-10**6*1.6022D-19*zbrac*3.D0/(1836.D0*apgasa)/
     >         tau_e*ene_d(j)*1.D-6*(te_d(j)-ti_d(j))*1.D3
         enddo                    
       endif     !watts/meter**3
       if(iexp_exch.eq.-1) then
        do j=1,nj_d
         qdelt_d(j)=0.
        enddo
       endif
c
c---power to ions from beam, watts/m**3
c   
           do j=1,nj_d
             qbeami_d(j)=u2d(j,12)        !watts/meter**3
           enddo
c    
c---qrfe, rf electron heating, watts/m**3
c   
c... JAK 6/20/95 read from fields QECHE + QICRHE
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(qrfe_d(j),j=1,nj_d)          !watts/meter**3 
           do j=1,nj_d
            qrfe_d(j)=u2d(j,33)*echconv + u2d(j,34)
           enddo
c    
c---qrfi, rf ion heating, watts/m**3
c
c... JAK 6/20/95 read from fields QECHI + QICRHI
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(qrfi_d(j),j=1,nj_d)          !watts/meter**3 
           do j=1,nj_d
            qrfi_d(j)=u2d(j,35) + u2d(j,36)
           enddo
c
c... JEK calculate integrated RF heating
c
       prfesum=0.D0
       prfisum=0.D0
       do j=2,nj_d
         drm=r_d(j)-r_d(j-1)
         dvoldr_p=2.D0*pi_m*r_d(j)*2.D0*pi_m*rmajor_d*hcap_d(j)
         dvoldr_m=2.D0*pi_m*r_d(j-1)*2.D0*pi_m*rmajor_d*hcap_d(j-1)
         ds=2.D0*pi_m*(r_d(j)*hcap_d(j)+r_d(j-1)*
     >      hcap_d(j-1))/2.D0*drm
         prfesum=prfesum+
     >     zhalf*(dvoldr_p*qrfe_d(j)+dvoldr_m*qrfe_d(j-1))*drm
         prfisum=prfisum+
     >     zhalf*(dvoldr_p*qrfi_d(j)+dvoldr_m*qrfi_d(j-1))*drm
       enddo
c       if(i_proc.eq.0) then
c         write(6,'(a27,2F10.3)') 'Total integrated Prfe,Prfi [MW]:',
c     >        prfesum*1.D-6,prfisum*1.D-6
c       endif
c
c---qione, recombination and impact ionization, watts/m**3
c
c... JAK 6/20/95 read from field QWALLE
c... --- is this correct??? QWALLE = ''thermal electron heat loss due
c... to the ionization of wall neutralsa W/m3''
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(qione_d(j),j=1,nj_d)         !watts/meter**3
           do j=1,nj_d
            qione_d(j)=u2d(j,37)
           enddo
c      
c---qioni, recombination and impact ionization, watts/m**3
c
c... JAK 6/20/95 read from field QWALLI
c... Note that QWALLI is the ''main thermal ion heat loss due to
c... ionization and charge exchange with wall neutrals in W/m3''
c... SO: QWALLI=qioni_d + qcx_d;
c... Put value in qioni_d       
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(qioni_d(j),j=1,nj_d)         !watts/meter**3
           do j=1,nj_d
            qioni_d(j)=u2d(j,38)
           enddo
c      
c---qcx,  neutral-ion charge exchange, watts/m**3
c
c... JAK 6/20/95 see comments above aat qioni_d    
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(qcx_d(j),j=1,nj_d)           !watts/meter**3
           do j=1,nj_d
            qcx_d(j)=0.
           enddo
c     
c---2d MHD equil. electron heating, watts/m**3
c      
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qe2d_d(j),j=1,nj_d)        !watts/meter**3
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c         
c---2d MHD equil.ion heating, watts/m**3
c      
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qi2d_d(j),j=1,nj_d)         !watts/meter**3
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c         
c---fusion electron heating, watts/m**3
c   JAK 96/03/20 added

           do j=1,nj_d
            qfuse_d(j)=u2d(j,39)
           enddo
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qfuse_d(j),j=1,nj_d)        !watts/meter**3
c          read(niterdb,10)(blank_d(j),j=1,nj_d) 
c               
c---fusion ion heating, watts/m**3
c
c   JAK 96/03/20 added
           do j=1,nj_d
            qfusi_d(j)=u2d(j,40)
           enddo    
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qfusi(j),j=1,nj_d)        !watts/meter**3
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c         
c---beam fusion electron heating,watts/m**3
c   
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qbfuse_d(j),j=1,nj_d)       !watts/meter**3
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c                
c---beam fusion ion heating profile
c      
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qbfusi_d(j),j=1,nj_d)       !watts/meter**3
c         read(niterdb,10)(blank_d(j),j=1,nj_d)
c                 
c---qmag electron heating, watts/m**3
c      
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qmag_d(j),j=1,nj_d)         !watts/meter**3
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c               
c---sawtooth electron heating, watts/m**3
c      
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qsawe_d(j),j=1,nj_d)        !watts/meter**3
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c               
c---sawtooth ion  heating, watts/m**3
c       
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qsawi_d(j),j=1,nj_d)        !watts/meter**3
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c         
c---radiated power density, watts/m**3
c       
            do j=1,nj_d
             qrad_d(j)=u2d(j,14)        !watts/meter**3
            enddo
c
c---(electron) ohmic power density,watts/m**3
c    
         do j=1,nj_d
            qohm_d(j)=u2d(j,13)        !watts/meter**3
         enddo
c--- GMS 5/18/98 total pressure, Kev/m**3
       if(u2d(1,44) .eq. 0.0) then
cgms   total pressure not availible so use thermal pressure
         if(ipptot.gt.0) then
           if(i_proc.eq.0 .and. istep.eq.1) then
            write(6,'(a31)') 'Error: total pressure not found'
            write(6,'(a32)') 'Computing total P w/ fast ions'
           endif
         do j=1,nj_d
            ptot_d(j) = ene_d(j)*te_d(j) +
     .            (ene_d(j) - (zimp_exp-1.D0)*
     .            en_d(j,nprim_d+1)+enbeam_d(j))*ti_d(j)
         enddo
         else
         if(i_proc.eq.0 .and. istep.eq.1) then
          write(6,'(a31)') 'Error: total pressure not found'
          write(6,'(a32)') 'Computing total thermal pressure'
         endif
         do j=1,nj_d
            ptot_d(j) = ene_d(j)*te_d(j) +
     .            (ene_d(j) - (zimp_exp-1.D0)*
     .            en_d(j,nprim_d+1))*ti_d(j)
         enddo
         endif
       else
         do j=1,nj_d
            ptot_d(j)=u2d(j,44)/1.602D-16     !Kev/m**3
         enddo
       endif
c
c... JEK calculate integrated radiation loss
c
       pradsum=0.D0
       do j=2,nj_d
         drm=r_d(j)-r_d(j-1)
         dvoldr_p=2.D0*pi_m*r_d(j)*
     >            2.D0*pi_m*rmajor_d*hcap_d(j)
         dvoldr_m=2.D0*pi_m*r_d(j-1)*
     >            2.D0*pi_m*rmajor_d*hcap_d(j-1)
         ds=2.D0*pi_m*(r_d(j)*hcap_d(j)+r_d(j-1)*
     >      hcap_d(j-1))/2.D0*drm
         pradsum=pradsum+
     >     zhalf*(dvoldr_p*qrad_d(j)+dvoldr_m*qrad_d(j-1))*drm
       enddo
c       if(i_proc.eq.0) write(6,'(a27,F10.3)') 
c     >   'Total integrated Prad [MW]:', pradsum*1.D-6
c
c... JEK calculate integrated ohmic heating power
c
       if (pohm_d.eq.0.) then
         pohsum=0.D0
         do j=2,nj_d
           drm=r_d(j)-r_d(j-1)
           dvoldr_p=2.D0*pi_m*r_d(j)*
     >              2.D0*pi_m*rmajor_d*hcap_d(j)
           dvoldr_m=2.D0*pi_m*r_d(j-1)*
     >              2.D0*pi_m*rmajor_d*hcap_d(j-1)
           pohsum=pohsum+
     >      zhalf*(dvoldr_p*qohm_d(j)+dvoldr_m*qohm_d(j-1))*drm
         enddo
c      if(i_proc.eq.0) write(6,'(a27,F10.3)') 
c    >   'Total integrated Pohm [MW]:', pohsum*1.D-6
c      else
c      if(i_proc.eq.0) write(6,'(a30,F7.3)') 
c    >   'Total Pohm [MW] from 1D ufile:',pohm_d*1.D-6
       endif
c
c ------------------------------------------------------------------
c...JAK 7/20/95 section placed here: calculate q-profile if not avbailable
c
c...JAK 960401 always read qa1_exp and q01_exp
c alfj_exp=(qa1_exp/q01_exp-1.) corresponds to ja=0 and shata_exp=2.
c
       qa1_exp=u1d(8)
       q01_exp=q0_exp
       alfj1_exp=(qa1_exp/q01_exp-1.)
      if(q_d(nj_d).eq.0.or.iexp_q.eq.-1) then
c
      do j=2,nj_d
       xxq=rho_xp(j)**2
       q_d(j)=1.D0/(1.D0/(q01_exp*xxq*(alfj1_exp+1.D0))*
     >               (1.D0-(1.D0+1.D-10-xxq)**(alfj1_exp+1)))
      enddo
c
      if(qmin_exp.ne.0.) then
      do j=2,nj_d
       xxq=rho_xp(j)**2
       q_d(j)=(rho_xp(j)-rho_qm_exp)**2/(1.D0-rho_qm_exp)**2
     >   *(qa_exp-qmin_exp)+qmin_exp
      enddo 
      endif
         q_d(1)=q_d(2)
      endif
c
c -JAK start 960401--------------------------------------------------
c
c...       In case nonexistent calculate j(r) from q-profile
c
           if (curden_d(1).EQ.0) then 
c
             if(i_proc.eq.0.and. istep.eq.1) 
     &          write(6,*) 'Calculating j(r) from q-profile'
c
c...         Check if total current is given
c
c...         Normalize to total current; if not available then use
c...         d3D formula to relate Ip (in MegAmps) to q_95, etc., (obtained from
c...         Jim Deboo in Feb. 1995., who says it usually agrees with efit
c...         to within 5-10%):
c...         1.e6 convert totohm_d from MA to A,
c
             if (tocur_d.EQ.0) then
               if(i_proc.eq.0.and. istep.eq.1) write(6,*) 'Total 
     &            current calculated from q(r)'
               tocur_d=(5.D0*r_d(nj_d)**2*btor_d/rmajor_d/q_d(nj_d))
     &         *(1.D0+1.5D0*(r_d(nj_d)/rmajor_d)**2)
     &         *(1.D0+kappa_d**2)/2.D0*1.D6
             endif
c
c...         Assume j(r)=j(0)*(1-rhohat**2)**(alpha)
c...         Note that by approximation (circular geometry)
c...         alpha = qa/qa0-1
c
             alpha = qa1_exp/q01_exp-1.D0
             sum=0.D0
             do j=2,nj_d
               drm=r_d(j)-r_d(j-1)
               dvoldr_p=2.D0*pi_m*r_d(j)*
     >                  2.D0*pi_m*rmajor_d*hcap_d(j)
               dvoldr_m=2.D0*pi_m*r_d(j-1)*
     >                  2.D0*pi_m*rmajor_d*hcap_d(j-1)
               ds=2.D0*pi_m*(r_d(j)*hcap_d(j)+
     >            r_d(j-1)*hcap_d(j-1))/2.D0*drm
               curden_d(j) = (1.D0-(r_d(j)/r_d(nj_d))**2)**alpha
               sum = sum + curden_d(j)*ds
c    >           zhalf*(dvoldr_p*curden_d(j)+dvoldr_m*curden_d(j-1))*drm
c    >           /(2.D0*pi_m*rmajor_d) 
             enddo
             cur0 = dabs(tocur_d)/sum
c
             do j=1,nj_d
               curden_d(j)=cur0*(1.D0-(r_d(j)/r_d(nj_d))**2)**alpha
             enddo
           endif 
c
         if (qohm_d(1).EQ.0) then
c
           if(i_proc.eq.0 .and. istep.eq.1) write(6,*)
     >       'Warning rep_iter: no ohmic power profile; making one up'
c
c...       calculate qohm_d from j(r) and loop voltage
c...       by: qohm = E.j
c...       Also calculate volume-integrated power 
c
           cursum=0.D0
           pohsum=0.D0
           do j=2,nj_d
             drm=r_d(j)-r_d(j-1)
             dvoldr_p=2.D0*pi_m*r_d(j)*
     >                2.D0*pi_m*rmajor_d*hcap_d(j)
             dvoldr_m=2.D0*pi_m*r_d(j-1)*
     >                2.D0*pi_m*rmajor_d*hcap_d(j-1)
             ds=2.D0*pi_m*(r_d(j)*hcap_d(j)+
     >          r_d(j-1)*hcap_d(j-1))/2.D0*drm
             cursum = cursum + curden_d(j)*ds
crew error fixed 7.1.96
crew             qohm_d(j)=ABS(vsurf_d*curden_d(j))/(2.*pi_m*rmajor_d)
             qohm_d(j)=(vsurf_d*curden_d(j))/(2.D0*pi_m*rmajor_d)
             pohsum=pohsum+
     >         zhalf*(dvoldr_p*qohm_d(j)+dvoldr_m*qohm_d(j-1))*drm
           enddo
           qohm_d(1)=(vsurf_d*curden_d(1))/(2.D0*pi_m*rmajor_d) 
crew....vsurf_d or curden_d may not be reliable...
crew     better renorm by pohm_d
crew thus in the end we use only a shape of current profile
crew and assume dE/dr=0 ie steady state current profile
           if(pohm_d.gt.0.) then
            do j=1,nj_d
              qohm_d(j)=pohm_d/pohsum*qohm_d(j)
            enddo
           endif
crew why calc?           pohm_d = pohsum
c        if(i_proc.eq.0) write(6,'(a27,F10.3,a10,F7.3)') 
c    >     'Total calculated Pohm [MW]:', pohsum*1.D-6,
c    >     'Vloop = ',vsurf_d
c        if(i_proc.eq.0) write(6,'(a27,F10.3)') 
c    >     'Total integrated Ip  [MA]:', cursum*1.D-6
c
         end if  !end making up ohmic profile
c ------------------------------------------------
c
c---avg major radius of each flux surface, meters, at elevation of mag. axis

c          read(niterdb,'(a)')stflg
          do j=1,nj_d
             rmajavnpsi_d(j) = u2d(j,15)   !meters
          enddo
c
c---avg minor radius of each flux surface,meters, at elevation of mag. axis
c
          do j=1,nj_d
             rminavnpsi_d(j) = u2d(j,16)   !meters
             if(rminavnpsi_d(1).lt.0) rminavnpsi_d(1)=1.D-6
          enddo  
c--- elongation of each flux surface
c
          do j=1,nj_d
           elongx_d(j)=u2d(j,26)
          enddo
c
         if( elongx_d(nj_d).eq.0.) then
          do j=1,nj_d
            elongx_d(j)= (r_d(j)/rminavnpsi_d(j))**2
          enddo
         endif  
c
c--- triangularity at each flux surface
          do j=1,nj_d
           deltax_d(j)=u2d(j,24)
          enddo
c -----------------------------------------------------------------------
c... Smoothing of profiles
c
      do j=1,nj_d
        rho_d(j)=rho(j-1)
        wt(j)=1.D0
      enddo
c
      eps=1.D-6
      sigmas=0.D0
      if (ismoo_ne.gt.0) then
        if (i_proc.eq.0 .and. istep.eq.1) write(*,91) ismoo_ne
        wt(1)=1.D-3     ! do not smooth central and edge values
        wt(nj_d)=1.D-3
        call curvss(nj_d,rho_d,ene_d*1.D-19,wt,0,ismoo_ne,eps,ys,ysp,
     >              sigmas,td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,ier)
        if (ier.gt.0 .and. i_proc.eq.0 .and. istep.eq.1) 
     >     write(*,*) 'ier = ',ier
        do j=1,nj_d
          if (lprint_smoo.gt.0 .and. istep.eq.1) write(*,53) j, 
     >        rho_d(j), ene_d(j), ys(j)*1.D19
          ene_d(j)=ys(j)*1.D19
        enddo
      endif
c
      if (ismoo_ni.gt.0) then
        if (i_proc.eq.0 .and. istep.eq.1) write(*,92) ismoo_ni
        wt(1)=1.D-3     ! do not smooth central and edge values  
        wt(nj_d)=1.D-3
        call curvss(nj_d,rho_d,en_d(1:nj_d,1)*1.D-19,wt,0,ismoo_ni,
     >       eps,ys,ysp,sigmas,td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,ier)
        if (ier.gt.0 .and. i_proc.eq.0 .and. istep.eq.1) 
     >      write(*,*) 'ier = ',ier
        do j=1,nj_d
          if (lprint_smoo.gt.0 .and. istep.eq.1) write(*,53) j, 
     >        rho_d(j), en_d(j,1), ys(j)*1.D19
          en_d(j,1)=ys(j)*1.D19
        enddo
      endif
c
      if (ismoo_zeff.gt.0) then
        if (i_proc.eq.0 .and. istep.eq.1) write(*,93) ismoo_zeff
        wt(1)=1.D-3     ! do not smooth central and edge values  
        wt(nj_d)=1.D-3
        call curvss(nj_d,rho_d,zeff_d,wt,0,ismoo_zeff,
     >       eps,ys,ysp,sigmas,td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,ier)
        if (ier.gt.0 .and. i_proc.eq.0 .and. istep.eq.1 ) 
     >      write(*,*) 'ier = ',ier
        do j=1,nj_d
          if (lprint_smoo.gt.0 .and. istep.eq.1) write(*,53) j,
     >        rho_d(j), zeff_d(j), ys(j)
          zeff_d(j)=ys(j)
        enddo
      endif
c
      if (ismoo_q.gt.0) then
        if (i_proc.eq.0 .and. istep.eq.1) write(*,94) ismoo_q
        wt(1)=1.D-3     ! do not smooth central and edge values  
        wt(nj_d)=1.D-3
        call curvss(nj_d,rho_d,q_d,wt,0,ismoo_q,
     >       eps,ys,ysp,sigmas,td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,ier)
        if (ier.gt.0 .and. i_proc.eq.0 .and. istep.eq.1) 
     >      write(*,*) 'ier = ',ier
        do j=1,nj_d
          if (lprint_smoo.gt.0 .and. istep.eq.1) write(*,53) j,
     >        rho_d(j), q_d(j), ys(j)
          q_d(j)=ys(j)
        enddo
      endif
c
      if (ismoo_vrot.gt.0) then
        if (i_proc.eq.0 .and. istep.eq.1) write(*,95) ismoo_vrot
        wt(1)=1.D-3     ! do not smooth central and edge values  
        wt(nj_d)=1.D-3
        call curvss(nj_d,rho_d,angrot_d*1.D-4,wt,0,ismoo_q,
     >       eps,ys,ysp,sigmas,td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,ier)
        if (ier.gt.0 .and. i_proc.eq.0 .and. istep.eq.1) 
     >      write(*,*) 'ier = ',ier
        do j=1,nj_d
          if (lprint_smoo.gt.0 .and. istep.eq.1) write(*,53) j,
     >        rho_d(j), angrot_d(j), ys(j)*1.D4
          angrot_d(j)=ys(j)*1.D4
        enddo
      endif
c
c -JAK end 960401--------------------------------------------------
c
c...    JAK 7/19/95 Neoclassical chii
c
        ng_nc=1                    !number of hyd. species
        numzones_nc=nj_d           !number of zones
        btf_nc=DABS(btor_d)         !tf (tesla) fld at cntr of outer flux
        drshaf_nc=0.D0              !shaf shift of outer zone bdr. (cm)
c        rminor_nc=r_d(nj_d)*1.e2   !plasma minor radius (cm)
c        rmajor_nc=rmajor_d*1.e2    !major radius (cntr of outer flux) (cm)
        rminor_nc=rminavnpsi_d(nj_d)*1.D2
        rmajor_nc=rmajavnpsi_d(nj_d)*1.D2
c
        do j=1,ng_nc
         aplasm_nc(j)=apgasa        !array of atomic masses of hyd. species
        enddo
c
        do j=1,nj_d
          rhoel_nc(j)=dabs(ene_d(j)*1.D-6)   !electron density (cm**-3)
          rhob_nc(j,1)=dabs(en_d(j,1)*1.D-6)   !hyd. spec. den s (cm**-3)
          rhi_nc(j)=en_d(j,2)*1.D-6    !z.c. array of av. impurity s density  
          rhoi_nc(j)=rhi_nc(j)+rhob_nc(j,1) !z.c. array of total ion density   
          te_nc(j)=dabs(te_d(j)*1.D3) !z.c. array of Te (ev)
          ti_nc(j)=dabs(ti_d(j)*1.D3) !z.c. array of Ti (ev)
          zeff_nc(j)=zeff_d(j)       !z.c. array of plasma zeff
          q_nc(j)=q_d(j)             !z.c. array of safety factor
          aimp_nc(j)=pimpa           !mass of impurity
          xzimp_nc(j)=pimpz          !charge of impurity
        enddo
c
c...    * means: default value given in input.mlt/mlt0in
c...    Note that the comments might be screwed up since some commentlines
c...    had to be deleted otherwise BASIS could not handle it (...) 
c...    See file neoclass.basf for the unscrewed comments
c
        call kapisn(
     >          nkimod_nc,         !*kapai model nr desired
     >          aimp_nc,           !*atomic mass of av. impurity
     >          xzimp_nc,          !*atomic number of av. impurity
     >          aplasm_nc,         !array of atomic masses of hyd. species
     >          ng_nc,             !number of hyd. species
     >          rhoel_nc,          !zone centered electron density (cm**-3)
     >          rhob_nc,           !z.c. array of hyd. spec. den s (cm**-3)
     >          rhi_nc,            !z.c. array of av. impurity s density  
     >          rhoi_nc,           !z.c. array of total ion density       
     >          te_nc,             !z.c. array of Te (ev)
     >          ti_nc,             !z.c. array of Ti (ev)
     >          zeff_nc,           !z.c. array of plasma zeff
     >          numzones_nc,       !number of zones 
     >          q_nc,              !z.c. array of safety factor
     >          btf_nc,            !tf (tesla) fld at cntr of outer flux 
     >          drshaf_nc,         !shaf shift of outer zone bdr. (cm)
     >          rminor_nc,         !plasma minor radius (cm)
     >          rmajor_nc,         !major radius (cntr of outer flux) (cm)
     >          istringer_nc,      !Stringer correctioon
     >          xkapi_nc,          !o Neo-Class Ion thermal diff. (cm**2/sec)
     >          xnstari_nc,        !o nu-star-ions, 
     >          xnstare_nc,        !o nu-star-elecs, 
     >          zfluxlim_nc)       !o flux lim flow max temp grad length
c
        do j=1,nj_d
          xkapi_nc(j)=xkapi_nc(j)*1.D-4*en_d(j,1)/elongx_d(j) ! 1/(m*s)
c error 8/1/96
        enddo
c
c...   calculate conductivity
c
       if ((xkineo_d(nj_d).EQ.0).OR.(imyneoclass.eq.1)) then
crew         write(6,*) 'Neoclassical xkineo calculated by kapisn'
         do j=1,nj_d
           xkineo_d(j)=xkapi_nc(j)
         enddo
       endif
c      
c---volume of each flux surface,  meters**3
c                
          do j=1,nj_d
           psivolp_d(j)=u2d(j,21)
          enddo
c                         
c---triangularity of each flux surface
c
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(triangnpsi_d(j),j=1,nj_d)
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c         
c---indentation of each flux surface
c          
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(pindentnpsi(j),j=1,nj_d)
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c           
c--- flux surface area,  meters**2 4.*pi*pi*R0*hcap*rho*<abs(grad rho)>
c         
           sfix=1.
           if(istk.eq.88) sfix=22.825
            do j=1,nj_d
             sfareanpsi_d(j)=sfix*u2d(j,19)
            enddo
c
c---cross sectional area of each flux surface, meters**2
c          
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(cxareanpsi(j),j=1,nj_d)  !meters**2
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c              
c---flux surface average absolute grad rho
c           
c
c FIX THIS CRAP --- o.k. crap fixed JAK
c            
         rfix=1.
         if(istk.eq.88) rfix=4.46
         do j=1,nj_d
	    grho1npsi_d(j)=u2d(j,5)*r_d(nj_d)/rfix
         enddo
c
c---flux surface average ( grad rho)**2
c           
         do j=1,nj_d
	    grho2npsi_d(j)=u2d(j,6)*r_d(nj_d)**2/rfix**2
        enddo        
c
c---plasma boundary:
c            nplasbdry : no. pts on plasma bdry
c            r pts for plasma boundary, meters
c            z pts for plasma boundary, meters
c           
c                 read(niterdb,'(a)')stflg
c                 read(niterdb,9)nplasbdry_d
c                
c                 read(niterdb,'(a)')stflg
cx                 read(niterdb,10)(rplasbdry_d(j), j=1,nplasbdry_d)
c                 read(niterdb,10)(bblank_d(j),j=1,nplasbdry_d)
c                
c                 read(niterdb,'(a)')stflg
cx                 read(niterdb,10)(zplasbdry_d(j), j=1,nplasbdry_d)
c                 read(niterdb,10)(bblank_d(j),j=1,nplasbdry_d)           
c
c       close(unit=niterdb,status='keep')       
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c      
cmnt   It is assumed that 0:jmaxm greater than or equal 1:nj_dd
cmnt   The volume source data is integrated to power and flows
cmnt   then all the data may be splined to a more dense jmaxm grid
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
       rhostar=bscale**(-z2thrd)*ascale**(-5.D0/6.D0)
c
       if(iscale.eq.1) pscale=ascale**(-z3qtr)*rhostar**(-1.5D0)
       if(iscale.eq.2) pscale=ascale**(-z3qtr)*rhostar**(-2.5D0)
       if(iscale.eq.3) pscale=ascale**(-z3qtr)*rhostar**(-3.D0)
       if(iscale.eq.4) pscale=ascale**(-z3qtr)*rhostar**(-3.5D0)
c
      if(pscale.eq.0) pscale=1.D0
c
c    global transfers                 
c
      arho_exp = r_d(nj_d)*ascale
      rmajor_exp = rmajor_d*ascale                             
      bt_exp = dabs(btor_d)*bscale                                 
      elonga_exp=elongx_d(nj_d)
      deltaa_exp=deltax_d(nj_d) 
c
c    direct profile transfers
c
      do j=1,nj_d
       te_exp(j-1)=te_d(j)*bscale**(z2thrd)*ascale**(zthird)
       ti_exp(j-1)=ti_d(j)*bscale**(z2thrd)*ascale**(zthird)
       ne_exp(j-1)=1.D-19*ene_d(j)*bscale**(z4thrd)/ascale**(zthird)
       ni_exp(j-1)=1.D-19*en_d(j,1)*bscale**(z4thrd)/ascale**(zthird)
       nz_exp(j-1)=1.D-19*en_d(j,2)*
     >      bscale**(z4thrd)/ascale**(zthird)
       nfast_exp(j-1)=1.D-19*enbeam_d(j)*bscale**(z4thrd)/
     >      ascale**(zthird)
       nitot_exp(j-1)=ni_exp(j-1)+nz_exp(j-1)
       fi_m(j-1)=ni_exp(j-1)/ne_exp(j-1)
       ptot_exp(j-1)=1.D-19*ptot_d(j)*bscale**2.D0
       te_d(j)=dabs(te_d(j))
       ti_d(j)=dabs(ti_d(j))
       ene_d(j)=dabs(ene_d(j))
       en_d(j,1)=dabs(en_d(j,1))
       en_d(j,2)=dabs(en_d(j,2))
       te_exp(j-1)=dabs(te_exp(j-1))
       ti_exp(j-1)=dabs(ti_exp(j-1))
       ne_exp(j-1)=dabs(ne_exp(j-1))
       ni_exp(j-1)=dabs(ni_exp(j-1))
       nz_exp(j-1)=dabs(nz_exp(j-1))
       nfast_exp(j-1)=dabs(nfast_exp(j-1))
       nitot_exp(j-1)=dabs(nitot_exp(j-1))
       ptot_exp(j-1)=dabs(ptot_exp(j-1))
c   we are assuming dln(ni)/dr=dln(ne)/dr
c   and will use ion plasma sources as electron
c   plasma sources. zeff enters only into collisionality.
c   ni_m is optained from ne_m assuming constant ni_exp/ne_exp
       zeff_exp(j-1)=zeff_d(j)
       q_exp(j-1)=q_d(j)
       rmin_exp(j-1)=rminavnpsi_d(j)*ascale
       rmaj_exp(j-1)=rmajavnpsi_d(j)*ascale
       angrot_exp(j-1)=angrot_d(j)*
     >     (bscale**(z2thrd)*ascale**(zthird))**(zhalf)/ascale    
       h_exp(j-1)=hcap_d(j)
       gradrhosq_exp(j-1)=grho2npsi_d(j)
       gradrho_exp(j-1)=grho1npsi_d(j)
       elong_exp(j-1)=elongx_d(j)
       delta_exp(j-1)=deltax_d(j)
       chiineo_exp(j-1)=xkineo_d(j)/ni_exp(j-1)/1.D19/
     >    (bscale*ascale**(zhalf))
     > /gradrhosq_exp(j-1)
       zptineomax(j-1)=arho_exp*100.D0/zfluxlim_nc(j)
       curden_exp(j-1) = curden_d(j)*bscale*ascale/ascale**2
c
      enddo
c 
        do j=1,jmaxm
         shat_exp(j)=(rho(j)+rho(j-1))/(q_exp(j)+q_exp(j-1))*
     >       (q_exp(j)-q_exp(j-1))/(rho(j)-rho(j-1))
c
c... JAK 7/06/95 zero check added
c
         if (DABS(shat_exp(j)).LT.1.D-6) shat_exp(j)=1.D-6
        enddo
        shat_exp(0)=0.D0
        rmin_exp(0)=1.D-6
        jstinv=1
c
        do j=2,jmaxm
         if (q_exp(j).gt.1..and.q_exp(j-1).le.1.) jstinv=j
        enddo
c
        jstmix=sqrt(2.D0)*jstinv
c
c  surface factor
c
       do j=0,jmaxm
        sfactor(j)=2.D0*pi_m*arho_exp*rho(j)*
     >             h_exp(j)*2.D0*pi_m*rmajor_exp
       enddo
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c    source integration
c
cmnt the neutral source is broaken into wall and volume
cmnt volume is from beams and recombination where as wall is from wall
cmnt recycled neutrals. The latter are estimated. To chamge the estimate
cmnt the parameter wallneut can be changed from 1.0. The ion and cx 
cmnt energy sources proportional to the neutral density will be corrected
cmnt JAK 960320 fusion power added
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
       flow_recom_exp(0)=0.D0
       flow_beam_exp(0)=0.D0
       flow_sdot_exp(0)=0.D0
       flow_exch_exp(0)=0.D0
       flow_exp(0)=0.D0
c       
       vol_exp(0)=0.D0
c       
       do j=2,nj_d 
        drm=r_d(j)-r_d(j-1)
        dvoldr_p=2.D0*pi_m*r_d(j)*2.D0*pi_m*rmajor_d*hcap_d(j)
        dvoldr_m=2.D0*pi_m*r_d(j-1)*2.D0*pi_m*rmajor_d*hcap_d(j-1)      
        vol_exp(j-1)=vol_exp(j-2)+
     >       zhalf*(dvoldr_p+dvoldr_m)*drm
c     
        powe_beam_exp(j-1)=powe_beam_exp(j-2)+
     >   1.D-6*zhalf*(dvoldr_p*qbeame_d(j)+dvoldr_m*qbeame_d(j-1))*drm
        powe_rf_exp(j-1)=powe_rf_exp(j-2)+
     >   1.D-6*zhalf*(dvoldr_p*qrfe_d(j)+dvoldr_m*qrfe_d(j-1))*drm
        powe_oh_exp(j-1)=powe_oh_exp(j-2)+
     >   1.D-6*zhalf*(dvoldr_p*qohm_d(j)+dvoldr_m*qohm_d(j-1))*drm
        powe_rad_exp(j-1)=powe_rad_exp(j-2)+
     >   1.D-6*zhalf*(dvoldr_p*qrad_d(j)+dvoldr_m*qrad_d(j-1))*drm
        powe_ion_exp(j-1)=powe_ion_exp(j-2)+wallneut*
     >   1.D-6*zhalf*(dvoldr_p*qione_d(j)+ dvoldr_m*qione_d(j-1))*drm
        powe_wdot_exp(j-1)=powe_wdot_exp(j-2)+
     >   1.D-6*zhalf*(dvoldr_p*dpedtc_d(j)+dvoldr_m*dpedtc_d(j-1))*drm
        powe_fus_exp(j-1)=powe_fus_exp(j-2)+
     >   1.D-6*zhalf*(dvoldr_p*qfuse_d(j)+dvoldr_m*qfuse_d(j-1))*drm
c    
c...JAK 6/20/95 typing error fixed in powe_rf_exp: should be powi_rf_exp
        powi_beam_exp(j-1)=powi_beam_exp(j-2)+
     >   1.D-6*zhalf*(dvoldr_p*qbeami_d(j)+dvoldr_m*qbeami_d(j-1))*drm
c       powe_rf_exp(j-1)=powe_rf_exp(j-2)+
        powi_rf_exp(j-1)=powi_rf_exp(j-2)+
     >   1.D-6*zhalf*(dvoldr_p*qrfi_d(j)+dvoldr_m*qrfi_d(j-1))*drm
        powi_ion_exp(j-1)=powi_ion_exp(j-2)+wallneut*
     >   1.D-6*zhalf*(dvoldr_p*qioni_d(j)+dvoldr_m*qioni_d(j-1))*drm
        powi_cx_exp(j-1)=powi_cx_exp(j-2)+wallneut*
     >   1.D-6*zhalf*(dvoldr_p*qcx_d(j)+dvoldr_m*qcx_d(j-1))*drm
        powi_wdot_exp(j-1)=powi_wdot_exp(j-2)+
     >   1.D-6*zhalf*(dvoldr_p*dpidtc_d(j)+dvoldr_m*dpidtc_d(j-1))*drm
        powi_fus_exp(j-1)=powi_fus_exp(j-2)+
     >   1.D-6*zhalf*(dvoldr_p*qfusi_d(j)+dvoldr_m*qfusi_d(j-1))*drm
c    
c note sign change to keep pow_ei_exp positive for te>ti
c (JAK 6/21/95; as in reading iterdb files)
c    
        pow_ei_exp(j-1)=pow_ei_exp(j-2)-
     >   1.D-6*zhalf*(dvoldr_p*qdelt_d(j)+dvoldr_m*qdelt_d(j-1))*drm
c
c... JAK 6/21/95 Sign changed for qrad, qione, qioni
c    
        powe_exp(j-1)=
     >    powe_beam_exp(j-1)+powe_rf_exp(j-1)+xoh_exp*powe_oh_exp(j-1)
     >   -xrad_exp*powe_rad_exp(j-1)-powe_ion_exp(j-1)-
     >    xwdot*powe_wdot_exp(j-1)
     >   -pow_ei_exp(j-1)
     >   +xfus_exp*powe_fus_exp(j-1)
c  
        powi_exp(j-1)=
     >    powi_beam_exp(j-1)+powi_rf_exp(j-1)
     >   -powi_ion_exp(j-1)+powi_cx_exp(j-1)-
     >    xwdot*powi_wdot_exp(j-1)
     >   +pow_ei_exp(j-1)
     >   +xfus_exp*powi_fus_exp(j-1)
c    
        flow_wall_exp(j-1)=flow_wall_exp(j-2)+wallneut*
     >       kevdsecpmw*zhalf*(dvoldr_p*(sion_d(j,1))+
     >       dvoldr_m*(sion_d(j-1,1)))*drm
        flow_recom_exp(j-1)=flow_recom_exp(j-2)+
     >       kevdsecpmw*0.5D0*(dvoldr_p*(srecom_d(j,1))+
     >       dvoldr_m*(srecom_d(j-1,1)))*drm 
        flow_beam_exp(j-1)=flow_beam_exp(j-2)+kevdsecpmw*0.5D0*
     >       (dvoldr_p*(snbscale*sbeam_d(j)+sbcx_d(j,1))+
     >       dvoldr_m*(snbscale*sbeam_d(j-1)+sbcx_d(j-1,1)))*drm
c
        flow_sdot_exp(j-1)=flow_sdot_exp(j-2)+
     >       kevdsecpmw*zhalf*(dvoldr_p*dudtsv_d(j,1)+
     >       dvoldr_m*dudtsv_d(j-1,1))*drm
c 
        flow_exp(j-1)=flow_wall_exp(j-1)+flow_recom_exp(j-1)+
     >                nfscale*flow_beam_exp(j-1)-
     >                xsdot*flow_sdot_exp(j-1)
c     
       enddo      
c
       do j=2,nj_d
        vol_exp(j-1)=vol_exp(j-1)*ascale**3.D0
        powi_exp(j-1)=powi_exp(j-1)*pscale
        powe_exp(j-1)=powe_exp(j-1)*pscale
        pow_ei_exp(j-1)=pow_ei_exp(j-1)*pscale
        powi_beam_exp(j-1)=powi_beam_exp(j-1)*pscale
c note if iexp_exch=-1 pow_ei_exp=0 and pow_ei_cor is the actual
c corrected exchange
c
        flow_exp(j-1)=flow_exp(j-1)*
     >     pscale/(bscale**(z2thrd)*ascale**(zthird))
       enddo
c
c      if (i_proc.eq.0) then
c        write(6,'(a27,2F10.3)') 'Total integrated Pnbe,Pnbi [MW]:',
c    >    powe_beam_exp(nj_d-1),powi_beam_exp(nj_d-1)
c        write(6,'(a27,2F10.3)') 'Total integrated Prfe,Prfi [MW]:',
c    >    powe_rf_exp(nj_d-1),powi_rf_exp(nj_d-1)
c      endif
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
     >  flow_exp(j)*ti_exp(j)*(dlog(dabs(ti_exp(j))*dabs(ni_exp(j)))
     >      -dlog(dabs(ti_exp(j-1))*dabs(ni_exp(j-1))))
       enddo
c
c..Diagnostic printout of power flows
c
      if (lprint_pflow .eq. 1) then
      open (40,file='pflow.out',status='unknown')
      write(40,50) nj_d-1
      do j=1,nj_d-1
        write(40,100) j, rho(j), powi_beam_exp(j), -powi_ion_exp(j),
     &                powi_cx_exp(j), pow_ei_exp(j), 
     &                -powi_wdot_exp(j), powi_exp(j)
      enddo
      write(40,55) nj_d-1
      do j=1,nj_d-1
        write(40,100) j, rho(j), powe_beam_exp(j), powe_oh_exp(j),
     &                -powe_ion_exp(j), -powe_rad_exp(j),  
     &                pow_ei_exp(j), -powe_wdot_exp(j), powe_exp(j)
      enddo
      close(40)
  50  format(i2,2x,'rho',5x,'powibeam',5x,'powiion',6x,'powicx',7x,
     &       'powei',8x,'powiwdot',5x,'powitot')
  55  format(i2,2x,'rho',5x,'powebeam',5x,'poweohm',6x,'poweion',6x,
     &       'powerad',6x,'powei',8x,'powewdot',5x,'powetot')
 100  format(i2,2x,0p1f4.2,1p8e13.5)
      endif
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
     >  *(vol_exp(j)/(2.D0*pi_m**2*rmajor_exp))**zhalf
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
        geoalpha(j)=1.D0/rmajor_exp*(vol_exp(j)-vol_exp(j-1))
     >   /(rho(j)-rho(j-1))/arho_exp/(2.D0*pi_m*rho(j)*arho_exp)**2
     >  *(vol_exp(j)/(2.D0*pi_m**2*rmajor_exp))**zhalf
        bteff_exp(j)=bt_exp*rho(j)*arho_exp/rmin_exp(j)*drhodr(j)
       enddo
       georotrate(0)=georotrate(1)
       geoalpha(0)=geoalpha(1)
       geofac(0)=geofac(1)
       drhodr(0)=drhodr(1)
       drhodrrrho(0)=drhodrrrho(1)
c 
c smoothing delta
c
c      call lspolysmooth(3,jmaxm,jmaxm,rho,delta_exp)
c
      if(iexb.eq.1) then
      open(unit=5,status='unknown',access='sequential',
     >             file='omega.dat')
       do j=0,jmaxm
         csda_exp(j)=9.79D5*(te_exp(j)*1.D3)**zhalf/
     >               (arho_exp*100.D0)/(amassgas_exp)**zhalf
         read(5,*) dummy1,aomega
         megamma_exp(j)=1.D3*aomega/csda_exp(j)
       enddo
       close(5)
       endif
c
      call derivedexpprofiles(amin_d,tocur_d)
c
      call ptorupdate      
c
      deallocate(rhob_nc,zfluxlim_nc, xnstari_nc, xnstare_nc, 
     &           xkapi_nc)
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
 53   format(i2,2x,0p1f4.2,0p6e13.5)
 57   format('   ierr = ',i3,',ufile=',a50)
 60   format(' Reading time-dependent 2D data ...')
 65   format(' Reading time-dependent 1D data ...')
 91   format(' Smoothing electron density profile ...',0p1f4.2)
 92   format(' Smoothing ion density profile ...',0p1f4.2)
 93   format(' Smoothing Zeff profile ...',0p1f4.2)
 94   format(' Smoothing q-profile ...',0p1f4.2)
 95   format(' Smoothing angrot profile ...',0p1f4.2)
 500  format(6x,' *** Warning: simulation time has exceeded ',
     >       'max data time ***')
c
 999  return
      end
