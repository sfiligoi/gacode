c reads full 2d time series from u-file.c@readiterdb.f
c jek 19-Jan-11 version 2.0
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c... Reads in data from Onetwo ascii iterdb file
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c  
      subroutine readiterdb(tok,shot,cudir,iptot,itorque,
     >           ncl_flag)
c
      use data_interface
c
      implicit none
cc      include 'input.m'
cc      include 'data_d.m'
cc      include 'data_exp.m'
cc      include 'glf.m'
c
      integer iptot,itorque,ncl_flag
      character(50) cudir
      character(40) shot
      character(6) phase
      character(10) tok
      integer iflag
      character*15 extension
      character*50 iterdbfile
      character cdfile*130
c
      integer niterdb, i, j, jj, jn
      integer jmaxm
      integer idchar, ifchar, ichar, ichmax
      character*132 headerline
      character*2 stflg
c
      real*8 pi_m, pgasa, pgasz, bgasa, bgasz
     & , rhostar, drm, dvoldr_p, dvoldr_m, dummy1, aomega
      real*8 kevdsecpmw
c
      parameter(niterdb=11)
cc      namelist/datafiles/nstk,xp_time,cudir,tok,shot,phase,ismooth
c      real, allocatable, dimension(:) :: blank_d, bblank_d
c      allocate (blank_d(nj), bblank_d(1000))
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      pi_m=3.1415926D0
c      xi=(0.D0,1.D0)
      kevdsecpmw=1.6022D-19*1.D3*1.D-6
c      jmaxm=mxgrid
c
c     Assume the worki and beam ions are same and the impurity
c     is carbon for the moment.
c
cc      pgasa=amassgas_exp
      pgasa = 2.0
      pgasz=1
      bgasa=pgasa
      bgasz=1
c     pimpa=12
c     pimpz=6
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... read ONETWO iterdb file from external directory
c    JEK 3/19/07
c
      extension = shot
      iterdbfile = 'iterdb.'//extension
c      write(*,*) 'iterdbfile = ',iterdbfile
c
c... count characters in directory name -> idchar
c
      ichmax = 60
      idchar = 0
      do j=1,ichmax
        if ( cudir(j:j)  .ne. ' ' ) then
          idchar = idchar + 1
        else
          go to 5
        endif
      enddo
 5    continue
c
c... count characters in iterdb filename -> ifchar
c
      ifchar = 0
      do j=1,ichmax
        if ( iterdbfile(j:j)  .ne. ' ' ) then
          ifchar = ifchar + 1
        else
          go to 6
        endif
      enddo
  6   continue
c
c       write (*,*) 'idchar = ',idchar
c       write (*,*) 'ifchar = ',ifchar
       if (    cudir(idchar:idchar) .ne. '/') then
          ichar = idchar + 1 + ifchar
          cdfile = cudir(1:idchar) // '/' // iterdbfile(1:ifchar)
       else
          ichar  = idchar + ifchar
          cdfile = cudir(1:idchar) // iterdbfile(1:ifchar)
       endif
       write (*,*) ' iterdb = ',cdfile(1:ichar)
c
       open(unit=niterdb,status='old',access='sequential',
     &      file=cdfile(1:ichar))
c old coding to read in same dir as xptor
c     &      file='iterdb.'//extension)
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
         read(niterdb,'(a)')headerline  

c ishot : shot number
         read(niterdb,'(a)') stflg
         read(niterdb,9) ishot_d

c nj : the size of the vectors printed in this file
         read(niterdb,'(a)') stflg
         read(niterdb,9) nj_d
         if(nj_d.gt.nj+1)then
           write(*,900) nj_d,nj+1
           stop
         else
           jmaxm=nj_d-1
cc           mxgrid=jmaxm
         endif

c nion : the number of ion species
         read(niterdb,'(a)')stflg
         read(niterdb,9)nion_d

c nprim : the number of primary ion species
         read(niterdb,'(a)')stflg
         read(niterdb,9)nprim_d

c nimp : the number of impurity ion species
         read(niterdb,'(a)')stflg
         read(niterdb,9)nimp_d

c nneu : the number of neutral ion species
         read(niterdb,'(a)')stflg
         read(niterdb,9)nneu_d

c ibion : index of beam species
         read(niterdb,'(a)')stflg
         read(niterdb,9)ibion_d

c namep : name(s) of primary ion species
         read(niterdb,'(a)')stflg
         read(niterdb,8)(namep_d(i),i=1,nprim_d)

c namei : name(s) of impurity ion species
         read(niterdb,'(a)')stflg
         read(niterdb,8)(namei_d(i),i=1,nimp_d)

c namen : name(s) of neutral ion species
         read(niterdb,'(a)')stflg
         read(niterdb,8)(namen_d(i),i=1,nneu_d)

c time :  time at which data is printed
         read(niterdb,'(a)')stflg
         read(niterdb,10)time_d

c Rgeom : major radius of geometric
         read(niterdb,'(a)')stflg
         read(niterdb,10)rgeom_d                   !meters

c Rmag :  major radius of mag axis, meters
         read(niterdb,'(a)')stflg
         read(niterdb,10)rmag_d                    !meters

c R0 : major radius of vaccuum btor ref
         read(niterdb,'(a)')stflg
         read(niterdb,10)rmajor_d                  !meters

c kappa : plasma elongation
         read(niterdb,'(a)')stflg
         read(niterdb,10)kappa_d

c delta : plasma triangularity
         read(niterdb,'(a)')stflg
         read(niterdb,10)deltao_d

c pindent : plasma indentation
         read(niterdb,'(a)')stflg
         read(niterdb,10)pindento_d

c volo : plasma volume,meters**3
         read(niterdb,'(a)')stflg
         read(niterdb,10)volo_d                    !meters**3

c cxareao :plasma cross sectional area, meters**2
         read(niterdb,'(a)')stflg
         read(niterdb,10)areao_d                   !meters**2

c Btor : vaccuum toroidal field at rmajor, tesla
         read(niterdb,'(a)')stflg
         read(niterdb,10)btor_d                    !tesla

c total,ohmic,bootstrap,beam,and rf currents, amps
         read(niterdb,'(a)')stflg
         read(niterdb,10)tocur_d,totohm_d,totboot_d,totbeam_d,totrf_d !amps

c betap : poloidal beta
         read(niterdb,'(a)')stflg
         read(niterdb,10)betap_d

c beta : toroidal beta
         read(niterdb,'(a)')stflg
         read(niterdb,10)beta_d

c ali : plasma inductance
         read(niterdb,'(a)')stflg
         read(niterdb,10)ali_d

c te0 : central electron temperature
         read(niterdb,'(a)')stflg
         read(niterdb,10)te0_d                        !kev

c ti0 : central ion temperature
         read(niterdb,'(a)')stflg
         read(niterdb,10)ti0_d                        !kev

c---psi on rho grid,volt*sec/rad

          read(niterdb,'(a)')stflg
          read(niterdb,10)(psir_d(j),j=1,nj_d)      !volt*se/rad
cx          read(niterdb,10)(blank_d(j),j=1,nj_d)

c---rho grid, meters

          read(niterdb,'(a)')stflg
          read(niterdb,10)(r_d(j),j=1,nj_d)          !meters

c---fcap, (ie f(psilim)/f(psi) )

          read(niterdb,'(a)')stflg
          read(niterdb,10)(fcap_d(j),j=1,nj_d)
cx          read(niterdb,10)(blank_d(j),j=1,nj_d)

c---gcap, (ie <(grad rho)**2*(R0/R)**2> )

          read(niterdb,'(a)')stflg
          read(niterdb,10)(gcap_d(j),j=1,nj_d)
cx          read(niterdb,10)(blank_d(j),j=1,nj_d)

c---hcap, (ie (dvolume/drho)/(4*pi*pi*R0*rho))

          read(niterdb,'(a)')stflg
          read(niterdb,10)(hcap_d(j),j=1,nj_d)

c*** read in NCLASS quantities if ncl_flag=1 ***
          if(ncl_flag.eq.1) then
            write(*,25)

c--xb2, <B^2>

          read(niterdb,'(a)')stflg
          read(niterdb,10)(xb2_d(j),j=1,nj_d)

c--xbm2, <1/B^2>

          read(niterdb,'(a)')stflg
          read(niterdb,10)(xbm2_d(j),j=1,nj_d)

c--xngrth, <n.grad theta>

          read(niterdb,'(a)')stflg
          read(niterdb,10)(xngrth_d(j),j=1,nj_d)

c--xgrbm2, <grad-rho^2/B^2>

          read(niterdb,'(a)')stflg
          read(niterdb,10)(xgrbm2_d(j),j=1,nj_d)

c--fm1, 1st poloidal moment

          read(niterdb,'(a)')stflg
          read(niterdb,10)(fm1_d(j),j=1,nj_d)

c--fm2, 2nd poloidal moment

          read(niterdb,'(a)')stflg
          read(niterdb,10)(fm2_d(j),j=1,nj_d)

c--fm3, 3rd poloidal moment

          read(niterdb,'(a)')stflg
          read(niterdb,10)(fm3_d(j),j=1,nj_d)

c--fhat, muo*F/dpsidr

          read(niterdb,'(a)')stflg
          read(niterdb,10)(fhat_d(j),j=1,nj_d)
c
          endif
c
c---electron temperature, kev

          read(niterdb,'(a)')stflg
          read(niterdb,10)(te_d(j),j=1,nj_d)           !kev

c---ion temperature, kev

          read(niterdb,'(a)')stflg
          read(niterdb,10)(ti_d(j),j=1,nj_d)           !kev
c
c---q (ie safety factor) profile

          read(niterdb,'(a)')stflg
          read(niterdb,10)(q_d(j),j=1,nj_d)
c
c---electron density,#/m**3 

          read(niterdb,'(a)')stflg
          read(niterdb,10)(ene_d(j),j=1,nj_d)           !#/meter**3

c---primary ion density,#/m**3,species
c---impurity ion density,#/m**3,species

       do jj=1,nion_d
         do j=1,nj_d
           en_d(j,jj)=0.D0
         enddo
       enddo
c      jp=0
c      ji=0
       do jj=1,nion_d
c          if(jj .le. nprim_d)jp=jp+1
c          if(jj .gt. nprim_d)ji=ji+1
               read(niterdb,'(a)')stflg
               read(niterdb,10)(en_d(j,jj),j=1,nj_d)    !#/meter**3
       enddo

c---sion   : source due to ionization,        #/(m**3*sec),species
c---srecom : source due to recombination,     #/(m**3*sec),species
c---scx    : source due to cx thermal neut.,  #/(m**3*sec),species
c---sbcx   : sink due to cx with beam neut.   #/(m**3*sec),species
c---s      : total source rate,               #/(m**3*sec),species
c---dudt   : s dot,                           #/(m**3*sec),species
 
       do jj=1,nprim_d
               read(niterdb,'(a)')stflg
               read(niterdb,10)(sion_d(j,jj),j=1,nj_d)    !#/meter**3/sec
               read(niterdb,'(a)')stflg
               read(niterdb,10)(srecom_d(j,jj),j=1,nj_d)  !#/meter**3/sec
               read(niterdb,'(a)')stflg
               read(niterdb,10)(scx_d(j,jj),j=1,nj_d)     !#/meter**3/sec
               read(niterdb,'(a)')stflg
               read(niterdb,10)(sbcx_d(j,jj),j=1,nj_d)    !#/meter**3/sec
               read(niterdb,'(a)')stflg
               read(niterdb,10)(s_d(j,jj),j=1,nj_d)       !#/meter**3/sec
               read(niterdb,'(a)')stflg
               read(niterdb,10)(dudtsv_d(j,jj),j=1,nj_d)  !#/meter**3/sec
          enddo

c---fast ion density, #/m**3, species
       
          read(niterdb,'(a)')stflg  
          read(niterdb,10)(enbeam_d(j),j=1,nj_d)          !#/meter**3
       
c---neutral density,  #/m**3,species

       do jn=1,nneu_d
          read(niterdb,'(a)')stflg
          read(niterdb,10)(enn_d(j,jn),j=1,nj_d)     !#/meter**3
       enddo

c---neutral density from wall source, #/m**3, species

       do jn=1,nneu_d
          read(niterdb,'(a)')stflg
          read(niterdb,10)(ennw_d(j,jn),j=1,nj_d)     !#/meter**3
       enddo

c---neutral density from volume source, #/m**3,species (recomb and beam cx)

       do jn=1,nneu_d
          read(niterdb,'(a)')stflg
          read(niterdb,10)(ennv_d(j,jn),j=1,nj_d)     !#/meter**3
       enddo

c---volume source of neutrals, #/(m**3*sec),species

       do jn=1,nneu_d
          read(niterdb,'(a)')stflg
          read(niterdb,10)(volsn_d(j,jn),j=1,nj_d)!#/(m**3*sec)  
       enddo

c---sbion : beam electron source,  #/(m**3*sec
           
          read(niterdb,'(a)')stflg
          read(niterdb,10)(sbion_d(j),j=1,nj_d)    !#/(m**3*sec)

c---sbeam : beam thermal ion source,  #/(m**3*sec)
           
          read(niterdb,'(a)')stflg   
          read(niterdb,10)(sbeam_d(j),j=1,nj_d)      !#/(m**3*sec)

c---total current density, amps/m**2
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(curden_d(j),j=1,nj_d)        !amps/meter**2
      
c---ohmic current density, amps/m**2
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(curohm_d(j),j=1,nj_d)        !amps/meter**2

c---bootstrap current density, amps/m**2
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(curboot_d(j),j=1,nj_d)       !amps/meter**2
       
c--- beam driven current density, amps/m**2
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(curbeam_d(j),j=1,nj_d)      !amps/meter**2
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
                
c---rf current density, amps/m**2
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(currf_d(j),j=1,nj_d)         !amps/meter**2

c---rho*bp0*fcap*gcap*hcap, tesla*meters
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(rbp_d(j),j=1,nj_d)
c          read(niterdb,10)(blank_d(j),j=1,nj_d)         !tesla*meters      

c---zeff profile:
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(zeff_d(j),j=1,nj_d)

c---angular rotation speed profile, rad/sec
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(angrot_d(j),j=1,nj_d)       !rad/sec

c---electron thermal diffusivity, meters*2*/sec on half grid
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(chieinv_d(j),j=1,nj_d)       !meters**2/sec
          
c---ion thermal diffusivity,meters*2*/sec  on half grid 
      
          read(niterdb,'(a)')stflg 
          read(niterdb,10)(chiinv_d(j),j=1,nj_d)        !meters**2/sec

c---ion neocl.  thermal conductivity, 1/(m*sec) on half grid
          
          read(niterdb,'(a)')stflg
          read(niterdb,10)(xkineo_d(j),j=1,nj_d)          !1/(m*sec)

c---wdot,electrons, watts/m**3 d(electron energy)/dt profile:
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(dpedtc_d(j),j=1,nj_d)        !watts/meter**3
     
c---wdot,ions,watts/m**3 d(ion energy)/dt profile:
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(dpidtc_d(j),j=1,nj_d)        !watts/meter**3

c---electron conduction, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qconde_d(j),j=1,nj_d)        !watts/meter**3
      
c---ion conduction, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qcondi_d(j),j=1,nj_d)        !watts/meter**3

c---electron convection,watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qconve_d(j),j=1,nj_d)        !watts/meter**3

c---ion convection, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qconvi_d(j),j=1,nj_d)        !watts/meter**3

c---power to elec.from beam, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qbeame_d(j),j=1,nj_d)        !watts/meter**3

c---qdelt : electron-ion equilibration, watts/m**3
      
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qdelt_d(j),j=1,nj_d)         !watts/meter**3
c
c---power to ions from beam, watts/m**3:
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qbeami_d(j),j=1,nj_d)        !watts/meter**3
     
c---qrfe, rf electron heating, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qrfe_d(j),j=1,nj_d)          !watts/meter**3
      
c---qrfi, rf ion heating, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qrfi_d(j),j=1,nj_d)          !watts/meter**3

c---qione, recombination and impact ionization, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qione_d(j),j=1,nj_d)         !watts/meter**3
       
c---qioni, recombination and impact ionization, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qioni_d(j),j=1,nj_d)         !watts/meter**3
       
c---qcx,  neutral-ion charge exchange, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qcx_d(j),j=1,nj_d)           !watts/meter**3

c---2d MHD equil. electron heating, watts/m**3
       
          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qe2d_d(j),j=1,nj_d)        !watts/meter**3
          read(niterdb,10)(blank_d(j),j=1,nj_d)
          
c---2d MHD equil.ion heating, watts/m**3
       
          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qi2d_d(j),j=1,nj_d)         !watts/meter**3
          read(niterdb,10)(blank_d(j),j=1,nj_d)
          
c---fusion electron heating, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qfuse_d(j),j=1,nj_d)        !watts/meter**3
                
c---fusion ion heating, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qfusi_d(j),j=1,nj_d)        !watts/meter**3

c---beam fusion electron heating,watts/m**3
    
          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qbfuse_d(j),j=1,nj_d)       !watts/meter**3
          read(niterdb,10)(blank_d(j),j=1,nj_d)
                 
c---beam fusion ion heating profile
       
          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qbfusi_d(j),j=1,nj_d)       !watts/meter**3
          read(niterdb,10)(blank_d(j),j=1,nj_d)
                 
c---qmag electron heating, watts/m**3
       
          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qmag_d(j),j=1,nj_d)         !watts/meter**3
          read(niterdb,10)(blank_d(j),j=1,nj_d)
                
c---sawtooth electron heating, watts/m**3
       
          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qsawe_d(j),j=1,nj_d)        !watts/meter**3
          read(niterdb,10)(blank_d(j),j=1,nj_d)
                
c---sawtooth ion  heating, watts/m**3
       
          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qsawi_d(j),j=1,nj_d)        !watts/meter**3
          read(niterdb,10)(blank_d(j),j=1,nj_d)

c---radiated power density, watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qrad_d(j),j=1,nj_d)           !watts/meter**3

c---(electron) ohmic power density,watts/m**3
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(qohm_d(j),j=1,nj_d)           !watts/meter**3

c---avg major radius of each flux surface, meters, at elevation of mag. axis
       
          read(niterdb,'(a)')stflg
          read(niterdb,10)(rmajavnpsi_d(j),j=1,nj_d)   !meters

c---avg minor radius of each flux surface,meters, at elevation of mag. axis

          read(niterdb,'(a)')stflg
          read(niterdb,10)(rminavnpsi_d(j),j=1,nj_d)   !meters
c make gradient smooth near axis so that drhodr will be smooth
          amin_d=rminavnpsi_d(nj_d)   !meters
          
c---volume of each flux surface,  meters**3
                 
          read(niterdb,'(a)')stflg
          read(niterdb,10)(psivolp_d(j),j=1,nj_d)      !meters**3
c
c--- elongation of each flux surface
c
          read(niterdb,'(a)')stflg
          read(niterdb,10)(elongx_d(j),j=1,nj_d)
c
c---triangularity of each flux surface
c
          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(triangnpsi_d(j),j=1,nj_d)
          read(niterdb,10)(deltax_d(j),j=1,nj_d)
c
c---indentation of each flux surface
c 
          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(pindentnpsi(j),j=1,nj_d)
          read(niterdb,10)(blank_d(j),j=1,nj_d)
c
c--- flux surface area,  meters**2 4.*pi*pi*R0*hcap*rho*<abs(grad rho)>
c
          read(niterdb,'(a)')stflg
          read(niterdb,'(a)')stflg
          read(niterdb,10)(sfareanpsi_d(j),j=1,nj_d)  !meters**2
c
c---cross sectional area of each flux surface, meters**2
c    
          read(niterdb,'(a)')stflg
          read(niterdb,10)(cxareanpsi_d(j),j=1,nj_d)  !meters**2
c      
c---flux surface average absolute grad rho
c     
         read(niterdb,'(a)')stflg
         read(niterdb,10)(grho1npsi_d(j),j=1,nj_d)
c
c---flux surface average ( grad rho)**2
c 
         read(niterdb,'(a)')stflg
         read(niterdb,10)(grho2npsi_d(j),j=1,nj_d)
c
c---plasma boundary:
c            nplasbdry : no. pts on plasma bdry
c            r pts for plasma boundary, meters
c            z pts for plasma boundary, meters
c
             read(niterdb,'(a)')stflg
             read(niterdb,9)nplasbdry_d
c
             read(niterdb,'(a)')stflg
cx           read(niterdb,10)(rplasbdry_d(j), j=1,nplasbdry_d)
             read(niterdb,10)(bblank_d(j),j=1,nplasbdry_d)
c
             read(niterdb,'(a)')stflg
cx           read(niterdb,10)(zplasbdry_d(j), j=1,nplasbdry_d)
             read(niterdb,10)(bblank_d(j),j=1,nplasbdry_d)
c
c---torque density  Nt-m/m^3 (old iterdb files were in dyne-cm/cm^3 = 10 Nt-m/m^3)
c
             if(itorque.ne.0)then
               read(niterdb,'(a)')stflg
               read(niterdb,10)(torque_d(j), j=1,nj_d)
c               write(*,10) (torque_d(j),j=1,nj_d)
c this rescaling is not needed anymore. 
c For old iterdb files in cgs units set s0(4)=0.1, s0(5)=0.1, s1(4)=0.1, s1(5)=0.1 
c to rescale the torque. 
c                if(idata.eq.1) then
c                 write(*,*) 'Rescaling torque density ...'
c                 do j=1,nj_d
c                  torque_d(j) = torque_d(j) / 10.D0  ! convert to N-M/M**3
c                 enddo
c                endif
             endif
c
c---total and fast ion pressure (keV/m^3)
c   Pa = N/m**2, keV/m**3=Pa/1.602e-16
c   If iptot_d=0, then total thermal pressure computed in pressure.f
c   Note: pfast_exp has 1.e19 factored out since alpha_exp,m
c   has 1.e19 factored out of densities
c
       if(iptot.eq.0)then
         write(*,'(a46)')
     >   'computing total pressure from thermal pressure'
         do j=1,nj_d
           ptot_d(j) = 1.6022D-16*ene_d(j)*(te_d(j)+ti_d(j))   ! Pascals
           pfast_d(j) = 0.0
         enddo
       elseif(iptot.eq.1) then
         write(*,'(a32)') 'Reading ptot only' 
         read(niterdb,'(a)')stflg
         read(niterdb,10)(ptot_d(j), j=1,nj_d)
         do j=1,nj_d
           pfast_d(j)=ptot_d(j)
     >     -1.6022D-16*ene_d(j)*(te_d(j)+ti_d(j))  !Pascals
         enddo
       elseif(iptot.eq.2) then
         write(*,'(a32)') 'Reading ptot and pfast'
         read(niterdb,'(a)')stflg
         read(niterdb,10)(pfast_d(j), j=1,nj_d)
         read(niterdb,'(a)')stflg
         read(niterdb,10)(ptot_d(j), j=1,nj_d)
       endif
c
       close(unit=niterdb,status='keep')
c      close(unit=niterdb)
c
      do j=1,nj_d
       bp0_d(j)=rbp_d(j)/( max(r_d(j),1.d-6)*fcap_d(j)*
     >          gcap_d(j)*hcap_d(j) )
c       write(*,50) j,rho_d(j), r_d(j), fcap_d(j),gcap_d(j),
c     >            hcap_d(j), bp0_d(j)
      enddo
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
  8   format(5(2x,a))      !common character variable write/read format
  9   format(5(2x,i6))     ! common integer write/read format
 10   format(5(2x,1pe14.4))! common write/read  format
 12   format(a9,i3)! common write/read  format
 25   format(' Reading geometric variables for NCLASS ...')
 50   format(i2,2x,0p1f4.2,0p6f10.5)
 51   format(i3,2x,0p1f6.4,1p6e12.4)
 52   format(i2,2x,i2,2x,0p1f4.2,0p6f10.5)
 53   format(i2,2x,0p1f4.2,0p6e13.5)
 54   format(i2,2x,0p6e13.5)
 85   format(' Smoothing electron density profile ...',0p1f4.2)
 86   format(' Smoothing ion density profile ...',0p1f4.2)
 87   format(' Smoothing Zeff profile ...',0p1f4.2)
 88   format(' Smoothing q-profile ...',0p1f4.2)
 89   format(' Smoothing angrot profile ...',0p1f4.2)
 90   format(' Smoothing fast ion density profile ...',0p1f4.2)
 91   format(' Smoothing triangularity profile ...',0p1f4.2)
 95   format(' Smoothing Qnb profile ...',0p1f4.2)
 100  format(' Error opening iterdb file')
 900  format(' Warning: nj_d larger than maximum pts, nj_d = ',i3,
     > 'max=',i3)
c
       end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
