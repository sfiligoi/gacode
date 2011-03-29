c@gks_driver.f
c 28-may-02 
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... This is the driver routine for GKS
c
c    Settings input:
c       1) files 'gksin' and 'gkso0in' - provide the default settings
c       2) file 'gksin1' -  switches
c    Data input:
c       1) idata=0 for reading ITER PDB ufiles with rep_iter.f
c       2) idata=1 for reading Onetwo iterdb files with readexpprofiles.f
c       3) idata=2 for creating test profiles with readtestprofiles.f
c       4) ipptot=1 to include fast ions in total pressure if PTOTR not given
c    Shafranov shift stabilization:
c       To include the effect of Shafranov shift stabilization in GKS
c       the variable 'xalpha' is used as follows:
c            xalpha
c              0         no alpha-stabilization
c             >0         use exp profiles for alpha
c             <0         compute alpha from model profiles w/ fast ion pressure
c                        if using ITER ufiles. Set ipptot=1 to compute ptot w/ 
c                        n_fast, t_fast=ti, otherwise compute ptot as thermal pressure
c    Updating sources:
c       The sources from ITER ufiles can be updated as a function of time.
c       1) set time_upd for simulation time when updating sources is to occur
c       2) set upsrc(1)=1 to update q-profile, upsrc(2)=1 to update power flows
c    Grid size:
c       The default grid size is 50 points
c       2) set jmaxm,mxgrid to desired size in file 'in' (defaulted to 50 in mltin)
c       3) set grid no. for BC in 'in' (default = 45 in mltin)
c       Note: Maximum number of time pts = 2500
c    Output:
c       Output switches:
c          lprint_cdf   =1 writes all output to a NetCDF file 'gksout.nc'
c       ASCII Output files: 
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      program gks_driver
c
      use gks_var
      use data_interface
c
c
      implicit none
c
      include 'input.m'
cc      include 'data_d.m'
      include 'data_exp.m'
      include 'data_t.m'
      include 'glf.m'
      include 'gks_out.m'
c
      character line*132
      integer j, jj, k, nsave
      integer njj,njjr,j1,j2,iproc,time_series,mxgrid
      real*8 dx, norm_bar, endtime, btscale
      real cputime, time_cp
      character*12 ZDATE
      character*10 UDATE
      character*30 foo1,foo2
      integer imon, ivals(8)
      character*3 mons(12)
      data mons/'jan','feb','mar','apr','may','jun','jul','aug',
     >          'sep','oct','nov','dec'/
c
c... namelists
c
      namelist /in/ ntheta,nperiod,ipar,ngauss,icv,negrid,
     * nspec,eps,shift,dbeam,shat,pk,epsa,epsl,width0,beta,
     * zeff,teti,zeff5,fprim1,fprim2,fprim3,fprim4,fprim5,
     * tprim1,tprim2,tprim3,tprim4,tprim5,vnewk1,vnewk2,vnewk3,
     * vnewk4,vnewk5,bakdif1,bakdif2,bakdif3,bakdif4,bakdif5,
     * amass2,amass3,amass4,amass5,temp2,temp3,temp4,temp5,
     * z2,z3,z4,z5,uprim1,uprim2,uprim3,uprim4,uprim5
c
      namelist /in0/ nscreen,nout,nstep,ngamstep,isvar,nwrite,iphi,
     * delt,ominst,phydif,power,tol,gamma,absom,fcv,fv,gridfac,fexp1,
     * fexp3,test1,test2,alr,ali,al1r,al1i,alar,alai,al1ar,al1ai,
     * ncspec1,ncspec2,ncspec3,ncspec4,ncspec5,nce1,nce2,nce3,nce4,
     * nce5,aky1,aky2,aky3,aky4,aky5,aky6,aky7,aky8,aky9,aky10,
     * aky11,aky12,aky13,aky14,aky15,aky16,naky,nt0,thetamin,thetamax,
     * nt0min,nt0max,ecut,adamp,nteststep
c
       namelist /in1/ igks_model,idata, bt_flag, tok, shot, delt,
     *  xp_time, cudir,i_solve,i_bpar,igeo,igeot,icontinue,
     *  cdebye,cbetae,xalpha,alpha_mhd_loc,ibranch,
     *  y_bpar,x_bpar,xy_bpar,y_mhd,signomega,
     *  igyro_fix,igyro_e,kys0,tol_f,xkymin_gf,xkymax_gf,iptot,
     *  ivar,jvar,irunmax,jrunmax,xvarmin,xvarmax,zvarmin,zvarmax,
     *  igraph,mprint,irot1,irot2,igeo_m,i_dengrad,
     *  jout_m,jin_m,echconv,igks,nstep,amassgas_exp,lprint_cdf,
     *  z2,z3,aky1,corot,igeo_print, itorque,
     *  cnewk2,cnewk3,ncspec2,nbasis_min,nbasis_max,amass2,amass3,
     *  alpha_mach, alpha_cur,alpha_e,alpha_p,alpha_p_cur,
     *  lprint_pflow, iexb, ipfst, idatzero,
     *  itest_ntcc, ialpha,save_tglf,overwrite_tglf,
     *  ismooth_data,tglf_defaults,gks_defaults
c
      version = "TGLF_1.93"
c
      do j=0,jmaxm
        ky_m(j) = 0.0
        anrate_m(j) = 0.0
        dnrate_m(j) = 0.0
        anfreq_m(j) = 0.0
        dnfreq_m(j) = 0.0
        gamma_ion(j)=0.0
        gamma_electron(j)=0.0
        freq_ion(j)=0.0
        freq_electron(j)=0.0
        te_bar_m(j)=0.0
        ne_bar_m(j)=0.0
        ne_te_phase_ion(j)=0.0
        ne_te_phase_electron(j)=0.0
        zpte_crit_m(j)=0.0
        zpti_crit_m(j)=0.0
        mhd_DR_m(j)=0.0
        peflx_m(j)=0.0
        qeflx_m(j)=0.0
        qiflx_m(j)=0.0
      enddo
c
c set defaults
c
      call set_gks_defaults
      iproc=0
      idatzero=0
      time_series=0
      endtime=1000.

c namelists ....................................................
c
c
      open (unit=11, file='gkslog',status='unknown')
      open (unit=7, file='gksin1', status='old' )
      read (7,nml=in1)
      close (7)
      iptotr=iptot
      write(11,*)"gksin1",shot
c
c
      ZDATE=' '
      call date_and_time(udate,foo1,foo2,ivals)
      read(udate(5:6),'(i2)') imon
      zdate=udate(7:8)//mons(imon)//udate(1:4)
      write(11,125)
      write(11,*) 'igks = ',igks
c
      if (igks.eq.0)then
c   single k run with no data read
         write(*,*) 'reading gksin ...'
         open (unit=7, file='gksin', status='old' )
         read (7,nml=in)
         close (7)
c         write(*,*)"gksin",ntheta
         open (unit=7, file='gks0in', status='old' )
         read (7,nml=in0)
         close (7)
c         write(*,*)"gks0in",nscreen
         if(igeo.eq.1)then
c run the test case for real geometry
            iptot=0
            call geogks_toq
            call geogks_profile(25,25,1)
         else
       an5= zeff5/(z5**2)
       an2= ( zeff -1. -z5*(z5-1.)*an5 )/( (z2-1.)*z2 )
       an1= 1. -an2*z2 -an5*z5 -dbeam
            call gstotal
         endif
           write(*,*)" agammas(1) = ",agammas(1)
           write(*,*)" dgammas(1) = ",dgammas(1)
           write(*,*)" afreqs(1)  = ",afreqs(1)
           write(*,*)" kys(1)     = ",kys(1)
           write(*,*)" peflux     = ",peflxa
           write(*,*)" qeflux     = ",eeflxa
           write(*,*)" qiflux     = ",eiflxa
c           write(*,*)"igeo = ",igeo
c           write(*,*)"cbetae = ",cbetae
c           write(*,*)"cdebye = ",cdebye
c           write(*,*)"xalpha = ",xalpha
c           write(*,*)"ipar = ",ipar
c           write(*,*)"iptot = ",iptot
c           write(*,*)"irot1,2 = ",irot1,irot2
c           write(*,*)"ncspec2 = ",ncspec2
c           write(*,*)"an(1)=",an(1),"an(2)=",an(2),"an(3)=",an(3)
       endif
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      if(igks.gt.0)then
c read in experimental data
        mxgrid=50
        jmaxm=50
        if(gks_defaults.eq.0)then
         open (unit=7, file='gksin', status='old' )
         read (7,nml=in)
         close (7)
c         write(*,*)"gksin",ntheta
         open (unit=7, file='gks0in', status='old' )
         read (7,nml=in0)
         close (7)
        endif
c
        if(jout_m.gt.jmaxm ) then
         write(11,35)
        endif
        if (idata.eq.2) then
          tok = 'test'
          shot= '999'
          phase='X'
         xp_time=99.D0
        endif 
        write(11,50) zdate, tok, shot
        call data_run(idata,shot,tok,cudir,xp_time,endtime,
     >  time_series,itorque,iptot,ncl_flag,mxgrid,ismooth,
     >  idatzero,iproc)        
cc        if (idata .eq. 0) then
cc            write(11,100)
cc          call readufiles(tok,shot,phase,cudir,mxgrid,ismooth,
cc     &             btscale,xp_time,endtime_pt,time_series,
cc     &             iexp_exch,iexp_q,
cc     &             p_glob_exp,gtauth_exp,gtautot_exp,iproc)
cc        elseif (idata .eq. 1) then
cc            write(11,105)
cc          call readiterdb(mxgrid)
cc          if(ismooth_data.ne.0) call datavg
cc        elseif (idata .eq. -1) then
cc            write(11,109)
cc        else
cc          write(11,*) 'Error: idata out of range',idata
cc		  stop
cc        endif
c
        call datmap
c
c fix ni_exp and nz_exp to insure charge balance
c
       do j=0,jmaxm
        nz_exp(j)=ne_exp(j)*(zeff_exp(j)-1.0)/(z2*(z2-1.0))
        ni_exp(j)=ne_exp(j)-z2*nz_exp(j)-nfst_exp(j)
       enddo      
c
        call expprofiles(1,mxgrid,arho_exp)
c
      endif
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c call to GKS gyro-kinetic stability
c
         if (igks.eq.1)then
c         lowk spectrum at fixed rho
           nsave=nstep
           nstep=1
           if(igeo.eq.0)call gks_profile(jin_m,jin_m)
           if(igeo.eq.1)call geogks_profile(jin_m,jin_m,1)
           nstep=nsave
           dx = -2.D0/20.D0
           call g_vs_logk(20,dx,0)
         endif
         if (igks.eq.2)then
c          highk spectrum at fixed rho
           delt=0.2
           nsave=nstep
           nstep=1
           if(igeo.eq.0)call gks_profile(jin_m,jin_m)
           if(igeo.eq.1)call geogks_profile(jin_m,jin_m,1)
           nstep=nsave
           dx = 2.D0/20.D0
           call g_vs_logk(20,dx,0)
         endif
         if (igks.eq.3)then
c    fixed k radial growthrate profile
           if(igeo.eq.0) call gks_profile(jin_m,jout_m)
           if(igeo.eq.1) call geogks_profile(jin_m,jout_m,1)
         endif
         if (igks.eq.4)then
c    maximum itg growthrate radial profile
           kys0=-DABS(kys0)
           xkymax_gf=2.0
           xkymin_gf=0.01
           if(igeo.eq.0) call gks_profile(jin_m,jout_m)
           if(igeo.eq.1) call geogks_profile(jin_m,jout_m,1)
         endif
         if (igks.eq.5)then
c          critical electron temperature gradient profile
c           xkymax_gf=200.D0
c           xkymin_gf=1.0D-4
           if(delt.lt.0.2D0)delt=0.2D0
           call gradte_crit(jin_m,jout_m,1)
          endif
         if (igks.eq.6)then
c          critical ion temperature gradient profile
c           xkymax_gf=200.D0
c           xkymin_gf=1.0D-4
           call gradti_crit(jin_m,jout_m,1)
          endif
         if (igks.eq.7)then
c          intermediate-k spectrum at fixed rho (ky=0.7-2.0)
           delt=0.15
           nsave=nstep
           nstep=1
           if(igeo.eq.0)call gks_profile(jin_m,jin_m)
           if(igeo.eq.1)call geogks_profile(jin_m,jin_m,1)
           nstep=nsave
           dx = 0.10D0
           call g_vs_logk(13,dx,1)
         endif
c
      if(igks.ne.0)then
       if(igks.le.2 .or. igks.eq.7)then
        if(igeo.eq.0)then
         norm_bar = rhosda_exp(jin_m)**2
         do k=0,jmaxm
          gamma_mks(k) = anrate_m(k)*csda_exp(jin_m)
          dgamma_mks(k)= dnrate_m(k)*csda_exp(jin_m)
          freq_mks(k)  = anfreq_m(k)*csda_exp(jin_m)
          ky_mks(k) = ky_m(k)/(rhosda_exp(jin_m)*arho_exp)
          if(igks_model.eq.1)then
            phi_bar_m(k) = phi_bar_m(k)*norm_bar
            ne_bar_m(k) = ne_bar_m(k)*norm_bar
            te_bar_m(k) = te_bar_m(k)*norm_bar
            ti_bar_m(k) = ti_bar_m(k)*norm_bar
          endif
         enddo
        else
         norm_bar = rhosda_loc_s(jin_m)**2
         do k=0,jmaxm
          gamma_mks(k) = anrate_m(k)*csda_loc_s(jin_m)
          dgamma_mks(k)= dnrate_m(k)*csda_loc_s(jin_m)
          freq_mks(k)  = anfreq_m(k)*csda_loc_s(jin_m)
          ky_mks(k) = ky_m(k)/(rhosda_loc_s(jin_m)*a_unit_exp)
          if(igks_model.eq.1)then
            phi_bar_m(k) = phi_bar_m(k)*norm_bar
            ne_bar_m(k) = ne_bar_m(k)*norm_bar
            te_bar_m(k) = te_bar_m(k)*norm_bar
            ti_bar_m(k) = ti_bar_m(k)*norm_bar
          endif
         enddo
        endif
       else   ! igks >= 3
        if(igeo.eq.0)then
         do k=0,jmaxm
          gamma_mks(k) = anrate_m(k)*csda_exp(k)
          dgamma_mks(k)= dnrate_m(k)*csda_exp(k)
          freq_mks(k)  = anfreq_m(k)*csda_exp(k)
          ky_mks(k) = ky_m(k)/(rhosda_exp(k)*arho_exp)
          if(igks_model.eq.1)then
            norm_bar = rhosda_exp(k)**2
            phi_bar_m(k) = phi_bar_m(k)*norm_bar
            ne_bar_m(k) = ne_bar_m(k)*norm_bar
            te_bar_m(k) = te_bar_m(k)*norm_bar
            ti_bar_m(k) = ti_bar_m(k)*norm_bar
          endif
         enddo
        else
         do k=0,jmaxm
          gamma_mks(k) = anrate_m(k)*csda_loc_s(k)
          dgamma_mks(k)= dnrate_m(k)*csda_loc_s(k)
          freq_mks(k)  = anfreq_m(k)*csda_loc_s(k)
          ky_mks(k) = ky_m(k)/(rhosda_loc_s(k)*a_unit_exp)
          if(igks_model.eq.1)then
            norm_bar = rhosda_loc_s(k)**2
            phi_bar_m(k) = phi_bar_m(k)*norm_bar
            ne_bar_m(k) = ne_bar_m(k)*norm_bar
            te_bar_m(k) = te_bar_m(k)*norm_bar
            ti_bar_m(k) = ti_bar_m(k)*norm_bar
          endif
         enddo
        endif
       endif
      endif
c
c
c
c NETCDF printout
c
      if ( lprint_cdf .gt. 0) call write_gksout
c
      if(lprint_cdf.eq.0)then
        write(11,*)
     * "tprim1 = ",tprim1,
     * "tprim2 = ",tprim2,
     * "tprim3 = ",tprim3,
     * "fprim1 = ",fprim1,
     * "fprim2 = ",fprim2,
     * "fprim3 = ",fprim3,
     * "shat = ",shat,
     * "shift = ",shift,
     * "zeff = ",zeff,
     * "aky1 = ",aky1,
     * "teti = ",teti,
     * "pk = ",pk,
     * "an(1) = ",an(1),
     * "an(2) = ",an(2),
     * "an(3) = ",an(3),
     * "an(4) = ",an(4),
     * "beta = ",beta,
     * "uprim1 = ",uprim1,
     * "arho_exp = ",arho_exp,
     * "drhodr(25) = ",drhodr(25),
     * "elong_exp(25) = ",elong_exp(25),
     * "zpte_exp(25) = ",zpte_exp(25),
     * "rho(25) = ",rho(25),
     * "psi_exp(25) = ",psi_exp(25),
     * "rmin_exp(25) = ",rmin_exp(25),
     * "rmaj_exp(25) = ",rmaj_exp(25),
     * "delta_exp(25) = ",delta_exp(25),
     * "q_exp(25) = ",q_exp(25),
     * "betat_exp(25) = ",betat_exp(25)
      endif
c
      close(11)
c
c
 20   format (a)
 35   format(' *** Warning: jout_m greater than grid size ***')
 50   format(a10,2x,a4,a6)
 100  format(' Using ITER/TRANSP ufiles ...')
 105  format(' Using ONETWO iterdb file ...')
 107  format(' Creating test profiles ...')
 109  format(' Using Std test case parameters ...')
 125  format(/,' Calling GKS code ...',/)
c
      end
c*********************************************************
c
c
      subroutine set_gks_defaults
c*******************************
c
c     sets the default values for all of the
c     namelist variables gksin,gks0in,gksin1
c
c*******************************
      use gks_var
      implicit none
c
      include 'input.m'
cc      include 'data_d.m'
      include 'data_exp.m'
      include 'glf.m'
c
c****  gksin 
      ntheta=32
      nperiod=2
      ipar=1
      ngauss=5
      icv=1
      negrid=10
      nspec=5
      eps=.16666666
      shift=0.0
      dbeam=0.0
      shat=1.0
      pk=0.3333333333
      epsa=.33333333
      epsl=0.666666666
      width0=3.
      beta=.00000000001
      zeff=1.0
      teti=1.0
      zeff5=0.0
      fprim1=1.0
      fprim2=1.0
      fprim3=1.0
      fprim4=1.0
      fprim5=1.0
      tprim1=3.
      tprim2=1.0
      tprim3=3.
      tprim4=10.
      tprim5=1.0
      vnewk1=0.000
      vnewk2=0.000
      vnewk3=0.000
      vnewk4=0.000
      vnewk5=0.000
      bakdif1=0.0
      bakdif2=0.0
      bakdif3=0.0
      bakdif4=0.0
      bakdif5=0.0
      amass2=6.0
      amass3=0.0002723   ! 5.446D-4/2.0 for deuterium
      amass4=1.00
      amass5=28.60
      temp2=1.0
      temp3=1.0
      temp4=5.952
      temp5=1.0
      z2=6.0
      z3=-1.0
      z4=1.0
      z5=27.0
      uprim1=0.0
      uprim2=0.0
      uprim3=0.0
      uprim4=0.0
      uprim5=0.0
      mach1 = 0.0
      mach2 = 0.0
      mach3 = 0.0
      mach4 = 0.0
      mach5 = 0.0
c
cc**** gks0in
c
      nscreen=25
      nout=100
      nstep=201
      ngamstep=1
      isvar=1
      nwrite=401
      iphi=0
      delt=0.1
      ominst=1.e09
      phydif=1.D50
      power=-0.5
      tol=-.005
      gamma=0.1
      absom=0.5
      fcv=0.0
      fv=0.0
      gridfac=4.e04
      fexp1=0.4
      fexp3=0.4
      test1=1.0
      test2=1.0
      alr=0.0
      ali=0.0
      al1r=0.0
      al1i=0.0
      alar=0.0
      alai=0.0
      al1ar=0.0
      al1ai=0.0
      ncspec1=1
      ncspec2=1
      ncspec3=1
      ncspec4=0
      ncspec5=0
      nce1=0
      nce2=0
      nce3=0
      nce4=0
      nce5=0
      aky1=0.42426
      aky2=0.28284
      aky3=0.42426
      aky4=0.56568
      aky5=0.7071
      aky6=0.8472
      aky7=0.98994
      aky8=1.1316
      aky9=1.27278
      aky10=1.4142
      aky11=0.3
      aky12=0.3
      aky13=0.3
      aky14=0.3
      aky15=0.3
      aky16=0.3
      naky=1
      nt0=0
      thetamin=0.0
      thetamax=0.0
      nt0min=0
      nt0max=0
c      ecut=2.5
      ecut = 5.0
      adamp=-0.25
      nteststep=20
c
c**** gksin1
c
      igks_model=1
      idata = 1
      bt_flag=0
      tok='d3d'
      shot='84736'
      delt = 0.1
      xp_time=1.5
      cudir="/u/staebler/GKS"
      i_solve=0
      i_bpar=1
      igeo=1
      igeot=0
      icontinue=0
      cdebye=0.0
      cbetae=1.0
      xalpha=1.0
      alpha_mhd_loc=0.0
      ibranch=0
      y_bpar=1.0
      x_bpar=1.
      xy_bpar=1.
      y_mhd=0.0
      signomega=-1
      igyro_fix=0
      igyro_e = 0
      kys0=0.3
      tol_f=0.2
      xkymin_gf=0.02
      xkymax_gf=2.0
      iptot=0
      ivar=0
      jvar=0
      irunmax=1
      jrunmax=1
      xvarmin=0.0
      xvarmax=0.0
      zvarmin=0.0
      zvarmax=0.0
      igraph=0
      mprint=0
      irot1=1
      irot2=1
      igeo_m=3      
      i_dengrad=2
      jout_m=50
      jin_m=0
      echconv=0.0
      igks=2
      nstep=301
      amassgas_exp=2.0
      lprint_cdf=1
      z2=6
      z3=-1
      aky1=0.42426
      corot=1.0
      igeo_print=0
      itorque=0
      cnewk2=0.0
      cnewk3=1.0
      nbasis_min=1
      nbasis_max=4
      alpha_mach=0.0
      alpha_cur=0.0
      alpha_e = 0.0
      alpha_p=1.0
      alpha_p_cur=0.0
      lprint_pflow=0
      iexb=0
      ipfst=0
      idatzero=0
      itest_ntcc=0
      ialpha=0
      save_tglf=.FALSE.
      overwrite_tglf=.FALSE.
      ismooth_data=0
      tglf_defaults=1
      gks_defaults=1
c
      return
      end
c**********************************
      subroutine geogks_toq

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cmnt  read_toq to setup geogks_driver
cmnt  gam_exp->1 is flat density profile
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       implicit none
       include 'data_exp.m'
       include 'input.m'
       real*8 nxt,gam_exp
       integer j
c
       call read_toq
c
       gam_exp=1.5D0
c       arho_exp=sqrt(torfO2pi_toq*2.D0/btvac_toq)/a_bnd_toq
crenorm field to 1. bt_vac_toq->1.
       bt_exp=1.D0
crmin_exp=a_unit_exp->1.
       elonga_exp=elong_exp(50)  
       do j=1,50
cnote beta_loc_0= 400.D0*(ne_exp(j)*te_exp(j)+ni_exp(j)*ti_exp(j))/
c               (1.e5*bt_mag_center**2*b_div**2)/ab_div 
       nxt=betat_exp(j)/400.D0/2.D0*(1.D5*bt_exp**2)
       ne_exp(j)=1.D-6*(nxt+1.D-10)**(1.-gam_exp)
       te_exp(j)=1.D+6*(nxt+1.D-10)**gam_exp
       ni_exp(j)=ne_exp(j)
       ti_exp(j)=te_exp(j)
c this insures correct beta but collisionality vanishing
     
       zeff_exp(j)=1.D0
       nz_exp(j)=1.D-6*ni_exp(j)
       enddo
       return
       end
cc***********************************************************************
       subroutine read_toq
c***********************************************************************
cmnt   reads profiles  from TOQ Miller's MHD quilibrium code
c***********************************************************************
       implicit none
       include 'input.m'
       include 'data_exp.m'
       integer nrho,nrhom1,i
  
       open(unit=75,file='toq.dat',status='old')
                      
       read(75,'(i5,4e13.5)') nrho,btvac_toq,btplas_toq,
     *     torfO2pi_toq,a_bnd_toq
       nrhom1=nrho-1
       read(75,'(6e13.5)')(rho(i),i=0,nrhom1)
       read(75,'(6e13.5)')(psi_exp(i),i=0,nrhom1)
       read(75,'(6e13.5)')(rmin_exp(i),i=0,nrhom1)
       read(75,'(6e13.5)')(rmaj_exp(i),i=0,nrhom1)
       read(75,'(6e13.5)')(elong_exp(i),i=0,nrhom1)
       read(75,'(6e13.5)')(delta_exp(i),i=0,nrhom1)
       read(75,'(6e13.5)')(q_exp(i),i=0,nrhom1)
       read(75,'(6e13.5)')(betat_exp(i),i=0,nrhom1)
       arho_exp=sqrt(torfO2pi_toq*2.D0/btvac_toq)/a_bnd_toq      
 
       close(75)
       return
       end
c***********************************************************************
 

