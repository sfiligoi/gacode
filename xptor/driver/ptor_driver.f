c@ptor_driver.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... This is the driver routine for XPTOR
c
c    Settings input:
c       1) files 'xpin' and 'xp0in' - provide the default settings
c       2) file 'in' - model and transport switches
c    Data input:
c       1) idata=0 for reading ITER PDB ufiles with readufiles.f
c       2) idata=1 for reading Onetwo iterdb files with readiterdb.f
c       3) idata=2 for reading Onetwo netcdf files with cdfread.f
c       4) idata=3 for creating test profiles with readtestprofiles.f
c       5) ipptot=1 to include fast ions in total pressure if PTOTR not given
c    Transport models:
c       1) Anomalous transport
c          The anomalous transport model is set using the variable 'imodel'
c
c            imodel      transport model
c              0             simple
c              2             IFS/PPPL
c              3             IIF (CDBM)
c              4             Culham
c              6             MM95
c              66            Weiland
c              81            GLF23    (iglf=0 for original version, 1 for new)
c              82            TGLF (xnu_model=1 for APS07, 2 for APS-09)
c              99            chi=1 test (chii=chie=1.0 m^2/s)
c       2) Neoclassical transport
c          There are 3 models for neoclassical transport which are selected
c          using the variable 'use_xneo_m'
c
c            use_xneo_m          model
c               -1         exponentially decaying artificial chi (ineo=1 for chie=chii)
c               0,1        KAPISN w/ experimental/model profiles
c               2,3        NCLASS w/ experimental/model profiles
c
c          Within the KAPISN routine, the neoclassical model is set
c          using the variable 'nkimod_nc'
c
c            nkimod_nc          transport model
c               1           Rutherford
c               2           modified Hazeltine-Hinton (default)
c               3           Bolton
c               4           original Chang-Hinton
c               5           modified Chang-Hinton for Zeff>1
c          To use chi-phi-neo=delta^3/2*chii-neo, set ineophi=1.
c
c       3) ExB shear stabilization
c             iexb
c              0         ExB shear using exp. vphi, neoclassical v_theta (Kim formula)
c              1         ExB shear using omega.dat file, stored in megamma_exp
c              3         NCLASS Hahm-Burrell formula, egamma_exp=egamma_ncl_exp
c
c       4) Shafranov shift stabilization:
c       To include the effect of Shafranov shift stabilization in GLF23
c       the variable 'xalpha' is used as follows:
c            xalpha
c              0         no alpha-stabilization
c             >0         use exp profiles for alpha
c             <0         compute alpha from model profiles w/ fast ion pressure
c                        if using ITER ufiles. Set ipptot=1 to compute ptot w/
c                        n_fast, t_fast=ti, otherwise compute ptot as thermal pressure
c    Artificial background diffusivities:
c       adiff_dv > 0 sets background chie, chii (m^2/s)
c       adiffphi_dv > 0 sets background chi-phi (m^2/s)
c    Updating sources:
c       The sources from ITER ufiles can be updated as a function of time.
c       1) set time_upd for simulation time when updating sources is to occur
c       2) set upsrc(1)=1 to update q-profile, upsrc(2)=1 to update power flows
c    Toroidal momentum source:
c       Setting itorque=1 uses torque density from ONETWO iterdb file.
c       Otherwise, momentum source computed from ion beam deposition
c    Radiation:
c       Bremmstrahlung and synchrotron radiation can be computed by
c       setting irad=-1 (=1 for no printout to screen). To compute
c       the synchrotron radiation, set the reflection coefficient
c       (radref) greater than zero. For example, radref=0.20 corresponds
c       to 20% lost from plasma and 80% reabsorbed.
c    Boundary conditions:
c       Experimental boundary conditions are enforced by setting jout_m
c       to the desired zone location. Pedestal scalings for T_ped may also
c       be used as follows:
c          ibound = 1 Eqn12 power independent MHD scaling (Thomsen,Cordey,IAEA02)
c          ibound = 2 Eqn1 power dependent scaling (Thomsen,Cordey)
c          ibound = 3 Wped1 power dependent scaling (Thomsen,Cordey,NF), Paux only
c          ibound = 4 Wped2 power dependent scaling w/o type III ELMy data (Thomsen,Cordey,NF)
c          ibound = 5 Wped1 power dependent scaling w/ Paux+Palpha
c          ibound = 6 Wped2 power dependent scaling w/o type III ELMy data, w/ Paux+Palpha
c          cped = arbitrary multiplier on pedestal scaling
c    Grid size:
c       The default grid size is 50 points
c       Currently, the grid size may only be changed when using ITER PDB ufiles.
c       Interpolation needs to be added for use with ONETWO iterdb files.
c       2) set jmaxm,mxgrid to desired size in file 'in' (defaulted to 50 in xpin)
c       3) set grid no. for BC in 'in' (default = 45 in xpin)
c       Note: Maximum number of time pts = 3300
c    Smoothing:
c       Profile smoothing of the experimental data is possible using the
c       FITPACK tension spline routine, curvss. The variables ismoo_* control
c       the level of smoothing where larger values correspond to more smoothing.
c          ismoo_ne > 0 for smoothing of the electron density
c          ismoo_ni > 0 for smoothing of the ion density
c          ismoo_nf > 0 for smoothing of the fast ion density
c          ismoo_zeff > 0 for smoothing of Zeff
c          ismoo_q > 0 for smoothing of the q-profile
c          ismoo_qnb > 0 for smoothing of the beam deposition profile (qbeame,i, sbeam)
c          ismoo_torque > 0 for smoothing of the torque density profile
c          ismoo_grho > 0 for smoothing of the gradrho, gradrho^2 profiles
c          ismoo_delta > 0 for smoothing of the triangularity
c          ismooth_all > 0 to smooth all experimental profiles
c       The model ion density and temperatures can also be smoothed, but not recommended.
c       It may be desirable, for example, to smooth the profiles just prior
c       to the end of the simulation in order to obtain smoother diffusivities
c       and growth rates. To do so, set istep_smoo to the time-step and
c          smoo_nim > 0 to smooth ni
c          smoo_tim > 0 to smooth both Te and Ti (a few percent is sufficient)
c          smoo_vphim > 0 to smooth vphi
c    Spatial averaging:
c       Both the fields and ExB shear velocities can be spatially averaged
c       over their neighboring pts where
c          x(j)=(1-ave)*x(j)+ave*(x(j-1)+x(j+1))
c       The averaging parameters are:
c          ave_dv > 0 for averaging of fields w/ Staebler DV solver
c          ave_field > 0 for averaging of fields 1-5
c          ave_ve > 0 for averaging of ve
c    Wdot terms:
c       The wdot and sdot terms are included by default. To exclude
c       them in the total power flows, set xwdot=1 for DWIR, DWER and xsdot=1 for DNER.
c    Fusion power:
c       Switches: ialpha, idt, ifusmodel,fuscale
c       XPTOR can compute the alpha power by setting ialpha=1. The default
c       calculation of the alpha power (ifusmodel=0) uses the NRL formula.
c       This was revised 3/15/01 to give a better fit for temperatures
c       between 7 and 100keV (jek). Setting idt=1 allows the user to compute
c       the alpha power using the deuterium and tritium densities provided
c       in the ufile data, for example, when the mix is not 50:50.
c       For ifusmodel=1, the alpha power is computed using the expression
c       by D. Boucher from the PRETOR code. The default maximum alpha power
c       is limited to 300 MW using the variable pfusion_max in subroutine setup.f.
c       The variable fuscale is used to scale the computed alpha power to
c       match the experimental value since code does not compute beam-thermal
c       or beam-beam interactions.
c       Note: E_crit assumes a pure 50:50 D:T plasma.
c       The total alpha power is printed to the screen in MW.
c    Output:
c       Output switches:
c          lprint_cdf   =1 writes all output to a NetCDF file 'xptor.cdf'
c                       >1 writes chi's, egamma, anrate, etc in addition to te,ti,vphi,vexb
c          lprint_prof  =1 writes profile output to ascii file 'out'
c                       =3 writes iterdb formatted profiles to iterdbprofs.dat
c          lprint_pulse =1 writes time-dependent ascii data to tepulse.dat, tipulse.dat
c          lprint_glf   =1 writes ascii data for input to GLF23 stand-alone code
c          lprint_glf2  =1 writes GLF23 diagnostic printout to file 'glf2.out'
c          lprint_smoo = 1 to print unsmoothed and smoothed values to screen
c          lprint_ufile = 1 to print ascii ufiles for profiles 
c                        (array: 1=Te only, 2=Ti, 3=ne)
c       ASCII Output files:
c         1) out - printout of profiles and diffusivities at final time
c         2) t0.dat - Te(0) and Ti(0) versus time
c         3) tepulse.dat, tipulse.dat - Te and Ti at various radii vs time
c         4) vphipulse.dat - toroidal velocity at various radii vs time
c         5) exb.dat - various contributions to ExB shear rate
c         6) vexb.dat - ExB shear velocity at various radii vs time
c         7) vetor.dat - vphi component of ExB shear velocity at various radii vs time
c         8) vepol.dat - vpol component of ExB shear velocity at various radii vs time
c         9) vstar.dat - vdia component of ExB shear velocity at various radii vs time
c        10) egamma.dat - gamma_ExB at various radii vs time
c        11) anrate.dat - gamma_max at various radii vs time
c        12) etaphi.dat - etaphi at various radii vs time
c        13) gammap.dat - parallel velocity shear at various radii vs time
c        14) egamma-vphi.dat, egamma-vpol.dat, egamma-vstar.dat - components of gamma_E
c        15) vout - vphi_exp(j), vexb_exp(j), vphi_m(j), vexb_m(j) vs radii to jout_m
c        16) vexbpulse.dat - vexb_t vs time
c        17) chie.dat - chie (m^2/s) at various radii vs time
c        18) chii.dat - chii (m^2/s) at various radii vs time
c       Diagnostic output:
c         Printout of the growth rates, frequencies and chi's from GLF23
c         at a particular time-step and grid point is possible by setting
c         lprint_glf2=1 and setting
c              lprint_glfstep = time-step desired (default=icallei)
c              lprint_glfgrid = grid point desired (default=25)
c         itest_ncl prints out diffusivities from NCLASS in readiterdb.f if set > 0
c    Restart capability:
c       To restart a simulation, set restart_pt.ge.0.D0
c       Currently, restarts are only operable if initial run had lprint_pulse=1
c       set to write out ascii *.dat files. Restart simulation must have lprint_pulse=0.
c       * Need to update this for NetCDF file read.
c    Running: 
c       Normally this is carried out using PBS batch submissions with the
c       batch.src file in xptor/bin.
c       To execute XPTOR interactively using MPI with two processors,
c       for example, ... mpirun -np 2 ./xptor
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      program ptor_driver
c
c      use tglf_pkg
      implicit none
c
      include 'mpif.h'
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/ptor.m'
      include '../inc/glf.m'
      include '../inc/vtrans.m'
      include '../inc/glf_dimensions.m'
      include '../inc/glf_hermite.m'
      include '../inc/gks.m'
c
      character line*132
      character*64, eqfile, msg
      integer iflag, ieqdsk, igrid
      integer i, j, k, lprint_cdf, lprint_prof, lprint_gks,
     &        lprint_pulse, lprint_glf, lprint_ufile(5), itgflag
      integer kmx, icallei, istkrun, iei0, ieicnt, j1, j2, jgeo_gks
      integer ishoot_hold,jout_m_hold,jin_m_hold
      integer jin_hold,jout_hold
      integer irun, irunmax, izrunmax, jrunmax
      integer iscans, iscanid, ixvar, izvar
      integer iglf2d_scan, iglf2d_jrun
      real xvarmin, xvarmax, xvar, xvarmax_h, rlti_gf_h, rlte_gf_h,
     &     cg_gf, gfac_gf, drhodr_gf, grhosq_exp_av,
     &     sfactor_av, flow_exp_av, powimm
      real*8 xvmin, xvmax, zvmin, zvmax
      parameter (kmx=3300,izrunmax=20,jrunmax=100)
c      real*8 xp_gk(jrunmax),yp_gk(izrunmax,jrunmax,3),
c     &       zp_gk(jrunmax)  ! put in gks.m
      real*8 xp_gf(jrunmax),yp_gf(izrunmax,jrunmax,3),zp_gf(jrunmax),
     &       tm_gf(izrunmax,jrunmax)
      real*8 dev_gkgf(3), stdev_gkgf(3),
     &       devsq_gkgf(3), expsq_gkgf(3),
     &       stdevtot_gkgf(3), stdevtot_smgkgf, rms_theta_max,
     &       dev_smgkgf, devsq_smgkgf, expsq_smgkgf, stdev_smgkgf,
     &       dev_smtgkgf, devsq_smtgkgf, expsq_smtgkgf, gamma_limit
      real*8 glf2d_zvar, glf2d_zvarmin,  glf2d_zvarmax
      real*8 stdev_zvar_gkgf(jrunmax,2)
      real*8 gkdat(jrunmax), gfdat(jrunmax),
     &       gkdat2(jrunmax), gfdat2(jrunmax)
      real*8 iota(0:60), stdev(4), relrms(4), offset(4), datmax(4)
      real*8 pow_ped,w_ped(8)
      real*8 score, rms_theta_save, chi_det
      double precision starttime, endtime,cputime
c
      real*8 rmin_gfm(mxgrd), rmaj_gfm(mxgrd), elong_gfm(mxgrd),
     &       rlne_gfm(mxgrd), rlni_gfm(mxgrd), rlte_gfm(mxgrd),
     &       rlti_gfm(mxgrd), taui_gfm(mxgrd), q_gfm(mxgrd),
     &       shat_gfm(mxgrd), alpha_gfm(mxgrd), xnu_gfm(mxgrd),
     &       xwell_gfm(mxgrd), park_gfm(mxgrd), ghat_gfm(mxgrd),
     &       gchat_gfm(mxgrd), adamp_gfm(mxgrd), alphae_gfm(mxgrd),
     &       gammae_gfm(mxgrd), alphap_gfm(mxgrd), gammap_gfm(mxgrd),
     &       powem_gfm(mxgrd), powim_gfm(mxgrd)
c
      character*12 ZDATE
      character*10 UDATE
      character*30 foo1,foo2
      integer imon, ivals(8)
      character*3 mons(12)
c
c namelists ....................................................
c
      include 'inpt.nml'
      include 'xpin.nml'
      include 'xp0in.nml'
c
      data mons/'jan','feb','mar','apr','may','jun','jul','aug',
     >          'sep','oct','nov','dec'/
c
c initialize MPI................................................
c
      call MPI_INIT(i_err)
      call MPI_COMM_RANK(MPI_COMM_WORLD,i_proc,i_err)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc,i_err)
      starttime = MPI_WTIME()
      if(i_proc.eq.0)write(*,*)"mpirun with ",n_proc," processors"
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      open (4,file='xpin')
      open (5,file='xp0in')
      open (3,file='in')
c      if (i_proc .eq. 0) then
c        open (30,file='vout',status='unknown')
c        open (31,file='vperpulse.dat',status='unknown')
c      endif
c
      ZDATE=' '
      call date_and_time(udate,foo1,foo2,ivals)
      read(udate(5:6),'(i2)') imon
      zdate=udate(7:8)//mons(imon)//udate(1:4)
c     call clock()
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c read input files xpin, xp0in
c
  10  continue
      read  (4,20,end=900,err=900) line
c
      if ( index ( line, '$xpin' ) .gt. 0
     &    .or. index ( line, '&xpin' ) .gt. 0 ) then
      backspace 4
      read(4,xpin)
      else
        go to 10
      endif
c
  30  continue
      read  (5,20,end=900,err=900) line
c
      if ( index ( line, '$xp0in' ) .gt. 0
     &    .or. index ( line, '&xp0in' ) .gt. 0 ) then
      backspace 5
      read(5,xp0in,end=32)
      else
        go to 30
      endif
c
 32   continue
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c setup for glf23, read 'in' input file
c
      if (imodel.eq.8 .or. imodel.eq.81) then
        call setup
      endif
c
      call init
c
      ieqdsk=0
      igrid=0
      iglf2d_jrun=1
      ipert_gf=0
      lprint_cdf=0
      lprint_time_cdf=0
      lprint_prof=0
      lprint_glfgrid=25
      lprint_ufile(1)=0
      lprint_ufile(2)=0
      lprint_ufile(3)=0
      tauscale(1)=1
      tauscale(2)=1
      tauscale(3)=1
      tauscale(4)=1
      taupin=0.25D0
      time_neut=999.D0
      irunmax=40
      xvarmin=1.D-5
      xvarmax=2.D0
c
c... read 'in' file
c
  38  continue
      read  (3,20,end=900,err=900) line
c
      if ( index ( line, '$inpt' ) .gt. 0
     &    .or. index ( line, '&inpt' ) .gt. 0 ) then
        backspace 3
        read(3,inpt)
      else
        go to 38
      endif
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      if(mxgrid.ne.50)then
        if(idata.eq.1)then
c          if(i_proc.eq.0)write(6,*)"error: mxgrid not 50 and idata=1"
c          stop
        endif
        if(mxgrid .gt. 300)then
          if(i_proc .eq. 0)write(6,*)"error: mxgrid must be <= 300"
          stop
        endif
      endif
      if(dvflag.eq.0)then
        nsteps_v=icallei
        lprint_glfstep=icallei
      endif
      jin_hold=jin_m
      jout_hold=jout_m
c
      if(jout_m.gt.jmaxm .and. i_proc .eq. 0) then
         write(6,35)
      endif
      if(i_proc .eq. 0 .and. igfac .gt. 0)
     >   write(6,*) '*** using Garys gfac ***'
      if (idata.eq.3) then
        tok = 'test'
        shot= '999'
        phase='X'
        xp_time=99.D0
      endif
c
      if (i_proc.eq.0 .and. iglf2d.ne.-1)
     >    write(*,50) zdate
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c read in experimental data
c
c     lprint_prof  =2  ! profile printout
c     lprint_glf   =1  ! printout GLF23 diagnostic data
c     lprint_pulse =1  ! printout temperature traces
c
      if(i_proc.eq.0 .and. iglf2d.ne.-1
     >   .and. idata.ne.-1) write(*,110) tok,shot,xp_time
      if(jmaxm.gt.nj .and. i_proc .eq. 0) then
         write(*,60) jmaxm, jout_m
      endif
c
c
c read eqdsk data
      if (ieqdsk.gt.0) then
        eqfile = 'eqdskin'
        call read_eqdsk(eqfile,iflag,msg)
      endif
c
      if (idata .eq. 0) then
        if (i_proc.eq.0) then
          write(*,100)
        endif
        call readufiles(tok,shot,phase,cudir,mxgrid,ismooth,btscale,
     >                xp_time,endtime_pt,time_series,iexp_exch,
     >                iexp_q,p_glob_exp,gtauth_exp,gtautot_exp,i_proc)
c        if(ismooth_all.ne.0) call datavg ! smoothing in readufiles for now
      elseif (idata .eq. 1) then
        if (i_proc.eq.0) then
          write(*,105)
        endif
        call readiterdb
        if(ismooth_all.ne.0) call datavg
        call datexch
        if(ineutp.eq.1) call neutxp
      elseif (idata .eq. 2) then
        if (i_proc.eq.0) then
          write(*,106)
        endif
        call cdfread(igrid)
      elseif (idata .eq. 3) then
        if (i_proc.eq.0) then
          write(*,107)
        endif
        call readtestprofiles
      elseif (idata .eq. -1) then
        if (i_proc.eq.0) then
          write(*,109)
        endif
      else
          write(*,*) 'Error: idata out of range'
                  stop
      endif
      if (iexb.eq.2) call read_er
c
c... interpolate to new grid
c
c      if(idata.ne.0) call gridsetup(igrid)
      call gridsetup(igrid)
c
c... tension spline smoothing
c
      if (ismoo_ne.gt.0 .or. ismoo_ni.gt.0 .or. ismoo_nf.gt.0 .or.
     &    ismoo_zeff.gt.0 .or. ismoo_q.gt.0 .or. ismoo_qnb.gt.0 .or.
     &    ismoo_vrot.gt.0 .or. ismoo_delta.gt.0 .or. 
     &    ismoo_grho.gt.0) call tenspline
      if (i_proc.eq.0)write(*,*) 'call to tenspline successful'
c
c... call neoclassical using exp data
c
      imyneoclass=0
      if (use_xneo_m.eq.0 .or. use_xneo_m.eq.1 .or.
     &    use_xneo_m.eq.-1) call datneo
      if(use_xneo_m.eq.2 .or. iexb.eq.3) call datnclass
c
c... map to _exp variables
c
      if(idata.ge.0 .and. idata.le.2) call datmap(ialpha)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c...rescale profiles (if desired)
c   first call derivedexpprofiles to get unscaled angrotp_exp
c   only if not call glf2d directly
c
      if(iglf2d.ne.-1.and.idata.ne.-1) then
        if(igyro.eq.1) then
          irot1=1
          irot2=1
          igeo_m=3
        endif
        call derivedexpprofiles(1,nj_d,nj,amin_d,tocur_d,
     &       rho_d,q_d,pow_ped,w_ped,ismoo_vrot)
      endif
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... use pedestal scaling to set boundary condition
c    set in subroutine ptorinit.f
c
      if(ibound.gt.0) then
        if(i_proc.eq.0) then
          write(*,'(/,a34)')
     &    'Using pedestal scaling for BC ...'
          write(*,'(a16,f8.4)') 'cped        = ',cped
          write(*,'(a16,f8.4)') 'nped        = ',ne_exp(jout_m)
          if(ibound.eq.1) then
            write(*,'(a16,f8.4)')   'Wped-MHD   = ',w_ped(1)
            write(*,'(a16,f8.4)')   'Tped-MHD   = ',t_ped(1)
          elseif(ibound.eq.2) then
            write(*,'(a16,f8.4)')   'pow_ped     = ',pow_ped
            write(*,'(a16,f8.4)')   'Wped1-iaea  = ',w_ped(2)
            write(*,'(a16,f8.4)')   'Tped1-iaea  = ',t_ped(2)
          elseif(ibound.eq.3 .or. ibound.eq.5) then
            write(*,'(a16,f8.4)')   'pow_ped     = ',pow_ped
            write(*,'(a16,f8.4)')   'Wped1-NF    = ',w_ped(3)
            write(*,'(a16,f8.4)')   'Tped1-NF    = ',t_ped(3)
          elseif(ibound.eq.4 .or. ibound.eq.6) then
            write(*,'(a16,f8.4)')   'pow_ped     = ',pow_ped
            write(*,'(a16,f8.4)')   'Wped2-NF    = ',w_ped(4)
            write(*,'(a16,f8.4)')   'Tped2-NF    = ',t_ped(4)
          endif
          write(*,'(a18,2f9.4,/)') 'Tped-exp-e,i  = ',te_exp(jout_m),
     &                           ti_exp(jout_m)
        endif
      endif
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c...write energy confinement scalings
c
      if (i_proc.eq.0.and.iglf2d.ne.-1.and.idata.ne.-1) then
        write(6,'(a26,F8.4)')   ' Pheat (mw)  :    exp   = ',
     &      pheat_exp
        write(6,'(a26,F8.4)')   ' wth (secs)  :    exp   = ',
     &      wstr_exp
        write(6,'(a26,F8.4,a10)')   ' taueth (secs)  : exp   = ',
     &        gtaut_exp,'(no prad)'
        write(6,'(a26,F8.4,a10)')   ' tauetot (secs) : exp   = ',
     &        gtautot_exp,'(no prad)'
        if(idata.eq.0) then
        write(6,'(a26,2F8.4)')   ' tauetot (secs): exp(0d)= ',
     &        gtautot_exp
        write(6,'(a26,2F8.4)')   ' taueth (secs) : exp(0d)= ',
     &        gtauth_exp
        endif
        if (tauscale(1).eq.1) then
          write(6,'(a34,2F8.4)') ' tauetot (secs) : ITER-89P       = ',
     &          tauescale_exp(1)
        endif
        if (tauscale(2).eq.1) then
          write(6,'(a34,2F8.4)') ' taueth (secs)  : ITER-97L       = ',
     &          tauescale_exp(2)
        endif
        if (tauscale(3).eq.1) then
          write(6,'(a34,2F8.4)') ' taueth (secs)  : ITER-98(y,2)   = ',
     &          tauescale_exp(3)
        endif
        if (tauscale(4).eq.1) then
          write(6,'(a34,2F8.4)') ' taueth (secs)  : IAEA04(DB3v13) = ',
     &          tauescale_exp(4)
        endif
      endif
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c create GYRO output 'expprofile'
c Note: set iexb=2 to read in exp Er data
c cub_spline is used in readiterdb for mapping
c
      if (igyro.gt.0) then
c
        write(*,*)
        write(*,*) 'Writing GYRO input file ...'
        write(*,*) 'setting irot2=1 for vphi, igeo_m=3'
        open (77,file='INPUT_profiles',status='unknown')
        write(77,'(a)') '#'
        write(77,'(a)') '# ** This file was generated by xptor **'
        write(77,'(a)') '#'
        write(77,'(a,a)') '#         INPUT FILE: ','iterdb.'//shot
        write(77,'(a,a)') '#        SHOT NUMBER: ',shot
        write(77,'(a,i3)') '#  RADIAL GRIDPOINTS: ',mxgrid+1
        write(77,'(a,/a,/a)') '#','#', '#'
        write(77,'(a)') '# Scalar data:'
        write(77,'(a)') '#'
        write(77,'(a,i3)') 'N_EXP=',mxgrid+1
        write(77,'(a,f7.5)') 'BT_EXP=',bt_exp
        write(77,'(a,f8.6)') 'ARHO_EXP=',arho_exp
        write(77,'(a)') '#'
        write(77,'(a,9x,a,8x,a,8x,a,11x,a)') '#rho(-)','rmin(m)',
     &            'rmaj(m)','q(-)','elong(-)'
        write(*,*) 'q(a) = ',q_exp(jmaxm)
        write(*,*) 'Ti(a) = ',ti_exp(jmaxm)
        write(*,*) 'Te(a) = ',te_exp(jmaxm)
        do j=0,jmaxm
          write(77,'(1pe13.6,1p4e15.6)')rho(j),rmin_exp(j),rmaj_exp(j),
     &          q_exp(j), elong_exp(j)
        enddo
        write(77,'(a)') '#'
        write(77,'(a,6x,a,8x,a,2x,a,7x,a)') '#delta0(-)','Te(keV)',
     &            'ne(10^19/m^3)','z_eff(-)','Er(kV/m)'
        do j=0,jmaxm
          write(77,'(1pe13.6,1p4e15.6)') delta_exp(j),te_exp(j),
     &          ne_exp(j), zeff_exp(j), er_exp(j)
        enddo
        write(77,'(a)') '#'
        write(77,'(a,3x,a,6x,a,6x,a,5x,a)') '#flow_mom(Nm)','pow_e(MW)',
     &            'pow_i(MW)','pow_ei(MW)','delta1(-)'
        do j=0,jmaxm
          write(77,'(1pe13.6,1p4e15.6)') delta_exp(j),powe_exp(j),
     &          powi_exp(j), pow_ei_exp(j), 0.0
        enddo
        write(77,'(a)') '#'
        write(77,'(a,6x,a,6x,a,12x,a,12x,a)') '#flow_beam',
     &         'flow_wall','[-]','[-]','[-]'
        do j=0,jmaxm
          write(77,'(1pe13.6,1p4e15.6)') flow_beam_exp(j),
     &          flow_wall_exp(j), 0.0, 0.0, 0.0
        enddo
        write(77,'(a)') '#'
        write(77,'(a,a,a,a,a)') '#ni_1(10^19/m^3)',
     &         'ni_2(10^19/m^3)','ni_3(10^19/m^3)',
     &         'ni_4(10^19/m^3)','ni_5(10^19/m^3)'
        do j=0,jmaxm
          write(77,'(1pe13.6,1p4e15.6)') ni_exp(j),
     &          nz_exp(j), 0.0, 0.0, 0.0
        enddo
        write(77,'(a)') '#'
        write(77,'(a,5x,a,5x,a,5x,a,5x,a)') '#Ti_1 (keV)',
     &         'Ti_2 (keV)','Ti_3 (keV)', 'Ti_4 (keV)','Ti_5 (keV)'
        do j=0,jmaxm
          write(77,'(1pe13.6,1p4e15.6)') ti_exp(j),
     &          ti_exp(j), 0.0, 0.0, 0.0
        enddo
        write(77,'(a)') '#'
        write(77,'(a,3x,a,3x,a,3x,a,3x,a)') '#vtor_1 (m/s)',
     &    'vtor_2 (m/s)','vtor_3 (m/s)', 'vtor_4 (m/s)','vtor_5 (m/s)'
        do j=0,jmaxm
          write(77,'(1pe13.6,1p4e15.6)') vphi_exp(j),
     &          0.0, 0.0, 0.0, 0.0
        enddo
        write(77,'(a)') '#'
        write(77,'(a,3x,a,3x,a,3x,a,3x,a)') '#vpol_1 (m/s)',
     &    'vpol_2 (m/s)','vpol_3 (m/s)', 'vpol_4 (m/s)','vpol_5 (m/s)'
        do j=0,jmaxm
          write(77,'(1pe13.6,1p4e15.6)') 0.0,
     &          0.0, 0.0, 0.0, 0.0
        enddo

        close(77)
        goto 900
      endif
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c restart option
c
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c write out central values at beginning of transport run
c
      if (i_proc .eq. 0) then
        if (itport_pt(1).gt.0) then
          write(*,70) ne_exp(0), ni_exp(0), te_exp(0), ti_exp(0),
     &                vphi_exp(1)
        else
          write(*,75) te_exp(0), ti_exp(0), vphi_exp(1)
        endif
      endif
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      do j=0,60
        iota(j)=j
      enddo
c
      if(xmult.gt.0) then
       if(temult.gt.0) then
        do j=0,jout_m
          te_exp_sav(j)=te_exp(j)
          te_exp(j)=temult*te_exp(j)*(1-(rho(j)/rho(jout_m)))**xmult+
     &         te_exp(jout_m)
        enddo
       endif
       if(timult.gt.0.) then
        do j=0,jout_m
          ti_exp_sav(j)=ti_exp(j)
          angrotp_exp_sav(j)=angrotp_exp(j)
          ti_exp(j)=timult*ti_exp(j)*(1-(rho(j)/rho(jout_m)))**xmult+
     &         ti_exp(jout_m)
          if(itport_pt(4).gt.0) then
            angrot_exp(j)=timult*angrot_exp(j)*
     &         (1-(rho(j)/rho(jout_m)))**xmult+angrot_exp(jout_m)
          endif
        enddo
       endif
      else
       if(temult.gt.0.) then
        do j=0,jmaxm
          te_exp_sav(j)=te_exp(j)
          te_exp(j)=temult*te_exp(j)
        enddo
       endif
       if(temult.lt.0.) then
        do j=0,jmaxm
          te_exp_sav(j)=te_exp(j)
          te_exp(j)=-temult*te_exp(jout_m)*
     &         (1.D0+1.D-6*(jout_m-iota(j)))
        enddo
       endif
       if(timult.gt.0.) then
        do j=0,jmaxm
          ti_exp_sav(j)=ti_exp(j)
          angrotp_exp_sav(j)=angrotp_exp(j)
          ti_exp(j)=timult*ti_exp(j)
          if(itport_pt(4).gt.0) then
            angrot_exp(j)=timult*angrot_exp(j)
          endif
        enddo
       endif
       if(timult.lt.0.) then
        do j=0,jmaxm
          ti_exp_sav(j)=ti_exp(j)
          ti_exp(j)=-timult*ti_exp(jout_m)*
     &          (1.D0+1.D-6*(jout_m-iota(j)))
        enddo
       endif
       do j=0,jmaxm
         zeff_exp(j)=zeff_exp(j)+zeff_hold
         vphip_exp_sav(j)=angrotp_exp_sav(j)*rmajor_exp
       enddo
      endif
c
      call derivedexpprofiles(2,nj_d,nj,amin_d,tocur_d,
     &            rho_d,q_d,pow_ped,w_ped,ismoo_vrot)
c
      if(dvflag.eq.0)then
       do j=0,jmaxm
        te_m(j)=te_exp(j)
        ti_m(j)=ti_exp(j)
        vphi_m(j)=rmajor_exp*angrotp_exp(j)
        ne_m(j)=ne_exp(j)
        ni_m(j)=ni_exp(j)
        nz_m(j)=nz_exp(j)
        ns_m(j)=nfast_exp(j)
       enddo
      endif
c
      do j=0,kmaxmt
        te_t(0,j)=0.D0
        ti_t(0,j)=0.D0
        ti_exp_t(0,j)=0.D0
        vphi_exp_t(0,j)=0.D0
        te_exp_t(0,j)=0.D0
        vphi_t(0,j)=0.D0
        vexb_t(0,j)=0.D0
        vexb_t(0,j)=0.D0
        vetor_t(0,j)=0.D0
        vepol_t(0,j)=0.D0
        vstar_t(0,j)=0.D0
        egamma_t(0,j)=0.D0
        egamma_vphi_t(0,j)=0.D0
        egamma_vpol_t(0,j)=0.D0
        egamma_vstar_t(0,j)=0.D0
        anratem_t(0,j)=0.D0
        etaphi_t(0,j)=0.D0
        gammap_t(0,j)=0.D0
        time_t(j)=0.D0
      enddo
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c run PTOR or PTOR_DV
c
      if (dvflag .gt. 0) then
        call ptor_dv
c        call derivedmodprofiles
      endif
c
      if(dvflag.eq.0)then
       te_m(0)=te_m(1)*1.001
       ti_m(0)=ti_m(1)*1.001
       ne_m(0)=ne_m(1)*1.001
       ni_m(0)=ni_m(1)*1.001
       nz_m(0)=nz_m(1)*1.001
       ns_m(0)=ns_m(1)*1.001
       te_m(0)=te_m(1)*1.001 
c
       call derivedmodprofiles
      endif
c
      jin_m=jin_hold
      jout_m=jout_hold
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c Printout results
c
c  lprint_cdf > 0 to write output to NetCDF file 'xptor.cdf'
c  lprint_prof > 0 to write profiles to ASCII files 'out' and 'out2'
c  lprint_pulse > 0 to printout full time histories to ASCII files
c  lprint_glf > 0 to printout namelist data for GLF23 stand-alone code
c
      if (i_proc .eq. 0) then
c
c NETCDF printout
c
c      if (lprint_cdf .gt. 0) call write_cdf(lprint_cdf,stdev,
c     &                                      relrms,offset)
      if (lprint_cdf .gt. 0) call write_netcdf
c
c ASCII printout
c
      if (lprint_prof .gt. 0) then
        call write_ascii(lprint_prof,lprint_pulse)
      endif
c
      if (lprint_glf .gt. 0) then
        call write_glf(zdate,j,j1,j2)
      endif
c
      if (lprint_ufile(1).gt.0 .or. lprint_ufile(2).gt.0 .or.
     &   lprint_ufile(3).gt.0) then
        if(i_proc.eq.0) write(*,*) 'Calling write_ufile ...'
        call write_ufile(lprint_ufile,imodel)
      endif
c
      endif
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c      if (lprint_glf .gt. 1) then
c        do j=1,10
c          write(*,350) j,xkyf_k_gf(j),gamma_k_j(j,25),freq_k_j(j,25),
c     &                 chie_k_j(j,25),chii_k_j(j,25),chie_e_k_j(j,25)
c        enddo
c      endif
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
 900  continue
c
c      call MPI_FINALIZE(i_err)
c
c      call mclock(time_cp)
c      cputime = 0.01 * time_cp
c      cputime = 0.01 * mclock()
      endtime=MPI_WTIME()
      if (i_proc .eq. 0)
     &  write(*,500) endtime-starttime
       call MPI_FINALIZE(i_err)
c
 20   format (a)
 35   format(' *** Warning: jout_m greater than grid size ***')
 50   format(a10,2x,a4,a6)
 60   format(/,'*** Increasing grid size ***',/,
     &      '    Grid size = ',i3,', BC zone = ',i3,/)
 70   format(/,1x,'ne-exp(0)=',f7.3,' ni-exp(0)=',f7.3,
     &      ' Te-exp(0)=',f7.3,' Ti-exp(0)=',
     &      f7.3,' vphi-exp(0)=',1pe12.4,/)
 75   format(/,1x,'Te-exp(0)=',f7.3,' Ti-exp(0)=',
     &      f7.3,' vphi-exp(0)=',1pe12.4,/)
 100  format('  Using ITER/TRANSP ufiles ...')
 105  format('  Using ONETWO iterdb file ...')
 106  format('  Using ONETWO netcdf file ...')
 107  format('  Creating test profiles ...')
 109  format('  Using Std test case parameters ...')
 110  format(2x,a10,/a11,/,2x,'time = ',0pf7.4,' secs')
 150  format(2i4)
 151  format(4x,'Te-max (keV)  = ',0pf8.3,6x,
     &       'Ti-max (keV)  = ',0pf8.3)
 152  format(4x,'Te-ped (keV)  = ',0pf8.3,6x,
     &       'Ti-ped (keV)  = ',0pf8.3)
 200  format(i2,2x,0p1f4.2,0p6f10.5)
 240  format(i2,2x,0p1f6.4,0p12f9.4)
 250  format(i2,2x,0p1f6.4,0p6f11.6)
 300  format(i4,2x,0p51f10.5)
 305  format(i4,2x,0p51f14.5)
 310  format(0pe16.7,2x,0p301e20.11)
 350  format(i4,2x,0pf10.5,1p6e12.4)
 500  format(/,'  ** CPU time (secs) = ',f8.3,' **')
c
      end
