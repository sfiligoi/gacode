 
  MODULE nbnamelist
 !--------------------------------------------------------------------------
 ! THIS NAMELIST MODULE IS USED FOR NUBEAM input prior 
 ! to the new nubeam version 201112
 ! For nubeam versions 201112 and P_Nfreya a separate namelist module
 ! is used because the post and prior namelists for nubeam are
 ! too different.
 ! the decision as to which module will be used is made through
 ! the inone input switch nubeam_version (which is defaulted to the
 ! older (2003-2006) version of nubeam
 !------------------------------------------------------------HSJ-2/14//2012
       USE nrtype,ONLY : Dp,I4B,SP
       USE param, ONLY : kj,kprim,kimp,kb
#ifdef ONETWO
       USE transp,ONLY : beam_data,nubeam_nclass,nlfhe3,nlfhe4,   &
                         nlfst,nlfsp, nlusf3,xdatsfa, xdatsf3,    &
                         xdatsft,xdatsfp,nlusfa,nlusfp,nlusft,    &
                         plfhe3,plfhe4,plfst,plfsp,dxbsmoo,wghta, &
                         nptcls,nptclf,ndep0,goocon,dtn_orbit,    &
                         nsigexc,nlminsv,nlbbcx,nlebei,dn0out,    &
                         nmsigx,xcfanbi,xdfanbi,xefanbi,xcfafus,  &
                         xdfafus,xefafus,nlbcde,nlbcoh,nlbcpa,    &
                         nlorbo,nlbflr,nlfbmflr,nlbgflr,gflr_min, &
                         nkdifb,fdifbe,edifbe,plhgt,nlcprb,       &
                         nmcurb,xdepmod,xp_nt1,inzri,             &
                         ndifbep,nubeam_restart,                  &
                         nubeam_dt
#else
!! gfortran wotn accept namelist variabless in use statments !!!
!       USE P_nfreya_interface,   ONLY :                           &
!                         beam_data,nubeam_nclass,nlfhe3,nlfhe4,   &
!                         nlfst,nlfsp, nlusf3,xdatsfa, xdatsf3,    &
!                         xdatsft,xdatsfp,nlusfa,nlusfp,nlusft,    &
!                         plfhe3,plfhe4,plfst,plfsp,dxbsmoo,wghta, &
!                         nptcls,nptclf,ndep0,goocon,dtn_orbit,    &
!                         nsigexc,nlminsv,nlbbcx,nlebei,dn0out,    &
!                         nmsigx,xcfanbi,xdfanbi,xefanbi,xcfafus,  &
!                         xdfafus,xefafus,nlbcde,nlbcoh,nlbcpa,    &
!                         nlorbo,nlbflr,nlfbmflr,nlbgflr,gflr_min, &
!                         nkdifb,fdifbe,edifbe,plhgt,nlcprb,       &
!                         nmcurb,xdepmod,xp_nt1,inzri,             &
!                         ndifbep,nubeam_restart,                  &
!                         nubeam_dt,nubeam0_dt,d_fast_ion
!! so use the following instead:
       USE P_nfreya_interface,   ONLY :                           &
                         beam_data,nubeam_nclass,nlfhe3a=>nlfhe3, &
                         nlfhe4a=>nlfhe4,                         &
                         nlfsta=>nlfst,nlfspa=>nlfsp,             &
                         nlusf3a=>nlusf3,xdatsfa, xdatsf3,                 &
                         xdatsft,xdatsfp,nlusfa,nlusfp,nlusft,    &
                         plfhe3,plfhe4,plfst,plfsp,dxbsmoo,wghta, &
                         nptcls,nptclf,ndep0,goocon,dtn_orbit,    &
                         nsigexc,nlminsv,nlbbcx,nlebei,dn0out,    &
                         nmsigx,xcfanbi,xdfanbi,xefanbi,xcfafus,  &
                         xdfafus,xefafus,nlbcde,nlbcoh,nlbcpa,    &
                         nlorbo,nlbflr,nlfbmflr,nlbgflr,gflr_min, &
                         nkdifb,fdifbe,edifbe,plhgt,nlcprb,       &
                         nmcurb,xdepmod,xp_nt1,inzri,             &
                         ndifbep,nubeam_restart,                  &
                         nubeam_dt,nubeam0_dt,d_fast_ion

#endif 

!
       USE numbrs,    ONLY   : nj,nprim,nimpc=>nimp
!
!
       IMPLICIT NONE
       SAVE
       INTEGER kj_check
!       INTEGER,PARAMETER :: nbeamx = 10
       INTEGER,PARAMETER  :: nbeamx = kb ! HSJ 8/10/11
       INTEGER,PARAMETER :: nrhixm = 4   !set to beam_data%nrhix ??
       INTEGER,PARAMETER :: nfracionsx = 1 ! dummy
       INTEGER,PARAMETER :: nfracimpsx = 1 !dummy 
!
!      nfracionsx,nfracimpsx and frac_ions, frac_imp are not
!      used to generate profiles in the Onetwo version.
!      thes e quantities are retained only so that the
!      namelist reads will work without having to modify the
!      namelists. See ...../nbdrive/nbdrive_prof.f90
!      on how these quantities are eliminated 
!
       INTEGER,PARAMETER :: nneutx = 10
       INTEGER,PARAMETER :: ngmaxx = 10
       REAL(DP) tbona(nbeamx),tboffa(nbeamx) !required for sub nblist_read only
!
!
  ! the PYTHON generated namelist values have a one-to-one correspondence
  ! with elements of NUBEAM input data structures -- units and meaning are
  ! documented in nubeam/nbspec.dat (** not repeated here **).
!
  ! only a subset of the NUBEAM controls are exposed in this NUBEAM
  ! test driver namelist.  This is an effort at simplification.  The
  ! omitted controls correspond to specialized "code developer"
  ! adjustments, or, to "advanced" features which (for the time
  ! being at least) have been omitted from the test driver.
!
  ! inputs needed specifically by the driver code (but not directly by
  ! NUBEAM itself), without this one-to-one correspondence, are declared
  ! and described in the "hand maintained section", below...
!
!
 !==========================================================================
  !-->PYTHON_START
  ! (nbigen.py: python generated code, do not edit)
 
  ! --> nubeam/nbspec.dat block: sys
    INTEGER :: nseed
    INTEGER :: nclass
  !****NAMELIST/nbdrive_naml/ nseed,nclass
 
  ! --> nubeam/nbspec.dat block: times
 
  ! --> nubeam/nbspec.dat block: grid
    INTEGER :: nzones
    INTEGER :: nznbma
    INTEGER :: nznbme
    REAL(DP) :: ebdmax
  !****NAMELIST/nbdrive_naml/ nzones,nznbma,nznbme,ebdmax
 
  ! --> nubeam/nbspec.dat block: beams
    INTEGER :: ngmax
    REAL(DP) :: backz(ngmaxx)  ! (1:ngmax)
    REAL(DP) :: aplasm(ngmaxx)  ! (1:ngmax)
    INTEGER :: nbeam
    REAL(DP) :: xzbeama(nbeamx)  ! (1:nbeam)
  !****NAMELIST/nbdrive_naml/ ngmax,backz,aplasm,nbeam,xzbeama
    REAL(DP) :: abeama(nbeamx)  ! (1:nbeam)
    LOGICAL :: nlco(nbeamx)  ! (1:nbeam)
    INTEGER :: ntrace(nbeamx)  ! (1:nbeam)
    REAL(DP) :: rtcena(nbeamx)  ! (1:nbeam)
    REAL(DP) :: xlbtna(nbeamx)  ! (1:nbeam)
  !****NAMELIST/nbdrive_naml/ abeama,nlco,ntrace,rtcena,xlbtna
    REAL(DP) :: xybsca(nbeamx)  ! (1:nbeam)
    REAL(DP) :: xbzeta(nbeamx)  ! (1:nbeam)
    INTEGER :: nbshapa(nbeamx)  ! (1:nbeam)
    REAL(DP) :: bmwidra(nbeamx)  ! (1:nbeam)
    REAL(DP) :: bmwidza(nbeamx)  ! (1:nbeam)
  !****NAMELIST/nbdrive_naml/ xybsca,xbzeta,nbshapa,bmwidra,bmwidza
    REAL(DP) :: divra(nbeamx)  ! (1:nbeam)
    REAL(DP) :: divza(nbeamx)  ! (1:nbeam)
    REAL(DP) :: foclra(nbeamx)  ! (1:nbeam)
    REAL(DP) :: foclza(nbeamx)  ! (1:nbeam)
    INTEGER :: nbapsha(nbeamx)  ! (1:nbeam)
  !****NAMELIST/nbdrive_naml/ divra,divza,foclra,foclza,nbapsha
    REAL(DP) :: rapedga(nbeamx)  ! (1:nbeam)
    REAL(DP) :: xzpedga(nbeamx)  ! (1:nbeam)
    REAL(DP) :: xlbapa(nbeamx)  ! (1:nbeam)
    REAL(DP) :: xybapa(nbeamx)  ! (1:nbeam)
    INTEGER :: nbapsh2(nbeamx)  ! (1:nbeam)
  !****NAMELIST/nbdrive_naml/ rapedga,xzpedga,xlbapa,xybapa,nbapsh2


    REAL(DP) :: rapedg2(nbeamx)  ! (1:nbeam)
    REAL(DP) :: xzpedg2(nbeamx)  ! (1:nbeam)
    REAL(DP) :: xlbapa2(nbeamx)  ! (1:nbeam)
  !****NAMELIST/nbdrive_naml/ rapedg2,xzpedg2,xlbapa2
 
  ! --> nubeam/nbspec.dat block: impurity
    INTEGER :: nrhix
    REAL(DP) :: xzimpx(nrhixm)  ! (1:nrhix)
    REAL(DP) :: aimpx(nrhixm)  ! (1:nrhix)
  !****NAMELIST/nbdrive_naml/ nrhix,xzimpx,aimpx
 
  ! --> nubeam/nbspec.dat block: minority
 
  ! --> nubeam/nbspec.dat block: powers
    REAL(DP) :: einja(nbeamx)  ! (1:nbeam)
    REAL(DP) :: pinja(nbeamx)  ! (1:nbeam)
    REAL(DP) :: ffulla(nbeamx)  ! (1:nbeam)
    REAL(DP) :: fhalfa(nbeamx)  ! (1:nbeam)
  !****NAMELIST/nbdrive_naml/ einja,pinja,ffulla,fhalfa
 
  ! --> nubeam/nbspec.dat block: fusion (from transp.mod )
   ! activate these for gfortran 2/25/13 HSJ
  ! NOTE gfortran is used for P_Nfreya on star cluster
#ifndef ONETWO 
    ! if compiling for P_Nfreya then need these
    ! for Onetwo they are defined in transp module: HSJ
    LOGICAL :: nlfhe3
    LOGICAL :: nlfhe4
    LOGICAL :: nlfst
    LOGICAL :: nlfsp
    LOGICAL :: nlusf3
#endif
  !****NAMELIST/nbdrive_naml/ nlfhe3,nlfhe4,nlfst,nlfsp,nlusf3
!    REAL(DP) :: xdatsf3
!    LOGICAL :: nlusfa
!    REAL(DP) :: xdatsfa
!   LOGICAL :: nlusft
!    REAL(DP) :: xdatsft
  !****NAMELIST/nbdrive_naml/ xdatsf3,nlusfa,xdatsfa,nlusft,xdatsft
!    LOGICAL :: nlusfp
!    REAL(DP) :: xdatsfp
!    REAL(DP) :: plfhe3
!    REAL(DP) :: plfhe4
!    REAL(DP) :: plfst
  !****NAMELIST/nbdrive_naml/ nlusfp,xdatsfp,plfhe3,plfhe4,plfst
!    REAL(DP) :: plfsp
    LOGICAL :: nlfatom
  !****NAMELIST/nbdrive_naml/ plfsp,nlfatom
 
  ! --> nubeam/nbspec.dat block: fpp
 
  ! --> nubeam/nbspec.dat block: num
  !  REAL(DP) :: dxbsmoo
!    INTEGER :: nptcls
!    INTEGER :: nptclf
!    INTEGER :: ndep0
    INTEGER :: nbbcal
  !****NAMELIST/nbdrive_naml/ dxbsmoo,nptcls,nptclf,ndep0,nbbcal
!    REAL(DP) :: wghta
!    REAL(DP) :: goocon       
    REAL(DP) :: dtn
  !****NAMELIST/nbdrive_naml/ wghta,goocon,dtn
 
  ! --> nubeam/nbspec.dat block: atomic
!    REAL(DP) :: xdepmod
!    INTEGER :: nsigexc
!    LOGICAL :: nlminsv
!    LOGICAL :: nlbbcx
!    LOGICAL :: nlebei
  !****NAMELIST/nbdrive_naml/ xdepmod,nsigexc,nlminsv,nlbbcx,nlebei
!    REAL(DP) :: dn0out
!    INTEGER :: nmsigx
    LOGICAL :: nlhvion
  !****NAMELIST/nbdrive_naml/ dn0out,nmsigx,nlhvion
 
  ! --> nubeam/nbspec.dat block: collid
!    REAL(DP) :: xcfanbi
!    REAL(DP) :: xdfanbi
!    REAL(DP) :: xefanbi
!    REAL(DP) :: xcfafus
!    REAL(DP) :: xdfafus
  !****NAMELIST/nbdrive_naml/ xcfanbi,xdfanbi,xefanbi,xcfafus,xdfafus
!    REAL(DP) :: xefafus
!    LOGICAL :: nlbcde
!    LOGICAL :: nlbcoh
!    LOGICAL :: nlbcpa
!    LOGICAL :: nlorbo
  !****NAMELIST/nbdrive_naml/ xefafus,nlbcde,nlbcoh,nlbcpa,nlorbo
 
  ! --> nubeam/nbspec.dat block: flr
!    LOGICAL :: nlbflr
!    LOGICAL :: nlfbmflr
!    LOGICAL :: nlbgflr
!    REAL(DP) :: gflr_min
  !****NAMELIST/nbdrive_naml/ nlbflr,nlfbmflr,nlbgflr,gflr_min
 
  ! --> nubeam/nbspec.dat block: saw
 
  ! --> nubeam/nbspec.dat block: adif
!    INTEGER :: nkdifb
     INTEGER :: ndifbe
  !****NAMELIST/nbdrive_naml/ nkdifb,ndifbe,fdifbe,edifbe
 
  ! --> nubeam/nbspec.dat block: ripple
 
  ! --> nubeam/nbspec.dat block: outcon
 
  ! --> nubeam/nbspec.dat block: fishbone
 
  ! --> nubeam/nbspec.dat block: box
 
  ! --> nubeam/nbspec.dat block: tube
 
  ! --> nubeam/nbspec.dat block: misc
!    REAL(DP) :: plhgt
!    LOGICAL :: nlcprb
!    INTEGER :: nmcurb
  !****NAMELIST/nbdrive_naml/ plhgt,nlcprb,nmcurb
 
  ! --> nubeam/nbspec.dat block: profiles
  !-->PYTHON_END
  !==========================================================================
  !-------------------> HAND MAINTAINED section <---------------------
!
  REAL(DP) time_start,time_stop   ! time range of simulation (sec)
  INTEGER :: nsteps             ! no. of steps to cover time range
!
  INTEGER :: nstep_start        ! restart at step (if .gt. 1)
  INTEGER :: nstep_stop         ! save restart data and stop before step
  !                               (if .gt. 1)
!
  LOGICAL cleanup               ! .TRUE. to clean up (delete) state files on 
  !                               successful completion (this is the default).
!
  INTEGER :: ntheta             ! no. of theta points for bicubic splines
  !                               (default value should be fine)
!
  REAL(DP) Zeff(kj)               ! kj must match kj in nbdrive routines
                                ! (dimensionless)
  !  Zeff equation:
  !                               ne*Zeff = sum{i}[n(i)*Z(i)^2]
  !                                 sum to include fast ion species
  !  quasineutrality equation:
  !                               ne = sum{i}[n(i)*Z(i)]
!
  !  *** ion species information ***
  !    the Zs and As of main plasma and impurity ions are specified 
  !    in the Python-generated section -- because there is a one-to-one 
  !    correspondence with actual NUBEAM input data structures for these
  !    items.  However, the relative concentrations of the thermal ion
  !    species (for Zeff & quasineutrality equations) must be specified
  !    as well.  The default is equal concentrations for all listed
  !    species; the fractions are dimenless but if the sum of the
  !    fractions is not one, nbdrive will renormalize.
!
  REAL (DP) te_kj(kj),ti_kj(kj),ne_kj(kj), omega_kj(kj)
  REAL (DP) ni_kj(kj,kprim),nimp_kj(kj,kimp),vloop_kj(kj)
  REAL (DP) n0w_kj(kj,2 ),n0v_kj(kj,2 )
  REAL(DP) frac_ions(nfracionsx)         ! actual number of active species << 100
                                ! (dimensionless)
!
  REAL(DP) frac_imps(nfracimpsx)         ! actual number of active species << 100
                                ! (dimensionless)
!
  !  *** profile information ***
!
  INTEGER :: nrho               ! no. of radial zone bdys -- nbdrive profiles
                                ! need not match NUBEAM nzones 
                                !   ==> equispaced radial grid x over [0,1]
                                !       with  h = 1/(nrho-1)
!
  !  for each item f, the values f0,fa,fxpin,fxpout are specified; 
  !    f(x) = f_a + (f_0-f_a)*(1-x**f_xpin)**f_xpout
!
  !  default values -- f_xpin=2, f_xpout=1, f_0=0, f_a=0
  !    f_0 and f_a in units of f; f_xpin & f_xpout are dimensionless.
  !
  !    generally, "real" values need to be provided at least for f_0 and
  !    f_a; for SOME profiles e.g. rotation, all-zero values are allowed.
!
  !  special treatment for wall-neutral-density profiles:
!
  !   n0w(x) = exp[ ln(n0w_a) + 
  !                   (ln(n0w_0)-ln(n0w_a))*(1-x**n0w_xpin)**n0w_xpout ]
!
  !  note that each specie has its own n0w_0 (central wall source neutral 
  !    density) and n0w_a (edge wall source neutral density), but that
  !    n0w_xpin and n0w_xpout are the same for *all* species.
!
  !  there is also a set of controls for volume-source-neutral-density
  !  profiles, which are fit with the usual formula (no ln, no exp).
  !    in real simulations, volume-source profiles are driven by the
  !  neutral beam charge-exchange internal neutral source, but for
  !  nbdrive, which has no neutral gas model, these must be input
  !  the defaults, zero volume source neutral density, can be used.
!
  !  simplifying assumptions:
  !    all ion species have the same temperature and rotation velocity
  !      (true of all NUBEAM runs).
  !    neutral specie temperatures = ion temperature 
  !      (not usually the case when NUBEAM is used with a transport code).
  !    neutral specie rotation velocities = ion rotation velocity
  !      (not usually the case when NUBEAM is used with a transport code).
!
  !  ** electron temperature (eV) **
!
  REAL(DP) te_0,te_a,te_xpin,te_xpout
!
  !  ** ion temperature (eV) **
!
  REAL(DP) ti_0,ti_a,ti_xpin,ti_xpout
!
  !  ** electron density (n/cm**3) **
!
  REAL(DP) ne_0,ne_a,ne_xpin,ne_xpout
!
  !  ** angular velocity (radians/second) **
  !     (local Vphi = R*omega(x))
!
  REAL(DP) omeg_0,omeg_a,omeg_xpin,omeg_xpout
!
  !  ** radial electrostatic potential (volts) **
!
  REAL(DP) epot_0,epot_a,epot_xpin,epot_xpout
!
  !  ** toroidal loop voltage (volts) **
!
  REAL(DP) vloop_0,vloop_a,vloop_xpin,vloop_xpout
!
  !  ** anomolous diffusivity (cm**2/sec) **
!
  REAL(DP) adiff_0,adiff_a,adiff_xpin,adiff_xpout
!
  !  ** wall-source-neutral densities  (n/cm**3)
  !     separate value for each main plasma ion species
  !       ** note exp / ln transformation in fit, described above **
!
  REAL(DP) n0w_0(nneutx),n0w_a(nneutx),n0w_xpin,n0w_xpout
!
  !  ** volume-source-neutral densities  (n/cm**3)
  !     separate value for each main plasma ion species
  !       (standard analytic fit formula used)
!
  REAL(DP) n0v_0(nneutx),n0v_a(nneutx),n0v_xpin,n0v_xpout
!
  !--------------------------
  !  radial resolution of fast ion distribution function
  !  spatial grid; no. of radial rows in grid extrapolation
  
  INTEGER nzone_fb  ! basic grid
!
  !  at present nzone_fb = 10 or nzone_fb = 20 are allowed, and, 
  !  nzones = [integer]*nzone_fb is enforced.  The grid used for
  !  computing distribution functions f(rho,theta,E,vpll/v) is
  !  irregular, with fewer poloidal angle zones in zone rows near
  !  the axis, and more out towards the edge; the cross-sectional
  !  areas of each 2d grid zone are similar.
  ! 
  !  For each fast ion specie, the size of the distribution function
  !  is O(nzone_fb**2*nznbme*nznbma); for updown symmetric equilibria
  !  only the top half of the poloidal variation is retained.
  ! 
  !  Number of poloidal angle zones in each radial zone row:
  !                nzone_fb=20                 nzone_fb=10
  !  updown:  symmetric   Asymmetric      symmetric    Asymmetric
  !  radial
  !  zone row:
  !     1         2            4              2            4
  !     2         4            8              4            8       
  !     3         6           12              6           12       
  !     4         8           16              8           16       
  !     5        10           20             10           20       
  !     6        12           24             12           24       
  !     7        14           28             14           28       
  !     8        16           32             16           32       
  !     9        18           36             18           36       
  !    10        20           40             20           40       
  !    11        22           44
  !    12        24           48
  !    13        26           52
  !    14        28           56
  !    15        30           60
  !    16        32           64
  !    17        34           68
  !    18        36           72
  !    19        38           76
  !    20        40           80
  ! total no.
  ! of spatial
  ! zones:      420          840            110          220
  !
  ! for default settings of the pitch & energy grid sizes, there are
  ! 5000 v-space zones in each spatial zone.
  !--------------------------
  !--------------------------
  !  MHD equilibrium selection
  !    an EFIT or TRANSP data source must be specified ***
  !
   CHARACTER*100 mhdpath
  !
  !  example *values* for mhdpath variable:
  !  EFIT data from G-eqdsk file:
  !    "EFIT:/u/user/efit_data/123456/geqdsk_8888.dat"
  !    "EFIT:geqdsk_8888.dat"     (file must be found in $cwd)
  !                               (relative paths can also be used)
  !  EFIT data from MDSplus:
  !    "EFIT:MDS+:EUROPA.PPPL.GOV:8501:EFIT01(103875;t=.25)"
  !       or in general
  !          EFIT:MDS+:<server-name[:port]>:<tree-name>(<shot-no>;t=<time>)
  !
  !  TRANSP data from TRANSP output files (usually <runid>.CDF):
  !    "TRANSP:/transp/results/TFTR.88/37065Z19"
  !    "TRANSP:37065Z19"          (37065Z19.CDF must be found in $cwd)
  !                               (relative paths can also be used)
  !  TRANSP data from MDSplus:
  !    "TRANSP:MDS+:CMODA.PSFC.MIT.EDU:TRANSP03(12345678)"
  !            (MIT-style TRANSP run identification)
  !    "TRANSP:MDS+:BIRCH.PPPL.GOV:TRANSP(TFTR.88,37065Y82)"
  !            (PPPL-style TRANSP run identification)
  !       or in general
  !          TRANSP:MDS+:<server-name[:port]>:<tree-name>(ident)
  !             where (ident) = (<shot-no>) for MIT-style
  !               and (ident) = (<tok.yy>,<runid>) for PPPL-style
  !                               (with <tree-name> = TRANSP in most cases)
  !
  !  if TRANSP is used, the time and +/- delta_time (for averaging) must
  !  also be given:
!
  REAL(DP) tr_time           ! TRANSP time, seconds
  REAL(DP) tr_delta_t        ! TRANSP time, seconds (default=0.0, no averaging).

!
  !  for TRANSP equilibrium only:  CCW = "counter-clockwise"
  !  set value = -1 for "clockwise" orientation of toroidal current or field.
!
  INTEGER nsnccwi_tr       ! =1: Iphi is CCW looking down on tokamak from above
  INTEGER nsnccwb_tr       ! =1: Bphi is CCW looking down on tokamak from above
!
  !  the following specifies the distance beyond the plasma to which orbits
  !  and FLR can extend without a fast ion orbit being considered "lost".
!
  REAL(DP) tr_sol_width      ! width (cm) of scrape-off layer
!
  !  comparison to EFIT:  EFIT (G-eqdsk) conventions define the field and
  !  current orientation; also, EFIT includes a piecewise linear limiter
  !  location specification, which defines the vacuum region beyond the
  !  plasma boundary.  (TRANSP also has a detailed description, but, it
  !  is not so easy to get at, and so it is omitted here for simplitity's
  !  sake).
!
  !  no. of points in (R,Z) grid, space into which vacuum field and 
  !  magnetic coordinate system will be extrapolated, beyond the plasma
  !  boundary in all directions:
!
  INTEGER nRgrid_tr,nZgrid_tr  ! (R,Z) grid sizes -- default: 100x100 OK.
  
  !  if d = tr_sol_width, e = a small additional distance, then:
  !  the (R,Z) grid will cover [Rmin-d-e,Rmax+d+e]x[Zmin-d-e,Zmax+d+e],
  !  where Rmin,Rmax, Zmin,Zmax are taken from the plasma boundary contour.
!
  !--------------------------
  
  ! added to this code by hsj 5/9/11 nlbdat 
  ! this appears to be a new variable in namelist which is written by
  ! C. Greenfield's IDL code )
   LOGICAL nlbdat
    REAL(DP) :: xrapoffa(nbeamx)  ! (1:nbeam)   ! added HSJ 8/9/11
    REAL(DP) :: xzapoffa(nbeamx)  ! (1:nbeam)   ! added HSJ 8/9/11
    REAL(DP) :: xrapoff2(nbeamx)  ! (1:nbeam)   ! added HSJ 8/9/11
    REAL(DP) :: xzapoff2(nbeamx)  ! (1:nbeam)   ! added HSJ 8/9/11
!
!
  !****NAMELIST/nbdrive_naml/                           &
!       time_start, time_stop, nsteps, nstep_start, &
!       nstep_stop, ntheta, Zeff, frac_ions,        &
!       frac_imps, nrho,                            &
!       te_0,te_a,te_xpin,te_xpout,                 &
!       ti_0,ti_a,ti_xpin,ti_xpout,                 &
!       ne_0,ne_a,ne_xpin,ne_xpout,                 &
!       te_kj,ti_kj,ne_kj,omega_kj,ni_kj,nimp_kj,   &
!       vloop_kj,n0w_kj,n0v_kj,                     &
!       omeg_0,omeg_a,omeg_xpin,omeg_xpout,         &
!       epot_0,epot_a,epot_xpin,epot_xpout,         &
!       vloop_0,vloop_a,vloop_xpin,vloop_xpout,     &
!       adiff_0,adiff_a,adiff_xpin,adiff_xpout,     &
!       n0w_0,n0w_a,n0w_xpin,n0w_xpout,             &
!       n0v_0,n0v_a,n0v_xpin,n0v_xpout,             &
!       mhdpath,tr_time,tr_delta_t,                 &
!       nsnccwi_tr,nsnccwb_tr,tr_sol_width,         &
!       nRgrid_tr,nZgrid_tr,                        &
!       nzone_fb,kj_check,nj,tbona,tboffa,          &
!       nubeam_dt,                                  &
!      xrapoffa,xzapoffa,xrapoff2,xzapoff2,nlbdat

! xrapoffa,xzapoffa,xrapoff2,xzapoff2,nlbdat
! temporarily  removed HSJ
! because nfreya  doesnt recognize these values
!(need new nfreya)
! need to add back in to use p_Nfreya ! 888888999999


!------------------------------------------------------------------------------
! the complete alphabetized namelist is given here :
! it is just the above added all together
!------------------------------------------------------------HSJ-2-14-2012-----
  NAMELIST/nbdrive_naml/                                                                      &    
      abeama          ,      adiff_0         ,      adiff_a         ,      adiff_xpin      ,  &
      adiff_xpout     ,      aimpx           ,      aplasm          ,      backz           ,  &
      bmwidra         ,      bmwidza         ,      divra           ,      divza           ,  &
      dn0out          ,      dtn             ,      dxbsmoo         ,      ebdmax          ,  &
      edifbe          ,      einja           ,      epot_0          ,      epot_a          ,  &
      epot_xpin       ,      epot_xpout      ,      fdifbe          ,      ffulla          ,  &
      fhalfa          ,      foclra          ,      foclza          ,      frac_imps       ,  &
      frac_ions       ,      gflr_min        ,      goocon          ,      kj_check        ,  &
      mhdpath         ,      n0v_0           ,      n0v_a           ,      n0v_kj          ,  &
      n0v_xpin        ,      n0v_xpout       ,      n0w_0           ,      n0w_a           ,  &
      n0w_kj          ,      n0w_xpin        ,      n0w_xpout       ,      nbapsh2         ,  &
      nbapsha         ,      nbbcal          ,      nbeam           ,      nbshapa         ,  &
      nclass          ,      ndep0           ,      ndifbe          ,      ne_0            ,  &
      ne_a            ,      ne_kj           ,      ne_xpin         ,      ne_xpout        ,  &
      ngmax           ,      nimp_kj         ,      ni_kj           ,      nj              ,  &
      nkdifb          ,      nlbbcx          ,      nlbcde          ,      nlbcoh          ,  &
      nlbcpa          ,      nlbflr          ,      nlbgflr         ,      nlco            ,  &
      nlcprb          ,      nlebei          ,      nlfatom         ,      nlfbmflr        ,  &
      nlfhe3          ,      nlfhe4          ,      nlfsp           ,      nlfst           ,  &
      nlhvion         ,      nlminsv         ,      nlorbo          ,      nlusf3          ,  &
      nlusfa          ,      nlusfp          ,      nlusft          ,      nmcurb          ,  &
      nmsigx          ,      nptclf          ,      nptcls          ,      nRgrid_tr       ,  &
      nrho            ,      nseed           ,      nsigexc         ,      nsnccwb_tr      ,  &
      nsnccwi_tr      ,      nsteps          ,      nstep_start     ,      nstep_stop      ,  &
      ntheta          ,      ntrace          ,      nubeam_dt       ,      nZgrid_tr       ,  &
      nznbma          ,      nznbme          ,      nzones          ,      nzone_fb        ,  &
      omega_kj        ,      omeg_0          ,      omeg_a          ,      omeg_xpin       ,  &
      omeg_xpout      ,      pinja           ,      plfhe3          ,      plfhe4          ,  &
      plfsp           ,      plfst           ,      plhgt           ,      rapedg2         ,  &
      rapedga         ,      nrhix           ,      rtcena          ,      tboffa          ,  &
      tbona           ,      te_0            ,      te_a            ,      te_kj           ,  &
      te_xpin         ,      te_xpout        ,      time_start      ,      time_stop       ,  &
      ti_0            ,      ti_a            ,      ti_kj           ,      ti_xpin         ,  &
      ti_xpout        ,      tr_delta_t      ,      tr_sol_width    ,      tr_time         ,  &
      vloop_0         ,      vloop_a         ,      vloop_kj        ,      vloop_xpin      ,  &
      vloop_xpout     ,      wghta           ,      xbzeta          ,      xcfafus         ,  &
      xcfanbi         ,      xdatsf3         ,      xdatsfa         ,      xdatsfp         ,  &
      xdatsft         ,      xdepmod         ,      xdfafus         ,      xdfanbi         ,  &
      xefafus         ,      xefanbi         ,      xlbapa          ,      xlbapa2         ,  &
      xlbtna          ,      xybapa          ,      xybsca          ,      xzbeama         ,  &
      xzimpx          ,      xzpedg2         ,      xzpedga         ,      Zeff



  END MODULE nbnamelist
