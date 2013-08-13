SUBROUTINE set_nbinputs

  USE nbi_types
  USE transp
  USE io,ONLY :           ncrt
  USE nbi_dimensions
  IMPLICIT NONE


  !---------------------
  ! initialize inputs to NBI calculation, 
  ! these are mostly items that get set ONCE and can not be reset; 
  ! for updatable items see "set_nbupdate".
  ! 
  ! method:  (a) call nbi_init_xxx to set defaults for block xxx
  !          (b) explicitly set xxx elements from transp module
  !          (c) call nbi_set_xxx to copy into MC code
  !              note this is only done ONCE at the beginning of a
  !              simulation; subsequent restarts use a state file.
  !
  ! note to non-TRANSP folks who might look at this as a template for
  ! interfacing to the fast ion code: the reason the type element names
  ! and the TRCOM names match so often, is that fast ion code used to be
  ! part of TRANSP and used to get its data from TRCOM directly...
  !
  ! summary comments are provided here on individual elements.  See also
  ! "nbspec.dat" or "nbi_com_mod.f90" for more detailed comments

  TYPE(nbitype_sys)    :: zsys         ! basic system
  TYPE(nbitype_grid)   :: zgrid        ! computational grids
  TYPE(nbitype_beams)  :: zbeams       ! list of neutral beams
  TYPE(nbitype_impurity) :: zimpurity  ! impurity info
  TYPE(nbitype_minority) :: zminority  ! minority species
  TYPE(nbitype_fusion) :: zfusion      ! fusion product species
  TYPE(nbitype_fpp)    :: zfpp         ! fpp code interface
  TYPE(nbitype_num)    :: znum         ! numerical controls 
  TYPE(nbitype_atomic) :: zatomic      ! atomic physics options
  TYPE(nbitype_collid) :: zcollid      ! collision operator adjustments
  TYPE(nbitype_flr)    :: zflr         ! finite larmor radius model
  TYPE(nbitype_saw)    :: zsaw         ! sawtooth controls
  TYPE(nbitype_adif)   :: zadif        ! anomolous diffusion
  TYPE(nbitype_ripple) :: zripple      ! TF rippple model controls
  TYPE(nbitype_fishbone) :: zfbone     ! "fishbone" model controls
  TYPE(nbitype_outcon) :: zoutcon      ! output control options
  TYPE(nbitype_box)    :: zbox         ! beam in box controls
  TYPE(nbitype_misc)   :: zmisc        ! miscellany

  !--------------------------------------------------------------

  INTEGER ierr,istat,nrhix

  INTEGER i,iptcls_max

  INTEGER id_ripple,id_mcgrid, nbeam_tr

  LOGICAL fpp,sawtooth,anom_diff,ripple,outcon ,fishbone
  LOGICAL beam_box
  fpp = .FALSE. ; sawtooth = .FALSE. ; anom_diff = .FALSE.
  ripple = .FALSE. ; outcon =.false. ; fishbone = .false.
  beam_box = .false.
  !--------------------------------------------------------------
  !  basic system

  CALL nbi_init_sys(zsys)     ! fetch defaults

  zsys%runid = runid          ! runid string
  zsys%nonlin = ncrt          ! LUN for messages
  zsys%nseed = beam_data%nseed          ! random number seed
  zsys%nlbout = nlbout        ! switch for "lost orbits" file
  zsys%lunnbx = lunnbx        ! LUN for lost orbits file
  zsys%lunres = lunres   ! LUN for lost orbits file backup (in case of restart)
  zsys%mrstrt = mrstrt   ! restart flag (set to zero to suppress restart capability)
  zsys%workpath = ' '         ! path to directory for debug state files
  ! (blank means, use current working directory)
  zsys%nclass= nubeam_nclass         ! activate NCLASS beamline-resolved outputs
  CALL nbi_set_sys(zsys,ierr) ! apply...
  IF(ierr.NE.0) THEN
     CALL errmsg_exit(' ?set_nbinputs:  nubeam/nbi_set_sys error code set.')
  ENDIF


  PRINT *,'in sub set_nbinputs ,done with nbi_set_sys'


  !--------------------------------------------------------------
  !  grid information

  CALL nbi_init_grid(zgrid)   ! defaults
  !set value for id_mcgrid :
  CALL xpload_mcgrid(id_mcgrid)

  zgrid%id_mcgrid = id_mcgrid ! xplasma id for 2d MC irregular grid
  zgrid%nzones = beam_data%nzone_nb     ! number of radial zones to represent plasma
  zgrid%nznbma = beam_data%nznbma       ! number of pitch angle zones for fbm
  zgrid%nznbme = beam_data%nznbme       ! number of energy zones (beams only) for fbm
  zgrid%ebdmax = beam_data%ebdmax       ! max energy (beams only) for fbm

  ! fbm: fast ion distribution function output information

  ! (note: the energy grid for fusion products is set to the maximum
  !  number of zones and Emax based on the fusion product birth energy)
  ! (set_nbdims sets the nubeam dimensions)

  CALL nbi_set_grid(zgrid,ierr)  ! apply...
  IF(ierr.NE.0) THEN
     CALL errmsg_exit(' ?set_nbinputs:  nubeam/nbi_set_grid error code set.')
     CALL STOP('set_nbinputs ',1)
  ENDIF
  PRINT *,'in sub set_nbinputs ,done  with nbi_set_grid'




  !  the thermal plasma species and all the neutral beams...
  !  see "nbspec.dat" for details...
  !
  !  generally: if a given isotope is injected by neutral beams, 
  !  the corresponding thermal gas must also be present in the model:
  !  this is what the beam fast ions thermalize to.
  ! 
  !  note: "gases" refers to "non-impurity" species, generally isotopes
  !  of H or He... there is a separate block, below, for impurities.
  ! aplasm,backz are set in  set_thermal_species   HSj


!  call set_thermal_species     !pass info through transp module
  

  CALL nbi_init_beams(zbeams)          ! defaults,source nbi_init.f90


  !  thermal plasma species 

  zbeams%ngmax =  beam_data%ngmax         ! no. of gases (eg thermal  species)
  zbeams%backz(1:beam_data%ngmax) = beam_data%backz(1:beam_data%ngmax)   ! Z (e.g. 1 for H, 2 for He)
  zbeams%aplasm(1:beam_data%ngmax) = beam_data%aplasm(1:beam_data%ngmax) ! A (AMU) (e.g. 3 for Tritium)


  !  neutral beams
  nbeam_tr =  beam_data%nbeam
  zbeams%nbeam = nbeam_tr   

  zbeams%xzbeama(1:nbeam_tr) = beam_data%xzbeama(1:nbeam_tr)   ! Z
  zbeams%abeama(1:nbeam_tr)  = beam_data%abeama(1:nbeam_tr)    ! A
  zbeams%nlco(1:nbeam_tr)    = beam_data%nlco(1:nbeam_tr)        ! co (with current) or ctr inj.
  zbeams%ntrace(1:nbeam_tr)  = beam_data%ntrace(1:nbeam_tr)    ! trace element option

  zbeams%rtcena(1:nbeam_tr)  = beam_data%rtcena(1:nbeam_tr)    ! geometry details...
  zbeams%xlbtna(1:nbeam_tr)  = beam_data%xlbtna(1:nbeam_tr)
  zbeams%xybsca(1:nbeam_tr)  = beam_data%xybsca(1:nbeam_tr)
  zbeams%nbshapa(1:nbeam_tr) = beam_data%nbshapa(1:nbeam_tr)  ! ion source grid...
  zbeams%bmwidra(1:nbeam_tr) = beam_data%bmwidra(1:nbeam_tr)
  zbeams%bmwidza(1:nbeam_tr) = beam_data%bmwidza(1:nbeam_tr)
  zbeams%divra(1:nbeam_tr)   = beam_data%divra(1:nbeam_tr)      ! divergence and focal length...
  zbeams%divza(1:nbeam_tr)   = beam_data%divza(1:nbeam_tr)
  zbeams%foclra(1:nbeam_tr)  = beam_data%foclra(1:nbeam_tr)
  zbeams%foclza(1:nbeam_tr)  = beam_data%foclza(1:nbeam_tr)
  zbeams%nbapsha(1:nbeam_tr) = beam_data%nbapsha(1:nbeam_tr)  ! aperture geometry & position
  zbeams%rapedga(1:nbeam_tr) = beam_data%rapedga(1:nbeam_tr)
  zbeams%xzpedga(1:nbeam_tr) = beam_data%xzpedga(1:nbeam_tr)
  zbeams%xlbapa(1:nbeam_tr)  = beam_data%xlbapa(1:nbeam_tr)
  zbeams%xybapa(1:nbeam_tr)  = beam_data%xybapa(1:nbeam_tr)
  zbeams%nbapsh2(1:nbeam_tr) = beam_data%nbapsh2(1:nbeam_tr)  ! optional 2nd aperture
  zbeams%rapedg2(1:nbeam_tr) = beam_data%rapedg2(1:nbeam_tr)
  zbeams%xzpedg2(1:nbeam_tr) = beam_data%xzpedg2(1:nbeam_tr)
  zbeams%xlbapa2(1:nbeam_tr) = beam_data%xlbapa2(1:nbeam_tr)
  PRINT *,'zbeams data dump:'
  PRINT *,'no. thermal  species, zbeams%ngmax =',zbeams%ngmax
  PRINT *,'Z of thermals  zbeams%backz ', zbeams%backz(1:zbeams%ngmax)
  PRINT *,'A of thermals zbeams%aplasm ',zbeams%aplasm(1:zbeams%ngmax)
  PRINT *,'no. of beams, zbeams%nbeam =  ',zbeams%nbeam
  PRINT *,'Z of beams  zbeams%xzbeama ', zbeams%xzbeama(1:zbeams%nbeam)
  PRINT *,'A of beams zbeams%abeama ',zbeams%abeama(1:zbeams%nbeam)
  PRINT *,'co or ctr injection, zbeams%nlco ',zbeams%nlco(1:zbeams%nbeam)
  PRINT *,'beam trace index zbeams%ntrace ',zbeams%ntrace(1:zbeams%nbeam)
  PRINT *,'tangency radius, zbeams%rtcena =',zbeams%rtcena(1:zbeams%nbeam)
  PRINT *,'dist source to tangncy pt zbeams%xlbtna =',zbeams%xlbtna(1:zbeams%nbeam)
  PRINT *,'aperture shape,zbeams%nbshapa =',zbeams%nbshapa(1:zbeams%nbeam)
  PRINT *,'source grid half width, zbeams%bmwidra ', zbeams%bmwidra(1:zbeams%nbeam)
  PRINT *,'source grid half height,zbeams%bmwidza ',zbeams%bmwidza(1:zbeams%nbeam)
  PRINT *,'horizontal divergence of beam (radians) zbeams%divra ',zbeams%divra(1:zbeams%nbeam)
  PRINT *,'vertical divergence of beam (radians),zbeams%divza ',zbeams%divza(1:zbeams%nbeam)
  PRINT *,'horizontal focal length of beam (cm),zbeams%foclra ',zbeams%foclra(1:zbeams%nbeam)
  PRINT *,'vertical focal length of beam (cm),zbeams%foclza ',zbeams%foclza(1:zbeams%nbeam)
  PRINT *,'shape of main aperture to vacuum vessel :zbeams%nbapsha ',zbeams%nbapsha(1:zbeams%nbeam)
  PRINT *,' aperture half-width,zbeams%rapedga ',zbeams%rapedga(1:zbeams%nbeam)
  PRINT *,'aperture half-height,zbeams%xzpedga ',zbeams%xzpedga(1:zbeams%nbeam)
  PRINT *,'distance, source grid to aperture, cm,zbeams%xlbapa ',zbeams%xlbapa(1:zbeams%nbeam)
  PRINT *,'shape of optional 2nd aperture,zbeams%nbapsh2 ',zbeams%nbapsh2(1:zbeams%nbeam)
  PRINT *,'aperture half-width ,zbeams%rapedg2 ',zbeams%rapedg2(1:zbeams%nbeam)
  PRINT *,'aperture half-height,zbeams%xzpedg2 ', zbeams%xzpedg2(1:zbeams%nbeam)
  PRINT *,'distance, source grid to aperture,zbeams%xlbapa2 ',zbeams%xlbapa2(1:zbeams%nbeam)

  CALL nbi_add_beams(zbeams,ierr)   ! apply...
  IF(ierr.NE.0) THEN
     CALL errmsg_exit(' ?set_nbinputs:  nubeam/nbi_add_beams error code set.')
  ENDIF
   PRINT *,'in sub set_nbinputs ,done  with nbi_add_beams'






  !--------------------------------------------------------------
  !  the impurity species...

  CALL nbi_init_impurity(zimpurity)  ! defaults
  nrhix = beam_data%nrhix
  zimpurity%nrhix = nrhix            ! number

  zimpurity%xzimpx(1:nrhix) = beam_data%xzimpx(1:nrhix)  ! Z (e.g. 6 for stripped Carbon)
  zimpurity%aimpx(1:nrhix) = beam_data%aimpx(1:nrhix)    ! A (AMU)

  CALL nbi_add_impurity(zimpurity,ierr)  ! apply...
  IF(ierr.NE.0) THEN
     CALL errmsg_exit(' ?set_nbinputs:  nubeam/nbi_add_impurity error code set.')
  ENDIF
   PRINT *,'in sub set_nbinputs ,done  with nbi_add_impurity'



  !--------------------------------------------------------------
  !  RF minority species... fast ion species computed by an RF code...
  !                         (or may be thermalized)

  CALL nbi_init_minority(zminority)

  zminority%nmini = nmini            ! number
  zminority%xzmini(1:nmini) = xzmini(1:nmini)  ! Z (e.g. 2 for He3)
  zminority%amini(1:nmini) = amini(1:nmini)    ! A (e.g. 3 for He3)

  CALL nbi_add_minority(zminority,ierr)  ! apply...
  IF(ierr.NE.0) THEN
     CALL errmsg_exit(' ?set_nbinputs:  nubeam/nbi_add_minority error code set.')
  ENDIF
   PRINT *,'in sub set_nbinputs ,done  with nbi_add_minority'




  !--------------------------------------------------------------
  !  Fusion product species... He3, He4(alphas), Tritons, Protons
  !    (note dmc 3/1/2002: proton source does not work; energy is 
  !     too high for collision operator currently installed; this
  !     could be fixed but it would be a project).

  CALL nbi_init_fusion(zfusion)

  zfusion%nlfhe3 = nlfhe3      ! switches to activate fusion product species
  zfusion%nlfhe4 = nlfhe4
  zfusion%nlfst = nlfst
  zfusion%nlfsp = nlfsp

  zfusion%nlusf3 = nlusf3      ! switches to supply birth rate data
  zfusion%nlusfa = nlusfa      ! otherwise use internal fusion rate calculation
  zfusion%nlusft = nlusft
  zfusion%nlusfp = nlusfp

  zfusion%plfhe3 = plfhe3      ! minimum power threshold for full MC statistics
  zfusion%plfhe4 = plfhe4
  zfusion%plfst = plfst
  zfusion%plfsp = plfsp

  zfusion%nlfatom = nlfatom    ! switch to compute atomic physics effects
                               ! on slowing down fusion product ions

  CALL nbi_set_fusion(zfusion,ierr)  ! apply...
  IF(ierr.NE.0) THEN
     CALL errmsg_exit(' ?set_nbinputs:  nubeam/nbi_set_fusion error code set.')
  ENDIF
   PRINT *,'in sub set_nbinputs ,done  with nbi_set_fusion'



  !--------------------------------------------------------------
  !  "fpp" block -- related to support of Fokker Planck (fpp) fast ion
  !                 module.  The Monte Carlo model can be operated
  !                 to calculate deposition, leaving slowing down to
  !                 a fokker planck model.
  IF(fpp) THEN
     CALL nbi_init_fpp(zfpp)

     zfpp%nlbfpp = nlbfpp           ! (set .FALSE. in sub init) if fpp 
     ! is to compute slowing down
     ! a deposition distribution function
     ! in the plasma frame will be computed.

     zfpp%nefpdep = mifpen          ! energy grid for deposition distribution 
     ! mifpen set in nbi_dimensions
     ALLOCATE (efpdep(1:mifpen),STAT = istat)
     IF(istat .NE. 0) &
          CALL allocate_error("efpdep,sub init",0,istat)
     efpdep(:) =0.0                 !0.0 for now, eventually control in inone
     zfpp%efpdep(1:mifpen) = efpdep(1:mifpen)


     zfpp%nxfpdep = mifpxi          ! pitch angle grid for deposition distribution
     ALLOCATE (xfpdep(1:mifpxi),STAT = istat)
     IF(istat .NE. 0) &
          CALL allocate_error("xfpdep,sub init",0,istat)
     xfpdep(:) =0.0               !0.0 for now, eventually control in inone
     zfpp%xfpdep(1:mifpxi) = xfpdep(1:mifpxi)



     zfpp%nxfpsh = nzone_fp         ! number of radial zones (norm.sqrt.tor.flux)


     CALL nbi_set_fpp(zfpp,ierr)  ! apply...
     IF(ierr.NE.0) THEN
        CALL errmsg_exit(' ?set_nbinputs:  nubeam/nbi_set_fpp error code set.')
     ENDIF
     PRINT *,'in sub set_nbinputs ,done  with nbi_set_fpp'

  ENDIF



  !--------------------------------------------------------------
  !  "num" block -- numerical controls and adjustments

  CALL nbi_init_num(znum)   ! set DEFAULTS (see nbspec.dat).

  znum%cxsplt = cxsplt      ! cx-splitting, default = 2.0 is generally OK
  znum%dxbsmoo = dxbsmoo    ! 1d output profile smoothing, default = 0.05 OK

  znum%nptcls = nptcls      ! *** number of MC model ions for BEAM SPECIES ***
  znum%nptclf = nptclf      ! *** number of MC ions for FUSION PRODUCTS ***
  znum%nchdvp = nchdvp      ! orbit integrator option  HSJ
  ! here TRANSP computes the max nptcls it will ever request...
  ! TRANSP "time zoom" feature can request a temporary increase...

  iptcls_max = MAX(nptcls,nptclf)
!  following noit implemented (no info on tinit,tzoom) HSJ
!  do i=1,8
!     if((tinit.le.tzoom(1,i)).and.(tzoom(1,i).le.ftime)) then
!        iptcls_max=max(iptcls_max,int(pzoom(8)+0.5))
!     endif
!  enddo

  znum%nptcl_max = iptcls_max

  znum%ndep0 = ndep0        ! minimum number of deposition tracks

  znum%nbbcal = beam_data%nbbcal ! number of evaluation points, fast ion - fast ion
                            ! fusion Monte Carlo integrals; default = 100000 OK

  znum%wghta = wghta        ! profile adjustment on MCs statistics
                            ! larger WGHTA -> better statistics at core
                            ! smaller WGHTA -> better statistics at edge

  znum%cxpcon = cxpcon      ! adjustment: frequency of atomic physics calls
  znum%fppcon = fppcon      ! adjustment: frequency of collision op. calls

  znum%dtn = dtn_orbit      ! *** ORBIT TIMESTEP PARAMETER *** CAUTION!!!
                            ! 0.0 default means automatic step adjustment

  CALL nbi_set_num(znum,ierr)  ! apply...
  IF(ierr.NE.0) THEN
     CALL errmsg_exit(' ?set_nbinputs:  nubeam/nbi_set_num error code set.')
  ENDIF
   PRINT *,'in sub set_nbinputs ,done  with nbi_set_num'




  !--------------------------------------------------------------
  !  Atomic physics options

  CALL nbi_init_atomic(zatomic)  ! defaults

  zatomic%xdepmod = xdepmod  ! anomolous deposition adjustment (default=1.0)
                             ! >1 to make plasma more opaque for deposition
                             ! <1 to make plasma more transparent

  zatomic%nsigexc = nsigexc  ! =1 for (simple) excited states correction
                             ! for neutral beam deposition calculation

  zatomic%nlminsv = nlminsv  ! .TRUE. -- loop over all impurities
                             ! .FALSE. -- use "average Z" approximation

  zatomic%nlbbcx = nlbbcx    ! .TRUE. to include beam-beam atomic physics

  zatomic%nlebei = nlebei    ! .TRUE. for beam energy correction to 
                             ! <sigma*v> electron impact ionization

  zatomic%dn0out = dn0out    ! neutral density (/cm3) assumed beyond plasma bdy

  zatomic%nmsigx = nmsigx    ! impurity stopping cross section option

  CALL nbi_set_atomic(zatomic,ierr)  ! apply...
  IF(ierr.NE.0) THEN
     CALL errmsg_exit(' ?set_nbinputs:  nubeam/nbi_set_atomic error code set.')
  ENDIF
   PRINT *,'in sub set_nbinputs ,done  with nbi_set_atomic'
  !--------------------------------------------------------------
  !  collision operator options...

  CALL nbi_init_collid(zcollid)  ! defaults

  !  the following anomolous adjustment factors all default to 1.0; they
  !  are rarely set to any other value...

  zcollid%xcfanbi = xcfanbi  ! electron drag & energy diffusion -- beam ions
  zcollid%xdfanbi = xdfanbi  ! ion drag & energy diffusion -- beam ions
  zcollid%xefanbi = xefanbi  ! pitch angle scattering -- beam ions

  zcollid%xcfafus = xcfafus  ! electron drag & energy diffusion -- fusion products
  zcollid%xdfafus = xdfafus  ! ion drag & energy diffusion -- fusion products
  zcollid%xefafus = xefafus  ! pitch angle scattering -- fusion products

  !  default .TRUE.:

  zcollid%nlbcde = nlbcde    ! .FALSE. to suppress energy diffusion
  zcollid%nlbcoh = nlbcoh    ! .FALSE. to suppress OH acceleration of fast ions
                             ! (due to toroidal electric field)
  zcollid%nlbcpa = nlbcpa    ! .FALSE. to suppress pitch angle scattering

  !  default .FALSE.:

  zcollid%nlorbo = nlorbo    ! .TRUE. to suppress all collisions(for debugging)

  CALL nbi_set_collid(zcollid,ierr)  ! apply...
  IF(ierr.NE.0) THEN
     CALL errmsg_exit(' ?set_nbinputs:  nubeam/nbi_set_collid error code set.')
  ENDIF
   PRINT *,'in sub set_nbinputs ,done  with nbi_set_collid'



  !--------------------------------------------------------------
  !  finite larmor radius corrections

  CALL nbi_init_flr(zflr)  ! defaults...

  zflr%nlbflr = nlbflr     ! .TRUE. for standard treatment (default)
  zflr%nlbgflr = nlbgflr   ! .TRUE. for enhanced treatment (expensive)


  zflr%gflr_op(15) = gflr_op(15)  ! (this is a TRANSP cobble)

  CALL nbi_set_flr(zflr,ierr)  ! apply...
  IF(ierr.NE.0) THEN
     CALL errmsg_exit(' ?set_nbinputs:  nubeam/nbi_set_flr error code set.')
  ENDIF
    PRINT *,'in sub set_nbinputs ,done  with nbi_set_flr'
  !--------------------------------------------------------------
  !  sawtooth options
  IF(sawtooth)THEN
     CALL nbi_init_saw(zsaw)  ! defaults...

     zsaw%nlsawb = .FALSE.        ! .TRUE. for sawtooth mixing of beam ions
     zsaw%nlsawf = .FALSE.        ! .TRUE. for sawtooth mixing of fusion products
     ! sawtooth controls may be set .TRUE. in set_nbupdate ... when sawteeth
     ! are expected.

     CALL nbi_set_saw(zsaw,ierr)  ! apply...
     IF(ierr.NE.0) THEN
        CALL errmsg_exit(' ?set_nbinputs:  nubeam/nbi_set_saw error code set.')
     ENDIF
     PRINT *,'in sub set_nbinputs ,done  with nbi_set_saw'
  ENDIF
 
  !--------------------------------------------------------------
  !  anomolous diffusion options (to activate the model, provide an
  !  anomolous diffusivity profile... see set_nbprofiles(...))
  IF(anom_diff)THEN
     CALL nbi_init_adif(zadif)  ! defaults...

     zadif%nkdifb = nkdifb      ! =3: apply to all Monte Carlo fast ions
     ! =1: beam ions only; =2: fusion products only.

     zadif%ndifbe = ndifbe      ! energy dependent adjustment (see nbspec.dat).
     zadif%fdifbe(1:ndifbe) = fdifbe(1:ndifbe)
     zadif%edifbe(1:ndifbe) = edifbe(1:ndifbe)

     CALL nbi_set_adif(zadif,ierr)  ! apply...
     IF(ierr.NE.0) THEN
        CALL errmsg_exit(' ?set_nbinputs:  nubeam/nbi_set_adif error code set.')
     ENDIF
     PRINT *,'in sub set_nbinputs ,done  with nbi_set_adif'

 ENDIF

  !--------------------------------------------------------------
  !  controls for TF ripple loss models (caution, these are crude...)
  !  see nbspec.dat for more details


  IF(ripple)THEN
     CALL nbi_init_ripple(zripple)  ! defaults...

     zripple%nrip = nrip            ! model selection switch (default=0, no loss).
     zripple%taurip = taurip
     CALL eq_gfnum('LOG_TF_RIPPLE',id_ripple)
     zripple%id_ripple = id_ripple
     zripple%ncoils = ncoils
     zripple%asrd = asrd
     zripple%bsrd = bsrd

     CALL nbi_set_ripple(zripple,ierr)  ! apply...
     IF(ierr.NE.0) THEN
        CALL errmsg_exit(' ?set_nbinputs:  nubeam/nbi_set_ripple error code set.')
     ENDIF
     PRINT *,'in sub set_nbinputs ,done  with nbi_set_ripple'

  ENDIF


  !--------------------------------------------------------------
  ! output options:  specify optional energy-range-restricted density and 
  ! trapping fraction profile output.
  if(outcon)then
     CALL nbi_init_outcon(zoutcon)  ! default:  nerngfi=0

     zoutcon%nerngfi = nerngfi
     zoutcon%erngfi(1:nerngfi) = erngfi(1:nerngfi)
     zoutcon%nsdbgb = nsdbgb

     CALL nbi_set_outcon(zoutcon,ierr)  ! apply...
     IF(ierr.NE.0) THEN
        CALL errmsg_exit(' ?set_nbinputs:  nubeam/nbi_set_outcon error code set.')
     ENDIF
     PRINT *,'in sub set_nbinputs ,done  with nbi_set_outcon'

  endif
  !--------------------------------------------------------------
  !  an ancient fishbone loss model... details: see nbspec.dat
  if(fishbone)then
     CALL nbi_init_fishbone(zfbone)  ! defaults...

     zfbone%nlfbon = nlfbon
     zfbone%fbemin = fbemin
     zfbone%fbemax = fbemax
     zfbone%fvpvmn = fvpvmn
     zfbone%fvpvmx = fvpvmx
     zfbone%fbltim = fbltim
     zfbone%fshper = fshper
     zfbone%fshwid = fshwid
     zfbone%tfshon = tfshon
     zfbone%tfshof = tfshof

     CALL nbi_set_fishbone(zfbone,ierr)  ! apply...
     IF(ierr.NE.0) THEN
        CALL errmsg_exit(' ?set_nbinputs:  nubeam/nbi_set_fishbone error code set.')
     ENDIF
     PRINT *,'in sub set_nbinputs ,done  with nbi_set_fishbones'


  endif

  !--------------------------------------------------------------
  !  optional "beam in box" 3d neutral density calculations.
  !   ...these cause a seldom used special output to be created.
  !   ...for details on the individual elements see nbspec.dat.
  if(beam_box)then
     CALL nbi_init_box(zbox)  ! defaults...

     zbox%nbbox = nbbox
     zbox%nbsbox(1:nbbox) = nbsbox(1:nbbox)
     zbox%nbebox(1:nbbox) = nbebox(1:nbbox)
     zbox%nxbox = nxbox
     zbox%nybox = nybox
     zbox%nlbox = nlbox
     zbox%xboxhw = xboxhw
     zbox%yboxhw = yboxhw
     zbox%xlbox1 = xlbox1
     zbox%xlbox2 = xlbox2
     zbox%ndepbox = ndepbox

     CALL nbi_set_box(zbox,ierr)  ! apply...

     IF(ierr.NE.0) THEN
        CALL errmsg_exit(' ?set_nbinputs:  nubeam/nbi_set_box error code set.')
     ENDIF
     PRINT *,'in sub set_nbinputs ,done  with nbi_set_box'

  endif


  !--------------------------------------------------------------
  ! code snippet:  block misc

  CALL nbi_init_misc(zmisc)  ! defaults...

  zmisc%plhgt = plhgt    ! updown symmetric runs: displacement of 
                         ! plasma midplane relative to vacuum vessel
                         ! Z=0 midplane.

  zmisc%nlcprb = nlcprb  ! set .FALSE. to suppress compression operator and
                         ! other effects of time varying MHD equilibrium.

  zmisc%edbfac = edbfac  ! debug tripwire: halt simulation if an ion with 
                         ! edbfac*(its birth energy) occurs.

  zmisc%nmcurb = nmcurb  ! =1 (default) for beam driven current with
                         ! aspect ratio approximate neoclassical shielding
                         ! factor correction.

  zmisc%nlfdep = nlfdep  ! .TRUE. to compute deposition distribution function
                                 ! in the lab frame on the "fbm" grid.

  CALL nbi_set_misc(zmisc,ierr)  ! apply...
  IF(ierr.NE.0) THEN
     CALL errmsg_exit(' ?set_nbinputs:  nubeam/nbi_set_misc error code set.')
  ENDIF
    PRINT *,'in sub set_nbinputs ,done  with nbi_set_misc'
  !---------------------
  !---------------------
END SUBROUTINE set_nbinputs
