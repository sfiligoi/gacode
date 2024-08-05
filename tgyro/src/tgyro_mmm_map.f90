!\***************************************************************
! pt_mmm_mod.f90: pt_solver implementation of mmm
! series turbulent model.
!/***************************************************************
MODULE tgyro_mmm_mod
  use modmmm, only : MMM_NCH
  implicit none

  integer, parameter, private :: r8 = selected_real_kind(8) 
  ! public available varaibles
  PUBLIC
  
  real(r8)   :: mmm_chi_i       !Ion thermal diffusivity
  real(r8)   :: mmm_chi_e       !Electron thermal diffusivity
  real(r8)   :: mmm_chi_ne      !Electron particle diffusivity
  real(r8)   :: mmm_chi_nx      !Impurity ion diffusivity
  real(r8)   :: mmm_chi_phi     !Toroidal momentum diffusivity
  real(r8)   :: mmm_chi_theta   !Poloidal momentum diffusivity

  real(r8)   :: mmm_vheat_i     !Electron thermal
  real(r8)   :: mmm_vheat_e     !Ion thermal
  real(r8)   :: mmm_vgx_ne      !Particle
  real(r8)   :: mmm_vgx_nx      !Impurity particle
  real(r8)   :: mmm_vgx_phi     !Toroidal momentum
  ! 
  real(r8), dimension(4,1) :: mmm_gamma_w20
  real(r8), dimension(4,1) :: mmm_omega_w20
  real(r8), dimension(1)   :: mmm_gamma_drbm
  real(r8), dimension(1)   :: mmm_omega_drbm
  real(r8), dimension(1)   :: mmm_gamma_mtm
  real(r8), dimension(1)   :: mmm_omega_mtm
  real(r8), dimension(1)   :: mmm_kyrhosmtm

  real(r8), dimension(1)   :: mmm_chii_w20
  real(r8), dimension(1)   :: mmm_chie_w20
  real(r8), dimension(1)   :: mmm_chini_w20

  real(r8), dimension(1)   :: mmm_chii_drbm
  real(r8), dimension(1)   :: mmm_chie_drbm
  real(r8), dimension(1)   :: mmm_chini_drbm

  real(r8), dimension(1)   :: mmm_chie_mtm
  real(r8), dimension(1)   :: mmm_chie_etg
  
  real(r8), dimension(1)   :: mmm_dbsqprf

  real(r8), dimension(MMM_NCH,1) :: mmm_velthi
  real(r8), dimension(MMM_NCH,1) :: mmm_flux

  !\
  ! private variables
  !/
  real(r8), private :: rmaj_exp
  real(r8), dimension(1), private :: rmajor_exp
  real(r8), dimension(1), private :: rminor_exp
  real(r8), dimension(1), private :: elong_exp
  real(r8), dimension(1), private :: te_exp
  real(r8), dimension(1), private :: ti_exp
  real(r8), dimension(1), private :: ne_exp
  real(r8), dimension(1), private :: ni_exp
  real(r8), dimension(1), private :: nt_exp
  real(r8), dimension(1), private :: nh_exp
  real(r8), dimension(1), private :: nz_exp
  real(r8), dimension(1), private :: nf_exp
  real(r8), dimension(1), private :: zi_exp
  real(r8), dimension(1), private :: mi_exp
  real(r8), dimension(1), private :: mh_exp
  real(r8), dimension(1), private :: zz_exp
  real(r8), dimension(1), private :: mz_exp
  real(r8), dimension(1), private :: zeff_exp
  real(r8), dimension(1), private :: q_exp
  real(r8), dimension(1), private :: btor_exp
  real(r8), dimension(1), private :: gradne_exp
  real(r8), dimension(1), private :: gradnh_exp
  real(r8), dimension(1), private :: gradnt_exp
  real(r8), dimension(1), private :: gradnz_exp
  real(r8), dimension(1), private :: gradte_exp
  real(r8), dimension(1), private :: gradti_exp
  real(r8), dimension(1), private :: gradq_exp
  real(r8), dimension(1), private :: gradelg_exp
  real(r8), dimension(1), private :: rlne_exp
  real(r8), dimension(1), private :: rlnh_exp
  real(r8), dimension(1), private :: rlnt_exp
  real(r8), dimension(1), private :: rlnz_exp
  real(r8), dimension(1), private :: rlte_exp
  real(r8), dimension(1), private :: rlti_exp
  real(r8), dimension(1), private :: shear_exp
  real(r8), dimension(1), private :: gamma_e_exp

  ! velocity for momentum prediction
  real(r8), dimension(1), private :: vtor_exp
  real(r8), dimension(1), private :: vpol_exp
  real(r8), dimension(1), private :: vpar_exp
  real(r8), dimension(1), private :: gradvpar_exp
  real(r8), dimension(1), private :: gradvtor_exp
  real(r8), dimension(1), private :: gradvpol_exp
  real(r8), dimension(1), private :: etanc_exp

  ! options & switches
  !
  !------- Input Variables -------------------------------------------
  ! The content of the input variables follows their corresponding arguments
  ! of the mmm subroutine. See modmmm.f90 for more details.
  integer, private, save :: lprint
  real(r8), private :: epslon = 1.0d-12

  real(r8), allocatable, dimension(:), save, private   :: cmodel  ! Internal model weights
  real(r8), allocatable, dimension(:,:), save, private :: cmmm07
  integer, allocatable, dimension(:,:), save, private  :: lmmm07

  ! function accessibility
  PUBLIC :: tgyro_mmm_map
  PUBLIC :: tgyro_mmm_run

CONTAINS
  !\------------------------------------------
  ! load data for mmm mode
  !/
  SUBROUTINE tgyro_mmm_map(ier)
    use modmmm, only : KW20, KDBM, KMTM, KETG, MAXNOPT, &
        MMM_NMODE, set_mmm_switches, KETGM
    use tgyro_globals
    implicit none

    ! inputs
    integer, intent(out)  :: ier

    ! local variables
    integer  :: i
    !
    integer, parameter :: hfIn = 34
    !
    real(r8) :: nimp_tot,nmain_tot
    logical :: ismoo
    logical, save :: initialized = .false.
    real(r8), parameter :: BADREAL = -1E10_r8
    integer, parameter  :: BADINT = -1000000
    real(r8), dimension(MAXNOPT) :: &
       cW20, cDBM, cETG, cMTM ! To be passed to cmmm
    character(len=9) :: mmm_namelist = 'input.mmm'
    integer, dimension(MAXNOPT) :: &
       lW20, lDBM, lETG, lMTM ! To be passed to lmmm

    NAMELIST /mmm_input/                         &
       cmodel, lprint, cW20, cDBM, cETG,         &
        cMTM, lW20, lDBM, lETG, lMTM
    
    ! allocate memory

    rmajor_exp(1) = 1d-2 * r_maj(i_r)       ! [m]
    if (r(i_r).eq.0d0) then
      rminor_exp(1) = epslon
    else  
      rminor_exp(1) = 1d-2 * r(i_r)         ! [m]
    endif  
    rmaj_exp = 1d-2 * r_maj(n_r)            ! [m]
    ! Btor_exp
    btor_exp(1) = 1d-4*b_ref(i_r)           ! [T] 
    ! te/ti, ne/ni/nh/nt/nf/nz
    te_exp(1) = 1d-3*te(i_r)                ! [keV]
    ne_exp(1) = 1d6*ne(i_r)                 ! [m^-3]
    ti_exp(1) = 0_r8
    zi_exp(1) = 0_r8
    mi_exp(1) = 0_r8
    mh_exp(1) = 0_r8
    zz_exp(1) = 0_r8
    mz_exp(1) = 0_r8
    ni_exp(1) = 0_r8
    nh_exp(1) = 0_r8
    nt_exp(1) = 0_r8
    nf_exp(1) = 0_r8
    nz_exp(1) = 0_r8
    do i = 1, loc_n_ion
      if (calc_flag(i) == 0) cycle
      ti_exp(1) = ti_exp(1) + ti(i,i_r)*ni(i,i_r)
      zi_exp(1) = zi_exp(1) + zi_vec(i)*ni(i,i_r)
      mi_exp(1) = mi_exp(1) + mi(i)*ni(i,i_r)
      ni_exp(1) = ni_exp(1) + ni(i,i_r)
      if (zi_vec(i).eq.1) then
        nh_exp(1) = nh_exp(1) + ni(i,i_r)
        mh_exp(1) = mh_exp(1) + mi(i)*ni(i,i_r)
      endif  
      if (therm_flag(i).eq.1) then
        nt_exp(1) = nt_exp(1) + ni(i,i_r)
      else  
        nf_exp(1) = nf_exp(1) + ni(i,i_r)
      endif
      if (zi_vec(i)>2.5) then
        nz_exp(1) = nz_exp(1) + ni(i,i_r)
        zz_exp(1) = zz_exp(1) + zi_vec(i)*ni(i,i_r)
        mz_exp(1) = mz_exp(1) + mi(i)*ni(i,i_r)
      endif  
    enddo
    ti_exp(1) = 1d-3*ti_exp(1)/ni_exp(1)       ! [keV]
    zi_exp(1) = zi_exp(1)/ni_exp(1)
    mi_exp(1) = mi_exp(1)/ni_exp(1)/md*2
    mh_exp(1) = mh_exp(1)/nh_exp(1)/md*2
    zz_exp(1) = zz_exp(1)/nz_exp(1)
    mz_exp(1) = mz_exp(1)/nz_exp(1)/md*2 
    ni_exp(1) = 1d6*ni_exp(1)                  ! [m^-3]
    nh_exp(1) = 1d6*nh_exp(1)                  ! [m^-3]
    nt_exp(1) = 1d6*nt_exp(1)                  ! [m^-3]
    nf_exp(1) = 1d6*nf_exp(1)                  ! [m^-3]
    nz_exp(1) = 1d6*nz_exp(1)                  ! [m^-3]
    ! zeff
    zeff_exp(1) = z_eff(i_r)
    ! q profile
    q_exp(1) = abs(q(i_r))
    ! elongation
    elong_exp(1) = kappa(i_r)
    ! resistivity
    !# ...
    !\-------------------------------
    ! calculate the gradients

    ! density gradients
    rlne_exp(1) = 1d2*rmajor_exp(1)*dlnnedr(i_r)
    rlte_exp(1) = 1d2*rmajor_exp(1)*dlntedr(i_r)
    rlti_exp(1) = 0_r8
    rlnh_exp(1) = 0_r8
    rlnt_exp(1) = 0_r8
    rlnz_exp(1) = 0_r8
    do i = 1, loc_n_ion
      if (calc_flag(i) == 0) cycle
      if (zi_vec(i).eq.1) then
        rlnh_exp(1) = rlnh_exp(1) + dlnnidr(i,i_r)*ni(i,i_r)*1d6
      endif  
      if (therm_flag(i).eq.1) then
        rlti_exp(1) = rlti_exp(1) + dlntidr(i,i_r)*ni(i,i_r)*ti(i,i_r)*1d3
        rlnt_exp(1) = rlnt_exp(1) + dlnnidr(i,i_r)*ni(i,i_r)*1d6
      endif
      if (zi_vec(i)>2.5) then
        rlnz_exp(1) = rlnz_exp(1) + dlnnidr(i,i_r)*ni(i,i_r)*1d6
      endif  
    enddo
    rlti_exp(1) = 1d2*rmajor_exp(1)*rlti_exp(1)/(ni_exp(1)+epslon)/(ti_exp(i_r)+epslon)
    rlnh_exp(1) = 1d2*rmajor_exp(1)*rlnh_exp(1)/(nh_exp(1)+epslon)
    rlnt_exp(1) = 1d2*rmajor_exp(1)*rlnt_exp(1)/(nt_exp(1)+epslon)
    rlnz_exp(1) = 1d2*rmajor_exp(1)*rlnz_exp(1)/(nz_exp(1)+epslon)
    shear_exp(1) = rmajor_exp(1)*s(i_r)/rminor_exp(1)
    gradelg_exp(1) = rmajor_exp(1)*s_kappa(i_r)/rminor_exp(1)

    !\
    ! velocities and their derivatives
    !/
    vtor_exp(1) = w0(i_r)*rmajor_exp(1)
    vpar_exp(1) = vtor_exp(i_r)
    vpol_exp(1) = 0.
    do i = 1, loc_n_ion
      if (calc_flag(i) == 0) cycle
      vpol_exp(1) = vpol_exp(1) + v_pol(i,i_r)*ni(i,i_r)
    enddo  
    vpol_exp(1) = vpol_exp(1)*1d-2/ni_exp(1)     ! m/s
    gradvtor_exp(1) = -1.0d2*vtor_p(i_r)*rmajor_exp(1)/(vtor_exp(1)+epslon)
    gradvpol_exp(1) = -1.0d2*vpol_p(i_r)*rmajor_exp(1)/(vtor_exp(1)+epslon)
    gradvpar_exp(1) = -1.0d2*vtor_p(i_r)*rmajor_exp(1)/(vpar_exp(1)+epslon)
    !\ 
    ! ExB flow shear
    !/ 
    gamma_e_exp(1)  = -rminor_exp(1)*f_rot(i_r)*w0_norm/(q_exp(1))
    
    !\
    ! Load MMM controls
    !/
    if (.not.(initialized)) then
      open( hfIn, file=mmm_namelist, &
         form='formatted', status='old', iostat=ier)
      if ( ier /= 0 ) then
        print *, 'Cannot open MMM input!'
        return
      endif
      
      if (.not.(allocated(cmodel))) allocate(cmodel(MMM_NMODE))
      if (.not.(allocated(cmmm07))) allocate(cmmm07(MAXNOPT, MMM_NMODE))
      if (.not.(allocated(lmmm07))) allocate(lmmm07(MAXNOPT, MMM_NMODE))

      lprint = 0
      cmodel = BADREAL
      cW20 = BADREAL; cDBM = BADREAL; cETG = BADREAL; cMTM = BADREAL
      lW20 = BADINT;  lDBM = BADINT;  lETG = BADINT;  lMTM = BADINT

      Read( hfIn, NML=mmm_input )
      Close( hfIn )
      !.. Fill parameter arrays with default values
      Call set_mmm_switches( cmmm = cmmm07, lmmm = lmmm07 )

      !.. Assign user specified parameters
      do i=1, MAXNOPT
        if ( abs( cW20(i) - BADREAL ) > epslon ) cmmm07(i,KW20)=cW20(i)
        if ( abs( cDBM(i) - BADREAL ) > epslon ) cmmm07(i,KDBM)=cDBM(i)
        if ( abs( cETG(i) - BADREAL ) > epslon ) cmmm07(i,KETG)=cETG(i)
        if ( abs( cMTM(i) - BADREAL ) > epslon ) cmmm07(i,KMTM)=cMTM(i)
        if ( lW20(i) /= BADINT ) lmmm07(i,KW20)=lW20(i)
        if ( lDBM(i) /= BADINT ) lmmm07(i,KDBM)=lDBM(i)
        if ( lETG(i) /= BADINT ) lmmm07(i,KETG)=lETG(i)
        if ( lMTM(i) /= BADINT ) lmmm07(i,KMTM)=lMTM(i)
      enddo

      do i=1, MMM_NMODE 
        if (cmodel(i) == BADINT ) cmodel(i) = 1
      enddo
      ! Disable the new ETG model for now
      cmodel(KETGM) = 0d0

      initialized = .true.
    endif

    return
  END SUBROUTINE tgyro_mmm_map

  !\
  ! run mmm model
  !/
  SUBROUTINE tgyro_mmm_run(ier)
    use modmmm
    use tgyro_globals
    implicit none

    ! output
    integer, intent(out) :: ier

    ! local variables
    integer :: j

    integer :: iprint = 0,  &! iprint > 0 for diagnostic printout
         iskip = 100, &! frequency of diagnostic printout
         initial = 0, &! = 0 the first time this routine is called
         inprout

    real(r8), dimension(1)   :: zthiig,zthdig,ztheig
    real(r8), dimension(1)   :: zthzig,zthtig,zthttig
    real(r8), dimension(6,1) :: zvelthi

    ier = 0

    inprout = 0

    ! calculate the characteristic length

    if (lmmm07(1,KDBM) /= 0 ) cmmm07(1,KDBM) = 0.0_r8 ! Disable ExB for DRIBM
    call mmm(                             &
         rmin     = rminor_exp,           &
         rmaj     = rmajor_exp,           &
         rmaj0    = rmaj_exp,             &
         elong    = elong_exp,            &
         ne       = ne_exp,               &
         nh       = nh_exp,               &
         nz       = nz_exp,               &
         nf       = nf_exp,               &
         zeff     = zeff_exp,             &
         te       = te_exp,               &
         ti       = ti_exp,               &
         q        = q_exp,                &
         btor     = btor_exp,             &
         zimp     = zz_exp,               &
         aimp     = mz_exp,               &
         ahyd     = mh_exp,               &
         aimass   = mi_exp,               &
         wexbs    = gamma_e_exp,          &
         gne      = rlne_exp,             &
         gni      = rlnt_exp,             &
         gnh      = rlnh_exp,             &
         gnz      = rlnz_exp,             &
         gte      = rlte_exp,             &
         gti      = rlti_exp,             &
         gq       = shear_exp,            &
         vtor     = vtor_exp,             &
         vpol     = vpol_exp,             &
         vpar     = vpar_exp,             &
         gvtor    = gradvtor_exp,         &
         gvpol    = gradvpol_exp,         &
         gvpar    = gradvpar_exp,         &
         gelong   = gradelg_exp,          &
         ! etanc    = etanc_exp,            &
         xti      = zthiig,               &
         xdi      = zthdig,               &
         xte      = ztheig,               &
         xdz      = zthzig,               &
         xvt      = zthtig,               &
         xvp      = zthttig,              &
         xtiW20   = mmm_chii_w20,         &
         xdiW20   = mmm_chini_w20,        &
         xteW20   = mmm_chie_w20,         &
         xtiDBM   = mmm_chii_drbm,        &
         xdiDBM   = mmm_chini_drbm,       &
         xteDBM   = mmm_chie_drbm,        &
         xteMTM   = mmm_chie_mtm,         &
         xteETG   = mmm_chie_etg,         &
         gammaW20 = mmm_gamma_w20(1:4,:), &
         omegaW20 = mmm_omega_w20(1:4,:), &
         gammaDBM = mmm_gamma_drbm,       &
         omegaDBM = mmm_omega_drbm,       &
         gammaMTM = mmm_gamma_mtm,        &
         omegaMTM = mmm_omega_mtm,        &
         kyrhosMTM = mmm_kyrhosmtm,       &
         dbsqprf  = mmm_dbsqprf,          &
         vconv    = mmm_velthi(1:6,:),    &
         vflux    = mmm_flux(1:6,:),      &
         npoints  = 1,                    &
         lprint   = iprint,               &
         nprout   = inprout,              &
         nerr     = ier,                  &
         cmodel   = cmodel,               &
         cswitch  = cmmm07,               &
         lswitch  = lmmm07)

    if ( ier /= 0 ) then 
      print *, 'Error in mmm subroutine!'
      return
    endif
    mmm_chi_i     = zthiig(1)*1d4      ! Ion thermal diffusivity
    mmm_chi_e     = ztheig(1)*1d4      ! Electron thermal diffusivity
    mmm_chi_ne    = zthdig(1)*1d4      ! Electron particle diffusivity
    mmm_chi_nx    = zthzig(1)*1d4      ! Impurity ion diffusivity
    mmm_chi_phi   = zthtig(1)*1d4      ! Toroidal momentum diffusivity
    mmm_chi_theta = zthttig(1)*1d4     ! Poloidal momentum diffusivity

    mmm_vheat_i = zvelthi(3,1)*1d2     ! Electron thermal
    mmm_vheat_e = zvelthi(1,1)*1d2     ! Ion thermal
    mmm_vgx_ne  = zvelthi(2,1)*1d2     ! Particle
    mmm_vgx_nx  = zvelthi(4,1)*1d2     ! Impurity particle
    mmm_vgx_phi = zvelthi(5,1)*1d2     ! Toroidal momentum
!     print *, 'mmm', i_r, mmm_chi_e*1d-4, mmm_chi_i*1d-4
!     print *, 'rmin     =',  rminor_exp
!     print *, 'rmaj     =',  rmajor_exp
!     print *, 'rmaj0    =',  rmaj_exp
!     print *, 'elong    =',  elong_exp
!     print *, 'ne       =',  ne_exp
!     print *, 'nh       =',  nh_exp
!     print *, 'nz       =',  nz_exp
!     print *, 'nf       =',  nf_exp
!     print *, 'zeff     =',  zeff_exp
!     print *, 'te       =',  te_exp
!     print *, 'ti       =',  ti_exp
!     print *, 'q        =',  q_exp
!     print *, 'btor     =',  btor_exp
!     print *, 'zimp     =',  zz_exp
!     print *, 'aimp     =',  mz_exp
!     print *, 'ahyd     =',  mh_exp
!     print *, 'aimass   =',  mi_exp
!     print *, 'wexbs    =',  gamma_e_exp
!     print *, 'gne      =',  rlne_exp, r(i_r), -r_maj(i_r)*(ne(i_r)-ne(i_r-1))/ne(i_r)/(r(i_r)-r(i_r-1)) 
!     print *, 'gni      =',  rlnt_exp
!     print *, 'gnh      =',  rlnh_exp
!     print *, 'gnz      =',  rlnz_exp, r(i_r), -r_maj(i_r)*(ni(2,i_r)-ni(2,i_r-1))/ni(2,i_r)/(r(i_r)-r(i_r-1))
!     print *, 'gte      =',  rlte_exp, r(i_r), -r_maj(i_r)*(te(i_r)-te(i_r-1))/te(i_r)/(r(i_r)-r(i_r-1))
!     print *, 'gti      =',  rlti_exp
!     print *, 'gq       =',  shear_exp
!     print *, 'vtor     =',  vtor_exp
!     print *, 'vpol     =',  vpol_exp
!     print *, 'vpar     =',  vpar_exp
!     print *, 'gvtor    =',  gradvtor_exp
!     print *, 'gvpol    =',  gradvpol_exp
!     print *, 'gvpar    =',  gradvpar_exp
!     print *, 'gelong   =',  gradelg_exp
    
    return
  END SUBROUTINE tgyro_mmm_run

END MODULE tgyro_mmm_mod
