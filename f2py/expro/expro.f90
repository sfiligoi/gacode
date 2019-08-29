module expro

  character(len=2) :: ident='# '
  double precision :: expro_mass_deuterium=3.34358e-24  ! md (g)
  

  !-----------------------------------------------------------
  ! DATA:
  !
  != Header entries

  integer :: expro_n_exp
  integer :: expro_n_ion
  character(len=10), dimension(20) :: expro_name
  character(len=10), dimension(20) :: expro_type

  double precision :: expro_masse=5.44887413e-4   ! me/m_H
  double precision, dimension(:), allocatable :: expro_mass
  double precision :: expro_ze=-1.0
  double precision, dimension(:), allocatable :: expro_z
  double precision :: expro_torfluxa=0.0
  double precision :: expro_rvbv=0.0
  double precision :: expro_ipa=0.0

  != 1D and 2D profile arrays contained in input.gacode

  double precision, dimension(:), allocatable :: expro_rho
  double precision, dimension(:), allocatable :: expro_rmin
  double precision, dimension(:), allocatable :: expro_polflux
  double precision, dimension(:), allocatable :: expro_q
  double precision, dimension(:), allocatable :: expro_w0
  double precision, dimension(:), allocatable :: expro_rmaj
  double precision, dimension(:), allocatable :: expro_zmag
  double precision, dimension(:), allocatable :: expro_kappa
  double precision, dimension(:), allocatable :: expro_delta
  double precision, dimension(:), allocatable :: expro_zeta
  double precision, dimension(:), allocatable :: expro_ne
  double precision, dimension(:,:), allocatable :: expro_ni
  double precision, dimension(:), allocatable :: expro_te
  double precision, dimension(:,:), allocatable :: expro_ti
  double precision, dimension(:), allocatable :: expro_ptot
  double precision, dimension(:), allocatable :: expro_johm
  double precision, dimension(:), allocatable :: expro_jbs
  double precision, dimension(:), allocatable :: expro_jrf
  double precision, dimension(:), allocatable :: expro_jnb
  double precision, dimension(:), allocatable :: expro_jbstor
  double precision, dimension(:), allocatable :: expro_sigmapar
  double precision, dimension(:), allocatable :: expro_z_eff
  double precision, dimension(:,:), allocatable :: expro_vpol
  double precision, dimension(:,:), allocatable :: expro_vtor
  double precision, dimension(:), allocatable :: expro_qohme
  double precision, dimension(:), allocatable :: expro_qbeame
  double precision, dimension(:), allocatable :: expro_qbeami
  double precision, dimension(:), allocatable :: expro_qrfe
  double precision, dimension(:), allocatable :: expro_qrfi
  double precision, dimension(:), allocatable :: expro_qfuse
  double precision, dimension(:), allocatable :: expro_qfusi
  double precision, dimension(:), allocatable :: expro_qbrem
  double precision, dimension(:), allocatable :: expro_qsync
  double precision, dimension(:), allocatable :: expro_qline
  double precision, dimension(:), allocatable :: expro_qei
  double precision, dimension(:), allocatable :: expro_qione
  double precision, dimension(:), allocatable :: expro_qioni
  double precision, dimension(:), allocatable :: expro_qcxi
  double precision, dimension(:), allocatable :: expro_qpar
  double precision, dimension(:), allocatable :: expro_qmom

  != 1D Derived quantities
  double precision, dimension(:), allocatable :: expro_bunit
  double precision, dimension(:), allocatable :: expro_gamma_e
  double precision, dimension(:), allocatable :: expro_gamma_p
  double precision, dimension(:), allocatable :: &
       expro_s,&
       expro_drmaj,&
       expro_dzmag,&
       expro_sdelta,&
       expro_skappa,&
       expro_szeta,&
       expro_dlnnedr,&
       expro_dlntedr,&
       expro_sdlnnedr,&
       expro_sdlntedr,&
       expro_dlnptotdr,&
       expro_w0p,&
       expro_vol,&
       expro_volp,&
       expro_cs,&
       expro_rhos,&
       expro_nuee,&
       expro_ni_new,&
       expro_dlnnidr_new,&
       expro_sdlnnidr_new,&
       expro_grad_r0,&
       expro_ave_grad_r,&
       expro_bp0,&
       expro_bt0,&
       expro_ip,&
       expro_mach,&
       expro_thetascale,&
       expro_flow_beam,&
       expro_flow_wall,&
       expro_flow_mom,&
       expro_pow_e,&
       expro_pow_i,&
       expro_pow_ei,&
       expro_pow_e_aux,&
       expro_pow_i_aux,&
       expro_pow_e_fus,&
       expro_pow_i_fus,&
       expro_pow_e_sync,&
       expro_pow_e_brem,&
       expro_pow_e_line

  != 2D Derived quantities
  double precision, dimension(:,:), allocatable :: &
       expro_dlnnidr,&
       expro_dlntidr,&
       expro_sdlnnidr,&
       expro_sdlntidr
  !-----------------------------------------------------------

  ! input.gacode.geo dimension and arrays

  integer :: expro_nfourier
  double precision, dimension(:,:,:), allocatable :: expro_geo
  double precision, dimension(:,:,:), allocatable :: expro_dgeo

  ! Field orientation parameters

  integer :: expro_signb
  integer :: expro_signq

  ! Control parameters

  integer :: expro_ctrl_n_ion = -1
  integer :: expro_ctrl_quasineutral_flag 
  integer :: expro_ctrl_numeq_flag
  integer :: expro_error=0

  ! Header information
  character(len=70) :: expro_head_original =  '#  *original : null'
  character(len=70) :: expro_head_statefile = '# *statefile : null'
  character(len=70) :: expro_head_gfile =     '#     *gfile : null'
  character(len=70) :: expro_head_cerfile =   '#   *cerfile : null'
  character(len=70) :: expro_head_vgen =      '#      *vgen : null'
  character(len=70) :: expro_head_tgyro =     '#     *tgyro : null'

contains

  subroutine expro_init(flag)

    implicit none

    integer :: nexp,nion
    integer, intent(in) :: flag

    if (flag == 1) then

       nexp = expro_n_exp
       nion = expro_n_ion

       allocate(expro_mass(nion))    ; expro_mass = 1.0
       allocate(expro_z(nion))       ; expro_z = 1.0

       allocate(expro_rho(nexp))     ; expro_rho = 0.0
       allocate(expro_rmin(nexp))    ; expro_rmin = 0.0
       allocate(expro_q(nexp))       ; expro_q = 0.0
       allocate(expro_polflux(nexp)) ; expro_polflux = 0.0
       allocate(expro_w0(nexp))      ; expro_w0 = 0.0
       allocate(expro_rmaj(nexp))    ; expro_rmaj = 0.0

       allocate(expro_zmag(nexp))    ; expro_zmag = 0.0
       allocate(expro_kappa(nexp))   ; expro_kappa = 0.0
       allocate(expro_delta(nexp))   ; expro_delta = 0.0
       allocate(expro_zeta(nexp))    ; expro_zeta = 0.0
       allocate(expro_ne(nexp))      ; expro_ne = 0.0
       allocate(expro_te(nexp))      ; expro_te = 0.0
       allocate(expro_ptot(nexp))    ; expro_ptot = 0.0
       allocate(expro_johm(nexp))    ; expro_johm = 0.0
       allocate(expro_jbs(nexp))     ; expro_jbs = 0.0
       allocate(expro_jrf(nexp))     ; expro_jrf = 0.0
       allocate(expro_jnb(nexp))     ; expro_jnb = 0.0
       allocate(expro_jbstor(nexp))  ; expro_jbstor = 0.0
       allocate(expro_sigmapar(nexp)); expro_sigmapar = 0.0
       allocate(expro_z_eff(nexp))   ; expro_z_eff = 0.0

       allocate(expro_qohme(nexp))  ; expro_qohme = 0.0
       allocate(expro_qbeame(nexp)) ; expro_qbeame = 0.0
       allocate(expro_qbeami(nexp)) ; expro_qbeami = 0.0
       allocate(expro_qrfe(nexp))   ; expro_qrfe = 0.0
       allocate(expro_qrfi(nexp))   ; expro_qrfi = 0.0
       allocate(expro_qfuse(nexp))  ; expro_qfuse = 0.0
       allocate(expro_qfusi(nexp))  ; expro_qfusi = 0.0
       allocate(expro_qbrem(nexp))  ; expro_qbrem = 0.0
       allocate(expro_qsync(nexp))  ; expro_qsync = 0.0
       allocate(expro_qline(nexp))  ; expro_qline = 0.0
       allocate(expro_qei(nexp))    ; expro_qei = 0.0
       allocate(expro_qione(nexp))  ; expro_qione = 0.0
       allocate(expro_qioni(nexp))  ; expro_qioni = 0.0
       allocate(expro_qcxi(nexp))   ; expro_qcxi = 0.0
       allocate(expro_qpar(nexp))   ; expro_qpar = 0.0
       allocate(expro_qmom(nexp))   ; expro_qmom = 0.0

       allocate(expro_ni(nion,nexp))   ; expro_ni = 0.0
       allocate(expro_ti(nion,nexp))   ; expro_ti = 0.0
       allocate(expro_vpol(nion,nexp)) ; expro_vpol = 0.0
       allocate(expro_vtor(nion,nexp)) ; expro_vtor = 0.0

       ! Derived quantities

       allocate(expro_bunit(nexp))        ; expro_bunit = 0.0
       allocate(expro_s(nexp))            ; expro_s = 0.0
       allocate(expro_drmaj(nexp))        ; expro_drmaj = 0.0
       allocate(expro_dzmag(nexp))        ; expro_dzmag = 0.0
       allocate(expro_sdelta(nexp))       ; expro_sdelta = 0.0
       allocate(expro_skappa(nexp))       ; expro_skappa = 0.0
       allocate(expro_szeta(nexp))        ; expro_szeta = 0.0
       allocate(expro_dlnnedr(nexp))      ; expro_dlnnedr = 0.0
       allocate(expro_dlntedr(nexp))      ; expro_dlntedr = 0.0
       allocate(expro_sdlnnedr(nexp))     ; expro_sdlnnedr = 0.0
       allocate(expro_sdlntedr(nexp))     ; expro_sdlntedr = 0.0
       allocate(expro_dlnptotdr(nexp))    ; expro_dlnptotdr = 0.0
       allocate(expro_w0p(nexp))          ; expro_w0p = 0.0
       allocate(expro_vol(nexp))          ; expro_vol = 0.0
       allocate(expro_volp(nexp))         ; expro_volp = 0.0
       allocate(expro_cs(nexp))           ; expro_cs = 0.0
       allocate(expro_rhos(nexp))         ; expro_rhos = 0.0
       allocate(expro_nuee(nexp))         ; expro_nuee = 0.0
       allocate(expro_ni_new(nexp))       ; expro_ni_new = 0.0
       allocate(expro_dlnnidr_new(nexp))  ; expro_dlnnidr_new = 0.0
       allocate(expro_sdlnnidr_new(nexp)) ; expro_sdlnnidr_new = 0.0
       allocate(expro_grad_r0(nexp))      ; expro_grad_r0 = 0.0
       allocate(expro_ave_grad_r(nexp))   ; expro_ave_grad_r = 0.0
       allocate(expro_bp0(nexp))          ; expro_bp0 = 0.0
       allocate(expro_bt0(nexp))          ; expro_bt0 = 0.0
       allocate(expro_ip(nexp))           ; expro_ip = 0.0
       allocate(expro_gamma_e(nexp))      ; expro_gamma_e = 0.0
       allocate(expro_gamma_p(nexp))      ; expro_gamma_p = 0.0
       allocate(expro_mach(nexp))         ; expro_mach = 0.0
       allocate(expro_thetascale(nexp))   ; expro_thetascale = 0.0
       allocate(expro_flow_beam(nexp))  ; expro_flow_beam = 0.0
       allocate(expro_flow_wall(nexp))  ; expro_flow_wall = 0.0
       allocate(expro_flow_mom(nexp))   ; expro_flow_mom = 0.0
       allocate(expro_pow_e(nexp))      ; expro_pow_e = 0.0
       allocate(expro_pow_i(nexp))      ; expro_pow_i = 0.0
       allocate(expro_pow_ei(nexp))     ; expro_pow_ei = 0.0
       allocate(expro_pow_e_aux(nexp))  ; expro_pow_e_aux = 0.0
       allocate(expro_pow_i_aux(nexp))  ; expro_pow_i_aux = 0.0
       allocate(expro_pow_e_fus(nexp))  ; expro_pow_e_fus = 0.0
       allocate(expro_pow_i_fus(nexp))  ; expro_pow_i_fus = 0.0
       allocate(expro_pow_e_sync(nexp)) ; expro_pow_e_sync = 0.0
       allocate(expro_pow_e_brem(nexp)) ; expro_pow_e_brem = 0.0
       allocate(expro_pow_e_line(nexp)) ; expro_pow_e_line = 0.0

       allocate(expro_dlnnidr(nion,nexp))  ; expro_dlnnidr = 0.0
       allocate(expro_dlntidr(nion,nexp))  ; expro_dlntidr = 0.0
       allocate(expro_sdlnnidr(nion,nexp)) ; expro_sdlnnidr = 0.0
       allocate(expro_sdlntidr(nion,nexp)) ; expro_sdlntidr = 0.0

    else

       deallocate(expro_mass) 
       deallocate(expro_z) 

       deallocate(expro_rho)
       deallocate(expro_rmin)
       deallocate(expro_q)
       deallocate(expro_polflux)
       deallocate(expro_w0)
       deallocate(expro_rmaj)
       deallocate(expro_zmag)
       deallocate(expro_kappa)
       deallocate(expro_delta)
       deallocate(expro_zeta)
       deallocate(expro_ne)
       deallocate(expro_te)
       deallocate(expro_ptot)
       deallocate(expro_johm)
       deallocate(expro_jbs)
       deallocate(expro_jrf)
       deallocate(expro_jnb)
       deallocate(expro_jbstor)
       deallocate(expro_sigmapar)
       deallocate(expro_z_eff)

       deallocate(expro_qohme)
       deallocate(expro_qbeame)
       deallocate(expro_qbeami)
       deallocate(expro_qrfe)
       deallocate(expro_qrfi)
       deallocate(expro_qfuse)
       deallocate(expro_qfusi)
       deallocate(expro_qbrem)
       deallocate(expro_qsync)
       deallocate(expro_qline)
       deallocate(expro_qei)
       deallocate(expro_qione)
       deallocate(expro_qioni)
       deallocate(expro_qcxi)
       deallocate(expro_qpar)
       deallocate(expro_qmom)

       deallocate(expro_ni,expro_ti,expro_vpol,expro_vtor)

       ! Derived

       deallocate(expro_bunit)
       deallocate(expro_s)
       deallocate(expro_drmaj)
       deallocate(expro_dzmag) 
       deallocate(expro_sdelta)       
       deallocate(expro_skappa)       
       deallocate(expro_szeta)        
       deallocate(expro_dlnnedr)      
       deallocate(expro_dlntedr)      
       deallocate(expro_sdlnnedr)      
       deallocate(expro_sdlntedr)      
       deallocate(expro_dlnptotdr)    
       deallocate(expro_w0p)          
       deallocate(expro_vol)          
       deallocate(expro_volp)         
       deallocate(expro_cs)           
       deallocate(expro_rhos)         
       deallocate(expro_nuee)
       deallocate(expro_ni_new)       
       deallocate(expro_dlnnidr_new)  
       deallocate(expro_sdlnnidr_new)
       deallocate(expro_grad_r0)      
       deallocate(expro_ave_grad_r)   
       deallocate(expro_bp0)          
       deallocate(expro_bt0)          
       deallocate(expro_ip)           
       deallocate(expro_gamma_e)      
       deallocate(expro_gamma_p)      
       deallocate(expro_mach)   
       deallocate(expro_thetascale)  
       deallocate(expro_flow_beam)
       deallocate(expro_flow_wall)
       deallocate(expro_flow_mom)
       deallocate(expro_pow_e)
       deallocate(expro_pow_i)
       deallocate(expro_pow_ei)
       deallocate(expro_pow_e_aux)
       deallocate(expro_pow_i_aux)
       deallocate(expro_pow_e_fus)
       deallocate(expro_pow_i_fus)
       deallocate(expro_pow_e_sync)
       deallocate(expro_pow_e_brem)
       deallocate(expro_pow_e_line)

       deallocate(expro_dlnnidr)
       deallocate(expro_dlntidr)  
       deallocate(expro_sdlnnidr) 
       deallocate(expro_sdlntidr)  

    endif

  end subroutine expro_init

  subroutine expro_read(thisinfile)

    implicit none

    character(len=*), intent(in) :: thisinfile 
    integer :: nexp,nion,ierr,i,nd
    character(len=70) :: ytag
    character(len=22) :: c

    ! ORDERING NOTE: nexp should appear before any profile arrays
    
    open(unit=1,file=trim(thisinfile),status='old')

    do

       read(1,'(a)',end=99) ytag

       if (index(ytag,'*original') > 0) then
          expro_head_original=ytag ; cycle
       else if (index(ytag,'*statefile') > 0) then
          expro_head_statefile=ytag ; cycle
       else if (index(ytag,'*gfile') > 0) then
          expro_head_gfile=ytag ; cycle
       else if (index(ytag,'*cerfile') > 0) then
          expro_head_cerfile=ytag ; cycle
       else if (index(ytag,'*vgen') > 0) then
          expro_head_vgen=ytag ; cycle
       else if (index(ytag,'*tgyro') > 0) then
          expro_head_tgyro=ytag ; cycle
       endif

       nd = scan(ytag,'|')
       if (nd == 0) then
          ! no units field so trim all whitespace
          c = trim(ytag(3:))
       else
          ! trim units and whitespace
          c = trim(ytag(3:nd-1))
       endif

       select case (c)
       case ('nexp')
          call expro_icomm(expro_n_exp) 
          nexp = expro_n_exp
       case ('nion')
          call expro_icomm(expro_n_ion)
          nion = expro_n_ion
          if (allocated(expro_rho)) call expro_init(0)
          call expro_init(1) 
       case ('name')
          call expro_tcomm(expro_name(1:nion),nion)
       case ('type')
          call expro_tcomm(expro_type(1:nion),nion)
       case ('masse')
          call expro_rcomm(expro_masse) 
       case ('mass')
          call expro_lcomm(expro_mass,nion)
       case ('ze')
          call expro_rcomm(expro_ze) 
       case ('z')
          call expro_lcomm(expro_z,nion)
       case ('torfluxa')
          call expro_rcomm(expro_torfluxa) 
       case ('rvbv')
          call expro_rcomm(expro_rvbv) 
       case ('ipa')
          call expro_rcomm(expro_ipa) 
       case ('rho')
          call expro_vcomm(expro_rho,nexp)  
       case ('rmin')
          call expro_vcomm(expro_rmin,nexp)  
       case ('polflux')
          call expro_vcomm(expro_polflux,nexp)  
       case ('q')
          call expro_vcomm(expro_q,nexp)  
       case ('w0')
          call expro_vcomm(expro_w0,nexp)  
       case ('rmaj')
          call expro_vcomm(expro_rmaj,nexp)  
       case ('zmag')
          call expro_vcomm(expro_zmag,nexp)  
       case ('kappa')
          call expro_vcomm(expro_kappa,nexp)  
       case ('delta')
          call expro_vcomm(expro_delta,nexp)  
       case ('zeta')
          call expro_vcomm(expro_zeta,nexp)  
       case ('ne')
          call expro_vcomm(expro_ne,nexp) 
       case ('te')
          call expro_vcomm(expro_te,nexp)  
       case ('ptot')
          call expro_vcomm(expro_ptot,nexp)  
       case ('johm')
          call expro_vcomm(expro_johm,nexp)  
       case ('jbs')
          call expro_vcomm(expro_jbs,nexp)  
       case ('jrf')
          call expro_vcomm(expro_jrf,nexp)  
       case ('jnb')
          call expro_vcomm(expro_jnb,nexp)  
       case ('jbstor')
          call expro_vcomm(expro_jbstor,nexp)  
       case ('sigmapar')
          call expro_vcomm(expro_sigmapar,nexp)  
       case ('z_eff')
          call expro_vcomm(expro_z_eff,nexp) 
       case ('ni')
          call expro_acomm(expro_ni(:,:),nion,nexp) 
       case ('ti')
          call expro_acomm(expro_ti(:,:),nion,nexp) 
       case ('vpol')
          call expro_acomm(expro_vpol(:,:),nion,nexp) 
       case ('vtor')
          call expro_acomm(expro_vtor(:,:),nion,nexp) 
       case ('qohme')
          call expro_vcomm(expro_qohme,nexp) 
       case ('qbeame')
          call expro_vcomm(expro_qbeame,nexp) 
       case ('qbeami')
          call expro_vcomm(expro_qbeami,nexp) 
       case ('qrfe')
          call expro_vcomm(expro_qrfe,nexp) 
       case ('qrfi')
          call expro_vcomm(expro_qrfi,nexp) 
       case ('qfuse')
          call expro_vcomm(expro_qfuse,nexp) 
       case ('qfusi')
          call expro_vcomm(expro_qfusi,nexp) 
       case ('qbrem')
          call expro_vcomm(expro_qbrem,nexp) 
       case ('qsync')
          call expro_vcomm(expro_qsync,nexp) 
       case ('qline')
          call expro_vcomm(expro_qline,nexp) 
       case ('qei')
          call expro_vcomm(expro_qei,nexp) 
       case ('qione')
          call expro_vcomm(expro_qione,nexp) 
       case ('qioni')
          call expro_vcomm(expro_qioni,nexp) 
       case ('qcxi')
          call expro_vcomm(expro_qcxi,nexp) 
       case ('qpar')
          call expro_vcomm(expro_qpar,nexp) 
       case ('qmom')
          call expro_vcomm(expro_qmom,nexp) 
       end select

    enddo

99  close(1)

    ! ** input.gacode.geo **

    nexp = expro_n_exp
    open(unit=1,file=trim(thisinfile)//'.geo',status='old',iostat=ierr)
    if (ierr == 0) then
       call expro_skip_header(1)
       call expro_icomm(expro_nfourier)
       if (allocated(expro_geo)) deallocate(expro_geo)
       if (allocated(expro_dgeo)) deallocate(expro_dgeo)
       allocate(expro_geo(4,0:expro_nfourier,nexp)) ; expro_geo(:,:,:)=0.0
       allocate(expro_dgeo(4,0:expro_nfourier,nexp)) ; expro_dgeo(:,:,:)=0.0
       do i=1,nexp
          call expro_scomm(expro_geo(:,:,i),4*(expro_nfourier+1))
       enddo
    else
       expro_nfourier = -1
    endif
    close(1)

    ! BCAST HERE

    if (expro_ctrl_n_ion <= expro_n_ion) then
       call expro_compute_derived
    else
       expro_error = 1
    endif

  end subroutine expro_read

  subroutine expro_write(thisinfile)

    implicit none

    integer :: i,nexp,nion
    character(len=*), intent(in) :: thisinfile 

    nexp = expro_n_exp
    nion = expro_n_ion

    ! Write header
    open(unit=1,file=trim(thisinfile),status='replace')
    write(1,'(a)') expro_head_original
    write(1,'(a)') expro_head_statefile 
    write(1,'(a)') expro_head_gfile
    write(1,'(a)') expro_head_cerfile
    write(1,'(a)') expro_head_vgen
    write(1,'(a)') expro_head_tgyro
    write(1,'(a)') '#'

    ! Write data
    write(1,'(a)') ident//'nexp' ; write(1,'(i0)') nexp
    write(1,'(a)') ident//'nion' ; write(1,'(i0)') nion
    write(1,'(a)') ident//'name' ; write(1,'(20(a,1x))') (trim(expro_name(i)),i=1,nion)
    write(1,'(a)') ident//'type' ; write(1,'(20(a,1x))') (trim(expro_type(i)),i=1,nion)
    write(1,'(a)') ident//'masse'; write(1,30) expro_masse
    write(1,'(a)') ident//'mass' ; write(1,40) expro_mass
    write(1,'(a)') ident//'ze'   ; write(1,30) expro_ze
    write(1,'(a)') ident//'z'    ; write(1,40) expro_z

    ! Write vector/array data, skipping objects that are 0.0
    call expro_writes(expro_torfluxa,'torfluxa','Wb/radian')
    call expro_writes(expro_rvbv,'rvbv','Tm')
    call expro_writes(expro_ipa,'ipa','MA')
    call expro_writev(expro_rho,nexp,'rho','-')
    call expro_writev(expro_rmin,nexp,'rmin','m')
    call expro_writev(expro_polflux,nexp,'polflux','Wb/radian')
    call expro_writev(expro_q,nexp,'q','-')
    call expro_writev(expro_w0,nexp,'w0','rad/s')
    call expro_writev(expro_rmaj,nexp,'rmaj','m')
    call expro_writev(expro_zmag,nexp,'zmag','m')
    call expro_writev(expro_kappa,nexp,'kappa','-')
    call expro_writev(expro_delta,nexp,'delta','-')
    call expro_writev(expro_zeta,nexp,'zeta','-')
    call expro_writev(expro_ne,nexp,'ne','10^19/m^3')
    call expro_writea(expro_ni(:,:),nion,nexp,'ni','10^19/m^3')
    call expro_writev(expro_te,nexp,'te','keV')
    call expro_writea(expro_ti(:,:),nion,nexp,'ti','keV')
    call expro_writev(expro_ptot,nexp,'ptot','Pa')
    call expro_writev(expro_johm,nexp,'johm','MA/m^2')
    call expro_writev(expro_jbs,nexp,'jbs','MA/m^2')
    call expro_writev(expro_jrf,nexp,'jrf','MA/m^2')
    call expro_writev(expro_jnb,nexp,'jnb','MA/m^2')
    call expro_writev(expro_jbstor,nexp,'jbstor','MA/m^2')
    call expro_writev(expro_sigmapar,nexp,'sigmapar','MS/m')
    call expro_writev(expro_z_eff,nexp,'z_eff','-')
    call expro_writea(expro_vpol(:,:),nion,nexp,'vpol','m/s')
    call expro_writea(expro_vtor(:,:),nion,nexp,'vtor','m/s')
    call expro_writev(expro_qohme,nexp,'qohme','MW/m^3')
    call expro_writev(expro_qbeame,nexp,'qbeame','MW/m^3')
    call expro_writev(expro_qbeami,nexp,'qbeami','MW/m^3')
    call expro_writev(expro_qrfe,nexp,'qrfe','MW/m^3')
    call expro_writev(expro_qrfi,nexp,'qrfi','MW/m^3')
    call expro_writev(expro_qfuse,nexp,'qfuse','MW/m^3')
    call expro_writev(expro_qfusi,nexp,'qfusi','MW/m^3')
    call expro_writev(expro_qbrem,nexp,'qbrem','MW/m^3')
    call expro_writev(expro_qsync,nexp,'qsync','MW/m^3')
    call expro_writev(expro_qline,nexp,'qline','MW/m^3')
    call expro_writev(expro_qei,nexp,'qei','MW/m^3')
    call expro_writev(expro_qione,nexp,'qione','MW/m^3')
    call expro_writev(expro_qioni,nexp,'qioni','MW/m^3')
    call expro_writev(expro_qcxi,nexp,'qcxi','MW/m^3')
    call expro_writev(expro_qpar,nexp,'qpar','MW/m^3')
    call expro_writev(expro_qmom,nexp,'qmom','MW/m^3')

    close(1)

30  format(1pe14.7)
40  format(10(1pe14.7))

  end subroutine expro_write

end module expro
