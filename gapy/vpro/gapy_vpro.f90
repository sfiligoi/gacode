module vpro

  ! Fundamental input
  integer, parameter :: ntag = 40
  character(len=12) :: input = 'input.gacode'
  
  character(len=2) :: ident='# '
  character(len=20), dimension(ntag) :: tag = (/&
       'nexp      ',& !1
       'nion      ',& !2
       'mass      ',& !3
       'z         ',& !4
       'bt_exp    ',& !5
       'arho_exp  ',& !6
       'rho       ',& !7
       'rmin      ',& !8
       'polflux   ',& !9
       'q         ',& !10
       'w0        ',& !11
       'rmaj      ',& !12
       'zmag      ',& !13
       'kappa     ',& !14
       'delta     ',& !15
       'zeta      ',& !16
       'ne        ',& !17
       'Te        ',& !18
       'ptot      ',& !19
       'z_eff     ',& !20
       'flow_beam ',& !21
       'flow_wall ',& !22
       'flow_mom  ',& !23
       'pow_e     ',& !24
       'pow_i     ',& !25
       'pow_ei    ',& !26
       'pow_e_aux ',& !27
       'pow_i_aux ',& !28
       'pow_e_fus ',& !29
       'pow_i_fus ',& !30
       'pow_e_sync',& !31
       'pow_e_brem',& !32
       'pow_e_line',& !33
       'sbeame    ',& !34
       'sbcx      ',& !35
       'sscxl     ',& !36
       'ni        ',& !37
       'ti        ',& !38
       'vpol      ',& !39
       'vtor      '&  !40
       /)

  integer :: &
       EXPRO_n_exp,&
       EXPRO_n_ion

  double precision, dimension(:), allocatable :: &
       EXPRO_mass,&
       EXPRO_z

  double precision :: &
       EXPRO_b_ref,&
       EXPRO_arho

  double precision, dimension(:), allocatable :: &
       EXPRO_rho,&
       EXPRO_rmin,&
       EXPRO_polflux,&
       EXPRO_q,&
       EXPRO_w0,&
       EXPRO_rmaj,&
       EXPRO_zmag,&
       EXPRO_kappa,&
       EXPRO_delta,&
       EXPRO_zeta,&
       EXPRO_ne,&
       EXPRO_te,&
       EXPRO_ptot,&
       EXPRO_z_eff,&
       EXPRO_flow_beam,&
       EXPRO_flow_wall,&
       EXPRO_flow_mom,&
       EXPRO_pow_e,&
       EXPRO_pow_i,&
       EXPRO_pow_ei,&
       EXPRO_pow_e_aux,&
       EXPRO_pow_i_aux,&
       EXPRO_pow_e_fus,&
       EXPRO_pow_i_fus,&
       EXPRO_pow_e_sync,&
       EXPRO_pow_e_brem,&
       EXPRO_pow_e_line,&
       EXPRO_sbeame,&
       EXPRO_sbcx,&
       EXPRO_sscxl

  double precision, dimension(:,:), allocatable :: &
       EXPRO_ni,&
       EXPRO_ti,&
       EXPRO_vpol,&
       EXPRO_vtor
  
  ! Derived quantities

  double precision, dimension(:), allocatable :: &
       EXPRO_bunit,&
       EXPRO_s,&
       EXPRO_drmaj,&
       EXPRO_dzmag,&
       EXPRO_sdelta,&
       EXPRO_skappa,&
       EXPRO_szeta,&
       EXPRO_dlnnedr,&
       EXPRO_dlntedr,&
       EXPRO_sdlnnedr,&
       EXPRO_sdlntedr,&
       EXPRO_dlnptotdr,&
       EXPRO_w0p,&
       EXPRO_vol,&
       EXPRO_volp,&
       EXPRO_cs,&
       EXPRO_rhos,&
       EXPRO_nuee,&
       EXPRO_ni_new,&
       EXPRO_dlnnidr_new,&
       EXPRO_sdlnnidr_new,&
       EXPRO_grad_r0,&
       EXPRO_ave_grad_r,&
       EXPRO_drdrho,&
       EXPRO_bp0,&
       EXPRO_bt0,&
       EXPRO_ip,&
       EXPRO_gamma_e,&
       EXPRO_gamma_p,&
       EXPRO_mach,&
       EXPRO_thetascale

  double precision, dimension(:,:), allocatable :: &
       EXPRO_dlnnidr,&
       EXPRO_dlntidr,&
       EXPRO_sdlnnidr,&
       EXPRO_sdlntidr

  ! input.profiles.geo dimension and arrays

  integer :: EXPRO_nfourier
  double precision, dimension(:,:,:),allocatable :: EXPRO_geo
  double precision, dimension(:,:,:),allocatable :: EXPRO_dgeo

  ! Field orientation parameters

  integer :: EXPRO_signb
  integer :: EXPRO_signq

  ! Control parameters (force nonsensical default -1 for usage check)

  integer :: EXPRO_ctrl_n_ion 
  integer :: EXPRO_ctrl_quasineutral_flag 
  integer :: EXPRO_ctrl_numeq_flag

  ! *** locsim variables ***
  
  integer :: n_species_exp
  
  double precision, parameter :: mass_deuterium  = 3.3452
  double precision, parameter :: temp_norm_fac   = 1602.2
  double precision, parameter :: charge_norm_fac = 1.6022

  double precision, dimension(:), allocatable :: rmin_exp

  double precision, dimension(:,:), allocatable :: temp_exp
  double precision, dimension(:,:), allocatable :: dens_exp
  double precision, dimension(:,:), allocatable :: dlntdr_exp
  double precision, dimension(:,:), allocatable :: dlnndr_exp
  double precision, dimension(:,:), allocatable :: sdlntdr_exp
  double precision, dimension(:,:), allocatable :: sdlnndr_exp

  double precision, dimension(:), allocatable :: gamma_e_exp
  double precision, dimension(:), allocatable :: gamma_p_exp
  double precision, dimension(:), allocatable :: mach_exp
  
  double precision, dimension(:,:,:), allocatable :: geo_yin_exp

  ! Local values

  double precision :: shift_loc
  double precision :: q_loc
  double precision :: s_loc
  double precision :: kappa_loc
  double precision :: delta_loc
  double precision :: zeta_loc
  double precision :: s_kappa_loc
  double precision :: s_delta_loc
  double precision :: s_zeta_loc
  double precision :: zmag_loc
  double precision :: dzmag_loc
  double precision :: gamma_e_loc
  double precision :: gamma_p_loc
  double precision :: mach_loc
  double precision :: rmin_loc
  double precision :: rmaj_loc
  double precision :: rhos_loc
  double precision :: z_eff_loc
  double precision :: b_unit_loc
  double precision :: rho_norm_loc
  double precision :: psi_norm_loc
  double precision :: psi_a_loc
  double precision :: cs_loc
  double precision :: betae_loc
  double precision :: beta_star_loc

  double precision, dimension(9) :: mass_loc
  double precision, dimension(9) :: z_loc
  double precision, dimension(9) :: dens_loc
  double precision, dimension(9) :: temp_loc
  double precision, dimension(9) :: dlnndr_loc
  double precision, dimension(9) :: dlntdr_loc
  double precision, dimension(9) :: sdlnndr_loc
  double precision, dimension(9) :: sdlntdr_loc

  integer :: geo_ny_loc
  double precision, dimension(:,:), allocatable :: geo_yin_loc
  
contains

  subroutine vpro_read

    implicit none

    character(len=22) :: ytag

    ! NOTE: nexp should appear before any profile arrays

    open(unit=1,file=input,status='old')

    do 

       read(1,'(a)',end=99) ytag

       select case (trim(ytag(3:22)))
       case('nexp')
          read(1,*) EXPRO_n_exp
       case('nion')
          read(1,*) EXPRO_n_ion
          if (allocated(EXPRO_rho)) call vpro_init(0)
          call vpro_init(1)
       case('mass')
          read(1,*) EXPRO_mass
       case('z')
          read(1,*) EXPRO_z
       case('bt_exp')
          read(1,30) EXPRO_b_ref
       case('arho_exp')
          read(1,30) EXPRO_arho
       case ('rho')
          read(1,30) EXPRO_rho 
       case ('rmin')
          read(1,30) EXPRO_rmin 
       case ('polflux')
          read(1,30) EXPRO_polflux 
       case ('q')
          read(1,30) EXPRO_q 
       case ('w0')
          read(1,30) EXPRO_w0 
       case ('rmaj')
          read(1,30) EXPRO_rmaj 
       case ('zmag')
          read(1,30) EXPRO_zmag 
       case ('kappa')
          read(1,30) EXPRO_kappa 
       case ('delta')
          read(1,30) EXPRO_delta 
       case ('zeta')
          read(1,30) EXPRO_zeta 
       case ('ne')
          read(1,30) EXPRO_ne
       case ('Te')
          read(1,30) EXPRO_te 
       case ('ptot')
          read(1,30) EXPRO_ptot 
       case ('z_eff')
          read(1,30) EXPRO_z_eff 
       case ('ni')
          read(1,30) EXPRO_ni(:,:)
       case ('ti')
          read(1,30) EXPRO_ti(:,:)
       case ('vpol')
          read(1,30) EXPRO_vpol(:,:)
       case ('vtor')
          read(1,30) EXPRO_vtor(:,:)
       case ('flow_beam')
          read(1,30) EXPRO_flow_beam 
       case ('flow_wall')
          read(1,30) EXPRO_flow_wall 
       case ('flow_mom')
          read(1,30) EXPRO_flow_mom 
       case ('pow_e')
          read(1,30) EXPRO_pow_e 
       case ('pow_i')
          read(1,30) EXPRO_pow_i 
       case ('pow_ei')
          read(1,30) EXPRO_pow_ei 
       case ('pow_e_aux')
          read(1,30) EXPRO_pow_e_aux 
       case ('pow_i_aux')
          read(1,30) EXPRO_pow_i_aux 
       case ('pow_e_fus')
          read(1,30) EXPRO_pow_e_fus 
       case ('pow_i_fus')
          read(1,30) EXPRO_pow_i_fus 
       case ('pow_e_sync')
          read(1,30) EXPRO_pow_e_sync
       case ('pow_e_brem')
          read(1,30) EXPRO_pow_e_brem 
       case ('pow_e_line')
          read(1,30) EXPRO_pow_e_line 
       case ('sbeame')
          read(1,30) EXPRO_sbeame 
       case ('sbcx')
          read(1,30) EXPRO_sbcx 
       case ('sscxl')
          read(1,30) EXPRO_sscxl 
       end select

    enddo

99  close(1)

30  format(1pe12.5)

  end subroutine vpro_read

  subroutine vpro_init(flag)

    implicit none

    integer :: nexp,nion
    integer, intent(in) :: flag

    if (flag == 1) then

       nexp = EXPRO_n_exp
       nion = EXPRO_n_ion

       allocate(EXPRO_mass(nion)) ; EXPRO_mass = 1.0
       allocate(EXPRO_z(nion))    ; EXPRO_z = 1.0
       
       allocate(EXPRO_rho(nexp))
       allocate(EXPRO_rmin(nexp))
       allocate(EXPRO_q(nexp))
       allocate(EXPRO_polflux(nexp))
       allocate(EXPRO_w0(nexp))
       allocate(EXPRO_rmaj(nexp))
       allocate(EXPRO_zmag(nexp))
       allocate(EXPRO_kappa(nexp))
       allocate(EXPRO_delta(nexp))
       allocate(EXPRO_zeta(nexp))
       allocate(EXPRO_ne(nexp))
       allocate(EXPRO_te(nexp))
       allocate(EXPRO_ptot(nexp))
       allocate(EXPRO_z_eff(nexp))

       allocate(EXPRO_flow_beam(nexp))
       allocate(EXPRO_flow_wall(nexp))
       allocate(EXPRO_flow_mom(nexp))
       allocate(EXPRO_pow_e(nexp))
       allocate(EXPRO_pow_i(nexp))
       allocate(EXPRO_pow_ei(nexp))
       allocate(EXPRO_pow_e_aux(nexp))
       allocate(EXPRO_pow_i_aux(nexp))
       allocate(EXPRO_pow_e_fus(nexp))
       allocate(EXPRO_pow_i_fus(nexp))
       allocate(EXPRO_pow_e_sync(nexp))
       allocate(EXPRO_pow_e_brem(nexp))
       allocate(EXPRO_pow_e_line(nexp))
       allocate(EXPRO_sbeame(nexp))
       allocate(EXPRO_sbcx(nexp))
       allocate(EXPRO_sscxl(nexp))

       EXPRO_rho     = 0.0
       EXPRO_rmin    = 0.0
       EXPRO_q       = 0.0
       EXPRO_polflux = 0.0
       EXPRO_w0      = 0.0
       EXPRO_rmaj    = 0.0
       EXPRO_zmag    = 0.0
       EXPRO_kappa   = 0.0
       EXPRO_delta   = 0.0
       EXPRO_zeta    = 0.0
       EXPRO_ne      = 0.0
       EXPRO_te      = 0.0
       EXPRO_ptot    = 0.0
       EXPRO_z_eff   = 0.0

       EXPRO_flow_beam = 0.0
       EXPRO_flow_wall = 0.0
       EXPRO_flow_mom  = 0.0
       EXPRO_pow_e     = 0.0
       EXPRO_pow_i     = 0.0
       EXPRO_pow_ei    = 0.0
       EXPRO_pow_e_aux = 0.0
       EXPRO_pow_i_aux = 0.0
       EXPRO_pow_e_fus = 0.0
       EXPRO_pow_i_fus = 0.0
       EXPRO_pow_e_sync = 0.0
       EXPRO_pow_e_brem = 0.0
       EXPRO_pow_e_line = 0.0
       EXPRO_sbeame = 0.0
       EXPRO_sbcx   = 0.0
       EXPRO_sscxl  = 0.0

       allocate(EXPRO_ni(nion,nexp))   ; EXPRO_ni = 0.0
       allocate(EXPRO_ti(nion,nexp))   ; EXPRO_ti = 0.0
       allocate(EXPRO_vpol(nion,nexp)) ; EXPRO_vpol = 0.0
       allocate(EXPRO_vtor(nion,nexp)) ; EXPRO_vtor = 0.0
       
       ! Derived quantities

       allocate(EXPRO_bunit(nexp))        ; EXPRO_bunit = 0.0
       allocate(EXPRO_s(nexp))            ; EXPRO_s = 0.0
       allocate(EXPRO_drmaj(nexp))        ; EXPRO_s = 0.0
       allocate(EXPRO_dzmag(nexp))        ; EXPRO_s = 0.0
       allocate(EXPRO_sdelta(nexp))       ; EXPRO_s = 0.0
       allocate(EXPRO_skappa(nexp))       ; EXPRO_s = 0.0
       allocate(EXPRO_szeta(nexp))        ; EXPRO_s = 0.0
       allocate(EXPRO_dlnnedr(nexp))      ; EXPRO_s = 0.0
       allocate(EXPRO_dlntedr(nexp))      ; EXPRO_s = 0.0
       allocate(EXPRO_sdlnnedr(nexp))     ; EXPRO_s = 0.0
       allocate(EXPRO_sdlntedr(nexp))     ; EXPRO_s = 0.0
       allocate(EXPRO_dlnptotdr(nexp))    ; EXPRO_s = 0.0
       allocate(EXPRO_w0p(nexp))          ; EXPRO_s = 0.0
       allocate(EXPRO_vol(nexp))          ; EXPRO_s = 0.0
       allocate(EXPRO_volp(nexp))         ; EXPRO_s = 0.0
       allocate(EXPRO_cs(nexp))           ; EXPRO_s = 0.0
       allocate(EXPRO_rhos(nexp))         ; EXPRO_s = 0.0
       allocate(EXPRO_nuee(nexp))         ; EXPRO_nuee = 0.0
       allocate(EXPRO_ni_new(nexp))       ; EXPRO_s = 0.0
       allocate(EXPRO_dlnnidr_new(nexp))  ; EXPRO_dlnnidr_new = 0.0
       allocate(EXPRO_sdlnnidr_new(nexp)) ; EXPRO_sdlnnidr_new = 0.0
       allocate(EXPRO_grad_r0(nexp))      ; EXPRO_s = 0.0
       allocate(EXPRO_ave_grad_r(nexp))   ; EXPRO_s = 0.0
       allocate(EXPRO_drdrho(nexp))       ; EXPRO_s = 0.0
       allocate(EXPRO_bp0(nexp))          ; EXPRO_s = 0.0
       allocate(EXPRO_bt0(nexp))          ; EXPRO_s = 0.0
       allocate(EXPRO_ip(nexp))           ; EXPRO_s = 0.0
       allocate(EXPRO_gamma_e(nexp))      ; EXPRO_s = 0.0
       allocate(EXPRO_gamma_p(nexp))      ; EXPRO_s = 0.0
       allocate(EXPRO_mach(nexp))         ; EXPRO_s = 0.0
       allocate(EXPRO_thetascale(nexp))   ; EXPRO_s = 0.0

       allocate(EXPRO_dlnnidr(nion,nexp))  ; EXPRO_dlnnidr = 0.0
       allocate(EXPRO_dlntidr(nion,nexp))  ; EXPRO_dlntidr = 0.0
       allocate(EXPRO_sdlnnidr(nion,nexp)) ; EXPRO_sdlnnidr = 0.0
       allocate(EXPRO_sdlntidr(nion,nexp)) ; EXPRO_sdlntidr = 0.0
   
    else

       deallocate(EXPRO_mass) 
       deallocate(EXPRO_z) 

       deallocate(EXPRO_rho)
       deallocate(EXPRO_rmin)
       deallocate(EXPRO_q)
       deallocate(EXPRO_polflux)
       deallocate(EXPRO_w0)
       deallocate(EXPRO_rmaj)
       deallocate(EXPRO_zmag)
       deallocate(EXPRO_kappa)
       deallocate(EXPRO_delta)
       deallocate(EXPRO_zeta)
       deallocate(EXPRO_ne)
       deallocate(EXPRO_te)
       deallocate(EXPRO_ptot)
       deallocate(EXPRO_z_eff)

       deallocate(EXPRO_flow_beam)
       deallocate(EXPRO_flow_wall)
       deallocate(EXPRO_flow_mom)
       deallocate(EXPRO_pow_e)
       deallocate(EXPRO_pow_i)
       deallocate(EXPRO_pow_ei)
       deallocate(EXPRO_pow_e_aux)
       deallocate(EXPRO_pow_i_aux)
       deallocate(EXPRO_pow_e_fus)
       deallocate(EXPRO_pow_i_fus)
       deallocate(EXPRO_pow_e_sync)
       deallocate(EXPRO_pow_e_brem)
       deallocate(EXPRO_pow_e_line)
       deallocate(EXPRO_sbeame)
       deallocate(EXPRO_sbcx)
       deallocate(EXPRO_sscxl)

       deallocate(EXPRO_ni,EXPRO_ti,EXPRO_vpol,EXPRO_vtor)

       ! Derived

       deallocate(EXPRO_bunit)
       deallocate(EXPRO_s)
       deallocate(EXPRO_drmaj)
       deallocate(EXPRO_dzmag) 
       deallocate(EXPRO_sdelta)       
       deallocate(EXPRO_skappa)       
       deallocate(EXPRO_szeta)        
       deallocate(EXPRO_dlnnedr)      
       deallocate(EXPRO_dlntedr)      
       deallocate(EXPRO_sdlnnedr)      
       deallocate(EXPRO_sdlntedr)      
       deallocate(EXPRO_dlnptotdr)    
       deallocate(EXPRO_w0p)          
       deallocate(EXPRO_vol)          
       deallocate(EXPRO_volp)         
       deallocate(EXPRO_cs)           
       deallocate(EXPRO_rhos)         
       deallocate(EXPRO_nuee)
       deallocate(EXPRO_ni_new)       
       deallocate(EXPRO_dlnnidr_new)  
       deallocate(EXPRO_sdlnnidr_new)
       deallocate(EXPRO_grad_r0)      
       deallocate(EXPRO_ave_grad_r)   
       deallocate(EXPRO_drdrho)       
       deallocate(EXPRO_bp0)          
       deallocate(EXPRO_bt0)          
       deallocate(EXPRO_ip)           
       deallocate(EXPRO_gamma_e)      
       deallocate(EXPRO_gamma_p)      
       deallocate(EXPRO_mach)   
       deallocate(EXPRO_thetascale)  

       deallocate(EXPRO_dlnnidr)
       deallocate(EXPRO_dlntidr)  
       deallocate(EXPRO_sdlnnidr) 
       deallocate(EXPRO_sdlntidr)  

    endif

  end subroutine vpro_init

  subroutine vpro_write

    implicit none

    integer :: nexp,nion

    nexp = EXPRO_n_exp
    nion = EXPRO_n_ion

    open(unit=1,file=input,status='replace')

    write(1,20) ident//tag(1)  ; write(1,'(i0)') nexp
    write(1,20) ident//tag(2)  ; write(1,'(i0)') nion
    write(1,20) ident//tag(3)  ; write(1,30) EXPRO_mass
    write(1,20) ident//tag(4)  ; write(1,30) EXPRO_z
    write(1,20) ident//tag(5)  ; write(1,30) EXPRO_b_ref
    write(1,20) ident//tag(6)  ; write(1,30) EXPRO_arho
    write(1,20) ident//tag(7)  ; call vpro_writev(EXPRO_rho,nexp)
    write(1,20) ident//tag(8)  ; call vpro_writev(EXPRO_rmin,nexp)
    write(1,20) ident//tag(9)  ; call vpro_writev(EXPRO_polflux,nexp)
    write(1,20) ident//tag(10) ; call vpro_writev(EXPRO_q,nexp)
    write(1,20) ident//tag(11) ; call vpro_writev(EXPRO_w0,nexp)
    write(1,20) ident//tag(12)
    write(1,30) EXPRO_rmaj

    write(1,20) ident//tag(13)
    write(1,30) EXPRO_zmag

    write(1,20) ident//tag(14)
    write(1,30) EXPRO_kappa

    write(1,20) ident//tag(15)
    write(1,30) EXPRO_delta

    write(1,20) ident//tag(16)
    write(1,30) EXPRO_zeta

    write(1,20) ident//tag(17)
    write(1,30) EXPRO_ne

    write(1,20) ident//tag(18)
    write(1,30) EXPRO_te

    write(1,20) ident//tag(19)
    write(1,30) EXPRO_ptot

    write(1,20) ident//tag(20)
    write(1,30) EXPRO_z_eff

    write(1,20) ident//tag(21)
    write(1,30) EXPRO_flow_beam

    write(1,20) ident//tag(22)
    write(1,30) EXPRO_flow_wall

    write(1,20) ident//tag(23)
    write(1,30) EXPRO_flow_mom

    write(1,20) ident//tag(24)
    write(1,30) EXPRO_pow_e

    write(1,20) ident//tag(25)
    write(1,30) EXPRO_pow_i

    write(1,20) ident//tag(26)
    write(1,30) EXPRO_pow_ei

    write(1,20) ident//tag(27)
    write(1,30) EXPRO_pow_e_aux

    write(1,20) ident//tag(28)
    write(1,30) EXPRO_pow_i_aux

    write(1,20) ident//tag(29)
    write(1,30) EXPRO_pow_e_fus

    write(1,20) ident//tag(30)
    write(1,30) EXPRO_pow_i_fus

    write(1,20) ident//tag(31)
    write(1,30) EXPRO_pow_e_sync

    write(1,20) ident//tag(32)
    write(1,30) EXPRO_pow_e_brem

    write(1,20) ident//tag(33)
    write(1,30) EXPRO_pow_e_line

    write(1,20) ident//tag(34)
    write(1,30) EXPRO_sbeame

    write(1,20) ident//tag(35)
    write(1,30) EXPRO_sbcx

    write(1,20) ident//tag(36)
    write(1,30) EXPRO_sscxl

    write(1,20) ident//tag(37)
    write(1,30) EXPRO_ni(:,:)

    write(1,20) ident//tag(38)
    write(1,30) EXPRO_ti(:,:)

    write(1,20) ident//tag(39)
    write(1,30) EXPRO_vpol(:,:)

    write(1,20) ident//tag(40)
    write(1,30) EXPRO_vtor(:,:)

    close(1)

20  format(a)
30  format(1pe12.5)

  end subroutine vpro_write

  subroutine vpro_read_legacy

    implicit none

    integer :: i
    integer :: nexp,nion
    character(len=99) :: line
    double precision :: x(5)

    open(unit=1,file='input.profiles',status='old')
    do while (line(1:2) /= '#r')
       read(1,'(a)') line
       if (line(1:5) == 'N_EXP') then
          read(line(7:),*) EXPRO_n_exp
       endif
       if (line(1:5) == 'N_ION') then
          read(line(7:),*) EXPRO_n_ion
       endif
       if (line(1:6) == 'BT_EXP') then
          read(line(8:),*) EXPRO_b_ref
       endif
       if (line(1:8) == 'ARHO_EXP') then
          read(line(10:),*) EXPRO_arho
       endif
    enddo

    call vpro_init(1)

    nexp = EXPRO_n_exp
    nion = EXPRO_n_ion

    ! 1
    do i=1,nexp
       read(1,*) x
       EXPRO_rho(i)     = x(1)
       EXPRO_rmin(i)    = x(2)
       EXPRO_polflux(i) = x(3)
       EXPRO_q(i)       = x(4)
       EXPRO_w0(i)      = x(5)
    enddo

    read(1,'(a)') line
    read(1,'(a)') line

    ! 2
    do i=1,nexp
       read(1,*) x
       EXPRO_rmaj(i)  = x(1)
       EXPRO_zmag(i)  = x(2)
       EXPRO_kappa(i) = x(3)
       EXPRO_delta(i) = x(4)
       EXPRO_zeta(i)  = x(5)
    enddo

    read(1,'(a)') line
    read(1,'(a)') line

    ! 3
    do i=1,nexp
       read(1,*) x
       EXPRO_ne(i)    = x(1)
       EXPRO_te(i)    = x(2)
       EXPRO_ptot(i)  = x(3)
       EXPRO_z_eff(i) = x(4)
    enddo

    read(1,'(a)') line
    read(1,'(a)') line

    ! 4
    do i=1,nexp
       read(1,*) x
       EXPRO_ni(1:nion,i) = x(1:nion)
    enddo

    read(1,'(a)') line
    read(1,'(a)') line

    ! 5 (assume < 6 ions, so skip)
    do i=1,nexp
       read(1,*) x
    enddo

    read(1,'(a)') line
    read(1,'(a)') line

    ! 6
    do i=1,nexp
       read(1,*) x
       EXPRO_ti(1:nion,i) = x(1:nion)
    enddo

    read(1,'(a)') line
    read(1,'(a)') line

    ! 7 (assume < 6 ions, so skip)
    do i=1,nexp
       read(1,*) x
    enddo
    close(1)

  end subroutine vpro_read_legacy

  subroutine vpro_writev(x,n)

    implicit none

    integer, intent(in) :: n
    double precision, intent(in), dimension(n) :: x
    integer :: i

    do i=1,n
       write(1,10) x(i),i
    enddo

10  format(1pe12.5,1x,i3)

  end subroutine vpro_writev

end module vpro
