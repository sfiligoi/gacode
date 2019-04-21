module vpro

  ! Fundamental input
  integer, parameter, private :: ntag = 40
  character*12, private :: infile = 'input.gacode'

  character*2, private :: ident='# '
  character*20, dimension(ntag), private :: tag = (/&
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
       expro_n_exp,&
       expro_n_ion

  double precision, dimension(:), allocatable :: &
       expro_mass,&
       expro_z

  character*7, dimension(:), allocatable :: &
       expro_name,&
       expro_type

  double precision :: &
       expro_b_ref,&
       expro_arho,&
       expro_rvbv,&
       expro_ip_exp

  double precision, dimension(:), allocatable :: &
       expro_rho,&
       expro_rmin,&
       expro_polflux,&
       expro_q,&
       expro_w0,&
       expro_rmaj,&
       expro_zmag,&
       expro_kappa,&
       expro_delta,&
       expro_zeta,&
       expro_ne,&
       expro_te,&
       expro_ptot,&
       expro_z_eff,&
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
       expro_pow_e_line,&
       expro_sbeame,&
       expro_sbcx,&
       expro_sscxl

  double precision, dimension(:,:), allocatable :: &
       expro_ni,&
       expro_ti,&
       expro_vpol,&
       expro_vtor

  ! Derived quantities

  double precision, dimension(:), allocatable :: &
       expro_bunit,&
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
       expro_drdrho,&
       expro_bp0,&
       expro_bt0,&
       expro_ip,&
       expro_gamma_e,&
       expro_gamma_p,&
       expro_mach,&
       expro_thetascale

  double precision, dimension(:,:), allocatable :: &
       expro_dlnnidr,&
       expro_dlntidr,&
       expro_sdlnnidr,&
       expro_sdlntidr

  ! input.profiles.geo dimension and arrays

  integer :: expro_nfourier
  double precision, dimension(:,:,:),allocatable :: expro_geo
  double precision, dimension(:,:,:),allocatable :: expro_dgeo

  ! Field orientation parameters

  integer :: expro_signb
  integer :: expro_signq

  ! Control parameters

  integer :: expro_ctrl_n_ion = -1
  integer :: expro_ctrl_quasineutral_flag 
  integer :: expro_ctrl_numeq_flag

contains

  subroutine vpro_init(flag)

    implicit none

    integer :: nexp,nion
    integer, intent(in) :: flag

    if (flag == 1) then

       nexp = expro_n_exp
       nion = expro_n_ion

       allocate(expro_mass(nion))    ; expro_mass = 1.0
       allocate(expro_z(nion))       ; expro_z = 1.0
       allocate(expro_type(nion))    ; expro_type = ' '
       allocate(expro_name(nion))    ; expro_name = ' '

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
       allocate(expro_z_eff(nexp))   ; expro_z_eff = 0.0

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
       allocate(expro_sbeame(nexp))     ; expro_sbeame = 0.0
       allocate(expro_sbcx(nexp))       ; expro_sbcx = 0.0
       allocate(expro_sscxl(nexp))      ; expro_sscxl = 0.0

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
       allocate(expro_drdrho(nexp))       ; expro_drdrho = 0.0
       allocate(expro_bp0(nexp))          ; expro_bp0 = 0.0
       allocate(expro_bt0(nexp))          ; expro_bt0 = 0.0
       allocate(expro_ip(nexp))           ; expro_ip = 0.0
       allocate(expro_gamma_e(nexp))      ; expro_gamma_e = 0.0
       allocate(expro_gamma_p(nexp))      ; expro_gamma_p = 0.0
       allocate(expro_mach(nexp))         ; expro_mach = 0.0
       allocate(expro_thetascale(nexp))   ; expro_thetascale = 0.0

       allocate(expro_dlnnidr(nion,nexp))  ; expro_dlnnidr = 0.0
       allocate(expro_dlntidr(nion,nexp))  ; expro_dlntidr = 0.0
       allocate(expro_sdlnnidr(nion,nexp)) ; expro_sdlnnidr = 0.0
       allocate(expro_sdlntidr(nion,nexp)) ; expro_sdlntidr = 0.0

    else

       deallocate(expro_mass) 
       deallocate(expro_z) 
       deallocate(expro_type) 
       deallocate(expro_name) 

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
       deallocate(expro_z_eff)

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
       deallocate(expro_sbeame)
       deallocate(expro_sbcx)
       deallocate(expro_sscxl)

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
       deallocate(expro_drdrho)       
       deallocate(expro_bp0)          
       deallocate(expro_bt0)          
       deallocate(expro_ip)           
       deallocate(expro_gamma_e)      
       deallocate(expro_gamma_p)      
       deallocate(expro_mach)   
       deallocate(expro_thetascale)  

       deallocate(expro_dlnnidr)
       deallocate(expro_dlntidr)  
       deallocate(expro_sdlnnidr) 
       deallocate(expro_sdlntidr)  

    endif

  end subroutine vpro_init

  subroutine vpro_read(path)

    implicit none

    character(len=*), intent(in) :: path 
    integer :: nexp,nion,ierr,i
    character*22 :: ytag,c

    ! ORDERING NOTE: nexp should appear before any profile arrays

    open(unit=1,file=trim(path)//infile,status='old')

    do 

       read(1,'(a)',end=99) ytag ; c = trim(ytag(3:22))

       select case (c)
       case ('nexp')
          call vpro_icomm(expro_n_exp) 
          nexp = expro_n_exp
       case ('nion')
          call vpro_icomm(expro_n_ion)
          nion = expro_n_ion
          if (allocated(expro_rho)) call vpro_init(0)
          call vpro_init(1) 
       case ('mass')
          call vpro_scomm(expro_mass,nion)
       case ('z')
          call vpro_scomm(expro_z,nion)
       case ('bt_exp')
          call vpro_rcomm(expro_b_ref) 
       case ('arho_exp')
          call vpro_rcomm(expro_arho) 
       case ('rho')
          call vpro_vcomm(expro_rho,nexp)  
       case ('rmin')
          call vpro_vcomm(expro_rmin,nexp)  
       case ('polflux')
          call vpro_vcomm(expro_polflux,nexp)  
       case ('q')
          call vpro_vcomm(expro_q,nexp)  
       case ('w0')
          call vpro_vcomm(expro_w0,nexp)  
       case ('rmaj')
          call vpro_vcomm(expro_rmaj,nexp)  
       case ('zmag')
          call vpro_vcomm(expro_zmag,nexp)  
       case ('kappa')
          call vpro_vcomm(expro_kappa,nexp)  
       case ('delta')
          call vpro_vcomm(expro_delta,nexp)  
       case ('zeta')
          call vpro_vcomm(expro_zeta,nexp)  
       case ('ne')
          call vpro_vcomm(expro_ne,nexp) 
       case ('Te')
          call vpro_vcomm(expro_te,nexp)  
       case ('ptot')
          call vpro_vcomm(expro_ptot,nexp)  
       case ('z_eff')
          call vpro_vcomm(expro_z_eff,nexp) 
       case ('ni')
          call vpro_acomm(expro_ni(:,:),nion,nexp) 
       case ('ti')
          call vpro_acomm(expro_ti(:,:),nion,nexp) 
       case ('vpol')
          call vpro_acomm(expro_vpol(:,:),nion,nexp) 
       case ('vtor')
          call vpro_acomm(expro_vtor(:,:),nion,nexp) 
       case ('flow_beam')
          call vpro_vcomm(expro_flow_beam,nexp)  
       case ('flow_wall')
          call vpro_vcomm(expro_flow_wall,nexp)  
       case ('flow_mom')
          call vpro_vcomm(expro_flow_mom,nexp)  
       case ('pow_e')
          call vpro_vcomm(expro_pow_e,nexp)  
       case ('pow_i')
          call vpro_vcomm(expro_pow_i,nexp)  
       case ('pow_ei')
          call vpro_vcomm(expro_pow_ei,nexp)
       case ('pow_e_aux')
          call vpro_vcomm(expro_pow_e_aux,nexp)  
       case ('pow_i_aux')
          call vpro_vcomm(expro_pow_i_aux,nexp)  
       case ('pow_e_fus')
          call vpro_vcomm(expro_pow_e_fus,nexp)  
       case ('pow_i_fus')
          call vpro_vcomm(expro_pow_i_fus,nexp)  
       case ('pow_e_sync')
          call vpro_vcomm(expro_pow_e_sync,nexp) 
       case ('pow_e_brem')
          call vpro_vcomm(expro_pow_e_brem,nexp) 
       case ('pow_e_line')
          call vpro_vcomm(expro_pow_e_line,nexp)  
       case ('sbeame')
          call vpro_vcomm(expro_sbeame,nexp) 
       case ('sbcx')
          call vpro_vcomm(expro_sbcx,nexp) 
       case ('sscxl')
          call vpro_vcomm(expro_sscxl,nexp) 
       end select

    enddo

99  close(1)

    ! ** input.profiles.geo **

    nexp = expro_n_exp
    open(unit=1,file=trim(path)//'input.profiles.geo',status='old',iostat=ierr)
    if (ierr == 0) then
       call vpro_skip_header(1)
       call vpro_icomm(expro_nfourier)
       if (allocated(expro_geo)) deallocate(expro_geo)
       if (allocated(expro_dgeo)) deallocate(expro_dgeo)
       allocate(expro_geo(4,0:expro_nfourier,nexp)) ; expro_geo(:,:,:)=0.0
       allocate(expro_dgeo(4,0:expro_nfourier,nexp)) ; expro_dgeo(:,:,:)=0.0
       do i=1,nexp
          call vpro_scomm(expro_geo(:,:,i),4*(expro_nfourier+1))
       enddo
    else
       expro_nfourier = -1
    endif
    close(1)

    ! BCAST HERE

    call vpro_compute_derived

30  format(1pe14.7)

  end subroutine vpro_read

  subroutine vpro_write

    implicit none

    integer :: nexp,nion

    nexp = expro_n_exp
    nion = expro_n_ion

    open(unit=1,file=infile,position='append')

    write(1,20) ident//tag(1)  ; write(1,'(i0)') nexp
    write(1,20) ident//tag(2)  ; write(1,'(i0)') nion
    write(1,20) ident//tag(3)  ; write(1,40) expro_mass
    write(1,20) ident//tag(4)  ; write(1,40) expro_z
    write(1,20) ident//tag(5)  ; write(1,30) expro_b_ref
    write(1,20) ident//tag(6)  ; write(1,30) expro_arho
    write(1,20) ident//tag(7)  ; call vpro_writev(expro_rho,nexp)
    write(1,20) ident//tag(8)  ; call vpro_writev(expro_rmin,nexp)
    write(1,20) ident//tag(9)  ; call vpro_writev(expro_polflux,nexp)
    write(1,20) ident//tag(10) ; call vpro_writev(expro_q,nexp)
    write(1,20) ident//tag(11) ; call vpro_writev(expro_w0,nexp)
    write(1,20) ident//tag(12) ; call vpro_writev(expro_rmaj,nexp)
    write(1,20) ident//tag(13) ; call vpro_writev(expro_zmag,nexp)
    write(1,20) ident//tag(14) ; call vpro_writev(expro_kappa,nexp)
    write(1,20) ident//tag(15) ; call vpro_writev(expro_delta,nexp)
    write(1,20) ident//tag(16) ; call vpro_writev(expro_zeta,nexp)
    write(1,20) ident//tag(17) ; call vpro_writev(expro_ne,nexp)
    write(1,20) ident//tag(18) ; call vpro_writev(expro_te,nexp)
    write(1,20) ident//tag(19) ; call vpro_writev(expro_ptot,nexp)
    write(1,20) ident//tag(20) ; call vpro_writev(expro_z_eff,nexp)
    write(1,20) ident//tag(21) ; call vpro_writev(expro_flow_beam,nexp)
    write(1,20) ident//tag(22) ; call vpro_writev(expro_flow_wall,nexp)
    write(1,20) ident//tag(23) ; call vpro_writev(expro_flow_mom,nexp)
    write(1,20) ident//tag(24) ; call vpro_writev(expro_pow_e,nexp)
    write(1,20) ident//tag(25) ; call vpro_writev(expro_pow_i,nexp)
    write(1,20) ident//tag(26) ; call vpro_writev(expro_pow_ei,nexp)
    write(1,20) ident//tag(27) ; call vpro_writev(expro_pow_e_aux,nexp)
    write(1,20) ident//tag(28) ; call vpro_writev(expro_pow_i_aux,nexp)
    write(1,20) ident//tag(29) ; call vpro_writev(expro_pow_e_fus,nexp)
    write(1,20) ident//tag(30) ; call vpro_writev(expro_pow_i_fus,nexp)
    write(1,20) ident//tag(31) ; call vpro_writev(expro_pow_e_sync,nexp)
    write(1,20) ident//tag(32) ; call vpro_writev(expro_pow_e_brem,nexp)
    write(1,20) ident//tag(33) ; call vpro_writev(expro_pow_e_line,nexp)
    write(1,20) ident//tag(34) ; call vpro_writev(expro_sbeame,nexp)
    write(1,20) ident//tag(35) ; call vpro_writev(expro_sbcx,nexp)
    write(1,20) ident//tag(36) ; call vpro_writev(expro_sscxl,nexp)
    write(1,20) ident//tag(37) ; call vpro_writea(expro_ni(:,:),nion,nexp)
    write(1,20) ident//tag(38) ; call vpro_writea(expro_ti(:,:),nion,nexp)
    write(1,20) ident//tag(39) ; call vpro_writea(expro_vpol(:,:),nion,nexp)
    write(1,20) ident//tag(40) ; call vpro_writea(expro_vtor(:,:),nion,nexp)

    close(1)

20  format(a)
30  format(1pe14.7)
40  format(10(1pe14.7))

  end subroutine vpro_write

  subroutine vpro_read_legacy

    implicit none

    integer :: i
    integer :: nexp,nion
    character*99 :: line
    double precision :: x(5)

    open(unit=1,file='input.profiles',status='old')
    do while (line(1:2) /= '#r')
       read(1,'(a)') line
       if (line(1:5) == 'N_EXP') then
          read(line(7:),*) expro_n_exp
       endif
       if (line(1:5) == 'N_ION') then
          read(line(7:),*) expro_n_ion
       endif
       if (line(1:6) == 'BT_EXP') then
          read(line(8:),*) expro_b_ref
       endif
       if (line(1:8) == 'ARHO_EXP') then
          read(line(10:),*) expro_arho
       endif
    enddo

    call vpro_init(1)

    nexp = expro_n_exp
    nion = expro_n_ion

    ! 1
    do i=1,nexp
       read(1,*) x
       expro_rho(i)     = x(1)
       expro_rmin(i)    = x(2)
       expro_polflux(i) = x(3)
       expro_q(i)       = x(4)
       expro_w0(i)      = x(5)
    enddo

    read(1,'(a)') line
    read(1,'(a)') line

    ! 2
    do i=1,nexp
       read(1,*) x
       expro_rmaj(i)  = x(1)
       expro_zmag(i)  = x(2)
       expro_kappa(i) = x(3)
       expro_delta(i) = x(4)
       expro_zeta(i)  = x(5)
    enddo

    read(1,'(a)') line
    read(1,'(a)') line

    ! 3
    do i=1,nexp
       read(1,*) x
       expro_ne(i)    = x(1)
       expro_te(i)    = x(2)
       expro_ptot(i)  = x(3)
       expro_z_eff(i) = x(4)
    enddo

    read(1,'(a)') line
    read(1,'(a)') line

    ! 4
    do i=1,nexp
       read(1,*) x
       expro_ni(1:nion,i) = x(1:nion)
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
       expro_ti(1:nion,i) = x(1:nion)
    enddo

    read(1,'(a)') line
    read(1,'(a)') line

    ! 7 (assume < 6 ions, so skip)
    do i=1,nexp
       read(1,*) x
    enddo

    read(1,'(a)') line
    read(1,'(a)') line

    ! 8
    do i=1,nexp
       read(1,*) x
       expro_vtor(1:nion,i) = x(1:nion)
    enddo

    read(1,'(a)') line
    read(1,'(a)') line

    ! 9 (assume < 6 ions, so skip)
    do i=1,nexp
       read(1,*) x
    enddo

    read(1,'(a)') line
    read(1,'(a)') line

    ! 10
    do i=1,nexp
       read(1,*) x
       expro_vpol(1:nion,i) = x(1:nion)
    enddo

    read(1,'(a)') line
    read(1,'(a)') line

    ! 11 (assume < 6 ions, so skip)
    do i=1,nexp
       read(1,*) x
    enddo
    close(1)

  end subroutine vpro_read_legacy

end module vpro
