module expro

  ! Fundamental input
  integer, parameter :: nextag = 42
  character(len=12) :: infile = 'input.gacode'

  character(len=2) :: ident='# '
  character(len=20), dimension(nextag) :: extag = (/&
       'nexp                ',& !1
       'nion                ',& !2
       'mass (m_H)          ',& !3
       'z (-)               ',& !4
       'torfluxa (Wb/radian)',& !5
       'rho (-)             ',& !6
       'rmin (m)            ',& !7
       'polflux (Wb/radian) ',& !8
       'q (-)               ',& !9
       'w0 (rad/s)          ',& !10
       'rmaj (m)            ',& !11
       'zmag (m)            ',& !12
       'kappa (-)           ',& !13
       'delta (-)           ',& !14
       'zeta (-)            ',& !15
       'ne (10^19/m^3)      ',& !16
       'ni (10^19/m^3)      ',& !17
       'Te (keV)            ',& !18
       'Ti (keV)            ',& !19
       'ptot (Pa)           ',& !20
       'jbs (MA/m^2)        ',& !21
       'jbstor (MA/m^2)     ',& !22
       'sigmapar (MS/m)     ',& !23
       'z_eff (-)           ',& !24
       'vpol (m/s)          ',& !25
       'vtor (m/s)          ',& !26
       'flow_beam (MW/keV)  ',& !27
       'flow_wall (MW/keV)  ',& !28
       'flow_mom (Nm)       ',& !29
       'pow_e (MW)          ',& !30
       'pow_i (MW)          ',& !31
       'pow_ei (MW)         ',& !32
       'pow_e_aux (MW)      ',& !33
       'pow_i_aux (MW)      ',& !34
       'pow_e_fus (MW)      ',& !35
       'pow_i_fus (MW)      ',& !36
       'pow_e_sync (MW)     ',& !37
       'pow_e_brem (MW)     ',& !38
       'pow_e_line (MW)     ',& !39
       'sbeame (1/m^3/s)    ',& !40
       'sbcx (1/m^3/s)      ',& !41
       'sscxl (1/m^3/s)     '&  !42
       /)

  integer :: &
       expro_n_exp,&
       expro_n_ion

  double precision, dimension(:), allocatable :: &
       expro_mass,&
       expro_z

  character(len=7), dimension(:), allocatable :: &
       expro_name,&
       expro_type

  double precision :: &
       expro_torfluxa,&
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
       expro_jbs,&
       expro_jbstor,&
       expro_sigmapar,&
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

  ! header information
  character(len=50), dimension(6) :: expro_header = (/&
       '#  original :                                     ',& !1
       '# statefile :                                     ',& !2
       '#     gfile :                                     ',& !3
       '#   cerfile :                                     ',& !4
       '#      vgen :                                     ',& !5
       '#     tgyro :                                     '& !6
       /)

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
       allocate(expro_type(nion))    ; expro_type = '       '
       allocate(expro_name(nion))    ; expro_name = '       '

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
       allocate(expro_jbs(nexp))     ; expro_jbs = 0.0
       allocate(expro_jbstor(nexp))  ; expro_jbstor = 0.0
       allocate(expro_sigmapar(nexp)); expro_sigmapar = 0.0
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
       deallocate(expro_jbs)
       deallocate(expro_jbstor)
       deallocate(expro_sigmapar)
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

  end subroutine expro_init

  subroutine expro_read(path)

    implicit none

    character(len=*), intent(in) :: path 
    integer :: nexp,nion,ierr,i,nd
    character(len=22) :: ytag,c

    ! ORDERING NOTE: nexp should appear before any profile arrays

    open(unit=1,file=trim(path)//infile,status='old')

    ! header
    do i=1,6
       read(1,'(a)') expro_header(i)
    enddo

    do 

       read(1,'(a)',end=99) ytag

       nd = scan(ytag,'(')
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
       case ('mass')
          call expro_lcomm(expro_mass,nion)
       case ('z')
          call expro_lcomm(expro_z,nion)
       case ('torfluxa')
          call expro_rcomm(expro_torfluxa) 
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
       case ('Te')
          call expro_vcomm(expro_te,nexp)  
       case ('ptot')
          call expro_vcomm(expro_ptot,nexp)  
       case ('jbs ')
          call expro_vcomm(expro_jbs,nexp)  
       case ('jbstor')
          call expro_vcomm(expro_jbstor,nexp)  
       case ('sigmapar')
          call expro_vcomm(expro_sigmapar,nexp)  
       case ('z_eff')
          call expro_vcomm(expro_z_eff,nexp) 
       case ('ni')
          call expro_acomm(expro_ni(:,:),nion,nexp) 
       case ('Ti')
          call expro_acomm(expro_ti(:,:),nion,nexp) 
       case ('vpol')
          call expro_acomm(expro_vpol(:,:),nion,nexp) 
       case ('vtor')
          call expro_acomm(expro_vtor(:,:),nion,nexp) 
       case ('flow_beam')
          call expro_vcomm(expro_flow_beam,nexp)  
       case ('flow_wall')
          call expro_vcomm(expro_flow_wall,nexp)  
       case ('flow_mom')
          call expro_vcomm(expro_flow_mom,nexp)  
       case ('pow_e')
          call expro_vcomm(expro_pow_e,nexp)  
       case ('pow_i')
          call expro_vcomm(expro_pow_i,nexp)  
       case ('pow_ei')
          call expro_vcomm(expro_pow_ei,nexp)
       case ('pow_e_aux')
          call expro_vcomm(expro_pow_e_aux,nexp)  
       case ('pow_i_aux')
          call expro_vcomm(expro_pow_i_aux,nexp)  
       case ('pow_e_fus')
          call expro_vcomm(expro_pow_e_fus,nexp)  
       case ('pow_i_fus')
          call expro_vcomm(expro_pow_i_fus,nexp)  
       case ('pow_e_sync')
          call expro_vcomm(expro_pow_e_sync,nexp) 
       case ('pow_e_brem')
          call expro_vcomm(expro_pow_e_brem,nexp) 
       case ('pow_e_line')
          call expro_vcomm(expro_pow_e_line,nexp)  
       case ('sbeame')
          call expro_vcomm(expro_sbeame,nexp) 
       case ('sbcx')
          call expro_vcomm(expro_sbcx,nexp) 
       case ('sscxl')
          call expro_vcomm(expro_sscxl,nexp) 
       end select

    enddo

99  close(1)

    ! ** input.gacode.geo **

    nexp = expro_n_exp
    open(unit=1,file=trim(path)//infile//'.geo',status='old',iostat=ierr)
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

    call expro_compute_derived

  end subroutine expro_read

  subroutine expro_write(thisinfile)

    implicit none

    integer :: i,nexp,nion
    character(len=*) :: thisinfile 

    nexp = expro_n_exp
    nion = expro_n_ion

    ! Write header
    open(unit=1,file=trim(thisinfile),status='replace')
    do i=1,6
       write(1,'(a)') expro_header(i)
    enddo
    write(1,'(a)') '#'

    ! Write data
    write(1,20) ident//extag(1)  ; write(1,'(i0)') nexp
    write(1,20) ident//extag(2)  ; write(1,'(i0)') nion
    write(1,20) ident//extag(3)  ; write(1,40) expro_mass
    write(1,20) ident//extag(4)  ; write(1,40) expro_z
    write(1,20) ident//extag(5)  ; write(1,30) expro_torfluxa
    write(1,20) ident//extag(6)  ; call expro_writev(expro_rho,nexp)
    write(1,20) ident//extag(7)  ; call expro_writev(expro_rmin,nexp)
    write(1,20) ident//extag(8)  ; call expro_writev(expro_polflux,nexp)
    write(1,20) ident//extag(9)  ; call expro_writev(expro_q,nexp)
    write(1,20) ident//extag(10) ; call expro_writev(expro_w0,nexp)
    write(1,20) ident//extag(11) ; call expro_writev(expro_rmaj,nexp)
    write(1,20) ident//extag(12) ; call expro_writev(expro_zmag,nexp)
    write(1,20) ident//extag(13) ; call expro_writev(expro_kappa,nexp)
    write(1,20) ident//extag(14) ; call expro_writev(expro_delta,nexp)
    write(1,20) ident//extag(15) ; call expro_writev(expro_zeta,nexp)
    write(1,20) ident//extag(16) ; call expro_writev(expro_ne,nexp)
    write(1,20) ident//extag(17) ; call expro_writea(expro_ni(:,:),nion,nexp)
    write(1,20) ident//extag(18) ; call expro_writev(expro_te,nexp)
    write(1,20) ident//extag(19) ; call expro_writea(expro_ti(:,:),nion,nexp)
    write(1,20) ident//extag(20) ; call expro_writev(expro_ptot,nexp)
    write(1,20) ident//extag(21) ; call expro_writev(expro_jbs,nexp)
    write(1,20) ident//extag(22) ; call expro_writev(expro_jbstor,nexp)
    write(1,20) ident//extag(23) ; call expro_writev(expro_sigmapar,nexp)
    write(1,20) ident//extag(24) ; call expro_writev(expro_z_eff,nexp)
    write(1,20) ident//extag(25) ; call expro_writea(expro_vpol(:,:),nion,nexp)
    write(1,20) ident//extag(26) ; call expro_writea(expro_vtor(:,:),nion,nexp)
    write(1,20) ident//extag(27) ; call expro_writev(expro_flow_beam,nexp)
    write(1,20) ident//extag(28) ; call expro_writev(expro_flow_wall,nexp)
    write(1,20) ident//extag(29) ; call expro_writev(expro_flow_mom,nexp)
    write(1,20) ident//extag(30) ; call expro_writev(expro_pow_e,nexp)
    write(1,20) ident//extag(31) ; call expro_writev(expro_pow_i,nexp)
    write(1,20) ident//extag(32) ; call expro_writev(expro_pow_ei,nexp)
    write(1,20) ident//extag(33) ; call expro_writev(expro_pow_e_aux,nexp)
    write(1,20) ident//extag(34) ; call expro_writev(expro_pow_i_aux,nexp)
    write(1,20) ident//extag(35) ; call expro_writev(expro_pow_e_fus,nexp)
    write(1,20) ident//extag(36) ; call expro_writev(expro_pow_i_fus,nexp)
    write(1,20) ident//extag(37) ; call expro_writev(expro_pow_e_sync,nexp)
    write(1,20) ident//extag(38) ; call expro_writev(expro_pow_e_brem,nexp)
    write(1,20) ident//extag(39) ; call expro_writev(expro_pow_e_line,nexp)
    write(1,20) ident//extag(40) ; call expro_writev(expro_sbeame,nexp)
    write(1,20) ident//extag(41) ; call expro_writev(expro_sbcx,nexp)
    write(1,20) ident//extag(42) ; call expro_writev(expro_sscxl,nexp)

    close(1)

20  format(a)
30  format(1pe14.7)
40  format(10(1pe14.7))

  end subroutine expro_write

  subroutine expro_read_legacy

    implicit none

    integer :: i
    integer :: nexp,nion
    character(len=99) :: line
    double precision :: x(5)
    double precision :: b_ref,arho

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
          read(line(8:),*) b_ref
       endif
       if (line(1:8) == 'ARHO_EXP') then
          read(line(10:),*) arho
       endif
    enddo
    expro_torfluxa = 0.5*b_ref*arho**2

    call expro_init(1)

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

  end subroutine expro_read_legacy

end module expro
