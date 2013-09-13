MODULE paleo
    IMPLICIT NONE
    REAL*8, ALLOCATABLE, DIMENSION(:) ::  &
    &  ne_p,    &   !electron density m^-3
    & te_p,     &       !electron temperature eV
    & zeff_p,   &       !Effective ion charge
    & rho_p,    &       !sqrt(Toroidal Flux) m
    & kappa_p,  &       !Elongation
    & hcap_p,   &       !Geometric Factor
    & q_p,      &       !Safety Factor
    & eta_nc_p, &       !Neoclassical parallel resistivity ohm-m
    & D_eta_p             !Paleoclassical diffusivity m^2/s
    REAL*8 :: &
    & R0_p              !Major Radius m
    REAL*8, PARAMETER :: &
    & pi=3.141592653589793238, &
    & mu_0=4*pi*1e-7,    &       !Permeability of free space [H/m]
    & m_e=9.1094e-31,   &       !Electron mass [kg]
    & m_D=2*1.6726e-27, &       !Deuteron mass [kg]
    & q_e=1.60218e-19,  &       !Electron charge [C]
    & eV2J=1.60218e-19          !Conversion from eV to Joules []

CONTAINS
    SUBROUTINE init_paleo_given_neoclassical(&
        & ne_in,    &   !electron density m^-3
        & te_in,    &   !electron temperature eV
        & zeff_in,  &   !Effective ion charge
        & rho_in,   &   !sqrt(Toroidal Flux) m
        & kappa_in, &   !Elongation
        & hcap_in,  &   !Geometric Factor
        & q_in,     &   !Safety Factor
        & R0_in,    &   !Major Radius m
        & eta_nc_in)    !Neoclassical parallel resistivity ohm-m
        
        REAL*8, DIMENSION(:), INTENT(IN) :: &
        & ne_in,    &   !electron density m^-3
        & te_in,    &   !electron temperature eV
        & zeff_in,  &   !Effective ion charge
        & rho_in,   &   !sqrt(Toroidal Flux) m
        & kappa_in, &   !Elongation
        & hcap_in,  &   !Geometric Factor
        & q_in,     &   !Safety Factor
        & eta_nc_in     !Neoclassical parallel resistivity ohm-m
        REAL*8, INTENT(IN) :: &
        & R0_in         !Major Radius m
        INTEGER :: &
        & n             !Size of arrays
        !Initialize module variables
        n=size(ne_in)
        IF (.NOT. ALLOCATED(ne_p)) THEN
            ALLOCATE(ne_p(n),te_p(n),zeff_p(n),rho_p(n),kappa_p(n),hcap_p(n),&
            &q_p(n),eta_nc_p(n),D_eta_p(n))
        ENDIF
        ne_p    =ne_in
        te_p    =te_in
        zeff_p  =zeff_in
        rho_p   =rho_in
        kappa_p =kappa_in
        hcap_p  =hcap_in
        q_p     =q_in
        eta_nc_p=eta_nc_in
        R0_p    =R0_in
        D_eta_p=eta_nc_p/mu_0
        
    ENDSUBROUTINE init_paleo_given_neoclassical
    FUNCTION dydx( x, y )
        REAL*8, INTENT(IN), DIMENSION(0:) :: x,y
        REAL*8, DIMENSION(0:size(x)-1) :: dydx
        REAL*8, DIMENSION(0:size(x)-3) :: dy1, dy2, dx1, dx2
        INTEGER :: n
        n=size(x)
        
        dy1=y(2:)-y(1:n-2)
        dy2=y(1:n-2)-y(0:n-3)
        dx1=x(2:)-x(1:n-2)
        dx2=x(1:n-2)-x(0:n-3)
        dydx(1:n-2)=(dy1/dx1+dy2/dx2)/2.
        dydx(0)=(-(y(2)-y(0))*dx2(0)**2+(x(2)-x(0))**2*dy2(0))/((x(2)-x(0))*dx1(0)*dx2(0))
        dydx(n-1)=(-(y(n-3)-y(n-1))*(x(n-2)-x(n-1))**2+(y(n-2)-y(n-1))*(x(n-3)-x(n-1))**2)
        dydx(n-1)=dydx(n-1)/((x(n-3)-x(n-1))*(x(n-2)-x(n-1))*(x(n-3)-x(n-2)))
    ENDFUNCTION dydx
    FUNCTION D_eta( eta_nc_in )
        REAL*8, DIMENSION(:), INTENT(IN) :: eta_nc_in ![ohm-m]
        REAL*8, DIMENSION(size(eta_nc_in)):: D_eta ![m^2/s]
        D_eta =  eta_nc_in/mu_0   !result
    ENDFUNCTION D_eta
    FUNCTION a_bar( rho, kappa )
        REAL*8, DIMENSION(:), INTENT(IN) :: rho, kappa ![m], []
        REAL*8, DIMENSION(size(rho)) :: a_bar ![m]
        a_bar=maxval(rho)*(2.*kappa**2/(1.+kappa**2))**0.5 !UW-CPTC 07-5 Eq. (11)
    ENDFUNCTION a_bar
    FUNCTION D_eta_bar( eta_nc_in, rho, kappa )
        REAL*8, DIMENSION(:), INTENT(IN) :: eta_nc_in, rho, kappa   ![ohm-m], [m], []
        REAL*8, DIMENSION(size(eta_nc_in)) :: D_eta_bar ![m^2/s]
        D_eta_bar = D_eta(eta_nc_in)*maxval(rho)**2/a_bar(rho,kappa)**2    !result
    ENDFUNCTION D_eta_bar 
    FUNCTION R_bar( R0 )
        REAL*8 :: R_bar ![m] 
        REAL*8, INTENT(IN) :: R0    ![m]
        R_bar = R0
    ENDFUNCTION
    FUNCTION delta_e( ne )
        REAL*8, INTENT(IN), DIMENSION(:) :: ne
        REAL*8, DIMENSION(size(ne)) :: delta_e ![m]
        delta_e = (m_e/(mu_0*q_e**2*ne))**.5 !NRL Formulary
    ENDFUNCTION delta_e
    FUNCTION delta_e_bar( ne, rho, kappa )
        REAL*8, INTENT(IN), DIMENSION(:) :: ne, rho, kappa ![m^-3], [m], []
        REAL*8, DIMENSION(size(ne)) :: delta_e_bar ![]
        delta_e_bar =  delta_e(ne)/a_bar(rho,kappa) !Eq (14)
    ENDFUNCTION delta_e_bar
    FUNCTION n_max( ne, q, rho, kappa )
        REAL*8, INTENT(IN), DIMENSION(:) :: ne, q, rho, kappa ![m^-3], [], [m], []
        REAL*8, DIMENSION(size(ne)) :: n_max, n_maxp ![], []
        n_max = (1./(pi*delta_e_bar(ne,rho,kappa)*abs(dydx(rho,q))))**0.5 !Eq (13)
        n_maxp  = 1. / ((pi*delta_e_bar(ne,rho,kappa))**2*&
            & ABS(dydx(rho,dydx(rho,q))))**0.3333 !UW-CPTC 07-5 Eq (16)
        WHERE (n_max>n_maxp) n_max=n_maxp
    ENDFUNCTION n_max
    FUNCTION l_max( ne, q, rho, kappa, R0 )
        REAL*8, DIMENSION(:), INTENT(IN) :: ne, q, rho, kappa ![m^-3], [], [m], []
        REAL*8, DIMENSION(size(ne)) :: l_max ![]
        REAL*8 :: R0 ![m]
        l_max = pi*R_bar(R0)*q*n_max(ne,q,rho,kappa) !Eq (12)
    ENDFUNCTION l_max
    FUNCTION ln_Lambda( ne, te )
        REAL*8, INTENT(IN), DIMENSION(:) :: ne, te ![m^-3], [eV]
        REAL*8, DIMENSION(size(ne)) :: ln_Lambda ![]
        ln_Lambda = 24.0-log((ne/1e6)**.5/te) !NRL Formulary 2002 pg 34 (b)-2
    ENDFUNCTION ln_Lambda
    FUNCTION lambda_e( ne, te, zeff )
        REAL*8, INTENT(IN), DIMENSION(:) :: ne, te, zeff ![m^-3], [eV], []
        REAL*8, DIMENSION(size(ne)) :: lambda_e ![]
        lambda_e =     1.2e16*te**2/(zeff*ne)*(17./ln_Lambda(ne,te)) !Eq (8)
    ENDFUNCTION lambda_e
    FUNCTION M( ne, te, zeff, q, rho, kappa, R0 )
        REAL*8, DIMENSION(:), INTENT(IN) :: ne,te,zeff,q,rho,kappa ! [m^-3], [eV], [], [], [m], []
        REAL*8, INTENT(IN) :: R0 ![m]
        REAL*8, DIMENSION(size(ne)) :: M, l_max_tmp, lambda_e_tmp ![]
        l_max_tmp=l_max(ne,q,rho,kappa,R0)
        lambda_e_tmp=lambda_e(ne,te,zeff)
        M=0
        WHERE(l_max_tmp+lambda_e_tmp>0)
            M = l_max_tmp*lambda_e_tmp/(pi*R0*q)/(l_max_tmp+lambda_e_tmp) !Eq (11)
        ENDWHERE
    ENDFUNCTION M
    FUNCTION V_Gamma_pc( ne, rho, kappa, hcap, eta_nc )
        REAL*8, INTENT(IN), DIMENSION(:) :: ne, rho, kappa, hcap, eta_nc ![m^-3], [m], [], [], [ohm-m]
        REAL*8, DIMENSION(size(ne)) :: V_Gamma_pc ![m/s]
        !*******(-1) Factor? *******
        V_Gamma_pc = dydx(rho,hcap*rho*D_eta_bar(eta_nc,rho,kappa))/(hcap*rho)
        !***************************
    ENDFUNCTION V_Gamma_pc
    FUNCTION Gamma_pc( ne, rho, kappa, hcap, eta_nc )
        REAL*8, INTENT(IN), DIMENSION(:) :: ne, rho, kappa, hcap, eta_nc
        REAL*8, DIMENSION(size(ne)) :: Gamma_pc
        Gamma_pc = -1./(hcap*rho)*dydx(rho,hcap*rho*D_eta_bar(eta_nc,rho,kappa)*ne)
    ENDFUNCTION Gamma_pc
    FUNCTION fsa_div_gamma_pc( ne, rho, kappa, hcap, eta_nc )
        REAL*8, INTENT(IN), DIMENSION(:) :: ne, rho, kappa, hcap, eta_nc
        REAL*8, DIMENSION(size(ne)) :: fsa_div_gamma_pc
        INTEGER :: i, iout
!        write(*,'(a)') 'Inside div_gamma_pc' 
!        write(*,'(a,i2)') 'len(ne)=',size(ne),'len(rho)=',size(rho),&
!        & 'len(kappa)=',size(kappa),'len(hcap)',size(hcap),'len(eta_nc)',size(eta_nc)
        fsa_div_gamma_pc = 1/(hcap*rho)*dydx(rho,hcap*rho*&
            & Gamma_pc(ne,rho,kappa,hcap,eta_nc))
!        OPEN(unit=iout,file='fsa_div_gamma_pc.dat',status='unknown')
!    write(iout,'(6a20)') '#ne','r','kappa','hcap','eta',&
!        &  'fsa_div_gamma_pc'
!    write(iout,'(6a20)') '#1/m^3','m','1','1','ohm-m',&
!        & '1/m^3/s'
!    DO i=1,size(ne)
!        write(iout,'(6g20.8)') ne(i),rho(i),kappa(i),hcap(i),eta_nc(i),&
!        &  fsa_div_gamma_pc(i)
!    ENDDO
!    close(iout)
    ENDFUNCTION fsa_div_gamma_pc
    FUNCTION fsa_div_qe_pc( ne, te, zeff, q, rho, kappa, R0, hcap, eta_nc )
        REAL*8, INTENT(IN), DIMENSION(:) :: ne, te, zeff, q, rho, kappa, hcap, eta_nc
        REAL*8, INTENT(IN) :: R0
        REAL*8, DIMENSION(size(ne)) :: fsa_div_qe_pc
        fsa_div_qe_pc = (M(ne,te,zeff,q,rho,kappa,R0)+1.)/(hcap*rho)*&
        & dydx(rho,dydx(rho,hcap*rho*3./2.*ne*te*eV2J*D_eta_bar(eta_nc,rho,kappa)))
    ENDFUNCTION fsa_div_qe_pc
    FUNCTION fsa_div_qi_pc( ni, ti, rho, kappa, hcap, eta_nc )
        REAL*8, INTENT(IN), DIMENSION(:) :: ni, ti, rho, kappa, hcap, eta_nc
        REAL*8, DIMENSION(size(ni)) :: fsa_div_qi_pc
        fsa_div_qi_pc = (1.)/(hcap*rho)*&
        & dydx(rho,dydx(rho,hcap*rho*3./2.*ni*ti*eV2J*D_eta_bar(eta_nc,rho,kappa)))
    ENDFUNCTION fsa_div_qi_pc
    FUNCTION D_qe_pc( ne, te ,zeff, q, rho, kappa, R0,  eta_nc )
        REAL*8, INTENT(IN), DIMENSION(:) :: ne, te ,zeff, q, rho, kappa, eta_nc
        REAL*8, INTENT(IN) :: R0
        REAL*8, DIMENSION(size(ne)) :: D_qe_pc
        D_qe_pc = (M(ne,te,zeff,q,rho,kappa,R0)+1.)*3./2.*ne*D_eta_bar(eta_nc,rho,kappa)
    ENDFUNCTION D_qe_pc
    FUNCTION V_qe_pc( ne, te, zeff, q, rho, kappa, R0, hcap, eta_nc )
        REAL*8, INTENT(IN), DIMENSION(:) :: ne, te, zeff, q, rho, kappa, hcap, eta_nc
        REAL*8, INTENT(IN) :: R0
        REAL*8, DIMENSION(size(ne)) :: V_qe_pc
        V_qe_pc = (M(ne,te,zeff,q,rho,kappa,R0)+1.)/(hcap*rho)*&
        & dydx(rho,hcap*rho*3./2.*ne*D_eta_bar(eta_nc,rho,kappa))
    ENDFUNCTION V_qe_pc
    FUNCTION W_qe_pc( ne, te, zeff, q, rho, kappa, R0, hcap, eta_nc )
        REAL*8, INTENT(IN), DIMENSION(:) :: ne, te, zeff, q, rho, kappa, hcap, eta_nc
        REAL*8, INTENT(IN) :: R0
        REAL*8, DIMENSION(size(ne)) :: W_qe_pc
        W_qe_pc = 1./(hcap*rho)*dydx(rho,(M(ne,te,zeff,q,rho,kappa,R0)+1.))*&
        & dydx(rho,hcap*rho*3./2.*ne*D_eta_bar(eta_nc,rho,kappa))
    ENDFUNCTION W_qe_pc
    FUNCTION S_qe_pc( ne, te, zeff, q, rho, kappa, R0, eta_nc )
        REAL*8, INTENT(IN), DIMENSION(:) :: ne, te, zeff, q, rho, kappa, eta_nc
        REAL*8, INTENT(IN) :: R0
        REAL*8, DIMENSION(size(ne)) :: S_qe_pc
        S_qe_pc = 3./2.*ne*D_eta_bar(eta_nc,rho,kappa)*&
        & dydx(rho,M(ne,te,zeff,q,rho,kappa,R0))*dydx(rho,te)*eV2J
    ENDFUNCTION S_qe_pc
    FUNCTION D_qi_pc( ni, rho, kappa, eta_nc )
        REAL*8, INTENT(IN), DIMENSION(:) :: ni, rho, kappa, eta_nc
        REAL*8, DIMENSION(size(ni)) :: D_qi_pc
        D_qi_pc = 3./2.*ni*D_eta_bar(eta_nc,rho,kappa)
    ENDFUNCTION D_qi_pc
    FUNCTION V_qi_pc( ni, rho, kappa, hcap, eta_nc )
        REAL*8, INTENT(IN), DIMENSION(:) :: ni, rho, kappa, hcap, eta_nc
        REAL*8, DIMENSION(size(ni)) :: V_qi_pc
        V_qi_pc = 1/(hcap*rho)*&
        & dydx(rho,hcap*rho*3./2.*ni*D_eta_bar(eta_nc,rho,kappa))
    ENDFUNCTION V_qi_pc
ENDMODULE paleo

