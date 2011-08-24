module gkcoll_equilibrium
  
  implicit none
  
  public :: EQUIL_alloc, EQUIL_do
  
  ! equilibrium parameters (theta)
  real :: d_theta                                 ! delta theta
  ! local (th)
  real, dimension(:), allocatable :: k_par        ! bhat dot grad/a
  real, dimension(:), allocatable :: v_drift_x    ! radial curvature drift vel
  real, dimension(:), allocatable :: v_drift_th   ! theta  curvature drift vel
  real, dimension(:), allocatable :: gradr        ! | grad r|
  real, dimension(:), allocatable :: w_theta    ! flux surface avg weights
  real, dimension(:), allocatable :: Btor       ! b_t / Bunit
  real, dimension(:), allocatable :: Bpol       ! b_p / Bunit
  real, dimension(:), allocatable :: Bmag       ! B/Bunit
  real, dimension(:), allocatable :: gradpar_Bmag ! bhat dot grad/aBmag/Bunit
  real, dimension(:), allocatable :: bigR       ! R/a (global)
  real, dimension(:,:), allocatable :: k_perp
  ! radial deriv of the volume enclosed by a flux surface
  real, dimension(:), allocatable :: v_prime_g  ! (ngr)
  
  ! Parameters needed for post-processing analysis
  real :: I_div_psip    ! I(psi)/psi_prime = q f / r
  
  logical, private :: initialized = .false.
  
contains
  
  subroutine EQUIL_alloc(flag)
    use gkcoll_globals
    use GEO_interface
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    integer :: it
    integer, parameter :: geo_ntheta=1001 ! num grid pts for Miller geo grid
    
    if(flag == 1) then
       if(initialized) return
       
       allocate(theta(n_theta))
       allocate(k_par(n_theta))
       allocate(v_drift_x(n_theta))
       allocate(v_drift_th(n_theta))
       allocate(gradr(n_theta))
       allocate(w_theta(n_theta))
       allocate(bigR(n_theta))
       allocate(Btor(n_theta))
       allocate(Bpol(n_theta))
       allocate(Bmag(n_theta))
       allocate(gradpar_Bmag(n_theta))
       allocate(v_prime_g(n_gr))
       allocate(k_perp(n_theta,n_radial))
       

       d_theta = 2*pi/n_theta
       do it=1,n_theta
          theta(it) = -pi+(it-1)*d_theta
       enddo
       
       if(equilibrium_model == 0) then
          GEO_model_in = 0
       else if (equilibrium_model == 2 .or. equilibrium_model == 3) then
          GEO_model_in    = geo_numeq_flag
       endif
       GEO_ntheta_in   = geo_ntheta
       GEO_nfourier_in = geo_ny
       call GEO_alloc(1)

       initialized = .true.
       
    else
       if(.NOT. initialized) return
       
       deallocate(theta)
       deallocate(k_par)
       deallocate(v_drift_x)
       deallocate(v_drift_th)
       deallocate(gradr)
       deallocate(w_theta)
       deallocate(bigR)
       deallocate(Btor)
       deallocate(Bpol)
       deallocate(Bmag)
       deallocate(gradpar_Bmag)
       deallocate(v_prime_g)
       deallocate(k_perp)

       call GEO_alloc(0)
        
       initialized = .false.

    endif
    
  end subroutine EQUIL_alloc
  
  subroutine EQUIL_do(ir)
    use gkcoll_globals
    use GEO_interface
    implicit none
    integer, intent(in) :: ir
    integer :: it, jt, id
    real :: sum
    ! parameters needed for Miller equilibrium

    sum=0.0

    ! geo_numeq_flag, geo_ny, and geo_yin already set
    
    GEO_signb_in     = sign_bunit
    GEO_rmin_in      = r(ir)
    GEO_rmaj_in      = rmaj(ir)
    GEO_drmaj_in     = shift(ir)
    GEO_zmag_in      = zmag(ir)
    GEO_dzmag_in     = s_zmag(ir)
    GEO_q_in         = q(ir)
    GEO_s_in         = shat(ir)
    GEO_kappa_in     = kappa(ir)
    GEO_s_kappa_in   = s_kappa(ir)
    GEO_delta_in     = delta(ir)
    GEO_s_delta_in   = s_delta(ir)
    GEO_zeta_in      = zeta(ir)
    GEO_s_zeta_in    = s_zeta(ir)
    GEO_beta_star_in = 0.0    ! EAB: set beta_star=0;t in first-order calc
    ! NOTE: it is t implemented in v_drift defs
    GEO_fourier_in(:,:) = geo_yin(:,:,ir)
    call GEO_do()  
    
    do it=1,n_theta
       call GEO_interp(theta(it))
       k_par(it) = 1.0 / (q(ir) * rmaj(ir) * GEO_g_theta)
       bigR(it) = GEO_bigr
       Bmag(it)  = GEO_b
       Btor(it)  = GEO_bt
       Bpol(it)  = GEO_bp
       gradpar_Bmag(it) = k_par(it) * GEO_dbdt
       gradr(it)        = GEO_grad_r
       v_drift_x(it)  = -rho(ir)/(rmaj(ir) * Bmag(it)) * &
            GEO_grad_r * GEO_gsin
       v_drift_th(it) = -rho(ir)/(rmaj(ir) * Bmag(it) * GEO_l_t) &
            * (GEO_gcos1 * GEO_bt / Bmag(it)  & 
            -  GEO_gsin * GEO_nsin * GEO_grad_r) &
            * r(ir)

       ! flux-surface average weights
       w_theta(it) = GEO_g_theta / Bmag(it)
       sum = sum + w_theta(it)

       do ir=1,n_radial
          k_perp(it,ir) = sqrt((2.0*pi*indx_r(ir)*GEO_grad_r/r_length &
               + k_theta*GEO_gq*GEO_captheta)**2 &
               + (k_theta*GEO_gq)**2)
       enddo
       
    enddo
    
    I_div_psip = GEO_f * q(ir) / r(ir)
    
    do it=1,n_theta
       w_theta(it) = w_theta(it) / sum
    enddo
    
    ! v_prime_g = d Volume of flux surface / dr
    v_prime_g(ir) = r(ir)*rmaj(ir)*sum/n_theta * 4 * pi * pi
    
  end subroutine EQUIL_DO


end module gkcoll_equilibrium
