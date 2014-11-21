module cgyro_equilibrium
  
  implicit none
  
  public :: EQUIL_alloc, EQUIL_do
  
  ! equilibrium parameters (theta)
  real :: d_theta
  real, dimension(:,:), allocatable   :: theta_B
  real, dimension(:), allocatable   :: w_theta
  real, dimension(:,:), allocatable :: k_perp    
  real, dimension(:), allocatable   :: Bmag
  real, dimension(:,:), allocatable :: omega_stream
  real, dimension(:,:), allocatable :: omega_trap
  real, dimension(:,:), allocatable :: omega_rdrift
  real, dimension(:,:), allocatable :: omega_adrift
  real, dimension(:,:), allocatable :: omega_aprdrift
  
  logical, private :: initialized = .false.
  
contains
  
  subroutine EQUIL_alloc(flag)
    use cgyro_globals
    use GEO_interface
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    integer :: it, ir
    integer, parameter :: geo_ntheta=1001 ! num grid pts for Miller geo grid
    
    if(flag == 1) then
       if(initialized) return
       
       allocate(theta(n_theta))
       allocate(theta_B(n_radial/box_size,n_theta))
       allocate(w_theta(n_theta))
       allocate(Bmag(n_theta))
       allocate(k_perp(n_theta,n_radial))
       allocate(omega_stream(n_theta,n_species))
       allocate(omega_trap(n_theta,n_species))
       allocate(omega_rdrift(n_theta,n_species))
       allocate(omega_adrift(n_theta,n_species))
       allocate(omega_aprdrift(n_theta,n_species))

       d_theta = (2*pi/n_theta)
       do it=1,n_theta
          theta(it) = -pi+(it-1)*d_theta
       enddo

       do ir=1,n_radial/box_size
          do it=1,n_theta
             theta_B(ir,it) = theta(it)+2*pi*(ir-1-n_radial/2/box_size)
          enddo
       enddo
       
       if(equilibrium_model == 0) then
          GEO_model_in = 0
       else if (equilibrium_model == 2 .or. equilibrium_model == 3) then
          GEO_model_in = geo_numeq_flag
       endif
       GEO_ntheta_in   = geo_ntheta
       GEO_nfourier_in = geo_ny
       call GEO_alloc(1)

       initialized = .true.
       
    else
       if(.NOT. initialized) return
       
       deallocate(theta)
       deallocate(theta_B)
       deallocate(w_theta)
       deallocate(Bmag)
       deallocate(k_perp)
       deallocate(omega_stream)
       deallocate(omega_trap)
       deallocate(omega_rdrift)
       deallocate(omega_adrift)
       deallocate(omega_aprdrift)

       call GEO_alloc(0)
        
       initialized = .false.

    endif
    
  end subroutine EQUIL_alloc
  
  subroutine EQUIL_do
    use cgyro_globals
    use GEO_interface
    implicit none
    integer :: it, ir, is
    real :: sum

    ! parameters needed for Miller equilibrium
    ! geo_numeq_flag, geo_ny, and geo_yin already set    
    GEO_rmin_in      = rmin
    GEO_rmaj_in      = rmaj
    GEO_drmaj_in     = shift
    GEO_zmag_in      = zmag
    GEO_dzmag_in     = s_zmag
    GEO_q_in         = q
    GEO_s_in         = s
    GEO_kappa_in     = kappa
    GEO_s_kappa_in   = s_kappa
    GEO_delta_in     = delta
    GEO_s_delta_in   = s_delta
    GEO_zeta_in      = zeta
    GEO_s_zeta_in    = s_zeta
    GEO_beta_star_in = beta_star
    GEO_fourier_in(:,:) = geo_yin(:,:)

    call GEO_do()  
    
    sum=0.0
    do it=1,n_theta

       call GEO_interp(theta(it))

       do is=1,n_species

          omega_stream(it,is) = sqrt(2.0) * vth(is) / (q * rmaj * GEO_g_theta)
          
          omega_trap(it,is) = -0.5*sqrt(2.0) * vth(is) &
               * (GEO_dbdt / GEO_b) / (q * rmaj * GEO_g_theta) 

          omega_rdrift(it,is) = -rho * vth(is)**2 * mass(is)/(Z(is)*GEO_b) &
               * GEO_grad_r / rmaj * GEO_gsin 

          omega_adrift(it,is) = -rho * vth(is)**2 * mass(is)/(Z(is)*GEO_b) &
               * GEO_gq / rmaj &
               * (GEO_gcos1 + GEO_gcos2 + GEO_captheta * GEO_gsin) 

          omega_aprdrift(it,is) = 2.0 * rho * vth(is)**2 &
               * mass(is)/(Z(is)*GEO_b) * GEO_gq / rmaj * GEO_gcos2

       enddo

       Bmag(it) = GEO_b

       ! flux-surface average weights
       w_theta(it) = GEO_g_theta / GEO_b
       sum = sum + w_theta(it)
       
       do ir=1,n_radial
          k_perp(it,ir) = sqrt((2.0*pi*px(ir)*GEO_grad_r/length &
               + k_theta*GEO_gq*GEO_captheta)**2 &
               + (k_theta*GEO_gq)**2) 
       enddo
       
    enddo
    
    do it=1,n_theta
       w_theta(it) = w_theta(it) / sum
    enddo
    
  end subroutine EQUIL_DO


end module cgyro_equilibrium
