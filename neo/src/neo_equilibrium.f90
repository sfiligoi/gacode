module neo_equilibrium
  
  implicit none
  
  public :: EQUIL_alloc, EQUIL_do
  
  ! equilibrium parameters (theta)
  real :: d_theta                                 ! delta theta
  ! local (th)
  real, dimension(:), allocatable :: k_par        ! bhat dot grad/a
  real, dimension(:), allocatable :: v_drift_x    ! radial curvature drift vel
  real, dimension(:), allocatable :: v_drift_th   ! theta  curvature drift vel
  real, dimension(:), allocatable :: gradr        ! | grad r|
  real, dimension(:), allocatable :: gradpar_gradr ! bhat dot grad/a | grad r|
  real, dimension(:), allocatable :: w_theta      ! flux surface avg weights
  real, dimension(:), allocatable :: Btor         ! B_t / Bunit
  real, dimension(:), allocatable :: Bpol         ! B_p / Bunit
  real, dimension(:), allocatable :: Bmag         ! B/Bunit
  real, dimension(:), allocatable :: Bmag_rderiv  ! (dB/dr) (a/Bunit)
  real, dimension(:), allocatable :: gradpar_Bmag ! bhat dot grad/a Bmag/Bunit
  real, dimension(:), allocatable :: bigR       ! R/a (global)
  real, dimension(:), allocatable :: bigR_rderiv       ! R/a (global)
  real, dimension(:), allocatable :: gradpar_bigR ! bhat dot grad/a R/a
  real, dimension(:), allocatable :: jacobln_rderiv ! 1/sqrt(g) dsqrt(g)/dr
  real                            :: bigR_th0    ! R/a at theta=0
  real                            :: bigR_th0_rderiv ! dR/dr at theta=0
  real                            :: gradr_th0   ! | grad r| at theta=0
  real                            :: Btor_th0    ! b_t/Bunit at theta=0
  real                            :: Bpol_th0    ! b_p/Bunit at theta=0
  real                            :: Bmag_th0    ! B/Bunit   at theta=0
  real                            :: Bmag_th0_rderiv  ! (dB/dr) (a/Bunit) th=0
  real                            :: ftrap       ! frac of trapped particles
  real, dimension(:), allocatable :: theta_nc    ! NCLASS theta grid
  
  ! radial deriv of the volume enclosed by a flux surface
  real, dimension(:), allocatable :: v_prime_g  ! (nr)
  
  ! Parameters needed for post-processing analysis
  real :: I_div_psip    ! I(psi)/psi_prime = q f / r
  real :: Bmag2_avg     ! <(B/Bunit)^2>
  real :: Bmag2inv_avg  ! <(Bunit/B)^2>
  real :: Btor2_avg     ! <(Btor/Bunit)^2>
  real :: bigRinv_avg   ! <a/R> 
  real :: gradpar_Bmag2_avg ! <(bhat dot grad/aBmag/Bunit)^2>
  real, dimension(:,:), allocatable :: geo_param
  
  logical, private :: initialized = .false.
  
contains
  
  subroutine EQUIL_alloc(flag)

    use neo_globals
  
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    integer :: it

    if(flag == 1) then
       if(initialized) return
       
       allocate(theta(n_theta))
       allocate(k_par(n_theta))
       allocate(v_drift_x(n_theta))
       allocate(v_drift_th(n_theta))
       allocate(gradr(n_theta))
       allocate(gradpar_gradr(n_theta))
       allocate(w_theta(n_theta))
       allocate(bigR(n_theta))
       allocate(bigR_rderiv(n_theta))
       allocate(gradpar_bigR(n_theta))
       allocate(Btor(n_theta))
       allocate(Bpol(n_theta))
       allocate(Bmag(n_theta))
       allocate(Bmag_rderiv(n_theta))
       allocate(gradpar_Bmag(n_theta))
       allocate(theta_nc(n_theta))
       allocate(jacobln_rderiv(n_theta))
       allocate(v_prime_g(n_radial))
       allocate(geo_param(n_radial,5))

       d_theta = (2*pi/n_theta)
       do it=1,n_theta
          theta(it) = -pi+(it-1)*d_theta
       enddo
       
       initialized = .true.
       
    else
       if(.NOT. initialized) return
       
       deallocate(theta)
       deallocate(k_par)
       deallocate(v_drift_x)
       deallocate(v_drift_th)
       deallocate(gradr)
       deallocate(gradpar_gradr)
       deallocate(w_theta)
       deallocate(bigR_rderiv)
       deallocate(gradpar_bigR)
       deallocate(bigR)
       deallocate(Btor)
       deallocate(Bpol)
       deallocate(Bmag)
       deallocate(Bmag_rderiv)
       deallocate(gradpar_Bmag)
       deallocate(theta_nc)
       deallocate(jacobln_rderiv)
       deallocate(v_prime_g)
       deallocate(geo_param)
              
       initialized = .false.

    endif
    
  end subroutine EQUIL_alloc
  
  subroutine EQUIL_do(ir)
    use neo_globals
    use geo
    implicit none
    integer, intent(in) :: ir
    integer :: it, jt, id, is
    real :: sum
    real, dimension(:), allocatable :: ttmp
    
    sum=0.0

    if (equilibrium_model == 2) then 

       GEO_signb_in   = sign_bunit
       GEO_rmin_in    = r(ir)
       GEO_rmaj_in    = rmaj(ir)
       GEO_drmaj_in   = shift(ir)
       GEO_zmag_in    = zmag(ir)
       GEO_dzmag_in   = s_zmag(ir)
       GEO_q_in       = q(ir)
       GEO_s_in       = shear(ir)
       GEO_kappa_in   = kappa(ir)
       GEO_s_kappa_in = s_kappa(ir)
       GEO_delta_in   = delta(ir)
       GEO_s_delta_in = s_delta(ir)
       GEO_zeta_in    = zeta(ir)
       GEO_s_zeta_in  = s_zeta(ir)
       GEO_shape_sin3_in    = shape_sin(3,ir)
       GEO_shape_s_sin3_in  = shape_s_sin(3,ir)
       GEO_shape_sin4_in    = shape_sin(4,ir)
       GEO_shape_s_sin4_in  = shape_s_sin(4,ir)
       GEO_shape_sin5_in    = shape_sin(5,ir)
       GEO_shape_s_sin5_in  = shape_s_sin(5,ir)
       GEO_shape_sin6_in    = shape_sin(6,ir)
       GEO_shape_s_sin6_in  = shape_s_sin(6,ir)
       GEO_shape_cos0_in    = shape_cos(0,ir)
       GEO_shape_s_cos0_in  = shape_s_cos(0,ir)
       GEO_shape_cos1_in    = shape_cos(1,ir)
       GEO_shape_s_cos1_in  = shape_s_cos(1,ir)
       GEO_shape_cos2_in    = shape_cos(2,ir)
       GEO_shape_s_cos2_in  = shape_s_cos(2,ir)
       GEO_shape_cos3_in    = shape_cos(3,ir)
       GEO_shape_s_cos3_in  = shape_s_cos(3,ir)
       GEO_shape_cos4_in    = shape_cos(4,ir)
       GEO_shape_s_cos4_in  = shape_s_cos(4,ir)
       GEO_shape_cos5_in    = shape_cos(5,ir)
       GEO_shape_s_cos5_in  = shape_s_cos(5,ir)
       GEO_shape_cos6_in    = shape_cos(6,ir)
       GEO_shape_s_cos6_in  = shape_s_cos(6,ir)
       
       !!!!! beta_star !!!!!
       ! EAB: does not enter in the first-order calc
       ! NOTE: it is not implemented in v_drift defs
       ! This input param is then only used in anisotropic species mode
       ! which currently only works in local profile mode
       GEO_beta_star_in = 0.0
       do is=1,n_species
          if(aniso_model(is) == 2) then
             GEO_beta_star_in = beta_star(ir)
          endif
       enddo

       ! GEO dimensions
       GEO_ntheta_in   = 2001
       GEO_model_in    = 0

       ! Get initial geo solution, then set geo params at theta=0
       allocate(ttmp(1))
       ttmp(1) = 0.0
       call geo_interp(1,ttmp,.true.)
       deallocate(ttmp)
       bigR_th0        = GEO_bigr(1)
       bigR_th0_rderiv = GEO_bigr_r(1)
       gradr_th0       = GEO_grad_r(1)
       Btor_th0        = GEO_bt(1)
       Bpol_th0        = GEO_bp(1)
       Bmag_th0        = GEO_b(1)
       Bmag_th0_rderiv = -GEO_b(1)/(rmaj(ir)*GEO_grad_r(1)) &
               * (GEO_gcos1(1) + GEO_gcos2(1))
            
       call geo_interp(n_theta,theta,.false.)

       do it=1,n_theta
          
          k_par(it)     = 1.0 / (q(ir) * rmaj(ir) * GEO_g_theta(it))
          w_theta(it)   = GEO_g_theta(it) / GEO_b(it)
          sum = sum + w_theta(it)
          
          bigR(it) = GEO_bigr(it)
          bigR_rderiv(it) = GEO_bigr_r(it)
          gradpar_bigR(it) = k_par(it) * GEO_bigr_t(it) 
          Bmag(it)  = GEO_b(it)
          Btor(it)  = GEO_bt(it)
          Bpol(it)  = GEO_bp(it)
          Bmag_rderiv(it)  = -GEO_b(it)/(rmaj(ir)*GEO_grad_r(it)) &
               * (GEO_gcos1(it) + GEO_gcos2(it))
          gradpar_Bmag(it) = k_par(it) * GEO_dbdt(it)
          gradr(it)        = GEO_grad_r(it)
          v_drift_x(it)  = -rho(ir)/(rmaj(ir) * Bmag(it)) * &
               GEO_grad_r(it) * GEO_gsin(it)
          v_drift_th(it) = -rho(ir)/(rmaj(ir) * Bmag(it) * GEO_l_t(it)) &
               * (GEO_gcos1(it) * GEO_bt(it) / Bmag(it)  & 
               -  GEO_gsin(it) * GEO_nsin(it) * GEO_grad_r(it)) &
               * r(ir)
          theta_nc(it) = GEO_theta_nc(it)
          jacobln_rderiv(it) = 2.0/bigR(it)*bigR_rderiv(it) &
               - 1.0/Btor(it) * r(ir)/(q(ir)*bigR(it)) * GEO_ffprime/GEO_f
         
       enddo

       I_div_psip = GEO_f * q(ir) / r(ir)

       do it=1,n_theta
          gradpar_gradr(it)  = 0.0
          do id=-2,2
             jt = thcyc(it+id)
             gradpar_gradr(it) = gradpar_gradr(it) &
                  + gradr(jt) * cderiv(id) / (12.0*d_theta)
          enddo
          gradpar_gradr(it) = gradpar_gradr(it) * k_par(it)
       enddo
          
       ! 1/J dJ/dr
       allocate(ttmp(n_theta))
       ttmp(:) = theta(:)
       GEO_rmin_in = r(ir) + 0.001
       call geo_interp(n_theta,ttmp,.true.)
       do it=1,n_theta
          ttmp(it) = GEO_g_theta(it)/GEO_b(it)
       enddo
       jacobln_rderiv(:) = 0.0
       do it=1,n_theta
          jacobln_rderiv(it) = (ttmp(it) - w_theta(it))/0.001
       enddo
       jacobln_rderiv(:) = jacobln_rderiv(:) / w_theta(:)
       deallocate(ttmp)
       
    else
       
       ! concentric circular geometry
       shift(ir)   = 0.0
       kappa(ir)   = 1.0
       s_kappa(ir) = 0.0
       delta(ir)   = 0.0
       s_delta(ir) = 0.0

       do it=1,n_theta
          k_par(it)       = 1.0 / (q(ir) * rmaj(ir)) * sign_bunit
          bigR(it)        = rmaj(ir) + r(ir) * cos(theta(it))
          bigR_rderiv(it) = cos(theta(it))
          gradpar_bigR(it) = -r(ir) * sin(theta(it)) * k_par(it)
          gradr(it)          = 1.0
          gradpar_gradr(it)  = 0.0
          theta_nc(it)    = theta(it)
          bigR_th0        = rmaj(ir) + r(ir)
          bigR_th0_rderiv = 1.0 
          gradr_th0       = 1.0

          if (equilibrium_model == 1) then
             ! large aspect ratio geometry
             Btor_th0 = 1.0 - (r(ir)/rmaj(ir))
             Bpol_th0 = r(ir) / (q(ir)*rmaj(ir)) * (1.0 - (r(ir)/rmaj(ir)))
             Bmag_th0 = (1.0 - (r(ir)/rmaj(ir))) * sign_bunit
             Bmag_th0_rderiv = - 1.0/rmaj(ir) * sign_bunit
             Bmag(it) = (1.0 - (r(ir)/rmaj(ir)) * cos(theta(it))) * sign_bunit
             Btor(it) = 1.0 - (r(ir)/rmaj(ir)) * cos(theta(it))
             Bpol(it) = r(ir) / (q(ir)*rmaj(ir)) &
                  * (1.0 - (r(ir)/rmaj(ir)) * cos(theta(it)))
             Bmag_rderiv(it) = - 1.0/rmaj(ir) * cos(theta(it)) * sign_bunit
             gradpar_Bmag(it) = k_par(it) * (r(ir)/rmaj(ir)) &
                  * sin(theta(it)) * sign_bunit
             v_drift_x(it) = -rho(ir)/rmaj(ir) * sin(theta(it)) &
                  / (1.0 - (r(ir)/rmaj(ir)) * cos(theta(it)))**2
             v_drift_th(it) = -rho(ir)/(rmaj(ir)) * cos(theta(it)) &
                  / (1.0 - (r(ir)/rmaj(ir)) * cos(theta(it)))**2

          else
             ! s-alpha geometry
             Btor_th0 = 1.0 / (1.0 + (r(ir)/rmaj(ir)))
             Bpol_th0 = r(ir) / (q(ir)*rmaj(ir)) &
                  / (1.0 + (r(ir)/rmaj(ir)))
             Bmag_th0 = (1.0 / (1.0 + (r(ir)/rmaj(ir)))) &
                  * sign_bunit
             Bmag_th0_rderiv = (1.0 / (1.0 + (r(ir)/rmaj(ir)))**2) &
                  * sign_bunit * (-1.0/rmaj(ir))
             Bmag(it) = (1.0 / (1.0 + (r(ir)/rmaj(ir)) * cos(theta(it)))) &
                  * sign_bunit
             Btor(it) = 1.0 / (1.0 + (r(ir)/rmaj(ir)) * cos(theta(it)))
             Bpol(it) = r(ir) / (q(ir)*rmaj(ir)) &
                  / (1.0 + (r(ir)/rmaj(ir)) * cos(theta(it)))
             Bmag_rderiv(it) = (1.0 / (1.0 &
                  + (r(ir)/rmaj(ir))* cos(theta(it)))**2) &
                  * sign_bunit * (-1.0/rmaj(ir) * cos(theta(it)))
             gradpar_Bmag(it) = k_par(it) * (r(ir)/rmaj(ir)) &
                  * sin(theta(it)) / (1.0 + (r(ir)/rmaj(ir)) &
                  * cos(theta(it)))**2 * sign_bunit
             v_drift_x(it)  = -rho(ir)/rmaj(ir) * sin(theta(it))
             v_drift_th(it) = -rho(ir)/rmaj(ir) * cos(theta(it))
          endif

          ! 1/J dJ/dr
          jacobln_rderiv(it) = -(1.0/Bmag(it))*Bmag_rderiv(it)
          
          ! flux-surface average weights
          w_theta(it) = 1.0 * sign_bunit / Bmag(it)
          sum = sum + w_theta(it)
       enddo

       I_div_psip = rmaj(ir) * q(ir) / r(ir)

    endif
    
    do it=1,n_theta
       w_theta(it) = w_theta(it) / sum
    enddo

    ! v_prime_g = d Volume of flux surface / dr
    v_prime_g(ir) = r(ir)*rmaj(ir)*sum/n_theta * 4 * pi * pi

    ! Flux surface averages of B^2 and 1/B^2
    Bmag2_avg         = 0.0
    Bmag2inv_avg      = 0.0
    gradpar_Bmag2_avg = 0.0
    Btor2_avg         = 0.0
    bigRinv_avg       = 0.0
    do it=1,n_theta
       Bmag2_avg    = Bmag2_avg    + w_theta(it)*Bmag(it)**2
       Bmag2inv_avg = Bmag2inv_avg + w_theta(it)/Bmag(it)**2
       gradpar_Bmag2_avg = gradpar_Bmag2_avg + w_theta(it)*gradpar_Bmag(it)**2
       Btor2_avg = Btor2_avg + w_theta(it)*Btor(it)**2
       bigRinv_avg =  bigRinv_avg + w_theta(it) / bigR(it)
    enddo
    
    call compute_fractrap(ftrap)

    if(silent_flag == 0 .and. i_proc == 0) then
       open(unit=1,file=trim(path)//'out.neo.diagnostic_geo',status='replace')
       write(1,'(a,1pe16.8)') "# I/psi'            = ",I_div_psip
       write(1,'(a,1pe16.8)') "# <1/B^2>-1/<B^2>   = ",Bmag2inv_avg-1.0/Bmag2_avg
       write(1,'(a,1pe16.8)') "# f_trap            = ",ftrap
       write(1,'(a,i3)')      "# n_theta           = ",n_theta
       write(1,'(a)') "# Functions:" 
       write(1,'(a)') "#   theta(:)" 
       write(1,'(a)') "#   v_drift_x(:)" 
       write(1,'(a)') "#   gradpar_Bmag(:)" 
       write(1,'(a)') "#   Bmag(:)" 
       write(1,'(a)') "#   w_theta(:)" 
       write(1,'(a)') "#   R(:)"
       write(1,'(a)') "#   R(theta=0)" 
       write(1,'(a)') "#   dR(theta=0)/dr" 
       write(1,'(1pe16.8)') theta(:)
       write(1,'(1pe16.8)') v_drift_x(:)
       write(1,'(1pe16.8)') gradpar_Bmag(:)
       write(1,'(1pe16.8)') Bmag(:)
       write(1,'(1pe16.8)') w_theta(:)
       write(1,'(1pe16.8)') bigR(:)
       write(1,'(1pe16.8)') bigR_th0
       write(1,'(1pe16.8)') bigR_th0_rderiv
       close(1)
    endif
    
    neo_geo_out(1)  = I_div_psip
    neo_geo_out(2)  = ftrap
    neo_geo_out(3)  = Bmag2_avg
    neo_geo_out(4)  = Bpol_th0
    neo_geo_out(5)  = (Bmag2inv_avg-1.0/Bmag2_avg)*Bmag2_avg

    if(silent_flag == 0 .and. i_proc == 0) then
       open(unit=1,file=trim(path)//'out.neo.diagnostic_geo2',status='replace')
       do it=1,5
          write(1,'(1pe16.8)') neo_geo_out(it)
       enddo
    endif

    geo_param(ir,:) = neo_geo_out(:)
    
  end subroutine EQUIL_DO
  

  ! Computes the fraction of trapped particles
  subroutine compute_fractrap(ft)
    use neo_globals
    implicit none
    real, intent(out) :: ft
    integer :: nlambda=500
    real, dimension(:), allocatable :: lambda
    real :: dlambda, fac_lambda, sum_th
    integer i,it
    real :: Bmax
    
    Bmax = Bmag(1) * sign_bunit
    do it=2, n_theta
       if(Bmag(it)*sign_bunit > Bmax) then
          Bmax = Bmag(it) * sign_bunit
       endif
    enddo

    allocate(lambda(nlambda))
    dlambda = 1.0/(nlambda-1)
    do i=1, nlambda
       lambda(i) = (i-1)*dlambda
    enddo

    ft = 0.0
    
    do i=2, nlambda-1
       
       ! open integration for lambda
       if(i==2 .or. i == nlambda-1)  then
          fac_lambda = 23.0/12.0
       else if(i==3 .or. i == nlambda-2) then
          fac_lambda = 7.0/12.0
       else
          fac_lambda = 1.0
       endif
       
       ! closed integration for th: <sqrt(1-lambda*B/Bmax)>
       sum_th = 0.0
       do it=1, n_theta
          sum_th = sum_th + w_theta(it) &
               * sqrt(1.0 - lambda(i) * Bmag(it) * sign_bunit / Bmax)
       enddo
       
       ft = ft + fac_lambda * lambda(i) / sum_th;
       
    end do
    

    ft = ft * dlambda * 0.75 * Bmag2_avg / Bmax**2
    ft = 1.0 - ft
    
    deallocate(lambda)
    
  end subroutine compute_fractrap


end module neo_equilibrium
