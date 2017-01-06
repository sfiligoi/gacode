subroutine cgyro_init_rotation

  use timer_lib
  use cgyro_globals
  use cgyro_io
  
  implicit none

  real, dimension(:), allocatable :: phi_rot, phi_rot_deriv
  real :: phi_rot_avg
  integer, dimension(:), allocatable :: thcyc
  real, dimension(:), allocatable :: thcderiv
  real :: x, x0, fac, sum_zn, dsum_zn
  integer :: is, it, j, id, jt
  integer, parameter :: jmax = 200

  allocate(phi_rot(n_theta))
  allocate(phi_rot_deriv(n_theta))

  if(cf_flag == 0) then
     dens_rot(:,:) = 1.0
     phi_rot(:)  = 0.0
     phi_rot_deriv(:)  = 0.0
     lambda_rot(:,:) = 0.0
     dlambda_rot(:,:) = 0.0
     omega_rot_trap(:,:) = 0.0
     omega_rot_u(:,:) = 0.0
     omega_rot_drift(:,:) = 0.0
     omega_rot_drift_r(:,:) = 0.0
     omega_rot_star(:,:) = 0.0
     omega_rot_drift_e(:,:) = 0.0
     return
  endif

  ! solve the quasi-neutrality relation for poloidal part of phi

  ! partial component of n/n(theta0) 
  ! -- will add the phi_rot component after the QN solve
  do is=1, n_species
     do it=1,n_theta
        dens_rot(it,is) = exp(0.5 * (mach/vth(is))**2 &
             * (bigR(it)**2 - bigR_th0**2) / rmaj**2) 
     enddo
  enddo
  
  phi_rot_avg = 0.0
  x = 0.05          ! initial guess for phi_rot(1)
  do it=1,n_theta
     
     j=1
     
     do
        ! use Newton's method to solve the quasi-neutrality relation
        
        sum_zn  = 0.0
        dsum_zn = 0.0
        do is=1, n_species
           fac = z(is) * dens(is) * dens_rot(it,is) &
                * exp(-x * z(is)/temp(is))
           sum_zn  = sum_zn  + fac
           dsum_zn = dsum_zn - fac * z(is)/temp(is) 
        enddo
        
        if(ae_flag == 1) then
           fac = -dens_ele * exp(x/temp_ele)
           sum_zn  = sum_zn  + fac
           dsum_zn = dsum_zn + fac/temp_ele 
        endif

        x0 = x
        x  = x0 - sum_zn / dsum_zn

        if (abs(x-x0) < 1.0e-12) exit
        
        j = j + 1
        if(j > jmax) exit
        
     enddo

     if(j > jmax) then
        call cgyro_error('ERROR: (CGYRO) Rotation poloidal density computation failed to converge')
        return
     endif
     
     ! phi_rot = [phi(theta) - phi(theta0)]
     phi_rot(it) = x

     ! <phi_rot>
     phi_rot_avg = phi_rot_avg + w_theta(it) * phi_rot(it)

  enddo

  ! d phi_rot/d_theta
  
  ! cyclic index (for theta-periodicity)
  allocate(thcyc(1-n_theta:2*n_theta))
  do it=1,n_theta
     thcyc(it-n_theta) = it
     thcyc(it) = it
     thcyc(it+n_theta) = it
  enddo

  ! 4th-order centered derivative
  allocate(thcderiv(-2:2))
  thcderiv(-2) =  1.0 / (12.0 * d_theta)
  thcderiv(-1) = -8.0 / (12.0 * d_theta)
  thcderiv(0)  =  0.0 / (12.0 * d_theta)
  thcderiv(1)  =  8.0 / (12.0 * d_theta)
  thcderiv(2)  = -1.0 / (12.0 * d_theta)
  
  do it=1,n_theta
     phi_rot_deriv(it) = 0.0
     do id=-nup_theta,nup_theta
        jt = thcyc(it+id)
        phi_rot_deriv(it) = phi_rot_deriv(it) &
             + phi_rot(jt) * thcderiv(id) 
     enddo
  enddo
  deallocate(thcyc)
  deallocate(thcderiv)
  
  do is=1,n_species
     do it=1,n_theta

        ! n(theta)/n(theta0)
        dens_rot(it,is) = dens_rot(it,is) * exp(-z(is)/temp(is)*phi_rot(it))
        
        lambda_rot(it,is) = z(is)/temp(is) * (phi_rot(it) - phi_rot_avg) &
             - 0.5 * (mach * bigR(it) / rmaj / vth(is))**2
        
        dlambda_rot(it,is) = dlambda_rot(it,is) &
             + z(is)/temp(is) * phi_rot_deriv(it) / g_theta(it) 
           
        omega_rot_trap(it,is) = -0.5*sqrt(2.0)*vth(is) &
             /(q*rmaj)*dlambda_rot(it,is)

        omega_rot_u(it,is) = -vth(is)/(sqrt(2.0)*q*rmaj) *dlambda_rot(it,is)

        omega_rot_star(it,is) = omega_rot_star(it,is) &
             + dlntdr(is)*z(is)/temp(is)*phi_rot(it) 

        omega_rot_drift_e(it,is) = 0.0  ! EAB: NOT YET SURE ABOUT THIS
        
        
     enddo
  enddo

  deallocate(phi_rot)
  deallocate(phi_rot_deriv)
  
end subroutine cgyro_init_rotation
