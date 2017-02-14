subroutine cgyro_init_rotation

  use timer_lib
  use cgyro_globals
  use cgyro_io
  
  implicit none

  real, dimension(:), allocatable :: phi_rot, phi_rot_tderiv, phi_rot_rderiv
  real, dimension(:), allocatable :: sum_pressure_t
  real :: phi_rot_avg
  integer, dimension(:), allocatable :: thcyc
  real, dimension(:), allocatable :: thcderiv
  real :: x, x0, fac, sum_zn, dsum_zn
  integer :: is, it, j, id, jt
  integer, parameter :: jmax = 200

  if(cf_model == 0) then
     dens_rot(:,:) = 1.0
     dens_ele_rot(:) = 1.0
     lambda_rot(:,:) = 0.0
     dlambda_rot(:,:) = 0.0
     omega_rot_trap(:,:) = 0.0
     omega_rot_u(:,:) = 0.0
     omega_rot_drift(:,:) = 0.0
     omega_rot_drift_r(:,:) = 0.0
     omega_rot_prdrift(:,:) = 0.0
     omega_rot_prdrift_r(:,:) = 0.0
     omega_rot_star(:,:) = 0.0
     omega_rot_edrift(:,:) = 0.0
     omega_rot_edrift_r(:,:) = 0.0
     omega_rot_edrift_0(:) = 0.0
     return
  endif

  allocate(phi_rot(n_theta))
  allocate(phi_rot_tderiv(n_theta))
  allocate(phi_rot_rderiv(n_theta))
  allocate(sum_pressure_t(n_theta))
  
  ! solve the quasi-neutrality relation for poloidal part of phi

  ! partial component of n/n(theta0) 
  ! -- will add the phi_rot component after the QN solve
  do is=1, n_species
     do it=1,n_theta
        dens_rot(it,is) = exp(0.5 * (mach/vth(is))**2 &
             * (bigR(it)**2 - bigR_th0**2) / rmaj**2) 
     enddo
  enddo

  if(ae_flag == 1) then
     dens_ele_rot(:) = 1.0
  else
     dens_ele_rot(:) = exp(0.5 * (mach/vth(is_ele))**2 &
          * (bigR(:)**2 - bigR_th0**2) / rmaj**2)
  endif
     
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
           fac = -dens_ele * dens_ele_rot(it) * exp(x/temp_ele)
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

  ! d phi_rot / d theta
  phi_rot_tderiv(:) = 0.0
  do it=1,n_theta
     do id=-2,2
        jt = thcyc(it+id)
        phi_rot_tderiv(it) = phi_rot_tderiv(it) &
             + phi_rot(jt) * thcderiv(id) 
     enddo
  enddo
  deallocate(thcyc)
  deallocate(thcderiv)
  
  ! n(theta)/n(theta0)
  do is=1,n_species
     do it=1,n_theta
        dens_rot(it,is) = dens_rot(it,is) * exp(-z(is)/temp(is)*phi_rot(it))
     enddo
  enddo
  dens_ele_rot(:) = dens_ele_rot(:) * exp(1.0/temp_ele*phi_rot(:))
  
  sum_pressure_t(:) = 0.0
  do is=1,n_species
     do it=1,n_theta

        lambda_rot(it,is) = z(is)/temp(is) * (phi_rot(it) - phi_rot_avg) &
             - 0.5 * (mach * bigR(it) / rmaj / vth(is))**2

        ! bhat dot grad lambda
        dlambda_rot(it,is) = dlambda_rot(it,is) &
             + z(is)/temp(is) * phi_rot_tderiv(it) / (q*rmaj*g_theta(it)) 
           
        omega_rot_trap(it,is) = -0.5*sqrt(2.0)*vth(is) *dlambda_rot(it,is)

        omega_rot_u(it,is) = -vth(is)/sqrt(2.0)*dlambda_rot(it,is)

        omega_rot_star(it,is) = omega_rot_star(it,is) &
             + dlntdr(is)*z(is)/temp(is)*phi_rot(it) 

        omega_rot_edrift(it,is) =  omega_rot_edrift(it,is) &
             * phi_rot_tderiv(it) / (q*rmaj*g_theta(it)) 

        omega_rot_edrift_r(it,is) = omega_rot_edrift_r(it,is) &
             * phi_rot_tderiv(it) / (q*rmaj*g_theta(it))
        
        ! bhat dot grad pressure
        sum_pressure_t(it) = sum_pressure_t(it) - dens(is)*temp(is) &
             * dlambda_rot(it,is)
     enddo
  enddo  
  
  ! 1/(ne(0)Te) bhat dot grad pressure
  sum_pressure_t(:) = sum_pressure_t(:)/(dens_ele*temp_ele)
  
  do is=1,n_species
     do it=1,n_theta
        omega_rot_prdrift(it,is) = omega_rot_prdrift(it,is) &
             * sum_pressure_t(it)
        
        omega_rot_prdrift_r(it,is) = omega_rot_prdrift_r(it,is) &
             * sum_pressure_t(it)
     enddo
  enddo

  ! Solve d/dr of QN to get d phi_rot / dr

  do it=1,n_theta
     
     ! sum_a z^2/T n(theta)
     sum_zn = 0.0
     do is=1,n_species
        sum_zn = sum_zn + z(is)*z(is)/temp(is)*dens(is)*dens_rot(it,is)
     enddo
     if(ae_flag == 1) then
        sum_zn = sum_zn + 1.0/temp_ele*dens_ele*dens_ele_rot(it)
     endif

     phi_rot_rderiv(it) = 0.0
     do is=1,n_species
        phi_rot_rderiv(it) = phi_rot_rderiv(it) &
             + z(is)*dens(is)*dens_rot(it,is)*(-dlnndr(is) &
             - dlntdr(is)*(z(is)/temp(is)*phi_rot(it) &
             - 0.5*(mach/rmaj/vth(is))**2 * (bigR(it)**2 - bigR_th0**2)) &
             + omega_rot_edrift_0(it)/vth(is)**2)
     enddo
     if(ae_flag == 1) then
        phi_rot_rderiv(it) = phi_rot_rderiv(it) &
             - dens_ele*dens_ele_rot(it)*(-dlnndr_ele &
             + dlntdr_ele/temp_ele*phi_rot(it))
     endif
     
     phi_rot_rderiv(it) = phi_rot_rderiv(it) / sum_zn
     
  enddo
  
  ! rho * dphi_*/dr
  do it=1,n_theta
     omega_rot_edrift_0(it) = rho*phi_rot_rderiv(it)
  enddo


  ! print out some diagnostics
  if (silent_flag == 0 .and. i_proc == 0) then
     open(unit=io,file=trim(path)//'out.cgyro.rotation',status='replace')
     do it=1,n_theta
        write(io,'(4(es16.9,1x))') theta(it), phi_rot(it), &
             phi_rot_tderiv(it), phi_rot_rderiv(it)
     enddo
     close(io)
  endif

  ! just cf trapping (no coriolis and no cf drift)
  if(cf_model == 2) then

     omega_cdrift(:,:) = 0.0
     omega_cdrift_r(:,:) = 0.0
     
     omega_rot_drift(:,:) = 0.0
     omega_rot_drift_r(:,:) = 0.0

  endif

  ! just cf drift (no coriolis and no cf trap)
  if(cf_model == 3) then

     omega_cdrift(:,:) = 0.0
     omega_cdrift_r(:,:) = 0.0
     
     dens_rot(:,:) = 1.0
     dens_ele_rot(:) = 1.0
     phi_rot(:)  = 0.0
     phi_rot_tderiv(:)  = 0.0
     phi_rot_rderiv(:)  = 0.0
     lambda_rot(:,:) = 0.0
     dlambda_rot(:,:) = 0.0
     omega_rot_trap(:,:) = 0.0
     omega_rot_u(:,:) = 0.0
     omega_rot_prdrift(:,:) = 0.0
     omega_rot_prdrift_r(:,:) = 0.0
     omega_rot_star(:,:) = 0.0
     omega_rot_edrift(:,:) = 0.0
     omega_rot_edrift_r(:,:) = 0.0
     omega_rot_edrift_0(:) = 0.0
     
  endif

  ! just cf mach^2 terms (no coriolis and no parallel velocity shear)
  if(cf_model == 4) then

     omega_cdrift(:,:) = 0.0
     omega_cdrift_r(:,:) = 0.0
     
     omega_gammap(:) = 0.0

  endif
  
  deallocate(phi_rot)
  deallocate(phi_rot_tderiv)
  deallocate(phi_rot_rderiv)
  deallocate(sum_pressure_t)
  
end subroutine cgyro_init_rotation
