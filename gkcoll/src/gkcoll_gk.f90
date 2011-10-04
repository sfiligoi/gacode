module gkcoll_gk
  implicit none

  public :: GK_alloc, GK_init, GK_do

  logical, private :: initialized = .false.

  complex, dimension(:,:,:,:,:,:), allocatable, private :: rhs
  complex, dimension(:,:,:,:,:), allocatable, private :: h0_x
  complex, dimension(:,:,:,:,:), allocatable, private :: cap_h_x_deriv

  ! theta derivative variables
  real, dimension(-2:2) :: uderiv
  real, dimension(-2:2) :: cderiv
  integer, dimension(:), allocatable :: thcyc

  ! xi conversion matrices
  complex, dimension(:,:), allocatable :: xi_mat, xi_mat_inv, xi_deriv_mat
  complex, private :: alpha = (1.0,0.0)
  complex, private :: beta  = (0.0,0.0)

contains

  subroutine GK_alloc(flag)
    use gkcoll_globals
    use gkcoll_legendre
    use gkcoll_equilibrium, only : d_theta
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    integer :: it, ix, jx
    real :: val
    ! LAPACK
    integer :: info
    integer, dimension(:), allocatable :: i_piv
    complex, dimension(:), allocatable :: work

    if(flag == 1) then
       if(initialized) return

       allocate(rhs(4,n_species,n_radial,n_theta,n_energy,n_xi))
       allocate(h0_x(n_species,n_radial,n_theta,n_energy,n_xi))
       allocate(cap_h_x_deriv(n_species,n_radial,n_theta,n_energy,n_xi))

       ! cyclic index (for theta-periodicity)
       allocate(thcyc(1-n_theta:2*n_theta))
       do it=1,n_theta
          thcyc(it-n_theta) = it
          thcyc(it) = it
          thcyc(it+n_theta) = it
       enddo
       ! coefficients for 4th order centered derivative
       cderiv(-2) =  1.0 / (12.0 * d_theta)
       cderiv(-1) = -8.0 / (12.0 * d_theta)
       cderiv(0)  =  0.0 / (12.0 * d_theta)
       cderiv(1)  =  8.0 / (12.0 * d_theta)
       cderiv(2)  = -1.0 / (12.0 * d_theta)
       ! coefficients for 4th order filter for 3rd order upwinded derivative
       uderiv(-2) =  1.0 / (12.0 * d_theta)
       uderiv(-1) = -4.0 / (12.0 * d_theta)
       uderiv(0)  =  6.0 / (12.0 * d_theta)
       uderiv(1)  = -4.0 / (12.0 * d_theta)
       uderiv(2)  =  1.0 / (12.0 * d_theta)

       ! xi grid to legendre conversion matrices
       allocate(xi_mat(n_xi,n_xi))
       allocate(xi_mat_inv(n_xi,n_xi))
       allocate(xi_deriv_mat(n_xi,n_xi))
       do ix=1,n_xi
          do jx=1,n_xi
             call legendre(indx_xi(jx),xi(ix),val)
             xi_mat(ix,jx) = val
             call legendre_deriv(indx_xi(jx),xi(ix),val)
             xi_deriv_mat(ix,jx) = val
          enddo
       enddo
       xi_mat_inv = xi_mat
       allocate(work(n_xi))
       allocate(i_piv(n_xi))
       call ZGETRF(n_xi,n_xi,xi_mat_inv,n_xi,i_piv,info)
       call ZGETRI(n_xi,xi_mat_inv,n_xi,i_piv,work,n_xi,info)
       deallocate(work)
       deallocate(i_piv)

       initialized = .true.

    else
       if(.NOT. initialized) return

       deallocate(rhs)
       deallocate(h0_x)
       deallocate(cap_h_x_deriv)
       deallocate(thcyc)
       deallocate(xi_mat)
       deallocate(xi_mat_inv)
       deallocate(xi_deriv_mat)

       initialized = .false.
    endif

  end subroutine GK_alloc
  
  subroutine GK_init
    use gkcoll_globals
    use gkcoll_poisson
    use gkcoll_gyro
    implicit none
    integer :: is,ir,it,ie,ix

    phi_old(:,:) = (0.0,0.0)
    h_x(:,:,:,:,:) = (0.0,0.0)
    do it=1,n_theta
       h_x(1,n_radial/2+1,it,:,:) = (1.0e-3) * (cos(theta(it)/2.0))**2
    enddo
    do ir=1,n_radial
       do it=1,n_theta
          call POISSONx_do(ir,it)
       enddo
    enddo
    do is=1,n_species
       do ir=1,n_radial
          do it=1,n_theta
             do ie=1,n_energy
                do ix=1,n_xi
                   cap_h_x(is,ir,it,ie,ix) = h_x(is,ir,it,ie,ix) &
                        + z(is)/temp(is) * gyrox_J0(is,ir,it,ie,ix) &
                        * phi(ir,it)
                enddo
                call ZGEMV('N',n_xi,n_xi,alpha,xi_mat_inv,n_xi,&
                     cap_h_x(is,ir,it,ie,:),1,beta,cap_h_p(is,ir,it,ie,:),1)
             enddo
          enddo
       enddo
    enddo

  end subroutine GK_init

  subroutine GK_do
    use gkcoll_globals
    use gkcoll_poisson
    use gkcoll_gyro
    implicit none
    integer :: is,ir,it,ie,ix
    
    ! compute new collisionless little_h: h = H - ze/T G phi
    ! RK4
    
    h0_x = h_x
    
    !!!!!!!!!!
    ! Stage 1
    !!!!!!!!!!

    do is=1,n_species
       do ir=1,n_radial
          do it=1,n_theta
             do ie=1,n_energy
                call ZGEMV('N',n_xi,n_xi,alpha,xi_mat,n_xi,&
                     cap_h_p(is,ir,it,ie,:),1,beta,cap_h_x(is,ir,it,ie,:),1)
                call ZGEMV('N',n_xi,n_xi,alpha,xi_deriv_mat,n_xi,&
                     cap_h_p(is,ir,it,ie,:),1,beta,&
                     cap_h_x_deriv(is,ir,it,ie,:),1)
             enddo
          enddo
       enddo
    enddo
    
    call get_gkRHS(1)
    h_x = h0_x + 0.5 * delta_t * rhs(1,:,:,:,:,:)
    do ir=1,n_radial
       do it=1,n_theta
          call POISSONx_do(ir,it)
       enddo
    enddo
    
    !!!!!!!!!!
    ! Stage 2
    !!!!!!!!!!

    do is=1,n_species
       do ir=1,n_radial
          do it=1,n_theta
             do ie=1,n_energy
                do ix=1,n_xi
                   cap_h_x(is,ir,it,ie,ix) = h_x(is,ir,it,ie,ix) &
                        + z(is)/temp(is) * gyrox_J0(is,ir,it,ie,ix) &
                        * phi(ir,it)
                enddo
                call ZGEMV('N',n_xi,n_xi,alpha,xi_mat_inv,n_xi,&
                     cap_h_x(is,ir,it,ie,:),1,beta,cap_h_p(is,ir,it,ie,:),1)
                call ZGEMV('N',n_xi,n_xi,alpha,xi_deriv_mat,n_xi,&
                     cap_h_p(is,ir,it,ie,:),1,beta,&
                     cap_h_x_deriv(is,ir,it,ie,:),1)
             enddo
          enddo
       enddo
    enddo
    
    call get_gkRHS(2)
    h_x = h0_x + 0.5 * delta_t * rhs(2,:,:,:,:,:)
    do ir=1,n_radial
       do it=1,n_theta
          call POISSONx_do(ir,it)
       enddo
    enddo
    
    !!!!!!!!!!
    ! Stage 3
    !!!!!!!!!!

    do is=1,n_species
       do ir=1,n_radial
          do it=1,n_theta
             do ie=1,n_energy
                do ix=1,n_xi
                   cap_h_x(is,ir,it,ie,ix) = h_x(is,ir,it,ie,ix) &
                        + z(is)/temp(is) * gyrox_J0(is,ir,it,ie,ix) &
                        * phi(ir,it)
                enddo
                call ZGEMV('N',n_xi,n_xi,alpha,xi_mat_inv,n_xi,&
                     cap_h_x(is,ir,it,ie,:),1,beta,cap_h_p(is,ir,it,ie,:),1)
                call ZGEMV('N',n_xi,n_xi,alpha,xi_deriv_mat,n_xi,&
                     cap_h_p(is,ir,it,ie,:),1,beta,&
                     cap_h_x_deriv(is,ir,it,ie,:),1)
             enddo
          enddo
       enddo
    enddo
    
    call get_gkRHS(3)
    h_x = h0_x + delta_t * rhs(3,:,:,:,:,:)
    do ir=1,n_radial
       do it=1,n_theta
          call POISSONx_do(ir,it)
       enddo
    enddo
    
    !!!!!!!!!!
    ! Stage 4
    !!!!!!!!!!

    do is=1,n_species
       do ir=1,n_radial
          do it=1,n_theta
             do ie=1,n_energy
                do ix=1,n_xi
                   cap_h_x(is,ir,it,ie,ix) = h_x(is,ir,it,ie,ix) &
                        + z(is)/temp(is) * gyrox_J0(is,ir,it,ie,ix) &
                        * phi(ir,it)
                enddo
                call ZGEMV('N',n_xi,n_xi,alpha,xi_mat_inv,n_xi,&
                     cap_h_x(is,ir,it,ie,:),1,beta,cap_h_p(is,ir,it,ie,:),1)
                call ZGEMV('N',n_xi,n_xi,alpha,xi_deriv_mat,n_xi,&
                     cap_h_p(is,ir,it,ie,:),1,beta,&
                     cap_h_x_deriv(is,ir,it,ie,:),1)
             enddo
          enddo
       enddo
    enddo
    
    call get_gkRHS(4)
    
    h_x = h0_x + delta_t/6.0 * (rhs(1,:,:,:,:,:) + 2.0*rhs(2,:,:,:,:,:) &
         + 2.0*rhs(3,:,:,:,:,:) + rhs(4,:,:,:,:,:))     
    do ir=1,n_radial
       do it=1,n_theta
          call POISSONx_do(ir,it)
       enddo
    enddo
    
    ! Re-map to get new cap_h
    do is=1,n_species
       do ir=1,n_radial
          do it=1,n_theta
             do ie=1,n_energy
                do ix=1,n_xi
                   cap_h_x(is,ir,it,ie,ix) = h_x(is,ir,it,ie,ix) &
                        + z(is)/temp(is) * gyrox_J0(is,ir,it,ie,ix) &
                        * phi(ir,it)
                enddo
                call ZGEMV('N',n_xi,n_xi,alpha,xi_mat_inv,n_xi,&
                     cap_h_x(is,ir,it,ie,:),1,beta,cap_h_p(is,ir,it,ie,:),1)
             enddo
          enddo
       enddo
    enddo
    
  end subroutine GK_do
  
  subroutine get_gkRHS(ij)
    use gkcoll_globals
    use gkcoll_equilibrium
    use gkcoll_gyro
    implicit none
    integer, intent(in) :: ij
    integer :: is, ir, it, ie, ix
    integer :: id, jt
    real    :: rval
    complex :: val
    complex :: thfac
    complex :: ir_fac

    ! Get the RHS collisionless GK little_h: h = H - ze/T G phi
    rhs(ij,:,:,:,:,:) = (0.0,0.0)

    do is=1,n_species
       do ir=1,n_radial
          do it=1,n_theta
             do ie=1,n_energy
                do ix=1,n_xi
                   
                   rval = omega_stream(it,is) * sqrt(energy(ie)) * xi(ix) 
                   do id=-2,2
                      jt = thcyc(it+id)
                      if((it+id) < 1) then
                         thfac = cos(2*pi*k_theta*rmin) &
                              + i_c * sin(2*pi*k_theta*rmin)
                         if(ir-1 >= 1) then
                            ir_fac = cap_h_x(is,ir-1,jt,ie,ix)
                         else
                            ir_fac = cap_h_x(is,n_radial,jt,ie,ix)
                         endif
                      else if((it+id) > n_theta) then
                         thfac = cos(2*pi*k_theta*rmin) &
                              - i_c * sin(2*pi*k_theta*rmin)
                         if(ir+1 <= n_radial) then
                            ir_fac = cap_h_x(is,ir+1,jt,ie,ix)
                         else
                            ir_fac = cap_h_x(is,1,jt,ie,ix)
                         endif
                      else
                         thfac  = 1.0
                         ir_fac =  cap_h_x(is,ir,jt,ie,ix)
                      endif
                      rhs(ij,is,ir,it,ie,ix) = rhs(ij,is,ir,it,ie,ix) &
                           - rval * cderiv(id) * thfac * ir_fac
                    enddo

                    val = omega_trap(it,is) &
                         * sqrt(energy(ie)) * (1.0 - xi(ix)**2)
                       rhs(ij,is,ir,it,ie,ix) = rhs(ij,is,ir,it,ie,ix) &
                            - val  * cap_h_x_deriv(is,ir,it,ie,ix)

                    ! omega_rdrift
                    val = omega_rdrift(it,is) &
                         * energy(ie) * (1.0 + xi(ix)**2) &
                         * (2.0*pi*i_c*indx_r(ir)/r_length) 
                    rhs(ij,is,ir,it,ie,ix) = rhs(ij,is,ir,it,ie,ix) &
                         - val * cap_h_x(is,ir,it,ie,ix)

                    ! omega_dalpha
                    val = omega_adrift(it,is) &
                         * energy(ie) * (1.0 + xi(ix)**2) &
                         * (i_c * k_theta)
                    rhs(ij,is,ir,it,ie,ix) = rhs(ij,is,ir,it,ie,ix) &
                         - val * cap_h_x(is,ir,it,ie,ix)
                         
                    ! omega_dalpha - pressure component
                    val = omega_aprdrift(it,is) &
                         * energy(ie) * xi(ix)**2 &
                         * (i_c * k_theta)
                    rhs(ij,is,ir,it,ie,ix) = rhs(ij,is,ir,it,ie,ix) &
                         - val * cap_h_x(is,ir,it,ie,ix)

                    ! omega_dalpha - pressure component with xi derivative
                    val = omega_xprdrift(it,is) &
                         * energy(ie) * xi(ix) * (1.0 - xi(ix)**2)
                    rhs(ij,is,ir,it,ie,ix) = rhs(ij,is,ir,it,ie,ix) &
                         - val * cap_h_x(is,ir,it,ie,ix)

                    ! omega_star
                    val = i_c * k_theta &
                         * rho *sqrt(temp(is)*mass(is))/(1.0*z(is)) * vth(is) &
                         * (dlnndr(is) + dlntdr(is) * (energy(ie)-1.5))
                    rhs(ij,is,ir,it,ie,ix) = rhs(ij,is,ir,it,ie,ix) &
                            - val * z(is)/temp(is) &
                            * gyrox_J0(is,ir,it,ie,ix) * phi(ir,it)


                 end do
              end do
           end do
        end do
     end do
     
   end subroutine get_gkRHS
   
   
 end module gkcoll_gk
