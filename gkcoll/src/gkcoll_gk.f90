module gkcoll_gk
  implicit none

  public :: GK_alloc, GK_do

  logical, private :: initialized = .false.

  complex, dimension(:,:,:,:,:,:), allocatable, private :: rhs
  complex, dimension(:,:,:,:,:), allocatable, private :: h0_x
  complex, dimension(:,:,:,:,:), allocatable, private :: cap_h_x_deriv

  ! theta derivative variables
  integer, dimension(-2:2) :: cderiv
  integer, dimension(:), allocatable :: thcyc

  ! xi conversion matrices
  real, dimension(:,:), allocatable :: xi_mat, xi_mat_inv, xi_deriv_mat

contains

  subroutine GK_alloc(flag)
    use gkcoll_globals
    use gkcoll_legendre
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    integer :: it, ix, jx
    ! LAPACK
    integer :: info
    integer, dimension(:), allocatable :: i_piv
    real, dimension(:), allocatable :: work

    if(flag == 1) then
       if(initialized) return

       allocate(rhs(4,n_species,n_radial,n_theta,n_energy,n_xi))
       allocate(h0_x(n_species,n_radial,n_theta,n_energy,n_xi))
       allocate(cap_h_x_deriv(n_species,n_radial,n_theta,n_energy,n_xi))

       ! cyclic index (for theta-periodicity)
       ! EAB: NOTE: NOTE CORRECT
       allocate(thcyc(1-n_theta:2*n_theta))
       do it=1,n_theta
          thcyc(it-n_theta) = it
          thcyc(it) = it
          thcyc(it+n_theta) = it
       enddo
       ! coefficients for 4th order centered derivative
       cderiv(-2) =  1
       cderiv(-1) = -8
       cderiv(0)  =  0
       cderiv(1)  =  8
       cderiv(2)  = -1

       ! xi grid to legendre conversion matrices
       allocate(xi_mat(n_xi,n_xi))
       allocate(xi_mat_inv(n_xi,n_xi))
       allocate(xi_deriv_mat(n_xi,n_xi))
       do ix=1,n_xi
          do jx=1,n_xi
             call legendre(indx_xi(jx),xi(ix),xi_mat(ix,jx))
             call legendre_deriv(indx_xi(jx),xi(ix),xi_deriv_mat(ix,jx))
          enddo
       enddo
       xi_mat_inv = xi_mat
       allocate(work(n_xi))
       allocate(i_piv(n_xi))
       call DGETRF(n_xi,n_xi,xi_mat_inv,n_xi,i_piv,info)
       call DGETRI(n_xi,xi_mat_inv,n_xi,i_piv,work,n_xi,info)
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
                call DGEMV('N',n_xi,n_xi,1.0,xi_mat,n_xi,&
                     cap_h_p(is,ir,it,ie,:),1,0.0,cap_h_x(is,ir,it,ie,:),1)
                call DGEMV('N',n_xi,n_xi,1.0,xi_deriv_mat,n_xi,&
                     cap_h_p(is,ir,it,ie,:),1,0.0,&
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
                call DGEMV('N',n_xi,n_xi,1.0,xi_mat_inv,n_xi,&
                     cap_h_x(is,ir,it,ie,:),1,0.0,cap_h_p(is,ir,it,ie,:),1)
                call DGEMV('N',n_xi,n_xi,1.0,xi_deriv_mat,n_xi,&
                     cap_h_p(is,ir,it,ie,:),1,0.0,&
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
                call DGEMV('N',n_xi,n_xi,1.0,xi_mat_inv,n_xi,&
                     cap_h_x(is,ir,it,ie,:),1,0.0,cap_h_p(is,ir,it,ie,:),1)
                call DGEMV('N',n_xi,n_xi,1.0,xi_deriv_mat,n_xi,&
                     cap_h_p(is,ir,it,ie,:),1,0.0,&
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
                call DGEMV('N',n_xi,n_xi,1.0,xi_mat_inv,n_xi,&
                     cap_h_x(is,ir,it,ie,:),1,0.0,cap_h_p(is,ir,it,ie,:),1)
                call DGEMV('N',n_xi,n_xi,1.0,xi_deriv_mat,n_xi,&
                     cap_h_p(is,ir,it,ie,:),1,0.0,&
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
                call DGEMV('N',n_xi,n_xi,1.0,xi_mat_inv,n_xi,&
                     cap_h_x(is,ir,it,ie,:),1,0.0,cap_h_p(is,ir,it,ie,:),1)
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
    real :: val
    
    ! Get the RHS collisionless GK little_h: h = H - ze/T G phi
    rhs(ij,:,:,:,:,:) = (0.0,0.0)

    do is=1,n_species
       do ir=1,n_radial
          do it=1,n_theta
             do ie=1,n_energy
                do ix=1,n_xi
                   
                   ! EAB note: need to fix boundary conditions
                   val = omega_stream(it,is) * sqrt(energy(ie)) * xi(ix) 
                   do id=-2,2
                      jt = thcyc(it+id)
                      rhs(ij,is,ir,it,ie,ix) = rhs(ij,is,ir,it,ie,ix) &
                           - val * cderiv(id) * cap_h_x(is,ir,jt,ie,ix)
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
