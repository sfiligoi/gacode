module gkcoll_implicit
  implicit none

  public :: GKimp_alloc, GKimp_do

  logical, private :: initialized = .false.

  complex, dimension(:,:,:,:,:), allocatable, private :: h0_x
  complex, dimension(:,:,:,:,:), allocatable, private :: gemat, gimat
  complex, dimension(:,:), allocatable, private :: pimat
  integer, dimension(:,:), allocatable, private :: indx_gmat
  integer,private :: msize
  ! parameters for matrix solve
  integer, private :: info
  integer, dimension(:), allocatable, private :: i_piv
  complex, dimension(:), allocatable, private :: work 

contains

  subroutine GKimp_alloc(flag)
    use gkcoll_globals
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    integer :: ir,it,jr,jt,p
    
    if(flag == 1) then
       if(initialized) return

       allocate(h0_x(n_species,n_radial,n_theta,n_energy,n_xi))

       msize = n_radial * n_theta
       allocate(gemat(n_species,n_energy,n_xi,msize,msize))
       allocate(gimat(n_species,n_energy,n_xi,msize,msize))
       allocate(pimat(msize,msize))
       allocate(indx_gmat(n_radial,n_theta))
       allocate(work(msize))
       allocate(i_piv(msize))
       p = 0
       do ir=1,n_radial
          do it=1,n_theta
             p = p + 1
             indx_gmat(ir,it) = p
          enddo
       enddo
       call GKimp_matinit

       initialized = .true.

    else
       if(.NOT. initialized) return

       deallocate(h0_x)
       deallocate(gemat)
       deallocate(gimat)
       deallocate(pimat)
       deallocate(indx_gmat)
       deallocate(work)
       deallocate(i_piv)

       initialized = .false.
    endif

  end subroutine GKimp_alloc

  subroutine GKimp_matinit
    use gkcoll_globals
    use gkcoll_gyro
    use gkcoll_equilibrium
    implicit none
    integer :: p, pp, ir, jr, it, jt, id, is, ie, ix
    complex :: thfac
    real :: rval
    real :: sum_den
    complex, dimension(:,:), allocatable :: amat
    complex :: alpha = (1.0,0.0)
    complex :: beta  = (0.0,0.0)
    real, dimension(-2:2) :: cderiv
    integer, dimension(:), allocatable :: thcyc

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

    allocate(amat(msize,msize))

    gemat(:,:,:,:,:) = (0.0,0.0)
    gimat(:,:,:,:,:) = (0.0,0.0)

    do is=1,n_species
       do ie=1,n_energy
          do ix=1,n_xi

             amat(:,:) = (0.0,0.0)
             do p=1,msize
                gemat(is,ie,ix,p,p) = gemat(is,ie,ix,p,p) + 1.0
             enddo

             do ir=1,n_radial
                do it=1,n_theta
                   p = indx_gmat(ir,it)
                   rval = omega_stream(it,is) * sqrt(energy(ie)) * xi(ix)

                   do id=-2,2
                      jt = thcyc(it+id)
                      if((it+id) < 1) then
                         thfac = cos(2*pi*k_theta*rmin) &
                              + i_c * sin(2*pi*k_theta*rmin)
                         if(ir-1 >= 1) then
                            jr = ir - 1
                         else
                            jr = n_radial
                         endif
                      else if((it+id) > n_theta) then
                         thfac = cos(2*pi*k_theta*rmin) &
                              - i_c * sin(2*pi*k_theta*rmin)
                         if(ir+1 <= n_radial) then
                            jr = ir+1
                         else
                            jr = 1
                         endif
                      else
                         thfac = 1.0
                         jr = ir
                      endif
                      pp = indx_gmat(jr,jt)
                      gemat(is,ie,ix,p,pp) = gemat(is,ie,ix,p,pp) &
                           + delta_t &
                           * rval * cderiv(id) * thfac
                      amat(p,pp) = amat(p,pp) - delta_t &
                           * rval * cderiv(id) * thfac
                      
                   enddo
                enddo
             enddo

             ! Lapack factorization and inverse
             call ZGETRF(msize,msize,gemat(is,ie,ix,:,:),msize,i_piv,info)
             call ZGETRI(msize,gemat(is,ie,ix,:,:),msize,i_piv,work,msize,info)
             ! Matrix multiply
             call ZGEMM('N','N',msize,msize,msize,alpha,gemat(is,ie,ix,:,:),&
                  msize,amat,msize,beta,gimat(is,ie,ix,:,:),msize)

          enddo
       enddo
    enddo


    deallocate(amat)
    deallocate(thcyc)

    pimat(:,:) = (0.0,0.0)
    do ir=1,n_radial
       do it=1,n_theta
          p = indx_gmat(ir,it)
          
          sum_den = 0.0
          do is=1,n_species
             do ie=1,n_energy
                do ix=1,n_xi
                   sum_den = sum_den &
                        + 0.5 * w_xi(ix) &
                        * (1.0 - gyrox_J0(is,ir,it,ie,ix)**2) &
                        * z(is)**2/temp(is) *dens(is) * w_e(ie)
                enddo
             enddo
          enddo
          if(adiabatic_ele_model == 1) then
             sum_den = sum_den + dens_ele / temp_ele
          endif
          sum_den = sum_den - k_perp(it,ir)**2 * lambda_debye**2 &
               * dens_ele / temp_ele

          pimat(p,p) = sum_den
       enddo
    enddo

    do ir=1,n_radial
       do it=1,n_theta
          p = indx_gmat(ir,it)
          do jr=1,n_radial
             do jt=1,n_theta
                pp = indx_gmat(jr,jt)
                do is=1,n_species
                   do ie=1,n_energy
                      do ix=1,n_xi
                         pimat(p,pp) = pimat(p,pp) - 0.5 * w_xi(ix) * w_e(ie) &
                              * z(is)**2 * dens(is) / temp(is) &
                              * gyrox_J0(is,ir,it,ie,ix) &
                              * gyrox_J0(is,jr,jt,ie,ix) &
                              * gimat(is,ie,ix,p,pp)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    ! Lapack factorization and inverse
    call ZGETRF(msize,msize,pimat,msize,i_piv,info)
    call ZGETRI(msize,pimat,msize,i_piv,work,msize,info)

  end subroutine GKimp_matinit

  subroutine GKimp_do
    use gkcoll_globals
    use gkcoll_gyro
    use gkcoll_equilibrium
    use gkcoll_gk, only : xi_mat, xi_mat_inv
    implicit none
    integer :: is,ir,it,ie,ix, jr, jt, p, pp
    complex :: val
    complex :: alpha = (1.0,0.0)
    complex :: beta  = (0.0,0.0)

    ! compute new collisionless little_h: h = H - ze/T G phi

    h0_x = h_x
    
    do is=1,n_species
       do ir=1,n_radial
          do it=1,n_theta
             do ie=1,n_energy
                call ZGEMV('N',n_xi,n_xi,alpha,xi_mat,n_xi,&
                     cap_h_p(is,ir,it,ie,:),1,beta,cap_h_x(is,ir,it,ie,:),1)
             enddo
          enddo
       enddo
    enddo

    ! explicit terms
    h_x(:,:,:,:,:) = (0.0,0.0)
    do is=1,n_species
       do ir=1,n_radial
          do it=1,n_theta
             do ie=1,n_energy
                do ix=1,n_xi

                   do jr=1,n_radial
                      do jt=1,n_theta
                         p  = indx_gmat(ir,it)
                         pp = indx_gmat(jr,jt)
                         
                         ! time derivative term
                         h_x(is,ir,it,ie,ix) = h_x(is,ir,it,ie,ix) &
                              + h0_x(is,jr,jt,ie,ix) * gemat(is,ie,ix,p,pp)

                         ! omega_rdrift
                         val = omega_rdrift(jt,is) &
                              * energy(ie) * (1.0 + xi(ix)**2) &
                              * (2.0*pi*i_c*indx_r(jr)*r_length_inv) 
                         h_x(is,ir,it,ie,ix) = h_x(is,ir,it,ie,ix) &
                              - val * cap_h_x(is,jr,jt,ie,ix) &
                              * delta_t * gemat(is,ie,ix,p,pp)
                         ! radial upwinding
                         val = omega_rdrift(jt,is) &
                              * energy(ie) * (1.0 + xi(ix)**2)
                         if(abs(val) /= 0.0) then
                            val = abs(val)/val
                         endif
                         val = val * rupwind_eps & 
                              * (1.0*indx_r(jr)/(1.0*n_radial))**rupwind_n &
                              * (2.0*pi*r_length_inv) 
                         h_x(is,ir,it,ie,ix) = h_x(is,ir,it,ie,ix) &
                              - val * h0_x(is,jr,jt,ie,ix)&
                              * delta_t * gemat(is,ie,ix,p,pp)

                         ! omega_dalpha
                         val = omega_adrift(jt,is) &
                              * energy(ie) * (1.0 + xi(ix)**2) &
                              * (i_c * k_theta)
                         h_x(is,ir,it,ie,ix) = h_x(is,ir,it,ie,ix) &
                              - val * cap_h_x(is,jr,jt,ie,ix)&
                              * delta_t * gemat(is,ie,ix,p,pp)
                         
                         ! omega_dalpha - pressure component
                         val = omega_aprdrift(jt,is) &
                              * energy(ie) * xi(ix)**2 &
                              * (i_c * k_theta)
                         h_x(is,ir,it,ie,ix) = h_x(is,ir,it,ie,ix) &
                              - val * cap_h_x(is,jr,jt,ie,ix)&
                              * delta_t * gemat(is,ie,ix,p,pp)

                         ! omega_star
                         val = i_c * k_theta &
                              * rho *sqrt(temp(is)*mass(is))/(1.0*z(is)) &
                              * vth(is) &
                              * (dlnndr(is) + dlntdr(is) * (energy(ie)-1.5))
                         h_x(is,ir,it,ie,ix) = h_x(is,ir,it,ie,ix) &
                              - val * z(is)/temp(is) &
                              * gyrox_J0(is,jr,jt,ie,ix) * phi(jr,jt)&
                              * delta_t * gemat(is,ie,ix,p,pp)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    
    ! compute the new phi
    call GKimp_POISSONx_do

    ! implicit streaming term
    do is=1,n_species
       do ir=1,n_radial
          do it=1,n_theta
             do ie=1,n_energy
                do ix=1,n_xi

                   do jr=1,n_radial
                      do jt=1,n_theta
                         p  = indx_gmat(ir,it)
                         pp = indx_gmat(jr,jt)
                         
                         h_x(is,ir,it,ie,ix) = h_x(is,ir,it,ie,ix) &
                              + gimat(is,ie,ix,p,pp) * z(is)/temp(is) &
                              * gyrox_J0(is,jr,jt,ie,ix) * phi(jr,jt)
                      enddo
                   enddo
                enddo
             enddo
          enddo
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

  end subroutine GKimp_do

  subroutine GKimp_POISSONx_do
    use gkcoll_globals
    use gkcoll_gyro
    implicit none
    integer :: is, ie, ix, ir, it, jr, jt, p, pp

    phi = (0.0,0.0)
    do ir=1,n_radial
       do it=1,n_theta
          p  = indx_gmat(ir,it)
          do jr=1,n_radial
             do jt=1,n_theta
                pp = indx_gmat(jr,jt)
                do is=1,n_species
                   do ie=1,n_energy
                      do ix=1,n_xi
                         phi(ir,it) = phi(ir,it) &
                              + 0.5 * w_xi(ix) &
                              * gyrox_J0(is,jr,jt,ie,ix) &
                              * z(is)*dens(is) * w_e(ie) &
                              * h_x(is,jr,jt,ie,ix) * pimat(p,pp)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

  end subroutine GKimp_POISSONx_do

end module gkcoll_implicit
