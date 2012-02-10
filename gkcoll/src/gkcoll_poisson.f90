module gkcoll_poisson
  
  implicit none
  
  public :: POISSON_alloc, POISSONp_do, POISSONx_do
  logical, private :: initialized = .false.
  complex, private :: sum_den_p, sum_den_x
  complex, dimension(:,:,:), allocatable, private :: pzf ! for n=0 test
  complex, dimension(:,:,:), allocatable, private :: dzf
  complex, dimension(:), allocatable, private :: ptemp
  
contains
  
  subroutine POISSON_alloc(flag)
    use gkcoll_globals
    use gkcoll_gyro
    use gkcoll_equilibrium
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    integer :: is, ie, ix, ir, it, jt
    integer :: info
    integer, dimension(:), allocatable :: i_piv
    complex, dimension(:), allocatable :: work
    
    if(flag == 1) then
       if(initialized) return
       
       sum_den_p = (0.0,0.0)
       sum_den_x = (0.0,0.0)
       do is=1,n_species
          do ie=1,n_energy
             do ix=1,n_xi
                sum_den_p = sum_den_p &
                     + 0.5 * w_xi(ix) &
                     * z(is)**2/temp(is) *dens(is) * w_e(ie)
                sum_den_x = sum_den_x &
                     + 0.5 * w_xi(ix) &
                     * (1.0 - gyrox_J0(is,ir,it,ie,ix)**2) &
                     * z(is)**2/temp(is) *dens(is) * w_e(ie)
             enddo
          enddo
       enddo
       if(adiabatic_ele_model == 1) then
          sum_den_p = sum_den_p + dens_ele / temp_ele
          sum_den_x = sum_den_x + dens_ele / temp_ele
       endif
       
       if(toroidal_model == 2 .and. adiabatic_ele_model == 1 &
            .and. neoclassical_model /= 1) then
          
          allocate(pzf(n_radial,n_theta,n_theta))
          pzf(:,:,:) = (0.0,0.0)      
          do ir=1,n_radial
             do it=1,n_theta
                pzf(ir,it,it) = -k_perp(it,ir)**2 * lambda_debye**2 &
                     * dens_ele / temp_ele + sum_den_p
                do jt=1,n_theta
                   pzf(ir,it,jt) = pzf(ir,it,jt) &
                        - dens_ele / temp_ele * w_theta(jt)
                enddo
             enddo
          enddo
          allocate(work(n_theta))
          allocate(i_piv(n_theta))
          do ir=1,n_radial
             call ZGETRF(n_theta,n_theta,pzf(ir,:,:),n_theta,i_piv,info)
             call ZGETRI(n_theta,pzf(ir,:,:),n_theta,i_piv,work,n_theta,info)
          enddo
          deallocate(i_piv)
          deallocate(work)
          
          allocate(dzf(n_radial,n_theta,n_theta))
          dzf(:,:,:) = (0.0,0.0)      
          do ir=1,n_radial
             do it=1,n_theta
                dzf(ir,it,it) = -k_perp(it,ir)**2 * lambda_debye**2 &
                     * dens_ele / temp_ele + sum_den_x
                do jt=1,n_theta
                   dzf(ir,it,jt) = dzf(ir,it,jt) &
                        - dens_ele / temp_ele * w_theta(jt)
                enddo
             enddo
          enddo
          allocate(work(n_theta))
          allocate(i_piv(n_theta))
          do ir=1,n_radial
             call ZGETRF(n_theta,n_theta,dzf(ir,:,:),n_theta,i_piv,info)
             call ZGETRI(n_theta,dzf(ir,:,:),n_theta,i_piv,work,n_theta,info)
          enddo
          deallocate(i_piv)
          deallocate(work)
          
          allocate(ptemp(n_theta))
       endif

       initialized = .true.
       
    else
       if(.NOT. initialized) return

       if(toroidal_model == 2 .and. adiabatic_ele_model == 1 &
            .and. neoclassical_model /= 1) then
          deallocate(pzf)
          deallocate(dzf)
          deallocate(ptemp)
       endif

       initialized = .false.
    endif

  end subroutine POISSON_alloc

  subroutine POISSONp_do
    use gkcoll_globals
    use gkcoll_gyro
    use gkcoll_equilibrium
    implicit none
    integer :: is, ie, ix, ir, it
    complex :: alpha = (1.0,0.0)
    complex :: beta  = (0.0,0.0)

    do ir=1,n_radial
       do it=1,n_theta
          phi(ir,it) = (0.0,0.0)
          do is=1,n_species
             do ie=1,n_energy
                do ix=1,n_xi
                   phi(ir,it) = phi(ir,it) &
                        + gyrop_J0(is,ir,it,ie,ix) &
                        * z(is)*dens(is) * w_e(ie) &
                        * cap_h_p(is,ir,it,ie,ix)
                enddo
             enddo
          enddo
       enddo

       if(toroidal_model == 2 .and. adiabatic_ele_model == 1 &
            .and. neoclassical_model /= 1) then
          call ZGEMV('N',n_theta,n_theta,alpha,pzf(ir,:,:),&
               n_theta,phi(ir,:),1,beta,ptemp,1)
          phi(ir,:) = ptemp(:)
       else
          do it=1,n_theta
             phi(ir,it) = phi(ir,it) / (-k_perp(it,ir)**2 * lambda_debye**2 &
                  * dens_ele / temp_ele + sum_den_p)
          enddo
       endif
    enddo
    
  end subroutine POISSONp_do
  
  subroutine POISSONx_do
    use gkcoll_globals
    use gkcoll_gyro
    use gkcoll_equilibrium
    implicit none
    integer :: is, ie, ix, ir, it
    complex :: alpha = (1.0,0.0)
    complex :: beta  = (0.0,0.0)
    
    do ir=1,n_radial
       do it=1,n_theta
          phi(ir,it) = (0.0,0.0)
          do is=1,n_species
             do ie=1,n_energy
                do ix=1,n_xi
                   phi(ir,it) = phi(ir,it) &
                        + 0.5 * w_xi(ix) &
                        * gyrox_J0(is,ir,it,ie,ix) &
                        * z(is)*dens(is) * w_e(ie) &
                        * h_x(is,ir,it,ie,ix)
                enddo
             enddo
          enddo
       enddo

       if(toroidal_model == 2 .and. adiabatic_ele_model == 1 &
            .and. neoclassical_model /= 1) then
          call ZGEMV('N',n_theta,n_theta,alpha,dzf(ir,:,:),&
               n_theta,phi(ir,:),1,beta,ptemp,1)
          phi(ir,:) = ptemp(:)
       else
          do it=1,n_theta
             phi(ir,it) = phi(ir,it) / (-k_perp(it,ir)**2 * lambda_debye**2 &
                  * dens_ele / temp_ele + sum_den_x)
          enddo
       endif
    enddo

  end subroutine POISSONx_do

end module gkcoll_poisson
