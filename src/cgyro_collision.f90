module cgyro_collision
  
  implicit none
  
  public :: COLLISION_alloc, COLLISION_do
  logical, private :: initialized = .false.

  real, dimension(:,:,:,:), allocatable, private :: cmat
  real, dimension(:,:), allocatable, private :: cvec,bvec
  integer, dimension(:,:,:), allocatable, private :: indx_coll 
  integer, private :: msize

contains

  subroutine COLLISION_alloc(flag)
    use cgyro_globals
    use cgyro_gyro
    use cgyro_equilibrium, only : omega_trap, k_perp
    !use cgyro_gk
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    real, dimension(:,:,:), allocatable :: nu_d, nu_s
    real, dimension(:,:), allocatable :: rs
    real, external :: derf
    real :: xa, xb, tauinv_ab
    real :: sum_nu, sum_nud, sum_den
    integer :: is,ir,it,ix,ie,js,je,jx,kx,ks,p,pp
    ! parameters for matrix solve
    integer :: info
    integer, dimension(:), allocatable :: i_piv
    real, dimension(:), allocatable :: work 
    real, dimension(:,:), allocatable :: amat, bmat
    real, dimension(:,:), allocatable :: rs_lor
    real, dimension(:), allocatable :: vecin_xi, vecout_xi

    if(collision_model == -1) return

    if(flag == 1) then
       if(initialized) return
       
       if(collision_model == 0 .and. adiabatic_ele_model == 1) then
          call cgyro_error('ERROR: (GKCOLL) collision_model=0 requires kinetic electrons')
          return
       endif

       allocate(nu_d(n_energy,n_species,n_species))
       allocate(nu_s(n_energy,n_species,n_species))
       allocate(rs(n_species,n_species))

       nu_d(:,:,:) = 0.0
       nu_s(:,:,:) = 0.0

       do ie=1,n_energy
          do is=1,n_species
             do js=1,n_species
                xa = sqrt(energy(ie))
                xb = xa * vth(is)**2 / vth(js)**2
                tauinv_ab = nu(is) * (1.0*Z(js))**2 / (1.0*Z(is))**2 &
                     * dens(js)/dens(is)
                
                if(collision_model == 0) then
                   ! Only ee,ei Connor-like Lorentz
                   if(is == is_ele) then
                      if(is == js) then
                         ! e-e
                         nu_d(ie,is,js) = tauinv_ab * (1.0/xa**3) &
                           * (exp(-xb*xb)/(xb*sqrt(pi)) &
                           + (1.0-1.0/(2.0*xb*xb)) * DERF(xb))
                      else
                         ! e-i
                         nu_d(ie,is,js) = tauinv_ab * (1.0/xa**3)
                      endif
                   endif

                else if(collision_model == 1) then
                   ! Connor model
                   if(is == js .or. &
                        (abs(mass(is) - mass(js)) < epsilon(0.))) then
                      ! case 1: like-species/same mass collisions
                      nu_d(ie,is,js) = tauinv_ab * (1.0/xa**3) &
                           * (exp(-xb*xb)/(xb*sqrt(pi)) &
                           + (1.0-1.0/(2.0*xb*xb)) * DERF(xb))
                      
                   else if(mass(is) < mass(js)) then
                      ! case 2: ele-ion and ion-imp(heavy) collisions
                      nu_d(ie,is,js) = tauinv_ab * (1.0/xa**3)
                      
                   else
                      ! case 3: ion-ele and imp(heavy)-ion collisions
                      
                      nu_d(ie,is,js) = tauinv_ab * 4.0/(3.0*sqrt(pi)) &
                           * sqrt(mass(js)/mass(is)) &
                           * (temp(is)/temp(js))**1.5
                   endif
                   nu_s(ie,is,js) = nu_d(ie,is,js)
                   
                else
                   ! Reduced Hirshman-Sigmar model
                   nu_d(ie,is,js) = tauinv_ab * (1.0/xa**3) &
                        * (exp(-xb*xb)/(xb*sqrt(pi)) &
                        + (1.0-1.0/(2.0*xb*xb)) * DERF(xb))
                   nu_s(ie,is,js) = tauinv_ab * (1.0/xa) &
                        * (-exp(-xb*xb)/(xb*sqrt(pi)) &
                        + (1.0/(2.0*xb*xb)) * DERF(xb)) &
                        * (2.0*temp(is)/temp(js))*(1.0+mass(js)/mass(is))
                endif
                
             enddo
          enddo
       enddo

       do is=1,n_species
          do js=1,n_species
             rs(is,js) = 0.0
             do ie=1,n_energy
                rs(is,js) = rs(is,js) + w_e(ie) * nu_s(ie,is,js) &
                     * energy(ie)
             enddo
          enddo
       enddo

       msize = n_species*n_energy*n_xi
       allocate(cmat(n_radial,n_theta,msize,msize))
       allocate(cvec(msize,2))
       allocate(bvec(msize,2))
       allocate(indx_coll(n_species,n_energy,n_xi))

       p = 0
       do is=1,n_species
          do ie=1,n_energy
             do ix=1,n_xi
                p = p+1
                indx_coll(is,ie,ix) = p
             enddo
          enddo
       enddo
       
       ! set-up the collision matrix
       sum_den = 0.0
       do is=1,n_species
          sum_den = sum_den + z(is)**2 * dens(is) / temp(is)
       enddo
       if(adiabatic_ele_model == 1) then
          sum_den = sum_den + dens_ele / temp_ele
       endif
       
       ! matrix solve parameters
       allocate(work(msize))
       allocate(i_piv(msize))
       allocate(amat(msize,msize))
       allocate(bmat(msize,msize))
       
       if(collision_model == 3) then
          allocate(rs_lor(n_species,n_species))
          allocate(vecin_xi(n_xi))
          allocate(vecout_xi(n_xi))
          do ix=1,n_xi
             vecin_xi(ix) = xi(ix)
          enddo
          call DGEMV('N',n_xi,n_xi,num1,xi_lor_mat(:,:),n_xi,&
               vecin_xi,1,num0,vecout_xi,1)
          rs_lor(:,:) = 0.0
          do is=1,n_species
             do js=1,n_species
                do ie=1,n_energy
                   do ix=1,n_xi
                      rs_lor(is,js) = rs_lor(is,js) &
                           + w_e(ie) * energy(ie) &
                           * w_xi(ix) * 0.5*vecout_xi(ix) &
                           * nu_d(ie,is,js)
                   enddo
                enddo
             enddo
          enddo
          call DGEMV('T',n_xi,n_xi,num1,xi_lor_mat(:,:),n_xi,&
               vecout_xi,1,num0,vecin_xi,1)
       endif

       do ir=1,n_radial
          do it=1,n_theta
             
             cmat(ir,it,:,:) = (0.0,0.0)
             amat(:,:)       = (0.0,0.0)

             do is=1,n_species     
                do ie=1,n_energy
                   do ix=1,n_xi
                      p  = indx_coll(is,ie,ix)

                      sum_nud = 0.0
                      do ks=1,n_species
                         sum_nud = sum_nud + nu_d(ie,is,ks)
                      enddo

                      do js=1,n_species
                         do je=1,n_energy
                            do jx=1,n_xi
                               pp = indx_coll(js,je,jx)
                               
                               ! constant part
                               if(is==js .and. ie==je .and. ix==jx) then
                                  cmat(ir,it,p,pp) =  cmat(ir,it,p,pp) + 1.0
                                  amat(p,pp) = amat(p,pp) + 1.0
                               endif

                               ! Trapping term
                               if(trap_method == 1) then
                                  if(is==js .and. ie==je) then
                                     cmat(ir,it,p,pp) =  cmat(ir,it,p,pp) &
                                          + (0.5*delta_t) * omega_trap(it,is) &
                                          * sqrt(energy(ie)) &
                                          * (1.0 - xi(ix)**2) &
                                          * xi_deriv_mat(ix,jx)
                                     amat(p,pp) =  amat(p,pp) &
                                          - (0.5*delta_t) * omega_trap(it,is) &
                                          * sqrt(energy(ie)) &
                                          * (1.0 - xi(ix)**2) &
                                          * xi_deriv_mat(ix,jx)
                                  endif
                               endif

                               ! Collision component: Lorentz
                               if(is==js .and. ie==je) then
                                 cmat(ir,it,p,pp)  &
                                       =  cmat(ir,it,p,pp) &
                                       - (0.5*delta_t) * xi_lor_mat(ix,jx) &
                                       *0.5*sum_nud
                                  
                                  amat(p,pp) &
                                       = amat(p,pp) &
                                       + (0.5*delta_t) * xi_lor_mat(ix,jx) &
                                       *0.5*sum_nud
                               endif

                               ! Collision component: Restoring
                               if(collision_model == 3) then
                                  ! EAB: NEED to recheck 0.5 with Lorentz
                                  cmat(ir,it,p,pp)  &
                                       =  cmat(ir,it,p,pp) &
                                       - (0.5*delta_t) &
                                       * (-mass(js)/mass(is)) &
                                       * (dens(js)/dens(is)) &
                                       * (vth(js)/vth(is)) &
                                       * nu_d(ie,is,js) * sqrt(energy(ie)) &
                                       * 0.5*vecout_xi(ix) &
                                       * nu_d(je,js,is) * sqrt(energy(je)) &
                                       * 0.5*vecin_xi(jx) * w_e(je)*w_xi(jx) &
                                       / rs_lor(is,js)
                                  amat(p,pp)  &
                                       =  amat(p,pp) &
                                       + (0.5*delta_t) &
                                       * (-mass(js)/mass(is)) &
                                       * (dens(js)/dens(is)) &
                                       * (vth(js)/vth(is)) &
                                       * nu_d(ie,is,js) * sqrt(energy(ie)) &
                                       * 0.5*vecout_xi(ix) &
                                       * nu_d(je,js,is) * sqrt(energy(je)) &
                                       * 0.5*vecin_xi(jx) * w_e(je)*w_xi(jx) &
                                       / rs_lor(is,js)

                               else if(abs(rs(is,js)) > epsilon(0.)) then
                                  cmat(ir,it,p,pp)  &
                                       =  cmat(ir,it,p,pp) &
                                       - (0.5*delta_t) &
                                       * 1.5 * (mass(js)/mass(is)) &
                                       * (dens(js)/dens(is)) &
                                       * (vth(js)/vth(is)) &
                                       * nu_s(ie,is,js) * sqrt(energy(ie)) &
                                       * xi(ix) &
                                       * nu_s(je,js,is) * sqrt(energy(je)) &
                                       * xi(jx) * w_e(je) * w_xi(jx) &
                                       / rs(is,js)
                                  amat(p,pp)  &
                                       =  amat(p,pp) &
                                       + (0.5*delta_t) &
                                       * 1.5 * (mass(js)/mass(is)) &
                                       * (dens(js)/dens(is)) &
                                       * (vth(js)/vth(is)) &
                                       * nu_s(ie,is,js) * sqrt(energy(ie)) &
                                       * xi(ix) &
                                       * nu_s(je,js,is) * sqrt(energy(je)) &
                                       * xi(jx) * w_e(je) * w_xi(jx) &
                                       / rs(is,js)
                               endif

                               if(collision_model == 2) then
                                  if(is==js .and. ie==je) then
                                     sum_nu = 0.0
                                     do ks=1,n_species
                                        sum_nu = sum_nu &
                                             + (nu_d(ie,is,ks) &
                                             -nu_s(ie,is,ks))
                                     enddo
                                     cmat(ir,it,p,pp)  &
                                          =  cmat(ir,it,p,pp) &
                                          - (0.5*delta_t) * sum_nu &
                                          * 1.5 * xi(ix) * xi(jx) * w_xi(jx)
                                     amat(p,pp)  &
                                          =  amat(p,pp) &
                                          + (0.5*delta_t) * sum_nu &
                                          * 1.5 * xi(ix) * xi(jx) * w_xi(jx)
                                  endif
                               endif

                               ! Poisson component 
                               cmat(ir,it,p,pp)  &
                                    =  cmat(ir,it,p,pp) &
                                    - z(is)/temp(is) / &
                                    (-k_perp(it,ir)**2 * lambda_debye**2 &
                                    * dens_ele / temp_ele &
                                     + sum_den) &
                                    * gyrox_J0(is,ir,it,ie,ix) &
                                    * z(js)*dens(js) &
                                    * gyrox_J0(js,ir,it,je,jx) * w_e(je) &
                                    * 0.5 * w_xi(jx)
                               amat(p,pp)  &
                                    =  amat(p,pp) &
                                    - z(is)/temp(is) / &
                                    (-k_perp(it,ir)**2 * lambda_debye**2 &
                                    * dens_ele / temp_ele &
                                     + sum_den) &
                                    * gyrox_J0(is,ir,it,ie,ix) &
                                    * z(js)*dens(js) &
                                    * gyrox_J0(js,ir,it,je,jx) * w_e(je) &
                                    * 0.5 * w_xi(jx)
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             enddo

             ! H_bar = (1 - dt/2 C - Poisson)^(-1) * (1 + dt/2 C + Poisson) H
             ! Lapack factorization and inverse of LHS
             call DGETRF(msize,msize,cmat(ir,it,:,:),msize,i_piv,info)
             call DGETRI(msize,cmat(ir,it,:,:),msize,i_piv,work,msize,info)
             ! Matrix multiply
             call DGEMM('N','N',msize,msize,msize,num1,cmat(ir,it,:,:),&
                  msize,amat,msize,num0,bmat,msize)
             cmat(ir,it,:,:) = bmat(:,:)
      
          enddo
       enddo

       deallocate(amat)
       deallocate(bmat)
       deallocate(work)
       deallocate(i_piv)
       deallocate(nu_d)
       deallocate(nu_s)
       deallocate(rs)
       if(collision_model == 3) then
          deallocate(rs_lor)
          deallocate(vecin_xi)
          deallocate(vecout_xi)
       endif

       initialized = .true.

    else
       if(.NOT. initialized) return
       deallocate(cmat)
       deallocate(cvec)
       deallocate(bvec)
       deallocate(indx_coll)
       initialized = .false.
    endif
    
  end subroutine COLLISION_alloc
  
  subroutine COLLISION_do

    use timer_lib

    use cgyro_globals
    use cgyro_poisson
    use cgyro_gyro

    implicit none

    integer :: is,ir,it,ie,ix
    integer :: p

    if (collision_model == -1) return

    ! compute new collisional cap_H: H = h + ze/T G phi
    ! assumes have cap_h_x

   call timer_lib_in('collision')
  
    do ir=1,n_radial
       do it=1,n_theta
          
          ! Set-up the RHS: H = f + ze/T G phi
          ! 
          ! Pack Re(H) and Im(H) into two-column matrix for use with DGEMM
          do is=1,n_species
             do ie=1,n_energy
                do ix=1,n_xi
                   p = indx_coll(is,ie,ix)
                   cvec(p,1) = real(cap_h_x(is,ir,it,ie,ix))
                   cvec(p,2) = imag(cap_h_x(is,ir,it,ie,ix))
                enddo
             enddo
          enddo
 
          
          ! Solve for H
           call DGEMM('N','N',msize,2,msize,num1,cmat(ir,it,:,:),&
               msize,cvec,msize,num0,bvec,msize)
          
          do is=1,n_species
             do ie=1,n_energy
                do ix=1,n_xi
                   p = indx_coll(is,ie,ix)
                   cap_h_x(is,ir,it,ie,ix) = bvec(p,1) + i_c * bvec(p,2)
                enddo
             enddo
          enddo

          ! Compute the new phi
          call POISSONh_do
          
          ! Compute the new h_x
          do is=1,n_species
             do ie=1,n_energy
                do ix=1,n_xi
                   h_x(is,ir,it,ie,ix) = cap_h_x(is,ir,it,ie,ix) &
                        - z(is)/temp(is) * gyrox_J0(is,ir,it,ie,ix) &
                        * phi(ir,it)
                enddo
             enddo
          enddo

       enddo
    enddo

  call timer_lib_out('collision')
   
  end subroutine COLLISION_do
  
end module cgyro_collision
