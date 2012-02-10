module gkcoll_collision
  
  implicit none
  
  public :: COLLISION_alloc, COLLISION_do
  logical, private :: initialized = .false.

  complex, dimension(:,:,:,:), allocatable, private :: cmat
  complex, dimension(:), allocatable, private :: cvec, bvec
  integer, dimension(:,:,:), allocatable, private :: indx_coll 
  integer, private :: msize

contains

  subroutine COLLISION_alloc(flag)
    use gkcoll_globals
    use gkcoll_gyro
    use gkcoll_equilibrium, only : omega_trap, k_perp
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    real, dimension(:,:,:), allocatable :: nu_d, nu_s
    real, dimension(:,:), allocatable :: rs
    real, external :: derf
    real :: xa, xb, tauinv_ab
    real :: sum_nu, sum_den
    integer :: is,ir,it,ix,ie,js,je,jx, ks,p, pp
    ! parameters for matrix solve
    integer :: info
    integer, dimension(:), allocatable :: i_piv
    complex, dimension(:), allocatable :: work 
    complex, dimension(:,:), allocatable :: amat, bmat
    complex :: alpha = (1.0,0.0)
    complex :: beta  = (0.0,0.0) 

    if(collision_model == -1) return
    if(imp_flag==2) return

    if(flag == 1) then
       if(initialized) return
       
       allocate(nu_d(n_energy,n_species,n_species))
       allocate(nu_s(n_energy,n_species,n_species))
       allocate(rs(n_species,n_species))

       do ie=1,n_energy
          do is=1,n_species
             do js=1,n_species
                xa = sqrt(energy(ie))
                xb = xa * vth(is)**2 / vth(js)**2
                tauinv_ab = nu(is) * (1.0*Z(js))**2 / (1.0*Z(is))**2 &
                     * dens(js)/dens(is)
                
                if(collision_model == 1) then
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
                rs(is,js) = rs(is,js) + w_e(ie) * nu_s(ie,is,js) * energy(ie)
             enddo
          enddo
       enddo

       msize = n_species*n_energy*n_xi
       allocate(cmat(n_radial,n_theta,msize,msize))
       allocate(cvec(msize))
       allocate(bvec(msize))
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

       do ir=1,n_radial
          do it=1,n_theta
             
             cmat(ir,it,:,:) = (0.0,0.0)
             amat(:,:)       = (0.0,0.0)

             do is=1,n_species     
                do ie=1,n_energy
                   do ix=1,n_xi
                      p  = indx_coll(is,ie,ix)
                      do js=1,n_species
                         do je=1,n_energy
                            do jx=1,n_xi
                               pp = indx_coll(js,je,jx)
                               
                               if(is==js .and. ie==je .and. ix==jx) then
                                  ! Collision component: Lorentz
                                  
                                  sum_nu = 0.0
                                  do ks=1,n_species
                                     sum_nu = sum_nu + nu_d(ie,is,ks)
                                  enddo

                                  cmat(ir,it,p,pp)  &
                                       =  cmat(ir,it,p,pp) &
                                       + 1.0 &
                                       - (0.5*delta_t) &
                                       * 0.5 * (-1.0*indx_xi(ix)) &
                                       * (indx_xi(ix)+1.0) * sum_nu
                                  
                                  amat(p,pp) &
                                       = amat(p,pp) + 1.0 &
                                       + (0.5*delta_t) * 0.5 &
                                       * (-1.0*indx_xi(ix)) &
                                       * (indx_xi(ix)+1.0) * sum_nu
                               endif

                               ! Collision component: Restoring
                               if(ix==1 .and. jx==1) then 
                                  if(abs(rs(is,js)) > epsilon(0.)) then
                                     cmat(ir,it,p,pp)  &
                                          =  cmat(ir,it,p,pp) &
                                          - (0.5*delta_t) &
                                          * (mass(js)/mass(is)) &
                                          *(dens(js)/dens(is)) &
                                          * (vth(js)/vth(is)) &
                                          * nu_s(ie,is,js) * sqrt(energy(ie)) &
                                          * nu_s(je,js,is) * sqrt(energy(je)) &
                                          * w_e(je)/ rs(is,js)
                                     amat(p,pp)  &
                                          =  amat(p,pp) &
                                          + (0.5*delta_t) &
                                          * (mass(js)/mass(is)) &
                                          *(dens(js)/dens(is)) &
                                          * (vth(js)/vth(is)) &
                                          * nu_s(ie,is,js) * sqrt(energy(ie)) &
                                          * nu_s(je,js,is) * sqrt(energy(je)) &
                                          * w_e(je)/ rs(is,js)
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
                                          - (0.5*delta_t) * sum_nu
                                        amat(p,pp)  &
                                          =  amat(p,pp) &
                                          + (0.5*delta_t) * sum_nu
                                     endif
                                  endif
                               endif

                               ! Poisson component 
                               cmat(ir,it,p,pp)  &
                                    =  cmat(ir,it,p,pp) &
                                    - z(is)/temp(is) / &
                                    (-k_perp(it,ir)**2 * lambda_debye**2 &
                                    * dens_ele / temp_ele &
                                     + sum_den) &
                                    * (2.0*indx_xi(ix) + 1) &
                                    * gyrop_J0(is,ir,it,ie,ix) &
                                    * z(js)*dens(js) &
                                    * gyrop_J0(js,ir,it,je,jx) * w_e(je)
                               amat(p,pp)  &
                                    =  amat(p,pp) &
                                    - z(is)/temp(is) / &
                                    (-k_perp(it,ir)**2 * lambda_debye**2 &
                                    * dens_ele / temp_ele &
                                     + sum_den) &
                                    * (2.0*indx_xi(ix) + 1) &
                                    * gyrop_J0(is,ir,it,ie,ix) &
                                    * z(js)*dens(js) &
                                    * gyrop_J0(js,ir,it,je,jx) * w_e(je)

                               ! Trapping component
                               if(trap_method /= 0) then
                                  if(is==js .and. ie==je) then
                                     if(indx_xi(jx) == indx_xi(ix+1)) then
                                        cmat(ir,it,p,pp)  &
                                             =  cmat(ir,it,p,pp) & 
                                             + (0.5*delta_t) &
                                             * omega_trap(it,is) &
                                             * sqrt(energy(ie)) &
                                             * (indx_xi(ix)+1.0) &
                                             * (indx_xi(ix)+2.0) &
                                             / (2*indx_xi(ix)+3.0)
                                        amat(p,pp)  &
                                             =  amat(p,pp) & 
                                             - (0.5*delta_t) &
                                             * omega_trap(it,is) &
                                             * sqrt(energy(ie)) &
                                             * (indx_xi(ix)+1.0) &
                                             * (indx_xi(ix)+2.0) &
                                             / (2*indx_xi(ix)+3.0)
                                     else if(indx_xi(jx) == indx_xi(ix-1)) then
                                        cmat(ir,it,p,pp)  &
                                             =  cmat(ir,it,p,pp) & 
                                             + (0.5*delta_t) &
                                             * omega_trap(it,is) &
                                             * sqrt(energy(ie)) &
                                             * (-1.0*indx_xi(ix)) &
                                             * (indx_xi(ix)-1.0)&
                                             / (2*indx_xi(ix)-1.0)
                                        amat(p,pp)  &
                                             =  amat(p,pp) & 
                                             - (0.5*delta_t) &
                                             * omega_trap(it,is) &
                                             * sqrt(energy(ie)) &
                                             * (-1.0*indx_xi(ix)) &
                                             * (indx_xi(ix)-1.0)&
                                             / (2*indx_xi(ix)-1.0)
                                     endif
                                  endif
                               endif

                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             enddo

             ! H_bar = (1 - dt/2 C - Poisson)^(-1) * (1 + dt/2 C + Poisson) H
             ! Lapack factorization and inverse of LHS
             call ZGETRF(msize,msize,cmat(ir,it,:,:),msize,i_piv,info)
             call ZGETRI(msize,cmat(ir,it,:,:),msize,i_piv,work,msize,info)
             ! Matrix multiply
             call ZGEMM('N','N',msize,msize,msize,alpha,cmat(ir,it,:,:),&
                  msize,amat,msize,beta,bmat,msize)
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
    use gkcoll_globals
    use gkcoll_poisson
    use gkcoll_gk, only : xi_mat
    use gkcoll_gyro
    implicit none
    integer :: is,ir,it,ie,ix
    integer :: p
    complex :: alpha = (1.0,0.0)
    complex :: beta  = (0.0,0.0)

    if(collision_model == -1) return
    if(imp_flag==2) return

    ! compute new collisional cap_H: H = h + ze/T G phi
    
    do ir=1,n_radial
       do it=1,n_theta
          
          ! set-up the RHS: H = f + ze/T G phi
          do is=1,n_species
             do ie=1,n_energy
                do ix=1,n_xi
                   p = indx_coll(is,ie,ix)
                   cvec(p) = cap_h_p(is,ir,it,ie,ix)
                enddo
             enddo
          enddo

          ! Solve for H
          call ZGEMV('N',msize,msize,alpha,cmat(ir,it,:,:),&
               msize,cvec,1,beta,bvec,1)
          
          do is=1,n_species
             do ie=1,n_energy
                do ix=1,n_xi
                   p = indx_coll(is,ie,ix)
                   cap_h_p(is,ir,it,ie,ix) = bvec(p)
                enddo
             enddo
          enddo
          
          ! Compute the new phi
          call POISSONp_do

       enddo
    enddo
    
    ! Compute the new h_x
    do is=1,n_species
       do ir=1,n_radial
          do it=1,n_theta
             do ie=1,n_energy
                call ZGEMV('N',n_xi,n_xi,alpha,xi_mat,n_xi,&
                     cap_h_p(is,ir,it,ie,:),1,beta,cap_h_x(is,ir,it,ie,:),1)
                do ix=1,n_xi
                   h_x(is,ir,it,ie,ix) = cap_h_x(is,ir,it,ie,ix) &
                        - z(is)/temp(is) * gyrox_J0(is,ir,it,ie,ix) &
                        * phi(ir,it)
                enddo
             enddo
          enddo
       enddo
    enddo
    

  end subroutine COLLISION_do
  
end module gkcoll_collision
