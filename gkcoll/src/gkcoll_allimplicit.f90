module gkcoll_allimplicit
  implicit none

  public :: GKallimp_alloc, GKallimp_do

  logical, private :: initialized = .false.

  complex, dimension(:,:), allocatable, private :: gmat
  complex, dimension(:), allocatable, private :: cvec, bvec, neovec
  integer, dimension(:,:,:,:,:), allocatable, private :: indx_gmat
  integer,private :: msize

contains

  subroutine GKallimp_alloc(flag)
    use gkcoll_globals
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    integer :: ir,it,is,ie,ix,p

    if(flag == 1) then
       if(initialized) return

       msize = n_radial * n_theta * n_species * n_energy * n_xi

       allocate(gmat(msize,msize))
       allocate(cvec(msize))
       allocate(bvec(msize))
       if(neoclassical_model == 1) then
          allocate(neovec(msize))
       endif
       allocate(indx_gmat(n_radial,n_theta,n_species,n_energy,n_xi))

       p = 0
       do ir=1,n_radial
          do it=1,n_theta
             do is=1,n_species
                do ie=1,n_energy
                   do ix=1,n_xi
                      p = p + 1
                      indx_gmat(ir,it,is,ie,ix) = p
                   enddo
                enddo
             enddo
          enddo
       enddo

       call GKallimp_matinit

       initialized = .true.
       
    else
       if(.NOT. initialized) return

       deallocate(gmat)
       deallocate(cvec)
       deallocate(bvec)
       if(neoclassical_model==1) then
          deallocate(neovec)
       endif
       deallocate(indx_gmat)

       initialized = .false.
    endif

  end subroutine GKallimp_alloc

  subroutine GKallimp_matinit
    use gkcoll_globals
    use gkcoll_gyro
    use gkcoll_equilibrium
    implicit none
    integer :: p, pp, ir, jr, it, jt, id, is, js, ie, je, ix, jx, ks
    real, external :: derf
    real :: xa, xb, tauinv_ab
    real, dimension(:,:,:), allocatable :: nu_d, nu_s
    real, dimension(:,:), allocatable :: rs
    real :: sum_nu, sum_den
    real, dimension(-2:2) :: cderiv
    integer, dimension(:), allocatable :: thcyc
    complex :: thfac, val
    real :: rval
    complex :: alpha = (1.0,0.0)
    complex :: beta  = (0.0,0.0) 
    complex, dimension(:,:), allocatable :: gexp, gtemp
    ! parameters for matrix solve
    integer :: info
    integer, dimension(:), allocatable :: i_piv
    complex, dimension(:), allocatable :: work
    ! n=0 test
    complex, dimension(:,:,:), allocatable :: pzf

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

    sum_den = 0.0
    do is=1,n_species
       sum_den = sum_den + z(is)**2 * dens(is) / temp(is)
    enddo
    if(adiabatic_ele_model == 1) then
       sum_den = sum_den + dens_ele / temp_ele
    endif

    if(toroidal_model == 2 .and. adiabatic_ele_model == 1 &
         .and. neoclassical_model /= 1) then
       allocate(pzf(n_radial,n_theta,n_theta))
       pzf(:,:,:) = (0.0,0.0)
       do ir=1,n_radial
          do it=1,n_theta
             pzf(ir,it,it) = -k_perp(it,ir)**2 * lambda_debye**2 &
                  * dens_ele / temp_ele + sum_den
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
    endif

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

    allocate(gexp(msize,msize))

    gmat(:,:) = (0.0,0.0)
    gexp(:,:) = (0.0,0.0)
    do ir=1,n_radial
       do it=1,n_theta
          do is=1,n_species
             do ie=1,n_energy
                do ix=1,n_xi
                   p = indx_gmat(ir,it,is,ie,ix)

                   ! identity
                   gmat(p,p) = gmat(p,p) + 1.0
                   gexp(p,p) = gexp(p,p) + 1.0

                   ! Collisions: Lorentz
                   js = is; je = ie; jx = ix
                   jr = ir; jt = it
                   pp = indx_gmat(jr,jt,js,je,jx)
                   sum_nu = 0.0
                   do ks=1,n_species
                      sum_nu = sum_nu + nu_d(ie,is,ks)
                   enddo
                   val = (0.5*delta_t) &
                        * 0.5 * (-1.0*indx_xi(ix)) &
                        * (indx_xi(ix)+1.0) * sum_nu
                   gmat(p,pp) = gmat(p,pp) - val
                   gexp(p,pp) = gexp(p,pp) + val

                   ! Collisions: Restoring terms
                   if(indx_xi(ix) == 1) then
                      jr = ir; jt = it
                      jx = ix
                      do je=1,n_energy
                         do js=1,n_species
                            pp = indx_gmat(jr,jt,js,je,jx)
                            if(abs(rs(is,js)) < epsilon(0.)) then
                               val = 0.0
                            else
                               val = (0.5*delta_t) &
                                    * (mass(js)/mass(is))*(dens(js)/dens(is)) &
                                    * (vth(js)/vth(is)) &
                                    * nu_s(ie,is,js) * sqrt(energy(ie)) &
                                    * nu_s(je,js,is) * sqrt(energy(je)) &
                                    * w_e(je)/ rs(is,js)
                            endif
                            gmat(p,pp) = gmat(p,pp) - val
                            gexp(p,pp) = gexp(p,pp) + val
                         enddo
                      enddo
                      if(collision_model == 2) then
                         je=ie; js=is
                         pp = indx_gmat(jr,jt,js,je,jx)
                         sum_nu = 0.0
                         do ks=1,n_species
                            sum_nu = sum_nu + (nu_d(ie,is,ks)-nu_s(ie,is,ks))
                         enddo
                         val = (0.5*delta_t) * sum_nu
                         gmat(p,pp) = gmat(p,pp) - val
                         gexp(p,pp) = gexp(p,pp) + val
                      endif
                   endif


                   ! Trapping
                   js = is; je = ie
                   jr = ir; jt = it
                   jx = ix-1
                   if(jx >= 1) then
                      pp = indx_gmat(jr,jt,js,je,jx)
                      val = (0.5*delta_t) &
                           * omega_trap(it,is) &
                           * sqrt(energy(ie)) &
                           * (-1.0*indx_xi(ix)) &
                           * (indx_xi(ix)-1.0)&
                           / (2*indx_xi(ix)-1.0)
                      gmat(p,pp) = gmat(p,pp) + val
                      gexp(p,pp) = gexp(p,pp) - val
                   endif
                   jx = ix+1
                   if(jx <= n_xi) then
                      pp = indx_gmat(jr,jt,js,je,jx)
                      val = (0.5*delta_t) &
                           * omega_trap(it,is) &
                           * sqrt(energy(ie)) &
                           * (indx_xi(ix)+1.0) &
                           * (indx_xi(ix)+2.0) &
                           / (2*indx_xi(ix)+3.0)
                      gmat(p,pp) = gmat(p,pp) + val
                      gexp(p,pp) = gexp(p,pp) - val
                   endif

                   ! Streaming
                   js = is; je = ie
                   rval = omega_stream(it,is) * sqrt(energy(ie))
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
                      jx=ix-1
                      if (jx >= 1) then
                         pp = indx_gmat(jr,jt,js,je,jx)
                         val = (0.5*delta_t) * rval * cderiv(id) * thfac &
                              * (1.0*indx_xi(ix)) / (2*indx_xi(ix)-1.0)
                         gmat(p,pp) = gmat(p,pp) + val
                         gexp(p,pp) = gexp(p,pp) - val
                      endif
                      jx=ix+1
                      if (jx <= n_xi) then
                         pp = indx_gmat(jr,jt,js,je,jx)
                         val = (0.5*delta_t) * rval * cderiv(id) * thfac &
                              * (indx_xi(ix)+1.0) / (2*indx_xi(ix)+3.0)
                         gmat(p,pp) = gmat(p,pp) + val
                         gexp(p,pp) = gexp(p,pp) - val     
                      endif
                   enddo

                   ! Radial and alpha drift
                   js = is; je = ie
                   jr = ir; jt = it
                   jx = ix
                   pp = indx_gmat(jr,jt,js,je,jx)
                   val = (0.5*delta_t) * energy(ie) &
                        * (omega_rdrift(it,is) &
                        * (2.0*pi*i_c*indx_r(ir)*r_length_inv) &
                        + omega_adrift(it,is) * i_c * k_theta) &
                        * (1.0 + (indx_xi(ix)+1.0)*(indx_xi(ix)+1.0) &
                        / ( (2*indx_xi(ix)+1.0) * (2*indx_xi(ix)+3.0) ) &
                        + (1.0*indx_xi(ix)) * (1.0*indx_xi(ix)) &
                        / ( (2*indx_xi(ix)+1.0) * (2*indx_xi(ix)-1.0) ) )
                   gmat(p,pp) = gmat(p,pp) + val
                   gexp(p,pp) = gexp(p,pp) - val
                   jx = ix-2
                   if(jx >= 1) then
                      pp = indx_gmat(jr,jt,js,je,jx)
                      val = (0.5*delta_t) * energy(ie) &
                           * (omega_rdrift(it,is) &
                           * (2.0*pi*i_c*indx_r(ir)*r_length_inv) &
                           + omega_adrift(it,is) * i_c * k_theta) &
                           * ( (1.0*indx_xi(ix)) * (1.0*indx_xi(ix)-1.0) &
                           / ( (2*indx_xi(ix)-3.0) * (2*indx_xi(ix)-1.0) ) )
                      gmat(p,pp) = gmat(p,pp) + val
                      gexp(p,pp) = gexp(p,pp) - val
                   endif
                   jx = ix+2
                   if(jx <= n_xi) then
                      pp = indx_gmat(jr,jt,js,je,jx)
                      val = (0.5*delta_t) * energy(ie) &
                           * (omega_rdrift(it,is) &
                           * (2.0*pi*i_c*indx_r(ir)*r_length_inv) &
                           + omega_adrift(it,is) * i_c * k_theta) &
                           * ( (1.0*indx_xi(ix)+1.0) * (1.0*indx_xi(ix)+2.0) &
                           / ( (2*indx_xi(ix)+3.0) * (2*indx_xi(ix)+5.0) ) )
                      gmat(p,pp) = gmat(p,pp) + val
                      gexp(p,pp) = gexp(p,pp) - val
                   endif

                   ! phi: (-ze/T) * G * (1-dt/2 * i*omega_star) * phi
                   if(toroidal_model == 2 .and. adiabatic_ele_model == 1 &
                        .and. neoclassical_model /= 1) then
                      val = (0.5*delta_t) * i_c * k_theta &
                           * rho *sqrt(temp(is)*mass(is))/(1.0*z(is)) &
                           * vth(is) &
                           * (dlnndr(is) + dlntdr(is) * (energy(ie)-1.5))
                      jr=ir
                      do js=1,n_species
                         do je=1,n_energy
                            do jx=1,n_xi
                               do jt=1,n_theta
                                  pp = indx_gmat(jr,jt,js,je,jx)
                                  gmat(p,pp)  &
                                       =  gmat(p,pp) &
                                       - z(is)/temp(is) &
                                       * pzf(ir,it,jt) &
                                       * (2.0*indx_xi(ix) + 1) &
                                       * gyrop_J0(is,ir,it,ie,ix) &
                                       * z(js)*dens(js) &
                                       * gyrop_J0(js,ir,jt,je,jx) * w_e(je) &
                                       * (1.0 - val)
                                  gexp(p,pp)  &
                                    =  gexp(p,pp) &
                                    - z(is)/temp(is) &
                                    * pzf(ir,it,jt) &
                                    * (2.0*indx_xi(ix) + 1) &
                                    * gyrop_J0(is,ir,it,ie,ix) &
                                    * z(js)*dens(js) &
                                    * gyrop_J0(js,ir,jt,je,jx) * w_e(je) &
                                    * (1.0 + val)
                               enddo
                            enddo
                         enddo
                      enddo
                   else
                      val = (0.5*delta_t) * i_c * k_theta &
                           * rho *sqrt(temp(is)*mass(is))/(1.0*z(is)) &
                           * vth(is) &
                           * (dlnndr(is) + dlntdr(is) * (energy(ie)-1.5))
                      jr=ir; jt=it
                      do js=1,n_species
                         do je=1,n_energy
                            do jx=1,n_xi
                               pp = indx_gmat(jr,jt,js,je,jx)
                               gmat(p,pp)  &
                                    =  gmat(p,pp) &
                                    - z(is)/temp(is) / &
                                    (-k_perp(it,ir)**2 * lambda_debye**2 &
                                    * dens_ele / temp_ele + sum_den) &
                                    * (2.0*indx_xi(ix) + 1) &
                                    * gyrop_J0(is,ir,it,ie,ix) &
                                    * z(js)*dens(js) &
                                    * gyrop_J0(js,ir,it,je,jx) * w_e(je) &
                                    * (1.0 - val)
                               gexp(p,pp)  &
                                    =  gexp(p,pp) &
                                    - z(is)/temp(is) / &
                                    (-k_perp(it,ir)**2 * lambda_debye**2 &
                                    * dens_ele / temp_ele + sum_den) &
                                    * (2.0*indx_xi(ix) + 1) &
                                    * gyrop_J0(is,ir,it,ie,ix) &
                                    * z(js)*dens(js) &
                                    * gyrop_J0(js,ir,it,je,jx) * w_e(je) &
                                    * (1.0 + val)
                            enddo
                         enddo
                      enddo
                   endif

                enddo
             enddo
          enddo
       enddo
    enddo
                                  
    deallocate(thcyc)
    deallocate(nu_d)
    deallocate(nu_s)
    deallocate(rs)
    if(toroidal_model == 2 .and. adiabatic_ele_model == 1 &
         .and. neoclassical_model /= 1) then
       deallocate(pzf)
    endif

    ! Lapack factorization and inverse
    allocate(work(msize))
    allocate(i_piv(msize))
    call ZGETRF(msize,msize,gmat,msize,i_piv,info)
    call ZGETRI(msize,gmat,msize,i_piv,work,msize,info)
    deallocate(i_piv)
    deallocate(work)
    ! Matrix multiply
    allocate(gtemp(msize,msize))
    call ZGEMM('N','N',msize,msize,msize,alpha,gmat,&
         msize,gexp,msize,beta,gtemp,msize)
    deallocate(gexp)
    
    if(neoclassical_model == 1) then
       cvec(:) = (0.0,0.0)
       do ir=1,n_radial
          do it=1,n_theta
             do is=1,n_species
                do ie=1,n_energy
                   do ix=1,n_xi
                      if(indx_xi(ix) == 0) then
                         p = indx_gmat(ir,it,is,ie,ix)
                         cvec(p) = (4.0/3.0) * delta_t * energy(ie) &
                              * omega_rdrift(it,is) &
                              * (dlnndr(is) + dlntdr(is) * (energy(ie)-1.5))
                      else if(indx_xi(ix) == 2) then
                         p = indx_gmat(ir,it,is,ie,ix)
                         cvec(p) = (2.0/3.0) * delta_t * energy(ie) &
                              * omega_rdrift(it,is) &
                              * (dlnndr(is) + dlntdr(is) * (energy(ie)-1.5))
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
       call ZGEMV('N',msize,msize,alpha,gmat,&
            msize,cvec,1,beta,neovec,1)
    endif

    gmat = gtemp
    deallocate(gtemp)

  end subroutine GKallimp_matinit

  subroutine GKallimp_do
    use gkcoll_globals
    use gkcoll_poisson
    use gkcoll_gk, only : xi_mat
    use gkcoll_gyro
    implicit none
    integer :: ir,it,is,ie,ix,p
    complex :: alpha = (1.0,0.0)
    complex :: beta  = (0.0,0.0)

    cvec(:) = (0.0,0.0)
    ! set-up the RHS
    do ir=1,n_radial
       do it=1,n_theta
          do is=1,n_species
             do ie=1,n_energy
                do ix=1,n_xi
                   p = indx_gmat(ir,it,is,ie,ix)
                   cvec(p) = cap_h_p(is,ir,it,ie,ix)
                enddo
             enddo
          enddo
       enddo
    enddo

    ! Solve for H
    call ZGEMV('N',msize,msize,alpha,gmat,&
         msize,cvec,1,beta,bvec,1)

    ! Neoclassical driving term
    if(neoclassical_model == 1) then
       do p=1,msize
          bvec(p) = bvec(p) + neovec(p)
       enddo
    endif

    do ir=1,n_radial
       do it=1,n_theta
          do is=1,n_species
             do ie=1,n_energy
                do ix=1,n_xi
                   p = indx_gmat(ir,it,is,ie,ix)
                   cap_h_p(is,ir,it,ie,ix) = bvec(p)
                enddo
             enddo
          enddo
       enddo
    enddo
    
    
    ! Compute the new phi
    call POISSONp_do

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
    
  end subroutine GKallimp_do
  
end module gkcoll_allimplicit
