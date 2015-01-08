module cgyro_collision

  implicit none

  public :: COLLISION_alloc, cgyro_step_collision
  logical, private :: initialized = .false.

  real, dimension(:,:,:), allocatable, private :: cmat
  complex, dimension(:), allocatable, private :: cvec,bvec

contains

  subroutine COLLISION_alloc(flag)

    use timer_lib

    use cgyro_globals
    use cgyro_equilibrium, only : omega_trap, k_perp, bmag

    implicit none

    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    real, dimension(:,:,:), allocatable :: nu_d, nu_s, nu_par, nu_par_deriv
    real, dimension(:,:), allocatable :: rs
    real, dimension(:,:,:,:), allocatable :: rsvec
    real, external :: derf
    real :: xa, xb, tauinv_ab
    real :: sum_den
    integer :: jv
    integer :: is,ir,it,ix,ie,js,je,jx,ks
    ! parameters for matrix solve
    real, dimension(:,:), allocatable :: amat, bmat
    real, dimension(:,:,:,:,:,:), allocatable :: ctest, cfield
    real, dimension(:,:,:,:,:,:,:), allocatable :: ctest_k, cfield_k
    real, external :: BESJ0
    real :: arg1, arg2
    real :: bessel(0:2)
    integer :: ierr

    if (collision_model == 0) return

    if (flag == 1) then

       if(initialized) return

       allocate(nu_d(n_energy,n_species,n_species))
       allocate(nu_s(n_energy,n_species,n_species))
       allocate(nu_par(n_energy,n_species,n_species))
       allocate(nu_par_deriv(n_energy,n_species,n_species))
       nu_d(:,:,:) = 0.0
       nu_s(:,:,:) = 0.0
       nu_par(:,:,:) = 0.0
       nu_par_deriv(:,:,:) = 0.0

       call timer_lib_in('coll_set1')

       do ie=1,n_energy
          do is=1,n_species
             do js=1,n_species

                xa = sqrt(energy(ie))
                xb = xa * vth(is) / vth(js)
                tauinv_ab = nu(is) * (1.0*Z(js))**2 / (1.0*Z(is))**2 &
                     * dens(js)/dens(is)

                select case (collision_model)

                case (1)

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

                case (2)

                   ! Connor model
                   if(is == js .or. &
                        (abs(mass(is) - mass(js)) < epsilon(0.))) then
                      ! case 1: like-species/same mass collisions
                      nu_d(ie,is,js) = tauinv_ab * (1.0/xa**3) &
                           * (exp(-xb*xb)/(xb*sqrt(pi)) &
                           + (1.0-1.0/(2.0*xb*xb)) * DERF(xb))

                   else if (mass(is) < mass(js)) then
                      ! case 2: ele-ion and ion-imp(heavy) collisions
                      nu_d(ie,is,js) = tauinv_ab * (1.0/xa**3)

                   else
                      ! case 3: ion-ele and imp(heavy)-ion collisions

                      nu_d(ie,is,js) = tauinv_ab * 4.0/(3.0*sqrt(pi)) &
                           * sqrt(mass(js)/mass(is)) &
                           * (temp(is)/temp(js))**1.5
                   endif
                   nu_s(ie,is,js) = nu_d(ie,is,js)

                case (3)

                   ! Reduced Hirshman-Sigmar model
                   nu_d(ie,is,js) = tauinv_ab * (1.0/xa**3) &
                        * (exp(-xb*xb)/(xb*sqrt(pi)) &
                        + (1.0-1.0/(2.0*xb*xb)) * DERF(xb))
                   nu_s(ie,is,js) = tauinv_ab * (1.0/xa) &
                        * (-exp(-xb*xb)/(xb*sqrt(pi)) &
                        + (1.0/(2.0*xb*xb)) * DERF(xb)) &
                        * (2.0*temp(is)/temp(js))*(1.0+mass(js)/mass(is))

                case(4)
                   
                   ! Ad hoc op
                   nu_d(ie,is,js) = tauinv_ab * (1.0/xa**3) &
                        * (exp(-xb*xb)/(xb*sqrt(pi)) &
                        + (1.0-1.0/(2.0*xb*xb)) * DERF(xb))
                   ! No i-e Lorentz
                   if(is /= is_ele .and. js == is_ele) then
                      nu_d(ie,is,js) = 0.0
                   endif

                   ! Only ii, ee Diffusion
                   if(is == js) then
                      nu_par(ie,is,js) = tauinv_ab * (2.0/xa**3) &
                           * (-exp(-xb*xb)/(xb*sqrt(pi)) &
                           + (1.0/(2.0*xb*xb)) * DERF(xb))
                      nu_par_deriv(ie,is,js) = tauinv_ab * (2.0/xa**3) &
                           * (-3/xa * (-exp(-xb*xb)/(xb*sqrt(pi)) &
                           + (1.0/(2.0*xb*xb)) * DERF(xb)) &
                           + vth(is)/vth(js) &
                           * (2.0*exp(-xb*xb)/(xb**2*sqrt(pi)) &
                           + 2.0*exp(-xb*xb)/sqrt(pi) - DERF(xb)/xb**3))
                   endif

                end select

             enddo
          enddo
       enddo

       allocate(ctest(n_species,n_species,n_xi,n_xi,n_energy,n_energy))
       allocate(cfield(n_species,n_species,n_xi,n_xi,n_energy,n_energy))
       allocate(ctest_k(n_species,n_species,n_xi,n_xi,n_energy,n_energy,&
            ic_loc))
       allocate(cfield_k(n_species,n_species,n_xi,n_xi,n_energy,n_energy,&
            ic_loc))
       allocate(rs(n_species,n_species))
       allocate(rsvec(n_species,n_species,n_xi,n_energy))

       ! Collision test particle component
       ctest   = 0.0
       ctest_k = 0.0

       ! Lorentz
       do is=1,n_species
          do ix=1,n_xi
             do ie=1,n_energy
                do js=1, n_species
                   do jx=1, n_xi
                      je = ie
                      ctest(is,js,ix,jx,ie,je) &
                           = ctest(is,js,ix,jx,ie,je) &
                           + xi_lor_mat(ix,jx) *0.5*nu_d(ie,is,js)
                   enddo
                enddo
             enddo
          enddo
       enddo

       ! U factor for HS0
       if (collision_model == 3 .and. collision_mom_restore == 1) then
          do is=1,n_species
             do ix=1,n_xi
                do ie=1,n_energy
                   do js=1, n_species
                      do jx=1, n_xi
                         je = ie
                         ctest(is,js,ix,jx,ie,je) &
                              = ctest(is,js,ix,jx,ie,je) &
                              + (nu_d(ie,is,js)-nu_s(ie,is,js)) &
                              * 1.5 * xi(ix) * xi(jx) * w_xi(jx)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       endif

       ! Diffusion
       if(collision_model == 4 .and. collision_ene_diffusion == 1) then
          do is=1,n_species 
             do ix=1,n_xi
                do ie=1,n_energy
                   do js=1, n_species
                      do je=1, n_energy
                         jx=ix                       
                         if(ie == je) then
                            ctest(is,js,ix,jx,ie,je) &
                                 = ctest(is,js,ix,jx,ie,je) &
                                 + (1.0-temp(is)/temp(js)) &
                                 * (-nu_par_deriv(ie,is,js) &
                                 * energy(ie)**1.5 & 
                                 + nu_par(ie,is,js) & 
                                 * (2.0*energy(ie)**2 &
                                 - 5.0*energy(ie)))
                         endif
                         ctest(is,js,ix,jx,ie,je) &
                              = ctest(is,js,ix,jx,ie,je) &
                              + nu_par(ie,is,js) * 0.5 *energy(ie) &
                              * e_deriv2_mat(ie,je) / e_max
                         ctest(is,js,ix,jx,ie,je) &
                              = ctest(is,js,ix,jx,ie,je) &
                              + e_deriv1_mat(ie,je)/sqrt(1.0*e_max) &
                              * (nu_par_deriv(ie,is,js) &
                              * 0.5*energy(ie) &
                              + nu_par(ie,is,js) &
                              * (2.0*sqrt(energy(ie)) &
                              + (temp(is)/temp(js)-2.0) &
                              * energy(ie)**1.5))
                      enddo
                   enddo
                enddo
             enddo
          enddo
       endif

       ! Finite-kperp corrections
       if(collision_model == 4 .and. collision_kperp == 1) then
          ic_loc = 0
          do ic=nc1,nc2
             ic_loc = ic_loc+1
             it = it_c(ic)
             ir = ir_c(ic)
             do is=1,n_species   
                do ix=1,n_xi
                   do ie=1,n_energy
                      do js=1, n_species
                         jx=ix
                         je=ie
                         ctest_k(is,js,ix,jx,ie,je,ic_loc) &
                              = ctest_k(is,js,ix,jx,ie,je,ic_loc) &
                              - 0.25*(k_perp(it,ir)*rho*vth(is)*mass(is) &
                              / (z(is)*Bmag(it)))**2 &
                              * 2.0*energy(ie) &
                              * (nu_d(ie,is,js) * (1+xi(ix)**2) &
                              + nu_par(ie,is,js) * (1-xi(ix)**2))
                      enddo
                   enddo
                enddo
             enddo
          enddo
       endif
       
       ! Collision field particle component
       cfield   = 0.0
       cfield_k = 0.0

       select case (collision_model)

       case(1,2,3)
          if(collision_mom_restore == 1) then
             do is=1,n_species
                do js=1,n_species
                   rs(is,js) = 0.0
                   do ie=1,n_energy
                      rs(is,js) = rs(is,js) &
                           + w_e(ie)*nu_s(ie,is,js)*energy(ie)
                   enddo
                enddo
             enddo
             do is=1,n_species
                do js=1, n_species   
                   do ix=1,n_xi
                      do jx=1, n_xi
                         do ie=1,n_energy
                            do je=1, n_energy
                               if (abs(rs(is,js)) > epsilon(0.)) then
                                  cfield(is,js,ix,jx,ie,je) = &
                                       cfield(is,js,ix,jx,ie,je) &
                                       + 1.5 * (mass(js)/mass(is)) &
                                       * (dens(js)/dens(is)) &
                                       * (vth(js)/vth(is)) &
                                       * nu_s(ie,is,js) * sqrt(energy(ie)) &
                                       * xi(ix) &
                                       * nu_s(je,js,is) * sqrt(energy(je)) &
                                       * xi(jx) * w_e(je) * w_xi(jx) &
                                       / rs(is,js)
                               endif
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          endif

       case(4)

          ! Momentum Restoring

          if(collision_mom_restore == 1) then
             ! C_test_ab(v_par f0a,f0b) / vth_a
             rsvec = 0.0
             do is=1,n_species
                do js=1,n_species
                   do ix=1,n_xi
                      do jx=1,n_xi
                         do ie=1,n_energy
                            do je=1,n_energy
                               rsvec(is,js,ix,ie) = rsvec(is,js,ix,ie) &
                                    + ctest(is,js,ix,jx,ie,je) &
                                    * sqrt(2.0*energy(je)) * xi(jx) 
                            enddo
                         enddo
                      enddo
                   enddo

                   ! int v_par C_test_ab(v_par f0a,f0b) / (n_0a vth_a^2)
                   rs(is,js) = 0.0
                   do ix=1,n_xi
                      do ie=1,n_energy
                         rs(is,js) = rs(is,js) + 0.5*w_e(ie)*w_xi(ix) &
                              * rsvec(is,js,ix,ie) &
                              * sqrt(2.0*energy(ie)) * xi(ix) 
                      enddo
                   enddo
                enddo
             enddo

             if(collision_kperp == 0) then
                do is=1,n_species
                   do js=1, n_species         
                      do ix=1,n_xi
                         do jx=1, n_xi
                            do ie=1,n_energy
                               do je=1, n_energy
                                  if (abs(rs(is,js))>epsilon(0.0)) then
                                     cfield(is,js,ix,jx,ie,je) &
                                          = cfield(is,js,ix,jx,ie,je) &
                                          - mass(js)/mass(is) &
                                          * (dens(js)/dens(is)) &
                                          * (vth(js)/vth(is)) &
                                          * rsvec(is,js,ix,ie) &
                                          / rs(is,js) &
                                          * rsvec(js,is,jx,je) &
                                          * 0.5*w_e(je)*w_xi(jx)
                                  endif
                               enddo
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
                
             else
                ic_loc = 0
                do ic=nc1,nc2
                   ic_loc = ic_loc+1
                   it = it_c(ic)
                   ir = ir_c(ic)
                   do is=1,n_species  
                      do js=1, n_species   
                         do ix=1,n_xi
                            do jx=1, n_xi
                               do ie=1,n_energy
                                  do je=1, n_energy
                                     if (abs(rs(is,js)) > epsilon(0.)) then
                                        arg1 = k_perp(it,ir)*rho*vth(is)*mass(is)&
                                             /(z(is)*Bmag(it)) &
                                             *sqrt(2.0*energy(ie))&
                                             *sqrt(1.0-xi(ix)**2)
                                        arg2 = k_perp(it,ir)*rho*vth(js)*mass(js)&
                                             /(z(js)*Bmag(it)) &
                                             *sqrt(2.0*energy(je))&
                                             *sqrt(1.0-xi(jx)**2)
                                        call RJBESL(abs(arg2),0.0,2,bessel,ierr)

                                        cfield_k(is,js,ix,jx,ie,je,ic_loc) &
                                             = cfield_k(is,js,ix,jx,ie,je,ic_loc) &
                                             - mass(js)/mass(is) &
                                             * (dens(js)/dens(is)) &
                                             * (vth(js)/vth(is))**2 &
                                             * rsvec(is,js,ix,ie) &
                                             * BESJ0(abs(arg1)) &
                                             / rs(is,js) &
                                             * rsvec(js,is,jx,je) &
                                             * bessel(1) &
                                             * 0.5*w_xi(jx)*w_e(je)
                                        cfield_k(is,js,ix,jx,ie,je,ic_loc) &
                                             = cfield_k(is,js,ix,jx,ie,je,ic_loc) &
                                             - mass(js)/mass(is) &
                                             * (dens(js)/dens(is)) &
                                             * (vth(js)/vth(is))**2 &
                                             * rsvec(is,js,ix,ie) &
                                             * BESJ1(abs(arg1)) &
                                             * sqrt(1.0-xi(ix)**2)/xi(ix) &
                                             / rs(is,js) &
                                             * rsvec(js,is,jx,je) &
                                             * bessel(1) &
                                             * sqrt(1.0-xi(jx)**2)/xi(jx) &
                                             * 0.5*w_xi(jx)*w_e(je)
                                     endif
                                  enddo
                               enddo
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
                
             endif

          endif

          ! Energy Restoring

          if(collision_ene_restore == 1) then
             ! C_test_ab(v^2 f0a,f0b) / vth_a^2
             rsvec = 0.0
             do is=1,n_species
                do js=1,n_species
                   do ix=1,n_xi
                      do jx=1,n_xi
                         do ie=1,n_energy
                            do je=1,n_energy
                               rsvec(is,js,ix,ie) = rsvec(is,js,ix,ie) &
                                    + ctest(is,js,ix,jx,ie,je) &
                                    * 2.0*energy(je)
                            enddo
                         enddo
                      enddo
                   enddo
                   
                   ! int v^2 C_test_ab(v^2 f0a,f0b) / (n_0a vth_a^4)
                   rs(is,js) = 0.0
                   do ix=1,n_xi
                      do ie=1,n_energy
                         rs(is,js) = rs(is,js) + 0.5*w_e(ie)*w_xi(ix) &
                              * rsvec(is,js,ix,ie) &
                              * 2.0*energy(ie) 
                      enddo
                   enddo                    
                enddo
             enddo

             if(collision_kperp == 0) then
                do is=1,n_species
                   do js=1, n_species   
                      do ix=1,n_xi
                         do jx=1, n_xi
                            do ie=1,n_energy
                               do je=1, n_energy
                                  if (abs(rs(is,js)) > epsilon(0.)) then
                                     cfield(is,js,ix,jx,ie,je) &
                                          = cfield(is,js,ix,jx,ie,je) &
                                          - mass(js)/mass(is) &
                                          * (dens(js)/dens(is)) &
                                          * (vth(js)/vth(is))**2 &
                                          * rsvec(is,js,ix,ie) &
                                          / rs(is,js) &
                                          * rsvec(js,is,jx,je) &
                                          * 0.5*w_xi(jx)*w_e(je)
                                  endif
                               enddo
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
                
             else
                ic_loc = 0
                do ic=nc1,nc2
                   ic_loc = ic_loc+1
                   it = it_c(ic)
                   ir = ir_c(ic)
                   do is=1,n_species  
                      do js=1, n_species   
                         do ix=1,n_xi
                            do jx=1, n_xi
                               do ie=1,n_energy
                                  do je=1, n_energy
                                     if (abs(rs(is,js)) > epsilon(0.)) then
                                        arg1 = k_perp(it,ir)*rho*vth(is)*mass(is)&
                                             /(z(is)*Bmag(it)) &
                                             *sqrt(2.0*energy(ie))&
                                             *sqrt(1.0-xi(ix)**2)
                                        arg2 = k_perp(it,ir)*rho*vth(js)*mass(js)&
                                             /(z(js)*Bmag(it)) &
                                             *sqrt(2.0*energy(je))&
                                             *sqrt(1.0-xi(jx)**2)
                                        call RJBESL(abs(arg2),0.0,2,bessel,ierr)
                                        cfield_k(is,js,ix,jx,ie,je,ic_loc) &
                                             = cfield_k(is,js,ix,jx,ie,je,ic_loc) &
                                             - mass(js)/mass(is) &
                                             * (dens(js)/dens(is)) &
                                             * (vth(js)/vth(is))**2 &
                                             * rsvec(is,js,ix,ie) &
                                             * BESJ0(abs(arg1)) &
                                             / rs(is,js) &
                                             * rsvec(js,is,jx,je) &
                                             * bessel(1) &
                                             * 0.5*w_xi(jx)*w_e(je)
                                     endif
                                  enddo
                               enddo
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
                
             endif
             
          endif

       end select
       
       allocate(cmat(nv,nv,nc_loc))
       allocate(cvec(nv))
       allocate(bvec(nv))

       sum_den = 0.0
       do is=1,n_species
          do ie=1,n_energy
             do ix=1,n_xi
                sum_den = sum_den + 0.5 * w_xi(ix) &
                     * z(is)**2/temp(is) *dens(is) * w_e(ie)
             enddo
          enddo
       enddo
       if (ae_flag == 1) then
          sum_den = sum_den + dens_ele / temp_ele
       endif

       call timer_lib_out('coll_set1')

       ! matrix solve parameters
       allocate(work(nv))
       allocate(i_piv(nv))
       allocate(amat(nv,nv))
       allocate(bmat(nv,nv))

       call timer_lib_in('coll_set2')

       ! set-up the collision matrix
       ic_loc = 0
       do ic=nc1,nc2
          ic_loc = ic_loc+1

          it = it_c(ic)
          ir = ir_c(ic)

          cmat(:,:,ic_loc) = 0.0
          amat(:,:)        = 0.0

          do iv=1,nv

             is = is_v(iv)
             ix = ix_v(iv)
             ie = ie_v(iv)

             do jv=1,nv

                js = is_v(jv)
                jx = ix_v(jv)
                je = ie_v(jv)

                ! constant part
                if (iv == jv) then
                   cmat(iv,jv,ic_loc) =  cmat(iv,jv,ic_loc) + 1.0
                   amat(iv,jv) = amat(iv,jv) + 1.0
                endif

                ! Trapping term
                if (is == js .and. ie == je) then
                   cmat(iv,jv,ic_loc) = cmat(iv,jv,ic_loc) &
                        + (0.5*delta_t) * omega_trap(it,is) &
                        * sqrt(energy(ie)) &
                        * (1.0 - xi(ix)**2) &
                        * xi_deriv_mat(ix,jx)
                   amat(iv,jv) = amat(iv,jv) &
                        - (0.5*delta_t) * omega_trap(it,is) &
                        * sqrt(energy(ie)) &
                        * (1.0 - xi(ix)**2) &
                        * xi_deriv_mat(ix,jx)
                endif

                ! Collision component: Test particle
                if(is == js) then
                   do ks=1,n_species
                      cmat(iv,jv,ic_loc) = cmat(iv,jv,ic_loc) &
                           - (0.5*delta_t) * (ctest(is,ks,ix,jx,ie,je) &
                           + ctest_k(is,ks,ix,jx,ie,je,ic_loc))
                      amat(iv,jv) = amat(iv,jv) &
                           + (0.5*delta_t) * (ctest(is,ks,ix,jx,ie,je) &
                           + ctest_k(is,ks,ix,jx,ie,je,ic_loc))
                   enddo
                endif

                ! Collision component: Field particle
                cmat(iv,jv,ic_loc) = cmat(iv,jv,ic_loc) &
                     - (0.5*delta_t) * (cfield(is,js,ix,jx,ie,je) &
                     + cfield_k(is,js,ix,jx,ie,je,ic_loc))
                amat(iv,jv) = amat(iv,jv) &
                     + (0.5*delta_t) * (cfield(is,js,ix,jx,ie,je) &
                     + cfield_k(is,js,ix,jx,ie,je,ic_loc))

                ! Poisson component 
                if (n == 0 .and. ae_flag == 1) then
                   ! Cannot include Poisson in collision matrix
                   ! for n=0 with ade because depends on theta
                   ! i.e. ne0 ~ phi - <phi>
                   cmat(iv,jv,ic_loc) = cmat(iv,jv,ic_loc) + 0.0
                   amat(iv,jv)        = amat(iv,jv) + 0.0
                else
                   cmat(iv,jv,ic_loc) = cmat(iv,jv,ic_loc) &
                        - z(is)/temp(is) / &
                        (k_perp(it,ir)**2 * lambda_debye**2 &
                        * dens_ele / temp_ele &
                        + sum_den) &
                        * j0_v(ic_loc,iv) &
                        * z(js)*dens(js) &
                        * j0_v(ic_loc,jv) * w_e(je) &
                        * 0.5 * w_xi(jx)
                   amat(iv,jv) = amat(iv,jv) &
                        - z(is)/temp(is) / &
                        (k_perp(it,ir)**2 * lambda_debye**2 &
                        * dens_ele / temp_ele &
                        + sum_den) &
                        * j0_v(ic_loc,iv) &
                        * z(js)*dens(js) &
                        * j0_v(ic_loc,jv) * w_e(je) &
                        * 0.5 * w_xi(jx)
                endif

                ! Ampere component
                if (n_field > 1) then
                   cmat(iv,jv,ic_loc) = cmat(iv,jv,ic_loc) &
                        - z(is)/temp(is) / &
                        (2.0*k_perp(it,ir)**2 * rho**2 / betae_unit & 
                        * dens_ele * temp_ele) &
                        * (-j0_v(ic_loc,iv)) &
                        * z(js)*dens(js) &
                        * xi(ix) * sqrt(2.0*energy(ie)) *vth(is) &
                        * j0_v(ic_loc,jv) * w_e(je) &
                        * 0.5 * w_xi(jx) &
                        * xi(jx) * sqrt(2.0*energy(je)) * vth(is)
                   amat(iv,jv) = amat(iv,jv) &
                        - z(is)/temp(is) / &
                        (2.0*k_perp(it,ir)**2 * rho**2 / betae_unit & 
                        * dens_ele * temp_ele) &
                        * (-j0_v(ic_loc,iv)) &
                        * z(js)*dens(js) &
                        * xi(ix) * sqrt(2.0*energy(ie)) * vth(is) &
                        * j0_v(ic_loc,jv) * w_e(je) &
                        * 0.5 * w_xi(jx) &
                        * xi(jx) * sqrt(2.0*energy(je)) * vth(is)
                endif

             enddo
          enddo

          ! H_bar = (1 - dt/2 C - Poisson)^(-1) * (1 + dt/2 C + Poisson) H
          ! Lapack factorization and inverse of LHS
          call DGETRF(nv,nv,cmat(:,:,ic_loc),nv,i_piv,info)
          call DGETRI(nv,cmat(:,:,ic_loc),nv,i_piv,work,nv,info)
          ! Matrix multiply
          call DGEMM('N','N',nv,nv,nv,num1,cmat(:,:,ic_loc),&
               nv,amat,nv,num0,bmat,nv)
          cmat(:,:,ic_loc) = bmat(:,:)

       enddo

       call timer_lib_out('coll_set2')

       deallocate(amat)
       deallocate(bmat)
       deallocate(work)
       deallocate(i_piv)
       deallocate(nu_d)
       deallocate(nu_s)
       deallocate(nu_par)
       deallocate(nu_par_deriv)
       deallocate(rs)
       deallocate(rsvec)
       deallocate(ctest)
       deallocate(cfield)
       deallocate(ctest_k)
       deallocate(cfield_k)

       initialized = .true.

    else

       if(.NOT. initialized) return
       deallocate(cmat)
       deallocate(cvec)
       deallocate(bvec)
       initialized = .false.

    endif

  end subroutine COLLISION_alloc

  subroutine cgyro_step_collision

    use parallel_lib
    use timer_lib

    use cgyro_globals

    implicit none

    integer :: is,ir,it,ie,ix
    integer :: ivp

    if (collision_model == 0) return

    ! compute new collisional cap_H: H = h + ze/T G phi
    ! assumes have cap_h_x

    call timer_lib_in('comm')
    call parallel_lib_r(transpose(cap_h_c),cap_h_v)
    call timer_lib_out('comm')

    call timer_lib_in('collision')

    ic_loc = 0
    do ic=nc1,nc2
       ic_loc = ic_loc+1

       ! Set-up the RHS: H = f + ze/T G phi

       cvec(:) = cap_h_v(ic_loc,:)

       ! This is a key loop for performance
       bvec = (0.0,0.0)
       do ivp=1,nv
          do iv=1,nv
             bvec(iv) = bvec(iv)+cmat(iv,ivp,ic_loc)*cvec(ivp)
          enddo
       enddo

       cap_h_v(ic_loc,:) = bvec(:)

    enddo

    call timer_lib_out('collision')

    ! Compute the new phi
    call cgyro_field_v

    call timer_lib_in('comm')
    call parallel_lib_f(cap_h_v,cap_h_ct)
    cap_h_c = transpose(cap_h_ct)
    call timer_lib_out('comm')

    call timer_lib_in('collision')
    ! Compute the new h_x
    iv_loc = 0
    do iv=nv1,nv2

       iv_loc = iv_loc+1

       is = is_v(iv)
       ix = ix_v(iv)
       ie = ie_v(iv)

       do ic=1,nc

          ir = ir_c(ic) 
          it = it_c(ic)

          h_x(ic,iv_loc) = cap_h_c(ic,iv_loc) &
               -z(is)/temp(is)*j0_c(ic,iv_loc)*field(ir,it,1)

          if (n_field > 1) then
             h_x(ic,iv_loc) = h_x(ic,iv_loc) &
                  +z(is)/temp(is)*j0_c(ic,iv_loc)*field(ir,it,2) &
                  *xi(ix)*sqrt(2.0*energy(ie))*vth(is)
          endif

       enddo
    enddo

    call timer_lib_out('collision')

  end subroutine cgyro_step_collision
  
end module cgyro_collision
