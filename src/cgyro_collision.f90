module cgyro_collision

  implicit none

  public :: COLLISION_alloc, cgyro_step_collision
  logical, private :: initialized = .false.

  real, dimension(:,:,:), allocatable, private :: cmat
  complex, dimension(:), allocatable, private :: cvec,bvec

contains

  subroutine COLLISION_alloc(flag)

    use cgyro_globals
    use cgyro_equilibrium, only : omega_trap, k_perp

    implicit none

    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    real, dimension(:,:,:), allocatable :: nu_d, nu_s
    real, dimension(:,:), allocatable :: rs
    real, external :: derf
    real :: xa, xb, tauinv_ab
    real :: sum_nu, sum_nud, sum_den
    integer :: jv
    integer :: is,ir,it,ix,ie,js,je,jx,ks
    ! parameters for matrix solve
    real, dimension(:,:), allocatable :: amat, bmat

    if (collision_model == 0) return

    if (flag == 1) then

       if(initialized) return

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

                end select

             enddo
          enddo
       enddo

       do is=1,n_species
          do js=1,n_species
             rs(is,js) = 0.0
             do ie=1,n_energy
                rs(is,js) = rs(is,js)+w_e(ie)*nu_s(ie,is,js)*energy(ie)
             enddo
          enddo
       enddo

       allocate(cmat(nv,nv,nc_loc))
       allocate(cvec(nv))
       allocate(bvec(nv))

       ! set-up the collision matrix
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

       ! matrix solve parameters
       allocate(work(nv))
       allocate(i_piv(nv))
       allocate(amat(nv,nv))
       allocate(bmat(nv,nv))

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

             sum_nud = 0.0
             do ks=1,n_species
                sum_nud = sum_nud + nu_d(ie,is,ks)
             enddo

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

                ! Collision component: Lorentz
                if (is==js .and. ie==je) then
                   cmat(iv,jv,ic_loc) = cmat(iv,jv,ic_loc) &
                        - (0.5*delta_t) * xi_lor_mat(ix,jx) &
                        *0.5*sum_nud

                   amat(iv,jv) = amat(iv,jv) &
                        + (0.5*delta_t) * xi_lor_mat(ix,jx) &
                        *0.5*sum_nud
                endif

                ! Collision component: Restoring
                if (abs(rs(is,js)) > epsilon(0.)) then
                   cmat(iv,jv,ic_loc) = cmat(iv,jv,ic_loc) &
                        - (0.5*delta_t) &
                        * 1.5 * (mass(js)/mass(is)) &
                        * (dens(js)/dens(is)) &
                        * (vth(js)/vth(is)) &
                        * nu_s(ie,is,js) * sqrt(energy(ie)) &
                        * xi(ix) &
                        * nu_s(je,js,is) * sqrt(energy(je)) &
                        * xi(jx) * w_e(je) * w_xi(jx) &
                        / rs(is,js)
                   amat(iv,jv) = amat(iv,jv) &
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

                if (collision_model == 3) then
                   if (is==js .and. ie==je) then
                      sum_nu = 0.0
                      do ks=1,n_species
                         sum_nu = sum_nu &
                              + (nu_d(ie,is,ks)-nu_s(ie,is,ks))
                      enddo
                      cmat(iv,jv,ic_loc) = cmat(iv,jv,ic_loc) &
                           - (0.5*delta_t) * sum_nu &
                           * 1.5 * xi(ix) * xi(jx) * w_xi(jx)
                      amat(iv,jv) = amat(iv,jv) &
                           + (0.5*delta_t) * sum_nu &
                           * 1.5 * xi(ix) * xi(jx) * w_xi(jx)
                   endif
                endif

                ! Poisson component 
                if (zf_test_flag == 1 .and. ae_flag == 1) then
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
