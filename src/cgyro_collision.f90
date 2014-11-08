module cgyro_collision
  
  implicit none
  
  public :: COLLISION_alloc, COLLISION_do
  logical, private :: initialized = .false.

  real, dimension(:,:,:), allocatable, private :: cmat
  real, dimension(:,:), allocatable, private :: cvec,bvec

contains

  subroutine COLLISION_alloc(flag)

    use cgyro_globals
    use cgyro_gyro
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
          call cgyro_error('ERROR: (CGYRO) collision_model=0 requires kinetic electrons')
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

       allocate(cmat(nv,nv,nc_loc))
       allocate(cvec(nv,2))
       allocate(bvec(nv,2))

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
       if(adiabatic_ele_model == 1) then
          sum_den = sum_den + dens_ele / temp_ele
       endif

       ! matrix solve parameters
       allocate(work(nv))
       allocate(i_piv(nv))
       allocate(amat(nv,nv))
       allocate(bmat(nv,nv))

       if (collision_model == 3) then

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

       ic_loc = 0
       do ic=nc1,nc2
          ic_loc = ic_loc+1

          it = it_c(ic)
          ir = ir_c(ic)

          cmat(:,:,ic_loc) = (0.0,0.0)
          amat(:,:)       = (0.0,0.0)

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
                if (collision_model == 3) then
                   ! EAB: NEED to recheck 0.5 with Lorentz
                   cmat(iv,jv,ic_loc) = cmat(iv,jv,ic_loc) &
                        - (0.5*delta_t) &
                        * (-mass(js)/mass(is)) &
                        * (dens(js)/dens(is)) &
                        * (vth(js)/vth(is)) &
                        * nu_d(ie,is,js) * sqrt(energy(ie)) &
                        * 0.5*vecout_xi(ix) &
                        * nu_d(je,js,is) * sqrt(energy(je)) &
                        * 0.5*vecin_xi(jx) * w_e(je)*w_xi(jx) &
                        / rs_lor(is,js)
                   amat(iv,jv) = amat(iv,jv) &
                        + (0.5*delta_t) &
                        * (-mass(js)/mass(is)) &
                        * (dens(js)/dens(is)) &
                        * (vth(js)/vth(is)) &
                        * nu_d(ie,is,js) * sqrt(energy(ie)) &
                        * 0.5*vecout_xi(ix) &
                        * nu_d(je,js,is) * sqrt(energy(je)) &
                        * 0.5*vecin_xi(jx) * w_e(je)*w_xi(jx) &
                        / rs_lor(is,js)

                else if (abs(rs(is,js)) > epsilon(0.)) then
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

                if (collision_model == 2) then
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
                if(toroidal_model == 2 .and. adiabatic_ele_model == 1) then
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
                        * gyrox_J0(is,ir,it,ie,ix) &
                        * z(js)*dens(js) &
                        * gyrox_J0(js,ir,it,je,jx) * w_e(je) &
                        * 0.5 * w_xi(jx)
                   amat(iv,jv) = amat(iv,jv) &
                        - z(is)/temp(is) / &
                        (k_perp(it,ir)**2 * lambda_debye**2 &
                        * dens_ele / temp_ele &
                        + sum_den) &
                        * gyrox_J0(is,ir,it,ie,ix) &
                        * z(js)*dens(js) &
                        * gyrox_J0(js,ir,it,je,jx) * w_e(je) &
                        * 0.5 * w_xi(jx)
                endif

                ! Ampere component
                if(n_field > 1) then
                   cmat(iv,jv,ic_loc) = cmat(iv,jv,ic_loc) &
                        - z(is)/temp(is) / &
                        (2.0*k_perp(it,ir)**2 * rho**2 / betae_unit & 
                        * dens_ele * temp_ele) &
                        * (-gyrox_J0(is,ir,it,ie,ix)) &
                        * z(js)*dens(js) &
                        * xi(ix) * sqrt(2.0*energy(ie)) &
                        * gyrox_J0(js,ir,it,je,jx) * w_e(je) &
                        * 0.5 * w_xi(jx) &
                        * xi(jx) * sqrt(2.0*energy(je))
                   amat(iv,jv) = amat(iv,jv) &
                        - z(is)/temp(is) / &
                        (2.0*k_perp(it,ir)**2 * rho**2 / betae_unit & 
                        * dens_ele * temp_ele) &
                        * (-gyrox_J0(is,ir,it,ie,ix)) &
                        * z(js)*dens(js) &
                        * xi(ix) * sqrt(2.0*energy(ie)) &
                        * gyrox_J0(js,ir,it,je,jx) * w_e(je) &
                        * 0.5 * w_xi(jx) &
                        * xi(jx) * sqrt(2.0*energy(je))
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
       initialized = .false.

    endif

  end subroutine COLLISION_alloc
  
  subroutine COLLISION_do

    use parallel_lib
    use timer_lib

    use cgyro_globals
    use cgyro_field
    use cgyro_gyro

    implicit none

    integer :: is,ir,it,ie,ix
    integer :: ivp

    if (collision_model == -1) return

    ! compute new collisional cap_H: H = h + ze/T G phi
    ! assumes have cap_h_x

    ! allocate(cap_h_c(nc,nv_loc))
    ! allocate(cap_h_ct(nv_loc,nc))
    ! allocate(cap_h_v(nc_loc,nv))

    call timer_lib_in('comm')
    call parallel_lib_r(transpose(cap_h_c),cap_h_v)
    call timer_lib_out('comm')


    call timer_lib_in('collision')

    ic_loc = 0
    do ic=nc1,nc2
       ic_loc = ic_loc+1

       it = it_c(ic)
       ir = ir_c(ic)

       ! Set-up the RHS: H = f + ze/T G phi
       ! 
       ! Pack Re(H) and Im(H) into two-column matrix for use with DGEMM
       do iv=1,nv
          cvec(iv,1) = real(cap_h_v(ic_loc,iv))
          cvec(iv,2) = imag(cap_h_v(ic_loc,iv))
       enddo

       ! Solve for H
       !call DGEMM('N','N',nv,2,nv,num1,cmat(ic_loc,:,:),&
       !      nv,cvec,nv,num0,bvec,nv)

       bvec = 0.0
       do iv=1,nv
          do ivp=1,nv
             bvec(iv,1) = bvec(iv,1)+cmat(iv,ivp,ic_loc)*cvec(ivp,1)
             bvec(iv,2) = bvec(iv,2)+cmat(iv,ivp,ic_loc)*cvec(ivp,2)
          enddo
       enddo


       do iv=1,nv
          cap_h_v(ic_loc,iv) = bvec(iv,1) + i_c * bvec(iv,2)
       enddo

    enddo

    ! Compute the new phi
    call FIELDh_do

    call timer_lib_out('collision')

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
               - z(is)/temp(is) * gyrox_J0(is,ir,it,ie,ix) &
               * field(ir,it,1)
          if(n_field > 1) then
             h_x(ic,iv_loc) = h_x(ic,iv_loc) &
                  + z(is)/temp(is) * gyrox_J0(is,ir,it,ie,ix) &
                  * field(ir,it,2) * xi(ix) * sqrt(2.0*energy(ie))
          endif
       enddo
    enddo

    call timer_lib_out('collision')

  end subroutine COLLISION_do
  
end module cgyro_collision
