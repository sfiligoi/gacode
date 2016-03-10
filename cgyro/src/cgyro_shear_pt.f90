!---------------------------------------------------------
! cgyro_shear_pt.f90
!
! PURPOSE:
!  Spectral shear algorithm.  Wavenumbers are shifted
!  to the left or right depending upon sign of omega_eb.
!
! NOTE:
!                       k_theta*length*gamma_e
!           omega_eb  = ---------------------- 
!                               2 pi
!---------------------------------------------------------

subroutine cgyro_shear_pt

  use cgyro_globals
  use timer_lib
  use GEO_interface

  implicit none

  integer :: ir,it,is,ie,ix
  real :: a,arg
  complex, dimension(n_theta,nv_loc) :: a1
  complex :: carg

  a     = omega_eb*delta_t
  gtime = gtime+a

  do iv=nv1,nv2
     do ic=1,nc

        iv_loc = iv-nv1+1
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
        ir = ir_c(ic) 
        it = it_c(ic)

        ! omega_rdrift
        omega_cap_h(ic,iv_loc) = omega_cap_h(ic,iv_loc) & 
             -omega_rdrift(it,is)*energy(ie)*&
             (1.0 + xi(ix)**2)*(2.0*pi*i_c*a/length) 

        k_perp(ic_c(ir,it)) = sqrt((2.0*pi*(px(ir)+gtime)*GEO_grad_r/length &
             + k_theta*GEO_gq*GEO_captheta)**2 &
             + (k_theta*GEO_gq)**2) 

        arg = k_perp(ic)*rho*vth(is)*mass(is)/(z(is)*bmag(it)) &
             *sqrt(2.0*energy(ie))*sqrt(1.0-xi(ix)**2)

        jvec_c(1,ic,iv_loc) = bessel_j0(abs(arg))

        carg = -i_c*k_theta*rho*(dlnndr(is)+dlntdr(is)*(energy(ie)-1.5)) &
             -i_c*k_theta*rho*(sqrt(2.0*energy(ie))*xi(ix)/vth(is) &
             *omega_gammap(it))

        omega_s(1,ic,iv_loc) = carg*jvec_c(1,ic,iv_loc)

     enddo
  enddo

  ! Forward shearing
  if (gtime > 1.0) then

     gtime = gtime-1.0

     a1 = h_x(ic_c(1,:),:)

     do ir=2,n_radial
        h_x(ic_c(ir-1,:),:) = h_x(ic_c(ir,:),:)
     enddo

     h_x(ic_c(n_radial,:),:) = a1*gamma_e_decay

     do iv=nv1,nv2
        do ic=1,nc

           iv_loc = iv-nv1+1
           is = is_v(iv)
           ix = ix_v(iv)
           ie = ie_v(iv)
           ir = ir_c(ic) 
           it = it_c(ic)

           ! omega_rdrift
           omega_cap_h(ic,iv_loc) = omega_cap_h(ic,iv_loc) & 
                -omega_rdrift(it,is)*energy(ie)*&
                (1.0 + xi(ix)**2)*(2.0*pi*i_c*(-1.0)/length) 

           k_perp(ic_c(ir,it)) = sqrt((2.0*pi*(px(ir)+gtime)*GEO_grad_r/length &
                + k_theta*GEO_gq*GEO_captheta)**2 &
                + (k_theta*GEO_gq)**2) 

           arg = k_perp(ic)*rho*vth(is)*mass(is)/(z(is)*bmag(it)) &
                *sqrt(2.0*energy(ie))*sqrt(1.0-xi(ix)**2)

           jvec_c(1,ic,iv_loc) = bessel_j0(abs(arg))

           carg = -i_c*k_theta*rho*(dlnndr(is)+dlntdr(is)*(energy(ie)-1.5)) &
                -i_c*k_theta*rho*(sqrt(2.0*energy(ie))*xi(ix)/vth(is) &
                *omega_gammap(it))

           omega_s(1,ic,iv_loc) = carg*jvec_c(1,ic,iv_loc)

        enddo
     enddo

     call cgyro_field_c

  endif

end subroutine cgyro_shear_pt
