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

  implicit none

  integer :: ir,it,is,ie,ix
  real :: a
  complex, dimension(n_theta,nv_loc) :: a1

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

        enddo
     enddo

     call cgyro_field_c

  endif

end subroutine cgyro_shear_pt
