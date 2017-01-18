!---------------------------------------------------------
! cgyro_shear_advect.f90
!
! NOTE:
!                       k_theta*length*gamma_e
!           omega_eb  = ---------------------- 
!                               2 pi
!---------------------------------------------------------

subroutine cgyro_shear_advect

  use cgyro_globals

  implicit none

  integer :: j,ir
  real :: eps
  complex, dimension(n_radial,n_radial) :: a
  complex, dimension(n_radial) :: h0
  complex, dimension(n_radial) :: s0,cwork
  integer, dimension(n_radial) :: ipiv

  eps = 0.25*omega_eb*delta_t

  do iv_loc=1,nv_loc
     do j=1,n_theta

        ! Solve A x = s

        h0(:) = h_x(ic_c(:,j),iv_loc)
        s0(:) = 0.0      
        a(:,:) = 0.0
        do ir=1,n_radial
           a(ir,ir) = 1.0
           if (ir > 1) a(ir,ir-1) = eps
           if (ir < n_radial) a(ir,ir+1) = -eps
           if (ir > 1) s0(ir) = s0(ir)-eps*h0(ir-1)
           if (ir < n_radial) s0(ir) = s0(ir)+eps*h0(ir+1)
        enddo

        call ZSYSV('U',n_radial,1,a,n_radial,ipiv,s0,n_radial,cwork,n_radial,info)

        h_x(ic_c(:,j),iv_loc) = s0(:)

     enddo
  enddo

  call cgyro_field_c

end subroutine cgyro_shear_advect
