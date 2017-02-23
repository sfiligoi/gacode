subroutine cgyro_advect_wavenumber(ij)

  use cgyro_globals

  implicit none

  integer, intent(in) :: ij
  integer :: ir,j
  complex, dimension(0:n_radial+1) :: h0

  ! Zero work array
  h0 = 0.0

  ! Wavenumber advection ExB shear
  if (shear_method == 2) then
     do iv=nv1,nv2
        iv_loc = iv-nv1+1
        do j=1,n_theta
           h0(1:n_radial) = h_x(ic_c(:,j),iv_loc)
           do ir=1,n_radial
              rhs(ic_c(ir,j),iv_loc,ij) = rhs(ic_c(ir,j),iv_loc,ij)+ &
                   omega_eb*0.5*(h0(ir+1)-h0(ir-1))
           enddo
        enddo
     enddo
  endif

  ! Wavenumber advection profile shear
  if (profile_shear_flag == 1) then
     do iv=nv1,nv2
        iv_loc = iv-nv1+1
        do j=1,n_theta
           do ir=1,n_radial
              h0(ir) = sum(omega_ss(:,ic_c(ir,j),iv_loc)*field(:,ic_c(ir,j)))
           enddo
           do ir=1,n_radial
              rhs(ic_c(ir,j),iv_loc,ij) = rhs(ic_c(ir,j),iv_loc,ij)+ &
                   0.5*(h0(ir+1)-h0(ir-1))
           enddo
        enddo
     enddo
  endif

end subroutine cgyro_advect_wavenumber
