!---------------------------------------------------------
! cgyro_advect_wavenumber.f90
!
! PURPOSE:
!  Manage shearing by wavenumber advection.
!---------------------------------------------------------

subroutine cgyro_advect_wavenumber(ij)

  use cgyro_globals
  use timer_lib

  implicit none

  integer, intent(in) :: ij
  integer :: ir,ip,j
  complex, dimension(1-nup_wave:n_radial+nup_wave) :: h0
  complex :: dh

  call timer_lib_in('shear')

  ! Zero work array including zero boundary regions
  h0 = 0.0

  ! Wavenumber advection ExB shear
  if (shear_method == 2) then
     do iv=nv1,nv2
        iv_loc = iv-nv1+1
        do j=1,n_theta
           h0(1:n_radial) = h_x(ic_c(:,j),iv_loc)
           do ir=1,n_radial
              dh = 0.0
              do ip=-nup_wave,nup_wave
                 dh = dh+der_wave(ip)*h0(ir+ip)
              enddo
              rhs(ic_c(ir,j),iv_loc,ij) = rhs(ic_c(ir,j),iv_loc,ij)+ &
                   omega_eb*dh
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
              dh = 0.0
              do ip=-nup_wave,nup_wave
                 dh = dh+der_wave(ip)*h0(ir+ip)
              enddo
              rhs(ic_c(ir,j),iv_loc,ij) = rhs(ic_c(ir,j),iv_loc,ij)+dh
           enddo
        enddo
     enddo
  endif

  call timer_lib_out('shear')

end subroutine cgyro_advect_wavenumber
