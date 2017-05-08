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
  integer :: ir,l,ll,j
  complex, dimension(:),allocatable :: h0
  complex :: dh
  
  call timer_lib_in('shear')

  ! Wavenumber advection ExB shear
  if (shear_method == 2) then
     allocate(h0(1-2*n_wave:n_radial+2*n_wave))
     ! Zero wavenumber (probably causes boundary reflection)
     h0 = 0.0
!$omp parallel do private(j,h0,ir,dh,l,ll,ic,c_wave)
     do iv_loc=1,nv_loc
        do j=1,n_theta
           h0(1:n_radial) = h_x(ic_c(:,j),iv_loc)
           do ir=1,n_radial
              dh = 0.0
              do l=1,n_wave
                 ll = 2*l-1
                 dh = dh+c_wave(l)*(h0(ir+ll)-h0(ir-ll))
              enddo
              ic = ic_c(ir,j)
              rhs(ic,iv_loc,ij) = rhs(ic,iv_loc,ij)+omega_eb*dh
           enddo
        enddo
     enddo
     deallocate(h0)
  endif

  !-------------------------------------------------------------------------

  ! Wavenumber advection profile shear
  if (profile_shear_flag == 1) then
     allocate(h0(1-2*n_wave:n_radial+2*n_wave))
     ! Zero wavenumber (probably causes boundary reflection)
     h0 = 0.0
!$omp parallel do private(j,h0,ir,dh,l,ll,ic,c_wave)
     do iv_loc=1,nv_loc
        do j=1,n_theta
           do ir=1,n_radial
              h0(ir) = sum(omega_ss(:,ic_c(ir,j),iv_loc)*field(:,ic_c(ir,j)))
           enddo
           do ir=1,n_radial
              dh = 0.0
              do l=1,n_wave
                 ll = 2*l-1
                 dh = dh+c_wave(l)*(h0(ir+ll)-h0(ir-ll))
              enddo
              ic = ic_c(ir,j)
              rhs(ic,iv_loc,ij) = rhs(ic,iv_loc,ij)+dh
           enddo
        enddo
     enddo
     deallocate(h0)
  endif

  call timer_lib_out('shear')

end subroutine cgyro_advect_wavenumber
