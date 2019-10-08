!---------------------------------------------------------
! cgyro_source.f90
!
! PURPOSE:
!  Time-delay source
!---------------------------------------------------------

subroutine cgyro_source

  use cgyro_globals
  use timer_lib

  implicit none

  integer :: p,k,ir,j,icc,in
  complex, dimension(:,:),allocatable :: he
  real :: nu_eff
  
  if (nonlinear_flag == 0) return

  if (profile_shear_flag == 1 .or. shear_method == 2) then

     nu_eff = nu_global*(abs(gamma_e)+abs(maxval(sdlnndr(:)))+abs(maxval(sdlntdr(:))))
     sa = 1.0+exp(-delta_t*nu_eff)*sa

     ! Time-delay source
     if (n == 0) then

        ir = 1+n_radial/2

        do p=-1,1,2
           icc = (ir-1+p)*n_theta
           k = (p+1)/2
           do j=1,n_theta
              ha(k,j,:) = ha(k,j,:)+(h_x(icc+j,:)-ha(k,j,:))/sa
              h_x(icc+j,:) = h_x(icc+j,:)-nu_eff*delta_t*ha(k,j,:)  
           enddo
        enddo

     endif
  endif

end subroutine cgyro_source
