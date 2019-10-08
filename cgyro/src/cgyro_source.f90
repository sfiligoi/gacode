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

  if (nonlinear_flag == 0) return

  sa = 1.0+exp(-delta_t*nu_global)*sa

  ! Time-delay source
  if (n == 0) then

     ir = 1+n_radial/2

     do p=-1,1,2
        icc = (ir-1+p)*n_theta
        k = (p+1)/2
        do j=1,n_theta
           ha(k,j,:) = (h_x(icc+j,:) + ha(k,j,:)*(sa-1))/sa
           h_x(icc+j,:) = h_x(icc+j,:)-ha(k,j,:)*delta_t
        enddo
     enddo

  endif

end subroutine cgyro_source
