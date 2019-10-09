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

  integer :: p,k,ir,j
  integer :: icm,icp
  real :: nu_eff
  real :: tau_ave = 100.0

  if (nonlinear_flag == 0) return

  if (profile_shear_flag == 1 .or. shear_method == 2) then

     nu_eff = nu_global*(abs(gamma_e)+maxval(abs(sdlnndr(:)))+maxval(abs(sdlntdr(:))))
     sa = 1.0+exp(-delta_t/tau_ave)*sa

     ! Time-delay source
     if (n == 0) then

        ir = 1+n_radial/2

        icm = (ir-1-1)*n_theta
        icp = (ir-1+1)*n_theta

        do j=1,n_theta
           ! Recursive update of p=+1 source 
           source(j,:) = source(j,:)+(h_x(icp+j,:)-source(j,:))/sa
           ! Subtract source from h(0,+1)
           h_x(icp+j,:) = h_x(icp+j,:)-nu_eff*delta_t*source(j,:)
           ! Subtract source from h(0,-1)
           h_x(icm+j,:) = h_x(icm+j,:)-nu_eff*delta_t*conjg(source(j,:))
        enddo

     endif
  endif

end subroutine cgyro_source
