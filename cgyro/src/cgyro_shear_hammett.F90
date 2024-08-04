!-------------------------------------------------------------------
! cgyro_shear_hammett.f90
!
! PURPOSE:
!  Spectral shear algorithm using the Hammett discrete shift method.
!  Wavenumbers are shifted to the left or right at selected time 
!  intervals depending upon sign of omega_eb.
!
! NOTE:
!                       k_theta*length*gamma_e
!           omega_eb  = ---------------------- 
!                               2 pi
!--------------------------------------------------------------------

subroutine cgyro_shear_hammett

  use cgyro_globals
  use timer_lib

  implicit none

  integer :: ir
  complex, dimension(n_theta) :: a1

  
  gtime = gtime+omega_eb*delta_t

  ! Forward shearing
  if (gtime > 0.5) then

     call timer_lib_in('shear')
     gtime = gtime-1.0

#ifdef _OPENACC
!$acc parallel loop independent gang private(a1,ir) present(h_x,ic_c)
#else
!$omp parallel do private(a1,ir)
#endif
     do iv_loc=1,nv_loc
       a1(:) = h_x(ic_c(1,:),iv_loc)

!acc loop seq
       do ir=2,n_radial
         h_x(ic_c(ir-1,:),iv_loc) = h_x(ic_c(ir,:),iv_loc)
       enddo

       h_x(ic_c(n_radial,:),iv_loc) = 0.0
     enddo

     call timer_lib_out('shear')

     call cgyro_field_c

  endif

  ! Backward shearing
  if (gtime < -0.5) then

     call timer_lib_in('shear')
     gtime = gtime+1.0

#ifdef _OPENACC
!$acc parallel loop independent gang private(a1,ir) present(h_x,ic_c)
#else
!$omp parallel do private(a1,ir)
#endif
     do iv_loc=1,nv_loc
       a1(:) = h_x(ic_c(n_radial,:),iv_loc) 

!acc loop seq
       do ir=n_radial-1,1,-1
         h_x(ic_c(ir+1,:),iv_loc) = h_x(ic_c(ir,:),iv_loc)
       enddo

       h_x(ic_c(1,:),iv_loc) = 0.0
     enddo
     call timer_lib_out('shear')

     call cgyro_field_c

  endif

end subroutine cgyro_shear_hammett
