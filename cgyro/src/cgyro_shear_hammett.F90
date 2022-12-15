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
  integer :: n_changes
  integer, dimension(nt1:nt2) :: change_sign
  integer :: my_sign

  
  n_changes = 0

  call timer_lib_in('shear')
  gtime(my_toroidal) = gtime(my_toroidal)+omega_eb_base*my_toroidal*delta_t

  ! Forward shearing
  if (gtime(my_toroidal) > 0.5) then

     gtime(my_toroidal) = gtime(my_toroidal)-1.0
     n_changes = n_changes + 1
     change_sign(my_toroidal) = 1

  ! Backward shearing
  else if (gtime(my_toroidal) < -0.5) then

     gtime(my_toroidal) = gtime(my_toroidal)+1.0
     n_changes = n_changes + 1
     change_sign(my_toroidal) = -1

  else
     change_sign(my_toroidal) = 0
  endif

  if (n_changes > 0) then
#ifdef _OPENACC
!$acc parallel loop independent gang vector_length(n_theta) &
!$acc&              private(ir,my_sign) present(h_x,ic_c) copyin(change_sign)
#else
!$omp parallel do private(ir,my_sign)
#endif
     do iv_loc=1,nv_loc

      my_sign = change_sign(my_toroidal)
      if (my_sign == 1) then
       ! Forward shearing
!acc loop seq
       do ir=2,n_radial
         h_x(ic_c(ir-1,:),iv_loc,my_toroidal) = h_x(ic_c(ir,:),iv_loc,my_toroidal)
       enddo

       h_x(ic_c(n_radial,:),iv_loc,my_toroidal) = 0.0
      else if (my_sign == -1) then
       ! Backward shearing
!acc loop seq
       do ir=n_radial-1,1,-1
         h_x(ic_c(ir+1,:),iv_loc,my_toroidal) = h_x(ic_c(ir,:),iv_loc,my_toroidal)
       enddo

       h_x(ic_c(1,:),iv_loc,my_toroidal) = 0.0
      endif
     enddo
  endif

  call timer_lib_out('shear')

  if (n_changes > 0) then
     call cgyro_field_c
  endif

end subroutine cgyro_shear_hammett
