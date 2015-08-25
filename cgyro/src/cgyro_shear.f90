subroutine cgyro_shear

  use cgyro_globals
  use timer_lib
  
  implicit none

  integer :: ir,it

  call timer_lib_in('stream')

  if (i_time == 1) gtime = 0.0

  gtime = gtime+omega_eb*delta_t

  if (gtime >= 0.5) then

     gtime = gtime-1.0

     do ir=2,n_radial
        do it=1,n_theta
           h_x(ic_c(ir-1,it),:) = h_x(ic_c(ir,it),:)
           cap_h_c(ic_c(ir-1,it),:) = cap_h_c(ic_c(ir,it),:)
           psi(ic_c(ir-1,it),:) = psi(ic_c(ir,it),:)
        enddo
     enddo
     h_x(ic_c(n_radial,:),:) = 0.0
     cap_h_c(ic_c(n_radial,:),:) = 0.0
     psi(ic_c(n_radial,:),:) = 0.0

     do ir=2,n_radial
        field(ir-1,:,:) = field(ir,:,:)
     enddo
     field(n_radial,:,:) = 0.0

  endif
  call timer_lib_out('stream')

end subroutine cgyro_shear
