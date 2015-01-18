subroutine cgyro_shear

  use cgyro_globals

  implicit none

  integer :: ir,it

  if (i_time == 1) gtime = 0.0

  gtime = gtime+omega_eb*delta_t

  if (gtime >= 1.0) then

     gtime = gtime-1.0

     do ir=2,n_radial
        do it=1,n_theta
           h_x(ic_c(ir-1,it),:) = h_x(ic_c(ir,it),:)
           cap_h_c(ic_c(ir-1,it),:) = cap_h_c(ic_c(ir,it),:)
           psi(ic_c(ir-1,it),:) = psi(ic_c(ir,it),:)
        enddo
     enddo

     do ir=2,n_radial
        field(ir-1,:,:) = field(ir,:,:)
     enddo
     field(n_radial,:,:) = 0.0

  endif

end subroutine cgyro_shear
