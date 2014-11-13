subroutine cgyro_init_arrays

  use cgyro_globals
  use cgyro_equilibrium

  implicit none

  real, external :: BESJ0
  real :: arg
  integer :: is,ir,it,ie,ix

  do is=1,n_species
     do ir=1,n_radial
        do it=1,n_theta
           do ie=1,n_energy
              do ix=1,n_xi
                 arg = k_perp(it,ir) * rho * vth(is) &
                      * mass(is) / (z(is) * Bmag(it)) &
                      * sqrt(2.0* energy(ie)) &
                      * sqrt(1-xi(ix)**2)
                 gyrox_J0(is,ir,it,ie,ix) = BESJ0(abs(arg))
              enddo
           enddo
        enddo
     enddo
  enddo

  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc

        ir = ir_c(ic) 
        it = it_c(ic)

        arg = k_perp(it,ir)*rho*vth(is)*mass(is)/(z(is)*Bmag(it)) &
             *sqrt(2.0* energy(ie))*sqrt(1.0-xi(ix)**2)

        j0_c(ic,iv_loc) = BESJ0(abs(arg))

     enddo
  enddo

  ic_loc = 0
  do ic=nc1,nc2
     ic_loc = ic_loc+1

     it = it_c(ic)
     ir = ir_c(ic)

     do iv=1,nv

        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)

        arg = k_perp(it,ir)*rho*vth(is)*mass(is)/(z(is)*Bmag(it)) &
             *sqrt(2.0* energy(ie))*sqrt(1.0-xi(ix)**2)

        j0_v(ic_loc,iv) = BESJ0(abs(arg))

     enddo
  enddo


end subroutine cgyro_init_arrays
