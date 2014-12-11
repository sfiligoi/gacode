subroutine cgyro_flux

  use mpi
  use cgyro_globals
  use cgyro_equilibrium

  implicit none

  integer :: ie,ix,is,it,ir
  real :: c_n,c_t

  flux_loc(:,:) = 0.0

  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     c_n = dens(is)*                    (w_e(ie)*0.5*w_xi(ix))*k_theta*rho
     c_t = dens(is)*temp(is)*energy(ie)*(w_e(ie)*0.5*w_xi(ix))*k_theta*rho

     do ic=1,nc

        it = it_c(ic)

        flux_loc(is,1) = flux_loc(is,1)-&
             c_n*aimag(2.0*cap_h_c(ic,iv_loc)*conjg(psi(ic,iv_loc)))*w_theta(it)

        flux_loc(is,2) = flux_loc(is,2)-&
             c_t*aimag(2.0*cap_h_c(ic,iv_loc)*conjg(psi(ic,iv_loc)))*w_theta(it)

     enddo
  enddo

  ! GyroBohm normalization
  flux_loc = flux_loc/rho**2

  ! Reduced flux, below, is still distributed over n

  call MPI_ALLREDUCE(flux_loc(:,:), &
       flux(:,:), &
       size(flux), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)


  ! Now capture flux-surface averaged field intensities
  ! No need for reduction since fields are not distributed within NEW_COMM_1.

  power(:,:) = 0.0
  do ic=1,nc

     ir = ir_c(ic)
     it = it_c(ic)

     power(ir,:) = power(ir,:)+w_theta(it)*abs(field(ir,it,:))**2

  enddo

end subroutine cgyro_flux

!==========================================================================================
