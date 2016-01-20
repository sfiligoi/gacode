!-----------------------------------------------------------------
! cgyro_flux.f90
!
! PURPOSE:
!  Compute flux and mean-square fluctuation amplitudes as
!  functions of (kx,ky).
!-----------------------------------------------------------------

subroutine cgyro_flux

  use mpi
  use cgyro_globals

  implicit none

  integer :: ie,ix,is,it,ir
  real :: dv,c_n,c_t

  flux_loc(:,:,:) = 0.0
  moment_loc(:,:) = 0.0

  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     dv  = w_xi(ix)*w_e(ie)
     c_n = dens(is)*                    dv*k_theta*rho
     c_t = dens(is)*temp(is)*energy(ie)*dv*k_theta*rho

     do ic=1,nc

        ir = ir_c(ic)
        it = it_c(ic)

        ! Density flux: Gamma_a
        flux_loc(ir,is,1) = flux_loc(ir,is,1)-&
             c_n*aimag(2.0*cap_h_c(ic,iv_loc)*conjg(psi(ic,iv_loc)))*w_theta(it)

        ! Energy flux: Q_a
        flux_loc(ir,is,2) = flux_loc(ir,is,2)-&
             c_t*aimag(2.0*cap_h_c(ic,iv_loc)*conjg(psi(ic,iv_loc)))*w_theta(it)

        ! Density moment: delta n_a
        if (it == it0) then
           moment_loc(ir,is) = moment_loc(ir,is)+dens(is)*dv*cap_h_c(ic,iv_loc)
        endif

     enddo

  enddo

  ! GyroBohm normalization
  flux_loc = flux_loc/rho**2
  moment_loc = moment_loc/rho
  
  ! Reduced real flux(kx,ky), below, is still distributed over n 

  call MPI_ALLREDUCE(flux_loc(:,:,:), &
       flux(:,:,:), &
       size(flux), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)

  ! Reduced complex moment(kx,ky), below, is still distributed over n 

  call MPI_ALLREDUCE(moment_loc(:,:), &
       moment(:,:), &
       size(moment), &
       MPI_DOUBLE_COMPLEX, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)

end subroutine cgyro_flux

!==========================================================================================
