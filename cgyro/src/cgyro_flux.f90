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

  integer :: ie,ix,is,it,ir,l
  real :: dv
  real :: c_n,c_n0
  real :: c_t,c_t0

  flux_loc(:,:,:) = 0.0
  moment_loc(:,:,:) = 0.0
  gflux_loc(:,:) = 0.0

  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     ! Integration weight
     dv  = w_xi(ix)*w_e(ie)

     ! Density moment weight
     c_n = dv*dens(is)

     ! Energy moment weight
     c_t = dv*dens(is)*temp(is)*energy(ie)

     ! Adiabatic coefficient
     c_n0 = z(is)*dens(is)/temp(is)
     c_t0 = 1.5*temp(is)*c_n0

     do ic=1,nc

        ir = ir_c(ic)
        it = it_c(ic)

        ! Density flux: Gamma_a
        flux_loc(ir,is,1) = flux_loc(ir,is,1) &
             -c_n*aimag(cap_h_c(ic,iv_loc)*conjg(psi(ic,iv_loc)))*w_theta(it)

        ! Energy flux : Q_a
        flux_loc(ir,is,2) = flux_loc(ir,is,2) &
             -c_t*aimag(cap_h_c(ic,iv_loc)*conjg(psi(ic,iv_loc)))*w_theta(it)

        if (it == it0) then
           ! Density moment: (delta n_a)/(n_norm rho_norm)
           moment_loc(ir,is,1) = moment_loc(ir,is,1)-c_n0*field(1,ic) &
                +c_n*cap_h_c(ic,iv_loc)*jvec_c(1,ic,iv_loc)
           ! Energy moment : (delta E_a)/(n_norm T_norm rho_norm)
           moment_loc(ir,is,2) = moment_loc(ir,is,2)-c_t0*field(1,ic) &
                +c_t*cap_h_c(ic,iv_loc)*jvec_c(1,ic,iv_loc)
        endif

        ! Global fluxes (complex)
        do l=0,4
           if (ir-l > 0) then
              gflux_loc(l,is) = gflux_loc(l,is) &
                   +c_t*i_c*cap_h_c(ic,iv_loc)*conjg(psi(ic_c(ir-l,it),iv_loc))*w_theta(it)
           endif
        enddo
 
     enddo

  enddo

  ! Complete definition of fluxes
  flux_loc  =  flux_loc*(2*k_theta*rho)
  gflux_loc = gflux_loc*(2*k_theta*rho)

  ! GyroBohm normalizations
  flux_loc   =  flux_loc/rho**2
  gflux_loc  = gflux_loc/rho**2
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

  call MPI_ALLREDUCE(moment_loc(:,:,:), &
       moment(:,:,:), &
       size(moment), &
       MPI_DOUBLE_COMPLEX, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)

  ! Global complex gflux(l,ky), below, is still distributed over n 

  call MPI_ALLREDUCE(gflux_loc(:,:), &
       gflux(:,:), &
       size(gflux), &
       MPI_DOUBLE_COMPLEX, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)

end subroutine cgyro_flux
