!-----------------------------------------------------------------
! cgyro_upwind.f90
!
! PURPOSE:
!  Compute corrected distribution used in conservative advection scheme:
!
!                  /
!              J0  | dv J0 |vp| g 
!                  /
!  |vp| g -  ----------------------
!                   /
!                   | dv J0^2 
!                   /
!
!-----------------------------------------------------------------

subroutine cgyro_upwind

  use cgyro_globals
  use mpi

  implicit none

  integer :: is,ie,ix
  complex :: fac
  complex, dimension(n_species,nc) :: res_loc 
  complex, dimension(n_species,nc) :: res

!$omp workshare
  res_loc(:,:) = (0.0,0.0)
!$omp end workshare

!$omp parallel private(ic,iv_loc,is,ix,ie,fac)
!$omp do reduction(+:res_loc)
  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc
        fac = w_e(ie)*w_xi(ix)*abs(xi(ix))*sqrt(energy(ie))*g_x(ic,iv_loc)
        res_loc(is,ic) = res_loc(is,ic)+jvec_c(1,ic,iv_loc)*fac
     enddo
  enddo
!$omp end do
!$omp end parallel

  call MPI_ALLREDUCE(res_loc(:,:),&
       res(:,:),&
       size(res(:,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  do iv=nv1,nv2
     do ic=1,nc

        iv_loc = iv-nv1+1
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)

        g_x(ic,iv_loc) = abs(xi(ix))*sqrt(energy(ie))*g_x(ic,iv_loc) &
             -jvec_c(1,ic,iv_loc)*res(is,ic)/res_norm(is,ic)

     enddo
  enddo

end subroutine cgyro_upwind
