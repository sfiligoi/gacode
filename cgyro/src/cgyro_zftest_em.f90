!-----------------------------------------------------------------
! cgyro_zftest_em.f90
!
! PURPOSE:
! initializes h and evaluates relevant integrals 
! for an electromagnetic zonal flow test 
! Corresp
!-----------------------------------------------------------------

subroutine cgyro_zftest_em

  use mpi
  use cgyro_globals
  use cgyro_io

  implicit none

  integer :: ir,it,is,ie,ix,itor
  real :: arg
  real, dimension(n_species,n_theta,nt1:nt2) :: ansum,adsum,alphah,sum_loc

  sum_loc(:,:,:) = 0.0

  ! Calculating the species and theta dependent function, alphah, that ensures
  !   no B|| fluctuation at t=0. It is a fraction of two integrals, 
  !   ansum and adsum, which are calculated first.  

!$omp parallel private(iv_loc,is,ix,ie,ic,it,arg)
!$omp do collapse(2) reduction(+:sum_loc)
  do itor=nt1,nt2
   do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     ix = ix_v(iv) 
     ie = ie_v(iv)
     do ic=1,nc
        it = it_c(ic)

        arg = k_perp(ic,itor)*rho*vth(is)*mass(is)/(z(is)*bmag(it)) &
             *vel2(ie)*sqrt(1.0-xi(ix)**2)

        sum_loc(is,it,itor) = sum_loc(is,it,itor) & 
             + w_exi(ie,ix)*energy(ie)*(1.0-xi(ix)**2) &
             *(2.0*(bessel_j1(arg)/arg)*(2.0-bessel_j0(arg)) -1.0)
     enddo
   enddo
  enddo
!$omp end do
!$omp end parallel

  call MPI_ALLREDUCE(sum_loc,&
       ansum,&
       size(ansum),&
       MPI_DOUBLE_PRECISION,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  sum_loc(:,:,:) = 0.0

!$omp parallel private(iv_loc,is,ix,ie,ic,it,arg)
!$omp do collapse(2) reduction(+:sum_loc)
  do itor=nt1,nt2
   do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     ix = ix_v(iv) 
     ie = ie_v(iv)
     do ic=1,nc
        it = it_c(ic)

        arg = k_perp(ic,itor)*rho*vth(is)*mass(is)/(z(is)*bmag(it)) &
             *vel2(ie)*sqrt(1.0-xi(ix)**2)

        sum_loc(is,it,itor) = sum_loc(is,it,itor) & 
             + w_exi(ie,ix)*energy(ie)*(1.0-xi(ix)**2) &
             *2.0*(bessel_j1(arg)/arg)*(2.0-bessel_j0(arg))*(energy(ie)-1.5)
     enddo
   enddo
  enddo
!$omp end do
!$omp end parallel

  call MPI_ALLREDUCE(sum_loc,&
       adsum,&
       size(adsum),&
       MPI_DOUBLE_PRECISION,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  alphah(:,:,:) = ansum(:,:,:) / adsum(:,:,:) 

  ! Constructing the h_x initial condition
!$omp parallel do private(ic,ir,it,iv,iv_loc,is,ix,ie,arg) shared(h_x)
  do itor=nt1,nt2
   do ic=1,nc

     ir = ir_c(ic) 
     it = it_c(ic)

     do iv=nv1,nv2

        iv_loc = iv-nv1+1
        is = is_v(iv)
        ix = ix_v(iv)       
        ie = ie_v(iv)

        if (px(ir) /= 0) then
           arg = k_perp(ic,itor)*rho*vth(is)*mass(is)/(z(is)*bmag(it)) &
                *vel2(ie)*sqrt(1.0-xi(ix)**2)           


           h_x(ic,iv_loc,itor) = z(is)/temp(is) * (2.0 - bessel_j0(arg)) &
                *(1.0 - alphah(is,it,itor)*(energy(ie) - 1.5)) &
                - z(is)/temp(is) * bessel_j0(arg)

        endif
     enddo
   enddo
  enddo

end subroutine cgyro_zftest_em
