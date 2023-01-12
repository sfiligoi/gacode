subroutine cgyro_rhs_comm_async(which)
  use cgyro_nl_comm
  use cgyro_globals

  implicit none

  integer, intent(in) :: which

  if (nonlinear_flag == 1) then
     if (which == 1) then
       call cgyro_nl_fftw_comm1_async
     else
       call cgyro_nl_fftw_comm2_async
     endif
  endif

end subroutine cgyro_rhs_comm_async

! Note: Calling test propagates the async operations in some MPI implementations
subroutine cgyro_rhs_comm_test(which)
  use cgyro_nl_comm
  use cgyro_globals

  implicit none

  integer, intent(in) :: which

  if (nonlinear_flag == 1) then
     if (which == 1) then
       call cgyro_nl_fftw_comm1_test
     else
       call cgyro_nl_fftw_comm2_test
     endif
  endif

end subroutine cgyro_rhs_comm_test

subroutine cgyro_rhs(ij)

  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer, intent(in) :: ij
  integer :: is, ix, ie, js, jx, je, jv, it, j, k, itor
  integer :: id,jc
  real :: rval,rval2
  complex :: rhs_stream,rhs_el
  complex, dimension(:,:), allocatable :: rhs_trap
  complex, dimension(:), allocatable   :: bvec_trap
  integer :: nj_loc

  ! h_x is not modified after this and before nl_fftw
  call cgyro_rhs_comm_async(1)

  ! Prepare suitable distribution (g, not h) for conservative upwind method
  if (n_field > 1) then
     call timer_lib_in('str')

!$omp parallel do collapse(2) private(iv_loc,is,ic)
     do itor=nt1,nt2
      do iv=nv1,nv2
        iv_loc = iv-nv1+1
        is = is_v(iv)
        do ic=1,nc
           g_x(ic,iv_loc,itor) = h_x(ic,iv_loc,itor)+ & 
                (z(is)/temp(is))*jvec_c(2,ic,iv_loc,itor)*field(2,ic,itor)
        enddo
      enddo
     enddo
     call timer_lib_out('str')
  else
     call timer_lib_in('str_mem')

     g_x(:,:,:) = h_x(:,:,:)

     call timer_lib_out('str_mem')
  endif
  call cgyro_rhs_comm_test(1)

  ! Correct g_x for number conservation
  call cgyro_upwind

  call cgyro_rhs_comm_test(1)
  call cgyro_rhs_comm_async(2)

  call timer_lib_in('str')

!$omp parallel do collapse(3) &
!$omp& private(iv_loc,is,rval,rval2,rhs_stream,id,jc,rhs_el) 
  do itor=nt1,nt2
   do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        ! Diagonal terms
        rhs_el = &
             omega_cap_h(ic,iv_loc,itor)*cap_h_c(ic,iv_loc,itor)+&
             omega_h(ic,iv_loc,itor)*h_x(ic,iv_loc,itor)

        is = is_v(iv)
        ! Parallel streaming with upwind dissipation 
        rval  = omega_stream(it_c(ic),is,itor)*vel(ie_v(iv))*xi(ix_v(iv))
        rval2 = abs(omega_stream(it_c(ic),is,itor))

        rhs_stream = 0.0
        do id=-nup_theta,nup_theta
           jc = icd_c(id, ic,itor)
           rhs_stream = rhs_stream &
                -rval*dtheta(id,ic,itor)*cap_h_c(jc,iv_loc,itor)  &
                -rval2*dtheta_up(id,ic,itor)*g_x(jc,iv_loc,itor) 
        enddo

        rhs_el = rhs_el + &
             sum(omega_s(:,ic,iv_loc,itor)*field(:,ic,itor))

        rhs(ic,iv_loc,itor,ij) = rhs_el+rhs_stream

     enddo
   enddo
  enddo
  call cgyro_rhs_comm_test(1)

  ! Explicit trapping term
  if (explicit_trap_flag == 1) then

     allocate(rhs_trap(nc,nv_loc))
     allocate(bvec_trap(nv))
     call parallel_lib_rtrans_pack(cap_h_c)
     call parallel_lib_r_do(cap_h_v)
     call parallel_lib_nj_loc(nj_loc)
     
     do itor=nt1,nt2
      do ic=nc1,nc2
        ic_loc = ic-nc1+1
        it = it_c(ic)
        
        bvec_trap(:) = (0.0,0.0)
        do iv=1,nv
           is = is_v(iv)
           ix = ix_v(iv)
           ie = ie_v(iv)
           do jv=1,nv
              js = is_v(jv)
              jx = ix_v(jv)
              je = ie_v(jv)
              if (is == js) then
                 if (ie == je) then
                    bvec_trap(iv) = bvec_trap(iv) - (omega_trap(it,is,itor) * vel(ie) &
                         + omega_rot_trap(it,is) / vel(ie)) &
                         * (1.0 - xi(ix)**2) * xi_deriv_mat(ix,jx) * cap_h_v(ic_loc,itor,jv)
                 endif
                 if (ix == jx) then
                    bvec_trap(iv) = bvec_trap(iv) - omega_rot_u(it,is) * xi(ix) &
                         * e_deriv1_rot_mat(ie,je)/sqrt(1.0*e_max) * cap_h_v(ic_loc,itor,jv)
                 endif
              endif
           enddo
        enddo

        do k=1,nproc
           do j=1,nj_loc
              fsendf(j,itor,ic_loc,k) = bvec_trap(j+(k-1)*nj_loc)
           enddo
        enddo
      enddo
     enddo
     
     call parallel_lib_f_i_do(cap_h_ct)
     do itor=nt1,nt2
      do iv=nv1,nv2
        iv_loc = iv-nv1+1
        do ic=1,nc
           rhs_trap(ic,iv_loc) = cap_h_ct(iv_loc,itor,ic)
        enddo
      enddo

      rhs(:,:,itor,ij) = rhs(:,:,itor,ij) +  rhs_trap(:,:)
     enddo
     
     deallocate(rhs_trap)
     deallocate(bvec_trap)
     
  endif
  
  call timer_lib_out('str')
  call cgyro_rhs_comm_test(1)

  ! Wavenumber advection shear terms
  call cgyro_advect_wavenumber(ij)

  ! Nonlinear evaluation [f,g]
  if (nonlinear_flag == 1) then     
        call cgyro_nl_fftw(ij)
  endif

end subroutine cgyro_rhs
