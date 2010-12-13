!---------------------------------------------------
! gyro_collision_kernel_ebelli.f90
!
! PURPOSE:
!  Take a collision step using the cross-species coupled
!  operator with the RBF method.  
!---------------------------------------------------

subroutine gyro_collision_kernel_ebelli

  use gyro_globals
  use gyro_collision_private
  use gyro_pointers
  use math_constants

  !--------------------------------------------------------------
  implicit none
  integer :: p, pp, js, ks, ia, ib, ja, jb
  complex, dimension(n_rbf,n_kinetic) :: fc_lor
  complex, dimension(n_rbf,n_kinetic,n_kinetic) :: fc_res
  complex, dimension(n_rbf,n_kinetic) :: fcp_tot, fcp_int
  complex, dimension(n_rbf) :: cphase
  complex, dimension(n_kinetic,n_kinetic,n_stack,n_lambda,n_ine_loc_1) :: f_int
  complex, dimension(n_kinetic,n_stack,n_lambda,n_ine_loc_1) :: h_C_int
  complex, dimension(n_kinetic,n_kinetic) :: f_temp
  !---------------------------------------------------------------

  !call proc_time(CPU_Ct_in)

  do ks=1,n_kinetic
     call energy_sum(h_C_all(ks,:,:,:), f_int(ks,:,:,:,:),ks)
  enddo

  p_ine_loc = 0
  do p_ine = 1+i_proc_1,n_ine_1,n_proc_1   
     p_ine_loc = p_ine_loc+1
     i  = ine_i(p_ine)
     ie = ine_e(p_ine)
     
     p = 0
     do k=1,n_lambda
        do m=1,n_stack
           p = p+1
           cphase(p) = exp(i_c*angp(i)*theta_t(i,k,m))
           do ks=1,n_kinetic
              fc_lor(p,ks) = h_C_all(ks,m,k,p_ine_loc) * cphase(p)
              do js=1,n_kinetic
                 fc_res(p,ks,js) = f_int(ks,js,m,k,p_ine_loc) * cphase(p)
              enddo
           enddo
        enddo
     enddo
     
     do p=1,n_rbf
        do ks = 1, n_kinetic
           fcp_int(p,ks)    = 0.0
           fcp_tot(p,ks)    = 0.0
           
           do pp=1,n_rbf  
              ! Lorentz term
              fcp_tot(p,ks) = fcp_tot(p,ks) &
                   + d_rbf_lorentz(pp,p,ks,p_ine_loc) * fc_lor(pp,ks)
              
              fcp_int(p,ks) = fcp_int(p,ks) &
                   + d_rbf_lorentz_int(pp,p,ks,p_ine_loc) * fc_lor(pp,ks)  
           enddo
        enddo
     enddo
              

     ! Restoring term Rs
     do p=1,n_rbf
        do js=1,n_kinetic 
           do ks = 1, n_kinetic
              do pp=1,n_rbf
                 fcp_tot(p,ks) = fcp_tot(p,ks) & 
                      + d_rbf_rs(pp,p,ks,js,p_ine_loc) * fc_res(pp,js,ks)
                 fcp_int(p,ks) = fcp_int(p,ks) &
                      + d_rbf_rs_int(pp,p,ks,js,p_ine_loc) * fc_res(pp,js,ks)
              enddo     
           enddo
        enddo
     enddo

     do ks = 1, n_kinetic
        p = 0
        do k=1,n_lambda
           do m=1,n_stack
              p = p+1
              h_C_all(ks,m,k,p_ine_loc) = fcp_tot(p,ks)/cphase(p)
              h_C_int(ks,m,k,p_ine_loc) = fcp_int(p,ks)/cphase(p)
           enddo
        enddo
     enddo

  enddo

  do ks=1,n_kinetic
     call energy_sum(h_C_int(ks,:,:,:), f_int(ks,:,:,:,:),ks)
  enddo

  p_ine_loc = 0
  do p_ine = 1+i_proc_1,n_ine_1,n_proc_1   
     p_ine_loc = p_ine_loc+1
     i  = ine_i(p_ine)
     ie = ine_e(p_ine)
     do k=1,n_lambda
        do m=1,n_stack

           ! solve for int d^3v sqrt(ene) * nu * xi * f  
           do ia=1,n_kinetic
              do ib=1,n_kinetic
                 p = indx_coll(ia,ib)
                 f_temp(ia,ib) = 0.0
                 do ja=1,n_kinetic
                    do jb=1,n_kinetic
                       pp = indx_coll(ja,jb)
                       f_temp(ia,ib) = f_temp(ia,ib) &
                            + d_rbf_velint(pp,p,p_ine_loc) &
                            * f_int(ja,jb,m,k,p_ine_loc)
                    enddo
                 enddo
              enddo
           enddo

           ! total solution
           do ks=1,n_kinetic
              do js=1,n_kinetic
                 h_C_all(ks,m,k,p_ine_loc) = h_C_all(ks,m,k,p_ine_loc) &
                      + rs_coll_const(i,ks,js) * 0.5 * dt * xi(i,k,m) &
                      * sqrt(energy(ie,1)) * nu_coll_d(i,ie,ks,js) &
                      * f_temp(js,ks)
              enddo
           enddo

        enddo
     enddo
  enddo

  !call proc_time(CPU_Ct_out)
  !if(i_proc == 0) then
  !   print *, 'time_total = ', CPU_Ct_out - CPU_Ct_in
  !endif

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[do_collision_advance done]'
  endif

end subroutine gyro_collision_kernel_ebelli


subroutine energy_sum(f_in,f_out,ks)
  use gyro_globals
  use gyro_pointers
  use gyro_collision_private
  implicit none
  complex, intent(in), dimension(n_stack,n_lambda,n_ine_loc_1) :: f_in
  complex, intent(out), dimension(n_kinetic,n_stack,n_lambda,n_ine_loc_1) &
       :: f_out
  integer, intent(in) :: ks
  integer :: n_i, n_j, n_k, n_d, js
  integer, external :: parallel_dim
  complex, dimension(:,:,:), allocatable :: f_tr
  complex, dimension(:,:,:), allocatable :: f_sum
  complex, dimension(:,:,:,:), allocatable :: f_tr_all
  real :: efac

  ! Permute p_ie,k -> p_ki,e
  n_i = n_x
  n_j = n_energy
  n_k = n_lambda
  n_d = parallel_dim(n_k*n_i,n_proc_1)

  allocate(f_tr(n_stack,n_energy,n_d))
  allocate(f_sum(n_stack,n_d,n_kinetic))
  allocate(f_tr_all(n_stack,n_energy,n_d,n_kinetic))

  ! f_tr(m,e,{ki}}
  call rTRANSP_INIT(n_i,n_j,n_k,NEW_COMM_1)
  do m=1,n_stack
     call rTRANSP_DO(f_in(m,:,:),f_tr(m,:,:))
  enddo
  call rTRANSP_CLEANUP

  f_sum(:,:,:) = 0.0

  p_ki_loc = 0
  do p_ki = 1+i_proc_1,n_ki_1,n_proc_1
     p_ki_loc = p_ki_loc+1
     k = ki_k(p_ki)
     i = ki_i(p_ki)

     do ie=1,n_energy
        efac = sqrt(energy(ie,1)) * w_energy(ie,1)
        do js=1,n_kinetic
           f_sum(:,p_ki_loc,js) = f_sum(:,p_ki_loc,js) + f_tr(:,ie,p_ki_loc) *efac &
                * nu_coll_d(i,ie,ks,js)
        enddo
     enddo
     do ie=1,n_energy
        do js=1,n_kinetic
           f_tr_all(:,ie,p_ki_loc,js) = f_sum(:,p_ki_loc,js)
        enddo
     enddo
  enddo

  ! Permute p_ki,e -> p_ie,k
  n_i = n_lambda
  n_j = n_x
  n_k = n_energy
  n_d = parallel_dim(n_j*n_k,n_proc_1)

  do js=1,n_kinetic
     call fTRANSP_INIT(n_i,n_j,n_k,NEW_COMM_1)
     do m=1,n_stack
        call fTRANSP_DO(f_tr_all(m,:,:,js),f_out(js,m,:,:))
     enddo
     call fTRANSP_CLEANUP
  enddo

  deallocate(f_tr)
  deallocate(f_sum)
  deallocate(f_tr_all)

end subroutine energy_sum
