!-----------------------------------------------------------------
! cgyro_kernel.f90
!
! PURPOSE:
!  Subroutinized main cgyro program.  
!
! NOTES:
!  This can be called directly using the driver routine cgyro 
!  (in which case input data will read from input.dat) or called 
!  as a subroutine using cgyro_sub.
!-----------------------------------------------------------------

subroutine cgyro_kernel

  use timer_lib
  use mpi

  use cgyro_globals
  use cgyro_io
  use cgyro_equilibrium
  use cgyro_collision

  implicit none


  if (silent_flag == 0 .and. i_proc == 0) then
     open(unit=io_run,file=trim(path)//runfile,status='replace')
     close(io_run)
  endif

  ! Timer initialization
  call timer_lib_init('init_arrays')
  call timer_lib_init('field_v')
  call timer_lib_init('field_c')
  call timer_lib_init('rhs')
  call timer_lib_init('rhs_nl')
  call timer_lib_init('collision')
  call timer_lib_init('comm')
  call timer_lib_init('comm_nl')

  ! 1. MPI setup
  call cgyro_mpi_grid
  if (error_status > 0) goto 100

  ! 2. Profile setup
  call cgyro_make_profiles
  if (error_status > 0) goto 100

  ! 3. Parameter consistency checks
  call cgyro_check
  if (error_status > 0) goto 100

  ! Construct energy nodes and weights
  allocate(energy(n_energy))
  allocate(w_e(n_energy))
  allocate(e_deriv1_mat(n_energy,n_energy))
  allocate(e_deriv2_mat(n_energy,n_energy))
  call pseudo_maxwell(n_energy,e_max,energy,w_e,e_deriv1_mat,e_deriv2_mat)

  ! Construct xi (pitch-angle) nodes and weights
  allocate(xi(n_xi))
  allocate(w_xi(n_xi))
  allocate(xi_lor_mat(n_xi,n_xi))
  allocate(xi_deriv_mat(n_xi,n_xi))
  call pseudo_legendre(n_xi,xi,w_xi,xi_deriv_mat,xi_lor_mat)

  ! Allocate distribution function and field arrays
  allocate(j0_c(nc,nv_loc))
  allocate(j0_v(nc_loc,nv))
  allocate(h_x(nc,nv_loc))
  allocate(psi(nc,nv_loc))
  allocate(f_nl(nc,nsplit,n_toroidal))
  allocate(g_nl(nc,nsplit,n_toroidal))
  allocate(rhs(4,nc,nv_loc))
  allocate(h0_x(nc,nv_loc))
  allocate(cap_h_c(nc,nv_loc))
  allocate(cap_h_ct(nv_loc,nc))
  allocate(cap_h_v(nc_loc,nv))
  allocate(cap_h_v_prime(nc_loc,nv))
  allocate(omega_cap_h(nc,nv_loc))
  allocate(omega_h(nc,nv_loc))
  allocate(omega_s(n_field,nc,nv_loc))
  allocate(field(n_radial,n_theta,n_field))
  allocate(field_loc(n_radial,n_theta,n_field))
  allocate(field_old(n_radial,n_theta,n_field))
  allocate(field_old2(n_radial,n_theta,n_field))
  allocate(field_old3(n_radial,n_theta,n_field))
  allocate(field_est(n_radial,n_theta,n_field))
  allocate(f_balloon(n_radial/box_size,n_theta))
  allocate(recv_status(MPI_STATUS_SIZE))

  allocate(thcyc(1-n_theta:2*n_theta))
  allocate(rcyc(n_radial,n_theta,-2:2))
  allocate(pcyc(-n_radial/2-n_radial:n_radial/2-1+n_radial))
  allocate(dtheta(n_radial,n_theta,-2:2))
  allocate(dtheta_up(n_radial,n_theta,-2:2))

  call EQUIL_alloc(1)
  call EQUIL_do

  !4. Array initialization
  call cgyro_init_arrays

  call COLLISION_alloc(1)

  if (silent_flag == 0 .and. i_proc == 0) then

     open(unit=io_data,file=trim(path)//runfile_grids,status='replace')
     write(io_data,'(i4)') n_species
     write(io_data,'(i4)') n_radial
     write(io_data,'(i4)') n_theta
     write(io_data,'(i4)') n_energy
     write(io_data,'(i4)') n_xi
     write(io_data,'(i4)') box_size
     write(io_data,'(i4)') px(:)
     write(io_data,'(1pe12.5)') theta(:)
     write(io_data,'(1pe12.5)') energy(:)
     write(io_data,'(1pe12.5)') xi(:)
     write(io_data,'(1pe12.5)') transpose(theta_B(:,:))
     close(io_data)

  endif

  ! Time-stepping
  n_time = nint(max_time/delta_t)

  i_time = 0

  io_control = 1*(1-silent_flag)
  call cgyro_write_timedata
  io_control = 2*(1-silent_flag)

  do i_time=1,n_time

     ! Collisionless step: returns new h_x, cap_h_x, fields 
     call cgyro_step_gk

     ! Collision step: returns new h_x, cap_h_x, fields
     call cgyro_step_collision

     ! Error estimate
     call cgyro_error_estimate

     ! Print results
     call cgyro_write_timedata

     if (abs(signal) == 1) exit

  enddo

  ! Print final distribution
  if (n_toroidal == 1) then
     io_control = 1*(1-silent_flag)
     call write_distribution(trim(path)//runfile_hx,io_data)
     io_control = 2*(1-silent_flag)
     call write_distribution(trim(path)//runfile_hx,io_data)
  endif

  ! Print timers
  if (i_proc == 0) then
     print *
     print '(a)', 'Timing Summary'
     print '(a,1x,1pe11.4)',' init_arrays ',timer_lib_time('init_arrays')
     print '(a,1x,1pe11.4)',' field_v     ',timer_lib_time('field_v')
     print '(a,1x,1pe11.4)',' field_c     ',timer_lib_time('field_c')
     print '(a,1x,1pe11.4)',' rhs         ',timer_lib_time('rhs')
     print '(a,1x,1pe11.4)',' rhs_nl      ',timer_lib_time('rhs_nl')
     print '(a,1x,1pe11.4)',' collision   ',timer_lib_time('collision')
     print '(a,1x,1pe11.4)',' comm        ',timer_lib_time('comm')
     print '(a,1x,1pe11.4)',' comm_nl     ',timer_lib_time('comm_nl')
  endif

100 continue

  call EQUIL_alloc(0)
  call COLLISION_alloc(0)

  if(allocated(indx_xi))       deallocate(indx_xi)
  if(allocated(px))        deallocate(px)
  if(allocated(energy))        deallocate(energy)
  if(allocated(w_e))           deallocate(w_e)
  if(allocated(e_deriv1_mat))  deallocate(e_deriv1_mat)
  if(allocated(e_deriv2_mat))  deallocate(e_deriv2_mat)
  if(allocated(xi))            deallocate(xi)
  if(allocated(w_xi))          deallocate(w_xi)
  if(allocated(xi_lor_mat))    deallocate(xi_lor_mat)
  if(allocated(xi_deriv_mat))  deallocate(xi_deriv_mat)
  if(allocated(h_x))           deallocate(h_x)
  if(allocated(cap_h_c))       deallocate(cap_h_c)
  if(allocated(cap_h_v))       deallocate(cap_h_v)
  if(allocated(field))         deallocate(field)
  if(allocated(field_loc))     deallocate(field_loc)
  if(allocated(field_old))     deallocate(field_old)
  if(allocated(field_est))     deallocate(field_est)
  if(allocated(f_balloon))     deallocate(f_balloon)

  if (zf_test_flag == 2 .and. ae_flag == 1) then
     deallocate(hzf)
     deallocate(xzf)
     deallocate(pvec_in)
     deallocate(pvec_outr)
     deallocate(pvec_outi)
  endif

end subroutine cgyro_kernel

!==================================================================================
! Provide integration error estimate via quadratic interpolation.
!==================================================================================

subroutine cgyro_error_estimate

  use cgyro_globals

  implicit none

  if (i_time == 1) then

     field_old2 = 0.0

  else 

     ! Estimate of field via quadratic interpolation
     field_est = 3.0*field_old-3.0*field_old2+field_old3

     field_error = sum(abs(field-field_est))/sum(abs(field))

  endif

  field_old3 = field_old2
  field_old2 = field_old
  field_old  = field

end subroutine cgyro_error_estimate
