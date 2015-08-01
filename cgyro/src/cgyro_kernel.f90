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
  use GEO_interface

  implicit none

  integer :: ir, it

  ! Need to initialize the info and error runfiles very early
  if (silent_flag == 0 .and. i_proc == 0) then
     open(unit=io,file=trim(path)//runfile_info,status='replace')
     open(unit=io,file=trim(path)//runfile_err,status='replace')
  endif

  ! Timer initialization
  call timer_lib_init('init_arrays')
  call timer_lib_init('field_v')
  call timer_lib_init('field_c')
  call timer_lib_init('rhs')
  call timer_lib_init('rhs_impgk')
  call timer_lib_init('rhs_impphi')
  call timer_lib_init('rhs_nl')
  call timer_lib_init('coll_set1')
  call timer_lib_init('coll_set2')
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

  if (test_flag  == 0) then

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
     allocate(xi_upderiv_mat(n_xi,n_xi))
     call pseudo_legendre(n_xi,xi,w_xi,xi_deriv_mat,xi_lor_mat,xi_upderiv_mat)

     ! Allocate distribution function and field arrays
     allocate(j0_c(nc,nv_loc))
     allocate(j0_v(nc_loc,nv))
     allocate(h_x(nc,nv_loc))
     allocate(h_xs(nc,nv_loc))
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
     allocate(f_balloon(n_radial/box_size,n_theta))
     allocate(    flux(n_species,2))
     allocate(flux_loc(n_species,2))
     allocate(power(n_radial,n_field))
     allocate(recv_status(MPI_STATUS_SIZE))

     allocate(thcyc(1-n_theta:2*n_theta))
     allocate(rcyc(n_radial,n_theta,-nup:nup))
     allocate(dtheta(n_radial,n_theta,-nup:nup))
     allocate(dtheta_up(n_radial,n_theta,-nup:nup))


     ! Equilibrium set-up
     allocate(theta(n_theta))
     allocate(thetab(n_radial/box_size,n_theta))
     allocate(w_theta(n_theta))
     allocate(Bmag(n_theta))
     allocate(k_perp(n_theta,n_radial))
     allocate(omega_stream(n_theta,n_species))
     allocate(omega_trap(n_theta,n_species))
     allocate(omega_rdrift(n_theta,n_species))
     allocate(omega_adrift(n_theta,n_species))
     allocate(omega_aprdrift(n_theta,n_species))
     allocate(omega_cdrift(n_theta,n_species))
     allocate(omega_gammap(n_theta))

     d_theta = (2*pi/n_theta)
     do it=1,n_theta
        theta(it) = -pi+(it-1)*d_theta
     enddo

     do ir=1,n_radial/box_size
        do it=1,n_theta
           thetab(ir,it) = theta(it)+2*pi*(ir-1-n_radial/2/box_size)
        enddo
     enddo

     GEO_model_in = geo_numeq_flag
     GEO_ntheta_in   = geo_ntheta
     GEO_nfourier_in = geo_ny
     call GEO_alloc(1)

     call cgyro_equilibrium

     ! 4. Array initialization
     call cgyro_init_arrays

     call cgyro_init_implicit_gk

     if (collision_model /= 0) then
        allocate(cmat(nv,nv,nc_loc))
        allocate(cvec(nv))
        allocate(bvec(nv))
     endif
     call cgyro_init_collision

  endif

  ! 5. Write initial data
  call cgyro_write_initdata

  if (test_flag == 1) return

  !---------------------------------------------------------------------------
  !
  ! Time-stepping
  n_time = nint(max_time/delta_t)

  ! Initialize h (via restart or analytic IC)
  call cgyro_init_h

  if (restart_flag == 0) then
     io_control = 1*(1-silent_flag)
  else
     io_control = 3*(1-silent_flag)
  endif
  call cgyro_write_timedata
  io_control = 2*(1-silent_flag)

  do i_time=1,n_time

     t_current = t_current+delta_t

     ! Collisionless step: returns new h_x, cap_h_x, fields 
     call cgyro_step_gk

     ! Spectral ExB shear
     if (abs(gamma_e) > 1e-10) call cgyro_shear

     ! Collisionless implicit streaming term step
     ! : returns new h_x, cap_h_x, fields 
     call cgyro_step_implicit_gk

     ! Collision step: returns new h_x, cap_h_x, fields
     call cgyro_step_collision

     ! Compute fluxes
     call cgyro_flux

     ! Error estimate
     call cgyro_error_estimate

     ! Print results
     call cgyro_write_timedata

     ! Print restart data
     call cgyro_write_restart

     if (abs(signal) == 1) exit

  enddo
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ! Print timers
  if (i_proc == 0) then
     print *
     print '(a)', 'Timing Summary'
     print '(a,1x,1pe11.4)',' init_arrays ',timer_lib_time('init_arrays')
     print '(a,1x,1pe11.4)',' field_v     ',timer_lib_time('field_v')
     print '(a,1x,1pe11.4)',' field_c     ',timer_lib_time('field_c')
     print '(a,1x,1pe11.4)',' rhs         ',timer_lib_time('rhs')
     print '(a,1x,1pe11.4)',' rhs_impgk   ',timer_lib_time('rhs_impgk')
     print '(a,1x,1pe11.4)',' rhs_impphi  ',timer_lib_time('rhs_impphi')
     print '(a,1x,1pe11.4)',' rhs_nl      ',timer_lib_time('rhs_nl')
     print '(a,1x,1pe11.4)',' coll_set1   ',timer_lib_time('coll_set1')
     print '(a,1x,1pe11.4)',' coll_set2   ',timer_lib_time('coll_set2')
     print '(a,1x,1pe11.4)',' collision   ',timer_lib_time('collision')
     print '(a,1x,1pe11.4)',' comm        ',timer_lib_time('comm')
     print '(a,1x,1pe11.4)',' comm_nl     ',timer_lib_time('comm_nl')
  endif
  !---------------------------------------------------------------------------

100 continue

  if(allocated(theta))          deallocate(theta)
  if(allocated(thetab))         deallocate(thetab)
  if(allocated(w_theta))        deallocate(w_theta)
  if(allocated(Bmag))           deallocate(Bmag)
  if(allocated(k_perp))         deallocate(k_perp)
  if(allocated(omega_stream))   deallocate(omega_stream)
  if(allocated(omega_trap))     deallocate(omega_trap)
  if(allocated(omega_rdrift))   deallocate(omega_rdrift)
  if(allocated(omega_adrift))   deallocate(omega_adrift)
  if(allocated(omega_aprdrift)) deallocate(omega_aprdrift)
  if(allocated(omega_cdrift))   deallocate(omega_cdrift)
  if(allocated(omega_gammap))   deallocate(omega_gammap)

  call GEO_alloc(0)

  if(allocated(indx_xi))       deallocate(indx_xi)
  if(allocated(px))            deallocate(px)
  if(allocated(energy))        deallocate(energy)
  if(allocated(w_e))           deallocate(w_e)
  if(allocated(e_deriv1_mat))  deallocate(e_deriv1_mat)
  if(allocated(e_deriv2_mat))  deallocate(e_deriv2_mat)
  if(allocated(xi))            deallocate(xi)
  if(allocated(w_xi))          deallocate(w_xi)
  if(allocated(xi_lor_mat))    deallocate(xi_lor_mat)
  if(allocated(xi_deriv_mat))  deallocate(xi_deriv_mat)
  if(allocated(xi_upderiv_mat)) deallocate(xi_upderiv_mat)
  if(allocated(h_x))           deallocate(h_x)
  if(allocated(cap_h_c))       deallocate(cap_h_c)
  if(allocated(cap_h_v))       deallocate(cap_h_v)
  if(allocated(field))         deallocate(field)
  if(allocated(field_loc))     deallocate(field_loc)
  if(allocated(field_old))     deallocate(field_old)
  if(allocated(f_balloon))     deallocate(f_balloon)
  if(allocated(hzf))           deallocate(hzf)
  if(allocated(xzf))           deallocate(xzf)
  if(allocated(pvec_in))       deallocate(pvec_in)
  if(allocated(pvec_outr))     deallocate(pvec_outr)
  if(allocated(pvec_outi))     deallocate(pvec_outi)

  if(collision_model /= 0) then
     if(allocated(cmat))       deallocate(cmat)
     if(allocated(cvec))       deallocate(cvec)
     if(allocated(bvec))       deallocate(bvec)
  endif

  call cgyro_clean_implicit_gk

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
     field_loc   = 3.0*field_old-3.0*field_old2+field_old3
     field_error = sum(abs(field-field_loc))/sum(abs(field))

  endif

  field_old3 = field_old2
  field_old2 = field_old
  field_old  = field

end subroutine cgyro_error_estimate
