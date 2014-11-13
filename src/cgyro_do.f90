!-----------------------------------------------------------------
! cgyro_do.f90
!
! PURPOSE:
!  Subroutinized main cgyro program.  
!
! NOTES:
!  This can be called directly using the driver routine cgyro 
!  (in which case input data will read from input.dat) or called 
!  as a subroutine using cgyro_sub.
!-----------------------------------------------------------------

subroutine cgyro_do

  use timer_lib
  use mpi

  use cgyro_globals
  use cgyro_io
  use cgyro_equilibrium
  use cgyro_gk
  use cgyro_field
  use cgyro_collision

  implicit none

  complex, dimension(:,:), allocatable :: h_x_glob


  if (silent_flag == 0 .and. i_proc == 0) then
     open(unit=io_run,file=trim(path)//runfile,status='replace')
     close(io_run)
  endif

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
  allocate(gyrox_J0(n_species,n_radial,n_theta,n_energy,n_xi))
  allocate(j0_c(nc,nv_loc))
  allocate(j0_v(nc_loc,nv))
  allocate(h_x(nc,nv_loc))
  allocate(cap_h_c(nc,nv_loc))
  allocate(cap_h_ct(nv_loc,nc))
  allocate(cap_h_v(nc_loc,nv))
  allocate(cap_h_v_prime(nc_loc,nv))
  allocate(field(n_radial,n_theta,n_field))
  allocate(field_loc(n_radial,n_theta,n_field))
  allocate(field_old(n_radial,n_theta,n_field))
  allocate(f_balloon(n_radial,n_theta))
  allocate(recv_status(MPI_STATUS_SIZE))

  call EQUIL_alloc(1)
  call EQUIL_do

  !4. Array initialization
  call cgyro_init_arrays

  call FIELD_alloc(1)
  call GK_alloc(1)
  call COLLISION_alloc(1)

  ! Timer initialization
  call timer_lib_init('gk_init')
  call timer_lib_init('fieldx')
  call timer_lib_init('gkrhs')
  call timer_lib_init('collision')
  call timer_lib_init('comm')

  call timer_lib_in('gk_init')
  call GK_init
  call timer_lib_out('gk_init')

  if (silent_flag == 0 .and. i_proc == 0) then


     open(unit=myio,file=trim(path)//runfile_hx,status='replace')
     close(myio)

     open(unit=myio,file=trim(path)//runfile_time,status='replace')
     close(myio)

     open(unit=myio,file=trim(path)//runfile_grids,status='replace')
     write(myio,'(i4)') n_species
     write(myio,'(i4)') n_radial
     write(myio,'(i4)') n_theta
     write(myio,'(i4)') n_energy
     write(myio,'(i4)') n_xi
     write(myio,'(i4)') indx_r(:)
     write(myio,'(1pe12.5)') theta(:)
     write(myio,'(1pe12.5)') energy(:)
     write(myio,'(1pe12.5)') xi(:)
     write(myio,'(1pe12.5)') transpose(theta_B(:,:))
     close(myio)

  endif

  ! Time-stepping
  nt_step = nint(max_time/delta_t)

  io_control = 1*(1-silent_flag)
  call cgyro_write_timedata
  io_control = 2*(1-silent_flag)

  do itime=1,nt_step

     ! Collisionless gyrokinetic equation
     ! Returns new h_x, cap_h_x, and fields 
     call GK_do

     ! Collision step
     ! Returns new h_x, cap_h_x, and fields
     call COLLISION_do

     if (mod(itime,print_step) == 0) then
        call cgyro_write_timedata
     endif

     if (abs(signal) == 1) exit

     field_old  = field

  enddo

  ! Print final distribution
  if (silent_flag == 0 .and. n_toroidal == 1) then

     if (i_proc == 0) then
        open(unit=myio,file=trim(path)//runfile_hx,status='old',&
             position='append')
     endif

     allocate(h_x_glob(nc,nv))

     ! Collect distribution onto process 0
     call MPI_GATHER(h_x(:,:),&
          size(h_x),&
          MPI_DOUBLE_COMPLEX,&
          h_x_glob(:,:),&
          size(h_x),&
          MPI_DOUBLE_COMPLEX,&
          0,&
          NEW_COMM_1,&
          i_err)

     if (i_proc == 0) then
        do iv=1,nv
           do ic=1,nc
              f_balloon(ir_c(ic),it_c(ic)) = h_x_glob(ic,iv) &
                   *exp(-2*pi*i_c*indx_r(ir_c(ic))*k_theta*rmin)
           enddo
           write(myio,'(1pe13.5e3)') transpose(f_balloon(:,:))
        enddo
     endif
     deallocate(h_x_glob)

     if (i_proc == 0) close(myio)

  endif

  !if (restart_write == 1) then
  !   open(unit=io_run,file=trim(path)//runfile_restart,status='replace')
  !   write(io_run,*) h_x
  !   close(io_run)
  !endif

  ! Print timers
  if (i_proc == 0) then
     print *
     print '(a)', 'Timing Summary'
     print '(a,1x,1pe11.4)',' gk_init   ',timer_lib_time('gk_init')
     print '(a,1x,1pe11.4)',' fieldx    ',timer_lib_time('fieldx')
     print '(a,1x,1pe11.4)',' gkrhs     ',timer_lib_time('gkrhs')
     print '(a,1x,1pe11.4)',' collision ',timer_lib_time('collision')
     print '(a,1x,1pe11.4)',' comm      ',timer_lib_time('comm')
  endif


100 continue

  call EQUIL_alloc(0)
  call GK_alloc(0)
  call FIELD_alloc(0)
  call COLLISION_alloc(0)

  if(allocated(indx_xi))       deallocate(indx_xi)
  if(allocated(indx_r))        deallocate(indx_r)
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
  if(allocated(f_balloon))     deallocate(f_balloon)

end subroutine cgyro_do
