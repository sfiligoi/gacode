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
  use cgyro_equilibrium
  use cgyro_gyro
  use cgyro_gk
  use cgyro_field
  use cgyro_collision
  use cgyro_freq

  implicit none

  integer :: ix, ir, it
  character(len=14), dimension(3)  :: runfile_field = &
       (/'out.cgyro.phi ','out.cgyro.apar','out.cgyro.bpar'/)
  character(len=15), dimension(3)  :: runfile_fieldb = &
       (/'out.cgyro.phiB ','out.cgyro.aparB','out.cgyro.bparB'/)
  character(len=80)  :: runfile_hx    = 'out.cgyro.hx'
  character(len=80)  :: runfile_grids = 'out.cgyro.grids'
  character(len=80)  :: runfile_time  = 'out.cgyro.time'
  integer :: myio = 20
  integer :: print_step=10
  complex, dimension(:,:), allocatable :: f_balloon
  complex, dimension(:,:), allocatable :: h_x_glob
  complex :: a_norm

  integer :: signal,i_field
  logical :: lfe
  integer, dimension(:), allocatable :: recv_status

  if (silent_flag == 0 .and. i_proc == 0) then
     open(unit=io_cgyroout,file=trim(path)//runfile,status='replace')
     close(io_cgyroout)
  endif

  ! MPI setup
  call cgyro_mpi_grid

  call cgyro_make_profiles
  if(error_status > 0) goto 100
  call cgyro_check
  if(error_status > 0) goto 100

  allocate(indx_xi(n_xi))
  do ix=1,n_xi
     indx_xi(ix) = ix-1
  enddo
  allocate(indx_r(n_radial))
  do ir=1,n_radial
     indx_r(ir) = -n_radial/2 + (ir-1)
  enddo
  if(toroidal_model == 2) then
     if(n_radial /= 1) then
        print *, 'Error: For zf test, n_radial must be 1'
        stop
     endif
     indx_r(1) = 1
  endif

  ! set-up energy grid and weights
  allocate(energy(n_energy))
  allocate(w_e(n_energy))
  allocate(e_deriv1_mat(n_energy,n_energy))
  allocate(e_deriv2_mat(n_energy,n_energy))
  call pseudo_maxwell(n_energy,e_max,energy,w_e,e_deriv1_mat,e_deriv2_mat)
  ! for old calls, e_max is real; now e_max is an int
  !call energy_integral(n_energy,e_max,energy,w_e)
  !call cgyro_energy_mesh(n_energy,e_max,energy,w_e)

  ! set-up xi-grid, weights
  allocate(xi(n_xi))
  allocate(w_xi(n_xi))
  allocate(xi_lor_mat(n_xi,n_xi))
  allocate(xi_deriv_mat(n_xi,n_xi))
  call pseudo_legendre(n_xi,xi,w_xi,xi_deriv_mat,xi_lor_mat)

  ! allocate distribution function and field arrays
  allocate(h_x(nc,nv_loc))
  allocate(cap_h_c(nc,nv_loc))
  allocate(cap_h_ct(nv_loc,nc))
  allocate(cap_h_v(nc_loc,nv))
  allocate(cap_h_v_prime(nc_loc,nv))
  allocate(field(n_radial,n_theta,n_field))
  allocate(field_loc(n_radial,n_theta,n_field))
  allocate(field_old(n_radial,n_theta,n_field))
  allocate(f_balloon(n_radial,n_theta))

  call EQUIL_alloc(1)
  call EQUIL_do
  call GYRO_alloc(1)
  call FIELD_alloc(1)
  call GK_alloc(1)
  call COLLISION_alloc(1)
  call FREQ_alloc(1)

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

     do i_field=1,n_field 

        open(unit=myio,file=trim(path)//runfile_field(i_field),status='replace')
        close(myio)

        open(unit=myio,file=trim(path)//runfile_fieldb(i_field),status='replace')
        close(myio)

     enddo

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

  do itime = 1, nt_step

     ! Collisionless gyrokinetic equation
     ! Returns new h_x, cap_h_x, and fields 
     call GK_do

     ! Collision step
     ! Returns new h_x, cap_h_x, and fields
     call COLLISION_do

     if (mod(itime,print_step) == 0) then

        ! Compute frequency and print
        call FREQ_do

        ! Print fields
        if (silent_flag == 0 .and. i_proc == 0) then

           ! time
           open(unit=myio,file=trim(path)//runfile_time,status='old',&
                position='append')
           write(myio,'(1pe13.5e3)') itime*delta_t
           close(myio)

           ! Fields: (phi,A_parallel)
           do i_field=1,n_field

              open(unit=myio,file=trim(path)//runfile_field(i_field),status='old',&
                   position='append')
              write(myio,'(1pe13.5e3)') transpose(field(:,:,i_field))
              close(myio)

              ! Construct ballooning-space form of field

              a_norm = field(n_radial/2+1,n_theta/2+1,1) 
              do ir=1,n_radial
                 do it=1,n_theta
                    f_balloon(ir,it) = field(ir,it,i_field) &
                         *exp(-2*pi*i_c*indx_r(ir)*k_theta*rmin)
                 enddo
              enddo

              open(unit=myio,file=trim(path)//runfile_fieldb(i_field),status='old',&
                   position='append')
              write(myio,'(1pe13.5e3)') transpose(f_balloon(:,:)/a_norm)
              close(myio)

           enddo

        endif

        ! Check for convergence
        if(abs(freq_err) < freq_tol) then
           if(silent_flag == 0 .and. i_proc == 0) then
              print *, 'Converged'
           endif
           exit
        endif

        ! Check for manual halt signal
        if(i_proc == 0) then
           inquire(file='halt',exist=lfe)
           if (lfe .eqv. .true.) then
              open(unit=1,file='halt',status='old')
              read(1,*) signal
              close(1)
           else
              signal = 0
           endif
        endif
        if (abs(signal) == 1) then
           exit
        endif

     endif

     field_old  = field

  enddo

  ! Print final distribution
  if (silent_flag == 0) then

     if (i_proc_1 == 0) then
        open(unit=myio,file=trim(path)//runfile_hx,status='old',&
             position='append')
     endif

     allocate(recv_status(MPI_STATUS_SIZE))
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

     if (i_proc_1 == 0) then
        do iv=1,nv
           do ic=1,nc
              f_balloon(ir_c(ic),it_c(ic)) = h_x_glob(ic,iv) &
                   *exp(-2*pi*i_c*indx_r(ir_c(ic))*k_theta*rmin)
           enddo
           write(myio,'(1pe13.5e3)') transpose(f_balloon(:,:))
        enddo
     endif
     deallocate(h_x_glob)

     if (i_proc_1 == 0) close(myio)

  endif

  if (restart_write == 1) then
     open(unit=io_cgyroout,file=trim(path)//runfile_restart,status='replace')
     write(io_cgyroout,*) h_x
     close(io_cgyroout)
  endif

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
  call GYRO_alloc(0)
  call GK_alloc(0)
  call FIELD_alloc(0)
  call COLLISION_alloc(0)
  call FREQ_alloc(0)

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
