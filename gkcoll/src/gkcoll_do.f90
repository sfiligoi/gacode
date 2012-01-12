!-----------------------------------------------------------------
! gkcoll_do.f90
!
! PURPOSE:
!  Subroutinized main gkcoll program.  
!
! NOTES:
!  This can be called directly using the driver routine gkcoll 
!  (in which case input data will read from input.dat) or called 
!  as a subroutine using gkcoll_sub.
!-----------------------------------------------------------------

subroutine gkcoll_do

  use gkcoll_globals
  use gkcoll_equilibrium
  use gkcoll_gyro
  use gkcoll_gk
  use gkcoll_collision
  use gkcoll_freq
  use gkcoll_allocate_profile
  use gkcoll_implicit
  use gkcoll_allimplicit
  implicit none

  integer :: ix, ir, it, jr, p, is, ie
  integer :: itime, nt_step
  character(len=80)  :: runfile_phi   = 'out.gkcoll.phi'
  character(len=80)  :: runfile_phiB  = 'out.gkcoll.phiB'
  character(len=80)  :: runfile_hx    = 'out.gkcoll.hx'
  character(len=80)  :: runfile_grids = 'out.gkcoll.grids'
  character(len=80)  :: runfile_time  = 'out.gkcoll.time'
  integer :: myio = 20
  integer :: print_step=10
  complex, dimension(:,:), allocatable :: f_balloon

  integer :: signal
  logical :: lfe

  if (silent_flag == 0 .and. i_proc == 0) then
     open(unit=io_gkcollout,file=trim(path)//runfile_gkcollout,status='replace')
     close(io_gkcollout)
  endif

  call gkcoll_make_profiles
  if(error_status > 0) goto 100
  call gkcoll_check
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
  call energy_integral(n_energy,e_max,1,energy,w_e)

  ! set-up xi-grid, weights
  allocate(xi(n_xi))
  allocate(w_xi(n_xi))
  call gauss_legendre(-1.0,1.0,xi,w_xi,n_xi)

  ! allocate distribution function and field arrays
  allocate(h_x(n_species,n_radial,n_theta,n_energy,n_xi))
  allocate(cap_h_x(n_species,n_radial,n_theta,n_energy,n_xi))
  allocate(cap_h_p(n_species,n_radial,n_theta,n_energy,n_xi))
  allocate(phi(n_radial,n_theta))
  allocate(phi_old(n_radial,n_theta))
  allocate(f_balloon(n_radial,n_theta))

  if(imp_flag == 1) then
     trap_method = 1
  else if(imp_flag == 2) then
     collision_model = -1
  endif

  call EQUIL_alloc(1)
  call EQUIL_do
  call GYRO_alloc(1)
  call GK_alloc(1)
  if(imp_flag == 1) then
     call GKimp_alloc(1)
  else if(imp_flag == 2) then
     call GKallimp_alloc(1)
  endif
  call COLLISION_alloc(1)
  call FREQ_alloc(1)

  ! Initialization
  call GK_init

  if(silent_flag == 0 .and. i_proc == 0) then
     open(unit=myio,file=trim(path)//runfile_phi,status='replace')
     close(myio)
     open(unit=myio,file=trim(path)//runfile_phiB,status='replace')
     close(myio)
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
     ! Returns new h_x, cap_h_x, cap_h_p, and phi 
     if(imp_flag == 1) then
        call GKimp_do
     else if(imp_flag == 2) then
        call GKallimp_do
     else
        call GK_do
     end if

     ! Collision step
     ! Returns new cap_h_p, h_x, and phi
     call COLLISION_do

     if(mod(itime,print_step) == 0) then

        ! Compute frequency and print
        call FREQ_do

        ! Print phi
        if(silent_flag == 0 .and. i_proc == 0) then
           open(unit=myio,file=trim(path)//runfile_time,status='old',&
                position='append')
           write(myio,'(1pe13.5e3)') (itime * delta_t)
           close(myio)
           open(unit=myio,file=trim(path)//runfile_phi,status='old',&
                position='append')
           write(myio,'(1pe13.5e3)') transpose(phi(:,:))
           close(myio)

           ! Construct ballooning-space form of phi
           do ir=1,n_radial
              do it=1,n_theta
                 f_balloon(ir,it) = phi(ir,it) &
                      *exp(-2*pi*i_c*indx_r(ir)*k_theta*rmin)
              enddo
           enddo
           open(unit=myio,file=trim(path)//runfile_phiB,status='old',&
                position='append')
           write(myio,'(1pe13.5e3)') transpose(f_balloon(:,:))
           close(myio)
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

     phi_old = phi

  enddo

  if(silent_flag == 0 .and. i_proc == 0) then
     open(unit=myio,file=trim(path)//runfile_hx,status='old',&
          position='append')

     do is=1,n_species
        do ie=1,n_energy
           do ix=1,n_xi

              ! Construct ballooning-space form of h_x
              do ir=1,n_radial
                 do it=1,n_theta
                    f_balloon(ir,it) = h_x(is,ir,it,ie,ix) &
                         *exp(-2*pi*i_c*indx_r(ir)*k_theta*rmin)
                 enddo
              enddo

              write(myio,'(1pe13.5e3)') transpose(f_balloon(:,:))

           enddo
        enddo
     enddo

     close(myio)
  endif

  if(restart_write == 1) then
     open(unit=io_gkcollout,file=trim(path)//runfile_restart,status='replace')
     write(io_gkcollout,*) h_x
     close(io_gkcollout)
  endif

100 continue
  call EQUIL_alloc(0)
  call GYRO_alloc(0)
  call GK_alloc(0)
  if(imp_flag == 1) then
     call GKimp_alloc(0)
  else if(imp_flag == 2) then
     call GKallimp_alloc(0)
  endif
  call COLLISION_alloc(0)
  call FREQ_alloc(0)

  if(allocated(indx_xi))       deallocate(indx_xi)
  if(allocated(indx_r))        deallocate(indx_r)
  if(allocated(energy))        deallocate(energy)
  if(allocated(w_e))           deallocate(w_e)
  if(allocated(xi))            deallocate(xi)
  if(allocated(w_xi))          deallocate(w_xi)
  if(allocated(h_x))           deallocate(h_x)
  if(allocated(cap_h_x))       deallocate(cap_h_x)
  if(allocated(cap_h_p))       deallocate(cap_h_p)
  if(allocated(phi))           deallocate(phi)
  if(allocated(phi_old))       deallocate(phi_old)
  if(allocated(f_balloon))     deallocate(f_balloon)

end subroutine gkcoll_do
