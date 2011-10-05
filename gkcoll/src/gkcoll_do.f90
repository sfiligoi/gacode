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
  implicit none
  integer :: ix, ir, it, jr, p
  integer :: itime, nt_step
  character(len=80)  :: runfile_phi   = 'out.gkcoll.phi'
  character(len=80)  :: runfile_phiB  = 'out.gkcoll.phiB'
  character(len=80)  :: runfile_hx    = 'out.gkcoll.hx'
  character(len=80)  :: runfile_grids = 'out.gkcoll.grids'
  character(len=80)  :: runfile_time  = 'out.gkcoll.time'
  integer :: myio = 20
  integer :: print_step=10
  complex, dimension(:,:), allocatable :: phi_B

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
  allocate(phi_B(n_radial,n_theta))

  call EQUIL_alloc(1)
  call EQUIL_do
  call GYRO_alloc(1)
  call GK_alloc(1)
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
     call GK_do

     ! Collision step
     ! Returns new cap_h_p and phi
     call COLLISION_do

     if(mod(itime,print_step) == 0) then

        ! Compute frequency and print
        call FREQ_do

        ! Print phi
        if(silent_flag == 0 .and. i_proc == 0) then
           open(unit=myio,file=trim(path)//runfile_time,status='old',&
                position='append')
           write(myio,'(1pe12.5)') (itime * delta_t)
           close(myio)
           open(unit=myio,file=trim(path)//runfile_phi,status='old',&
                position='append')
           write(myio,'(1pe12.5)') transpose(phi(:,:))
           close(myio)
           do ir=1,n_radial
              do it=1,n_theta
                 phi_B(ir,it) = phi(ir,it) &
                      *exp(-2*pi*i_c*indx_r(ir)*k_theta*rmin)
              enddo
           enddo
           open(unit=myio,file=trim(path)//runfile_phiB,status='old',&
                position='append')
           write(myio,'(1pe12.5)') transpose(phi_B(:,:))
           close(myio)
        end if
        
        ! Check for convergence
        if(abs(freq_err) < freq_tol) then
           if(silent_flag == 0 .and. i_proc == 0) then
              print *, 'Converged'
           endif
           exit
        endif
     endif
     
     phi_old = phi
     
  enddo

  if(silent_flag == 0 .and. i_proc == 0) then
     open(unit=myio,file=trim(path)//runfile_hx,status='old',&
          position='append')
     write(myio,'(1pe12.5)') h_x
     close(myio)
  endif
     

100 continue
  call EQUIL_alloc(0)
  call GYRO_alloc(0)
  call GK_alloc(0)
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
  if(allocated(phi_B))         deallocate(phi_B)

end subroutine gkcoll_do
