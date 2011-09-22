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
  integer :: ix, ir
  integer :: itime, nt_step
  real    :: max_time=1000.0

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

  call EQUIL_alloc(1)
  call GYRO_alloc(1)
  call GK_alloc(1)
  call COLLISION_alloc(1)
  call FREQ_alloc(1)

  h_x     = 0.0
  cap_h_p = (0.0,0.0)
  phi = 1.0e-3
  nt_step = nint(max_time/delta_t)
  do itime =1, nt_step

     ! Collisionless gyrokinetic equation
     ! Returns new h_x, cap_h_x, cap_h_p, and phi 
     call GK_do

     ! Collision step
     ! Returns new cap_h_p and phi
     call COLLISION_do

     ! Compute frequency and check for convergence
     call FREQ_do

     if(abs(freq_err) < freq_tol) exit

     phi_old = phi

  enddo

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


end subroutine gkcoll_do
