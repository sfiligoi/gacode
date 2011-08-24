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
  use gkcoll_allocate_profile
  implicit none

  integer :: is, ie, ix, it, js, je, jx, jt
  integer :: itime, nt_step
  real    :: dtime=0.01
  real    :: max_time=1000.0


  ! theta derivative variables
  integer, dimension(-2:2) :: cderiv
  integer, dimension(:), allocatable :: thcyc

  ! kinetic equation terms 
  real :: stream, trap, driftx
  
  ! energy grid
  real, dimension(:), allocatable :: e_grid, w_e
  real, dimension(:), allocatable :: xi_grid, w_xi
  real, dimension(:,:), allocatable :: xi_mat, xi_mat_inv

  ! LAPACK
  integer :: info
  integer, dimension(:), allocatable :: i_piv
  real, dimension(:), allocatable :: work

  integer, parameter :: io_gkcoll=10, io_f=11
  character(len=80)  :: runfile_f = 'out.gkcoll.f'

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

  call EQUIL_alloc(1)

  ! cyclic index (for theta-periodicity)
  allocate(thcyc(1-n_theta:2*n_theta))
  do it=1,n_theta
     thcyc(it-n_theta) = it
     thcyc(it) = it
     thcyc(it+n_theta) = it
  enddo
  ! coefficients for 4th order centered derivative
  cderiv(-2) =  1
  cderiv(-1) = -8
  cderiv(0)  =  0
  cderiv(1)  =  8
  cderiv(2)  = -1

  ! set-up energy grid and weights
  allocate(e_grid(n_energy))
  allocate(w_e(n_energy))
  call energy_integral(n_energy,e_max,1,e_grid,w_e)

  ! set-up xi-grid, weights, and conversion matrix
  allocate(xi_grid(n_xi))
  allocate(w_xi(n_xi))
  call gauss_legendre(-1.0,1.0,xi_grid,w_xi,n_xi)
  allocate(xi_mat(n_xi,n_xi))
  allocate(xi_mat_inv(n_xi,n_xi))
  do ix=1,n_xi
     do jx=1,n_xi
        call compute_legendre(ix-1,xi_grid(jx),xi_mat(ix,jx))
     enddo
  enddo
  xi_mat_inv = xi_mat
  allocate(work(n_xi))
  allocate(i_piv(n_xi))
  call DGETRF(n_xi,n_xi,xi_mat_inv,n_xi,i_piv,info)
  call DGETRI(n_xi,xi_mat_inv,n_xi,i_piv,work,n_xi,info)
  deallocate(work)
  deallocate(i_piv)

  ! allocate distribution function and field arrays
  allocate(f(n_species,n_radial,n_theta,n_xi,n_energy))
  allocate(phi(n_radial,n_theta))

  nt_step = nint(max_time/dtime)
  do itime =1, nt_step
     ! Solve the gyrokinetic equation
  enddo

100 continue
  if(allocated(indx_xi))    deallocate(indx_xi)
  if(allocated(indx_r))     deallocate(indx_r)
  if(allocated(thcyc))      deallocate(thcyc)
  if(allocated(e_grid))     deallocate(e_grid)
  if(allocated(w_e))        deallocate(w_e)
  if(allocated(xi_grid))    deallocate(xi_grid)
  if(allocated(w_xi))       deallocate(w_xi)
  if(allocated(xi_mat))     deallocate(xi_mat)
  if(allocated(xi_mat_inv)) deallocate(xi_mat_inv)
  if(allocated(f))          deallocate(f)
  if(allocated(phi))        deallocate(phi)

  call EQUIL_alloc(0)

contains

  subroutine compute_legendre(n,arg,val)
    integer, intent (in) :: n
    real, intent (in) :: arg
    real, intent(out) :: val
    real :: pmm, pmmp1, pnn
    integer :: k
    pmm=1.0
    if(n==0) then
       val = pmm
    else
       pmmp1 = arg*pmm;
       if(n==1) then
          val = pmmp1
       else
          do k=2, n
             pnn = (arg*(2*k-1)*pmmp1 - (k-1)*pmm)/(1.0*k)
             pmm=pmmp1
             pmmp1=pnn
          enddo
          val = pnn
       end if
    end if
  end subroutine compute_legendre


end subroutine gkcoll_do
