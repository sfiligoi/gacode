subroutine prgen_read_er

  use prgen_read_globals 

  implicit none

  integer :: n_in
  integer :: i
  integer :: ierr
  integer :: count
  real :: xa,xb
  real :: fa,fb
  real, dimension(2) :: x
  real, dimension(:), allocatable :: rho_in
  real, dimension(:), allocatable :: f_in
  real, dimension(:), allocatable :: f_temp

  !----------------------------------------------------
  ! Determine vector length
  !
  count = 0 
  open(unit=1,file=er_file,status='old')
  do 
     read(1,*,iostat=ierr) x(:)
     if (ierr < 0) exit
     count = count+1
  enddo
  close(1)
  !
  n_in = count
  !----------------------------------------------------

  !----------------------------------------------------
  ! Read data
  !
  allocate(rho_in(n_in))
  allocate(f_in(n_in))
  !
  open(unit=1,file=er_file,status='old')
  do i=1,n_in 
     read(1,*) x(:)
     rho_in(i) = x(1)
     ! Scale Er according to input units
     ! er_conv = 3 : kV/m, no scaling  
     f_in(i)   = x(2)*10.0**(er_conv-3)
  enddo

  select case (er_conv)

  case (0)

     print '(t2,a,a)',trim(er_file),': Expecting [rho,Er(V/m),...]'

  case (3)

     print '(t2,a,a)',trim(er_file),': Expecting [rho,Er(kV/m),...]'

  case (6)

     print '(t2,a,a)',trim(er_file),': Expecting [rho,Er(MV/m),...]'

  case default

     print *,'Error'
     stop

  end select

  close(1)
  !----------------------------------------------------

  !----------------------------------------------------
  ! Correct for R_major grid
  !
  if (rho_in(n_in) > 1.01) then

     allocate(f_temp(n_in))
     f_temp(:) = f_in(:)

     print *, 'WARNING: rho was likely the major radius grid.q'
     print *, '         ... generating rho grid to correct.'
     do i=1,n_in
        rho_in(i) = (i-1.0)/(n_in-1)
        f_in(i) = f_temp(n_in-i+1)
     enddo
     print *, 'Er data has R(1) < R(2) for rho(1) < rho(2), ie reverse order'
     print *, 'Er data corrected from reverse order'

     deallocate(f_temp)

  endif
  !----------------------------------------------------

  !----------------------------------------------------
  ! Map last gridpoint so that rho_er(n_er) = 1.
  !
  xa = rho_in(n_in-1)
  xb = rho_in(n_in)
  fa = f_in(n_in-1)
  fb = f_in(n_in)
  !
  rho_in(n_in) = 1.0
  f_in(n_in) = fa*(xb-1.0)/(xb-xa)+fb*(1.0-xa)/(xb-xa)
  !----------------------------------------------------

  !----------------------------------------------------
  ! Spline interpolation from rho_er to rho grid.
  !
  call cub_spline(rho_in,f_in,n_in,rho,er_exp,nx)
  !----------------------------------------------------

  deallocate(rho_in) 
  deallocate(f_in)

end subroutine prgen_read_er
