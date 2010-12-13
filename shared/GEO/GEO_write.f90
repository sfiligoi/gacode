subroutine GEO_write(datafile,io)

  use GEO_interface

  implicit none

  !-------------------------------------------------------
  integer, intent(in) :: io
  character (len=*), intent(in) :: datafile
  !
  integer :: i
  integer :: n_theta
  !
  character (len=12), dimension(6) :: tag
  real, parameter :: pi=3.141592653589793
  !-------------------------------------------------------

  n_theta = size(GEOV_b)

  open(unit=io,file=datafile,status='replace')

  write(io,*) '----------------------------------------------------------------'
  write(io,*) ' Sample output from GEO'
  write(io,*) '----------------------------------------------------------------'
  write(io,*)
  write(io,'(a,i4)') 'n_theta = ',n_theta
  write(io,*)
  write(io,10) 'GEO_fluxsurfave_grad_r',GEO_fluxsurfave_grad_r
  write(io,10) 'GEO_fluxsurfave_grad_r2',GEO_fluxsurfave_grad_r2
  write(io,10) 'GEO_volume',GEO_volume
  write(io,10) 'GEO_volume_prime',GEO_volume_prime
  write(io,10) 'GEO_volume/circle_vol',GEO_volume/(2*pi**2*GEO_rmin_in**2*GEO_rmaj_in)
  write(io,10) 'GEO_volume/circle_volp',GEO_volume_prime/(4*pi**2*GEO_rmin_in*GEO_rmaj_in)
  write(io,10) 'ffprime',GEO_ffprime
  write(io,10) 'f',GEO_f
  write(io,10) 'f/R0',GEO_f/GEO_rmaj_in
  write(io,10) 'beta_star (beta_unit*dlnpdr)',GEO_beta_star
  write(io,*)

  tag(1) = 'theta'
  tag(2) = 'b'
  tag(3) = 'dbdt'
  tag(4) = 'dbdt2'
  tag(5) = 'g_theta'
  tag(6) = 'grad_r'

  write(1,*)
  write(io,30) tag(:)
  write(1,*)

  do i=1,n_theta
     write(1,20) GEOV_theta(i),GEOV_b(i),GEOV_dbdt(i),GEOV_dbdt2(i),&
          GEOV_g_theta(i),GEOV_grad_r(i)
  enddo

  tag(1) = 'theta'
  tag(2) = 'gsin'
  tag(3) = 'gcos1'
  tag(4) = 'gcos2'
  tag(5) = 'gq'
  tag(6) = 'captheta'

  write(1,*)
  write(1,*)  tag(:)
  write(1,*)
  do i=1,n_theta
     write(1,20) GEOV_theta(i),GEOV_gsin(i),GEOV_gcos1(i),GEOV_gcos2(i),&
          GEOV_gq(i),GEOV_captheta(i)
  enddo

  tag(1) = 'theta'
  tag(2) = 'theta_nc'
  tag(3) = '(null)'
  tag(4) = '(null)'
  tag(5) = '(null)'
  tag(6) = '(null)'

  write(1,*)
  write(1,*)  tag(:)
  write(1,*)
  do i=1,n_theta
     write(1,20) GEOV_theta(i),GEOV_theta_nc(i),0.0,0.0,0.0,0.0
  enddo

  close(io)

10 format(t2,a,': ',f10.6)
20 format(10(es11.4,1x))
30 format(t2,6(a))

end subroutine GEO_write
