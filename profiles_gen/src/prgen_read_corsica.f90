!--------------------------------------------------------------
! prgen_read_corsica.f90
!
! PURPOSE:
!  Extract corsica data from native TEXT file.
!--------------------------------------------------------------

subroutine prgen_read_corsica

  use prgen_globals

  implicit none

  integer, parameter :: corsica_nvals = 81

  integer :: i
  real :: index
  character(len=16) :: buf

  ! profile lengths are always the same (I think)
  ! so arrays can be allocated ahead of time
  nx = corsica_nvals

  call allocate_corsica_vars
  call prgen_allocate

  open(unit=1,file=file_state,status='old')

  ! skip header
  do 
     read(1,'(a)',err=10) buf
     if (buf == '#New time slice ') goto 99
  enddo
10 continue

  print '(a)','ERROR: (prgen_read_corsica) Could not find time slice data'
  close(1)
  return

99 continue

  ! read scalar data
  read(1,100) corsica_time
  read(1,100) corsica_current
  read(1,100) corsica_loop_voltage
  read(1,100) corsica_fusion_power
  read(1,100) corsica_bootstrap_fraction
  read(1,100) corsica_betat
  read(1,100) corsica_betan
  read(1,100) corsica_taue
  read(1,100) corsica_q95
  read(1,100) corsica_r
  read(1,100) corsica_a
  read(1,100) corsica_kappa
  read(1,100) corsica_delta
  read(1,100) corsica_fusion_gain
  read(1,100) corsica_betap
  read(1,100) corsica_li3

  ! skip comments
  do i=1,4
     read(1,*)
  enddo

  ! read profiles
  do i=1,nx
     read(1,101) index, corsica_rho(i), corsica_r_a(i), dpsi(i), &
          corsica_vl(i), te_kev(i), ti_kev(i), corsica_ne(i), &
          corsica_ndt(i), corsica_nz(i), corsica_nalpha(i), &
          corsica_zeff(i), q(i), corsica_j(i), corsica_jbs(i)
  enddo

  close(1)
  
  rmin(:)   = 0.0
  rmaj(:)   = 0.0
  rho(:)    = 0.0
  omega0(:) = 0.0

100 format (1e12.5)
101 format (15e13.5)

end subroutine prgen_read_corsica
