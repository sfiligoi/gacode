!--------------------------------------------------------------
! EXPRO_read_driver.f90
!
! PURPOSE:
!  Low-level read control for EXPRO.
!--------------------------------------------------------------

subroutine EXPRO_read_driver

  use EXPRO_globals
  use EXPRO_interface

  implicit none

  integer, parameter :: io=1
  integer :: ierr

  real, dimension(:), allocatable :: dummy

  !--------------------------------------------------------------
  ! READ generated (stripped) version of input.profiles:
  !
  open(unit=io,&
       file=trim(path)//'input.profiles.gen',&
       status='old',&
       iostat=ierr)

  read(io,*) EXPRO_ncol
  read(io,*) EXPRO_nblock
  read(io,*) EXPRO_n_exp
  read(io,*) EXPRO_b_ref
  read(io,*) EXPRO_arho

  allocate(dummy(EXPRO_n_exp))

  ! 1-5
  read(io,*) EXPRO_rho(:)
  read(io,*) EXPRO_rmin(:)
  read(io,*) EXPRO_rmaj(:)
  read(io,*) EXPRO_q(:)     ! |q|
  read(io,*) EXPRO_kappa(:)

  ! 6-10
  read(io,*) EXPRO_delta(:)
  read(io,*) EXPRO_te(:)
  read(io,*) EXPRO_ne(:)
  read(io,*) EXPRO_z_eff(:)
  read(io,*) EXPRO_w0(:)     ! Note that EXPRO_w0 has GYRO/NEO sign

  ! 11-15
  read(io,*) EXPRO_flow_mom(:)
  read(io,*) EXPRO_pow_e(:)
  read(io,*) EXPRO_pow_i(:)
  read(io,*) EXPRO_pow_ei(:)
  read(io,*) EXPRO_zeta(:)

  ! 16-20
  read(io,*) EXPRO_flow_beam(:)
  read(io,*) EXPRO_flow_wall(:)
  read(io,*) EXPRO_zmag(:)
  read(io,*) EXPRO_ptot(:)
  read(io,*) EXPRO_poloidalfluxover2pi(:)

  ! 21-25
  read(io,*) EXPRO_ni(1,:)
  read(io,*) EXPRO_ni(2,:)
  read(io,*) EXPRO_ni(3,:)
  read(io,*) EXPRO_ni(4,:)
  read(io,*) EXPRO_ni(5,:)

  ! 26-30
  read(io,*) EXPRO_ti(1,:)
  read(io,*) EXPRO_ti(2,:)
  read(io,*) EXPRO_ti(3,:)
  read(io,*) EXPRO_ti(4,:)
  read(io,*) EXPRO_ti(5,:)

  ! 31-35
  read(io,*) EXPRO_vtor(1,:)
  read(io,*) EXPRO_vtor(2,:)
  read(io,*) EXPRO_vtor(3,:)
  read(io,*) EXPRO_vtor(4,:)
  read(io,*) EXPRO_vtor(5,:)

  ! 36-40
  read(io,*) EXPRO_vpol(1,:)
  read(io,*) EXPRO_vpol(2,:)
  read(io,*) EXPRO_vpol(3,:)
  read(io,*) EXPRO_vpol(4,:)
  read(io,*) EXPRO_vpol(5,:)

  ! 41-42 
  read(io,*,iostat=ierr) EXPRO_pow_e_fus(:)
  if (ierr == 0) then
     read(io,*) EXPRO_pow_i_fus(:)
  else
     close(io)
     open(unit=io,file=trim(path)//trim(runfile),status='replace')
     print '(a)', 'INFO: (EXPRO_read_driver) Old input.profiles detected.  Please regenerate with profiles_gen.'
     EXPRO_pow_e_fus(:) = 0.0
     EXPRO_pow_i_fus(:) = 0.0
  endif

  close(io)
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  ! READ general shape coefficients if they exist: 
  !
  if (EXPRO_nfourier > 0) then
     call EXPRO_read_geo
  endif
  !--------------------------------------------------------------

  deallocate(dummy)

end subroutine EXPRO_read_driver
 
 
