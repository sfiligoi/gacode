!------------------------------------------------------------
! prgen_map_ufile.f90
!
! PURPOSE:
!  Map UFILE data onto input.profiles standard.  
!  
! NOTES:
!  - See UFILE documentation at
!    http://tokamak-profiledb.ccfe.ac.uk/DOCS/PR08MAN/pdbman.html
!------------------------------------------------------------

subroutine prgen_map_ufile

  use prgen_globals

  implicit none
  real, dimension(:), allocatable :: powd_i
  real, dimension(:), allocatable :: powd_e
  real, dimension(:), allocatable :: powd_i_aux
  real, dimension(:), allocatable :: powd_e_aux
  real, dimension(:), allocatable :: powd_i_fus
  real, dimension(:), allocatable :: powd_e_fus
  real, dimension(:), allocatable :: powd_e_rad

  integer :: i,j

  allocate(powd_i(nx))
  allocate(powd_e(nx))
  allocate(powd_i_aux(nx))
  allocate(powd_e_aux(nx))
  allocate(powd_i_fus(nx))
  allocate(powd_e_fus(nx))
  allocate(powd_e_rad(nx))

  print '(a,i2)','INFO: (prgen) Number of ions: ',ufile_nion

  powd_i_aux(:) = ufile_qnbii(:) &
       +ufile_qicrhi(:) &
       +ufile_qechi(:) &
       +ufile_qlhi(:) 

  powd_e_aux(:) = ufile_qnbie(:) &
       +ufile_qicrhe(:)  &
       +ufile_qeche(:) &
       +ufile_qlhe(:) &
       +ufile_qohm(:) 

  powd_e_rad(:) = ufile_qrad(:)

  powd_e_fus(:) = ufile_qfuse(:)
  powd_i_fus(:) = ufile_qfusi(:)

  powd_i(:) = powd_i_aux(:) &
       +ufile_qei(:) &
       -ufile_qwalli(:)

  print '(a)','INFO: (prgen) i-power: (QNBII+QICHRI+QLHI+QECHI)+QEI-QWALLI'

  powd_e(:) = powd_e_aux(:) &
       -ufile_qei(:) &
       -ufile_qrad(:) &
       -ufile_qwalle(:) 

  print '(a)','INFO: (prgen) e-power: (QNBIE+QICHRE+QLHE+QECHE+QOHM)-QE-QRAD-QWALLE'

  ! Integrate these profiles (W/m^3) times 1e-6 to obtain MW
  call ufile_volint(1e-6*powd_i,pow_i,ufile_volume,nx)
  call ufile_volint(1e-6*powd_e,pow_e,ufile_volume,nx)
  call ufile_volint(1e-6*powd_i_aux,pow_i_aux,ufile_volume,nx)
  call ufile_volint(1e-6*powd_e_aux,pow_e_aux,ufile_volume,nx)
  call ufile_volint(1e-6*powd_i_fus,pow_i_fus,ufile_volume,nx)
  call ufile_volint(1e-6*powd_e_fus,pow_e_fus,ufile_volume,nx)
  call ufile_volint(1e-6*powd_e_rad,pow_e_rad,ufile_volume,nx)
  call ufile_volint(1e-6*ufile_qei,pow_ei,ufile_volume,nx)

  !---------------------------------------------------------
  ! Map profile data onto single array:
  !
  allocate(vec(n_indx,nx))
  vec(:,:) = 0.0
  !
  vec(1,:)  = rho(:)
  vec(2,:)  = rmin(:)
  vec(3,:)  = rmaj(:)
  vec(4,:)  = q(:)
  vec(5,:)  = kappa(:)
  vec(6,:)  = delta(:)
  vec(7,:)  = ufile_te(:)*1e-3
  vec(8,:)  = ufile_ne(:)*1e-19
  vec(9,:)  = ufile_zeff(:)
  vec(11,:) = 0.0
  vec(12,:) = pow_e(:)
  vec(13,:) = pow_i(:)
  vec(14,:) = pow_ei(:)
  vec(15,:) = 0.0
  vec(16,:) = 0.0
  vec(17,:) = 0.0
  vec(18,:) = 0.0
  vec(19,:) = ufile_pres(:)
  vec(20,:) = dpsi(:)

  !-----------------------------------------------------------------
  ! Construct ion densities and temperatures with reordering
  ! in general case.  Use vphi and vpol as temporary arrays.
  !
  do i=1,ufile_nion
     vec(30+i,:) = ufile_ni(:,i)*1e-19
     vec(35+i,:) = ufile_ti(:,i)*1e-3
  enddo

  ! reorder
  do i=1,5 
     vec(20+i,:) = vec(30+reorder_vec(i),:)
     vec(25+i,:) = vec(35+reorder_vec(i),:)
  enddo

  ! vphi
  vec(31:35,:) = 0.0

  ! Insert carbon toroidal velocity
  do i=1,5
     j = reorder_vec(i)
     if (nint(ufile_m(j)) == 12) then
        print '(a)', 'INFO: (prgen) Assuming VROT is the carbon toroidal rotation.'
        vec(30+i,:) = -ufile_vrot(:)*(rmaj(:)+rmin(:))
     endif
  enddo

  ! vpol
  vec(36:40,:) = 0.0

  ! Additional powers (fusion and radiation)
  ! * for iterdb, put all radiated power in pow_e_line
  vec(41,:) = pow_e_fus(:)
  vec(42,:) = pow_i_fus(:)
  vec(43,:) = 0.0
  vec(44,:) = 0.0
  vec(45,:) = pow_e_rad(:)

  ! Additional powers (external heating)
  vec(46,:) = pow_e_aux(:)
  vec(47,:) = pow_i_aux(:)

end subroutine prgen_map_ufile

!---------------------------------------------------------
! Simple routine to obtain volume integral of the UFILE 
! power densities by trapezoidal integration.
!---------------------------------------------------------
subroutine ufile_volint(f,fi,v,n)

  implicit none 

  integer, intent(in) :: n 
  real, intent(in) :: f(n),v(n)
  real, intent(inout) :: fi(n)
  integer :: i

  fi(1) = 0.0
  do i=2,n
     fi(i) = fi(i-1)+0.5*(f(i-1)+f(i))*(v(i)-v(i-1))
  enddo

end subroutine ufile_volint

