!------------------------------------------------------------
! prgen_map_iterdb.f90
!
! PURPOSE:
!  Map native iterdb data onto input.profiles standard.  
!  
! NOTES:
!  - This routine is common to both text and NetCDF formats.
!  - See map_plasmastate.f90 for analogous routine for 
!    plasmastate data.
!------------------------------------------------------------

subroutine prgen_map_ufile

  use prgen_globals

  implicit none
  real, dimension(:), allocatable :: powd_i
  real, dimension(:), allocatable :: powd_e

  integer :: i


  allocate(powd_i(nx))
  allocate(powd_e(nx))

  powd_i(:) = ufile_qnbii(:) &
       +ufile_qicrhi(:) &
       +ufile_qei(:) &
       +ufile_qechi(:) &
       -ufile_qwalli(:)

  powd_e(:) = ufile_qnbie(:) &
       +ufile_qicrhe(:) &
       -ufile_qei(:) &
       +ufile_qrad(:) &
       +ufile_qeche(:) &
       +ufile_qohm(:) &
       -ufile_qwalle(:) 

  ! Convert W to MW (1e-6):
  call ufile_volint(rho,1e-6*powd_i,pow_i,ufile_volume,nx)
  call ufile_volint(rho,1e-6*powd_e,pow_e,ufile_volume,nx)
  call ufile_volint(rho,1e-6*ufile_qei,pow_ei_exp,ufile_volume,nx)

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
  vec(14,:) = pow_ei_exp(:)
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
  vec(31,:) = ufile_nm1(:)*1e-19
  if (ufile_nion > 1) vec(32,:) = ufile_nm2(:)*1e-19
  if (ufile_nion > 2) vec(33,:) = ufile_nm3(:)*1e-19

  vec(36,:) = ufile_ti(:)
  if (ufile_nion > 1) vec(37,:) = ufile_ti(:)
  if (ufile_nion > 2) vec(38,:) = ufile_ti(:)

  ! Beam ions: one only
  !do i=1,ufile_ibion
  ! ni
  !   vec(31+i+ufile_nion-1,:) = ufile_enbeam(:)*1e-19
  ! Ti: T[keV] = (p/n)[J]/1.6022e-16[J/eV]
  !   vec(36+i+ufile_nion-1,:) = ufile_pfast(:)/ufile_enbeam(:)/1.6022e-16
  !enddo

  ! reorder
  do i=1,5 
     vec(20+i,:) = vec(30+reorder_vec(i),:)
     vec(25+i,:) = vec(35+reorder_vec(i),:)
  enddo

  ! vphi
  vec(31:35,:) = 0.0

  ! Insert carbon toroidal velocity
  do i=1,5
     if (reorder_vec(i) == 2) then
        vec(30+i,:) = -ufile_vrot(:)*(rmaj(:)+rmin(:))
     endif
  enddo

  ! vpol
  vec(36:40,:) = 0.0

end subroutine prgen_map_ufile

subroutine ufile_volint(x,f,fi,v,n)

  implicit none 

  integer, intent(in) :: n 
  real, intent(in) :: x(n),f(n),v(n)
  real, intent(inout) :: fi(n)
  real, dimension(n) :: vp
  integer :: i

  ! dv/dx
  call bound_deriv(vp,v,x,n)

  fi(1) = 0.0
  do i=2,n
     fi(i) = fi(i-1)+0.5*(vp(i-1)*f(i-1)+vp(i)*f(i))*(x(i)-x(i-1))
  enddo

end subroutine ufile_volint

