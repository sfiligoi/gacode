subroutine prgen_swap

  use vpro
  
  implicit none

  ! Perform reordering (may just be identity)

  call swapc(expro_type,expro_n_ion)
  call swapc(expro_name,expro_n_ion)
  call swap1(expro_z,expro_n_ion)
  call swap1(expro_mass,expro_n_ion)
  call swap2(expro_ni,expro_n_ion,expro_n_exp)
  call swap2(expro_ti,expro_n_ion,expro_n_exp)
  call swap2(expro_vpol,expro_n_ion,expro_n_exp)
  call swap2(expro_vtor,expro_n_ion,expro_n_exp)

end subroutine prgen_swap

subroutine swap1(f,nion)

  use prgen_globals

  implicit none

  integer, intent(in) :: nion
  real, intent(inout), dimension(nion) :: f

  integer :: i,ip
  real, dimension(nion) :: ftmp

  ftmp = f
  do i=1,nion
     ip   = reorder_vec(i)
     f(i) = ftmp(ip)
  enddo

end subroutine swap1

subroutine swap2(f,nion,nexp)

  use prgen_globals

  implicit none

  integer, intent(in) :: nion,nexp
  real, intent(inout), dimension(nion,nexp) :: f

  integer :: i,ip
  real, dimension(nion,nexp) :: ftmp

  ftmp = f
  do i=1,nion
     ip     = reorder_vec(i)
     f(i,:) = ftmp(ip,:)
  enddo

end subroutine swap2

subroutine swapc(f,nion)

  use prgen_globals

  implicit none

  integer, intent(in) :: nion
  character(len=7), intent(inout), dimension(nion) :: f

  integer :: i,ip
  character(len=7), dimension(nion) :: ftmp

  ftmp = f
  do i=1,nion
     ip   = reorder_vec(i)
     f(i) = ftmp(ip)
  enddo

end subroutine swapc
