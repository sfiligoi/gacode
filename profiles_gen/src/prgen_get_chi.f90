subroutine prgen_get_chi(n,q,psi,rho,torfluxa)

  implicit none

  integer, intent(in) :: n
  real, intent(in), dimension(n) :: q
  real, intent(in), dimension(n) :: psi
  real, intent(inout), dimension(n) :: rho
  real, intent(inout) :: torfluxa

  real, dimension(:), allocatable :: chi_t
  integer :: i 

  !---------------------------------------------------------
  ! Compute rho by integrating d(chi_t) = q d(psi_p)
  ! using the trapezoidal rule
  !
  allocate(chi_t(n))  
  !
  chi_t(1) = 0.0
  do i=2,n
     chi_t(i) = chi_t(i-1) + 0.5*(q(i)+q(i-1))*(psi(i)-psi(i-1))
  enddo
  ! 
  torfluxa = chi_t(n)
  !
  ! Normalized root of chi_t:
  !
  rho(:) = sqrt(chi_t(:)/chi_t(n))
  ! Prevent any roundoff
  rho(n) = 1.0
  deallocate(chi_t)
  !---------------------------------------------------------

end subroutine prgen_get_chi
