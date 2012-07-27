!--------------------------------------------------------------
! Dep_allocate.f90
!
! PURPOSE:
!  Allocates all variable size arrays from Dep_globals and various
!  from_*_to_* modules.
!
!---------------------------------------------------------------

subroutine Dep_allocate

  use Dep_global
  use from_nubeam_to_Dep

  implicit none

  allocate(aoLf_EP(ie_max,k_max,2))

end subroutine Dep_allocate
