!-------------------------------------------------------------------------
! Dep_to_nubeam.f90
!
! PURPOSE:
!  Provides interface NUBEAM to Dep
!
!-------------------------------------------------------------------------

module  Dep_to_nubeam

  implicit none

    ! integer,parameter :: ie_max=8
    ! integer,parameter :: k_max=8
    ! real :: D_EP_starOchi_i_TGLF(ie_max,k_max,2)
    real :: D_EP_starOchi_i_TGLF(8,8,2)
     ! ie_max = 8 & k_max = 8  are aready specified?
!rew added 8.20.12
    real :: D_EP_rEOchi_i_TGLF(8,8,2)
    real :: D_EP_ErOchi_i_TGLF(8,8,2)
    real :: D_EP_EEOchi_i_TGLF(8,8,2)

end       module Dep_to_nubeam
