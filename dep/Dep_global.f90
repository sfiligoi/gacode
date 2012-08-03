!-------------------------------------------------------------------------
! Dep_global.f90
!
! PURPOSE:
!  Provides interface within Dep
!-------------------------------------------------------------------------

module Dep_global 

!    use Dep_global          !e_hat(ie),lambda(k),sig(isig)
                             !e_wts(ie),lambda_wts(k),Fmax(ie)
                             !ie_max,ie,k_max,n_max,nb_max
                             !D_EP_starOchi_i_kernal(ie,k,isig,n,nb)
                             !A_EP(ie,k,isig,n)


  implicit none

  ! PRIVATE PARAMETERS
  ! Number  of n-modes:  n_max  = 10  !  tglf  low-k default
     !! integer, parameter :: n_max = 10
         integer, parameter :: n_max = 10
  ! Maximum number of unstable branches  nb_max = 2  !  tglf default
     !! integer, parameter :: nb_max = 2
         integer, parameter :: nb_max = 2

   integer, parameter :: ie_max=8
   integer, parameter :: k_max=8


  ! 


  real :: e_hat(ie_max)
  real :: lambda(k_max)
  real :: e_wts(ie_max)
  real :: lambda_wts(k_max)
  real :: Fmax(ie_max)

  real :: D_EP_starOchi_i_kernal(ie_max,k_max,2,n_max,nb_max)
  real :: D_EP_starOchi_i_kernal_den(n_max,nb_max)
  real :: A_EP(ie_max,k_max,2,n_max)

end module Dep_global
