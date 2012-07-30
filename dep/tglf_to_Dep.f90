!-------------------------------------------------------------------------
! tglf_to_Dep.f90
!
! PURPOSE:
!  Provides interface TGLF to Dep
!
! CALLING SEQUENCE:
!
! NOTES:
!
!-------------------------------------------------------------------------

  module  tglf_to_Dep
!
  implicit none

  ! PRIVATE PARAMETERS
  ! Number  of n-modes:  n_max  = 10  ! same as low-k tglf default
    integer, parameter :: nk_max = 10
  ! Maximum number of unstable modes at each ky  nb_max = 2! tglf default 
    integer, parameter :: nmode_max =2

  ! CONTROL PARAMETERS
    integer              :: tglf_to_Dep_flag = 1

  ! INPUT PARAMETERS

   !local plasma parameters passed TGLF
   
    real :: r_hat
    real :: rmaj_hat
    real :: q_saf
    real :: Ti_hat  !  local Ti/Te
    real :: aoLn_i
    real :: aoLT_i
    real :: ni_hat  ! ni/ne


  !tglf mode values

    real ::  chi_i_tglf  !  = Q_i_hat/[ni_hat*Ti_hat*aoLT_i]
       ! refers to hydrogenic only ion energy transport  NOT total ions

    real ::  krho(nk_max)
       ! likely tglf default:  krho(1) = 0.1, krho(2) = 0.2, ...krho(10) = 1.0

    real ::  chi_i_tglf_wt(nk_max,nmode_max)
       ! chi_i_tglf_wt(n,nb) = chi_i_tglf(n,nb)/chi_i= Q_i_tglf(n,nb)/Q_i_tgld

       ! Sum_n Sum_nb chi_i_tglf_wt(n,nb)  should be 1.0
       !  for stable n,nb modes  chi_i_tglf_wt(n,nb) = 0.
     
       ! this should tbe fraction of the hydrogenic ion energy flux

    complex :: omega_tglf(nk_max,nmode_max)   !(omega_re,omega_i)_tglf

    complex :: Rem(nk_max,nmode_max) ! = delta_A_par_hat/delta_phi_hat

       !usual GYRO norms  delta_phi_hat = e*delta_phi/Te 
       !                   delta_A_par_hat = e*(c_s/c)*delta_A_par/Te

    end module tglf_to_Dep
!
    SUBROUTINE map_tglf_to_Dep
! PURPOSE: collects data from TGLF for use in Dep.
! It is assumed that TGLF was called and that the species index 2 
! in the put_species call is the main ions. 
! The integration assumes that the ky spectrum is equally spaced for the first n_max ky's. 
! This will be true if n_max <= 10 for kygrid_model_in = 1 in TGLF.
!
    USE tglf_pkg
    USE tglf_to_Dep
! 
    IMPLICIT NONE
!
    INTEGER :: i,j,k
    REAL :: dky
    COMPLEX :: xi=(0.0,1.0)
!
    CALL get_DEP_parameters(r_hat,rmaj_hat,q_saf,Ti_hat,aoLn_i,aoLT_i,ni_hat)
!    write(*,*)"DEP",r_hat,rmaj_hat,q_saf,Ti_hat,aoLn_i,aoLt_i,ni_hat
!
    chi_i_tglf = 0.0
    dky = get_ky_spectrum_out(1)
    do i=1,nk_max
      krho(i) = get_ky_spectrum_out(i)
!      write(*,*)i,"krho=",krho(i)
      do j=1,nmode_max
        omega_tglf(i,j) = get_eigenvalue_spectrum_out(2,i,j) &
        + xi*get_eigenvalue_spectrum_out(1,i,j)
        chi_i_tglf_wt(i,j) = dky*get_flux_spectrum_out(2,2,1,i,j) 
        chi_i_tglf = chi_i_tglf + chi_i_tglf_wt(i,j)
!        write(*,*)i,j,"chi_i_tglf_wt=",chi_i_tglf_wt(i,j)
!        write(*,*)"omega_tglf=",omega_tglf(i,j) 
      enddo
    enddo
!    write(*,*)"chi_i_tglf=",chi_i_tglf
! renormalize the weights 
    if(chi_i_tglf.ne.0.0)chi_i_tglf_wt(:,:) = chi_i_tglf_wt(:,:)/chi_i_tglf
!
    END SUBROUTINE map_tglf_to_Dep
