!------------------------------------------------------------------
! Dep_driver.f90
!
! PURPOSE:
!  Given omega_tglf(n,nb), calls Dep_kernal.f90 to output
!  D_EP_starOchi_i_kernal(ie,k,isig,n,nb) for testing.
!------------------------------------------------------------------

program Dep_driver

  use Dep_from_tglf        !krho(n), omega_tglf(n,nb)
                           !chi_i_tglf_wt(n,nb), chi_i_tglf
                           !r_hat,rmaj_hat,q_saf, Ti_hat, aoLT_i, aoLn_i

  use Dep_global           !e_hat(ie),lambda(k),sig(isig)
                           !e_wts(ie),lambda_wts(k),Fmax(ie)
                           !ie_max,ie,k_max,n_max,nb_max
                           !D_EP_starOchi_i_kernal(ie,k,isig,n,nb)
                           !A_EP(ie,k,isig,n)

  use Dep_from_nubeam      !T_EP_hat, aoLf_EP

  use tglf_pkg             !tglf package interface



  !---------------------------------------------------------------
  implicit none
  !---------------------------------------------------------------

      integer :: n
      integer :: nb

      real, external :: BESJ0

   print *, 'start start j0=',BESJ0(0.1)

  ! Allocate arrays.
  call Dep_allocate
  !---------------------------------------------------------------
  ! Set all input parameters.
  call Dep_read_input
  !---------------------------------------------------------------


!!!    print *, 'before call to TGLF ---------------------------------------'
!!!    print *, 'chi_i_tglf=',chi_i_tglf
!!!    print *, 'r_hat=',r_hat
!!!    print *, 'rmaj_hat=',rmaj_hat
!!!    print *, 'q_saf=',q_saf
!!!    print *, 'Ti_hat=',Ti_hat
!!!    print *, 'aoLT_i=',aoLT_i
!!!    print *, 'aoLn_i=',aoLn_i
!!!    print *, 'ni_hat=',ni_hat
!!!   do n=1,n_max
!!!    print *,'krho=',krho(n)
!!!   enddo
!!!   do n=1,n_max
!!!    do nb=1,nb_max
!!!    print *, 'omega_tglf=',omega_tglf(n,nb)
!!!    print *, 'chi_i_tglf_wt=',chi_i_tglf_wt(n,nb)
!!!    enddo
!!!   enddo
!!!    print *,'T_EP_hat=',T_EP_hat
!!!    print *, '----------------------------------------'


  !---------------------------------------------------------------
  ! call tglf transport model and get Dep inputs 
  call tglf_TM
  call map_tglf_to_Dep
!
    print *, 'after call to TGLF ---------------------------------------'
   !print *, 'chi_i_tglf=',chi_i_tglf
    print *, 'chi_i_tglf=',chi_i_tglf/ni_hat/aoLT_i
       !may be a problem that ni_hat = ni/ne not what used in tglf call?
       !also....chi_i_tglf refers to only low-k total
    print *, 'r_hat=',r_hat
    print *, 'rmaj_hat=',rmaj_hat
    print *, 'q_saf=',q_saf
    print *, 'Ti_hat=',Ti_hat
    print *, 'aoLT_i=',aoLT_i
    print *, 'aoLn_i=',aoLn_i
    print *, 'ni_hat=',ni_hat
   do n=1,n_max
    print *,'krho=',krho(n)
   enddo
   do n=1,n_max
    do nb=1,nb_max
    print *, n,nb
    print *, 'omega_tglf=',omega_tglf(n,nb)
    print *, 'chi_i_tglf_wt=',chi_i_tglf_wt(n,nb)
    enddo
   enddo
    print *,'T_EP_hat=',T_EP_hat
    print *, '----------------------------------------'

  !---------------------------------------------------------------

  ! Set the passive diffusivity variable D_EP_starOchi_i_kernal(:,:,:,:,:)
!!!!!  call Dep_kernel
  !---------------------------------------------------------------

  ! test Dep_mainsub
   call Dep_mainsub  ! which will call Dep_kernel again

end program Dep_driver
