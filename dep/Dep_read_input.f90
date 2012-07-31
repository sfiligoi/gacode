!---------------------------------------------------------
! Dep_read_input.f90
!
! PURPOSE:
!  Reads all input parameters from input.Dep for testing
!  use by Dep_driver.f90
!  Inputs are read in order without pre-processing, so
!  the order of inputs must match the order here.
!
!---------------------------------------------------------

subroutine Dep_read_input

  use from_tglf_to_Dep     !krho(n), omega_tglf(n,nb)
                           !chi_i_tglf_wt(n,nb), chi_i_tglf
                           !r_hat,rmaj_hat,q_saf, Ti_hat, aoLT_i, aoLn_i

  use Dep_global           !e_hat(ie),lambda(k),sig(isig)
                           !e_wts(ie),lambda_wts(k),Fmax(ie)
                           !ie_max,ie,k_max,n_max,nb_max
                           !D_EP_starOchi_i_kernal(ie,k,isig,n,nb)
                           !A_EP(ie,k,isig,n)

  use from_nubeam_to_Dep   !T_EP_hat, aoLf_EP

  !--------------------------------------
  implicit none
  !
  real :: n_max_input
  real :: nb_max_input
  integer :: i_n
  integer :: i_nb
  character(len=80) :: comment
  !--------------------------------------

  open(unit=1,file='input.Dep',status='old')

  !----------------------------------------------------------
  ! Order and variable format in input.Dep must match here.
  !
  read(1,*) n_max_input
  read(1,*) nb_max_input
   print *, 'n_max_input=',n_max_input, '  nb_max_input=',nb_max_input

  do i_n = 1, n_max_input      ! Loop for  krho
   read(1,*) krho(i_n)
   print *, 'krho=',krho(i_n)
  enddo ! i_n

     print *, '----------------------------------------'
  do i_n = 1, n_max_input      ! Loop for complex omega_tglf
   do i_nb = 1, nb_max_input
     read(1,*) omega_tglf(i_n,i_nb)
     print *, 'omega_tglf=',omega_tglf(i_n,i_nb)
   enddo ! i_nb
  enddo ! i_n

     print *, '----------------------------------------'
  do i_n = 1, n_max_input      ! Loop for chi_i_tglf_wt QL weights
   do i_nb = 1, nb_max_input
     read(1,*) chi_i_tglf_wt(i_n,i_nb)
     print *, 'chi_i_tglf_wt=',chi_i_tglf_wt(i_n,i_nb)
   enddo ! i_nb
  enddo ! i_n
     print *, '----------------------------------------'

  read(1,*) chi_i_tglf
    print *, 'chi_i_tglf=',chi_i_tglf
  read(1,*) r_hat
    print *, 'r_hat=',r_hat
  read(1,*) rmaj_hat
    print *, 'rmaj_hat=',rmaj_hat
  read(1,*) q_saf
    print *, 'q_saf=',q_saf
  read(1,*) Ti_hat
    print *, 'Ti_hat=',Ti_hat
  read(1,*) aoLT_i
    print *, 'aoLT_i=',aoLT_i
  read(1,*) aoLn_i
    print *, 'aoLn_i=',aoLn_i

  read(1,*) T_EP_hat
    print *, 'T_EP_hat=',T_EP_hat

    print *, 'n_max=',n_max
    print *, 'nb_max',nb_max


end subroutine Dep_read_input
