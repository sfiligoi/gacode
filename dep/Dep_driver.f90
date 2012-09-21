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
                           !A_EP(ie,k,isig,n,nb)

  use Dep_from_nubeam      !T_EP_hat, aoLf_EP

  use Dep_to_nubeam        !D_EP_starOchi_i_TGLF(ie,k,isig)

  use tglf_pkg             !tglf package interface



  !---------------------------------------------------------------
  implicit none
  !---------------------------------------------------------------

      integer :: n
      integer :: nb
      integer :: ie
      integer :: k
      integer :: isig
 
      integer :: i_plot_output
      integer :: i_GYRO_bypass

      real :: D_EP_starOchi_i_pass_plus(ie_max) 
      real :: D_EP_starOchi_i_pass_minus(ie_max) 
      real :: D_EP_starOchi_i_trap_plus_minus(ie_max) 

      real :: D_EP_pass_plus_sum
      real :: D_EP_pass_minus_sum
      real :: D_EP_trap_plus_minus_sum

      real :: D_EP_starOchi_i_pass_plus_n(ie_max,n_max)
      real :: D_EP_starOchi_i_pass_minus_n(ie_max,n_max)
      real :: D_EP_starOchi_i_trap_plus_minus_n(ie_max,n_max)


      real, external :: BESJ0
      real :: ps_sum

      real :: chi_i_tglf_wt_sum

   print *, 'start start j0=',BESJ0(0.1)

   i_plot_output = 1

  ! Allocate arrays.
  !!!call Dep_allocate
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
  i_GYRO_bypass =0 
!  i_GYRO_bypass =1
  if(i_GYRO_bypass .eq. 1) then
!GYRO GA-std case run with ony leading mode

    chi_i_tglf_wt(:,:)=0.
 
    chi_i_tglf_wt(1,1) = 20.5873
    chi_i_tglf_wt(2,1) = 8.39682
    chi_i_tglf_wt(3,1) = 3.19437
    chi_i_tglf_wt(4,1) = 1.17799
    chi_i_tglf_wt(5,1) = 0.594771
    chi_i_tglf_wt(6,1) = 0.342773
    chi_i_tglf_wt(7,1) = 0.224239
    chi_i_tglf_wt(8,1) = 0.171076
    chi_i_tglf_wt(9,1) = 0.138766
    chi_i_tglf_wt(10,1) = 0.122791 

    chi_i_tglf_wt_sum = 0.
    do n=1,n_max
     chi_i_tglf_wt_sum = chi_i_tglf_wt_sum + chi_i_tglf_wt(n,1)
    enddo
     chi_i_tglf_wt(:,1)=chi_i_tglf_wt(:,1)/chi_i_tglf_wt_sum

    omega_tglf(:,:)= 0.0

    omega_tglf(1,1) = (-0.069,0.1009)
    omega_tglf(2,1) = (-0.185,0.2177)
    omega_tglf(3,1) = (-0.317,0.2560)
    omega_tglf(4,1) = (-0.438,0.2398)
    omega_tglf(5,1) = (-0.229,0.2018)
    omega_tglf(6,1) = (0.4269,0.2847)
    omega_tglf(7,1) = (0.4844,0.3681)
    omega_tglf(8,1) = (0.5415,0.4469)
    omega_tglf(9,1) = (0.5973,0.5305)
    omega_tglf(10,1) = (0.6555,0.6275)

    omega_tglf(:,2) = omega_tglf(:,1)


   endif ! GYRO by pass

    print *, 'Dep_driver print'
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
    print *, '----------------------------------------'
    print *, 'Dep_driver calling Dep_mainsub'
    print *, '----------------------------------------'
   call Dep_mainsub  ! which will call Dep_kernel again

  if(i_plot_output .eq. 1) then

   open(unit=2,file='chi_i_wts.out',status='replace')
     write(2,*) 'krho(n)'
   do n=1,n_max
        write(2,*)  krho(n)
   enddo
     write(2,*)  'chi_i_tglf_wt(n,1)'
   do n=1,n_max
       write(2,*) chi_i_tglf_wt(n,1)
   enddo
     write(2,*)  'chi_i_tglf_wt(n,2)'
   do n=1,n_max
       write(2,*) chi_i_tglf_wt(n,2)
   enddo
   close(2)

   open(unit=2,file='D_EP_starOchi_i_TGLF.out',status='replace')
     write(2,*) 'e_hat(ie)'
   do ie=1,ie_max
        write(2,*)  e_hat(ie)
   enddo
     write(2,*)  'lambda=',lambda(3)
     write(2,*)  'D_EP_starOchi_i_TGLF sample pass sig plus'
   do ie=1,ie_max
       write(2,*) D_EP_starOchi_i_TGLF(ie,3,1)
   enddo
     write(2,*)  'D_EP_starOchi_i_TGLF sample pass sig minus'
   do ie=1,ie_max
       write(2,*) D_EP_starOchi_i_TGLF(ie,3,2)
   enddo
     write(2,*)  'lambda=',lambda(6)
     write(2,*)  'D_EP_starOchi_i_TGLF sample trap sig plus'
   do ie=1,ie_max
       write(2,*) D_EP_starOchi_i_TGLF(ie,6,1)
   enddo
     write(2,*)  'D_EP_starOchi_i_TGLF sample trap sig minus'
   do ie=1,ie_max
       write(2,*) D_EP_starOchi_i_TGLF(ie,6,2)
   enddo
   close(2)

   open(unit=2,file='D_EP_starOchi_i_TGLF_lambda.out',status='replace')
     write(2,*) 'e_hat(ie)'
   do ie=1,ie_max
        write(2,*)  e_hat(ie)
   enddo
    do k=1,k_max
     write(2,*)  'k=',k
     write(2,*)  'lambda=',lambda(k)
     write(2,*)  'D_EP_starOchi_i_TGLF  sig plus'
   do ie=1,ie_max
       write(2,*) D_EP_starOchi_i_TGLF(ie,k,1)
   enddo
     write(2,*)  'D_EP_starOchi_i_TGLF  sig minus'
   do ie=1,ie_max
       write(2,*) D_EP_starOchi_i_TGLF(ie,k,2)
   enddo
    enddo
   close(2)

   open(unit=2,file='D_EP_rEOchi_i_TGLF_lambda_p_m.out',status='replace')
     write(2,*) 'e_hat(ie)'
   do ie=1,ie_max
        write(2,*)  e_hat(ie)
   enddo
    do k=1,k_max
     write(2,*)  'k=',k
     write(2,*)  'lambda=',lambda(k)
     write(2,*)  'D_EP_rEOchi_i_TGLF  sig plus'
   do ie=1,ie_max
       write(2,*) D_EP_rEOchi_i_TGLF(ie,k,1)
   enddo
     write(2,*)  'D_EP_rEOchi_i_TGLF  sig minus'
   do ie=1,ie_max
       write(2,*) D_EP_rEOchi_i_TGLF(ie,k,2)
   enddo
    enddo
   close(2)


   open(unit=2,file='D_EP_starOchi_i_TGLF_lambda_ave_sig.out',status='replace')
     write(2,*) 'e_hat(ie)'
   do ie=1,ie_max
        write(2,*)  e_hat(ie)
   enddo
    do k=1,k_max
     write(2,*)  'k=',k
     write(2,*)  'lambda=',lambda(k)
     write(2,*)  'D_EP_starOchi_i_TGLF  sig ave'
   do ie=1,ie_max
       write(2,*) (D_EP_starOchi_i_TGLF(ie,k,1)+D_EP_starOchi_i_TGLF(ie,k,2))/2.
   enddo
    enddo
   close(2)


!!  open(unit=2,file='D_EP_starOchi_i_TGLF_lambda_ave_Gsig.out',status='replace')
!     write(2,*) 'e_hat(ie)'
!   do ie=1,ie_max
!        write(2,*)  e_hat(ie)
!   enddo
!    do k=1,k_max
!     write(2,*)  'k=',k
!     write(2,*)  'lambda=',lambda(k)
!     write(2,*)  'D_EP_starOchi_i_TGLF  sig ave'
!   do ie=1,ie_max
!     write(2,*) sqrt(abs(D_EP_starOchi_i_TGLF(ie,k,1)*D_EP_starOchi_i_TGLF(ie,k,2)))
!   enddo
!    enddo
!   close(2)

   open(unit=2,file='D_EP_rEOchi_i_TGLF_lambda_ave_sig.out',status='replace')
     write(2,*) 'e_hat(ie)'
   do ie=1,ie_max
        write(2,*)  e_hat(ie)
   enddo
    do k=1,k_max
     write(2,*)  'k=',k
     write(2,*)  'lambda=',lambda(k)
     write(2,*)  'D_EP_rEOchi_i_TGLF  sig ave'
   do ie=1,ie_max
       write(2,*) (D_EP_rEOchi_i_TGLF(ie,k,1)+D_EP_rEOchi_i_TGLF(ie,k,2))/2.
   enddo
    enddo
   close(2)

!!   open(unit=2,file='D_EP_rEOchi_i_TGLF_lambda_Gave_sig.out',status='replace')
!     write(2,*) 'e_hat(ie)'
!   do ie=1,ie_max
!        write(2,*)  e_hat(ie)
!   enddo
!    do k=1,k_max
!     write(2,*)  'k=',k
!     write(2,*)  'lambda=',lambda(k)
!     write(2,*)  'D_EP_rEOchi_i_TGLF  sig ave'
!   do ie=1,ie_max
!       write(2,*) -sqrt(abs(D_EP_rEOchi_i_TGLF(ie,k,1)*D_EP_rEOchi_i_TGLF(ie,k,2)))
!   enddo
!    enddo
!!   close(2)



   open(unit=2,file='D_EP_EEOchi_i_TGLF_lambda_ave_sig.out',status='replace')
     write(2,*) 'e_hat(ie)'
   do ie=1,ie_max
        write(2,*)  e_hat(ie)
   enddo
    do k=1,k_max
     write(2,*)  'k=',k
     write(2,*)  'lambda=',lambda(k)
     write(2,*)  'D_EP_EEOchi_i_TGLF  sig ave'
   do ie=1,ie_max
       write(2,*) (D_EP_EEOchi_i_TGLF(ie,k,1)+D_EP_EEOchi_i_TGLF(ie,k,2))/2.
   enddo
    enddo
   close(2)

   open(unit=2,file='D_EP_starOchi_i_TGLF_max.out',status='replace')
     write(2,*) 'e_hat(ie)'
   do ie=1,ie_max
        write(2,*)  e_hat(ie)
   enddo
    D_EP_starOchi_i_pass_plus(:) = 0.
    D_EP_starOchi_i_pass_minus(:) = 0.
    D_EP_starOchi_i_trap_plus_minus(:) =0.
    do ie=1,ie_max
      do k=1,4
       D_EP_starOchi_i_pass_plus(ie) = D_EP_starOchi_i_pass_plus(ie)+&
         D_EP_starOchi_i_TGLF(ie,k,1)*lambda_wts(k)*e_wts(ie)
      D_EP_starOchi_i_pass_minus(ie) = D_EP_starOchi_i_pass_minus(ie)+&
         D_EP_starOchi_i_TGLF(ie,k,2)*lambda_wts(k)*e_wts(ie)
      enddo
       do k=5,8
       D_EP_starOchi_i_trap_plus_minus(ie) = D_EP_starOchi_i_trap_plus_minus(ie)+&
         D_EP_starOchi_i_TGLF(ie,k,1)*lambda_wts(k)*e_wts(ie)+&
         D_EP_starOchi_i_TGLF(ie,k,2)*lambda_wts(k)*e_wts(ie)
      enddo
    enddo

    D_EP_pass_plus_sum=0.
    write(2,*) 'D_EP_starOchi_i_pass_plus(ie)'
    do ie=1,ie_max
     write(2,*) D_EP_starOchi_i_pass_plus(ie)
     D_EP_pass_plus_sum=D_EP_pass_plus_sum+D_EP_starOchi_i_pass_plus(ie)
    enddo

    D_EP_pass_minus_sum=0.
    write(2,*) 'D_EP_starOchi_i_pass_minus(ie)'
    do ie=1,ie_max
     write(2,*) D_EP_starOchi_i_pass_minus(ie)
     D_EP_pass_minus_sum=D_EP_pass_minus_sum+D_EP_starOchi_i_pass_minus(ie)
    enddo

    D_EP_trap_plus_minus_sum=0.
    write(2,*) 'D_EP_starOchi_i_trap_plus_minus(ie)'
    do ie=1,ie_max
     write(2,*) D_EP_starOchi_i_trap_plus_minus(ie)
    D_EP_trap_plus_minus_sum=D_EP_trap_plus_minus_sum+D_EP_starOchi_i_trap_plus_minus(ie)
    enddo

    write(2,*) 'D_EP_starOchi_i_max_tot' 
    do ie=1,ie_max
    write(2,*) D_EP_starOchi_i_pass_plus(ie)+ &
         D_EP_starOchi_i_pass_minus(ie)+D_EP_starOchi_i_trap_plus_minus(ie)
    enddo

    write(2,*) 'D_EP_pass_plus_sum=',D_EP_pass_plus_sum
    write(2,*) 'D_EP_pass_minus_sum=',D_EP_pass_minus_sum
    write(2,*) 'D_EP_trap_plus_minus_sum=',D_EP_trap_plus_minus_sum
     write(2,*) 'D_EP_starOchi_i_TGLF_sum=', &
            D_EP_pass_plus_sum+D_EP_pass_minus_sum+D_EP_trap_plus_minus_sum
   close(2)

   open(unit=2,file='D_EP_starOchi_i_TGLF_max_n.out',status='replace')
     write(2,*) 'e_hat(ie)'
   do ie=1,ie_max
        write(2,*)  e_hat(ie)
   enddo
   do n=1,n_max
    D_EP_starOchi_i_pass_plus_n(:,n) = 0.
    D_EP_starOchi_i_pass_minus_n(:,n) = 0.
    D_EP_starOchi_i_trap_plus_minus_n(:,n) =0.
    do ie=1,ie_max
     do nb=1,nb_max
      do k=1,4
       D_EP_starOchi_i_pass_plus_n(ie,n) = D_EP_starOchi_i_pass_plus_n(ie,n)+&
         D_EP_starOchi_i_kernal(ie,k,1,n,nb)*lambda_wts(k)*e_wts(ie)*chi_i_tglf_wt(n,nb)
       D_EP_starOchi_i_pass_minus_n(ie,n) = D_EP_starOchi_i_pass_minus_n(ie,n)+&
         D_EP_starOchi_i_kernal(ie,k,2,n,nb)*lambda_wts(k)*e_wts(ie)*chi_i_tglf_wt(n,nb)
      enddo
       do k=5,8
       D_EP_starOchi_i_trap_plus_minus_n(ie,n) = D_EP_starOchi_i_trap_plus_minus_n(ie,n)+&
         D_EP_starOchi_i_kernal(ie,k,1,n,nb)*lambda_wts(k)*e_wts(ie)*chi_i_tglf_wt(n,nb)+&
         D_EP_starOchi_i_kernal(ie,k,2,n,nb)*lambda_wts(k)*e_wts(ie)*chi_i_tglf_wt(n,nb)
       enddo
      enddo !nb
     enddo !ie
    enddo !n

    do n=1,n_max
      write(2,*) 'n=',n,' krho=',krho(n)
      write(2,*) '--------------------------------------------------------'
      write(2,*) 'D_EP_starOchi_i_pass_plus_n'
      do ie=1,ie_max
       write(2,*) D_EP_starOchi_i_pass_plus_n(ie,n)
      enddo
      write(2,*) 'D_EP_starOchi_i_pass_minus_n'
      do ie=1,ie_max
       write(2,*) D_EP_starOchi_i_pass_minus_n(ie,n)
      enddo
      write(2,*) 'D_EP_starOchi_i_trap_plus_minus_n'
      do ie=1,ie_max
       write(2,*) D_EP_starOchi_i_trap_plus_minus_n(ie,n)
      enddo
      write(2,*) '--------------------------------------------------------'
    enddo !n

   close(2)


  endif

   


end program Dep_driver
