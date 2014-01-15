!-----------------------------------------------------------
! gyro_nonlinear_transfer.f90
!
! PURPOSE:
!  Compute nonlinear transfer (n,p) spectrum 
!  Compute turbulent entropy (n,p) spectrum S_k
!  for both periodic boundary conditions
!  
!  d(S_k+W_k)/dt = Tr_k + 2.*gammaS_k*S_k
!
!  k = (kx,ky)
!
!  Sum_k Tr_k = 0
!
!  Sum_k 2.*gamma_k*S_k = F + D's
!
!  d(S+W)/dt = F+D's   
! see Eq (41) amd Eqs(42)-(47) Candy & Waltz PoP 13 (2006) 032310
! 
!  on time average  d(S_k+W_k)/dt = 0  and d(S+W)/dt = 0
!
!  so -Tr_k/2./S_k = gammaS_k = (F_k+D's_k)/2./S_k can be regarded as 
!   1/2 the rate of entropy generagtion
!  at the (n,p) = (ky,kx) mode
!
!  "two * gammaS_k"  because S_k is a bilinear object and we want to compare
!  gammaS_k to gamma_k (amplitude growth rate)  
!
!  At steady state -Sum_k Tr_k =  F + D's = 0, i.e. F > 0 = -D's
!
!  tpyically S >> W the small polarization part and we need not worky about it
!  S is the "bare" entropy of the entropy in the gyrocenters
!  actually we will track only the ion entropy
!
!  Essentially Tr_k = (g_k_star*NL_k + g_k*NL_k_star)/4 = (h_k_star*NL_k + h_k*NL_k_star)/4
!  with a flux surface average included
!
!   S_k ----> Eng_k   is GYRO lable used.
!-----------------------------------------------------------

subroutine gyro_nonlinear_transfer

  use mpi
  use gyro_globals
  use gyro_pointers

  !---------------------------------------------------------------
  implicit none
  !
  complex, dimension(n_stack,n_x,n_nek_loc_1) :: rp
  complex, dimension(n_stack,n_x,n_nek_loc_1) :: hp

  integer :: p

  real, dimension(n_stack,n_x,n_nek_loc_1) :: nl1
  real, dimension(n_stack,n_x,n_nek_loc_1) :: nl2

  real, dimension(n_x,2) :: nl_transfer_loc
  real :: dw
  !---------------------------------------------------------------

  nl_transfer_loc(:,:) = 0.0

  do is=1,n_kinetic

     ! get p for each n

     rp(:,:,:) = (0.0,0.0)
     hp(:,:,:) = (0.0,0.0)

     do p=1,n_x
        do i=1,n_x
           rp(:,p,:) = rp(:,p,:) + rhs(:,i,:,is)*cri(p,i)
           hp(:,p,:) = hp(:,p,:) + h(:,i,:,is)*cri(p,i)
        enddo
     enddo

     nl1(:,:,:) = 0.5*real(conjg(hp(:,:,:))*rp(:,:,:))
     nl2(:,:,:) = 0.5*real(conjg(hp(:,:,:))*hp(:,:,:)) 

     ! The 1/2 factor corresponds to entropy (S) definition used in Candy & Waltz 
     ! PoP 13 (2006) 032310 Eq. (42).  Since Eq. (42) is a volume integration, a 
     ! flux surface average is implicit.  That is included in the w_p(ie,i,k,is) 
     ! factor embedded in the loop below. Since this routine is to be used in  
     ! periodic runs only, i = ir_norm

     p_nek_loc = 0

     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        ie = nek_e(p_nek)  
        k  = nek_k(p_nek)   

        dw = w_p(ie,ir_norm,k)

        do m=1,n_stack
           nl_transfer_loc(:,1) = nl_transfer_loc(:,1)+nl1(m,:,p_nek_loc)*dw
           nl_transfer_loc(:,2) = nl_transfer_loc(:,2)+nl2(m,:,p_nek_loc)*dw
        enddo ! m

     enddo ! p_nek

  enddo ! is

  nl_transfer_loc(:,:) = nl_transfer_loc(:,:)/rhos_norm**2 

  if (n_1(in_1) /= 0) then
     nl_transfer_loc(:,:) = 2.0*nl_transfer_loc(:,:)
  endif

  call MPI_ALLREDUCE(nl_transfer_loc, &
       nl_transfer, &
       size(nl_transfer), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'*[gyro_nonlinear_transfer done]'
  endif

end subroutine gyro_nonlinear_transfer
