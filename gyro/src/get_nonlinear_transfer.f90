!-----------------------------------------------------------
! get_nonlinear_transfer.f90 [caller gyro_rhs_total.f90]
!
! PURPOSE:
!  Compute  nonlinear transfer (n,p) spectrum 
!  Compute  turbulent entropy  (n,p) spectrum S_k
!  for both periodic boundary  conditions
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

subroutine get_nonlinear_transfer

  use mpi
  use gyro_globals
  use gyro_pointers

  !---------------------------------------------------------------
  implicit none
  !
  complex, dimension(n_stack,n_x,n_nek_loc_1) :: NL_p
  complex, dimension(n_stack,n_x,n_nek_loc_1) :: h_p

  integer :: p

  real, dimension(n_stack,n_x,n_nek_loc_1) :: Tr_m_nek
  real, dimension(n_stack,n_x,n_nek_loc_1) :: Eng_m_nek

  real, dimension(n_x) :: Tr_p_loc
  real, dimension(n_x) :: Eng_p_loc


  Tr_p_loc(:)  = 0.0
  Eng_p_loc(:) = 0.0


  do is = 1,n_kinetic

     ! get p for each n

     NL_p(:,:,:) = 0.0
     h_p(:,:,:) = 0.0

     do p=1,n_x
        do i=1,n_x
           NL_p(:,p,:) = NL_p(:,p,:) + RHS(:,i,:,is)*cRi(p,i)
           h_p(:,p,:) = h_p(:,p,:) + h(:,i,:,is)*cRi(p,i)
        enddo
     enddo



     do p=1,n_x
        Tr_m_nek(:,p,:) = (conjg(h_p(:,p,:))*NL_p(:,p,:) + h_p(:,p,:)*conjg(NL_p(:,p,:)))/2.
        Eng_m_nek(:,p,:) = (conjg(h_p(:,p,:))*h_p(:,p,:) + h_p(:,p,:)*conjg(h_p(:,p,:)))/2. 
     enddo
     Tr_m_nek(:,:,:)  = Tr_m_nek(:,:,:)/2.
     Eng_m_nek(:,:,:)  = Eng_m_nek(:,:,:)/2.

     ! the 1/2 factor corresponds to entropy (S) definition used in Candy & Waltz PoP 13 (2006) 032310
     ! Eq. (42). Since Eq. (42) is a volume integration, a  flux surface average is implicit.
     ! That is included in the w_p(ie,i,k,is) factor imbedded in the loop below. Since this routine 
     ! is to be used in a cyclic flux tube without profile variation i = n_x/2

     ! Essentially Eng (p,n) is the p & n (kx,ky) transform of S also called "entopy(2)" in 
     ! get_entropy.f90


     p_nek_loc = 0

     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        ie = nek_e(p_nek)  
        k  = nek_k(p_nek)   

        ck = class(k)

        do m=1,n_stack

           m0 = m_phys(ck,m)

           Tr_p_loc(:) = Tr_p_loc(:)+Tr_m_nek(m,:,p_nek_loc)*w_p(ie,n_x/2,k,is)
           Eng_p_loc(:) = Eng_p_loc(:)+Eng_m_nek(m,:,p_nek_loc)*w_p(ie,n_x/2,k,is)

        enddo ! m

     enddo ! p_nek

  enddo ! is

  Tr_p_loc(:)= Tr_p_loc(:)/rhos_norm**2 
  Eng_p_loc(:)= Eng_p_loc(:)/rhos_norm**2
  if (n_1(in_1) /= 0) then
     Tr_p_loc(:)= 2.0*Tr_p_loc(:)
     Eng_p_loc(:) = 2.0*Eng_p_loc(:)
  endif


  call MPI_ALLREDUCE(Tr_p_loc, &
       Tr_p, &
       n_x, &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)


  call MPI_ALLREDUCE(Eng_p_loc, &
       Eng_p, &
       n_x, &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)


  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'*[get_nonlinear_transfer done]'
  endif

end subroutine get_nonlinear_transfer
