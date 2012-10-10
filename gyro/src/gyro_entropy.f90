!------------------------------------------------------------------
! gyro_entropy.f90
!
! PURPOSE:
!  Calculate individual terms from entropy balance equation for 
!  each species:
!
!   sigma - H* dH/dt + H* C(H) + H* Diss_tau + H* Diss_r = 0 
!
! NOTES:
!  entropy(1) : sigma
!  entropy(2) : - H* dH/dt
!  entropy(3) : H* C(H)
!  entropy(4) : H* Diss_tau
!  entropy(5) : H* Diss_r
!
!  See GYRO Technical manual for more details.
!------------------------------------------------------------------

subroutine gyro_entropy

  use mpi
  use gyro_globals
  use gyro_pointers
  use math_constants

  !-------------------------------------------------------------
  implicit none
  !
  integer :: ii
  !
  real, dimension(n_kinetic,n_entro) :: entropy_loc
  complex :: h_cap_dot_0
  !-------------------------------------------------------------


  entropy_loc(:,:) = 0.0

  !------------------------------------------------------------
  ! The first element is Sugama's sigma function (the surface-
  ! averaged anomalous entropy production rate):
  !
  do is=1,n_kinetic

     ! Sum entropy (sigma) over fields

     entropy(is,1) = 0.0
     do ix=1,n_field
        entropy(is,1) = entropy(is,1)+& 
             (dlnndr_s(is,ir_norm)-1.5*dlntdr_s(is,ir_norm))*gbflux(is,ix,1) &
             +dlntdr_s(is,ir_norm)/tem_s(is,ir_norm)*gbflux(is,ix,2) &
             +q_s(ir_norm)/r_s(ir_norm)*gamma_e_s(ir_norm)*gbflux(is,ix,3) &
             +1.0/tem_s(is,ir_norm)*gbflux(is,ix,4)
     enddo ! ix

  enddo ! is
  !------------------------------------------------------------

  !------------------------------------------------------------
  p_nek_loc = 0
  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek) 
     k  = nek_k(p_nek)

     ck = class(k)

     do is=1,n_kinetic

        do i=1,n_x
           do m=1,n_stack

              ! 2: - H* dH/dt

              ! Remove advective part from time derivative
              h_cap_dot_0 = h_cap_dot(m,i,p_nek_loc,is) &
                   -i_c*omega_eb_s(i)*h_cap(m,i,p_nek_loc,is)

              entropy_loc(is,2) = entropy_loc(is,2)-w_p(ie,i,k,is) &
                   *real(conjg(h_cap(m,i,p_nek_loc,is))*h_cap_dot_0)

              ! 3: collisional dissipation

              if (collision_flag == 1 .and. is == indx_e) then
                 entropy_loc(is,3) = entropy_loc(is,3)+w_p(ie,i,k,is) &
                      *real(conjg(h_cap(m,i,p_nek_loc,is))* &
                      (fb_coll(m,i,p_nek_loc)-f_coll(m,i,p_nek_loc))/dt)
              endif

              ! 4: orbit dissipation

              entropy_loc(is,4) = entropy_loc(is,4)+&
                   w_p(ie,i,k,is)*RHS_dt(m,i,p_nek_loc,is) 

              ! 5: radial dissipation

              entropy_loc(is,5) = entropy_loc(is,5)+&
                   w_p(ie,i,k,is)*RHS_dr(m,i,p_nek_loc,is) 

           enddo ! m
        enddo ! i
     enddo ! is

  enddo ! p_nek
  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Normalize local results: 
  !
  do ii=2,n_entro

     if (n_1(in_1) /= 0) then
        entropy_loc(:,ii) = 2.0*entropy_loc(:,ii)
     endif

     entropy_loc(:,ii) = entropy_loc(:,ii)/n_x
     entropy_loc(:,ii) = entropy_loc(:,ii)/rhos_norm**2

  enddo ! ii 
  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Sum over all processors: 
  !
  call MPI_REDUCE(entropy_loc(:,2:n_entro), &
       entropy(:,2:n_entro), &
       (n_entro-1)*n_kinetic, &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       0,&
       GYRO_COMM_WORLD, &
       i_err)
  !------------------------------------------------------------

  if (step == 0) entropy = 0.0

  if (i_proc == 0) then
     if (debug_flag == 1) print *,'[get_entropy called]'
  endif

end subroutine gyro_entropy
