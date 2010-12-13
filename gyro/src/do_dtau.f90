!------------------------------------------------
! do_dtau.f90 [caller gyro_rhs_total]
!
! PURPOSE:
!  Compute the tau-derivative terms in the RHS 
!  (i.e., derivatives in the theta directions)
!  using 5-point rules.
!
! NOTES:
!  Upwind dissipation is necessary for Landau
!  damping of zonal flows.  No dissipation will
!  result in "bouncing" instead of damping.
!------------------------------------------------

subroutine do_dtau

  use gyro_globals
  use gyro_pointers

  !---------------------------------------
  implicit none
  !
  integer :: p
  !
  integer :: mff,mf
  integer :: mc
  integer :: mb,mbb
  !
  real, dimension(n_x) :: v0
  real :: upwind
  !
  complex :: pff,pf
  complex :: pc
  complex :: pb,pbb
  !
  complex, dimension(n_x,n_stack) :: fh0
  complex, dimension(n_x,n_stack) :: fh
  !
  complex :: temp1
  complex :: temp2
  !---------------------------------------

  RHS_dt = 0.0

  ! We should be calling this for gyrokinetic species only

  do is=1,n_gk

     upwind = orbit_upwind_vec(is)*abs(q0)/q0

     p_nek_loc = 0
     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        ie = nek_e(p_nek)  
        k  = nek_k(p_nek)   

        ck = class(k)

        do i=1,n_x

           fh0(i,:) = h(:,i,p_nek_loc,is) 
           fh(i,:)  = fh0(i,:)+ &
                z(is)*alpha_s(is,i)*gyro_u(:,i,p_nek_loc,is)

           v0(i) = mu_root(i,is)*v_theta(i,ie,k,is)/d_tau(ck)

           do m=1,n_stack

              p = p_phys(ck,m)

              mff = m_cyc(ck,m+2,p)
              mf  = m_cyc(ck,m+1,p)
              mc  = m_cyc(ck,m,p)
              mb  = m_cyc(ck,m-1,p)
              mbb = m_cyc(ck,m-2,p)

              pff = p_cyc(ck,i,m+2,p)
              pf  = p_cyc(ck,i,m+1,p)
              pc  = p_cyc(ck,i,m,p)
              pb  = p_cyc(ck,i,m-1,p)
              pbb = p_cyc(ck,i,m-2,p)

              temp1 = (-1.0/12.0)*pff*fh(i,mff) &
                   +(8.0/12.0)*pf*fh(i,mf) &
                   +(-8.0/12.0)*pb*fh(i,mb) &
                   +(1.0/12.0)*pbb*fh(i,mbb) 

              temp2 = (1.0/12.0)*pff*fh0(i,mff) &
                   +(-4.0/12.0)*pf*fh0(i,mf) &
                   +(6.0/12.0)*pc*fh0(i,mc) &
                   +(-4.0/12.0)*pb*fh0(i,mb) &
                   +(1.0/12.0)*pbb*fh0(i,mbb) 


              !------------------------------------------------
              ! Entropy tracking:
              !
              RHS_dt(m,i,p_nek_loc,is) =  RHS_dt(m,i,p_nek_loc,is) &
                   +real(v0(i)*(-upwind*temp2)*conjg(fh(i,m)))
              !------------------------------------------------

              RHS(m,i,p_nek_loc,is) = RHS(m,i,p_nek_loc,is)+ &
                   v0(i)*(sigma_tau(ck,p)*temp1-upwind*temp2)

           enddo ! i

        enddo ! m

     enddo ! p_nek

  enddo ! is

end subroutine do_dtau
