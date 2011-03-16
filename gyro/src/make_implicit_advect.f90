!------------------------------------------------
! make_implicit_advect.f90
!
! PURPOSE:
!  Generate matrices used for implicit electron 
!  advection, together with some associated 
!  operators used for time-stepping.
!------------------------------------------------

subroutine make_implicit_advect(i_print)

  use mpi
  use gyro_globals
  use gyro_pointers
  use math_constants

  !--------------------------------------------------
  implicit none
  !
  integer, intent(in) :: i_print
  !
  integer :: m_diff
  integer :: mp
  integer :: mpc
  integer :: n_s2
  integer :: m0p
  !
  complex :: pff,pf
  complex :: pc
  complex :: pb,pbb
  !
  real :: a_dt_v
  real :: x_s
  !
  complex, dimension(n_x,n_blend,8) :: imp_loc
  complex, dimension(n_x,n_blend,8) :: imp_glob
  complex, dimension(n_stack,n_stack) :: o_temp
  complex, dimension(2,n_x,n_stack,-2:2) :: deriv_4
  complex, dimension(:,:), allocatable :: o_work
  complex, dimension(:), allocatable :: work
  !--------------------------------------------------


  deriv_4 = (0.0,0.0)

  n_s2 = n_stack/2

  !-----
  ck = 1
  !-----

  do m=1,n_s2

     do i=1,n_x

        pff = p_cyc(ck,i,m+2,1)
        pf  = p_cyc(ck,i,m+1,1)
        pc  = p_cyc(ck,i,m,1)
        pb  = p_cyc(ck,i,m-1,1)
        pbb = p_cyc(ck,i,m-2,1)

        deriv_4(ck,i,m,2) = (-1.0/12.0)*pff
        deriv_4(ck,i,m,1) = (8.0/12.0)*pf
        deriv_4(ck,i,m,-1)  = (-8.0/12.0)*pb
        deriv_4(ck,i,m,-2)  = (1.0/12.0)*pbb

     enddo ! i

  enddo ! m

  !-----
  ck = 2
  !-----

  do m=1,n_stack

     do i=1,n_x

        deriv_4(ck,i,m,2) = (-1.0/12.0)
        deriv_4(ck,i,m,1) = (8.0/12.0)
        deriv_4(ck,i,m,-1)  = (-8.0/12.0)
        deriv_4(ck,i,m,-2)  = (1.0/12.0)

     enddo ! i

  enddo ! m

  ! Renormalize smoother:

  !-----------------------------------
  ! Compute orbit-blending matrices 
  !
  !  [OF]  -> o_f
  !  [OFv] -> o_fv
  !
  !-----------------------------------

  o_f(:,:,:,:)  = (0.0,0.0)
  o_fv(:,:,:,:) = (0.0,0.0)

  do i=1,n_x

     p_nek_loc = 0

     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        ie = nek_e(p_nek)  
        k  = nek_k(p_nek)   

        ck = class(k)

        x_s = mu(n_spec)*sqrt(tem_s(n_spec,i))

        a_dt_v = a_SDIRK*dt*x_s*v_theta(i,ie,k,indx_e)/d_tau(ck)

        o_temp = (0.0,0.0)

        if (ck == 1) then

           ! PASSING

           do m=1,n_s2

              o_temp(m,m) = 1.0
              o_temp(m+n_s2,m+n_s2) = 1.0

              do m_diff=-2,2

                 mpc = m_cyc(ck,m+m_diff,1)

                 ! First bunch with v<0

                 o_temp(m,mpc) = o_temp(m,mpc) & 
                      -a_dt_v*deriv_4(ck,i,m,m_diff) 

                 ! Second bunch with v>0

                 o_temp(m+n_s2,mpc+n_s2) = o_temp(m+n_s2,mpc+n_s2) & 
                      +a_dt_v*deriv_4(ck,i,m,m_diff)

              enddo ! m_diff

           enddo ! m

           allocate(i_piv(n_s2))
           allocate(work(n_s2))
           allocate(o_work(n_s2,n_s2))

           o_work(:,:) = o_temp(1:n_s2,1:n_s2) 
           call ZGETRF(n_s2,n_s2,o_work,n_s2,i_piv,info)
           call ZGETRI(n_s2,o_work,n_s2,i_piv,work,n_s2,info)
           o_temp(1:n_s2,1:n_s2) = o_work(:,:)
           o_advect(1:n_s2,1:n_s2,i,p_nek_loc) = o_work(:,:)

           o_work(:,:) = o_temp(n_s2+1:n_stack,n_s2+1:n_stack) 
           call ZGETRF(n_s2,n_s2,o_work,n_s2,i_piv,info)
           call ZGETRI(n_s2,o_work,n_s2,i_piv,work,n_s2,info)
           o_temp(n_s2+1:n_stack,n_s2+1:n_stack) = o_work(:,:)
           o_advect(n_s2+1:n_stack,n_s2+1:n_stack,i,p_nek_loc) = o_work(:,:)

           do m=1,n_s2
              do mp=1,n_s2
                 do j=1,n_blend

                    o_f(j,m,i,p_nek_loc) = o_f(j,m,i,p_nek_loc) + &
                         o_temp(m,mp)*c_blend(j,mp,i,p_nek_loc)

                    o_fv(j,m,i,p_nek_loc) = o_fv(j,m,i,p_nek_loc) + &
                         o_temp(m,mp)*c_blend(j,mp,i,p_nek_loc)* &
                         v_para(mp,i,p_nek_loc,n_spec) 

                 enddo ! j
              enddo ! mp
           enddo ! m

           do m=n_s2+1,n_stack
              do mp=n_s2+1,n_stack
                 m0p = m_phys(1,mp)
                 do j=1,n_blend

                    o_f(j,m,i,p_nek_loc) = o_f(j,m,i,p_nek_loc) + &
                         o_temp(m,mp)*c_blend(j,m0p,i,p_nek_loc)

                    o_fv(j,m,i,p_nek_loc) = o_fv(j,m,i,p_nek_loc) + &
                         o_temp(m,mp)*c_blend(j,m0p,i,p_nek_loc)*&
                         v_para(mp,i,p_nek_loc,n_spec)

                 enddo ! j 
              enddo ! mp
           enddo ! m 

           deallocate(i_piv)
           deallocate(work)
           deallocate(o_work)

        else

           ! TRAPPED

           do m=1,n_stack

              o_temp(m,m) = 1.0

              do m_diff=-2,2

                 mpc = m_cyc(ck,m+m_diff,1)

                 ! Begin with v<0 and loop around

                 o_temp(m,mpc) = o_temp(m,mpc) & 
                      -a_dt_v*deriv_4(ck,i,m,m_diff)

              enddo ! m_diff

           enddo ! m

           allocate(i_piv(n_stack))
           allocate(work(n_stack))
           allocate(o_work(n_stack,n_stack))

           o_work(:,:) = o_temp(1:n_stack,1:n_stack) 
           call ZGETRF(n_stack,n_stack,o_work,n_stack,i_piv,info)
           call ZGETRI(n_stack,o_work,n_stack,i_piv,work,n_stack,info)
           o_temp(1:n_stack,1:n_stack) = o_work(:,:)
           o_advect(1:n_stack,1:n_stack,i,p_nek_loc) = o_work(:,:)

           do m=1,n_stack
              do mp=1,n_stack
                 m0p = m_phys(2,mp)
                 do j=1,n_blend

                    o_f(j,m,i,p_nek_loc) = o_f(j,m,i,p_nek_loc) + &
                         o_temp(m,mp)*c_blend(j,m0p,i,p_nek_loc)

                    o_fv(j,m,i,p_nek_loc) = o_fv(j,m,i,p_nek_loc) + &
                         o_temp(m,mp)*c_blend(j,m0p,i,p_nek_loc)*&
                         v_para(mp,i,p_nek_loc,n_spec) 

                 enddo ! j
              enddo ! mp
           enddo ! m 

           deallocate(i_piv)
           deallocate(work)
           deallocate(o_work)

        endif

     enddo ! p_nek

  enddo ! i

  if (i_print == 1) then
     do i=1,n_x
        print *,'i = ',i
        do m=1,n_stack
           print 10,real(o_f(:,m,i,1))
        enddo
        print *
        do m=1,n_stack
           print 10,aimag(o_f(:,m,i,1))
        enddo
     enddo
  endif

  !-------------------------------------------
  ! Using orbit-blending matrices o_f, o_fv, 
  ! compute the "electron implicit advection" 
  ! blending matrices imp(1:8)
  !-------------------------------------------

  do jp=1,n_blend

     imp_loc(:,:,:) = (0.0,0.0)

     do i=1,n_x

        p_nek_loc = 0

        do p_nek=1+i_proc_1,n_nek_1,n_proc_1

           p_nek_loc = p_nek_loc+1

           ie = nek_e(p_nek)
           k  = nek_k(p_nek)

           ck = class(k)

           do m=1,n_stack

              m0 = m_phys(ck,m)

              !--------------------------------------------------------
              ! L_P phi     = imp(1) phi + imp(2) A_par + (-2)*imp(6) B_par
              ! L_A A_par   = imp(3) phi + imp(4) A_par + imp(5) B_par
              ! L_B B_par   = imp(6) phi + imp(7) A_par + imp(8) B_par
              !--------------------------------------------------------

              ! Note that elements 2,3,5,6 carry a relative
              ! minus sign.  This follows from the minus sign 
              ! in G(phi-(vp/c)*Ap) + G_perp(m*v_perp^2)/(z*e*B)*B_par.

              imp_loc(i,:,1) = imp_loc(i,:,1)+&
                   alpha_s(n_spec,i)*cs_blend(:,m0,i,p_nek_loc)*&
                   o_f(jp,m,i,p_nek_loc)
              
              imp_loc(i,:,2) = imp_loc(i,:,2)-&
                   alpha_s(n_spec,i)*cs_blend(:,m0,i,p_nek_loc)*&
                   o_fv(jp,m,i,p_nek_loc)
              
              imp_loc(i,:,3) = imp_loc(i,:,3)+&
                   alpha_s(n_spec,i)*cs_blend(:,m0,i,p_nek_loc)*&
                   v_para(m,i,p_nek_loc,n_spec)*o_f(jp,m,i,p_nek_loc)
              
              imp_loc(i,:,4) = imp_loc(i,:,4)-&
                   alpha_s(n_spec,i)*cs_blend(:,m0,i,p_nek_loc)*&
                   v_para(m,i,p_nek_loc,n_spec)*o_fv(jp,m,i,p_nek_loc)
              
              imp_loc(i,:,5) = imp_loc(i,:,5)-&
                   alpha_s(n_spec,i)*cs_blend(:,m0,i,p_nek_loc)*&
                   tem_s(n_spec,i)*o_f(jp,m,i,p_nek_loc)*&
                   v_para(m,i,p_nek_loc,n_spec)*energy(ie,indx_e)*lambda(i,k)
              
              imp_loc(i,:,6) = imp_loc(i,:,6)+&
                   0.5*alpha_s(n_spec,i)*cs_blend(:,m0,i,p_nek_loc)*&
                   tem_s(n_spec,i)*o_f(jp,m,i,p_nek_loc)*&
                   energy(ie,indx_e)*lambda(i,k)
              
              imp_loc(i,:,7) = imp_loc(i,:,7)-&
                   0.5*alpha_s(n_spec,i)*cs_blend(:,m0,i,p_nek_loc)*&
                   tem_s(n_spec,i)*o_fv(jp,m,i,p_nek_loc)*&
                   energy(ie,indx_e)*lambda(i,k)
              
              imp_loc(i,:,8) = imp_loc(i,:,8)-&
                   0.5*alpha_s(n_spec,i)*cs_blend(:,m0,i,p_nek_loc)*&
                   tem_s(n_spec,i)**2 *o_f(jp,m,i,p_nek_loc)*&
                   energy(ie,indx_e)**2 * lambda(i,k)**2
              

           enddo ! m

        enddo ! p_nek_loc

     enddo ! i

     call MPI_ALLREDUCE(imp_loc,&
          imp_glob,&
          n_x*n_blend*8,&
          MPI_DOUBLE_COMPLEX,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

     do i=1,n_x
        do j=1,n_blend
           imp(i,j,jp,:) = imp_glob(i,j,:)
        enddo ! j
     enddo ! i

  enddo ! jp

  ! Print samples:

  if (i_print == 1) then
     do i=1,n_x
        print *,'i = ',i
        do j=1,n_blend
           print 10,real(imp(i,j,:,4))
        enddo
        print *
        do j=1,n_blend
           print 10,aimag(imp(i,j,:,4))
        enddo
     enddo
  endif

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[make_implicit_advect done]'
  endif

10 format(32(f11.5,1x))

end subroutine make_implicit_advect
