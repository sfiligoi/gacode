!------------------------------------------------------------
! make_ampere_blend.f90
!
! PURPOSE:
!  Generate matrix L_A of paper by pushing the factor 
!  [-nabla_perp A] into a velocity integral. 
!------------------------------------------------------------

subroutine make_ampere_blend

  use mpi
  use gyro_globals
  use gyro_pointers
  use math_constants

  !---------------------------------------------------
  implicit none
  !
  real :: betae_eff
  !
  real :: v_perp0
  real :: x_s
  real :: rho_gyro
  real :: a_gyro
  real :: u_gyro
  real :: v_gyro
  complex :: inqr
  complex :: temp
  !
  complex, dimension(-mg_dx:mg_dx-ig_dx,n_blend,n_blend) :: vel_sum_loc
  complex, dimension(-mg_dx:mg_dx-ig_dx,n_blend,n_blend) :: vel_sum_glob
  complex, dimension(n_gk,-mg_dx:mg_dx-ig_dx) :: f_x
  complex, dimension(-m_dx:m_dx-i_dx) :: grad_perp_ap
  complex, dimension(n_gk,-mg_dx:mg_dx-ig_dx) :: ion_current
  !---------------------------------------------------


  betae_eff = betae_unit_norm*ampere_scale

  do i=1,n_x

     vel_sum_loc = (0.0,0.0)

     p_nek_loc = 0

     inqr = i_c*n_1(in_1)*q_s(i)/r_s(i) 

     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        ie = nek_e(p_nek)
        k  = nek_k(p_nek)

        ck = class(k)

        ! Now, compute m-projections of the RHS:

        do m=1,n_stack

           m0 = m_phys(ck,m)

           grad_perp_ap(:) =-2.0*rhos_norm**2/betae_eff*( &
                inqr**2*qrat_t(i,k,m)**2* &
                (1.0+captheta_t(i,k,m)**2)*w_d0(:)+ &
                2.0*inqr*qrat_t(i,k,m)*captheta_t(i,k,m)* &
                grad_r_t(i,k,m)*dr_eodr(i)*w_d1(:)+ &
                (grad_r_t(i,k,m)*dr_eodr(i))**2*w_d2(:))

           do is=1,n_gk       

              if (gyro_method == 1) then 

                 ion_current(is,:) = &
                      v_para(m,i,p_nek_loc,is)**2* &
                      alpha_s(is,i)*z(is)**2*w_gd0(:)

              else

                 !--------------------------------------------
                 ! Parameters which depend on equilibrium:
                 !
                 v_perp0 = sqrt(2.0*energy(ie,is)*lambda(i,k)*b0_t(i,k,m0))
                 x_s = sqrt(tem_s(is,i))/mu(is)/abs(z(is))
                 !
                 if (kill_gyro_b_flag == 0) then
                    rho_gyro = rhos_norm/b_unit_s(i)*v_perp0*x_s/b0_t(i,k,m0)
                 else
                    rho_gyro = rhos_norm/b_unit_s(i)*v_perp0*x_s
                 endif
                 !
                 a_gyro = grad_r_t(i,k,m0)/x_length*dr_eodr(i)
                 u_gyro = qrat_t(i,k,m0)*n_1(in_1)*q_s(i)/r_s(i)*captheta_t(i,k,m0)
                 v_gyro = qrat_t(i,k,m0)*n_1(in_1)*q_s(i)/r_s(i)
                 !--------------------------------------------

                 f_x(is,:) = (0.0,0.0) 

                 ! J_0^2
                 call gyro_bessel_operator(rho_gyro,&
                      a_gyro,&
                      u_gyro,&
                      v_gyro,&
                      f_x(is,:),&
                      2)

                 if (n_1(in_1) == 0) then   

                    if (i_gyro /= 1) then 

                       !! JC
                       ! Correct truncated gyroaverage                 

                       temp = sum(f_x(is,:))-1.0
                       f_x(is,0) = f_x(is,0)-temp

                    endif

                    ! Enforce EXACT reality 

                    f_x(is,:) = real(f_x(is,:))

                 endif

                 ion_current(is,:) = &
                      v_para(m,i,p_nek_loc,is)**2* &
                      alpha_s(is,i)*z(is)**2*f_x(is,:)

              endif

           enddo ! is

           do j=1,n_blend 
              do jp=1,n_blend

                 do i_diff=-m_dx,m_dx-i_dx

                    !----------------------------
                    ! This will stagnate at i=n_x  
                    ! for i+i_diff > n_x
                    !
                    ip = i_cyc(i+i_diff)
                    !----------------------------

                    vel_sum_loc(i_diff,j,jp) = vel_sum_loc(i_diff,j,jp)+&
                         grad_perp_ap(i_diff)* &
                         cs_blend(j,m0,i,p_nek_loc)* &
                         c_blend(jp,m0,ip,p_nek_loc)

                 enddo ! i_diff

                 do i_diff=-mg_dx,mg_dx-ig_dx

                    !----------------------------
                    ! This will stagnate at i=n_x  
                    ! for i+i_diff > n_x
                    !
                    ip = i_cyc(i+i_diff)
                    !----------------------------

                    vel_sum_loc(i_diff,j,jp) = vel_sum_loc(i_diff,j,jp)+&
                         fakefield_flag*sum(ion_current(:,i_diff))* &
                         cs_blend(j,m0,i,p_nek_loc)* &
                         c_blend(jp,m0,ip,p_nek_loc)

                 enddo ! i_diff

              enddo ! jp
           enddo ! j

        enddo ! m

     enddo ! p_nek_loc

     call MPI_ALLREDUCE(vel_sum_loc,&
          vel_sum_glob,&
          n_blend*n_blend*(2*mg_dx-ig_dx+1),&
          MPI_DOUBLE_COMPLEX,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

     aa_mm(i,:,:,:) = vel_sum_glob(:,:,:)

  enddo ! i

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[make_ampere_blend done]'
  endif

end subroutine make_ampere_blend
