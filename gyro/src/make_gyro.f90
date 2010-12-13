!-----------------------------------------------------
! make_gyro.f90 [caller: make_radial_operators]
!
! PURPOSE:
!  Driver for creation of (tau-space) gyroaverage 
!  stencils w_gyro and w_rot_gyro (optionally).
!---------------------------------------------

subroutine make_gyro

  use gyro_globals
  use gyro_pointers
  use math_constants

  !----------------------------------------------
  implicit none
  !
  integer :: p
  !
  real :: rho_gyro
  real :: a_gyro
  real :: u_gyro
  real :: v_gyro
  !
  real :: omega_c
  !
  complex :: temp
  !
  complex, dimension(n_theta(2),-m_gyro:m_gyro-i_gyro) :: w_temp
  complex, dimension(n_theta(2),-m_gyro:m_gyro-i_gyro) :: w_temp_rot
  complex, dimension(n_theta(2),-m_gyro:m_gyro-i_gyro) :: w_temp_aperp
  complex, dimension(-m_gyro:m_gyro-i_gyro) :: g
  !----------------------------------------------

  !---------------------------------------------
  ! Set dimensions and functions
  !
  do p=-n_x/2,n_x/2-1
     do m=-m_gyro,m_gyro-i_gyro
        z_gyro(m,p) = exp(-i_c*p*m*pi/(n_x/2))
     enddo
  enddo
  !---------------------------------------------

  w_gyro(:,:,:,:,:) = (0.0,0.0)

  !---------------------------------------------------------
  ! LOCAL gyroaverage:
  !
  do is=1,n_gk
     p_nek_loc = 0
     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        ie = nek_e(p_nek)  
        k  = nek_k(p_nek)   

        ck = class(k)

        do i=1,n_x

           do m0=1,n_theta(ck)

              !---------------------------------------------------------------
              ! Prepare argument of Bessel function
              !  
              omega_c = abs(z(is))*b_unit_s(i)*mu(is)**2

              if (kill_gyro_b_flag == 0) then
                 omega_c = omega_c*b0_t(i,k,m0)
              endif

              rho_gyro = rhos_norm*v_perp(m0,i,p_nek_loc,is)/omega_c
              !
              a_gyro = grad_r_t(i,k,m0)/x_length*dr_eodr(i)
              u_gyro = qrat_t(i,k,m0)*n_1(in_1)*q_s(i)/r_s(i)*captheta_t(i,k,m0)
              v_gyro = qrat_t(i,k,m0)*n_1(in_1)*q_s(i)/r_s(i)
              !----------------------------------------------------------------

              !---------------------------
              ! Standard gyroaverage: G
              !---------------------------

              call gyro_ave(rho_gyro,&
                   a_gyro,&
                   u_gyro,&
                   v_gyro,&
                   g,&
                   1)

              if (n_1(in_1) == 0) then

                 if (i_gyro /= 1) then 

                    ! Renormalize truncated gyroaverage
                    temp = sum(g(:))-1.0
                    g(0) = g(0)-temp

                 endif

                 ! Enforce EXACT reality 

                 g(:) = real(g(:))

              endif

              w_temp(m0,:) = g(:)

              !-------------------------------
              ! Momentum gyroaverage: cos(a) G
              !-------------------------------

              call gyro_ave_rot(rho_gyro,&
                   a_gyro,&
                   u_gyro,&
                   v_gyro,&
                   g)

              w_temp_rot(m0,:) = g(:)

              !-------------------------------
              ! A_perp gyroaverage: (1/2)(J0+J2)
              !-------------------------------

              if (n_field == 3) then

                 call gyro_ave_aperp(rho_gyro,&
                      a_gyro,&
                      u_gyro,&
                      v_gyro,&
                      g,1)

                 if (n_1(in_1) == 0) then

                    ! Enforce EXACT reality 

                    g(:) = real(g(:))

                 endif

                 w_temp_aperp(m0,:) = g(:)

              endif

           enddo ! m0

           do m=1,n_stack
              m0 = m_phys(ck,m)
              w_gyro(m,:,i,p_nek_loc,is) = w_temp(m0,:)
              w_gyro_rot(m,:,i,p_nek_loc,is) = w_temp_rot(m0,:)
              if (n_field == 3) w_gyro_aperp(m,:,i,p_nek_loc,is) = w_temp_aperp(m0,:)
           enddo

        enddo ! i
     enddo ! p_nek
  enddo ! is
  !
  !----------------------------------------------------------



  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[make_gyro done]'
  endif

end subroutine make_gyro
