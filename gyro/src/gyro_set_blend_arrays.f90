!-----------------------------------------------------
! gyro_set_blend_arrays.f90 [caller BigScience]
!
! PURPOSE:
!  Precomputed arrays required for field-solve.
!-----------------------------------------------------

subroutine gyro_set_blend_arrays

  use gyro_globals
  use gyro_pointers
  use math_constants

  !---------------------------------------------------
  implicit none
  !
  integer :: j_plot
  !
  real :: x
  real, dimension(n_blend,n_x) :: c_fluxave_loc
  !
  complex, dimension(n_blend,n_blend,2) :: ff_sum_loc
  complex, dimension(n_blend,n_blend,2) :: ff_sum_glob
  complex, dimension(2*n_blend,2*n_blend) :: ff2_sum_loc
  complex, dimension(2*n_blend,2*n_blend) :: ff2_sum_glob
  !
  complex :: f_b
  complex :: f_bp
  complex :: fmm
  !
  complex, external :: BLEND_F
  complex, external :: BLEND_Fp
  !---------------------------------------------------

  include 'mpif.h'

  call BLEND_init(blend_fit_order,n_blend)

  !--------------------------------------------------
  ! Initialize n=0 flux-surface-average coefficients:
  !
  c_fluxave(:,:) = 0.0
  c_fluxave_loc(:,:) = 0.0
  !--------------------------------------------------

  !--------------------------------------------------
  ! Do the ubiquitous sum over velocity space:
  !
  p_nek_loc = 0

  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek)  
     k  = nek_k(p_nek)   

     ck = class(k)

     do i=1,n_x

        do m=1,n_stack

           m0 = m_phys(ck,m)

           ! Get the j-th theta value on orbit 
           ! with lambda=lambda(k)

           x = theta_t(i,k,m)/pi

           do j=1,n_blend

              f_b = BLEND_F(j,x,phase(in_1,i))

              c_blend(j,m0,i,p_nek_loc)  = f_b
              cs_blend(j,m0,i,p_nek_loc) = conjg(f_b)*w_p(ie,i,k,1)

              f_bp = BLEND_Fp(j,x,phase(in_1,i))/pi

              cs_blend_prime(j,m0,i,p_nek_loc) = &
                   conjg(f_bp)*w_p(ie,i,k,1)            

              ! This is meant only for n=0, so take real part:

              if (n_1(in_1) == 0) then
                 c_fluxave_loc(j,i) = &
                      c_fluxave_loc(j,i)+real(f_b*w_p(ie,i,k,1))
              endif

           enddo ! j

        enddo ! m0

     enddo ! i

  enddo ! p_nek
  !--------------------------------------------------

  !--------------------------------------------------
  ! Need to reduce this sum on our subgroup:
  !
  call MPI_ALLREDUCE(c_fluxave_loc,&
       c_fluxave,&
       n_blend*n_x,&
       MPI_DOUBLE_PRECISION,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)
  !--------------------------------------------------

  !------------------------------------------------------------------
  ! Compute blending matrices:
  !
  !    ff_mm(j,jp,1) = FV[ F*_j F_jp ]
  !    ff_mm(j,jp,2) = FV[ F*_j v_para^2 F_jp ]
  !
  !    ff2_mm(j,jp) = FV[ F*_j {1,e,e,e^2} F_jp ]
  !
  is = 1
  !
  do i=1,n_x

     ff_sum_loc  = (0.0,0.0)
     ff2_sum_loc = (0.0,0.0)

     p_nek_loc = 0

     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        ie = nek_e(p_nek)
        k  = nek_k(p_nek)

        ck = class(k)

        ! Now, compute m-projections of the RHS:

        do m=1,n_stack

           m0 = m_phys(ck,m)

           do j=1,n_blend
              do jp=1,n_blend

                 !
                 ! fmm = (F*_j) (F_jp) dv
                 !

                 fmm = cs_blend(j,m0,i,p_nek_loc)* &
                      c_blend(jp,m0,i,p_nek_loc)

                 ff_sum_loc(j,jp,1) = ff_sum_loc(j,jp,1)+fmm

                 ff_sum_loc(j,jp,2) = ff_sum_loc(j,jp,2)+fmm*&
                      v_para(m,i,p_nek_loc,is)**2

                 ff2_sum_loc(j,jp) = &
                      ff2_sum_loc(j,jp)+&
                      fmm

                 ff2_sum_loc(j,jp+n_blend) = &
                      ff2_sum_loc(j,jp+n_blend)+&
                      fmm*energy(ie,is)

                 ff2_sum_loc(j+n_blend,jp) = &
                      ff2_sum_loc(j+n_blend,jp)+&
                      fmm*energy(ie,is)

                 ff2_sum_loc(j+n_blend,jp+n_blend) = &
                      ff2_sum_loc(j+n_blend,jp+n_blend)+&
                      fmm*energy(ie,is)**2


              enddo ! jp
           enddo ! j

        enddo ! m 

     enddo ! p_nek_loc 

     call MPI_ALLREDUCE(ff_sum_loc,&
          ff_sum_glob,&
          size(ff_sum_glob),&
          MPI_DOUBLE_COMPLEX,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

     call MPI_ALLREDUCE(ff2_sum_loc,&
          ff2_sum_glob,&
          size(ff2_sum_glob),&
          MPI_DOUBLE_COMPLEX,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

     ! LU-decomposition of ff's

     call ZGETRF(n_blend,&
          n_blend,&
          ff_sum_glob(:,:,1),&
          n_blend,&
          ff_mm_piv(:,i,1),&
          info)

     ff_mm(:,:,i,1) = ff_sum_glob(:,:,1)

     call ZGETRF(n_blend,&
          n_blend,&
          ff_sum_glob(:,:,2),&
          n_blend,&
          ff_mm_piv(:,i,2),&
          info)

     ff_mm(:,:,i,2) = ff_sum_glob(:,:,2)

     ! LU-decomposition of ff2

     call ZGETRF(2*n_blend,&
          2*n_blend,&
          ff2_sum_glob(:,:),&
          2*n_blend,&
          ff2_mm_piv(:,i),&
          info)

     ff2_mm(:,:,i) = ff2_sum_glob(:,:)

  enddo ! i 
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Store blending plot coefficients
  ! 
  do i=1,n_x
     do j_plot=1,n_theta_plot

        x = theta_plot(j_plot)/pi

        do j=1,n_blend
           blend_plot(j,j_plot,i) = BLEND_F(j,x,phase(in_1,i))
           blend_prime_plot(j,j_plot,i) = BLEND_Fp(j,x,phase(in_1,i))/pi
        enddo ! j

     enddo ! j_plot
  enddo ! i
  if (iohdf5out == 1) then
    do i=1,n_x
       do j_plot=1,n_theta_plot*n_theta_mult
          !x = -pi+REAL(j-1)*pi_2/REAL(n_theta_plot*n_theta_mult)
          x=theta_fine_start+REAL(j_plot-1)*theta_fine_angle/REAL(n_theta_plot*n_theta_mult)
          x=x/pi
          do j=1,n_blend
             blend_fine(j,j_plot,i) = BLEND_F(j,x,phase(in_1,i))
             blend_prime_fine(j,j_plot,i) = BLEND_Fp(j,x,phase(in_1,i))/pi
          enddo ! j

       enddo ! j_plot
    enddo ! i
  endif
  !---------------------------------------------------------------  

  !---------------------------------------------------------------
  ! Store blending plot coefficients
  ! 
  i = ir_norm
  !
  do j_plot=1,field_r0_grid

     x = theta_r0_plot(j_plot)/pi

     do j=1,n_blend
        blend_r0_plot(j,j_plot) = BLEND_F(j,x,phase(in_1,i))
     enddo ! j

  enddo ! j_plot
  !---------------------------------------------------------------  

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_set_blend_arrays done]'
  endif

end subroutine gyro_set_blend_arrays
