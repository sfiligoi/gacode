!-----------------------------------------------------------------
! gyro_fieldeigen_kernel.f90
!
! PURPOSE:
!  Kernel of the Maxwell dispersion matrix eigenvalue problem.
!  Compute the determinant of the dispersion matrix.
!------------------------------------------------------------------

subroutine gyro_fieldeigen_kernel

  use gyro_globals
  use gyro_pointers
  use math_constants
  use gyro_collision_private, only : nu_total
  use gyro_fieldeigen_private

  !---------------------------------------------------------------
  implicit none
  !---------------------------------------------------

  include 'mpif.h'

  a_eigen_loc(:,:) = (0.0,0.0)

  do is=1,n_kinetic
     p_nek_loc = 0
     do ie2=1,n_energy/n_proc

        ie = ie2+i_proc_1*n_energy/n_proc

        do k=1,n_lambda

           p_nek_loc = p_nek_loc+1

           ck = class(k)

           !-------------------------------------------------------------
           ! v_omegastar = -omega_star
           !
           do i=1,n_x

              do m=1,n_stack

                 mk = m + (k-1)*n_stack 

                 v_omegastar(i,mk) = -i_c*omega_star(m,i,p_nek_loc,is)

              enddo ! m

           enddo ! i
           !-------------------------------------------------------------

           !-------------------------------------------------------------
           ! fg = conjg(F_blend) * Gyroavgop
           !
           do m=1,n_stack

              m0 = m_phys(ck,m)
              mk = m+(k-1)*n_stack 

              do i=1,n_x
                 do i_diff=-m_gyro,m_gyro-i_gyro

                    fg(:,i_diff,i,mk,1) = w_gyro(m0,i_diff,i,p_nek_loc,is) &
                         *conjg(c_blend(:,m0,i,p_nek_loc))

                    if (n_field > 1) then
                       fg(:,i_diff,i,mk,2) = w_gyro(m0,i_diff,i,p_nek_loc,is) &
                            *(-v_para(m,i,p_nek_loc,is))  &
                            *conjg(c_blend(:,m0,i,p_nek_loc))
                    endif

                    if (n_field > 2) then
                       fg(:,i_diff,i,mk,3) = 2.0*energy(ie,is)*lambda(i,k) &
                            *tem_s(is,i)/z(is) &
                            *w_gyro_aperp(m0,i_diff,i,p_nek_loc,is) &
                            *conjg(c_blend(:,m0,i,p_nek_loc))

                    endif
                 enddo
              enddo
           enddo ! m
           !-------------------------------------------------------------

        enddo ! k

        if (collision_flag == 1 .and. is == indx_e) then

           propinvk(:,:) = (0.0,0.0)

           ! Dense-in-k (collisions)

           !------------------------------------------------------------
           ! Set propagator elements
           !
           ! Diagonal in k
           !
           do i=1,n_x
              do m=1,n_stack
                 do k=1,n_lambda
                    im = i+(m-1)*n_x
                    imk = i+(m-1)*n_x+(k-1)*n_x*n_stack
                    do ip=1,n_x
                       do mp=1,n_stack
                          impr = ip+(mp-1)*n_x
                          imkp = ip+(mp-1)*n_x+(k-1)*n_x*n_stack

                          propinvk(imk,imkp) = vdotgrad(im,impr,k,ie2,is)

                       enddo
                    enddo
                 enddo
              enddo
           enddo
           !------------------------------------------------------------

           !------------------------------------------------------------
           ! Diagonal in i
           !
           do i=1,n_x
              do m=1,n_stack
                 do k=1,n_lambda
                    mk = m+(k-1)*n_stack
                    imk = i+(m-1)*n_x+(k-1)*n_x*n_stack
                    do mp=1,n_stack
                       do kp=1,n_lambda
                          mkp = mp+(kp-1)*n_stack
                          imkp = i+(mp-1)*n_x+(kp-1)*n_x*n_stack

                          propinvk(imk,imkp) =  propinvk(imk,imkp) &
                               -nu_op(mkp,mk,i,ie2)

                       enddo
                    enddo
                 enddo
              enddo
           enddo
           !------------------------------------------------------------

           !------------------------------------------------------------
           ! Diagonal in imk
           !
           do imk=1,n_imk
              propinvk(imk,imk) = propinvk(imk,imk)-i_c*omega_eigen
           enddo
           !-----------------------------------------------------------

           !-----------------------------------------------------------
           ! Invert the propagator: final result is "propinv"
           ! (these are by far the most expensive operations)
           !
           call ZGETRF(n_imk,n_imk,propinvk,n_imk,i_piv_imk,info_imk)
           call ZGETRI(n_imk,propinvk,n_imk,i_piv_imk,work_imk,n_imk,info_imk)
           !-----------------------------------------------------------

           gk_rightk(:,:) = 0.0
           do ixp = 1,n_field
              do ip=1,n_x
                 do jp=1,n_blend
                    ijp = ip + (jp-1)*n_x + (ixp-1)*n_x*n_blend
                    do mp=1,n_stack
                       do kp=1,n_lambda
                          mkp = mp + (kp-1)*n_stack
                          do idiff_ipp=-m_gyro,m_gyro-i_gyro

                             ipp  = i_cyc(ip+idiff_ipp)
                             imkp = ipp+n_x*(mp-1)+n_x*n_stack*(kp-1)

                             gk_rightk(imkp,ijp) = gk_rightk(imkp,ijp) &
                                  +(-i_c*omega_eigen*z(is)*alpha_s(is,ipp) &
                                  + v_omegastar(ipp,mkp)) &
                                  * conjg(fg(jp,idiff_ipp,ip,mkp,ixp))

                          enddo
                       enddo
                    enddo
                 enddo
              enddo
           enddo ! ixp


           gk_leftk(:,:) = 0.0
           do ix=1,n_field
              do i=1,n_x
                 do j=1,n_blend
                    ij = i + (j-1)*n_x + (ix-1)*n_x*n_blend
                    do m=1,n_stack
                       do k=1,n_lambda
                          mk = m + (k-1)*n_stack
                          do idiff_ippp=-m_gyro,m_gyro-i_gyro

                             ippp = i_cyc(i+idiff_ippp)
                             imk  = ippp+n_x*(m-1)+n_x*n_stack*(k-1)

                             gk_leftk(ij,imk) = gk_leftk(ij,imk) &
                                  +z(is)*w_p(ie,i,k,is)  &
                                  *fg(j,idiff_ippp,i,mk,ix) 

                          enddo
                       enddo
                    enddo
                 enddo
              enddo
           enddo ! ix

           !-----------------------------------------------------------
           ! Matrix triple-product:  gk_leftk * propinvk * gk_rightk
           !
           call ZGEMM('N','N',n_eigen,n_imk,n_imk,cmplx_1,&
                gk_leftk,n_eigen,propinvk,n_imk,cmplx_0,prodk,n_eigen)
           call ZGEMM('N','N',n_eigen,n_eigen,n_imk,cmplx_1,&
                prodk,n_eigen,gk_rightk,n_imk,cmplx_0,a_eigen,n_eigen)
           !-----------------------------------------------------------

           a_eigen_loc(:,:) = a_eigen_loc(:,:)+a_eigen(:,:)

        else

           ! Diagonal-in-k (no collisions)

           do k=1,n_lambda

              !-----------------------------------------------------------
              ! Finalize definition of propagator
              !
              propinv(:,:) = vdotgrad(:,:,k,ie2,is)

              do im=1,n_im
                 propinv(im,im) = propinv(im,im)-i_c*omega_eigen
              enddo
              !-----------------------------------------------------------

              !-----------------------------------------------------------
              ! Invert the propagator: final result is "propinv"
              !
              call ZGETRF(n_im,n_im,propinv,n_im,i_piv_im,info_im)
              call ZGETRI(n_im,propinv,n_im,i_piv_im,work_im,n_im,info_im)
              !-----------------------------------------------------------

              gk_right(:,:) = 0.0
              do ixp = 1,n_field
                 do ip=1,n_x
                    do jp=1,n_blend
                       ijp = ip + (jp-1)*n_x + (ixp-1)*n_x*n_blend
                       do mp=1,n_stack
                          mkp = mp + (k-1)*n_stack
                          do idiff_ipp=-m_gyro,m_gyro-i_gyro

                             ipp  = i_cyc(ip+idiff_ipp)
                             impr = ipp + n_x*(mp-1)

                             gk_right(impr,ijp) = gk_right(impr,ijp) &
                                  +(-i_c*omega_eigen*z(is)*alpha_s(is,ipp) &
                                  + v_omegastar(ipp,mkp)) &
                                  * conjg(fg(jp,idiff_ipp,ip,mkp,ixp))

                          enddo
                       enddo
                    enddo
                 enddo
              enddo ! ixp

              gk_left(:,:) = 0.0
              do ix=1,n_field
                 do i=1,n_x
                    do j=1,n_blend
                       ij = i + (j-1)*n_x + (ix-1)*n_x*n_blend
                       do m=1,n_stack
                          mk = m + (k-1)*n_stack
                          do idiff_ippp=-m_gyro,m_gyro-i_gyro

                             ippp = i_cyc(i+idiff_ippp)
                             im  = ippp + n_x*(m-1)

                             gk_left(ij,im) = gk_left(ij,im) &
                                  +z(is)*w_p(ie,i,k,is)  &
                                  *fg(j,idiff_ippp,i,mk,ix) 

                          enddo
                       enddo
                    enddo
                 enddo
              enddo ! ix

              !-----------------------------------------------------------
              ! Matrix triple-product:  gk_left * propinv * gk_right
              !
              call ZGEMM('N','N',n_eigen,n_im,n_im,cmplx_1,&
                   gk_left,n_eigen,propinv,n_im,cmplx_0,prod,n_eigen)
              call ZGEMM('N','N',n_eigen,n_eigen,n_im,cmplx_1,&
                   prod,n_eigen,gk_right,n_im,cmplx_0,a_eigen,n_eigen)
              !-----------------------------------------------------------

              a_eigen_loc(:,:) = a_eigen_loc(:,:)+a_eigen(:,:)

           enddo ! k

        endif ! collision_flag test

     enddo  ! ie2

  enddo ! is

  ! Final velocity-space sum on subgroup
  !
  call MPI_ALLREDUCE(a_eigen_loc(:,:),&
       a_eigen(:,:),&
       size(a_eigen),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  ! Sign flip
  a_eigen = -a_eigen

  do i=1,n_x
     do j=1,n_blend
        do jp=1,n_blend

           ! h-independent part of field equations

           ! Poisson eqn
           ix = 1
           ixp = 1
           do i_diff=-m_gyro,m_gyro-i_gyro
              ip = i_cyc(i+i_diff)

              ijp = ip + (jp-1)*n_x + (ixp-1)*n_x*n_blend 
              ij  = i + (j-1)*n_x + (ix-1)*n_x*n_blend 

              a_eigen(ij,ijp) = a_eigen(ij,ijp) + ap_mm(i,i_diff,j,jp)

           enddo

           ! Ampere eqn
           if (n_field > 1) then
              ix = 2
              ixp = 2
              do i_diff=-mg_dx,mg_dx-ig_dx
                 ip = i_cyc(i+i_diff)

                 ij  = i + (j-1)*n_x + (ix-1)*n_x*n_blend 
                 ijp = ip + (jp-1)*n_x + (ixp-1)*n_x*n_blend 

                 a_eigen(ij,ijp) = a_eigen(ij,ijp) - aa_mm(i,i_diff,j,jp)

              enddo
           endif

           ! Ampere Perp eqn
           if (n_field > 2) then
              ix = 3
              ixp = 3
              do i_diff=-m_gyro,m_gyro-i_gyro
                 ip = i_cyc(i+i_diff)

                 ij  = i + (j-1)*n_x + (ix-1)*n_x*n_blend 
                 ijp = ip + (jp-1)*n_x + (ixp-1)*n_x*n_blend 

                 a_eigen(ij,ijp) = a_eigen(ij,ijp) - 2.0*ab_mm(i,i_diff,j,jp)

              enddo
           endif

        enddo
     enddo
  enddo

  !-------------------------------------------------------------
  ! Determinant calculation
  !
  ! (1) LAPACK factorization
  !
  a_eigen_loc = a_eigen
  call ZGETRF(n_eigen,n_eigen,a_eigen,n_eigen,i_piv_eigen,info_eigen)

  !
  ! (2) Compute the determinant by diagonal product 
  !
  !          det(a_eigen) = det(L) * det(U) * det(P) 
  !
  ! where det(L) = 1 (L has 1's on its diagonal)
  !       det(U) = product of elements on diagonal of LU factorized matrix
  !       det(P) = {+1 if #row permutations even, 
  !                 -1 if #row permutations odd}

  if (diag_scale(1) < 0.0) then
     do i=1,n_eigen
        diag_scale(i) = abs(a_eigen(i,i))
     enddo
  endif

  det = (1.0,0.0)
  do i=1,n_eigen
     if (i_piv_eigen(i) /= i) then
        sgn = -1.0
     else
        sgn = 1.0
     endif
     det = det*a_eigen(i,i)/diag_scale(i)*sgn
  enddo

end subroutine gyro_fieldeigen_kernel
