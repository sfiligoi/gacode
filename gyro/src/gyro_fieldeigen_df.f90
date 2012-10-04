!-----------------------------------------------------------------
! gyro_fieldeigen_df.f90
!
! PURPOSE:
!  Compute distribution function given field eigenmode.
!------------------------------------------------------------------

subroutine gyro_fieldeigen_df

  use mpi
  use gyro_globals
  use gyro_pointers
  use math_constants
  use gyro_collision_private, only : nu_total
  use gyro_fieldeigen_private

  !---------------------------------------------------------------
  implicit none
  !---------------------------------------------------

  complex, dimension(:), allocatable :: right_vec
  complex, dimension(:), allocatable :: h_vec

  integer, dimension(n_energy/n_proc,n_lambda) :: p_nek_vec

  a_eigen_loc(:,:) = (0.0,0.0)

  do is=1,n_kinetic
     p_nek_loc = 0
     do ie2=1,n_energy/n_proc

        ie = ie2+i_proc_1*n_energy/n_proc

        do k=1,n_lambda

           p_nek_loc = p_nek_loc+1
           p_nek_vec(ie2,k) = p_nek_loc

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

           allocate(right_vec(n_imk))
           allocate(h_vec(n_imk))

           call ZGEMV('N',n_imk,n_eigen,cmplx_1,gk_rightk,n_imk,&
                b_eigen,1,cmplx_0,right_vec,1)
           call ZGEMV('N',n_imk,n_imk,cmplx_1,propinvk,n_imk,&
                right_vec,1,cmplx_0,h_vec,1)

           do i=1,n_x
              do m=1,n_stack
                 do k=1,n_lambda
                    imk = i+(m-1)*n_x+(k-1)*n_x*n_stack
                    h(m,i,p_nek_vec(ie2,k),is) = h_vec(imk)
                 enddo
              enddo
           enddo

           deallocate(right_vec)
           deallocate(h_vec)

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

              allocate(right_vec(n_im))
              allocate(h_vec(n_im))

              call ZGEMV('N',n_im,n_eigen,cmplx_1,gk_right,n_im,&
                   b_eigen,1,cmplx_0,right_vec,1)
              call ZGEMV('N',n_im,n_im,cmplx_1,propinv,n_im,&
                   right_vec,1,cmplx_0,h_vec,1)

              do i=1,n_x
                 do m=1,n_stack
                    im = i+(m-1)*n_x
                    ! This is really H = cap_h
                    h(m,i,p_nek_vec(ie2,k),is) = h_vec(im)
                 enddo
              enddo

              deallocate(right_vec)
              deallocate(h_vec)

           enddo ! k

        endif ! collision_flag test

     enddo  ! ie2

  enddo ! is

  call gyro_field_interpolation

  p_nek_loc = 0
  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     do m=1,n_stack
        do is=1,n_kinetic
           do i=1,n_x
              ! h = H-z*alpha*<psi>
              h(m,i,p_nek_loc,is) = h(m,i,p_nek_loc,is)-&
                   z(is)*alpha_s(is,i)*gyro_u(m,i,p_nek_loc,is)
           enddo
        enddo
     enddo

  enddo

  h_old = h

end subroutine gyro_fieldeigen_df
