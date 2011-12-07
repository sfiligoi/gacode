!---------------------------------
! i_solve=0: factor
! i_solve=1: solve
!---------------------------------

subroutine sparse_solve_umfpack(n_elem,n_row,matnum,i_solve)

  use gyro_globals

  !---------------------------------------------------
  implicit none
  !
  integer, intent(in) :: n_elem
  integer, intent(in) :: n_row
  integer, intent(in) :: matnum
  integer, intent(in) :: i_solve
  !
  integer :: iter
  integer :: ij
  !
  complex, dimension(:), allocatable :: b_UMF
  complex, dimension(:), allocatable :: x_UMF
  complex, dimension(:), allocatable :: w_UMF
  !
  !---------------------------------------------------

  select case(matnum)

  case(1) 
     m_sparse => m_poisson
     indx_sparse => indx_poisson

  case(2)
     m_sparse => m_ampere
     indx_sparse => indx_ampere

  case(3)
     m_sparse => m_maxwell
     indx_sparse => indx_maxwell

  case(4)
     m_sparse => m_poissonaperp
     indx_sparse => indx_poissonaperp

  end select

  if (i_solve == 0) then

     !==============================================
     ! FACTORIZE
     !==============================================

     iter = 0

20   continue

     ! Define sparse field matrices

     select case(matnum)

     case(1) 

        call gyro_build_sparse_poisson

     case(2)

        call gyro_build_sparse_ampere

     case(3)

        call gyro_build_sparse_maxwell

     case(4)

        call gyro_build_sparse_poissonaperp

     end select

     call UMZ21I(keep(matnum,:),cntl(matnum,:),icntl(matnum,:))

     !--------------------------------------------
     ! Error and diagnostic messages from UMFPACK,
     ! to unit=5.
     !
     icntl(matnum,1) = 6
     icntl(matnum,2) = 6
     !
     ! Use icntl(3) = 2 for more verbosity.
     icntl(matnum,3) = 0
     !
     ! Use icntl(8)=n for n levels of iterative 
     ! refinement (set job=1).
     icntl(matnum,8) = 0
     !--------------------------------------------

     ! Factor Maxwell matrices
     call UMZ2FA(n_row,&
          n_elem,&
          0,&
          .false.,&
          lvalue(matnum),&
          lindx(matnum),&
          m_sparse,&
          indx_sparse,&
          keep(matnum,:),&
          cntl(matnum,:),&
          icntl(matnum,:),&
          uinfo,&
          rinfo)

     lindx(matnum)  = uinfo(19)+n_elem
     lvalue(matnum) = uinfo(21)+n_elem

     if (uinfo(1) < 0) then

        iter = iter+1

        select case(matnum)

        case(1) 

           deallocate(m_poisson)
           deallocate(indx_poisson)

           allocate(m_poisson(lvalue(matnum)))
           allocate(indx_poisson(lindx(matnum)))

           m_sparse => m_poisson
           indx_sparse => indx_poisson

        case(2)

           deallocate(m_ampere)
           deallocate(indx_ampere)

           allocate(m_ampere(lvalue(matnum)))
           allocate(indx_ampere(lindx(matnum)))

           m_sparse => m_ampere
           indx_sparse => indx_ampere

        case(3)

           deallocate(m_maxwell)
           deallocate(indx_maxwell)

           allocate(m_maxwell(lvalue(matnum)))
           allocate(indx_maxwell(lindx(matnum)))

           m_sparse => m_maxwell
           indx_sparse => indx_maxwell

        case(4)

           deallocate(m_poissonaperp)
           deallocate(indx_poissonaperp)

           allocate(m_poissonaperp(lvalue(matnum)))
           allocate(indx_poissonaperp(lindx(matnum)))

           m_sparse => m_poissonaperp
           indx_sparse => indx_poissonaperp

        end select

        goto 20

     endif

     if (uinfo(1) > 0) then
        if (i_proc == 0) then 
           print *,'matnum   =',matnum
           print *,'uinfo(1) =',uinfo(1)
        endif
        call catch_error('Error: UMZ2FA')
     endif

     call write_matrix_stat(n_elem,lvalue(matnum),lindx(matnum),iter,matnum)

  else

     !----------
     ! SOLVE
     !----------

     allocate(b_UMF(n_row))
     allocate(x_UMF(n_row))
     allocate(w_UMF(4*n_row))

     select case(matnum)

     case(1)
        do i=1,n_x
           do j=1,n_blend
              ij  = i+(j-1)*n_x
              b_UMF(ij) = vel_sum_p(j,i)
           enddo ! j
        enddo ! i

     case(2) 
        do i=1,n_x
           do j=1,n_blend
              ij = i+(j-1)*n_x
              b_UMF(ij) = vel_sum_a(j,i)
           enddo ! j
        enddo ! i

     case(3) 
        do i=1,n_x
           do j=1,n_blend
              ij = i+(j-1)*n_x
              b_UMF(ij) = vel_sum_p(j,i)
           enddo ! j
        enddo ! i

        if (n_field > 1) then

           do i=1,n_x
              do j=1,n_blend
                 ij = i+(j-1)*n_x + n_x*n_blend
                 b_UMF(ij) = vel_sum_a(j,i)
              enddo ! j
           enddo ! i

        endif

        if (n_field > 2) then

           do i=1,n_x
              do j=1,n_blend
                 ij = i+(j-1)*n_x + 2*n_x*n_blend
                 b_UMF(ij) = vel_sum_aperp(j,i)
              enddo ! j
           enddo ! i

        endif

     case(4)
        do i=1,n_x
           do j=1,n_blend
              ij = i+(j-1)*n_x
              b_UMF(ij) = vel_sum_p(j,i)
           enddo ! j
        enddo ! i        

        do i=1,n_x
           do j=1,n_blend
              ij = i+(j-1)*n_x + n_x*n_blend
              b_UMF(ij) = vel_sum_aperp(j,i)
           enddo ! j
        enddo ! i        

     end select

     call UMZ2SO(n_row,&
          0,&
          .false.,&
          lvalue(matnum),&
          lindx(matnum),&
          m_sparse,&
          indx_sparse,&
          keep(matnum,:),&
          b_UMF,&
          x_UMF,&
          w_UMF,&
          cntl(matnum,:),&
          icntl(matnum,:),&
          uinfo,&
          rinfo)

     if (uinfo(1) /= 0) then
        if (i_proc == 0) then 
           print *,'matnum   =',matnum
           print *,'uinfo(1) =',uinfo(1)
        endif
        call catch_error('Error: UMZ2SO')
     endif

     select case(matnum)

     case(1)
        do i=1,n_x
           do j=1,n_blend
              ij = i+(j-1)*n_x
              field_blend(j,i,1) = x_UMF(ij)
           enddo ! j
        enddo ! i

     case(2) 
        do i=1,n_x
           do j=1,n_blend
              ij = i+(j-1)*n_x
              field_blend(j,i,2) = x_UMF(ij)
           enddo
        enddo

     case(3) 
        do i=1,n_x
           do j=1,n_blend
              ij = i+(j-1)*n_x
              field_blend(j,i,1) = x_UMF(ij)
           enddo
        enddo

        if (n_field > 1) then
           do i=1,n_x
              do j=1,n_blend
                 ij = i+(j-1)*n_x + n_x*n_blend
                 field_blend(j,i,2) = x_UMF(ij)
              enddo
           enddo
        endif

        if (n_field > 2) then
           do i=1,n_x
              do j=1,n_blend
                 ij = i+(j-1)*n_x + 2*n_x*n_blend
                 field_blend(j,i,3) = x_UMF(ij)
              enddo
           enddo
        endif

     case(4)
        do i=1,n_x
           do j=1,n_blend
              ij = i+(j-1)*n_x
              field_blend(j,i,1) = x_UMF(ij)
           enddo
        enddo

        do i=1,n_x
           do j=1,n_blend
              ij = i+(j-1)*n_x + n_x*n_blend
              field_blend(j,i,3) = x_UMF(ij)
           enddo
        enddo

     end select

     deallocate(b_UMF)
     deallocate(x_UMF)
     deallocate(w_UMF)

  endif

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[sparse_solve done]'
  endif

end subroutine sparse_solve_umfpack
