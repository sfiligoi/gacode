!------------------------------------------------------------------
! gyro_sparse_solve_umfpack.f90
!
! PURPOSE:
!  Manage factorize or solve for various field matrix combinations 
!  using UMFPACK.
!
! NOTES:
!  i_solve=0: factor
!  i_solve=1: solve
!------------------------------------------------------------------

subroutine gyro_sparse_solve_umfpack(n_elem,n_row,matnum,i_solve)

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
        call catch_error('ERROR: (GYRO) UMZ2SO failed.')
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
     print *,'[sparse_solve_umfpack done]'
  endif

end subroutine gyro_sparse_solve_umfpack


subroutine write_matrix_stat(nelem,nval,nindx,niter,tag)

  use gyro_globals

  !------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nelem
  integer, intent(in) :: nval
  integer, intent(in) :: nindx
  integer, intent(in) :: niter
  !
  integer :: nelem_c(n_n)
  integer :: nval_c(n_n)
  integer :: nindx_c(n_n)
  integer :: niter_c(n_n)
  integer :: p
  !
  integer, intent(in) :: tag
  !
  character (len=36) :: lab1
  character (len=36) :: lab2
  character (len=36) :: lab3
  character (len=36) :: lab4
  !------------------------------------------------

  select case (output_flag)

  case (1)

     call collect_integer(nelem,nelem_c)
     call collect_integer(nval,nval_c)
     call collect_integer(nindx,nindx_c)
     call collect_integer(niter,niter_c)

     if (n(1) == 0 .and. n_n > 1) then
        p = 2
        lab1 = 'EXPLICIT POISSON:     n=0        n>0'
        lab2 = ' EXPLICIT AMPERE:     n=0        n>0'
        lab3 = '     TOTAL FIELD:     n=0        n>0'
        lab4 = '  POISSON-AMPERE:     n=0        n>0'
     endif

     if (n(1) == 0 .and. n_n == 1) then
        p = 1
        lab1 = 'EXPLICIT POISSON:     n=0'
        lab2 = ' EXPLICIT AMPERE:     n=0'
        lab3 = '     TOTAL FIELD:     n=0'
        lab4 = '  POISSON-AMPERE:     n=0'
     endif

     if (n(1) /= 0) then
        p = 1
        lab1 = 'EXPLICIT POISSON:     n>0'
        lab2 = ' EXPLICIT AMPERE:     n>0'
        lab3 = '     TOTAL FIELD:     n>0'
        lab4 = '  POISSON-AMPERE:     n>0'
     endif

     if ((i_proc == 0) .and. (gkeigen_j_set == 0)) then
        open(unit=1,file=trim(runfile),status='old',position='append')

        select case (tag)

        case (1)

           write(1,*) '----------- SPARSE MATRIX STATS ---------------'
           write(1,*) lab1

        case (2)

           write(1,*) lab2

        case (3)

           write(1,*) lab3

        case (4)

           write(1,*) '----------- SPARSE MATRIX STATS ---------------'
           write(1,*) lab4

        end select

        write(1,40) '        nonzeros:',nelem_c(1:p)
        write(1,40) '          values:',nval_c(1:p)
        write(1,40) '         indices:',nindx_c(1:p)
        write(1,40) '      iterations:',niter_c(1:p)
        write(1,*)

        close(1)

     endif

     return

  end select

40 format(t2,a,t20,4(i8,2x))

end subroutine write_matrix_stat
