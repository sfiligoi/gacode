!---------------------------------
! i_solve=0: factor
! i_solve=1: solve
!---------------------------------

subroutine sparse_solve_mumps(n_elem,n_row,matnum,i_solve)

  use mpi
  use gyro_globals
  use gyro_mumps_private

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
  !---------------------------------------------------

  call proc_time(CPU_field2_in)

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

     !----------
     ! FACTORIZE
     !----------

     iter = 1

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

     ! initialize MUMPS
     m_mumps(matnum)%COMM = MUMPS_COMM
     m_mumps(matnum)%JOB = -1  ! only initialize on call to ZMUMPS
     m_mumps(matnum)%SYM = 0
     m_mumps(matnum)%PAR = 1

     call ZMUMPS(m_mumps(matnum))

     ! allocate space for MUMPS arrays only on MYID = 0
     if (m_mumps(matnum)%MYID == 0) then

        m_mumps(matnum)%N = n_row
        m_mumps(matnum)%NZ = n_elem
        allocate(m_mumps(matnum)%A(n_elem))
        allocate(m_mumps(matnum)%IRN(n_elem))
        allocate(m_mumps(matnum)%JCN(n_elem))
        allocate(m_mumps(matnum)%RHS(n_row))

        ! copy data from general pointer into mumps structure
        m_mumps(matnum)%A = m_sparse(1:n_elem)
        m_mumps(matnum)%IRN = indx_sparse(1:n_elem)
        m_mumps(matnum)%JCN = indx_sparse(n_elem+1:2*n_elem)

     endif

     !-----------------------------------------
     ! Turn off diagnostic messages from MUMPS
     !
     m_mumps(matnum)%ICNTL(2) = 0
     m_mumps(matnum)%ICNTL(3) = 0
     m_mumps(matnum)%ICNTL(4) = 0
     !
     !--------------------------------------------

     ! Factor Maxwell matrices
     m_mumps(matnum)%JOB = 4  
     call ZMUMPS(m_mumps(matnum))

     if (m_mumps(matnum)%MYID == 0) then
        if (m_mumps(matnum)%INFOG(1) < 0 ) then
           print *,m_mumps(matnum)%INFOG(1)
           call catch_error('Exited ZMUMPS factorization with errors')
        endif
     endif

  else

     !----------
     ! SOLVE
     !----------

     call proc_time(CPU_field2_in)

     if (m_mumps(matnum)%MYID == 0) then

        select case(matnum)

        case(1)
           do i=1,n_x
              do j=1,n_blend
                 ij  = i+(j-1)*n_x
                 m_mumps(matnum)%RHS(ij) = vel_sum_p(j,i)
              enddo ! j
           enddo ! i

        case(2) 
           do i=1,n_x
              do j=1,n_blend
                 ij = i+(j-1)*n_x
                 m_mumps(matnum)%RHS(ij) = vel_sum_a(j,i)
              enddo ! j
           enddo ! i

        case(3) 
           do i=1,n_x
              do j=1,n_blend
                 ij = i+(j-1)*n_x
                 m_mumps(matnum)%RHS(ij) = vel_sum_p(j,i)
              enddo ! j
           enddo ! i

           if (n_field > 1) then
              do i=1,n_x
                 do j=1,n_blend
                    ij = i+(j-1)*n_x + n_x*n_blend
                    m_mumps(matnum)%RHS(ij) = vel_sum_a(j,i)
                 enddo ! j
              enddo ! i
           endif

           if (n_field > 2) then
              do i=1,n_x
                 do j=1,n_blend
                    ij = i+(j-1)*n_x + 2*n_x*n_blend
                    m_mumps(matnum)%RHS(ij) = vel_sum_aperp(j,i)
                 enddo ! j
              enddo ! i
           endif

        case(4)
           do i=1,n_x
              do j=1,n_blend
                 ij = i+(j-1)*n_x
                 m_mumps(matnum)%RHS(ij) = vel_sum_p(j,i)
              enddo ! j
           enddo ! i        

           do i=1,n_x
              do j=1,n_blend
                 ij = i+(j-1)*n_x + n_x*n_blend
                 m_mumps(matnum)%RHS(ij) = vel_sum_aperp(j,i)
              enddo ! j
           enddo ! i     

        end select

     endif

     m_mumps(matnum)%JOB = 3  ! solve only
     call ZMUMPS(m_mumps(matnum))

     if (m_mumps(matnum)%MYID == 0) then
        if (m_mumps(matnum)%INFOG(1) < 0 ) then
           print *,m_mumps(matnum)%INFOG(1)
           call catch_error('error in ZMUMPS solve')
        endif
     endif

     select case(matnum)

     case(1)
        if (m_mumps(matnum)%MYID == 0) then
           do i=1,n_x
              do j=1,n_blend
                 ij = i+(j-1)*n_x
                 field_blend(j,i,1) = m_mumps(matnum)%RHS(ij)
              enddo ! j
           enddo ! i
        endif
        call MPI_BCAST(field_blend(:,:,1),n_x*n_blend,MPI_DOUBLE_COMPLEX,0,NEW_COMM_1,i_err)

     case(2) 
        if (m_mumps(matnum)%MYID == 0) then
           do i=1,n_x
              do j=1,n_blend
                 ij = i+(j-1)*n_x
                 field_blend(j,i,2) = m_mumps(matnum)%RHS(ij)
              enddo
           enddo
        endif
        call MPI_BCAST(field_blend(:,:,2),n_x*n_blend,MPI_DOUBLE_COMPLEX,0,NEW_COMM_1,i_err)

     case(3) 
        if (m_mumps(matnum)%MYID == 0) then
           do i=1,n_x
              do j=1,n_blend
                 ij = i+(j-1)*n_x
                 field_blend(j,i,1) = m_mumps(matnum)%RHS(ij)
              enddo
           enddo

           if (n_field > 1) then
              do i=1,n_x
                 do j=1,n_blend
                    ij = i+(j-1)*n_x + n_x*n_blend
                    field_blend(j,i,2) = m_mumps(matnum)%RHS(ij)
                 enddo
              enddo
           endif

           if (n_field > 2) then
              do i=1,n_x
                 do j=1,n_blend
                    ij = i+(j-1)*n_x + 2*n_x*n_blend
                    field_blend(j,i,3) = m_mumps(matnum)%RHS(ij)
                 enddo
              enddo
           endif

        endif
        call MPI_BCAST(field_blend,size(field_blend),MPI_DOUBLE_COMPLEX,0,NEW_COMM_1,i_err)

     case(4)
        if (m_mumps(matnum)%MYID == 0) then
           do i=1,n_x
              do j=1,n_blend
                 ij = i+(j-1)*n_x
                 field_blend(j,i,1) = m_mumps(matnum)%RHS(ij)
              enddo
           enddo

           do i=1,n_x
              do j=1,n_blend
                 ij = i+(j-1)*n_x + n_x*n_blend
                 field_blend(j,i,3) = m_mumps(matnum)%RHS(ij)
              enddo
           enddo
        endif

        call MPI_BCAST(field_blend(:,:,1),n_x*n_blend,MPI_DOUBLE_COMPLEX,0,NEW_COMM_1,i_err)
        call MPI_BCAST(field_blend(:,:,3),n_x*n_blend,MPI_DOUBLE_COMPLEX,0,NEW_COMM_1,i_err)

     end select


  endif

  call proc_time(CPU_field2_out)
  CPU_field2 = CPU_field2 + (CPU_field2_out - CPU_field2_in)

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[sparse_solve done]'
  endif

end subroutine sparse_solve_mumps
