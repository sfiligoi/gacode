subroutine SOLVE_superlu(n_elem, n_size, a, a_iindx, a_jindx)
    use neo_globals
    implicit none
    integer, intent (in) :: n_elem, n_size
    real, dimension(:), intent (in) :: a
    integer, dimension(:), intent (in) :: a_iindx
    integer, dimension(:), intent (in) :: a_jindx

    real, dimension(:), allocatable :: acc
    real, dimension(:), allocatable :: b
    integer, dimension(:), allocatable :: rowindx
    integer, dimension(:), allocatable :: colindx

    integer :: n, nnz, nrhs=1
    integer :: num_threads

     num_threads = n_omp

     !Get the number of non-zero elements in compressed column format
     call st_to_cc_size ( n_elem, a_iindx, a_jindx, nnz )

     !Allocate arrays to hold compressed column format
     n = n_row
     allocate(acc(1:nnz))
     allocate(rowindx(1:nnz))
     allocate(colindx(1:n+1))
     allocate(b(1:n))
     b(1:n) = g(1:n)

     !Convert sparse triplet to compressed column format 
     call st_to_cc_index ( n_elem, a_iindx, a_jindx, nnz, n_row, rowindx, colindx )
     call st_to_cc_values ( n_elem, a_iindx, a_jindx, a, nnz, n_row, rowindx, colindx, acc )

     ! Factor the Matrix -- uses acc(:) and [row/col]indx(:)
     if(silent_flag == 0 .and. i_proc == 0) then
        open(unit=io_neoout,file=trim(path)//runfile_neoout,&
             status='old',position='append')
        write(io_neoout,'(t2,a)') 'Begin matrix factor and solve -- SuperLU_MT'
        close(io_neoout)
     endif

     call c_bridge_pdgssv(num_threads, n, nnz, nrhs, acc, rowindx, colindx, b, n, error_status)
     if (error_status .ne. 0) return

     if(silent_flag == 0 .and. i_proc == 0) then
        open(unit=io_neoout,file=trim(path)//runfile_neoout,&
             status='old',position='append')
        write(io_neoout,'(t2,a)') 'Done matrix factor and solve -- SuperLU_MT'
        close(io_neoout)
     endif

     g(1:n) = b(1:n)

     !Last, free the storage allocated inside SuperLU
     deallocate(acc)
     deallocate(b)
     deallocate(rowindx)
     deallocate(colindx)

end subroutine SOLVE_superlu

