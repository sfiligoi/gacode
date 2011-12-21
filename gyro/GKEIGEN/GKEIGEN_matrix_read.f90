!-----------------------------------------------------
! GKEIGEN_matrix_read.f90
!
! PURPOSE:
!  Read all columns of the linear time-derivative
!  matrix prevously written with write_eigensolve_matrix.f90.
!  Applicable only in eigensolve mode.
!-----------------------------------------------------

subroutine GKEIGEN_matrix_read(n_last_col)

  use gyro_globals
  use gyro_pointers
  use GKEIGEN_globals

  !---------------------------------
  implicit none
  !
  integer :: n_last_col
  integer :: ierr
  !
  real :: real_entry
  real :: imag_entry
  !
  complex, dimension(h_length_loc) :: one_col_r
  !
  character :: dummy
  !---------------------------------

  include 'mpif.h'

10 format(1x,'(',g14.6,',',g14.6,')')
20 format('----------------- Column: ',I6,'  ----------------------')

  n_last_col = 0

  open(unit=1,file=trim(path)//file_eigen_restart,status='old')

  read (1,*) dummy

  do istate_1=0,h_length-1

     read (1,20,end=100) n_last_col
     do i_proc_w = 0, n_proc_1-1

        if (i_proc == i_proc_w) Then

           do jstate=1,h_length_loc
              read (1,10,advance='no') real_entry, imag_entry
              one_col_r(jstate) = real_entry*(1.0,0.0) + imag_entry*(0.0,1.0)
           enddo
           eigensolve_matrix(:,istate_1+1) = one_col_r(:)
           read (1,*) dummy

        else

           read (1,*) dummy

        endif

        call MPI_BARRIER(GYRO_COMM_WORLD,ierr)

     enddo ! i_proc_w

  enddo ! istate_1

100 close (unit=1)

end subroutine GKEIGEN_matrix_read

