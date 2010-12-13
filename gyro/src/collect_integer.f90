!-----------------------------------------------------------
! collect_inetger.f90
!
! PURPOSE:
!  Merge an INTEGER array which is distributed over 
!  toroidal subgroups.
!
! NOTES:
!
! input:  v_distrib(n_n_1,n_v)
! output: v_collect(n_n,n_v)
!
! n_v is the length of the array for each n.
!-----------------------------------------------------------

subroutine collect_integer(v_distrib,v_collect)

  use gyro_globals

  !-------------------------------------
  implicit none
  !
  integer, intent(in) :: v_distrib
  integer, intent(inout) :: v_collect(n_n)
  !
  integer :: i_send
  integer :: i_group_send
  !
  integer :: temp_n
  !-------------------------------------

  include 'mpif.h'

  do in=1,n_n

     !-----------------------------------------
     ! Subgroup collector:
     !
     i_group_send = (in-1)/n_n_1

     if (i_group_send /= 0) then

        i_send = i_group_send*n_proc_1

        if (i_proc == 0) then

           call MPI_RECV(temp_n,&
                1,&
                MPI_INTEGER,&
                i_send,&
                in,&
                GYRO_COMM_WORLD,&
                recv_status,&
                i_err)

        else if (i_proc == i_send) then

           call MPI_SEND(v_distrib,&
                1,&
                MPI_INTEGER,&
                0,&
                in,&
                GYRO_COMM_WORLD,&
                i_err)

        endif

     else

        temp_n = v_distrib

     endif
     !-----------------------------------------

     v_collect(in) = temp_n

  enddo ! in

end subroutine collect_integer
