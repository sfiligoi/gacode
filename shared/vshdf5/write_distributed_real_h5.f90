!------------------------------------------------------
! write_distributed_real.f90
!
! PURPOSE:
!  Control merged output of distributed real array.
!------------------------------------------------------

subroutine write_distributed_real_h5(varName,rGid,n_fn,fn,h5in,h5err)

  use hdf5_api
  use gyro_globals, only : &
       n_n,&
       n_n_1,&
       n_proc_1,&
       recv_status,&
       data_step,&
       GYRO_COMM_WORLD,&
       i_proc,&
       i_err

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: varName
  integer(HID_T), intent(in) :: rGid
  integer, intent(in) :: n_fn
  complex, intent(in) :: fn(n_fn)
  type(hdf5InOpts), intent(inout) :: h5in
  type(hdf5ErrorType), intent(inout) :: h5err
  !
  integer :: data_loop
  integer :: i_group_send
  integer :: i_send
  integer :: in
  !
  real :: fn_recv(n_fn)
  character(128) :: n_varName
  character(3) :: n_name
  !------------------------------------------------------


  include 'mpif.h'

     do in=1,n_n

        !-----------------------------------------
        ! Subgroup collector:
        !
        i_group_send = (in-1)/n_n_1

        if (i_group_send /= 0) then

           i_send = i_group_send*n_proc_1

           if (i_proc == 0) then

              call MPI_RECV(fn_recv,&
                   n_fn,&
                   MPI_DOUBLE_PRECISION,&
                   i_send,&
                   in,&
                   GYRO_COMM_WORLD,&
                   recv_status,&
                   i_err)

           else if (i_proc == i_send) then

              call MPI_SEND(fn,&
                   n_fn,&
                   MPI_DOUBLE_PRECISION,&
                   0,&
                   in,&
                   GYRO_COMM_WORLD,&
                   i_err)

           endif

        else

           fn_recv(:) = fn(:)

        endif
        !
        !-----------------------------------------

        WRITE(n_name,fmt='(i3.3)') in
        n_varName=trim(varName)//"_"//n_name
        if (i_proc == 0) call dump_h5(rGid,n_varName,fn_recv,h5in,h5err)

     enddo ! in

end subroutine write_distributed_real_h5
