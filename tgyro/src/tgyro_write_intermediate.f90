subroutine tgyro_write_intermediate(index,res_in)

  use mpi
  use tgyro_globals

  implicit none

  integer, intent(in) :: index
  real, intent(in) :: res_in(p_max)
  integer :: ip


  if (i_proc_global == 0) then
     open(unit=1,file='out.tgyro.iterate',status='old',position='append')

     if (index == 0) write(1,'(t2,a,i3,(2x,1pe12.5))') 'ITERATION: ',i_tran

     write(1,10) &
          index,sum(res_in)/size(res_in),(sum(res_in(pmap(:,ip)))/size(res_in),ip=1,n_evolve)

     close(1)
  endif

  ! Reset the reset character
  b_flag(:) = ' ' 

  ! Halt TGYRO on GYRO failure
  if (maxval(gyro_exit_status) == 1) then
     call MPI_FINALIZE(ierr)
     stop
  endif

10 format(t2,i2,2x,1pe12.6,2x,4(2x,1pe11.5))

end subroutine tgyro_write_intermediate
