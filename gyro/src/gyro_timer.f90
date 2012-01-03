subroutine gyro_timer_init(tag)

  use gyro_globals

  implicit none
  character(len=*), intent(in) :: tag

  cpu_maxindx = cpu_maxindx+1
  cpu_tag(cpu_maxindx) = trim(tag)

end subroutine gyro_timer_init

subroutine gyro_timer_in(tag)

  use mpi
  use gyro_globals

  implicit none
  character(len=*), intent(in) :: tag
  integer :: indx

  do indx=1,cpu_maxindx
     if (trim(tag) == trim(cpu_tag(indx))) then
        cpu_in(indx) = MPI_Wtime()
     endif
  enddo

  cpu_in(indx) = MPI_Wtime()

end subroutine gyro_timer_in

subroutine gyro_timer_out(tag)

  use mpi
  use gyro_globals

  implicit none
  character(len=*), intent(in) :: tag
  integer :: indx

  do indx=1,cpu_maxindx
     if (trim(tag) == trim(cpu_tag(indx))) then
        cpu(indx) = cpu(indx)+MPI_Wtime()-cpu_in(indx)
     endif
  enddo

end subroutine gyro_timer_out

