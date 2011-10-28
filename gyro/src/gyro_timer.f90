subroutine gyro_timer(indx,tag)

  use mpi
  use gyro_globals

  implicit none

  integer, intent(in) :: indx
  character(len=*), intent(in) :: tag

  if (cpu(indx) < 0.0) then

     cpu(indx)     = 0.0
     cpu_in(indx)  = MPI_Wtime()
     cpu_tag(indx) = tag
     if (indx > cpu_maxindx) cpu_maxindx=indx

  else if (trim(tag) == 'out') then

     cpu(indx) = cpu(indx)+MPI_Wtime()-cpu_in(indx)

  else

     cpu_in(indx) = MPI_Wtime()

  endif

end subroutine gyro_timer
