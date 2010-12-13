      subroutine proc_time(seconds)

      implicit none
      include 'mpif.h'

      double precision seconds
 
      seconds = MPI_WTIME()
 
      end subroutine proc_time
