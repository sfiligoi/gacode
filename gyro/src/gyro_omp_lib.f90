module omp_lib

  implicit none

contains

  integer function omp_get_thread_num()

    implicit none

    omp_get_thread_num = 0

  end function omp_get_thread_num

  integer function omp_get_max_threads()

    implicit none

    omp_get_max_threads = 1

  end function omp_get_max_threads

end module omp_lib
