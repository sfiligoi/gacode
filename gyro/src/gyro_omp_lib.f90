module omp_lib

contains

  integer function omp_get_thread_num()

    omp_get_thread_num = 0

  end function omp_get_thread_num

  integer function omp_get_max_threads()

    omp_get_max_threads = 1

  end function omp_get_max_threads

end module omp_lib
