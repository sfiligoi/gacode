module omp_lib

contains

  integer function omp_get_thread_num()

    omp_get_thread_num = -1

  end function omp_get_thread_num

end module omp_lib
