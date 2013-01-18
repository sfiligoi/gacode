subroutine gyro_init_ompdata()

  use gyro_globals
  use ompdata

  implicit none
  integer, external :: omp_get_max_threads, omp_get_thread_num

  n_omp = omp_get_max_threads()
!$omp parallel 
  i_omp = omp_get_thread_num()
  call looplimits(i_omp,n_omp,1,n_x,ibeg,iend)
!$omp end parallel

end subroutine gyro_init_ompdata

!------------------------------------------------
! block partition do-loop indices i1 to i2 
! for thread : mythd = 0, ..., nthds - 1
! using a method that is as balanced as possible
!------------------------------------------------

subroutine looplimits(mythd,nthds,i1,i2,ibeg,iend)

  implicit none
  integer, intent(in) :: mythd, nthds, i1, i2
  integer, intent(out) :: ibeg, iend
  integer :: nwork, chunk, extra, ntcut  

  nwork = i2 - i1 + 1
  chunk = nwork/nthds
  extra = nwork - nthds*chunk
  ntcut = nthds - extra
  if (mythd < ntcut) then
    ibeg = i1 + mythd*chunk
    iend = ibeg + chunk - 1
  else
    ibeg = i1 + ntcut*chunk + (mythd - ntcut)*(chunk + 1)
    iend = ibeg + chunk
  end if
  if (iend > i2) iend = i2

end subroutine looplimits

