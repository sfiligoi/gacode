!-----------------------------------------------------------------
! cgyro_math.F90
!
! PURPOSE:
!  Helper functions for scalable math operations.
!-----------------------------------------------------------------

module cgyro_math

    implicit none 

contains

  !=========================================================
  ! Multiple-add of array, updating left in place
  !=========================================================

subroutine cgyro_cmpl_inplace_fma1(sz, left, c1, r1)
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: sz
    complex, intent(inout), dimension(*) :: left
    real, intent(in) :: c1
    complex, intent(in), dimension(*) :: r1
    !
    integer :: i
    !-------------------------------------------------------
#ifdef _OPENACC
!$acc parallel loop independent present(left,r1)
    do i=1,sz
       left(i) = left(i) + c1 * r1(i)
    enddo
#else
!$omp parallel do 
    do i=1,sz
       left(i) = left(i) + c1 * r1(i)
    enddo
#endif
end subroutine cgyro_cmpl_inplace_fma1

subroutine cgyro_cmpl_inplace_fma2(sz, left, c1, r1, c2, r2)
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: sz
    complex, intent(inout), dimension(*) :: left
    real, intent(in) :: c1
    complex, intent(in), dimension(*) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(*) :: r2
    !
    integer :: i
    !-------------------------------------------------------
#ifdef _OPENACC
!$acc parallel loop independent present(left,r1,r2)
    do i=1,sz
       left(i) = left(i) + c1 * r1(i) + c2 * r2(i)
    enddo
#else
!$omp parallel do 
    do i=1,sz
       left(i) = left(i) + c1 * r1(i) + c2 * r2(i)
    enddo
#endif
end subroutine cgyro_cmpl_inplace_fma2

subroutine cgyro_cmpl_inplace_fma3(sz, left, c1, r1, c2, r2, c3, r3)
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: sz
    complex, intent(inout), dimension(*) :: left
    real, intent(in) :: c1
    complex, intent(in), dimension(*) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(*) :: r2
    real, intent(in) :: c3
    complex, intent(in), dimension(*) :: r3
    !
    integer :: i
    !-------------------------------------------------------
#ifdef _OPENACC
!$acc parallel loop independent present(left,r1,r2,r3)
    do i=1,sz
       left(i) = left(i) + c1 * r1(i) + c2 * r2(i) + c3 * r3(i)
    enddo
#else
!$omp parallel do 
    do i=1,sz
       left(i) = left(i) + c1 * r1(i) + c2 * r2(i) + c3 * r3(i)
    enddo
#endif
end subroutine cgyro_cmpl_inplace_fma3


subroutine cgyro_cmpl_inplace_fma4(sz, left, c1, r1, c2, r2, c3, r3, c4, r4)
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: sz
    complex, intent(inout), dimension(*) :: left
    real, intent(in) :: c1
    complex, intent(in), dimension(*) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(*) :: r2
    real, intent(in) :: c3
    complex, intent(in), dimension(*) :: r3
    real, intent(in) :: c4
    complex, intent(in), dimension(*) :: r4
    !
    integer :: i
    !-------------------------------------------------------
#ifdef _OPENACC
!$acc parallel loop independent present(left,r1,r2,r3,r4)
    do i=1,sz
       left(i) = left(i) + c1 * r1(i) + c2 * r2(i) + c3 * r3(i)+ c4 * r4(i)
    enddo
#else
!$omp parallel do 
    do i=1,sz
       left(i) = left(i) + c1 * r1(i) + c2 * r2(i) + c3 * r3(i)+ c4 * r4(i)
    enddo
#endif
end subroutine cgyro_cmpl_inplace_fma4


  !=========================================================
  ! Copy one or more arrays
  !=========================================================

subroutine cgyro_cmpl_copy(sz, left, r1)
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: sz
    complex, intent(out), dimension(*) :: left
    complex, intent(in), dimension(*) :: r1
    !
    integer :: i
    !-------------------------------------------------------
#ifdef _OPENACC
!$acc parallel loop independent present(left,r1)
    do i=1,sz
       left(i) = r1(i)
    enddo
#else
!$omp parallel do 
    do i=1,sz
       left(i) = r1(i)
    enddo
#endif
end subroutine cgyro_cmpl_copy

subroutine cgyro_cmpl_copy2(sz, left1, left2, r1)
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: sz
    complex, intent(out), dimension(*) :: left1
    complex, intent(out), dimension(*) :: left2
    complex, intent(in), dimension(*) :: r1
    !
    integer :: i
    complex :: tmp
    !-------------------------------------------------------
#ifdef _OPENACC
!$acc parallel loop independent present(left1,left2,r1) private(tmp)
    do i=1,sz
       tmp = r1(i)
       left1(i) = tmp
       left2(i) = tmp
    enddo
#else
!$omp parallel do private(tmp)
    do i=1,sz
       tmp = r1(i)
       left1(i) = tmp
       left2(i) = tmp
    enddo
#endif
end subroutine cgyro_cmpl_copy2

  !=========================================================
  ! Multiple-add of array without using the value of left
  !=========================================================

subroutine cgyro_cmpl_fma2(sz, left, r1, c2, r2, abssum)
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: sz
    complex, intent(out), dimension(*) :: left
    complex, intent(in), dimension(*) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(*) :: r2
    real, intent(inout), optional :: abssum
    !
    integer :: i
    complex :: tmp
    real :: s
    !-------------------------------------------------------
    if (present(abssum)) then
      s = 0.0
#ifdef _OPENACC
!$acc parallel loop independent present(left,r1,r2) private(tmp) reduction(+:s)
#else
!$omp parallel do private(tmp) reduction(+:s)
#endif
      do i=1,sz
        tmp = r1(i) + c2 * r2(i)
        left(i) = tmp
        s = s + abs(tmp)
      enddo
      abssum = s
    else
#ifdef _OPENACC
!$acc parallel loop independent present(left,r1,r2)
#else
!$omp parallel do
#endif
      do i=1,sz
        left(i) = r1(i) + c2 * r2(i)
      enddo
    endif
end subroutine cgyro_cmpl_fma2

subroutine cgyro_cmpl_fma3(sz, left, r1, c2, r2, c3, r3, abssum)
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: sz
    complex, intent(out), dimension(*) :: left
    complex, intent(in), dimension(*) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(*) :: r2
    real, intent(in) :: c3
    complex, intent(in), dimension(*) :: r3
    real, intent(inout), optional :: abssum
    !
    integer :: i
    complex :: tmp
    real :: s
    !-------------------------------------------------------
    if (present(abssum)) then
      s = 0.0
#ifdef _OPENACC
!$acc parallel loop independent present(left,r1,r2,r3) private(tmp) reduction(+:s)
#else
!$omp parallel do private(tmp) reduction(+:s)
#endif
      do i=1,sz
        tmp = r1(i) + c2 * r2(i) + c3 * r3(i)
        left(i) = tmp
        s = s + abs(tmp)
      enddo
      abssum = s
    else
#ifdef _OPENACC
!$acc parallel loop independent present(left,r1,r2,r3)
#else
!$omp parallel do
#endif
      do i=1,sz
        left(i) = r1(i) + c2 * r2(i) + c3 * r3(i)
      enddo
    endif
end subroutine cgyro_cmpl_fma3

subroutine cgyro_cmpl_fma4(sz, left, r1, c2, r2, c3, r3, c4, r4, abssum)
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: sz
    complex, intent(out), dimension(*) :: left
    complex, intent(in), dimension(*) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(*) :: r2
    real, intent(in) :: c3
    complex, intent(in), dimension(*) :: r3
    real, intent(in) :: c4
    complex, intent(in), dimension(*) :: r4
    real, intent(inout), optional :: abssum
    !
    integer :: i
    complex :: tmp
    real :: s
    !-------------------------------------------------------
    if (present(abssum)) then
      s = 0.0
#ifdef _OPENACC
!$acc parallel loop independent present(left,r1,r2,r3,r4) private(tmp) reduction(+:s)
#else
!$omp parallel do private(tmp) reduction(+:s)
#endif
      do i=1,sz
        tmp = r1(i) + c2 * r2(i) + c3 * r3(i)+ c4 * r4(i)
        left(i) = tmp
        s = s + abs(tmp)
      enddo
      abssum = s
    else
#ifdef _OPENACC
!$acc parallel loop independent present(left,r1,r2,r3,r4)
#else
!$omp parallel do
#endif
      do i=1,sz
        left(i) = r1(i) + c2 * r2(i) + c3 * r3(i)+ c4 * r4(i)
      enddo
    endif
end subroutine cgyro_cmpl_fma4

subroutine cgyro_cmpl_fma5(sz, left, r1, c2, r2, c3, r3, c4, r4, c5, r5, abssum)
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: sz
    complex, intent(out), dimension(*) :: left
    complex, intent(in), dimension(*) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(*) :: r2
    real, intent(in) :: c3
    complex, intent(in), dimension(*) :: r3
    real, intent(in) :: c4
    complex, intent(in), dimension(*) :: r4
    real, intent(in) :: c5
    complex, intent(in), dimension(*) :: r5
    real, intent(inout), optional :: abssum
    !
    integer :: i
    complex :: tmp
    real :: s
    !-------------------------------------------------------
    if (present(abssum)) then
      s = 0.0
#ifdef _OPENACC
!$acc parallel loop independent present(left,r1,r2,r3,r4,r5) private(tmp) reduction(+:s)
#else
!$omp parallel do private(tmp) reduction(+:s)
#endif
      do i=1,sz
        tmp = r1(i) + c2 * r2(i) + c3 * r3(i)+ c4 * r4(i) + c5 * r5(i)
        left(i) = tmp
        s = s + abs(tmp)
      enddo
      abssum = s
    else
#ifdef _OPENACC
!$acc parallel loop independent present(left,r1,r2,r3,r4,r5)
#else
!$omp parallel do
#endif
      do i=1,sz
        left(i) = r1(i) + c2 * r2(i) + c3 * r3(i)+ c4 * r4(i) + c5 * r5(i)
      enddo
    endif
end subroutine cgyro_cmpl_fma5

subroutine cgyro_cmpl_fma6(sz, left, r1, c2, r2, c3, r3, c4, r4, c5, r5, c6, r6, abssum)
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: sz
    complex, intent(out), dimension(*) :: left
    complex, intent(in), dimension(*) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(*) :: r2
    real, intent(in) :: c3
    complex, intent(in), dimension(*) :: r3
    real, intent(in) :: c4
    complex, intent(in), dimension(*) :: r4
    real, intent(in) :: c5
    complex, intent(in), dimension(*) :: r5
    real, intent(in) :: c6
    complex, intent(in), dimension(*) :: r6
    real, intent(inout), optional :: abssum
    !
    integer :: i
    complex :: tmp
    real :: s
    !-------------------------------------------------------
    if (present(abssum)) then
      s = 0.0
#ifdef _OPENACC
!$acc parallel loop independent present(left,r1,r2,r3,r4,r5,r6) private(tmp) reduction(+:s)
#else
!$omp parallel do private(tmp) reduction(+:s)
#endif
      do i=1,sz
        tmp = r1(i) + c2 * r2(i) + c3 * r3(i)+ c4 * r4(i) + c5 * r5(i) + c6 * r6(i)
        left(i) = tmp
        s = s + abs(tmp)
      enddo
      abssum = s
    else
#ifdef _OPENACC
!$acc parallel loop independent present(left,r1,r2,r3,r4,r5,r6)
#else
!$omp parallel do
#endif
      do i=1,sz
        left(i) = r1(i) + c2 * r2(i) + c3 * r3(i)+ c4 * r4(i) + c5 * r5(i) + c6 * r6(i)
      enddo
    endif
end subroutine cgyro_cmpl_fma6

end module cgyro_math

