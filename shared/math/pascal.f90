!--------------------------------------------------------------
! pascal.f90 
!
! PURPOSE:
!  Compute coefficients in Pascal's triangle, with added 
!  alternating sign.
!
! NOTES:
!  Equivalently, these are the coefficients of the 
!  centered n-th derivative operator.  Below, we have 
!  supressed the alternating sign.
!
!  n    s(1)  s(2)  s(3) ...
! ---   ---------------------------------------------------------
!  3    1.0   2.0   1.0
!  4    1.0   3.0   3.0   1.0
!  5    1.0   4.0   6.0   4.0   1.0
!  6    1.0   5.0  10.0  10.0   5.0   1.0
!  7    1.0   6.0  15.0  20.0  15.0   6.0   1.0
!  8    1.0   7.0  21.0  35.0  35.0  21.0   7.0   1.0
!  9    1.0   8.0  28.0  56.0  70.0  56.0  28.0   8.0   1.0
! 
! REVISIONS
! 06 Feb 01: jeff.candy@gat.com
!  Created. 
!--------------------------------------------------------------

subroutine pascal(n,s)

  implicit none

  integer, intent(in) :: n
  integer :: i
  integer :: m
  real, dimension(n) :: s
  real, dimension(n) :: s0

  s = 0 

  if (n < 3) then
     print *,'PASCAL: error in n'
     stop
  endif

  if (n == 3) then

     s(1) = 1.0  
     s(2) = -2.0
     s(3) = 1.0

  else

     s0(1) = 1.0  
     s0(2) = 2.0
     s0(3) = 1.0

     do m=4,n

        s(1) = 1.0
        do i=2,m-1
           s(i) = s0(i)+s0(i-1)
        enddo
        s(m) = 1.0

        s0 = s

     enddo

     ! Ensure middle element is negative for n odd.
 
     do m=1,n
        s(m) = (-1)**((n+1)/2-m-1)*s(m)
     enddo

  endif

end subroutine pascal
