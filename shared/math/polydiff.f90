!-----------------------------------------------
! polydiff.f90 
!
! PURPOSE:
!  Compute coefficients for arbitrary-order, 
!  equally-spaced, uncentered finite-difference 
!  derivatives.
!
! NOTES:
!  The coefficients computed are:
!
!           1   n
!  f'(j) = --- Sum c(p) f(p)
!           h  p=0    
!
!  where h = x(j) - x(j-1) (and f(j) means f[x(j)]).
!       
! 
!
!                   /   Prod   (j-m) , if j /= p
!          1  f(p)  |  m /= p,j
!  c(p) = --- ----  |        
!          h   [p]  |   Sum     Prod   (p-m) , if j=p
!                   \  i /= j  m /= p,i
!
!
! [p] =  Prod  (p-i)
!       i /= p
!
! REVISIONS
! 11 Jan 01: jeff.candy@gat.com
!  Added this revision info; code originally written
!  in fall 2000.
!--------------------------------------------------------------

subroutine polydiff(n,j,c,denom)

  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: j

  real, intent(inout), dimension(0:n) :: c
  real, intent(inout) :: denom

  integer :: p,k,i,m

  real :: temp
  real :: add


  do p=0,n

     temp = 1.0
     do k=0,n
        if (k /= p) temp = temp*(p-k)
     enddo
     denom = temp

     temp = 1.0
     if (j /= p) then

        do m=0,n
           if (m /= p .and. m /= j) temp = temp*(j-m)
        enddo

        c(p) = temp/denom

     else

        add = 0.0
        do i=0,n
           if (i /= j) then
              temp = 1.0
              do m=0,n
                 if (m /= p .and. m /= i) temp = temp*(p-m)
              enddo
              add = add+temp
           endif
        enddo

        c(p) = add/denom

     endif

  enddo

end subroutine polydiff
