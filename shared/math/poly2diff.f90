!-----------------------------------------------
! poly2diff.f90 
!
! PURPOSE:
!  Compute coefficients for arbitrary-order, 
!  centered 2nd finite-difference derivative.
!  
! NOTES:
!  The coefficients computed are:
!
!                 n
!  f''(0) = ---  Sum c(i) f(i)
!           h^2  i=-n    
!
!  where f(j) means f[x(j)].
!      
!                1           1           -j
! c(i) = Sum_p ----- Sum_q ----- Prod_j -----
!        p/=i   i-p  q/=i   i-q   j/=i   i-j
!                    q/=p         j/=p
!                                 j/=q
!
! REVISIONS
! 18 July 01: jeff.candy@gat.com
!  New.
!--------------------------------------------------------------

subroutine poly2diff(n,c)

  implicit none

  integer, intent(in) :: n

  real, intent(inout), dimension(-n:n) :: c

  integer :: p
  integer :: q
  integer :: i
  integer :: j

  real :: prod


  do i=-n,n
 
     c(i) = 0.0

     do p=-n,n

        if (p /= i) then

           do q=-n,n

              if (q /= i .and. q /= p) then

                 prod = 1.0

                 do j=-n,n

                    if (j /= i .and. j /= q .and. j /= p) then

                    prod = -j*prod/(i-j)

                    endif

                 enddo ! j

                 c(i) = c(i)+prod/((i-p)*(i-q))

              endif

           enddo ! q

        endif

     enddo ! p

  enddo ! i
 
end subroutine poly2diff
