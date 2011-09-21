module gkcoll_legendre
  implicit none

  public :: legendre, legendre_deriv

contains

  subroutine legendre(n,arg,val)
    integer, intent (in) :: n
    real, intent (in) :: arg
    real, intent(out) :: val
    real :: pmm, pmmp1, pnn
    integer :: k
    
    pmm=1.0
    if(n==0) then
       val = pmm
    else
       pmmp1 = arg*pmm;
       if(n==1) then
          val = pmmp1
       else
          do k=2, n
             pnn = (arg*(2*k-1)*pmmp1 - (k-1)*pmm)/(1.0*k)
             pmm=pmmp1
             pmmp1=pnn
          enddo
          val = pnn
       end if
    end if
    
  end subroutine legendre

  subroutine legendre_deriv(n,arg,val)
    integer, intent (in) :: n
    real, intent (in) :: arg
    real, intent(out) :: val
    real :: pmm1, pmm2, pnn
    integer :: k
    
    pmm2=0.0
    pmm1=1.0
    if(n==0) then
       val = pmm2
    else if(n==1) then
       val = pmm1
    else
       do k=2, n
          call legendre(k-1,arg,pnn)
          val = pmm2 + (2.0*k-1.0) * pnn
          pmm2 = pmm1
          pmm1 = val
       enddo
    endif
 
  end subroutine legendre_deriv

end module gkcoll_legendre
