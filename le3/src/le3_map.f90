subroutine le3_map(x,a,b,c,d,nps,nts,direction)

  implicit none
  character(len=4), intent(in) :: direction 
  integer, intent(in) :: nps,nts
  real, intent(inout), dimension(0:nts,0:nps) :: a,b,c,d
  real, intent(inout), dimension(4*nts*nps+2*(nts+nps)) :: x
  integer :: ips,its,ix

  if (direction == 'setx') then

     ix=0
     do ips=0,nps
        do its=0,nts
           if (its > 0) then 
              ix = ix+1           
              x(ix) = a(its,ips)
           endif
           if (ips > 0 .and. its > 0) then
              ix = ix+1
              x(ix) = b(its,ips)
           endif
           if (ips + its > 0) then
              ix = ix+1           
              x(ix) = c(its,ips)
           endif
           if (ips > 0) then
              ix = ix+1
              x(ix) = d(its,ips)
           endif
        enddo
     enddo

  else if (direction == 'setc') then

     a = 0.0
     b = 0.0  
     c = 0.0
     d = 0.0

     ix=0
     do ips=0,nps
        do its=0,nts
           if (its > 0) then 
              ix = ix+1           
              a(its,ips) = x(ix)
           endif
           if (ips > 0 .and. its > 0) then
              ix = ix+1
              b(its,ips) = x(ix)
           endif
           if (ips + its > 0) then
              ix = ix+1           
              c(its,ips) = x(ix)
           endif
           if (ips > 0) then
              ix = ix+1
              d(its,ips) = x(ix)
           endif
        enddo
     enddo

  endif

end subroutine le3_map
