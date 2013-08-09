subroutine le3_func2(xsize,x,fvec,iflag)

  use le3_globals

  implicit none

  integer :: ix,i,j,its,ips
  real :: ang
  real :: iota
  integer, intent(in) :: xsize
  integer, intent(inout) :: iflag
  real, dimension(xsize), intent(inout) :: fvec
  real, dimension(xsize), intent(in) :: x

  iota = 1.0/q

  ix=0
  do ips=0,nps
     do its=1,nts 
        ix = ix+1           
        as(its,ips) = x(ix)
        if (ips > 0) then
           ix = ix+1
           bs(its,ips) = x(ix)
        endif
        ix = ix+1           
        cs(its,ips) = x(ix)
        if (ips > 0) then
           ix = ix+1
           ds(its,ips) = x(ix)
        endif
     enddo
  enddo
  bs(:,0) = 0.0
  ds(:,0) = 0.0

  !--------------------------------------------------------
  ! Use spectral form to compute tb, d(tb)/dt and d(tp)/dp 
  ! on finite grids over [0,2pi) 
  !
  do i=1,nt
     do j=1,np

        tb(i,j)    = t(i)
        dtbdt(i,j) = 1.0
        dtbdp(i,j) = 0.0

        do its=1,nts
           do ips=0,nps
              tb(i,j) = tb(i,j) + &
                   sin(its*t(i))*(bs(its,ips)*sin(ips*p(j))+as(its,ips)*cos(ips*p(j))) &
                   +(cos(its*t(i))-1)*(ds(its,ips)*sin(ips*p(j))+cs(its,ips)*cos(ips*p(j)))
              dtbdt(i,j) = dtbdt(i,j) + &
                   its*cos(its*t(i))*(bs(its,ips)*sin(ips*p(j))+as(its,ips)*cos(ips*p(j))) &
                   -its*sin(its*t(i))*(ds(its,ips)*sin(ips*p(j))+cs(its,ips)*cos(ips*p(j)))
              dtbdp(i,j) = dtbdp(i,j) + &
                   ips*sin(its*t(i))*(bs(its,ips)*cos(ips*p(j))-as(its,ips)*sin(ips*p(j))) &
                   +ips*(cos(its*t(i))-1)*(ds(its,ips)*cos(ips*p(j))-cs(its,ips)*sin(ips*p(j)))
           enddo
        enddo

     enddo
  enddo
  !--------------------------------------------------------

  do i=1,nt
     do j=1,np

        ang = m*tb(i,j)+n*p(j)

        ! R,Z
        r(i,j) = rmaj + rmin*cos(tb(i,j)) + hmin*cos(ang)
        z(i,j) =         rmin*sin(tb(i,j)) + hmin*sin(ang)

        ! dR/d(tb)
        drdtb(i,j) = -rmin*sin(tb(i,j)) - m*hmin*sin(ang) 

        ! dR/d(pb)
        drdpb(i,j) =                   -n*hmin*sin(ang)

        ! dZ/d(tb)
        dzdtb(i,j) =  rmin*cos(tb(i,j)) + m*hmin*cos(ang) 

        ! dZ/d(pb)
        dzdpb(i,j) =                  n*hmin*cos(ang)

        ! dR/dr
        drdr(i,j) = cos(tb(i,j))

        ! dZ/dr
        dzdr(i,j) = sin(tb(i,j))

        ! J [dtb/dr terms vanish]
        jac(i,j) = r(i,j)*(drdr(i,j)*dzdtb(i,j)-drdtb(i,j)*dzdr(i,j))

     enddo
  enddo

  bp(:,:) = r(:,:)
  br(:,:) = drdpb(:,:)+drdtb(:,:)*(dtbdp(:,:)+iota*dtbdt(:,:))
  bz(:,:) = dzdpb(:,:)+dzdtb(:,:)*(dtbdp(:,:)+iota*dtbdt(:,:))

  rp(:,:) = drdpb(:,:)+dtbdp(:,:)*drdtb(:,:)
  zp(:,:) = dzdpb(:,:)+dtbdp(:,:)*dzdtb(:,:)
  rt(:,:) = drdtb(:,:)*dtbdt(:,:)
  zt(:,:) = dzdtb(:,:)*dtbdt(:,:)

  fp(:,:) = (r(:,:)**2+br(:,:)*rp(:,:)+bz(:,:)*zp(:,:))/jac(:,:)/dtbdt(:,:)
  ft(:,:) = (br(:,:)*rt(:,:)+bz(:,:)*zt(:,:))/jac(:,:)/dtbdt(:,:)


  ix = 0
  fvec(:) = 0.0
  do ips=0,nps
     do its=1,nts

        ! Projections

        ! A 
        ix = ix+1           
        do j=1,np
           do i=1,nt
              fvec(ix) = fvec(ix)  &
                   -its*cos(its*t(i))*cos(ips*p(j))*fp(i,j) &
                   -ips*sin(its*t(i))*sin(ips*p(j))*ft(i,j) 
           enddo
        enddo
        if (ix == xsize) exit

        ! B 
        if (ips > 0) then
           ix = ix+1
           do j=1,np
              do i=1,nt
                 fvec(ix) = fvec(ix) &
                      -its*cos(its*t(i))*sin(ips*p(j))*fp(i,j) &
                      +ips*sin(its*t(i))*cos(ips*p(j))*ft(i,j) 
              enddo
           enddo
        endif
        if (ix == xsize) exit

        ! C
        ix = ix+1
        do j=1,np
           do i=1,nt
              fvec(ix) = fvec(ix) &
                   +its*sin(its*t(i))*cos(ips*p(j))*fp(i,j) &
                   -ips*cos(its*t(i))*sin(ips*p(j))*ft(i,j) 
           enddo
        enddo
        if (ix == xsize) exit

        ! D
        if (ips > 0) then
           ix = ix+1
           do j=1,np
              do i=1,nt
                 fvec(ix) = fvec(ix) &
                      +its*sin(its*t(i))*sin(ips*p(j))*fp(i,j) &
                      +ips*cos(its*t(i))*cos(ips*p(j))*ft(i,j) 
              enddo
           enddo
        endif
        if (ix == xsize) exit

     enddo
  enddo

end subroutine le3_func2
