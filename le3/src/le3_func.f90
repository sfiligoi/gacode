subroutine le3_func(xsize,x,fvec,iflag)

  use le3_globals

  implicit none

  integer :: id, k, i, j
  real :: ang
  real :: iota
  integer, intent(in) :: xsize
  integer, intent(inout) :: iflag
  real, dimension(xsize), intent(inout) :: fvec
  real, dimension(xsize), intent(in) :: x
  real :: tb2
  integer, parameter :: fix=1

  iota = 1.0/q

  ! tb is really the periodic function thetabar-theta

  k=1
  do j=1,np
     do i=2,nt
        tb(i,j) = x(k)
        k = k+1
     enddo
  enddo

  !-------------------------------------------------
  ! Compute d(tb)/dt and d(tp)/dp over [0,2pi) with 

  ! d(tb)/dt
  dtbdt(:,:) = 1.0
  do i=1,nt
     do id=-2,2
        k = tcyc(i+id)
        dtbdt(i,:) = dtbdt(i,:) + tb(k,:) * cderiv(id)/dt
     enddo
  enddo

  ! d(tb)/dp
  dtbdp(:,:) = 0.0
  do j=1,np
     do id=-2,2
        k = pcyc(j+id)
        dtbdp(:,j) = dtbdp(:,j) + tb(:,k) * cderiv(id)/dp
     enddo
  enddo

  !-------------------------------------------------

  do i=1,nt
     do j=1,np

        ang = m*(tb(i,j)+t(i))+n*p(j)

        ! R,Z
        r(i,j) = rmaj + rmin*cos(tb(i,j)+t(i)) + hmin*cos(ang)
        z(i,j) =         rmin*sin(tb(i,j)+t(i)) + hmin*sin(ang)

        ! dR/d(tb)
        drdtb(i,j) = -rmin*sin(tb(i,j)+t(i)) - m*hmin*sin(ang) 

        ! dR/d(pb)
        drdpb(i,j) =                   -n*hmin*sin(ang)

        ! dZ/d(tb)
        dzdtb(i,j) =  rmin*cos(tb(i,j)+t(i)) + m*hmin*cos(ang) 

        ! dZ/d(pb)
        dzdpb(i,j) =                  n*hmin*cos(ang)

        ! dR/dr
        drdr(i,j) = cos(tb(i,j)+t(i))

        ! dZ/dr
        dzdr(i,j) = sin(tb(i,j)+t(i))

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

  fp(:,:) = (br(:,:)*rp(:,:)+bz(:,:)*zp(:,:))/jac(:,:)/dtbdt(:,:)
  ft(:,:) = (br(:,:)*rt(:,:)+bz(:,:)*zt(:,:))/jac(:,:)/dtbdt(:,:)

  ! d(fp)/dt
  fpt(:,:) = 0.0
  do i=1,nt
     do id=-2,2
        k = tcyc(i+id)
        fpt(i,:) = fpt(i,:) + fp(k,:) * cderiv(id)/dt
     enddo
  enddo

  ! d(ft)/dp
  ftp(:,:) = 0.0
  do j=1,np
     do id=-2,2
        k = pcyc(j+id)
        ftp(:,j) = ftp(:,j) + ft(:,k) * cderiv(id)/dp
     enddo
  enddo

  fp(:,:) = bp(:,:)*r(:,:)/jac(:,:)
  do j=1,np
     do i=1,nt

        ! d^2(theta_bar)/dt^2
        tb2 = (tb(tcyc(i+1),j)-2.0*tb(i,j)+tb(tcyc(i-1),j))/dt**2

        ! Add "axisymmetric term": d/dtheta ( fp / dtbdt ) 
        fpt(i,j) = fpt(i,j) &
             +(fp(tcyc(i+1),j)-fp(tcyc(i-1),j))/(2*dt)/dtbdt(i,j) &
             -tb2*fp(i,j)/(dtbdt(i,j)**2)
     enddo
  enddo

  k=1
  do j=1,np
     do i=2,nt
        fvec(k) = fpt(i,j)-ftp(i,j)
        k = k+1
     enddo
  enddo

end subroutine le3_func
