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

  call le3_map(x,as,bs,cs,ds,nps,nts,'setc')

  !--------------------------------------------------------
  ! Use spectral form to compute tb, d(tb)/dt and d(tp)/dp 
  ! on finite grids over [0,2pi) 
  !
  do i=1,nt
     tb(i,:) = t(i)
  enddo
  dtbdt(:,:) = 1.0
  dtbdp(:,:) = 0.0

  do ips=0,nps
     do its=0,nts
        do j=1,np
           do i=1,nt

              tb(i,j) = tb(i,j) + &
                   sinm(i,its)*(bs(its,ips)*sinn(j,ips)+as(its,ips)*cosn(j,ips)) &
                   +cosm(i,its)*(ds(its,ips)*sinn(j,ips)+cs(its,ips)*cosn(j,ips))

              dtbdt(i,j) = dtbdt(i,j) + &
                   its*cosm(i,its)*(bs(its,ips)*sinn(j,ips)+as(its,ips)*cosn(j,ips)) &
                   -its*sinm(i,its)*(ds(its,ips)*sinn(j,ips)+cs(its,ips)*cosn(j,ips))

              dtbdp(i,j) = dtbdp(i,j) + &
                   ips*sinm(i,its)*(bs(its,ips)*cosn(j,ips)-as(its,ips)*sinn(j,ips)) &
                   +ips*cosm(i,its)*(ds(its,ips)*cosn(j,ips)-cs(its,ips)*sinn(j,ips))
           enddo
        enddo

     enddo
  enddo
  !--------------------------------------------------------

  do i=1,nt
     do j=1,np

        ang = m*tb(i,j)+n*p(j)

        ! R,Z
        r(i,j) = rmaj + rmin * cos(tb(i,j) + asin(delta) * sin(tb(i,j))) &
             + hmin*cos(ang)

        z(i,j) = zmag + kappa * rmin * sin(tb(i,j) + zeta * sin(2*tb(i,j))) &
             + hmin*sin(ang)

        ! dR/d(tb)
        drdtb(i,j) = -rmin * sin(tb(i,j)+ asin(delta) * sin(tb(i,j))) &
             * (1.0 + asin(delta) * cos(tb(i,j))) &
             - m*hmin*sin(ang) 

        ! dR/d(pb)
        drdpb(i,j) =  -n*hmin*sin(ang)                 
             

        ! dZ/d(tb)
        dzdtb(i,j) =  kappa * rmin * cos(tb(i,j) + zeta * sin(2*tb(i,j))) &
             * (1.0 + 2*zeta * cos(2*tb(i,j))) &
             + m*hmin*cos(ang) 

        ! dZ/d(pb)
        dzdpb(i,j) =  n*hmin*cos(ang)

        ! dR/dr
        drdr(i,j) = shift + cos(tb(i,j) + asin(delta) * sin(tb(i,j))) &
             - sin(tb(i,j) + asin(delta) * sin(tb(i,j))) * sin(tb(i,j)) &
             * 1.0/sqrt(1.0-delta**2) * s_delta

        ! dZ/dr
        dzdr(i,j) = dzmag + (1.0 + s_kappa) * kappa &
             * sin(tb(i,j) + zeta * sin(2*tb(i,j))) &
             + kappa * cos(tb(i,j) + zeta * sin(2*tb(i,j))) &
             * sin(2*tb(i,j)) * s_zeta

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

  fp(:,:) = (bp(:,:)*r(:,:)+br(:,:)*rp(:,:)+bz(:,:)*zp(:,:))/jac(:,:)/dtbdt(:,:)
  ft(:,:) = (br(:,:)*rt(:,:)+bz(:,:)*zt(:,:))/jac(:,:)/dtbdt(:,:)


  ix = 0
  fvec(:) = 0.0
  do ips=0,nps
     do its=0,nts

        !-------------------------------------------------------------
        ! Projections
        !
        ! NOTE: each projection is a nonlinear equation to be solved 
        !       as a nonlinear system by MINPACK.
        !-------------------------------------------------------------

        ! A: sin m cos n
        if (its > 0) then 
           ix = ix+1           
           do j=1,np
              do i=1,nt
                 fvec(ix) = fvec(ix)  &
                      -its*cosm(i,its)*cosn(j,ips)*fp(i,j) &
                      -ips*sinm(i,its)*sinn(j,ips)*ft(i,j) 
              enddo
           enddo
           if (ix == xsize) exit
        endif

        ! B: sin m sin n
        if (ips > 0 .and. its > 0) then
           ix = ix+1
           do j=1,np
              do i=1,nt
                 fvec(ix) = fvec(ix) &
                      -its*cosm(i,its)*sinn(j,ips)*fp(i,j) &
                      +ips*sinm(i,its)*cosn(j,ips)*ft(i,j) 
              enddo
           enddo
           if (ix == xsize) exit
        endif

        ! C: cos m cos n
        if (ips + its > 0) then
           ix = ix+1
           do j=1,np
              do i=1,nt
                 fvec(ix) = fvec(ix) &
                      +its*sinm(i,its)*cosn(j,ips)*fp(i,j) &
                      -ips*cosm(i,its)*sinn(j,ips)*ft(i,j) 
              enddo
           enddo
           if (ix == xsize) exit
        endif

        ! D: cos m sin n
        if (ips > 0) then
           ix = ix+1
           do j=1,np
              do i=1,nt
                 fvec(ix) = fvec(ix) &
                      +its*sinm(i,its)*sinn(j,ips)*fp(i,j) &
                      +ips*cosm(i,its)*cosn(j,ips)*ft(i,j) 
              enddo
           enddo
           if (ix == xsize) exit
        endif

     enddo
  enddo

end subroutine le3_func2
