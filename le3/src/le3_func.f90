subroutine le3_func(xsize,x,fvec,iflag)

  use le3_globals

  implicit none

  integer :: ix,i,j,its,ips
  integer, intent(in) :: xsize
  integer, intent(inout) :: iflag
  real, dimension(xsize), intent(inout) :: fvec
  real, dimension(xsize), intent(in) :: x
  real, dimension(:,:), allocatable :: bp,br,bz

  call le3_map(x,as,bs,cs,ds,nps,nts,'setc')

  allocate(bp(nt,np))
  allocate(br(nt,np))
  allocate(bz(nt,np))

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

        call le3_rz(tb(i,j),&
             p(j),&
             r(i,j),&
             z(i,j),&
             drdtb(i,j),&
             drdpb(i,j),&
             dzdtb(i,j),&
             dzdpb(i,j),&
             jac(i,j))

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

  deallocate(bp)
  deallocate(br)
  deallocate(bz)

end subroutine le3_func
