subroutine le3_alloc(flag)

  use le3_globals

  implicit none
  integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
  integer :: i,j,its,ips

  if (flag == 1) then

     if (initialized) return

     allocate(t(nt))
     allocate(p(np))
     allocate(tb(nt,np))
     allocate(dtbdt(nt,np))
     allocate(dtbdp(nt,np))
     allocate(r(nt,np))
     allocate(z(nt,np))
     allocate(drdtb(nt,np))
     allocate(drdpb(nt,np))
     allocate(dzdtb(nt,np))
     allocate(dzdpb(nt,np))
     allocate(jac(nt,np))
     allocate(rp(nt,np))
     allocate(rt(nt,np))
     allocate(zp(nt,np))
     allocate(zt(nt,np))
     allocate(fp(nt,np))
     allocate(ft(nt,np))

     do i=1,nt
        t(i) = 2*(i-1)*pi/(nt)
     enddo
     dt = t(2)-t(1)
     do j=1,np
        p(j) = 2*(j-1)*pi/(np)
     enddo
     dp = p(2)-p(1)

     ! Spectral
     msize = 4*nts*nps+2*(nts+nps)

     ! Fourier coefficients
     allocate(as(0:nts,0:nps))
     allocate(bs(0:nts,0:nps))
     allocate(cs(0:nts,0:nps))
     allocate(ds(0:nts,0:nps))

     ! Storage to eliminate direct sin/cos evaluation
     allocate(sinm(nt,0:nts))
     allocate(cosm(nt,0:nts))
     allocate(sinn(np,0:nps))
     allocate(cosn(np,0:nps))

     do its=0,nts
        do i=1,nt
           sinm(i,its) = sin(its*t(i))
           cosm(i,its) = cos(its*t(i))
        enddo
     enddo
     do ips=0,nps
        do j=1,np
           sinn(j,ips) = sin(ips*p(j))
           cosn(j,ips) = cos(ips*p(j))
        enddo
     enddo

     allocate(yfunc(msize))
     allocate(xfunc(msize))
     nwork = (msize*(3*msize+13))/2*10
     allocate(work(nwork))

     initialized = .true.

  else

     if (.NOT. initialized) return

     deallocate(t)
     deallocate(p)
     deallocate(tb)
     deallocate(dtbdt)
     deallocate(dtbdp)
     deallocate(r)
     deallocate(z)
     deallocate(drdtb)
     deallocate(drdpb)
     deallocate(dzdtb)
     deallocate(dzdpb)
     deallocate(jac)
     deallocate(rp)
     deallocate(rt)
     deallocate(zp)
     deallocate(zt)
     deallocate(fp)
     deallocate(ft)
     deallocate(xfunc)
     deallocate(yfunc)
     deallocate(work)

     initialized = .false.

  endif

end subroutine le3_alloc
