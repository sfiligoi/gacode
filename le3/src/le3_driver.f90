module le3_driver

  implicit none

  public :: le3_alloc, le3_solver

  real, dimension(:), allocatable, private :: mhdfunc
  real, dimension(:), allocatable, private :: xfunc
  integer, private :: info
  integer, private :: nwork
  integer, private :: msize
  real, dimension(:), allocatable, private :: work
  logical, private :: initialized = .false.

contains

  subroutine le3_alloc(flag)
    use le3_globals
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    integer :: i,j

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
       allocate(drdr(nt,np))
       allocate(dzdtb(nt,np))
       allocate(dzdpb(nt,np))
       allocate(dzdr(nt,np))
       allocate(jac(nt,np))
       allocate(bp(nt,np))
       allocate(br(nt,np))
       allocate(bz(nt,np))
       allocate(rp(nt,np))
       allocate(rt(nt,np))
       allocate(zp(nt,np))
       allocate(zt(nt,np))
       allocate(fp(nt,np))
       allocate(ft(nt,np))
       allocate(fpt(nt,np))
       allocate(ftp(nt,np))
       allocate(tcyc(1-nt:2*nt))
       allocate(pcyc(1-np:2*np))

       do i=1,nt
          t(i) = 2*(i-1)*pi/(nt)
       enddo
       dt = t(2)-t(1)
       do j=1,np
          p(j) = 2*(j-1)*pi/(np)
       enddo
       dp = p(2)-p(1)

       do i=1,nt
          tcyc(i-nt) = i
          tcyc(i) = i
          tcyc(i+nt) = i
       enddo
       do j=1,np
          pcyc(j-np) = j
          pcyc(j) = j
          pcyc(j+np) = j
       enddo
       ! coefficients for 1st derivative
       cderiv(-2) =  1.0/12.0
       cderiv(-1) = -8.0/12.0
       cderiv(0)  =  0
       cderiv(1)  =  8.0/12.0
       cderiv(2)  = -1.0/12.0

       ! matrix allocations
       if (solve_method == 1) then

          ! FD
          msize = nt*np-1

       else

          ! Spectral
          msize = 4*nts*nps+2*(nts+nps)

          allocate(as(0:nts,0:nps))
          allocate(bs(0:nts,0:nps))
          allocate(cs(0:nts,0:nps))
          allocate(ds(0:nts,0:nps))

       endif

       allocate(mhdfunc(msize))
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
       deallocate(drdr)
       deallocate(dzdtb)
       deallocate(dzdpb)
       deallocate(dzdr)
       deallocate(jac)
       deallocate(bp)
       deallocate(br)
       deallocate(bz)
       deallocate(rp)
       deallocate(rt)
       deallocate(zp)
       deallocate(zt)
       deallocate(fp)
       deallocate(ft)
       deallocate(fpt)
       deallocate(ftp)
       deallocate(mhdfunc)
       deallocate(xfunc)
       deallocate(work)
       deallocate(tcyc)
       deallocate(pcyc)

       initialized = .false.

    endif

  end subroutine le3_alloc

  subroutine le3_solver

    use le3_globals
    implicit none
    integer :: ix,i,j,its,ips
    external :: le3_func,le3_func2

    if (solve_method == 1) then

       do i=1,nt
          do j=1,np
             tb(i,j) = tb(i,j)-t(i)
          enddo
       enddo

       ix=0
       do j=1,np
          do i=1,nt
             if (i+j > 2) then
                ix = ix+1
                xfunc(ix) = tb(i,j)
             endif
          enddo
       enddo

       call hybrd1(le3_func,msize,xfunc,mhdfunc,tol,info,work,nwork)

       !print *, xfunc
       print '(a,1pe12.5)','INFO: (le3) Root accuracy ->',sum(abs(mhdfunc))/size(mhdfunc)

       do i=1,nt
          do j=1,np
             tb(i,j) = tb(i,j)+t(i)
          enddo
       enddo

       open(unit=1,file='out.le3.t',status='replace')
       do i=1,nt
          write(1,10) t(i)
       enddo
       close(1)
       open(unit=1,file='out.le3.p',status='replace')
       do j=1,np
          write(1,10) p(j)
       enddo
       close(1)
       open(unit=1,file='out.le3.tb',status='replace')
       do i=1,nt
          write(1,10) tb(i,:)
       enddo
       close(1)
       open(unit=1,file='out.le3.r',status='replace')
       do i=1,nt
          write(1,10) r(i,:)
       enddo
       close(1)
       open(unit=1,file='out.le3.z',status='replace')
       do i=1,nt
          write(1,10) z(i,:)
       enddo
       close(1)

    else

       as(:,:) = 0.0
       bs(:,:) = 0.0
       cs(:,:) = 0.0
       ds(:,:) = 0.0
       as(1,0) = rmin/rmaj

       ix=0
       do ips=0,nps
          do its=0,nts
             if (its > 0) then 
                ix = ix+1           
                xfunc(ix) = as(its,ips)
             endif
             if (ips > 0 .and. its > 0) then
                ix = ix+1
                xfunc(ix) = bs(its,ips)
             endif
             if (ips + its > 0) then
                ix = ix+1           
                xfunc(ix) = cs(its,ips)
             endif
             if (ips > 0) then
                ix = ix+1
                xfunc(ix) = ds(its,ips)
             endif
          enddo
       enddo

       call hybrd1(le3_func2,msize,xfunc,mhdfunc,tol,info,work,nwork)

       print '(a,1pe12.5)','INFO: (le3) Root accuracy ->',sum(abs(mhdfunc))/size(mhdfunc)

       ix=0
       do ips=0,nps
          do its=0,nts
             if (its > 0) then 
                ix = ix+1           
                print *,'a',ips,its,xfunc(ix)
             endif
             if (ips > 0 .and. its > 0) then
                ix = ix+1
                print *,'b',ips,its,xfunc(ix)
             endif
             if (ips + its > 0) then
                ix = ix+1           
                print *,'c',ips,its,xfunc(ix)
             endif
             if (ips > 0) then
                ix = ix+1
                print *,'d',ips,its,xfunc(ix)
             endif
          enddo
       enddo

       open(unit=1,file='out.le3.t',status='replace')
       do i=1,nt
          write(1,10) t(i)
       enddo
       close(1)
       open(unit=1,file='out.le3.p',status='replace')
       do j=1,np
          write(1,10) p(j)
       enddo
       close(1)
       open(unit=1,file='out.le3.tb',status='replace')
       do i=1,nt
          write(1,10) tb(i,:)
       enddo
       close(1)
       open(unit=1,file='out.le3.r',status='replace')
       do i=1,nt
          write(1,10) r(i,:)
       enddo
       close(1)
       open(unit=1,file='out.le3.z',status='replace')
       do i=1,nt
          write(1,10) z(i,:)
       enddo
       close(1)

    endif

10  format(200(1pe12.5,1x))

  end subroutine le3_solver
  
end module le3_driver
