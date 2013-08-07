module le3_driver

  implicit none

  public :: le3_alloc, le3_solver

  real, dimension(:), allocatable, private :: mhdfunc
  real, dimension(:), allocatable, private :: xfunc
  integer, private :: msize
  integer, private :: info
  integer, private :: nwork
  real, dimension(:), allocatable, private :: work
  logical, private :: initialized = .false.

contains

  subroutine le3_alloc(flag)
    use le3_globals
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    integer :: it, ip

    if (flag == 1) then
       if(initialized) return

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

       do it=1,nt
          t(it) = 2*(it-1)*pi/(nt)
       enddo
       dt = t(2)-t(1)
       do ip=1,np
          p(ip) = 2*(ip-1)*pi/(np)
       enddo
       dp = p(2)-p(1)

       do it=1,nt
          tcyc(it-nt) = it
          tcyc(it) = it
          tcyc(it+nt) = it
       enddo
       do ip=1,np
          pcyc(ip-np) = ip
          pcyc(ip) = ip
          pcyc(ip+np) = ip
       enddo
       ! coefficients for 1st derivative
       cderiv(-2) =  1.0/12.0
       cderiv(-1) = -8.0/12.0
       cderiv(0)  =  0
       cderiv(1)  =  8.0/12.0
       cderiv(2)  = -1.0/12.0

       ! matrix allocations
       msize = (nt-1)*np
       allocate(mhdfunc(msize))
       allocate(xfunc(msize))
       nwork = (msize*(3*msize+13))/2 * 10
       allocate(work(nwork))

       initialized = .true.

    else
       if(.NOT. initialized) return

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

    end if

  end subroutine le3_alloc

  subroutine le3_solver
    use le3_globals
    implicit none
    integer :: ir,it,ip
    external :: le3_func

    do it=1,nt
       do ip=1,np
          tb(it,ip) = tb(it,ip) - t(it)
       enddo
    enddo

    ir=1
    do ip=1,np
       do it=2,nt
          xfunc(ir) = tb(it,ip)
          ir = ir+1
       enddo
    enddo
    
    call hybrd1(le3_func,msize,xfunc,mhdfunc,tol,info,work,nwork)
    
    !print *, xfunc
    print '(a,1pe12.5)','INFO: (le3) Root accuracy ->',sum(abs(mhdfunc))/size(mhdfunc)
    
    do it=1,nt
       do ip=1,np
          tb(it,ip) = tb(it,ip) + t(it)
       enddo
    enddo

    open(unit=1,file='out.le3.t',status='replace')
    do it=1,nt
       write(1,10) t(it)
    enddo
    close(1)
    open(unit=1,file='out.le3.p',status='replace')
    do ip=1,np
       write(1,10) p(ip)
    enddo
    close(1)
    open(unit=1,file='out.le3.tb',status='replace')
    do it=1,nt
       write(1,10) tb(it,:)
    enddo
    close(1)
    open(unit=1,file='out.le3.r',status='replace')
    do it=1,nt
       write(1,10) r(it,:)
    enddo
    close(1)
    open(unit=1,file='out.le3.z',status='replace')
    do it=1,nt
       write(1,10) z(it,:)
    enddo
    close(1)
    
10  format(200(1pe12.5,1x))
    
  end subroutine le3_solver
  
end module le3_driver
