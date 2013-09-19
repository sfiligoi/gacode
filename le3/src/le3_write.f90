subroutine le3_write

  use le3_globals

  implicit none

  integer :: i,j,k,its,ips, ip, jp
  real :: jacs
  real :: drdtbs,dzdtbs
  real :: drdpbs,dzdpbs
  real, dimension(:), allocatable :: derivvec
  real, dimension(:,:), allocatable :: g, gpp, gtt, gpt
  real, dimension(:,:), allocatable :: rs,zs
  real, dimension(:,:), allocatable :: bpol, btor, bmag
  real, dimension(:,:), allocatable :: dbdt, dbdp
  real, dimension(:,:), allocatable :: bdotgrad, bdotgradB_overB, vdrift_x
 
  print '(a,1pe12.5)','INFO: (le3) Root accuracy ->',sum(abs(yfunc))/size(yfunc)

  print 30,'a_{mn}','b_{mn}','c_{mn}','d_{mn}'
  do ips=0,nps
     do its=0,nts
        print 20,its,ips,as(its,ips),bs(its,ips),cs(its,ips),ds(its,ips)
     enddo
  enddo

  deallocate(t)
  deallocate(tb)
  deallocate(dtbdt)
  deallocate(p)
  deallocate(sinm)
  deallocate(sinn)
  deallocate(cosm)
  deallocate(cosn)

  allocate(t(ntp))  
  allocate(tb(ntp,npp))  
  allocate(dtbdt(ntp,npp))  
  allocate(p(npp))  
  allocate(sinm(ntp,0:nts))  
  allocate(cosm(ntp,0:nts))  
  allocate(sinn(npp,0:nts))  
  allocate(cosn(npp,0:nts)) 

  do i=1,ntp
     t(i) = 2*(i-1)*pi/(ntp)
  enddo
  dt = t(2)-t(1)
  do j=1,npp
     p(j) = 2*(j-1)*pi/(npp)
  enddo
  dp = p(2)-p(1)

  do i=1,ntp
     tb(i,:) = t(i)
  enddo
  dtbdt(:,:) = 1.0

  do its=0,nts
     do i=1,ntp
        sinm(i,its) = sin(its*t(i))
        cosm(i,its) = cos(its*t(i))
     enddo
  enddo
  do ips=0,nps
     do j=1,npp
        sinn(j,ips) = sin(ips*p(j))
        cosn(j,ips) = cos(ips*p(j))
     enddo
  enddo

  allocate(rs(ntp,npp))
  allocate(zs(ntp,npp))
  allocate(g(ntp,npp))
  allocate(gpp(ntp,npp))
  allocate(gtt(ntp,npp))
  allocate(gpt(ntp,npp))

  allocate(btor(ntp,npp))
  allocate(bpol(ntp,npp))
  allocate(bmag(ntp,npp))
  allocate(dbdt(ntp,npp))
  allocate(dbdp(ntp,npp))
  allocate(bdotgrad(ntp,npp))
  allocate(bdotgradB_overB(ntp,npp))
  allocate(vdrift_x(ntp,npp))

  do j=1,npp
     do i=1,ntp

        do ips=0,nps
           do its=0,nts

              tb(i,j) = tb(i,j) + &
                   sinm(i,its)*(bs(its,ips)*sinn(j,ips)+as(its,ips)*cosn(j,ips)) &
                   +cosm(i,its)*(ds(its,ips)*sinn(j,ips)+cs(its,ips)*cosn(j,ips))

              dtbdt(i,j) = dtbdt(i,j) + &
                   its*cosm(i,its)*(bs(its,ips)*sinn(j,ips)+as(its,ips)*cosn(j,ips)) &
                   -its*sinm(i,its)*(ds(its,ips)*sinn(j,ips)+cs(its,ips)*cosn(j,ips))

           enddo
        enddo

        call le3_rz(tb(i,j),&
             p(j),&
             rs(i,j),&
             zs(i,j),&
             drdtbs,&
             drdpbs,&
             dzdtbs,&
             dzdpbs,&
             jacs)

        g(i,j)   = 1.0/rmin * dtbdt(i,j) * jacs

        gpp(i,j) = rs(i,j)**2 + (drdpbs + drdtbs * dtbdp(i,j))**2 &
             + (dzdpbs + dzdtbs * dtbdp(i,j))**2

        gtt(i,j) = (drdtbs**2 + dzdtbs**2) * dtbdt(i,j)**2

        gpt(i,j) = drdtbs * dtbdt(i,j) * (drdpbs + drdtbs * dtbdp(i,j)) &
             + dzdtbs * dtbdt(i,j) * (dzdpbs + dzdtbs * dtbdp(i,j)) 

     enddo
  enddo
 
  btor(:,:) = 1.0/(rs * g)
  
  bpol(:,:) = 1.0/g * (gpt/sqrt(gtt) + iota * sqrt(gtt))
  
  bmag(:,:) = 1.0/g * sqrt(gpp + 2.0*iota*gpt + iota**2 * gtt)

  ! db/dtheta
  allocate(derivvec(0:ntp-1))
  do i=1,ntp-1
     derivvec(i) = -0.5*(-1)**i/tan(0.5*t(i+1))
  enddo
  dbdt(:,:) = 0.0
  do j=1,npp
     do i=1,ntp
        do ip=1,ntp
           k = ip-i
           if(k < 0) then
              k = k + ntp
           endif
           dbdt(i,j) = dbdt(i,j) + derivvec(k) * bmag(ip,j)
        enddo
     enddo
  enddo
  deallocate(derivvec)
  
  ! db/dphi
  allocate(derivvec(0:ntp-1))
  do i=1,npp-1
     derivvec(i) = -0.5*(-1)**i/tan(0.5*p(i+1))
  enddo
  dbdp(:,:) = 0.0
  do j=1,ntp
     do i=1,npp
        do ip=1,npp
           k = ip-i
           if(k < 0) then
              k = k + npp
           endif
           dbdp(j,i) = dbdp(j,i) + derivvec(k) * bmag(j,ip)
        enddo
     enddo
  enddo
  deallocate(derivvec)
           
  
  ! bhat dot grad = bdotgrad * (iota d/dt + d/dp)  
  bdotgrad(:,:) = 1.0/(bmag * g)
  
  ! (bhat dot grad B)/B
  bdotgradB_overB(:,:) = bdotgrad * (iota * dbdt + dbdp) / bmag
  
  ! bhat cross grad B dot grad psi
  vdrift_x(:,:) = iota/(bmag * g**2) &
       * (-dbdt * (gpp + iota * gpt) + dbdp * (gpt + iota*gtt))

  open(unit=1,file='out.le3.t',status='replace')
  write(1,40) t(:)
  close(1)

  open(unit=1,file='out.le3.p',status='replace')
  write(1,40) p(:)
  close(1)

  open(unit=1,file='out.le3.tb',status='replace')
  do j=1,npp
     write(1,40) tb(:,j)-t(:)
  enddo
  close(1)
  open(unit=1,file='out.le3.r',status='replace')
  do i=1,ntp
     write(1,10) rs(i,:)
  enddo
  close(1)
  open(unit=1,file='out.le3.z',status='replace')
  do i=1,ntp
     write(1,10) zs(i,:)
  enddo
  close(1)
  open(unit=1,file='out.le3.b',status='replace')
  write(1,40) bmag(:,:)
  close(1)

  deallocate(rs)
  deallocate(zs)
  deallocate(g)
  deallocate(gpp)
  deallocate(gtt)
  deallocate(gpt)
  deallocate(btor)
  deallocate(bpol)
  deallocate(bmag)
  deallocate(dbdt)
  deallocate(dbdp)
  deallocate(bdotgrad)
  deallocate(bdotgradB_overB)
  deallocate(vdrift_x)
  deallocate(t)
  deallocate(tb)
  deallocate(dtbdt)
  deallocate(p)
  deallocate(sinm)
  deallocate(sinn)
  deallocate(cosm)
  deallocate(cosn)

10 format(200(1pe12.5,1x))
20 format('(',i2,',',i2,'):',2x,4(1pe14.7,1x))
30 format(t15,4(a,9x))
40 format(1pe12.5)

end subroutine le3_write
