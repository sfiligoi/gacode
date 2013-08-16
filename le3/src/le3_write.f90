subroutine le3_write

  use le3_globals

  implicit none

  integer :: i,j,its,ips
  real :: jacs
  real :: drdtbs,dzdtbs
  real :: drdpbs,dzdpbs
  real, dimension(:,:), allocatable :: bp
  real, dimension(:,:), allocatable :: bt
  real, dimension(:,:), allocatable :: rs
  real, dimension(:,:), allocatable :: zs

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

  allocate(bt(ntp,npp))
  allocate(bp(ntp,npp))
  allocate(rs(ntp,npp))
  allocate(zs(ntp,npp))

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

        bt(i,j) = rmin*dtbdt(i,j)/jacs
        bp(i,j) = iota*rmin/(rs(i,j)*jacs)*sqrt(drdtbs**2+dzdtbs**2)

     enddo
  enddo

  open(unit=1,file='out.le3.t',status='replace')
  do i=1,ntp
     write(1,10) t(i)
  enddo
  close(1)
  open(unit=1,file='out.le3.p',status='replace')
  do j=1,npp
     write(1,10) p(j)
  enddo
  close(1)
  open(unit=1,file='out.le3.tb',status='replace')
  do i=1,ntp
     write(1,10) tb(i,:)
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
  do i=1,ntp
     write(1,10) sqrt(bt(i,:)**2+bp(i,:)**2)
  enddo
  close(1)

10 format(200(1pe12.5,1x))
20 format('(',i2,',',i2,'):',2x,4(1pe14.7,1x))
30 format(t15,4(a,9x))

end subroutine le3_write
