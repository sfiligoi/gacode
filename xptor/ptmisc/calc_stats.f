      subroutine calc_stats(iflag,ndata,r,data,model,stdev,relrms,
     &                      offset,datmax)
c
      integer iflag, j, ndata
      real*8 dev, devtot, devsq, expsq, stdev, offset
      real*8 datmax, reldev, reldevsq, relrms
      real*8 r(0:ndata), data(0:ndata), model(0:ndata)
c
      dev=0.D0
      reldev=0.D0
      devtot=0.D0
      devsq=0.D0
      reldevsq=0.D0
      expsq=0.D0
      datmax=0.D0
c
c datmax = maximum experimental data value
c
      do j=0,ndata
        datmax = max(datmax,data(j))
c        write(*,*) j, data(j), model(j)
      enddo
c
      do j=0,ndata
        dev = data(j) - model(j)
        if(iflag.eq.1 .and. (abs(data(j)).lt.1.d-6 .or.
     >     abs(model(j)).lt.1.d-6) ) dev=0.D0
c        if(iflag.eq.1) then
c         write(*,100) j,r(j),data(j),model(j),dev
c        endif
        reldev = dev / datmax
        devtot = devtot + dev
        devsq = devsq  + dev**2.D0
        reldevsq = reldevsq  + reldev**2.D0
        expsq = expsq + data(j)**2.D0
c       write(*,100) j, r(j), data(j), model(j), devsq, expsq
      enddo
c
      stdev = dsqrt( devsq/expsq )
      relrms = dsqrt ( reldevsq / ndata )
      offset = devtot/expsq
c
c      write(*,150) stdev
c      write(*,200) offset
c      write(*,999)
c
 100  format(2x,i2,2x,0p1f4.2,0p6f10.5)
 150  format(/,2x,'sigma-Ti = ',0pf10.5)
 200  format(2x,'f-Ti     = ',0pf10.5)
 999  format(' -----------------------')
c
      return
      end
