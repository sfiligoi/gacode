


    subroutine plot_dump(x,nzones,y,name)
      real, intent(in) :: x(nzones),y(nzones)
      character*(*) name
      integer,intent(IN) :: nzones
      integer iounit
      iounit =12
      call getioun (iounit,iounit)
      open  (unit = iounit  , file = name , status = 'UNKNOWN')
      write(iounit,FMT='("#",a)')name
      do j=1,nzones
         write(iounit,1)x(j),y(j)
      enddo
      write(iounit,FMT='(2x,"  ")')
 1 format(2x,2(1pe14.3,2x))
      call giveupus(iounit)
      close(iounit) 
    end subroutine plot_dump

    subroutine plot_dumpd(x,nzones,y,name)
      real*8, intent(in) :: x(nzones),y(nzones)
      character*(*) name
      integer,intent(IN) :: nzones
      integer iounit
      iounit =12
      call getioun (iounit,iounit)
      open  (unit = iounit  , file = name , status = 'UNKNOWN')
      write(iounit,FMT='("#",a)')name
      do j=1,nzones
         write(iounit,1)x(j),y(j)
      enddo
      write(iounit,FMT='(2x,"  ")')
 1    format(2x,2(1pe14.3,2x))
      call giveupus(iounit)
      close(iounit)  
    end subroutine plot_dumpd



