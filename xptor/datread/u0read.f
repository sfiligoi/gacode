c@u0read.f
c  Read in data in ITER/ITPA PDB 0d ufiles
c  idat = 0 ITER PDB 0D ufile format
c  idat = 1 new IPTA PDB format (163 entries) 3/27/07
c           Note: must remove any slashes
c  from show0d.f from Stan Attenberger, ORNL, June 95.
c  adapted by Joop Konings 04/95, J. Kinsey 9/06
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine u0read(idat,udfile,name,value,ifail)
c
      implicit none
c
      integer ifail, nun, mxshot, mxrow, mxcol, nnames, i, ishot
      integer jdone, inames, nline, jcol
      integer idat, mrow, ju, imax
      character*100 udfile
      character*67 name(200)
      character*11 value(200)
      character*80 line
      character*132 line2
      data nun/10/, mrow/25/
c     data mxshot/50/,mxrow/7/
      data mxshot/1/,mxrow/16/,mxcol/12/,nnames/200/
      include 'mpif.h'
      include '../inc/glf.m'
c
c     open(unit=nun,file='readnames_0d.dat',status='old',err=888)
c     read(nun,'(i5)') nnames
c     read(unit=nun,fmt='(a67)',end=100,err=888) (name(i),i=1,nnames)
c     close(nun)
      if(i_proc.eq.0) write(*,500) udfile
      open(unit=nun,file=udfile,status='old',err=888)
c
      if(idat.eq.1.or. idat.eq.3) then
      imax=163
      do ishot=1,mxshot
        do inames=1,mrow
          read(nun, *,end=101 ) (value(ju),ju=1,imax)
c          write(*,*) ( value(ju),ju=1,imax )
        enddo
      enddo
c
      elseif(idat.eq.2) then
c      write(*,*) 'mxshot = ',mxshot,mrow
      imax=164
      do ishot=1,mxshot
        do inames=1,mrow
          read(nun, *,end=101 ) (value(ju),ju=1,imax)
c          write(*,*) ( value(ju),ju=1,imax )
        enddo
      enddo
c
      else
c
      do ishot=1,mxshot
        jdone=0
c       do inames=1,nnames,mxrow
        do inames=1,mxrow
          read(unit=nun,fmt='(a80)',end=101) line
          nline=80
          do jcol=1,mxcol
            jdone=jdone+1
            call part(line,nline,value(jdone))
c            write(*,*) inames,jdone,'value = ',value(jdone)
            if(jdone.eq.nnames)go to 50
          enddo
        enddo
 50     continue
      enddo
      endif
c
c*****normal exit
c
 101  return
c
c*****error exit
c
 888  continue
      ifail=1
      return
c
 100  write(*,*) 'fatal: expected ',nnames,' lines in names_0d'
      goto 888
 500  format('0D file =',a50)
      end
