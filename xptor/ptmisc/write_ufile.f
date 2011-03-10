      subroutine write_ufile(lprint_ufile,imodel)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c     This routine writes out results in ASCII ufile format
c     lprint_ufile(1)=1 for Te, (2) for Ti, (3) for ne
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      implicit none
c
      include 'mpif.h'
      include '../inc/glf.m'
      include '../inc/tport.m'
      include '../inc/input.m'
c
      character xname*8,unitx*8,fname*8,unitf*8
      character yname*8,unity*8
      character cid*20, desc*75
      character*26 cidass(1)
      integer imodel
      integer j, idchar, idchart, idchars, ichar0, ichar1, ichar2
      integer iufout, ier, iufmsg,inass, iass(1)
      integer lprint_ufile(5)
      real*8 f_var(jmaxm,1)
c
      xname='RHO'
      unitx='none'
      yname='time'
      unity='sec'
      inass=0
      iufout=50
      iufmsg=6
      if(imodel.eq.82) then
        desc='XPTOR simulation using TGLF-APS09 model'
      elseif(imodel.eq.81) then
        desc='XPTOR simulation using GLF23 model'
      endif
c
      idchar=0
      do j=1,10
       if (tok(j:j) .ne. ' ') then
         idchar = idchar + 1
       else
         goto 10
       endif
      enddo
 10   continue
      idchart=idchar
c
      idchar=0
      do j=1,6
        if (shot(j:j) .ne. ' ')then
          idchar = idchar + 1
        else
          goto 20
        endif
      enddo
 20   continue
      idchars=idchar
c
      cid(1:idchart)=tok(1:idchart)
      ichar0=idchart+1
      ichar1=idchart+1+idchars
      ichar2=ichar1+1
      cid(ichar0:ichar1)=shot(1:idchars)
      cid(ichar1:20)='_2d.sim'
      write(*,*) 'cid = ',cid
c
      open(iufout,file=tok(1:idchart)//shot(1:idchars)//
     &    '_2d.sim',status='unknown')
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      if(lprint_ufile(1).gt.0) then
        fname='TE'
        unitf='eV' 
        do j=1,jmaxm
          f_var(j,1)=1.D3*te_m(j)
        enddo
c
        call putu2d(cid,desc,inass,iass,cidass,
     &            xname,unitx,yname,unity,fname,unitf,
     &            jmaxm,1,rho,xp_time,f_var,jmaxm,
     &            iufout,iufmsg,ier)
      endif
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      if(lprint_ufile(2).gt.0) then
        fname='TI'
        unitf='eV' 
        do j=1,jmaxm
          f_var(j,1)=1.D3*ti_m(j)
        enddo
c
        call putu2d(cid,desc,inass,iass,cidass,
     &            xname,unitx,yname,unity,fname,unitf,
     &            jmaxm,1,rho,xp_time,f_var,jmaxm,
     &            iufout,iufmsg,ier)
      endif



c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      close(iufout)
c
      return
      end
c**********************************************************************
      subroutine putu2d(cid,desc,inass,iass,cidass, 
     >                  xname,unitx,yname,unity,fname,unitf, 
     >                  nx,ny,x,y,f,n1d,                       
     >                  iufout,iufmsg,ier)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c     This subroutine writes ITER PDB ufiles for 2d data at 1 time slice
c
c               CID       -  SOURCE IDENTIFIER,  C*20           (INPUT) 
c               CODE      -  SHORT  DESCRIPTION  C*8            (INPUT)
c               DESC      -  COMMENT             C*75           (INPUT)
c               NASS      -  NUMBER OF ASS.SCALARS              (INPUT)
c               IASS()    -  ASS.SCALARS                        (INPUT)
c               CIDASS()  -  LABEL FOR SCALARS   C*26           (INPUT)
c               XNAME     -  INDEPENDENT 1 NAME  C*8            (INPUT)
c               YNAME     -  INDEPENDENT 2 NAME  C*8            (INPUT)
c               FNAME     -  FUNCTION NAME       C*8            (INPUT)
c               UNITX     -  INDEP 1 UNITS       C*8            (INPUT)
c               UNITY     -  INDEP 2 UNITS       C*8            (INPUT) 
c               UNITF     -  FUNCTION UNITS      C*8            (INPUT)
c               NX        -  NUMBER OF POINTS: INDEP X          (INPUT) 
c               NY        -  NUMBER OF POINTS: INDEP Y          (INPUT)
c               X(NX)     -  INDEPENDENT LABEL X (=RHO)         (INPUT)
c               Y(NY)     -  INDEPENDENT LABEL Y (=TIME)        (INPUT)
c               F(N1D,NY) -  FUNCTION POOL                      (INPUT)
c               N1D       -  1ST DIMENSION OF F (>=NX)          (INPUT) 
c               IUFOUT    -  UNIT NUMBER FOR 2D DATABASE FILE   (INPUT)
c               IUFMSG    -  UNIT NUMBER FOR MESSAGES           (INPUT)
c               IER       -  RETURN CODE                       (OUTPUT)
c                             = 0 :  SUCCESSFUL                  
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c       
      character   cid*20, desc*75   
      character   xname*8,unitx*8, fname*8,unitf*8
      character   yname*8,unity*8          
      dimension   x(nx), y(ny), f(n1d,*)  
      dimension   iass(*)           
      character*26 cidass(*)
      parameter (nheadl=9)
      character*44 headl(nheadl)
      data        (headl(l),l=1,nheadl ) / 
     1             ';-SHOT #- SHOT IDENTIFICATION               ', 
     2             ';-SHOT DATE-                                ',
     3             ';-NUMBER OF ASSOCIATED SCALAR QUANTITIES-   ', 
     4             ';-INDEPENDENT VARIABLE LABEL: X-            ',
     5             ';-INDEPENDENT VARIABLE LABEL: Y-            ', 
     6             ';-DEPENDENT VARIABLE LABEL-                 ',
     7             ';-STATUS FLAG                               ',
     8             ';-# OF X PTS-                               ', 
     9             ';-# OF Y PTS- X,Y,F(X,Y) DATA FOLLOW:       ' /
      character*44 headls(2)
      data        (headls(l),l=1,2      ) / 
     1             ';-SCALAR, LABEL FOLLOWS                     ',
     2             '                                            ' /
c
      ier=0
      write(iufout,7002) cid,headl(1),headl(2),inass,headl(3) 
      if(inass.gt.0) then
        do l=1,inass
          write(iufout,7003) iass(l),headls(1),cidass(l),headls(2)
        enddo
      endif
      write(iufout,7004) xname,unitx,headl(4),yname,unity,headl(5),
     &                   fname,unitf,headl(6),headl(7),
     &                   nx,headl(8),ny,headl(9)
c
      write(iufout,7010) x
      write(iufout,7010) y
      write(iufout,7010) ((f(j,l),j=1,nx),l=1,ny)
      write(iufout,*)
     &     ';----END-OF-DATA-----------------COMMENTS:----------'
      write(iufout,'(1x,a)') desc
      write(iufout,7000)
c
 7000 format('****************************************'/  
     >       '****************************************') 
 7002 format(1x,a,t32,a/t32,a/1x,i4,t32,a)
 7003 format(4x,i4,t32,a/
     &       4x,a,t32,a)
 7004 format(1x,a,1x,t22,a,t32,a/
     &       1x,a,1x,t22,a,t32,a/
     &       1x,a,1x,t22,a,t32,a/
     &       1x,'0  ',t32,a/
     &       1x,i10,t32,a/1x,i10,t32,a)
 7010 format(1x,1p,6e13.6)
c
      return
      end
