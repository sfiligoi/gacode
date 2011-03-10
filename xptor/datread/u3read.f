c@u3read.f   17 Nov 1994   Glenn Bateman, PPPL
c Adapted by J.A.Konings to read experimental data without extrapolation
c to 51 pts grid 7/07/95
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine u3read ( cdir, cfile, kufile,
     & ppxu, kkxu, ppyu, pdata, kerr )
c
c   This routine reads 2-D ASCII U-files.
c
c  cdir      = directory name where U-file is found             (input)
c  cfile     = U-file name                                      (input)
c  kufile    = input unit number to read U-file (default 70)    (input)
c  kkxu      = number of elements in x array                    (output)
c  pxu(jx)   = 1-D independent variable x array from U-file     (output)
c  ppyu       = time of desired profile				(input)
c  parray(jx,jy) = 2-D array from U-file                        (output)
c  
c  kerr      = error indicator                            (input/output)
c
c..On input, set kerr < 0 if you want the program to stop on error.
c
c..If kerr .ge. 0 on input,
c    the routine will return with the following error flag set:
c
c  kerr      = 0 for normal output
c            = 1 if the U-file cannot be opened or read
c            = 2 if there are too many elements in the arrays
c                (that is, if kxu > kxdim)
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
      implicit none
c
      integer jmaxm,kmaxm
      parameter (jmaxm=301,kmaxm=2000)
      real*8 parray(jmaxm,kmaxm), pxu(jmaxm), pyu(kmaxm), zu(6)
      real*8 pdata(jmaxm), ppxu(jmaxm), ppyu
      real*8 pxu1, pxu0, parray1, parray0
c
      integer kxdim, kxu, kydim, kyu, kufile, kdchar, kfchar, kerr
c
      integer  idchar, ifchar, ichar, ichmax, iprint, ierr
     &  , iufile, iscalar, iu, iumax
     &  , j, jx, jy, ju, i, ii, kkxu
c
c
c  iprint  = integer to control amount of diagnostic output
c  ichmax  = maximum allowed length of file and directory names
c
c  idchar  = length of directory name cdir
c  ifchar  = length of U-file name cfile
c  iufile  = unit number for reading UFILE
c  iscalar = number of scalar labels
c  ierr    = input value of kerr
c
      character cdir*(*), cfile*(*), cdfile*130, cline*72
c
c
c..set defaults
c
      ierr   = kerr
c
      kerr   = 0
c
      iprint = 0
c
      ichmax = 60
c
c kxdim is the maximum number of radial points in the 2D ufiles.
c
      kxdim=jmaxm
c
c kydim is the maximum number of time points in the 2D ufiles.
c
      kydim=kmaxm
c
c..clear arrays
c
      do jx=1,kxdim
        pxu(jx)    = 0.D0
      enddo
c
      do jy=1,kydim
        pyu(jy)    = 0.D0
      enddo
c
      do jy=1,kydim
        do jx=1,kxdim
          parray(jx,jy) = 0.D0
        enddo
      enddo
c
c..determine the length of the directory and file names
c  Note:  These names must have no imbedded blanks.
c
      idchar = 0
      do j=1,ichmax
        if ( cdir(j:j)  .ne. ' ' ) then
          idchar = idchar + 1
        else
          go to 10
        endif
      enddo
  10  continue
c
      ifchar = 0
      do j=1,ichmax
        if ( cfile(j:j)  .ne. ' ' ) then
          ifchar = ifchar + 1
        else
          go to 11
        endif
      enddo
  11  continue
c
c..protect against no file name given
c
      if ( ifchar .lt. 1 ) then
        write (*,*)
        write (*,*) 'No file name is given for 2-D U-file'
        write (*,*) ifchar,' = ifchar .lt. 1 in sbrtn u3read'
        write (*,*) idchar,' = idchar'
        write (*,*) cdir,' = cdir'
        write (*,*) cfile,' = cfile'
        write (*,*) ' abort from sbrtn u3read'
        kerr = 1
        if ( ierr .lt. 0 ) stop
        return
      endif
c
      if ( idchar .lt. 1 ) then
c
c..if no directory name is given, then just use the file name
c
        ichar  = ifchar
        cdfile = cfile
c
      else
c
c..concatenate the directory and file name into the character string cdfile
c  Note:
c  On a UNIX system, the directory name must end with '/'
c  On a VMS system, the directory name must end with ']'
c  If the directory name is given but does not end with '/' or ']'
c  then use the UNIX convention and add '/' to the directory name
c
        if (    cdir(idchar:idchar) .ne. '/'
     &    .and. cdir(idchar:idchar) .ne. ']' ) then
c
          ichar = idchar + 1 + ifchar
          cdfile = cdir(1:idchar) // '/' // cfile(1:ifchar)
c
        else
c
          ichar  = idchar + ifchar
c
          cdfile = cdir(1:idchar) // cfile(1:ifchar)
c
        endif
c
      endif
c
      if ( iprint .gt. 0 ) then
        write (*,*) cdfile(1:ichar),' = cdfile'
        write (*,*) ichar,' = idchar'
      endif
c
c..protect against names too long
c
        if ( ichar .gt. 130 ) then
          write (*,*)
          write (*,*) 'directory and file names too long'
          write (*,*) ichar
     &      ,' = idchar + ifchar .gt. 130 in sbrtn u3read'
          write (*,*) ifchar,' = ifchar'
          write (*,*) idchar,' = idchar'
          write (*,*) cdir,' = cdir'
          write (*,*) cfile,' = cfile'
          write (*,*) ' abort from sbrtn u3read'
          kerr = 1
          if ( ierr .lt. 0 ) stop
          return
        endif
c
c..open UFILE
c
      iufile = kufile
      if ( kufile .lt. 1  .or.  kufile .gt. 99 ) iufile = 70
c
      open ( iufile, file=cdfile(1:ichar), status='old', err=90 )
c
c
c..read header
c
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
 100  format (a)
c
c..read number of scalars
c
      read (iufile,102) iscalar
 102  format (i4)
c
c..diagnostic printout
c
      if ( iprint .gt. 8 ) then
        write (*,*)
        write (*,*) cdfile(1:ichar),' = name of dir/file'
        write (*,*) cfile,' = name of UFILE 2'
        write (*,*) iscalar,' = iscalar'
      endif
c
c..read lines of scalars
c
      if ( iscalar .gt. 0 ) then
        do j=1,iscalar*2
          read (iufile,100,err=92) cline
        enddo
      endif
c
c..read variable labels
c
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
c
c..read number of x and y points
c
      read (iufile,104,err=94) kxu, kyu
 104  format (i11)
c
      if ( iprint .gt. 8 ) then
        write (*,*) kxu,' = kxu'
        write (*,*) kyu,' = kyu'
      endif
c
      if ( kxu .gt. kxdim  .or.  kyu .gt. kydim ) then
        write (*,*) 'Too many elements in 2-D U-file arrays'
        write (*,*) kxu, kxdim, kyu, kydim
     &   ,' = kxu, kxdim, kyu, kydim'
        write (*,*)
     &   ' kxu .gt. kxdim .or. kyu .gt. kydim in sbrtn u3read'
        kerr = 1
        if ( ierr .lt. 0 ) stop
        return
      endif
c
c
c..read x and y 1-D arrays
c
      iu = 7
      do jx=1,kxu
        iu = iu + 1
        if ( iu .gt. 6 ) then
          iumax = min ( 6, kxu + 1 - jx )
          read ( iufile, *, err=92 ) ( zu(ju),ju=1,iumax )
          iu = 1
        endif
        pxu(jx) = zu(iu)
      enddo
c
      iu = 7
      do jy=1,kyu
        iu = iu + 1
        if ( iu .gt. 6 ) then
          iumax = min ( 6, kyu + 1 - jy )
          read ( iufile, *, err=92 ) ( zu(ju),ju=1,iumax )
          iu = 1
        endif
        pyu(jy) = zu(iu)
      enddo
c
c
c..read 2-D profile
c
 
      iu = 7
      do jy=1,kyu
        do jx=1,kxu
          iu = iu + 1
          if ( iu .gt. 6 ) then
            iumax = min ( 6, kyu*kxu + 1 - (jy-1)*kxu - jx )
            read ( iufile, 111, err=92 ) ( zu(ju),ju=1,iumax )
            iu = 1
          endif
        parray(jx,jy) = zu(iu)
        enddo
      enddo
  
c
c
c..end of routine
c
      close ( iufile, err=90 )
c
c
c select timeslice in case more measurements are available
c In case of only one timeslice: this one is selected.
c
	kkxu=kxu
		
	i=0
	do ii=1,kyu
	  i=i+1
	  if(pyu(ii).ge.ppyu) goto 95
	enddo
  95    continue
  
         do j=1,kkxu
          ppxu(j)=pxu(j)
          pdata(j)=parray(j,i)
         enddo
        
      return
c
c
c..error conditions
c
  90  continue
c
c
c...JAK more compact output
crew      write (6,'(a22,a20)') 'cannot find 2d UFILE:',cfile
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  92  continue
c
      write (*,*)
      write (*,*) cfile,' = name of UFILE'
      write (*,*) cline
      write (*,*)
     & 'error reading line in UFILE in sbrtn u3read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  94  continue
c
      write (*,*)
      write (*,*) cfile,' = name of UFILE'
      write (*,*) kxu,' = kxu'
      write (*,*) kyu,' = kyu'
      write (*,*)
     & 'error reading kxu or kyu in UFILE in sbrtn u3read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  96  continue
c
      write (*,*)
      write (*,*) cfile,' = name of UFILE'
      write (*,*)
     & 'error reading pxu, pyu, or parray in UFILE in sbrtn u3read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
c
 111  format (1x,6e13.6)
 114  format (1p6e12.4)
 116  format (1p6e13.6)
c
      end
