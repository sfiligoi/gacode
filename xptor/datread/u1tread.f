c@u1tread.f   derived from u1read.f by G.M.Staebler 6/28/2000
c reads full time seris from 1d u-files.
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine u1tread ( cdir, cfile, kufile, 
     &  kxu ,pxu, parray, nsmooth, kerr )
c
c   This routine reads 1-D ASCII U-files.
c
c  cdir      = directory name where U-file is found             (input)
c  cfile     = U-file name                                      (input)
c  kufile    = input unit number to read U-file (default 70)    (input)
c  kxdim     = maximum number of elements allowed in arrays     (input)
c
c  kxu       = number of time points            		(input)
c  parray     = data time series         			(output)
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
      integer nrmax, ntmax, nsmooth
      parameter(nrmax=51) !max number of radial points in experimental grid
c      parameter(ntmax=801) !max number of time points in experimental grid
      parameter(ntmax=3801) !max number of time points in experimental grid
c
      real*8 parray(ntmax), pxu(ntmax), zu(6)
c
      integer kxdim, kxu, kufile, kerr
c
      integer  idchar, ifchar, ichar, ichmax, iprint, ierr
     &  , iufile, iscalar, iu, iumax
     &  , j, jx, ju, i, ii
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
c kxdim is the maximum number of time slices we can read
c
      kxdim=ntmax

c
c..clear arrays
c
      do j=1,kxdim
        parray(j) = 0.0
        pxu(j)    = 0.0
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
        write (*,*) 'No file name is given for 1-D U-file'
        write (*,*) ifchar,' = ifchar .lt. 1 in sbrtn u1read'
        write (*,*) idchar,' = idchar'
        write (*,*) cdir,' = cdir'
        write (*,*) cfile,' = cfile'
        write (*,*) ' abort from sbrtn u1read'
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
        write (*,*) ichar,' = ichar'
      endif
c
c..protect against names too long
c
        if ( ichar .gt. 130 ) then
          write (*,*)
          write (*,*) 'directory and file names too long'
          write (*,*) ichar
     &      ,' = idchar + ifchar .gt. 130 in sbrtn u1read'
          write (*,*) ifchar,' = ifchar'
          write (*,*) idchar,' = idchar'
          write (*,*) cdir,' = cdir'
          write (*,*) cfile,' = cfile'
          write (*,*) ' abort from sbrtn u1read'
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
c         write(*,*) 'file = ',cdfile(1:ichar)
      open ( iufile, file=cdfile(1:ichar), status='old', err=90 )
c
c..read header
c
      read (iufile,100,err=92) cline
      read (iufile,100,err=92) cline
 100  format (a)
c
c..read number of scalars
c..JAK 96/03/25 error check added
c
      read (iufile,102,err=92) iscalar
 102  format (i4)
c
c..diagnostic printout
c
      if ( iprint .gt. 8 ) then
        write (*,*)
        write (*,*) cdfile(1:ichar),' = name of dir/file'
        write (*,*) cfile,' = name of UFILE 1'
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
c
c..read number of x and y points
c
      read (iufile,104,err=94) kxu
 104  format (i11)
c
      if ( iprint .gt. 8 ) then
        write (*,*) kxu,' = kxu'
      endif
c
      if ( kxu .gt. kxdim  .or.  kxu .lt. 1 ) then
        write (*,*) 'Too many elements in 1-D U-file array'
        write (*,*) kxu, kxdim
     &   ,' = kxu, kxdim'
        write (*,*)
     &   ' kxu .gt. kxdim .or. kxu .lt. 1 in sbrtn u1read'
        kerr = 1
        if ( ierr .lt. 0 ) stop
        return
      endif
c
c
c..read x  1-D array
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
c       write(*,*) jx, pxu(jx), ' pxu-u1tread'
      enddo
c
c
c..read 1-D profile
c
  
cfix      if(cfile(1:3).ne.'d3d') then 

      iu = 7
        do jx=1,kxu
          iu = iu + 1
          if ( iu .gt. 6 ) then
            iumax = min ( 6, kxu + 1 - jx )
            read ( iufile, 111, err=92 ) ( zu(ju),ju=1,iumax )
            iu = 1
          endif
        parray(jx) = zu(iu)
        enddo

cfix      else 
cfix
cfix      iu = 7
cfix        do jx=1,kxu
cfix          iu = iu + 1
cfix          if ( iu .gt. 6 ) then
cfix            iumax = min ( 6, kxu + 1 - jx )
cfix            read ( iufile, *, err=92 ) ( zu(ju),ju=1,iumax )
cfix            iu = 1
cfix          endif
cfix        parray(jx) = zu(iu)
cfix        enddo
  

cfix      endif
      if(nsmooth.ne.0)then
c smooth in time
        do i=2,kxu-1
            parray(i) = 2.D0*parray(i)/3.D0 +
     >      (parray(i-1)+parray(i+1))/6.D0
        enddo
      endif       
c
c
c
c..end of routine
c
      close ( iufile, err=90 )
c
      return
c
c
c..error conditions
c
  90  continue
c
c...JAK more compact output
c      write (6,'(a22,a24)') 'cannot find 1d UFILE:',cfile
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
     & 'error reading line in UFILE in sbrtn u1read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  94  continue
c
      write (*,*)
      write (*,*) cfile,' = name of UFILE'
      write (*,*) kxu,' = kxu'
      write (*,*)
     & 'error reading kxu in UFILE in sbrtn u1read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  96  continue
c
      write (*,*)
      write (*,*) cfile,' = name of UFILE'
      write (*,*)
     & 'error reading pxu, pyu, or parray in UFILE in sbrtn u1read'
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
