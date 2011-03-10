c@u2read.f   17 Nov 1994   Glenn Bateman, PPPL
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine u2read ( idat, cdir, cfile, kufile,
     & ppxu, kkxu, ppyu, pdata, kerr, mxtime,
     & xtime, ytime, nty, time_flag, ishift, nsmooth, islice,it,iproc)
c
c   This routine reads 2-D ASCII U-files.
c
c JEK Use either format 111 or 116 for reading ufiles on line 309
c
c  idat      = flag for ufile format                            (input)
c              0 for 111 format, 1 for 116 format
c  cdir      = directory name where U-file is found             (input)
c  cfile     = U-file name                                      (input)
c  kufile    = input unit number to read U-file (default 70)    (input)
c  pxu(jx)   = 1-D independent variable x array from U-file     (output)
c  kkxu      = number of elements in interpolated x array       (output)
c  mxtime    = max number of time pts                           (input)
c  ppyu       = time of desired profile				(input)
c  kyu       = number of elements in y array from U-file        (output)
c  parray(jx,jy) = 2-D array from U-file                        
c  pdata(jjx) = interpolated radial profile from parray 
c  at time index nt_i 						(output)
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
c...adapted by J.A.Konings 95/10/25 to also supply time information 
c...added parameters;
c
c   ytime(51,mxtime) data array with max. mxtime timeslices
c   nty            number of actually read timeslices
c   xtime(mxtime)    array containing the time values
c   time_flag      = 0 don t read all timeslices
c                  = 1 fill arrays xtime, ytime 
c   ishift         number of sample points shifted with respect to ppyu
c...JAK 96/05/16 islice added
c   islice         index of array xtime giving the time of interest
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
      implicit none
c
      integer it, jmaxm, idat, iproc
      integer mxtime, kkxu
      parameter (jmaxm=301)
      real*8 parray(jmaxm,mxtime), pxu(jmaxm), pyu(mxtime), zu(6)
      real*8 pdata(kkxu), ppxu(kkxu), ppyu, der(jmaxm)
      real*8 pxu1, pxu0, parray1, parray0
c
      real*8 xtime(mxtime), ytime(kkxu,mxtime)
      real*8 ytmp(kkxu), ydata(kkxu), sum
      integer time_flag,nty, islice, ishift, nsmooth, ism
c
      integer kxdim, kxu, kydim, kyu, kufile, kdchar, kfchar, kerr
c
      integer  idchar, ifchar, ichar, ichmax, iprint, ierr
     &  , iufile, iscalar, iu, iumax
     &  , j, jx, jy, ju, i, ii
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
c      write(*,*) 'inside u2read ...'
c      write(*,*) 'u2read-ufile = ',cfile
c
c..set defaults
c
      ierr   = kerr
c
      kerr   = 0
c
      iprint = 0
c      iprint=10
c
c      write(*,*) 'mxtime = ',mxtime
c      if(it.eq.29) iprint = 10  ! for diagnostic purposes on a given 2d file
c     7 = NE, 18 = SNBII, 23 = Te
c
      ichmax = 60
c
c kxdim is the maximum number of radial points in the 2D ufiles.
c
      kxdim=jmaxm
c
c kydim is the maximum number of time points in the 2D ufiles.
c
      kydim=mxtime
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
        write (*,*) ifchar,' = ifchar .lt. 1 in sbrtn u2read'
        write (*,*) idchar,' = idchar'
        write (*,*) cdir,' = cdir'
        write (*,*) cfile,' = cfile'
        write (*,*) ' abort from sbrtn u2read'
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
        write (*,*) cdfile(1:ichar),' = cdfile',it
        write (*,*) ichar,' = idchar'
      endif
c
c..protect against names too long
c
        if ( ichar .gt. 130 ) then
          write (*,*)
          write (*,*) 'directory and file names too long'
          write (*,*) ichar
     &      ,' = idchar + ifchar .gt. 130 in sbrtn u2read'
          write (*,*) ifchar,' = ifchar'
          write (*,*) idchar,' = idchar'
          write (*,*) cdir,' = cdir'
          write (*,*) cfile,' = cfile'
          write (*,*) ' abort from sbrtn u2read'
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
        write (*,*) kkxu,' = kkxu'
        write (*,*) kxu,' = kxu'
        write (*,*) kyu,' = kyu'
      endif
c
      if ( kxu .gt. kxdim  .or.  kyu .gt. kydim ) then
        write (*,*) 'Too many elements in 2-D U-file arrays'
        write (*,*) kxu, kxdim, kyu, kydim
     &   ,' = kxu, kxdim, kyu, kydim'
        write (*,*)
     &   ' kxu .gt. kxdim .or. kyu .gt. kydim in sbrtn u2read'
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
        if(iprint.ge.10) write(*,*) jx, pxu(jx),' pxu'
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
        if(iprint.ge.11) write(*,*) jy, pyu(jy),' pyu'
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
            if (idat.eq.1) then
               read(iufile,116,err=92) ( zu(ju),ju=1,iumax )
            else
               read(iufile,111,err=92) ( zu(ju),ju=1,iumax )
            endif
            iu = 1
          endif
        parray(jx,jy) = zu(iu)
        if(iprint.ge.12) write(*,*) jy, jx, parray(jx,jy),' parray'
        enddo
      enddo
c
c
c..end of routine
c
      close ( iufile, err=90 )
c
c
c interpolate the iter data to a kxdim point grid and return 
c only the interpolated quantities.
c
c
c   kkxu=51
c

	do i=1,kkxu
	  ppxu(i)=dfloat(i-1)/dfloat(kkxu-1)
	enddo
	ppxu(1)=ppxu(1)+1.D-10
	ppxu(kkxu)=ppxu(kkxu)-1.D-10
	if(iprint.ge.10) then
	  write(*,*) 'time_flag = ',time_flag
	  write(*,*) 'kkxu = ',kkxu
	  do j=1,kkxu
	    write(*,*) j, ppxu(j), ' ppxu'
	  enddo
	endif
c
	i=0
	do ii=1,kyu
	  i=i+1
          if(iprint.ge.10) write(*,*) 'ppyu = ',ppyu,pyu(ii)
	  if(pyu(ii).ge.ppyu) goto 95
	enddo
  95    continue
        if(pyu(ii).gt.ppyu) then
          i=i-1
          pyu(ii)=pyu(ii-1) ! use time-slice where xp_time >= u2d-time
        endif
        if(iprint.ge.10) write(*,*) 'using data at t = ',pyu(ii)
c
c...JAK shift this timeslice and keep thi i-value of the actual timeslice
c...    Do not shift if timeslice is not available
c...JAK dont use ishift
c
ccc     i = i + ishift
ccc     if ((i.lt.1).or.(i.gt.kyu)) i = i - ishift
        islice = i  
c
c...JAK smooth data (only if no time data is read further on)
c
        if (time_flag.eq.0) then 
          nsmooth = min(nsmooth, min(i-1,kyu-i))
          do j=1,kxu
            sum=0.D0
            do ism=-nsmooth,nsmooth
              sum = sum + parray(j,i+ism)
            end do
            parray(j,i) = sum/(2.D0*dfloat(nsmooth)+1.D0)
          end do
        endif
c
c...JAK structure changed: don t return but continue for reading timedaata
c
        if(kxu.ge.kkxu) then
         do j=1,kkxu
          ppxu(j)=pxu(j)
          pdata(j)=parray(j,i)
         enddo
        if(iprint.ge.10) then
          write(*,*)
          do j=1,kxu
            write(*,150) j, pxu(j), parray(j,i)
          enddo
        endif
c
c could be problems here if pxu(1) not 0 or pxu(kxu) not 1.0
cjak     return
        endif
c 
      if(kxu.lt.kkxu) then
c-->
      if(iprint.ge.10) then
        write(*,*)
        do j=1,kxu
          write(*,150) j, pxu(j), parray(j,i)
        enddo
      endif
c
c      pxu0=0.D0     ! jek: disabled 5/22/00
      pxu0=pxu(1)
      if (pxu(1).gt.0.025) pxu0=1.D-3 ! set to zero if pxu(1) too large
      parray0=parray(1,i)+(parray(1,i)-parray(2,i))
     >        /(pxu(1)-pxu(2))*(pxu0-pxu(1))
c      pxu(1)=pxu0   ! jek: disabled 5/22/00
      if (pxu(1).gt.0.025) pxu(1)=pxu0
      parray(1,i)=parray0
  
c     pxu1=1.D0
      pxu1=pxu(kxu)
      parray1=parray(kxu,i)+(parray(kxu,i)-parray(kxu-1,i))
     >        /(pxu(kxu)-pxu(kxu-1))*(pxu1-pxu(kxu))
c      pxu(kxu)=pxu1   ! jek: disabled 5/22/00
      parray(kxu,i)=parray1
      if(iprint.ge.10) then
        do j=1,kxu
          write(*,155) j, pxu(j), parray(j,i)
        enddo
      endif
c
      call INTER_CSPL(kxu,pxu,parray(1,i),kkxu,ppxu,pdata,der)
c
      if(iprint.ge.10) then
        write(*,*)
        do j=1,kkxu
          write(*,66) j, ppxu(j), der(j), pdata(j)
        enddo
      endif
c      write(*,*) 'u2, okay so far-11c...'
c<--
      endif
  66  format(i2,2x,1p3e12.4,2x,'interpolated pdata')
c--------------------------------------------  
c
      if (time_flag.eq.1) then   !fill timearrays
c
c...    reset arrays
c
c     do jy=1,6
c       write(6,*) (parray(jx,jy),jx=1,6)
c     enddo

      do jy=1,mxtime
        xtime(jy) = 0.D0
        do jx=1,kkxu
          ytime(jx,jy) = 0.D0
        enddo
      enddo
c
        nty=kyu  
c
        if(iproc.eq.0) write(6,*) 'Number of timesamples:', kyu
c
        do i=1,kyu           !start time loop
c
          xtime(i)=pyu(i)
c
c
c...      note that ppxu is already calculated at the time ppyu of
c...      interest.
c
          if(kxu.lt.kkxu) then
c
c...        don t adjust the boundary condition for timeslice that
c...        has already been done
c
            if (i.ne.islice) then
              pxu0=0.D0
              parray0=parray(1,i)+(parray(1,i)-parray(2,i))
     >                /(pxu(1)-pxu(2))*(pxu0-pxu(1))
              pxu(1)=pxu0
              parray(1,i)=parray0
 
              pxu1=1.D0
              parray1=parray(kxu,i)+(parray(kxu,i)-parray(kxu-1,i))
     >                /(pxu(kxu)-pxu(kxu-1))*(pxu1-pxu(kxu))
              pxu(kxu)=pxu1
              parray(kxu,i)=parray1
            endif
c
c...        copy to local array
c
            do j=1,kxu
              ytmp(j) = parray(j,i)
            end do
            call INTER_CSPL(kxu,pxu,ytmp(1),kkxu,ppxu,ydata,der)
            do j=1,kkxu
              ytime(j,i) = ydata(j)
            end do
          else
            do j=1,kkxu
              ytime(j,i) = parray(j,i)
            end do
          endif 
        end do
c
      endif
c--------------------------------------------  
c
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
      write(*,*) 'zu =',zu(ju)
      write (*,*) cline
      write (*,*)
     & 'error reading line in UFILE in sbrtn u2read'
      kerr = 1
c
c...  JAK 96/04/04 close statement added
c
      close(iufile)
c
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
     & 'error reading kxu or kyu in UFILE in sbrtn u2read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  96  continue
c
      write (*,*)
      write (*,*) cfile,' = name of UFILE'
      write (*,*)
     & 'error reading pxu, pyu, or parray in UFILE in sbrtn u2read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
c
 111  format (1x,6e13.6)
 114  format (1p6e12.4)
 115  format (1pe14.6,1p5e13.6)
 116  format (1p6e13.6)
 150  format (2x,i2,2x,1p2e13.6,2x'parray-before')
 155  format (2x,i2,2x,1p2e13.6,2x'parray-after')
c
      end
