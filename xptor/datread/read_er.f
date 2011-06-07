c@read_er.f
c jek 17-Jan-11
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c... Reads in Er data
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c  
      subroutine read_er
c
      implicit none
      include '../inc/tport.m'
      include '../inc/data.m'
c
      character*50 erfile
      character cdfile*130
      integer j, mxnr_r
      integer idchar, ifchar, ichar, ichmax, nr_er
      parameter (ichmax=60,mxnr_r=300)
      real*8 rho_er_d(mxnr_r), er_raw_d(mxnr_r)
c
      write(*,*) '** inside read_er **'
c
c...  count characters in Er filename -> ifchar
c
      ifchar = 0
      do j=1,ichmax
        if ( shot(j:j)  .ne. ' ' ) then
          ifchar = ifchar + 1
        else
          go to 10
        endif
      enddo
 10   continue
c
      erfile = shot(1:ifchar)//'.er'
c
      write(*,*) 'cudir = ',cudir
      idchar = 0
      do j=1,ichmax
        if ( cudir(j:j)  .ne. ' ' ) then
          idchar = idchar + 1
        else
          go to 15
        endif
      enddo
 15   continue
c
      ichar = 0
      if (cudir(idchar:idchar) .ne. '/') then
        ichar = idchar + 1 + ifchar
        cdfile = cudir(1:idchar) // '/' // erfile
      else
        ichar  = idchar + ifchar
        cdfile = cudir(1:idchar) // erfile
      endif
      write(*,*) ' Opening Er file = ',erfile
c
c
      open(unit=5,status='old',access='sequential',
     >     file=cdfile)
c
c... count number of lines
c
      nr_er=0
      do j=1,mxnr_r
        read(5,*,end=20)
        nr_er=nr_er+1
      enddo
 20   continue
c
      rewind(unit=5)
      do j=1,nr_er
         read(5,*) rho_er_d(j), er_raw_d(j)
c         write(*,'(i3,0p6f10.5)') j, rho_er_d(j), er_raw_d(j)
      enddo
      rho_er_d(j)=1.0 + 1.d-6
      close(5)
c
      call cub_spline(rho_er_d,er_raw_d,nr_er,rho_d,er_d,nj_d)
      do j=1,nj_d
        write(*,'(i3,0p6f10.5)') j, rho_d(j), rho(j), er_d(j)
      enddo
c
      end

