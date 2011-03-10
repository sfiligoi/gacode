      subroutine readufiles(tok,shot,phase,cudir,mxgrid,ismooth,
     &                    btscale,xp_time,endtime_pt,time_series,
     &                    iexp_exch,iexp_q,
     &                    p_glob_exp,gtauth_exp,gtautot_exp,iproc)
c************************************************************************
c Reads experimental profiles and splines from ITER PDB 
c 1D and 2D Ufiles.   
c
c The following variables denote the mass and charge of the
c main, beam, and impurity ion species:
c     pgasa      plasma working gas atomic number
c     pgasz      plasma working gas charge
c     bgasa      beam species atomic number
c     bgasz      beam species charge
c     pimpa      impurity species atomic number
c     pimpz      impurity species charge
c 
c Global quantities vs. time in 1D ufiles:
c
c     amin       horizontal minor radius (m)
c     bt         vacuum toroidal field (T) at RGEO
c     delta      triangularity 
c     indent     indentation 
c     ip         plasma current (A)
c     kappa      elongation
c     li         induction 
c     q95        plasma safety factor evaluated at the flux surface that 
c                encloses 95% of the total poloidal flux
c     rgeo       plasma geometrical major radius (m)
c     wtot       total stored energy
c     zeff       ave zeff
c
c Profile quantities vs. rho, time in 2D ufiles:
c
c     chie         experimentally inferred chi_e
c     chii         experimentally inferred chi_i
c     curnbi       NBI currents??
c     curtot       Total current??
c     grho1        < abs(grad rho)>
c     grho2        <(grad rho)**2>
c     ne           electron density
c     nfast        fast ion density
c     nimp         impurity density
c     nm1          main ion density
c     qnbie        power from NBI to electrons
c     qnbii        power from NBI to ions
c     qohm         ohmic heating power
c     qrad         radiated power
c     rmajor       avg major radius of flux surface
c     rminor       avg minor radius of flux surface
c     snbie        electron particle source from NBI
c     snbii        ion particle source from NBI
c     swall        thermal ion particle source due to ionisation
c     surf         flux surface surface area
c     te           T_e
c     ti           T_i
c     volume       flux surface volume
c     zeffr        Z_eff profile
c     ptot         total pressure profile
c
c The data is assumed to be in directory cudir with a name of the form:
c machineNdshotnum.FIELD
c
c The ufiles are assumed to have names beginning with the tokamak,
c then 1d or 2d then the shot number. For example,
c
c /u/usernumber/iterpdb/d3d/69627/d3d1d69627.AMIN
c /u/usernumber/iterpdb/d3d/69627/d3d2d69627.TI
c
c A number of global quantites (e.g. Wtot) are often provided in a
c single 0D file. For example,
c
c /u/usernumber/iterpdb/d3d/69627/d3d_69627_0d.dat
c
c The 0D file is not in ufile format. It is in either in column format
c (idatzero=0) or using a column delimited format (idatzero=1, 2, or 3).
c
c If a file is not present, the relevant variables are either not
c set or some default value (e.g. zeff=1.0)
c
c nrmax = max number of radial grid points in experimental data
c         changed array declarations using nrmax -> nj
c nr   = actual number of radial grid points in experimental data
c rho_xp  = radial variable read from experimental data

c u3dx = radial positions of measured quantities like TEXP, TEXBR, etc.
c u3dy = values of measured quantities like TEXP, TEXBR, etc.
c nx3  = number of points in u3dx and u3dy.
c
c************************************************************************
      implicit none
c
c      include 'mpif.h'
      include '../inc/input.m'
      include '../inc/model.m'
      include '../inc/data.m'
c
      integer iproc,mxgrid,idat, nrmax, ntmax, n1d, n2d, n3d, nur
c
      parameter(nrmax=301) !max number of radial points in experimental grid
      parameter(ntmax=3800) !max number of time points in experimental grid
      parameter(n1d=13)   !number of 1d fields read;
      parameter(n2d=48) 
      parameter(n3d=6)   !number of fields read for experimental data
c      include '../inc/ut.m'  ! need ntmax defined first
c
      character*6 u0phase
      character*7 fields1d(n1d), fields2d(n2d), fields3d(n3d)
      character*50 ufile, erfile
      character*100 udfile
      character cdfile*130
      character*67 u0names(200)
      character*11 u0values(200)
      character*(40) shot
      character*(10) tok
      character*(6) phase
      character*50 cudir
c
      integer i, ig, is, i2d, ierr, time_channel
     &  , j, i1d, id, ifchar, ichar, ifail, k, jstart, jend
     &  , nnexp_d, ntexp_d, ntixp_d, jj, jn, niterdb, linblnk, ier
     &  , ichmax, nr_er, iread, nsmooth, time_flag, itt, itime_shift
      integer mxnr_r
      integer nx3(n3d)
      parameter(mxnr_r=300)
c
      real*8 rfix, sfix, q01_exp, qa1_exp, zeff_t, zbrac, alamda
     &  , tau_e, ds, z2, alfj1_exp, xxq, alpha, sum, cur0, cursum
     &  , rhostar, shift_shaf, dummy1, aomega
     &  , zhalf, zthird, z2thrd, z4thrd, z3qtr
      real*8 xp_time, endtime_pt, btscale, gtauth_exp,gtautot_exp
c
      real*8 rgeo_0d, rmin_0d, bt_0d, ip_0d, nebar_0d, pnbi_0d
     &  , ti0_0d, te0_0d, pich_0d, pech_0d, q0_0d, q95_0d
     &  , kap_0d, delt_0d, wth_0d, wtot_0d
     &  , rho_norm_loc, aj
      real*8 pi_m, kevdsecpmw, p_glob_exp
      real*8 lnlam,eta,pohsum,pradsum,drm,dvoldr_p,dvoldr_m
     &  , prfesum, prfisum, conv_torq
c
      real*8 xnexp_d(nrmax), ynexp_d(nrmax)
     &  , xnexpeb_d(nrmax), ynexpeb_d(nrmax)
     &  , xtexp_d(nrmax), ytexp_d(nrmax)
     &  , xtexpeb_d(nrmax), ytexpeb_d(nrmax)
     &  , xtixp_d(nrmax), ytixp_d(nrmax)
     &  , xtixpeb_d(nrmax), ytixpeb_d(nrmax)
c
      integer nr, nut, time_series, ismooth, iexp_exch, iexp_q
      real*8 rho_xp(nrmax+1), ut(ntmax)
      real*8 u1d_t(ntmax,n1d), u2d_t(nrmax,ntmax,n2d)
      real*8 rhox(0:nrmax-1)
c
      real*8 xtime_d(ntmax), ytime_d(nrmax,ntmax)
      real*8 rho_er_d(nrmax), er_raw_d(nrmax)
      real*8 ur(nrmax),rho_u(nrmax)
      real*8 pxp(nrmax),par(nrmax,501)
      real*8 u1d(n1d),u2d(nrmax,n2d),u1d_int(ntmax)
      real*8 u3dx(nrmax,n3d),u3dy(nrmax,n3d)
c
      data fields1d/'AMIN   ','BT     ','DELTA  ','INDENT ',
     .	'IP     ','KAPPA  ','LI     ','Q95    ','RGEO   ',
     .  'WTOT   ','ZEFF   ','POHM   ','VSURF  '/
c
      data fields2d/'CHIE   ','CHII   ','CURNBI ','CURTOT ',
     .	'GRHO1  ','GRHO2  ','NE     ','NFAST  ','NIMP   ','NM1    ',
     .	'QNBIE  ','QNBII  ','QOHM   ','QRAD   ','RMAJOR ','RMINOR ',
     .	'SNBIE  ','SNBII  ','SURF   ','ZEFFR  ','VOLUME ','TI     ',
     .	'TE     ','DELTAR ','INDENTR','KAPPAR ','Q      ','QEI    ',
     .	'VROT   ','SWALL  ','DWER   ','DWIR   ','QECHE ','QICRHE ',
     .	'QECHI  ','QICRHI ','QWALLE ','QWALLI ','QFUSE  ','QFUSI  ',
     .	'DNER','NM2 ','NM3 ','PTOTR','QLHE ','TORQ',
     .  'CHIENEO','CHIINEO'/
c
      data fields3d/'NEXP   ','NEXPEB ',
     >  'TEXP   ','TEXPEB ','TIXP   ','TIXPEB '/
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
       do i=1,mxgrid+1
         do j=1,n2d
            u2d(i,j)=0.D0
         enddo
       enddo
c
c...Misc switches
c
      idat=0
      rho_norm_loc=1.D0
      conv_torq=1.0
      nsmooth = ismooth
      jmaxm=mxgrid
      pi_m=3.1415926D0
      kevdsecpmw=1.6022D-19*1.D3*1.D-6
      zhalf = 1.D0/2.D0
      zthird = 1.D0/3.D0
      z2thrd = 2.D0/3.D0
      z4thrd = 4.D0/3.D0
      z3qtr = 3.D0/4.D0
c
c... Setup grids and constants
c
      nr=mxgrid+1
      do j=0,jmaxm
        aj=dfloat(j)
        rhox(j)=aj/dfloat(jmaxm)
        rho_d(j+1)=rhox(j)
      enddo
      rhox(0)=1.D-6
      rhox(jmaxm)=rhox(jmaxm)-1.D-6
      rho_d(1)=rhox(0)
      rho_d(jmaxm+1)=rhox(jmaxm)
c
c... Special shots needing different format reads
c
      is=linblnk(shot,10)
      if(shot(1:is).eq.'35156') idat=0
      if(shot(1:is).eq.'35171') idat=0
      if(shot(1:is).eq.'37379') idat=0
      if(shot(1:is).eq.'37944') idat=0
      if(shot(1:is).eq.'38285') idat=0
      if(shot(1:is).eq.'38287') idat=0
      if(shot(1:is).eq.'46664') idat=0
      if(shot(1:is).eq.'52009') idat=0
      if(shot(1:is).eq.'52014') idat=1
      if(shot(1:is).eq.'52022') idat=0
      if(shot(1:is).eq.'52025') idat=1
      if(shot(1:is).eq.'52096') idat=0
      if(shot(1:is).eq.'50844') idat=1
      if(shot(1:is).eq.'53299') idat=1
      if(shot(1:is).eq.'52014') idat=0
      if(shot(1:is).eq.'52015') idat=0
      if(shot(1:is).eq.'52979' .and. idatzero.eq.0) idat=1
      if(shot(1:is).eq.'55935') idat=1
      if(shot(1:is).eq.'57987') idat=0
      if(shot(1:is).eq.'58159') idat=0
      if(shot(1:is).eq.'58323') idat=0
      if(shot(1:is).eq.'60933') idat=0
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... Loop over the 2d files first:
c
      do i=1,n2d
c
c       compose the file name:
         ig=len_trim(tok)
         is=len_trim(shot)         
c         i2d=linblnk(fields2d(i),7)
         i2d=len_trim(fields2d(i))
         ufile=tok(1:ig)//'2d'//shot(1:is)//'.'//fields2d(i)(1:i2d)
c         write(*,*) i, ufile
c
         ierr=0
         time_flag=0
         if (i.eq.time_channel) then
           time_flag=1
           if(iproc.eq.0) write(6,*) 'Reading timedata for channel=',i
         endif
         if(time_series.eq.0)then
c          write(*,*) 'ufile = ',ufile
c          write(*,*) 'u2d = ',u2d(1,i)
          call u2read(idat,cudir,ufile,88,rho_xp,nr,xp_time,u2d(1,i),
     &               ierr,ntmax,xtime_d, ytime_d, ntime_d, time_flag,
     &               itime_shift,nsmooth,islice_d,i,iproc)
         else
          if(i.eq.1 .and. iproc.eq.0) write(6,60)
          do k=1,mxgrid+1
            pxp(k) = 0.D0
            do j=1,ntmax
              u2d_t(k,j,i) = 0.D0
            enddo
          enddo
          call u2tread(cudir,ufile,88,nur,nut,ur,ut,par,nsmooth,i,ierr)
          if(ierr.eq.0)then
           if(nut.lt.2)then
             if(iproc.eq.0)write(*,*)"error nut < 2 ",nut
             return
           endif
           if(ut(1).gt.endtime_pt.or.ut(nut).lt.xp_time
     >      .or.ut(nut).lt.endtime_pt.or.ut(1).gt.xp_time)then
             if(iproc.eq.0)
     >       write(*,*)"error: data time range is ",ut(1),ut(nut)
             stop
           endif
           do k=1,nur
             rho_u(k) = ur(k)
           enddo
c          rho_u(1) = 0.D-10
c          rho_u(nur) = 1.D0-1.D-10
           jstart = nut+1
           jend = 0
c
           do j=1,nut
c select times in the range xp_time to endtime_pt
             if(j.lt.nut.and.ut(j+1).lt.xp_time)go to 9
             if(j.gt.1.and.ut(j-1).gt.endtime_pt)go to 9
               if(jstart.gt.j)jstart=j
               if(jend.lt.j)jend=j
c fix the boundaries to 0,1
c              if(ur(1).gt.0.025)then
c                par(1,j) = par(1,j) - 
c     >        (ur(1)-rho_u(1))*(par(2,j)-par(1,j))/(ur(2)-ur(1))
                par(nur,j) = par(nur,j) + 
     >             (ur(nur)-rho_u(nur))*(par(nur,j)-par(nur-1,j))/
     >             (ur(nur)-ur(nur-1))
c              endif
c interpolate onto ptor data grid
              call INTER_CSPL(nur,rho_u,par(1,j),nr,rho_d,pxp)
c save timeseries
              do k=1,nr
               u2d_t(k,j,i) = pxp(k)
              enddo
 9            continue
           enddo
c
           do j=jstart,jend
             jj=j
             if(j.eq.1)jj=j+1    
             if(ut(jj).ge.xp_time.and.ut(jj-1).lt.xp_time)then
c interpolate to the time of interest
               do k=1,nr
                 u2d(k,i) = u2d_t(k,jj-1,i)+
     >           (xp_time-ut(jj-1))*(u2d_t(k,jj,i)-u2d_t(k,jj-1,i))/
     >           (ut(jj)-ut(jj-1))
c                u2d(k,i) = u2d_t(k,jj,i)
               enddo
             endif
           enddo
         endif
        endif
        if(ierr.eq.1) then
           do k=1,mxgrid+1
              u2d(k,i)=0.D0
           enddo
        else if(ierr.eq.2) then
           if(iproc.eq.0) write(*,*) 'array dimensions 
     &                     incompatible, i = ',i
           stop
        endif
      enddo
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... Smooth data
c
      if(ismooth_all.ne.0)call average7_2d(u2d,nrmax,nr,n2d)
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... Read expdata files
c
c     loop over the 2d files that contain experimental data
c     do not interpolate to the 51 pts grid 
c
      do i=1,n3d
c
c     compose the file name:
         ig=linblnk(tok,10)
         is=linblnk(shot,10)
         i2d=linblnk(fields3d(i),7)
         ufile=tok(1:ig)//'2d'//shot(1:is)//'.'//fields3d(i)(1:i2d)
c
c Read the 3D ufile:
c
         ierr=0
         nx3(i)=0
         call u3read(cudir,ufile,88,u3dx(1,i),nx3(i),
     >               xp_time,u3dy(1,i),ierr)
c
c...Reset nx3(i) to zero if an error
c
         if(ierr.eq.1) then
            if (iproc.eq.0) write(*,57) ierr,ufile
            nx3(i)=0
            do j=1,mxgrid+1
                  u3dx(j,i)=0.D0
                  u3dy(j,i)=0.D0
            enddo
         else if(ierr.eq.2) then
            if(iproc.eq.0) write(*,*) 'array dimensions 
     &                      incompatible, i = ',i
            stop
         endif
      enddo
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... Loop over the 1d files next:
c
      do i=1,n1d
c
c     compose the file name:
         ig=linblnk(tok,10)
         is=linblnk(shot,10)
         i1d=linblnk(fields1d(i),7)
         ufile=tok(1:ig)//'1d'//shot(1:is)//'.'//fields1d(i)(1:i1d)
c
c Read the 1D ufile:
c
         ierr=0
         if(time_series.eq.0)then
c          call u1read(cudir,ufile,88,xp_time,u1d(i),ierr)
           call u1tread(cudir,ufile,88,nut,ut,u1d_t(1,i),nsmooth,ierr)
           if(ierr.eq.0)then
            do k=1,nut
             j=k
             if(xp_time.le.ut(1)) then
               xp_time=ut(1)+1.D-6
               write(*,35) ut(1)
             endif
c            if(i.eq.1) write(*,*) 'ut,xp_time = ',ut(j),xp_time
             if(k.eq.1)j=k+1
             if(ut(j).ge.xp_time.and.ut(j-1).lt.xp_time)then
c interpolate to the time of interest
                u1d(i) = u1d_t(j-1,i) + (u1d_t(j,i)-u1d_t(j-1,i))*
     >            (xp_time-ut(j-1))/(ut(j)-ut(j-1))
             endif
            enddo 
          endif
         else
           if(i.eq.1 .and. iproc.eq.0) write(6,65)
           do j=1,ntmax
             u1d_t(j,i) = 0.D0
           enddo
           call u1tread(cudir,ufile,88,nut,ut,u1d_t(1,i),nsmooth,ierr)
c
           if(ierr.eq.0)then
            do k=1,nut
             j=k
             if(k.eq.1)j=k+1
             if(ut(j).ge.xp_time.and.ut(j-1).lt.xp_time)then
c interpolate to the time of interest
c                 u1d(i) = u1d_t(j,i)
                u1d(i) = u1d_t(j-1,i) + (u1d_t(j,i)-u1d_t(j-1,i))*
     >            (xp_time-ut(j-1))/(ut(j)-ut(j-1))
             endif
            enddo 
          endif 
         endif
         if(ierr.eq.1) then
	    continue
            if (iproc.eq.0) write(6,57) ierr,ufile
            u1d(i)=0.D0
         else if(ierr.eq.2) then
            if(iproc.eq.0) write(*,*) 'array dimensions 
     &                      incompatible, i = ',i
            stop
         endif
      enddo
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c Finally, read the 0D data:
c
         ig=linblnk(tok,10)
         is=linblnk(shot,10)
c
         ufile=tok(1:ig)//'_'//shot(1:is)//'_0d.dat'

         id=linblnk(cudir,50)
c        ifchar=ig+is+2
         ifchar=ig+is+8
         if ( id .lt. 1 ) then
c
c..if no directory name is given, then just use the file name
c
            ichar  = ifchar            
            udfile = ufile
         else
c
c..concatenate the directory and file name into the character string cdfile
c  Note:
c  On a UNIX system, the directory name must end with '/'
c  On a VMS system, the directory name must end with ']'
c  If the directory name is given but does not end with '/' or ']'
c  then use the UNIX convention and add '/' to the directory name
c
        if (    cudir(id:id) .ne. '/'
     &    .and. cudir(id:id) .ne. ']' ) then
          ichar = id + 1 + ifchar
          udfile = cudir(1:id) // '/' // ufile(1:ifchar)
        else
          ichar  = id + ifchar
          udfile = cudir(1:id) // ufile(1:ifchar)
        endif
c
      endif
c
c..protect against names too long
c
      if ( ichar .gt. 130 ) then
       if (iproc .eq. 0 ) then
         write (*,*)
         write (*,*) 'directory and file names too long'
         write (*,*) ichar
     &        ,' = id + ifchar .gt. 130 in sbrtn read'
         write (*,*) ifchar,' = ifchar'
         write (*,*) id,' = idchar'
         write (*,*) cudir,' = dir'
         write (*,*) ufile,' = ufile'
       endif
      endif
c
c...  Read in 0d ufile data
c     u0read reads each element of a 0d file into a string
c     idatzero=0 old ITER PDB 0d file format
c     idatzero=1 for new 0D file formats
c     idatzero=3 for new 0D file format without DWDIA,DWMHD
c     Note that phase = OHM, L, LHLHL, H, HSELM,
c                       HGELM, HGELMH, VH, PEP
      ifail=0
      call u0read(idatzero,udfile,u0names,u0values,ifail)
c     open(unit=51,file=ufile(1:ichar), err=90)
      if (ifail.ne.0 .and. iproc .eq. 0) then
         write(*,*) 'Failure opening 0D ufile'
      else
         if(iproc.eq.0) then
           write(*,*) ' 0D file read successful, ifail=',ifail
         endif
      endif
      if(iproc.eq.0) write(*,*) ' idatzero = ',idatzero
c
c      do j=1,165
c        write(*,*) j, u0values(j)
c      enddo
c
      p_glob_exp=0.D0
      phase=' '
      if (ifail.eq.0) then
        if(idatzero.eq.0) then
          read(u0values(14),1010) apgasa
          read(u0values(15),1010) apgasz
          read(u0values(16),1010) abgasa
          read(u0values(17),1010) abgasz
          read(u0values(25),1010) apimpa
          read(u0values(26),1010) apimpz
          read(u0values(61),*) bt_0d
          read(u0values(62),1010) ip_0d
          read(u0values(64),1010) q95_0d
c Estimated thermal power loss (not corrected for cx and orbit losses)
          read(u0values(133),1010,err=1011) p_glob_exp
          read(u0values(135),1010,err=1011) gtauth_exp
          read(u0values(134),1010,err=1011) gtautot_exp
          p_glob_exp = p_glob_exp*1.D-6  ! [MW]
          u0phase = u0values(7)
c        write(*,*) 'q95_0d = ',q95_0d
        endif
        if(idatzero.eq.1) then
          read(u0values(94),*) apgasa
          read(u0values(95),*) apgasz
          read(u0values(99),1010) apimpa
          read(u0values(100),1010) apimpz
          read(u0values(8),*) abgasa
          read(u0values(10),*) abgasz
          u0phase = u0values(96)
          read(u0values(115),*) rgeo_0d
          read(u0values(1),*) rmin_0d
          read(u0values(15),*) bt_0d
          read(u0values(51),*) ip_0d
          read(u0values(56),*) kap_0d
          read(u0values(19),*) delt_0d
          read(u0values(64),*) nebar_0d
          read(u0values(106),*) pnbi_0d
          read(u0values(98),*) pich_0d
          read(u0values(91),*) pech_0d
          read(u0values(139),*) ti0_0d
          read(u0values(134),*) te0_0d
          read(u0values(111),*) q0_0d
          read(u0values(110),*) q95_0d
          read(u0values(131),*) gtauth_exp
          read(u0values(133),*) gtautot_exp
          read(u0values(159),*) wth_0d
          read(u0values(160),*) wtot_0d
        endif
        if(idatzero.eq.2) then
          read(u0values(37),*) apgasa
          read(u0values(38),*) apgasz
          read(u0values(73),1010) apimpa
          read(u0values(74),1010) apimpz
          read(u0values(39),*) abgasa
          read(u0values(40),*) abgasz
          u0phase = u0values(8)
          read(u0values(25),*) rgeo_0d
          read(u0values(27),*) rmin_0d
          read(u0values(10),*) bt_0d
          read(u0values(11),*) ip_0d
          read(u0values(32),*) kap_0d
          read(u0values(35),*) delt_0d
          read(u0values(64),*) nebar_0d
          read(u0values(43),*) pnbi_0d
          read(u0values(84),*) pich_0d
          read(u0values(91),*) pech_0d
          read(u0values(56),*) ti0_0d
          read(u0values(53),*) te0_0d
          read(u0values(59),*) q0_0d
          read(u0values(60),*) q95_0d
          read(u0values(20),*) gtauth_exp
          read(u0values(21),*) gtautot_exp
          read(u0values(14),*) wth_0d
          read(u0values(15),*) wtot_0d
        endif
        if(idatzero.eq.3) then
          read(u0values(92),*) apgasa
          read(u0values(93),*) apgasz
          read(u0values(97),1010) apimpa
          read(u0values(98),1010) apimpz
          read(u0values(8),*) abgasa
          read(u0values(10),*) abgasz
          u0phase = u0values(94)
          read(u0values(113),*) rgeo_0d
          read(u0values(1),*) rmin_0d
          read(u0values(15),*) bt_0d
          read(u0values(49),*) ip_0d
          read(u0values(54),*) kap_0d
          read(u0values(19),*) delt_0d
          read(u0values(62),*) nebar_0d
          read(u0values(104),*) pnbi_0d
          read(u0values(96),*) pich_0d
          read(u0values(89),*) pech_0d
          read(u0values(137),*) ti0_0d
          read(u0values(132),*) te0_0d
          read(u0values(109),*) q0_0d
          read(u0values(108),*) q95_0d
          read(u0values(129),*) gtauth_exp
          read(u0values(131),*) gtautot_exp
          read(u0values(157),*) wth_0d
          read(u0values(158),*) wtot_0d
        endif
        ip_0d=ip_0d/1.D6
        nebar_0d=nebar_0d/1.D19
        pnbi_0d=pnbi_0d/1.D6
        pich_0d=pich_0d/1.D6
        pech_0d=pech_0d/1.D6
        ti0_0d=ti0_0d/1.D3
        te0_0d=te0_0d/1.D3
        wth_0d=wth_0d/1.D6
        wtot_0d=wtot_0d/1.D6
c
        if (u0phase(1:1).eq.'o'.or.u0phase(1:1).eq.'O') phase='O'
        if (u0phase(1:1).eq.'l'.or.u0phase(1:1).eq.'L') phase='L'
        if (u0phase(1:1).eq.'h'.or.u0phase(1:1).eq.'H') phase='Hnoelm'
        if (u0phase(2:2).eq.'s'.or.u0phase(2:2).eq.'S') phase='Helm'
        if (u0phase(2:2).eq.'g'.or.u0phase(2:2).eq.'G') phase='Helm'
      endif
c
c     Assume the working and beam ions are the same and the impurity
c     is set by pimpa,pimpz. In xpin, the defaults are the following:
c       amassgas_exp=2.0     amassimp_exp=12.0
c       zgas_exp=1.0         zimp_exp=6.0         
c       pimpa=12.0           pimpz=6.0
c     All five variables may also be set in the input file 
c     if the 0D file information is not present or desired.
c
      if (apgasa.lt.1.D-3) apgasa=amassgas_exp
      if (apgasz.lt.1.D-3) apgasz=zgas_exp
      if (abgasa.lt.1.D-2) abgasa=amassgas_exp
      if (abgasz.lt.1.D-2) abgasz=1.D0
      if (apimpa.lt.1.D-2) apimpa=pimpa
      if (apimpz.lt.1.D-2) apimpz=pimpz
c      if (apimpa.lt.1.D-2) apimpa=amassimp_exp
c      if (apimpz.lt.1.D-2) apimpz=zimp_exp
c
      if (iproc.eq.0) then
        write(*,11)
        if(idatzero.ge.1)then
          write(*,13) rgeo_0d,rmin_0d,bt_0d,ip_0d,nebar_0d,pnbi_0d,
     &                pich_0D,pech_0d,ti0_0d,te0_0d,q0_0d,q95_0d,
     &                kap_0d,delt_0d,wth_0d,wtot_0d
          write(*,12) apgasa,apgasz,abgasa,abgasz,apimpa,apimpz
        else
          write(*,12) apgasa,apgasz,abgasa,abgasz,apimpa,apimpz
        endif
      endif
c
c...PLTH second best option:
c...  Read total input power corrected for cx and orbit losses
c...  Note that here p_glob_exp can still be .le.0
c...  In that case it is assigned a value (according
c...  to 1D data) just before ITER89P is calculated.
c
      if (p_glob_exp.le.0.0 .and. idatzero.eq.0) then
        if (ifail.eq.0) then
          read(u0values(76),1010,err=1012) p_glob_exp
          p_glob_exp = p_glob_exp*1.D-6                  ![MW]
        endif
      endif
c
 1010 format(1pe10.3)
 1011 continue
 1012 continue
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c nj : the size of the vectors printed in this file
          nj_d=nr
c
c nion : the number of ion species          
          nion_d=2
c
c nprim : the number of primary ion species
          nprim_d=1
c
c nimp : the number of impurity ion species          
          nimp_d=1
c
c nneu : the number of neutral ion species
          nneu_d=1
c
c ibion : index of beam species          
          ibion_d=1
c
c namep : name(s) of primary ion species          
c         (namep_d(i),i=1,nprim_d)
c
c namei : name(s) of impurity ion species          
c         (namei_d(i),i=1,nimp_d)
c
c namen : name(s) of neutral ion species         
c         (namen_d(i),i=1,nneu_d)
c              
c time :  time at which data is printed          
          time_d=xp_time
c
c Rgeom : major radius of geometric          
         rgeom_d=u1d(9)                  !meters
c 
c Rmag :  major radius of mag axis, meters        
         rmag_d=u2d(1,15)                    !meters
         if (rmag_d.eq.0.) rmag_d=u1d(9)
c
c R0 : major radius of vaccuum btor ref            
         rmajor_d=u1d(9)                !meters
         amin_d=u1d(1)                  !meters
c
c JAK 960325 warning output added
         if (rmajor_d.eq.0.0)
     >     write(6,*) 'Warning rep_iter: rmajor zero'
c
c kappa : plasma elongation         
         if(u1d(6).gt.0.) then
           kappa_d=u1d(6)
         else
           u1d(6)=1.D0
           kappa_d=1.D0
         endif
c delta : plasma triangularity         
         deltao_d=u1d(3)
c  
c pindent : plasma indentation         
         pindento_d=u1d(4)
c
c volo : plasma volume,meters**3         
         volo_d=u2d(nr,21)                    !meters**3
c
c Btor : vaccuum toroidal field at rmajor, tesla         
         btor_d=u1d(2)                    !tesla
         if(btscale.gt.1.) btor_d=btor_d*btscale
c
c total, ohmic, bootstrap, beam, and rf currents, amps       
c        tocur_d
c        totohm_d
c        totboot_d
c        totbeam_d
c        totrf_d
c
         tocur_d=u1d(5)
         curtot=tocur_d/1.D6
c
c betap : poloidal beta         
c         betap_d
c
c beta : toroidal beta         
c         beta_d
c
c ali : plasma inductance         
         ali_d=u1d(7)
c
c te0 : central electron temperature         
         te0_d=u2d(1,23)/1.D3                        !kev
c
c ti0 : central ion temperature         
         ti0_d=u2d(1,22)/1.D3                        !kev
c
c pohm: total ohmic input power                      !Wa  (JAK)
         pohm_d=u1d(12)
c
c vsurf: loop voltage at the plasma boundary         !V   (JAK)
         vsurf_d=u1d(13)
c
c---psi on rho grid,volt*sec/rad 
c              
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(psir_d(j),j=1,nj_d)      !volt*se/rad        
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c   
c---rho grid, meters
c
c get rho_a from edge kappa
c r_d = rho*a*sqrt(kappa)
c
          do j=1,nj_d
             r_d(j)=rho_d(j)*u1d(1)*sqrt(u1d(6))    !meters
             r_d(j)=r_d(j)*rho_norm_loc
          enddo
          if (r_d(nj_d).eq.0.0 .and. iproc.eq.0) then
            write(6,*) 'Warning rep_iter: r_d zero'
          endif
c        stop
c
c---flux surface average absolute grad rho and
c   flux surface average ( grad rho)**2
c
         do j=1,nj_d
           grho1npsi_d(j)=u2d(j,5)*r_d(nj_d)
           grho2npsi_d(j)=u2d(j,6)*r_d(nj_d)**2
         enddo
c
c---fcap, (ie f(psilim)/f(psi) )
c   place holder for fcap until can get from U-file  GMS 1-12-2007
         do j=1,nj_d
           fcap_d(j) = 1.D0
         enddo
c
c---gcap, (ie <(grad rho)**2*(R0/R)**2> )
c   place holder for gcap until can get from U-file  GMS 1-12-2007
         do j=1,nj_d
           gcap_d(j) = grho2npsi_d(j)/(1.D0-(r_d(j)/rmajor_d)**2)**1.5
         enddo
c       
c---hcap, (surface area/(4*pi*pi*R0*gradrho*rho))
c         (vprime/(4*pi*pi*R0*gradrho*rho*<abs(gradrho)>))
c         with surface area = vprime*<abs(gradrho)> and
c         vprime = dV/drho
c
           rfix=1.D0
           if(istk.eq.88) rfix=4.46
           sfix=1.D0
           if(istk.eq.88) sfix=22.825
           if(istk.eq.88) write(6,*) 'sfix=22.825'
c
          do j=2,nj_d
            hcap_d(j)=sfix*u2d(j,19)/(u2d(j,5)*
     >      r_d(nj_d)/rfix)/(4.D0*pi_m*pi_m*rmajor_d*r_d(j))
c           write(*,100) j, r_d(j), u2d(j,19), hcap_d(j)
          enddo
          hcap_d(1)=hcap_d(2)
c
c---volume of each flux surface,  meters**3
c
        do j=1,nj_d
          psivolp_d(j)=u2d(j,21)
        enddo
c                         
c--- flux surface area,  meters**2 4.*pi*pi*R0*hcap*rho*<abs(grad rho)>
c
c        u2d(2,19)=0.D0
        if ( u2d(2,19).eq.0 ) then
          if(iproc.eq.0)
     >     write(6,'(a33)') 'Computing flux surface area, hcap'
          do j=2,nj_d-1
           sfareanpsi_d(j)=
     >       (psivolp_d(j+1)-psivolp_d(j-1))/2.D0 /
     >       (rho_d(j+1)-rho_d(j)) / r_d(nj_d) * grho1npsi_d(j)
           hcap_d(j)=
     >       (psivolp_d(j+1)-psivolp_d(j-1)) /2.D0 /
     >       (rho_d(j+1)-rho_d(j)) / r_d(nj_d) /
     >       (4.D0*pi_m*pi_m*rmajor_d*r_d(j))
          enddo
          sfareanpsi_d(1)=1.D-6
          sfareanpsi_d(nj_d)=
     >       (psivolp_d(nj_d)-psivolp_d(nj_d-1)) /
     >       (rho_d(nj_d)-rho_d(nj_d-1)) / r_d(nj_d) * 
     >       grho1npsi_d(nj_d)
          hcap_d(5)=hcap_d(6)
          hcap_d(4)=hcap_d(5)
          hcap_d(3)=hcap_d(4)
          hcap_d(2)=hcap_d(3)
          hcap_d(1)=hcap_d(2)
        else
          do j=1,nj_d
            sfareanpsi_d(j)=sfix*u2d(j,19)
          enddo
        endif
c        write(*,*) 'sfareanpsi = ',sfareanpsi_d(nj_d)
c        write(*,*) 'area = ',3.14*amin_d**2*kappa_d
c        write(*,*) 'kappa = ',kappa_d
c
c       do j=2,nj_d-1
c         write(*,54) j, r_d(j), sfareanpsi_d(j), 
c    >       (psivolp_d(j+1)-psivolp_d(j-1)) /2.D0 /
c    >       (rho_d(j+1)-rho_d(j)) / r_d(nj_d) * grho1npsi_d(j),
c    >        hcap_d(j),
c    >       (psivolp_d(j+1)-psivolp_d(j-1)) /2.D0 /
c    >       (rho_d(j+1)-rho_d(j)) / r_d(nj_d) /
c    >       (4.D0*pi_m*pi_m*rmajor_d*r_d(j))
c       enddo
c
c---read original experimental data for ne, te and ti ---JAK
c
c ne and errorbar
c
          nnexp_d = nx3(1)
          do j=1,nnexp_d
           xnexp_d(j)=u3dx(j,1)
           ynexp_d(j)=u3dy(j,1)*1.D-19
           xnexpeb_d(j)=u3dx(j,2)
           ynexpeb_d(j)=u3dy(j,2)*1.D-19
          enddo
c
c te and errorbar
c
          ntexp_d = nx3(3)
          do j=1,ntexp_d
           xtexp_d(j)=u3dx(j,3)
           ytexp_d(j)=u3dy(j,3)*1.D-3        !keV
           xtexpeb_d(j)=u3dx(j,4)
           ytexpeb_d(j)=u3dy(j,4)*1.D-3      !keV
          enddo
c
c ti and errorbar
c
          ntixp_d = nx3(5)
          do j=1,ntixp_d
           xtixp_d(j)=u3dx(j,5)
           ytixp_d(j)=u3dy(j,5)*1.D-3        !keV
           xtixpeb_d(j)=u3dx(j,6)
           ytixpeb_d(j)=u3dy(j,6)*1.D-3      !keV
          enddo
c
c---electron temperature, kev
c
          do j=1,nj_d
           te_d(j)=u2d(j,23)/1.D3            !kev
          enddo
c
c---ion temperature, kev
c
          do j=1,nj_d
           ti_d(j)=u2d(j,22)/1.D3          !kev
          enddo
  
       if( iexp_exch.eq.2) then
        do j=1,nj_d
         ti_d(j)=te_d(j)-.001
        enddo
       endif    
c
c---q (ie safety factor) profile
c
          do j=1,nj_d
           q_d(j)=u2d(j,27)
c           write(*,*) j, r_d(j), q_d(j), 'q_d'
          enddo
          q01_exp=q_d(1)
          qa1_exp=q_d(nj_d)
c
c---electron density,#/m**3
c
          do j=1,nj_d
           ene_d(j)=u2d(j,7)           !#/meter**3
c           write(*,54) j, rho_d(j), ene_d(j)
          enddo                         
c      
c---primary ion density,#/m**3,species (NM1,NM2,NM3)
c---impurity ion density,#/m**3,species (NIMP)
c   assuming abgasz=apgasz
             do j=1,nj_d
               en_d(j,1)=u2d(j,10) + u2d(j,42) + u2d(j,43) !#/meter**3
               en_nm1_d(j)=u2d(j,10)
               en_nm2_d(j)=u2d(j,42)
               en_nm3_d(j)=u2d(j,43)
c              if (imodel.lt.0) then
c                en_d(j,1)=u2d(j,10) + u2d(j,42) + 
c    >                     u2d(j,43) + u2d(j,9) ! chi=1 test
c              endif
               en_d(j,2)=u2d(j,9)    !#/meter**3
c               write(*,100) j, rho_d(j),u2d(j,20)
               zeff_t=u2d(j,20)
               if (zeff_t.eq.0.) zeff_t=u1d(11)
               if (zeff_t.eq.0.) then
                 zeff_t=1.
                 if(iproc.eq.0) write(6,*) 'zeff_t=1.'
               endif
c
               if(iproc.eq.0 .and. apimpz.lt.2.D0)
     >           write(6,*) '*** WARNING: apimpz looks fishy ***'
               if(en_d(j,1).eq.0.) 
     >          en_d(j,1)=(apimpz-zeff_t)/(apimpz-apgasz)/apgasz
     >                    *u2d(j,7)-u2d(j,8)
               if(en_d(j,2).eq.0.) 
     >          en_d(j,2)=(zeff_t-apgasz)/(apimpz-apgasz)/apimpz
     >		              *u2d(j,7)
               if(en_d(j,1).lt.0..or.en_d(j,2).lt.0.) then
                 if(iproc.eq.0) write(6,*) 'negative density fixed'
                 en_d(j,1)=dabs(en_d(j,1))
                 en_d(j,2)=dabs(en_d(j,2))
               endif
c               write(*,100) j, rho_d(j), ene_d(j), en_d(j,1), en_d(j,2)
             enddo
c
c---sion   : source due to ionization,        #/(m**3*sec),species
c---srecom : source due to recombination,     #/(m**3*sec),species
c---scx    : source due to cx thermal neut.,  #/(m**3*sec),species
c---sbcx   : sink due to cx with beam neut.   #/(m**3*sec),species
c---s      : total source rate,               #/(m**3*sec),species
c---dudt   : s dot,                           #/(m**3*sec),species
c
       do jj=1,nprim_d       
c               read(niterdb,'(a)')stflg
c               read(niterdb,10)(sion_d(j,jj),j=1,nj_d)    !#/meter**3/sec
c               
c               read(niterdb,'(a)')stflg
c              read(niterdb,10)(srecom_d(j,jj),j=1,nj_d)  !#/meter**3/sec
c               
c               read(niterdb,'(a)')stflg
c               read(niterdb,10)(scx_d(j,jj),j=1,nj_d)     !#/meter**3/sec
c               
c               read(niterdb,'(a)')stflg
c               read(niterdb,10)(sbcx_d(j,jj),j=1,nj_d)    !#/meter**3/sec
c               
c               read(niterdb,'(a)')stflg
c               read(niterdb,10)(s_d(j,jj),j=1,nj_d)       !#/meter**3/sec
c               
c               read(niterdb,'(a)')stflg
c               read(niterdb,10)(dudtsv_d(j,jj),j=1,nj_d)  !#/meter**3/sec
c
            do j=1,nj_d
              sion_d(j,jj)=0.
              srecom_d(j,jj)=0.
              scx_d(j,jj)=0.
              sbcx_d(j,jj)=0.
              s_d(j,jj)=0.
              dudtsv_d(j,jj)=0.
             enddo
          enddo
c
c---fast ion density, #/m**3, species
c
          do j=1,nj_d
             enbeam_d(j)=u2d(j,8)         !#/meter**3
          enddo
c      
c---neutral density,  #/m**3,species
c
       do jn=1,nneu_d
           do j=1,nj_d
              enn_d(j,jn)=0.     !#/meter**3
          enddo
       enddo
c
c---neutral density from  wall source, #/m**3, species
c
       do jn=1,nneu_d        
c               read(niterdb,'(a)')stflg
c               read(niterdb,10)(ennw_d(j,jn),j=1,nj_d)     !#/meter**3     
         do j=1,nj_d
          ennw_d(j,jn)=0.
         enddo
       enddo
c
c---neutral density from volume source, #/m**3,species (recomb and beam cx)
c
       do jn=1,nneu_d
c               read(niterdb,'(a)')stflg
c               read(niterdb,10)(ennv_d(j,jn),j=1,nj_d)     !#/meter**3
        do j=1,nj_d
         ennv_d(j,jn)=0.
        enddo
       enddo
c
c---volume source of neutrals, #/(m**3*sec),species
c
c       do jn=1,nneu_d   
c               read(niterdb,'(a)')stflg
cx               read(niterdb,10)(volsn_d(j,jn),j=1,nj_d)!#/(m**3*sec)
c               read(niterdb,10)(blank_d(j),j=1,nj_d) 
c       enddo
c
c---sbion : beam electron source,  #/(m**3*sec
c        
c               read(niterdb,'(a)')stflg
cx               read(niterdb,10)(sbion_d(j),j=1,nj_d)    !#/(m**3*sec)
c               read(niterdb,10)(blank_d(j),j=1,nj_d)          
c
c---sbeam : beam thermal ion source,  #/(m**3*sec)
c
c               read(niterdb,'(a)')stflg   
c               read(niterdb,10)(sbeam_d(j),j=1,nj_d)      !#/(m**3*sec)
c
c...Note that SNBII as read contains all sources and sinks due to
c...beam thermal ion source and charge exchange processes; 
c...so: field SNBII = sbeam + scx + sbcx + srecom    (note sbcx<0 sink)
c...Read the source in sbeam_d and let all other source terms be 0.          
c
             do j=1,nj_d
              sbeam_d(j)=u2d(j,18)                       !#/(m**3*sec)
             enddo
c
c...thermal ion particle source due to ionisation; field SWALL
c
             do j=1,nj_d
              sion_d(j,1)=u2d(j,30)                       !#/(m**3*sec)
             enddo
             do j=1,nj_d
               dudtsv_d(j,1)=u2d(j,41)                    ! dn/dt
             enddo
c
c---total current density, amps/m**2
c   (avoids negative currenty density)
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(curden_d(j),j=1,nj_d)        !amps/meter**2
             do j=1,nj_d
              curden_d(j)=DABS(u2d(j,4))
             enddo
c    
c---ohmic current density, amps/m**2
c     
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(curohm_d(j),j=1,nj_d)        !amps/meter**2
c          read(niterdb,10)(blank_d(j),j=1,nj_d) 
c          
c---bootstrap current density, amps/m**2
c      
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(curboot_d(j),j=1,nj_d)       !amps/meter**2
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c      
c--- beam driven current density, amps/m**2
c     
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(curdbeam_d(j),j=1,nj_d)      !amps/meter**2
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c               
c---rf current density, amps/m**2
c    
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(currf_d(j),j=1,nj_d)         !amps/meter**2
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c          
c---rho*bp0*fcap*gcap*hcap, tesla*meters
c     
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(rbp_d(j),j=1,nj_d)
c          read(niterdb,10)(blank_d(j),j=1,nj_d)         !tesla*meters      
c
c---zeff profile:
c
          do j=1,nj_d
           zeff_d(j)=u2d(j,20)
           if (zeff_d(j).eq.0.) zeff_d(j)=u1d(11)
c          write(*,54) j, rho_d(j), zeff_d(j)
          enddo
c
c---angular rotation speed profile, rad/sec
c
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(angrot_d(j),j=1,nj_d)       !rad/sec
          do j=1,nj_d
            angrot_d(j)=u2d(j,29)
          enddo
c
c---torque density, N/M**2
c Note: TRANSP ufiles have N-m/cm**3
c
          if(itorque.ne.0)then
            if(tok.eq.'jet') conv_torq=1.D0
            do j=1,nj_d
              torque_d(j)=u2d(j,46)*conv_torq
            enddo
          endif
c
c---electron thermal diffusivity, meters*2*/sec on half grid
c
          do j=1,nj_d
           chieinv_d(j)=u2d(j,1)       !meters**2/sec
          enddo
c                 
c---ion thermal diffusivity,meters*2*/sec  on half grid
c     
          do j=1,nj_d
           chiinv_d(j)=u2d(j,2)       !meters**2/sec
          enddo
c         
c---ion neocl.  thermal conductivity, 1/(m*sec) on half grid 
c        
          do j=1,nj_d
           xkineo_d(j)=0.D0      !meters**2/sec
          enddo
c
c---wdot,electrons, watts/m**3 d(electron energy)/dt profile (DWER):
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(dpedtc_d(j),j=1,nj_d)        !watts/meter**3      
          do j=1,nj_d
           dpedtc_d(j)=u2d(j,31)
          enddo
c    
c---wdot,ions,watts/m**3 d(ion energy)/dt profile (DWIR):
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(dpidtc_d(j),j=1,nj_d)        !watts/meter**3
          do j=1,nj_d
           dpidtc_d(j)=u2d(j,32)
          enddo
c
c---electron conduction, watts/m**3
cc      
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(qconde_d(j),j=1,nj_d)        !watts/meter**3 
          do j=1,nj_d
           qconde_d(j)=0.D0
          enddo
c   
c---ion conduction, watts/m**3
c
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(qcondi_d(j),j=1,nj_d)        !watts/meter**3
          do j=1,nj_d
           qcondi_d(j)=0.D0
          enddo
c
c---electron convection,watts/m**3 
c     
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(qconve_d(j),j=1,nj_d)        !watts/meter**3
          do j=1,nj_d
           qconve_d(j)=0.D0
          enddo
c
c---ion convection, watts/m**3 
c     
c         read(niterdb,'(a)')stflg
c         read(niterdb,10)(qconvi_d(j),j=1,nj_d)        !watts/meter**3
          do j=1,nj_d
           qconvi_d(j)=0.D0
          enddo
c
c---power to elec.from beam, watts/m**3 
c     
           do j=1,nj_d
             qbeame_d(j)=u2d(j,11)        !watts/meter**3
           enddo
c          
c---qdelt : electron-ion equilibration, watts/m**3
c
         do j=1,nj_d
          qdelt_d(j)=u2d(j,28)
         enddo
c
c NRL formulary :
c
       if(qdelt_d(nj_d).eq.0.or.iexp_exch.ge.1) then       
         do j=1,nj_d
           zbrac=en_d(j,1)/ene_d(j)+apgasa/
     >         apimpa*en_d(j,2)/ene_d(j)*apimpz**2
           alamda=24.D0-dlog((dabs(ene_d(j))*1.D-6)**zhalf/
     >         (dabs(te_d(j))*1.D3))
           tau_e=3.44D5*(dabs(te_d(j))*1.D3)**1.5D0/
     >         (ene_d(j)*1.D-6)/alamda
           qdelt_d(j)=-10**6*1.6022D-19*zbrac*3.D0/(1836.D0*apgasa)/
     >         tau_e*ene_d(j)*1.D-6*(te_d(j)-ti_d(j))*1.D3
         enddo                    
       endif     !watts/meter**3
       if(iexp_exch.eq.-1) then
        do j=1,nj_d
         qdelt_d(j)=0.D0
        enddo
       endif
c
c---power to ions from beam, watts/m**3
c   
           do j=1,nj_d
             qbeami_d(j)=u2d(j,12)        !watts/meter**3
           enddo
c    
c---qrfe, rf electron heating, watts/m**3 (QECHE + QICRHE)
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(qrfe_d(j),j=1,nj_d)          !watts/meter**3 
           do j=1,nj_d
            qrfe_d(j)=u2d(j,33)*echconv + u2d(j,34)
           enddo
c    
c---qrfi, rf ion heating, watts/m**3 (QECHI + QICRHI)
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(qrfi_d(j),j=1,nj_d)          !watts/meter**3 
           do j=1,nj_d
            qrfi_d(j)=u2d(j,35) + u2d(j,36)
           enddo
c
c... calculate integrated RF heating
c
       prfesum=0.D0
       prfisum=0.D0
       do j=2,nj_d
         drm=r_d(j)-r_d(j-1)
         dvoldr_p=2.D0*pi_m*r_d(j)*2.D0*pi_m*rmajor_d*hcap_d(j)
         dvoldr_m=2.D0*pi_m*r_d(j-1)*2.D0*pi_m*rmajor_d*hcap_d(j-1)
         ds=2.D0*pi_m*(r_d(j)*hcap_d(j)+r_d(j-1)*
     >      hcap_d(j-1))/2.D0*drm
         prfesum=prfesum+
     >     zhalf*(dvoldr_p*qrfe_d(j)+dvoldr_m*qrfe_d(j-1))*drm
         prfisum=prfisum+
     >     zhalf*(dvoldr_p*qrfi_d(j)+dvoldr_m*qrfi_d(j-1))*drm
       enddo
c      if(iproc.eq.0) then
c        write(6,'(a27,2F10.3)') 'Total integrated Prfe,Prfi [MW]:',
c    >        prfesum*1.D-6,prfisum*1.D-6
c      endif
c
c
c---qlhe, lower hybrid heating to electrons, watts/m**3
c   
        do j=1,nj_d
          qlhe_d(j)=u2d(j,45)
        enddo
c
c---qione, recombination and impact ionization, watts/m**3
c   use QWALLE = ''thermal electron heat loss due
c   to the ionization of wall neutrals W/m3''
c
        do j=1,nj_d
          qione_d(j)=u2d(j,37)
        enddo
c      
c---qioni, recombination and impact ionization, watts/m**3
c   use QWALLI = ''main thermal ion heat loss due to
c   ionization and charge exchange with wall neutrals in W/m3''
c   SO: QWALLI=qioni_d + qcx_d
c   
        do j=1,nj_d
          qioni_d(j)=u2d(j,38)
        enddo
c      
c---qcx,  neutral-ion charge exchange, watts/m**3
c          read(niterdb,'(a)')stflg
c          read(niterdb,10)(qcx_d(j),j=1,nj_d)           !watts/meter**3
        do j=1,nj_d
          qcx_d(j)=0.D0
        enddo
c     
c---2d MHD equil. electron heating, watts/m**3
c      
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qe2d_d(j),j=1,nj_d)        !watts/meter**3
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c         
c---2d MHD equil.ion heating, watts/m**3
c      
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qi2d_d(j),j=1,nj_d)         !watts/meter**3
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c         
c---fusion electron heating, watts/m**3
c
        do j=1,nj_d
          qfuse_d(j)=u2d(j,39)
        enddo
c               
c---fusion ion heating, watts/m**3
c
        do j=1,nj_d
          qfusi_d(j)=u2d(j,40)
        enddo    
c         
c---beam fusion electron heating,watts/m**3
c   
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qbfuse_d(j),j=1,nj_d)       !watts/meter**3
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c                
c---beam fusion ion heating profile
c      
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qbfusi_d(j),j=1,nj_d)       !watts/meter**3
c         read(niterdb,10)(blank_d(j),j=1,nj_d)
c                 
c---qmag electron heating, watts/m**3
c      
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qmag_d(j),j=1,nj_d)         !watts/meter**3
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c               
c---sawtooth electron heating, watts/m**3
c      
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qsawe_d(j),j=1,nj_d)        !watts/meter**3
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c               
c---sawtooth ion  heating, watts/m**3
c       
c          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(qsawi_d(j),j=1,nj_d)        !watts/meter**3
c          read(niterdb,10)(blank_d(j),j=1,nj_d)
c         
c---radiated power density, watts/m**3
c       
        do j=1,nj_d
          qrad_d(j)=u2d(j,14)
        enddo
c
c---(electron) ohmic power density,watts/m**3
c    
        do j=1,nj_d
          qohm_d(j)=u2d(j,13)
        enddo
c
c--- Total plasma pressure, Kev/m**3
c    Pa = N/m**2, keV/m**3=Pa/1.602e-16
c    PTOTR has units of Pascals
c    Use thermal pressure if total pressure not provided
c
       if(u2d(1,44) .eq. 0.0) then
         if(iproc.eq.0) then
          write(6,'(a31)') 'Error: total pressure not found'
          write(6,'(a32)') 'Computing total thermal pressure'
         endif
         do j=1,nj_d
            ptot_d(j) = ene_d(j)*(te_d(j) + ti_d(j))
            pfast_d(j) = 0.0
c            write(*,53) j, r_d(j)/r_d(nj_d),ptot_d(j)
         enddo
       else
         do j=1,nj_d
            ptot_d(j)=u2d(j,44)/1.602D-16     !Kev/m**3
            pfast_d(j)=ptot_d(j) - (ene_d(j)*te_d(j) +
     >            (ene_d(j) - (zimp_exp-1.D0)*
     >            en_d(j,nprim_d+1))*ti_d(j))
            if(ipfst.eq.1) pfast_d(j)=0.D0
c           write(*,53) j, r_d(j)/r_d(nj_d),ptot_d(j),pfast_d(j)
         enddo
       endif
c
c... calculate integrated radiation loss
c
       pradsum=0.D0
       do j=2,nj_d
         drm=r_d(j)-r_d(j-1)
         dvoldr_p=2.D0*pi_m*r_d(j)*
     >            2.D0*pi_m*rmajor_d*hcap_d(j)
         dvoldr_m=2.D0*pi_m*r_d(j-1)*
     >            2.D0*pi_m*rmajor_d*hcap_d(j-1)
         ds=2.D0*pi_m*(r_d(j)*hcap_d(j)+r_d(j-1)*
     >      hcap_d(j-1))/2.D0*drm
         pradsum=pradsum+
     >     zhalf*(dvoldr_p*qrad_d(j)+dvoldr_m*qrad_d(j-1))*drm
       enddo
c      if(iproc.eq.0) write(6,'(a27,F10.3)') 
c    >   'Total integrated Prad [MW]:', pradsum*1.D-6
c
c... calculate integrated ohmic heating power
c
       if (pohm_d.eq.0.) then
         pohsum=0.D0
         do j=2,nj_d
           drm=r_d(j)-r_d(j-1)
           dvoldr_p=2.D0*pi_m*r_d(j)*
     >              2.D0*pi_m*rmajor_d*hcap_d(j)
           dvoldr_m=2.D0*pi_m*r_d(j-1)*
     >              2.D0*pi_m*rmajor_d*hcap_d(j-1)
           pohsum=pohsum+
     >      zhalf*(dvoldr_p*qohm_d(j)+dvoldr_m*qohm_d(j-1))*drm
         enddo
       if(iproc.eq.0) write(6,'(a27,F10.3)') 
     >   'Total integrated Pohm [MW]:', pohsum*1.D-6
       else
       if(iproc.eq.0) write(6,'(a30,F7.3)') 
     >   'Total Pohm [MW] from 1D ufile:',pohm_d*1.D-6
       endif
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c...Calculate q-profile if not available
c   using qa1_exp and q01_exp values
c   alfj_exp=(qa1_exp/q01_exp-1.) corresponds to ja=0 and shata_exp=2.
c
       qa1_exp=u1d(8)
       q01_exp=q0_exp
       alfj1_exp=(qa1_exp/q01_exp-1.D0)
      if(q_d(nj_d).eq.0.or.iexp_q.eq.-1) then
        if(iproc.eq.0) write(6,*) 'Calculating q(r) ...'
c
      do j=2,nj_d
       xxq=rho_xp(j)**2
       q_d(j)=1.D0/(1.D0/(q01_exp*xxq*(alfj1_exp+1.D0))*
     >               (1.D0-(1.D0+1.D-10-xxq)**(alfj1_exp+1)))
      enddo
c
      if(qmin_exp.ne.0.) then
      do j=2,nj_d
       xxq=rho_xp(j)**2
       q_d(j)=(rho_xp(j)-rho_qm_exp)**2/(1.D0-rho_qm_exp)**2
     >   *(qa_exp-qmin_exp)+qmin_exp
      enddo 
      endif
         q_d(1)=q_d(2)
      endif
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c...Calculate J(r) from q-profile if not available
c
           if (curden_d(1).EQ.0) then 
c
             if(iproc.eq.0) 
     &          write(6,*) 'Calculating j(r) from q-profile'
c
c...         Check if total current is given
c
c...         Normalize to total current; if not available then use
c...         d3D formula to relate Ip (in MegAmps) to q_95, etc., (obtained from
c...         Jim Deboo in Feb. 1995., who says it usually agrees with efit
c...         to within 5-10%):
c...         1.e6 convert totohm_d from MA to A,
c
             if (tocur_d.EQ.0) then
               if(iproc.eq.0) write(6,*) 'Total 
     &            current calculated from q(r)'
               tocur_d=(5.D0*r_d(nj_d)**2*btor_d/rmajor_d/q_d(nj_d))
     &         *(1.D0+1.5D0*(r_d(nj_d)/rmajor_d)**2)
     &         *(1.D0+kappa_d**2)/2.D0*1.D6
             endif
c
c...         Assume j(r)=j(0)*(1-rhohat**2)**(alpha)
c...         Note that by approximation (circular geometry)
c...         alpha = qa/qa0-1
c
             alpha = qa1_exp/q01_exp-1.D0
             sum=0.D0
             do j=2,nj_d
               drm=r_d(j)-r_d(j-1)
               dvoldr_p=2.D0*pi_m*r_d(j)*
     >                  2.D0*pi_m*rmajor_d*hcap_d(j)
               dvoldr_m=2.D0*pi_m*r_d(j-1)*
     >                  2.D0*pi_m*rmajor_d*hcap_d(j-1)
               ds=2.D0*pi_m*(r_d(j)*hcap_d(j)+
     >            r_d(j-1)*hcap_d(j-1))/2.D0*drm
               curden_d(j) = (1.D0-(r_d(j)/r_d(nj_d))**2)**alpha
               sum = sum + curden_d(j)*ds
c    >           zhalf*(dvoldr_p*curden_d(j)+dvoldr_m*curden_d(j-1))*drm
c    >           /(2.D0*pi_m*rmajor_d) 
             enddo
             cur0 = dabs(tocur_d)/sum
c
             do j=1,nj_d
               curden_d(j)=cur0*(1.D0-(r_d(j)/r_d(nj_d))**2)**alpha
             enddo
           endif 
c
         if (qohm_d(1).EQ.0) then
c
           if(iproc.eq.0) write(6,*)
     >       'Warning rep_iter: no ohmic power profile; making one up'
c
c...       calculate qohm_d from j(r) and loop voltage
c...       by: qohm = E.j
c...       Also calculate volume-integrated power 
c
           cursum=0.D0
           pohsum=0.D0
           do j=2,nj_d
             drm=r_d(j)-r_d(j-1)
             dvoldr_p=2.D0*pi_m*r_d(j)*
     >                2.D0*pi_m*rmajor_d*hcap_d(j)
             dvoldr_m=2.D0*pi_m*r_d(j-1)*
     >                2.D0*pi_m*rmajor_d*hcap_d(j-1)
             ds=2.D0*pi_m*(r_d(j)*hcap_d(j)+
     >          r_d(j-1)*hcap_d(j-1))/2.D0*drm
             cursum = cursum + curden_d(j)*ds
crew error fixed 7.1.96
crew             qohm_d(j)=ABS(vsurf_d*curden_d(j))/(2.*pi_m*rmajor_d)
             qohm_d(j)=(vsurf_d*curden_d(j))/(2.D0*pi_m*rmajor_d)
             pohsum=pohsum+
     >         zhalf*(dvoldr_p*qohm_d(j)+dvoldr_m*qohm_d(j-1))*drm
           enddo
           qohm_d(1)=(vsurf_d*curden_d(1))/(2.D0*pi_m*rmajor_d) 
crew....vsurf_d or curden_d may not be reliable...
crew     better renorm by pohm_d
crew thus in the end we use only a shape of current profile
crew and assume dE/dr=0 ie steady state current profile
           if(pohm_d.gt.0.) then
            do j=1,nj_d
              qohm_d(j)=pohm_d/pohsum*qohm_d(j)
            enddo
           endif
crew why calc?           pohm_d = pohsum
         if(iproc.eq.0) write(6,'(a27,F10.3,a10,F7.3)') 
     >     'Total calculated Pohm [MW]:', pohsum*1.D-6,
     >     'Vloop = ',vsurf_d
         if(iproc.eq.0) write(6,'(a27,F10.3)') 
     >     'Total integrated Ip  [MA]:', cursum*1.D-6
c
         end if  !end making up ohmic profile
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c---avg major radius of each flux surface, meters, at elevation of mag. axis

c          read(niterdb,'(a)')stflg
          do j=1,nj_d
             rmajavnpsi_d(j) = u2d(j,15)   !meters
          enddo
c
c---avg minor radius of each flux surface,meters, at elevation of mag. axis
c
          do j=1,nj_d
             rminavnpsi_d(j) = u2d(j,16)   !meters
             if(rminavnpsi_d(1).lt.0) rminavnpsi_d(1)=1.D-6
          enddo  
c--- elongation of each flux surface
c
          do j=1,nj_d
           elongx_d(j)=u2d(j,26)
          enddo
c
c use 1D value for NTCC benchmarking if itest_ntcc.gt.0
c
         if( elongx_d(nj_d).lt.1.e-3) then
          if(itest_ntcc.gt.0) then
            write(*,32)
            do j=1,nj_d
              elongx_d(j)=u1d(6)
            enddo
          else
            do j=1,nj_d
              elongx_d(j)= (r_d(j)/rminavnpsi_d(j))**2
            enddo
          endif
         endif
c
c--- triangularity at each flux surface
          do j=1,nj_d
           deltax_d(j)=u2d(j,24)
          enddo 
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
 999  return
c
 89   continue
      if(iproc.eq.0) 
     >   write(*,*) 'Error reading iter namelist'
      close(50)
      goto 999
c     
 90   continue
      if(iproc.eq.0) 
     >   write(*,*) 'Error reading the 0d file'
      stop
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
 11   format(/,'0D quantities:',/)
 12   format(/,2x,'A-plasma = ',0p1f10.3,4x,'Z-plasma = ',0p1f10.3,/,
     >         2x,'A-beam   = ',0p1f10.3,4x,'Z-beam   = ',0p1f10.3,/,
     >         2x,'A-imp    = ',0p1f10.3,4x,'Z-imp    = ',0p1f10.3,/)
 13   format(2x,'R0  = ',0p1f10.3,4x,'a   = ',0p1f10.3,/,
     >       2x,'BT  = ',0p1f10.3,4x,'Ip  = ',0p1f10.3,/,
     >       2x,'ne  = ',0p1f10.3,4x,'Pnb = ',0p1f10.3,/,
     >       2x,'Pic = ',0p1f10.3,4x,'Pec = ',0p1f10.3,/,
     >       2x,'Ti0 = ',0p1f10.3,4x,'Te0 = ',0p1f10.3,/,
     >       2x,'q0  = ',0p1f10.3,4x,'q95 = ',0p1f10.3,/,
     >       2x,'kappa = ',0p1f8.3,4x,'delta = ',0p1f8.3,/,
     >       2x,'Wth = ',0p1f10.3,4x,'Wtot= ',0p1f10.3)
 32   format(' Using 1D ufile value for elongation')
 35   format(' Warning: xp_time <= min exp. time, setting xp_time = ',
     >         0pf10.5)
 50   format(i2,2x,'rho',5x,'powibeam',5x,'powiion',6x,'powicx',7x,
     &       'powei',8x,'powiwdot',5x,'powitot')
 53   format(i2,2x,0p1f4.2,0p6e13.5)
 54   format(i2,2x,0p1f5.3,1p6e13.5)
 55   format(i2,2x,'rho',5x,'powebeam',5x,'poweohm',6x,'poweion',6x,
     &       'powerad',6x,'powei',8x,'powewdot',5x,'powetot')
 57   format('   ierr = ',i3,',ufile=',a50)
 58   format(2x,i2,2x,0p1f10.6,1p6e15.7)
 60   format(' Reading time-dependent 2D data ...')
 65   format(' Reading time-dependent 1D data ...')
 91   format(' Smoothing electron density profile ...',0p1f4.2)
 92   format(' Smoothing ion density profile ...',0p1f4.2)
 93   format(' Smoothing Zeff profile ...',0p1f4.2)
 94   format(' Smoothing q-profile ...',0p1f4.2)
 95   format(' Smoothing Qnb profile ...',0p1f4.2)
 96   format(' Smoothing fast ion density profile ...',0p1f4.2)
 97   format(' Smoothing vrot profile ...',0p1f4.2)
 98   format(' Smoothing torque profile ...',0p1f4.2)
101   format(' Smoothing grho profile ...',0p1f4.2)
100   format(i2,2x,0p1f4.2,1p8e13.5)
150   format(2x,i2,2x,0p1f10.6,1p6e15.7)
c
      end
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      integer function linblnk(string,max)
      character string*(*)
      integer i,j,max
      i=0
      do j=1,max
         if ( string(j:j)  .ne. ' ' ) then
            i = i + 1
         else
            go to 10
         endif
      enddo
 10   continue 
      linblnk=i
      return
      end

