      Module toq_12_interface

!    NOTE:             subroutine get_toq (which gets the path for a
!                      toq executeable) is kept in ext_prog_info.f90



        implicit none 

        !type toq_input: !cant use type because variables are in namelist
           integer, parameter :: kffp =10
           integer, parameter :: kppcoef=10, kpcoef=kppcoef+1,kcdcoeff =20
           integer,dimension(20) :: nbug       !20 hard wired in toq
           integer ieqdtoq,npsi_toq,nthet_toq,imislo, modelbnd,            &
                   ishape , modelf ,nffp, modelp,modelq,nppcoef,           &
                   ieqdsk_toq,minmg, maxmg,loopinmg,nhior, nbndry,         &
                   iteqmx,npts,nppspline, ieqdsk_toqp
           logical fixcur,prtbal,prtboot,prtgato, prtdcon,                 &
                   prteqdata,prtfixb,prtflux,igimblett,prtcamino,          &
                   prtmars,prtgato2,prtbalr,prtmercy,prtonetwo
           real *8 alpsi,eshape,xshape, betiso,  bavmg,tolmg,psifactr_toq, &
                   dsrat,cdx,cddelx,cdlocbump,cdwidthbump,zeff_toq,        &
                   qedge_toq
           real *8,dimension (kffp) :: ffpcoef,ppcoef
          integer                                                &
                  io12,iterillu,nthet,isym_toq,        &
                  iqval,jlimzero,islwst,nbet,      &
                  nbavxz,nrotcoef,   &
                  nppcutoff,modelcd,      &
                  ncdcoeff,ncdfreeze,ifillboot,nqspline,    &
                  nffpspline,ncdspline

                  
          real *8                                                &
                  spack,xpack,dxpack,ashape,bshape,cshape, &
                  rzero, dchi,qax,qaxp,qlim,  &
                  qlimp,alphaq,gammaq,betis1,bavq,bava,   &
                  bavxz,bavboot,bavcd,qaxwant, bug1,bug2,  &
                  bug3,ncent,nedge,nexpin,nexpout, tcent,tedge,  &
                  texpin,texpout,rcent,rexpin,rexpout,     &
                  betapwant,baxis0,cdnorm,cdexpin,        &
                  cdexpout,cdmerc,cboot,gamboot,jbjratio,   &
                  qpaxspline,qplimspline,paxis,xphalf,widthp,    &
                  nedge13,tedgeEV,rmax_toq,toleq_toq
   
          real *8, dimension(:) :: bug(20)
          real *8, dimension(:) :: rotcoef(kpcoef)
          real *8, dimension(:) :: cdcoeff(kcdcoeff)
          integer, parameter :: kqspline=50                   
          real *8, dimension(:) :: qspline(kqspline)
          real *8, dimension(:) :: xqspline(kqspline)
          integer, parameter :: kppspline=50                   
          real*8 , dimension(:) :: ppspline(kppspline)
          real*8 , dimension(:) :: xppspline(kppspline)
          integer, parameter :: kffpspline=50                    
          real*8 , dimension(:) :: ffpspline(kffpspline)
          real*8 , dimension(:) :: xffpspline(kffpspline)
          integer, parameter :: kcdspline=50                     
          real*8 , dimension(:) :: cdspline(kcdspline)
          real*8 , dimension(:) :: xcdspline(kcdspline)
          real*8 , dimension(:),allocatable ::psinarray,ffparray,pparray
          logical tridirect,cdfreeze,fixcurwop,jbs
          character*10 jbtype            !len 10 from toq +2 for quotes
          character(len = 9) equiltype  !len =7 from toq +2 for quotes
          character(len=4) betafix      !len =2          +2
          character(len=64) fneqdsk
          character(len=3) updownsym    !len =1          +2
          character(len = 256) :: runid_toq  ='Toq input file created by Onetwo'
          character(len = 256) :: toq_input_file ='intoq'
          character(len = 256) :: toq_base_out


               common /dskeqdat/ modelcd, modelbnd



           data zeff_toq,npsi_toq,nthet_toq,qedge_toq,psifactr_toq,    &
                isym_toq,ieqdsk_toq                                    &
              /1.5,65,67,0.,0.99999,1,0 /
        !end type toq_input

        contains  


        subroutine  wrt_toq_input(totcur)
!-----------------------------------------------HSJ-8/01/03-------

          USE mhdpar, only: nw
!          USE ename,  only: eqdskfilename
          USE io,     only: nitre
          integer npsi,isym,ieqdsk,npts
          real *8 rmax,toleq,zeff,totcur,redge,qedge,psifactr
!         common /dskeqdat/ modelcd, modelbnd
       namelist /all/                                             &
         npsi,nthet,alpsi,spack,xpack,dxpack,isym,ieqdtoq ,       &
         ishape,ashape,bshape,cshape,eshape,xshape,tridirect,     &
         rzero,rmax, dchi,                                        &
         qax,qaxp,qlim,qlimp,iqval,alphaq,gammaq,jlimzero,        &
         islwst,nbet,betis1,betiso,                               &
         ieqdsk,iteqmx,bavq,bava,bavxz,toleq,                     &
         imislo,nbavxz,bavboot,bavcd,qaxwant,                     &
         nhior,nbndry,                                            &
         bug1,bug2,bug3,bug,nbug,                                 &
         modelp,nppcoef,ppcoef,nrotcoef,rotcoef,                  &
         nppcutoff,                                               &
         modelq,modelf,modelcd,modelbnd,                          &
         ncent,nedge,nexpin,nexpout,                              &
         tcent,tedge,texpin,texpout,                              &
         rcent,redge,rexpin,rexpout,                              &
         betapwant,betafix,                                       &
         ffpcoef,totcur,                                          &
         nffp,equiltype,baxis0,                                   &
         prtmars,prtboot,prtgato,prtgato2,prtfixb,prtbalr,        &
         prtbal,prtmercy,prtdcon,prteqdata,prtonetwo,             &
         prtcamino,prtflux,                                       &
         ncdcoeff,cdcoeff,cdx,cddelx,cdlocbump,cdwidthbump,       &
         ncdfreeze,cdfreeze,cdnorm,cdexpin,cdexpout,cdmerc,       &
         cboot,zeff,gamboot,ifillboot,jbjratio,jbtype,fixcur,     &
         fixcurwop,jbs,                                           &
         qspline,xqspline,qpaxspline,qplimspline,nqspline,        &
         ppspline,xppspline,nppspline,                            &
         ffpspline,xffpspline,nffpspline,                         &
         paxis,xphalf,widthp,                                     &
         cdspline,xcdspline,ncdspline,                            &
         loopinmg,bavmg,minmg,maxmg,tolmg,nedge13,tedgeEV         


         namelist /input/fneqdsk,psifactr,dsrat,npts,updownsym,   &
                   qedge,igimblett

          call getioun(io12,42)
          open  (unit = io12, file = toq_input_file, status = 'UNKNOWN')

!          fneqdsk = '"'//fneqdsk(1:LEN_TRIM(eqdskfilename))//'"'
!          fneqdsk ='"'//eqdskfilename(1:LEN_TRIM(eqdskfilename))//'"'
!          fneqdsk = ADJUSTL(fneqdsk)
!          print *,'Toq input eqdsk is =',fneqdsk(1:len_trim(fneqdsk))
          write(nitre,FMT='("Toq input eqdsk =",a)') &
              fneqdsk(1:len_trim(fneqdsk))
          !set those values that have name clashes with Onetwo:
          zeff = zeff_toq
          rmax = rmax_toq
          npsi = npsi_toq
          isym = isym_toq
          toleq = toleq_toq
          nthet = nthet_toq
          ieqdsk =ieqdsk_toq
          redge = -1.              !need to define
          prtonetwo = .true.


         qedge = qedge_toq
         psifactr = psifactr_toq
         npts = nw


          write (unit = io12, fmt = '(3x, a)') runid_toq
          write (unit = io12, nml = all)
          write (unit = io12, nml = input)
          close (unit = io12)
          call giveupus(io12)


        end subroutine wrt_toq_input







        subroutine  set_toq_default_input
!-----------------------------------------------HSJ-8/05/03-------
!     code copied from init.f (init.f is a toq source file)
!     totcur,isym,npsi,rmax,redge,ieqdsk,toleq,zeff  are 
!     variables that have the same name in Onetwo and
!     consequently need special treatment
!     (cant use f90 type because it doenst work with namelists)
!     set namelist defaults and then read in namelist
!------------------------------------------------------------------------
      integer i
!     mesh parms
!      isym=0
      ieqdtoq=0
!      npsi=67
      nthet=65
      alpsi=0.
      spack=0.
      xpack=0.5
      dxpack=0.075
!     geometry parameters
      ishape=0
      ashape=16.649
      bshape=0.6
      cshape=1.111
      eshape=0.780
      xshape=5.007
      tridirect=.false.
!      rmax=10.
      rzero=169.5
!     poloidal flux
      dchi=3.e6
!     safety factor parameters
      qax=1.05
      qlim=4.
      qaxp=0.8
      qlimp=9.
      iqval=0                   !see qfunc;iqval<0 uses gammaq 
                                !to determine alphaq
      alphaq=3.0
      jlimzero=0                !if=1 alphaq adjusted to make edge current = 0
!     density parms
      ncent=1.
      nedge=0.
      nexpin=1.
      nexpout=1.
!     temperature parms
      tcent=1.
      tedge=0.
      texpin=1.
      texpout=1.
!     rotation temperature parms
      rcent=0.
!      redge=0.
      rexpin=1.
      rexpout=1.
!     slow start variables (see sub slwstr)
      islwst=0
      nbet=5
      betis1=0.
      betiso=0.
      betapwant=0.
      betafix="bt"              ! fix beta_tor; "bp" fix beta_pol; "no"
!     equilibrium control parms
!      ieqdsk=0
      iteqmx=25
      bavq=0.
      bava=0.
      bavxz=0.
      bavboot=0
      bavcd=0
!      toleq=1.e-6              !set in wrt_toq_input
      imislo=3
      nbavxz=1000
!     interpolation parameters
      nbndry=199
      nhior=15
!     debugging parameters
      bug1=0.
      bug2=0.
      bug3=0.
      do i=1,20
         bug(i)=0.
         nbug(i)=0
      end do
!     equilibrium solution type
      equiltype='qsolver'
!     ffprime stuff
      nffp=2
      ffpcoef(1)=-1
      ffpcoef(2)=1
!      totcur=1.e6
      baxis0=1.e4
      iterillu=20
!     printing switches
      prtmars=.false.
      prtboot=.false.
      prtgato=.false.
      prtgato2=.false.
      prtfixb=.false.
      prtbalr=.false.
      prtbal=.false.
      prtmercy=.false.
      prtdcon=.false.
      prteqdata=.false.
      prtcamino=.false.
      prtflux=.false.
!    
      modelp=1
      modelq=1
      modelf=1
      modelcd=3                 !use polynomial model as default
      modelbnd=1
      nppcoef=0
      nppcutoff=0               !yr
      do i=1,kppcoef
         ppcoef(i)=0.           !modelp=3
      end do
      nrotcoef=0
      do i=1,kppcoef
         rotcoef(i)=0.          !modelp=3
      end do
      paxis=1.0                 !modelp=4
      xphalf=0.8                !modelp=4
      widthp=0.5                !modelp=4
      jbs=.false.
!      cdone=1.                  !modelcd=1--not a namelist parameter
      cdexpin=2.                !modelcd=1
      cdexpout=3.               !modelcd=1
      cdmerc=0.                 !modelcd=1
      ncdcoeff=2                !modelcd=3
      cdcoeff(1)=1.
      cdcoeff(2)=0.
      do i=3,kcdcoeff
         cdcoeff(i)=0.
      end do
!yr      cdx=0.5
!yr      cddelx=0.02
      cdx=5.0                  !yr
      cddelx=0.04              !yr
      cdlocbump=1.e10
      cdwidthbump=1.
      ncdfreeze=0
      cdfreeze=.false.
      qaxwant=2.0
      cdnorm=-1.
      cboot=1.0
      gamboot=0.5
!      zeff=1.5                !set in write_toq_input
      nedge13=0.
      tedgeEV=200.
!      nustr_edge=0.
      ifillboot=0
      jbjratio=1.
      jbtype='separate'         !other option is 'combined'
      fixcur=.true.             !hold total current fixed
      fixcurwop=.false.         !if true fixcur WithOutPressure
                                !i.e., without adjusting pressure
! q values for modelq=3 spline option
      qspline(1)=1.
      qspline(2)=1.
      xqspline(1)=0.
      xqspline(2)=1.
      qpaxspline=0.
      qplimspline=0.
      nqspline=2
! pp values for modelp=5 spline option
      ppspline(1)=1.
      ppspline(2)=1.
      xppspline(1)=0.
      xppspline(2)=1.
      nppspline=2
! pp values for modelffp=5 spline option
      ffpspline(1)=1.
      ffpspline(2)=1.
      xffpspline(1)=0.
      xffpspline(2)=1.
      nffpspline=2
!      initf=1                   ! this is not namelist parameter, it is a
                                ! switch to make sure fppspline used only
                                ! on first call to fsetup
! cd values for modelcd=5 spline option
      cdspline(1)=1.
      cdspline(2)=1.
      xcdspline(1)=0.
      xcdspline(2)=1.
      ncdspline=2
!     multi-grid parameters
      loopinmg=2
      bavmg=0.5
      minmg=3
      maxmg=10
      tolmg=1.e-8

!     second toq namelist,input:
!      qedge=0.                    !set in write_toq_input
      igimblett=.false.
      updownsym = 'a'
      dsrat =0.005
      end subroutine set_toq_default_input




      subroutine toq_drive(toq_current)
!------------------------------------------------------------------------
!
!-- driver to spawn Toq fixed boundary inverse equilibirum solver:
!
!------------------------------------------------------------------------

       !USE param
       !USE mhdpar 
       USE  ext_prog_info,     only: get_toq, toq_path
!       USE  toq_12_interface,  only: toq_input_file, wrt_toq_input
       USE io,only : ncrt,nout
       implicit none 
       real *8 toq_current
       integer len_str,ishell
       character(len = 256),save :: toq_path_out
       character(len =256)       :: command
       logical,save ::  first_time 
       data first_time / .true./




      if(first_time)then
         !get fully qualified name of toq to run:
          call get_toq(ncrt,nout,toq_path_out,len_str)
          first_time = .false.
      endif


      call wrt_toq_input(toq_current)


      command = ADJUSTL(toq_path_out(1:LEN_TRIM( toq_path_out)))
      command = command(1:LEN_TRIM(command))//' '//              &
         toq_input_file(1:LEN_TRIM(toq_input_file))


!     execute inverse equilibrium code TOQ:
      write  (6, '(/ '' ---- TOQ started'')')
      write(6, FMT ='(" running : ",a)')command
!
      if (ISHELL(command) .lt. 0)                                   &
       call STOP ('subroutine TOQ_DRIVE: failure of spawned TOQ', 67)
!
      write  (6, '(/ '' ---- TOQ  finished'')')

      return
      end subroutine toq_drive


      subroutine toqfile(rcfixb,raxtoq,zaxtoq,btorfixb,totcurtoq,    &
        psinarray,ffparray,pparray,npsitoq,rboundary,zboundary,nbound)
!=======================================================================
!     read file dskfixb  created by toq (taken from toq distribution HSJ)
!=======================================================================
      implicit none
      integer,intent(out):: npsitoq,nbound
      integer  nfixb,i
      real*8, dimension(:), intent(out) :: rboundary,zboundary
      real*8, dimension(:), intent(out) :: psinarray,ffparray,pparray
      real*8, intent(out)::  btorfixb,rcfixb,raxtoq,zaxtoq,totcurtoq
      character*80 line
!     
!     open output file
      nfixb = 42
      call getioun(nfixb,nfixb)
      open(unit=nfixb,file="dskfixb",status="old")
      read(nfixb,'(80a)') line
      read(nfixb,'(80a)') line
      read(nfixb,100) npsitoq, nbound
 100  format(2i5)
      write(6,'(2i5)') npsitoq,nbound
      read(nfixb,'(80a)') line
      read(nfixb,200) rcfixb,raxtoq,zaxtoq,btorfixb
      read(nfixb,200) totcurtoq
      read(nfixb,'(80a)') line
      read(nfixb,201) (psinarray(i),i=1,npsitoq) ! normalized psi 0-1
      write(6,'("read psi")')
      read(nfixb,'(80a)') line
      read(nfixb,201) (ffparray(i),i=1,npsitoq) ! ff'
      write(6,'("read ffp")')
      read(nfixb,'(80a)') line
      read(nfixb,201) (pparray(i),i=1,npsitoq) ! p'
      write(6,'("read pprime")')
      read(nfixb,'(80a)') line
      read(nfixb,200) (rboundary(i),i=1,nbound)
      write(6,'("read rboundary")')
      read(nfixb,'(80a)') line
      read(nfixb,200) (zboundary(i),i=1,nbound)
      write(6,'("read zboundary")')

      close(unit = nfixb)
      call giveupus(nfixb)


 200  format(4e19.12)
 201  format(4e19.0)
      return
      end subroutine toqfile

 
      end module toq_12_interface
