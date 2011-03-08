      subroutine cdfread (igrid)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c  07-Aug-08 J. Kinsey, GA
c  This subroutine reads Netcdf files from ONETWO
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      implicit none
c
      include 'mpif.h'
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/glf.m'
      include '../inc/ptor.m'
c
      include 'netcdf.inc'
c
      integer igrid, nfmax, nfmaxatts, dim_rho, dim_npsi, dim_ion
      parameter (nfmax=500,nfmaxatts=500,dim_rho=301,dim_npsi=501)
      parameter (dim_ion=4)
      integer ncid, status, ndims, nvars, natts, recdim, rcode
      integer i, j, k, dimid, numrecs, varid, lenstr
      integer nctype, ndsize, nvdim, iflag
      integer ichmax, idchar, ifchar, ichar, nex, npsi_d
      character*15 extension
      character*50 iterdbfile
      character cdfile*130
      character*120 msg
      character*1024 input_file, strbuf
      integer dimsiz(nfmax), vartyp(nfmax),
     &        nvdims(nfmax), vdims(nfmax), nvatts(nfmaxatts)
      integer vvdims(nfmax,nfmax), attype(nfmax,nfmaxatts),
     &        attlen(nfmax,nfmaxatts)
      integer start(2), count(2)
      integer*4 dnamsz(nfmax)
      character*128 dimnam(nfmax), varnam(nfmax)
      character*128 attnam(nfmax,nfmaxatts)
      character*11 avt(nfmax), avartyp(6)
      real*8 aj, temp(dim_npsi)
      real*8 rho_mhd(dim_npsi)
      real*8 rminavnpsi_mhd(dim_npsi),
     &       rmajavnpsi_mhd(dim_npsi), elong_mhd(dim_npsi),
     &       deltau_mhd(dim_npsi), deltal_mhd(dim_npsi), 
     &       delta_mhd(dim_npsi), psivolp_mhd(dim_npsi),
     &       sfareanpsi_mhd(dim_npsi), cxarea_mhd(dim_npsi),
     &       grho1_mhd(dim_npsi), grho2_mhd(dim_npsi)
      real*8 rminavnpsi_ex(dim_npsi),
     &       rmajavnpsi_ex(dim_npsi), elong_ex(dim_npsi),
     &       deltau_ex(dim_npsi), deltal_ex(dim_npsi), 
     &       delta_ex(dim_npsi), psivolp_ex(dim_npsi),
     &       sfareanpsi_ex(dim_npsi), cxarea_ex(dim_npsi),
     &       grho1_ex(dim_npsi), grho2_ex(dim_npsi)
      real*8 rho_ex(dim_rho), te_ex(dim_rho), ti_ex(dim_rho),
     &   ne_ex(dim_rho), nem_ex(dim_rho), nfast_ex(dim_rho),
     &   en_ex(dim_rho,dim_ion), enm_ex(dim_rho,2),
     &   z_ex(dim_rho,dim_ion),
     &   zeff_ex(dim_rho),rmin_ex(dim_rho),rmaj_ex(dim_rho),
     &   q_ex(dim_rho), ptot_ex(dim_rho), pfast_ex(dim_rho),
     &   angrot_ex(dim_rho),curden_ex(dim_rho),torque_ex(dim_rho),
     &   qbeame_ex(dim_rho),qbeami_ex(dim_rho), 
     &   qrfe_ex(dim_rho), qrfi_ex(dim_rho),
     &   qione_ex(dim_rho), qioni_ex(dim_rho), qcx_ex(dim_rho),
     &   qrad_ex(dim_rho), qohm_ex(dim_rho), qdelt_ex(dim_rho),
     &   qfuse_ex(dim_rho), qfusi_ex(dim_rho),
     &   qpedtc_ex(dim_rho), qpidtc_ex(dim_rho), 
     &   fcap_ex(dim_rho), gcap_ex(dim_rho), hcap_ex(dim_rho),
     &   enn_ex(dim_rho), ennw_ex(dim_rho), ennv_ex(dim_rho),
     &   sbeam_ex(dim_rho), sion_ex(dim_rho), srecom_ex(dim_rho),
     &   scx_ex(dim_rho), sbcx_ex(dim_rho), s_ex(dim_rho),
     &   dudtsv_ex(dim_rho), dpedtc_ex(dim_rho), dpidtc_ex(dim_rho)
c
      integer varnnj, varnnion, varnnprim, varnnimp, varnrmajor,
     &   varnkappa, varndelta, varnvolume, varnarea, varnbt, varnip,
     &   varnteo, varntio
      integer varnrhog, varnte, varnti, varnne, varnnf, varnzeff,
     &   varnz, varnni, varnq, varnptot, varnpfst, varnang, varncur,
     &   varntoq, varnneu, varnnw, varnnv, varnsion, varnsrec,
     &   varnscx, varnsbcx, varnstot, varnsdot, varnsbeam, varnfcap, 
     &   varngcap, varnhcap, varnrmhd, varnrmaj, varnrmin, varnkap,
     &   varndelu, varndell, varnvol, varnsfa, varncxa, varngr1,
     &   varngr2, varnqbe, varnqbi, varnqrfe, varnqrfi,
     &   varnqione, varnqioni, varnqcx, varnqrad, varnqohm,
     &   varnqfuse, varnqfusi, varnwde, varnwdi, varnqdelt
c
      real, allocatable, dimension(:,:) :: rhob_nc
      real, allocatable, dimension(:) :: zfluxlim_nc, xnstari_nc,
     &    xnstare_nc,xkapi_nc,zrhoi_nc,ztaui_nc,delta_nc,xkmneo_d

      allocate ( rhob_nc(nj,3), zfluxlim_nc(nj), xnstari_nc(nj), 
     &   xnstare_nc(nj), xkapi_nc(nj), zrhoi_nc(nj), ztaui_nc(nj),
     &   delta_nc(nj), xkmneo_d(nj) )
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      extension = shot
      input_file = 'iterdbnc.'//extension
c
c... count characters in directory name -> idchar
c
      ichmax = 60
      idchar = 0
      do j=1,ichmax
        if ( cudir(j:j)  .ne. ' ' ) then
          idchar = idchar + 1
        else
          go to 5
        endif
      enddo
 5    continue
c
c
c... count characters in iterdb filename -> ifchar
c
      ifchar = 0
      do j=1,ichmax
        if ( input_file(j:j)  .ne. ' ' ) then
          ifchar = ifchar + 1
        else
          go to 6
        endif
      enddo
  6   continue
c
c      write(*,*) 'file = ',input_file
c      write (*,*) 'idchar = ',idchar
c      write (*,*) 'ifchar = ',ifchar
      if (    cudir(idchar:idchar) .ne. '/') then
         ichar = idchar + 1 + ifchar
         cdfile = cudir(1:idchar) // '/' // input_file(1:ifchar)
      else
         ichar  = idchar + ifchar
         cdfile = cudir(1:idchar) // input_file(1:ifchar)
      endif
      input_file=cdfile(1:ichar)
      if(i_proc.eq.0) write (*,*) ' iterdbnc = ',input_file
c
c... open netcdf file
c
      status = nf_open(input_file,0,ncid)
      if (status.ne.nf_noerr) call handle_err(status)
c
c... inquire about dimensions, variables, attributes
c
      call ncinq(ncid,ndims,nvars,natts,recdim,rcode)
c      write(*,*) 'nvars = ',nvars,ndims,rcode,ncid
c
      do i=1,ndims
        dimid=i
        call ncdinq(ncid,dimid,dimnam(i),dimsiz(i),rcode)
        dnamsz(i)=index(dimnam(i),' ')-1
        if(recdim.ne.-1) numrecs=dimsiz(recdim)
c        write(*,50) i,dimnam(i),dimsiz(i)
 50    format(' dimension id= ',i3,' name= ',a31,' size= ',i3)
      enddo
c
c... read in variables
c
      do i=1,nvars
        varid=i
        call ncvinq(ncid,varid,varnam(i),vartyp(i),nvdims(i),
     &         vdims,nvatts(i),rcode)
        if(nvdims(i).ne.0) then
          do k=1,nvdims(i)
            vvdims(i,k)=vdims(k)
          enddo
        endif
        write(*,100) varid,varnam(i),vartyp(i),nvdims(i),
     &               nvatts(i),vdims(i)
 100   format(' var id= ',i3,' varnam= ',a10,' vartyp= ',i1,
     &     ' nvdims= ',i1,' num atts=',i2,' vdims= ',i3)
      enddo
c
c... get info on variable attributes
c
      do i=1,nvars
        varid=i
        do k=1,nvatts(i)
          call ncanam(ncid,varid,k,attnam(i,k),rcode)
          call ncainq(ncid,varid,attnam(i,k),attype(i,k),
     &                attlen(i,k),rcode)
c        write(*,150) varnam(i), nvatts(i), attnam(i,k),
c     &               attype(i,k),attlen(i,k)
        enddo
 150    format(' variable= ',a31,' attribute #',i2,' is: ',
     &        a31,'type= ',i3,' len= ',i3)
      enddo
c
c... store character name of vartyp in avt
c
      do i=1,nvars
        varid=i
        avt(i)=avartyp(vartyp(varid))
      enddo
c
c... store varid number for each 1D varname
c
      do i=1,nvars
        if(varnam(i).eq.'nj') varnnj=i
        if(varnam(i).eq.'nion') varnnion=i
        if(varnam(i).eq.'nprim') varnnprim=i
        if(varnam(i).eq.'nimp') varnnimp=i
        if(varnam(i).eq.'rmajor') varnrmajor=i
        if(varnam(i).eq.'kappa') varnkappa=i
        if(varnam(i).eq.'deltao') varndelta=i
        if(varnam(i).eq.'volume') varnvolume=i
        if(varnam(i).eq.'areao') varnarea=i
        if(varnam(i).eq.'btor') varnbt=i
        if(varnam(i).eq.'tot_cur') varnip=i
        if(varnam(i).eq.'te0') varnteo=i
        if(varnam(i).eq.'ti0') varntio=i
      enddo
c
c... store varid number for each 2D varname
c
      do i=1,nvars
        if(varnam(i).eq.'rho_grid') varnrhog=i
        if(varnam(i).eq.'Te') varnte=i
        if(varnam(i).eq.'Ti') varnti=i
        if(varnam(i).eq.'ene') varnne=i
        if(varnam(i).eq.'enbeam') varnnf=i
        if(varnam(i).eq.'zeff') varnzeff=i
        if(varnam(i).eq.'z') varnz=i
        if(varnam(i).eq.'enion') varnni=i
        if(varnam(i).eq.'q_value') varnq=i
        if(varnam(i).eq.'press') varnptot=i
        if(varnam(i).eq.'pressb') varnpfst=i
        if(varnam(i).eq.'angrot') varnang=i
        if(varnam(i).eq.'curden') varncur=i
        if(varnam(i).eq.'storqueb') varntoq=i
        if(varnam(i).eq.'enn') varnneu=i
        if(varnam(i).eq.'ennw') varnnw=i
        if(varnam(i).eq.'ennv') varnnv=i
        if(varnam(i).eq.'sion') varnsion=i
        if(varnam(i).eq.'srecom') varnsrec=i
        if(varnam(i).eq.'scx') varnscx=i
        if(varnam(i).eq.'sbcx') varnsbcx=i
        if(varnam(i).eq.'stsource') varnstot=i
        if(varnam(i).eq.'dudtsv') varnsdot=i
        if(varnam(i).eq.'sbeam') varnsbeam=i
        if(varnam(i).eq.'fcap') varnfcap=i
        if(varnam(i).eq.'gcap') varngcap=i
        if(varnam(i).eq.'hcap') varnhcap=i
      enddo
c
      do i=1,nvars
        if(varnam(i).eq.'rho_mhd_gridnpsi') varnrmhd=i
        if(varnam(i).eq.'rmajavnpsi') varnrmaj=i
        if(varnam(i).eq.'rminavnpsi') varnrmin=i
        if(varnam(i).eq.'elongxnpsi') varnkap=i
        if(varnam(i).eq.'triangnpsi_u') varndelu=i
        if(varnam(i).eq.'triangnpsi_l') varndell=i
        if(varnam(i).eq.'psivolpnpsi') varnvol=i
        if(varnam(i).eq.'sfareanpsi') varnsfa=i
        if(varnam(i).eq.'cxareanpsi') varncxa=i
        if(varnam(i).eq.'grho1npsi') varngr1=i
        if(varnam(i).eq.'grho2npsi') varngr2=i
      enddo
      do i=1,nvars
        if(varnam(i).eq.'qbeame') varnqbe=i
        if(varnam(i).eq.'qbeami') varnqbi=i
        if(varnam(i).eq.'qrfe') varnqrfe=i
        if(varnam(i).eq.'qrfi') varnqrfi=i
        if(varnam(i).eq.'qione') varnqione=i
        if(varnam(i).eq.'qioni') varnqioni=i
        if(varnam(i).eq.'qcx') varnqcx=i
        if(varnam(i).eq.'qrad') varnqrad=i
        if(varnam(i).eq.'qohm') varnqohm=i
        if(varnam(i).eq.'qfuse') varnqfuse=i
        if(varnam(i).eq.'qfusi') varnqfusi=i
        if(varnam(i).eq.'dpedtc') varnwde=i
        if(varnam(i).eq.'dpidtc') varnwdi=i
        if(varnam(i).eq.'dpedt') varnwde=i
        if(varnam(i).eq.'dpidt') varnwdi=i
        if(varnam(i).eq.'qdelt') varnqdelt=i
      enddo
c
      write(*,*) 'varnnj = ',varnnj
      write(*,*) 'varnnion = ',varnnion
      write(*,*) 'varnnprim = ',varnnprim
      write(*,*) 'varnnimp = ',varnnimp
      write(*,*) 'varnrmajor = ',varnrmajor
      write(*,*) 'varnkappa = ',varnkappa
      write(*,*) 'varndelta = ',varndelta
      write(*,*) 'varnvolume = ',varnvolume
      write(*,*) 'varnarea = ',varnarea
      write(*,*) 'varnbt = ',varnbt
      write(*,*) 'varnip = ',varnip
      write(*,*) 'varnte0 = ',varnteo
      write(*,*) 'varnti0 = ',varntio
c
      write(*,*) 'varnrhog = ',varnrhog
      write(*,*) 'varnte = ',varnte
      write(*,*) 'varnti = ',varnti
      write(*,*) 'varnne = ',varnne
      write(*,*) 'varnnf = ',varnnf
      write(*,*) 'varnzeff = ',varnzeff
      write(*,*) 'varnz = ',varnz
      write(*,*) 'varnni = ',varnni
      write(*,*) 'varnq = ',varnq
      write(*,*) 'varnptot = ',varnptot
      write(*,*) 'varnpfst = ',varnpfst
      write(*,*) 'varnang = ',varnang
      write(*,*) 'varncur = ',varncur
      write(*,*) 'varntoq = ',varntoq
      write(*,*) 'varnneu = ',varnneu
      write(*,*) 'varnnw = ',varnnw
      write(*,*) 'varnnv = ',varnnv
      write(*,*) 'varnsion = ',varnsion
      write(*,*) 'varnsrec = ',varnsrec
      write(*,*) 'varnscx = ',varnscx
      write(*,*) 'varnsbcx = ',varnsbcx
      write(*,*) 'varnstot = ',varnstot
      write(*,*) 'varnsdot = ',varnsdot
      write(*,*) 'varnsbeam = ',varnsbeam
      write(*,*) 'varnfcap = ',varnfcap
      write(*,*) 'varngcap = ',varngcap
      write(*,*) 'varnhcap = ',varnhcap
      write(*,*) 'varnrmhd = ',varnrmhd
      write(*,*) 'varnrmaj = ',varnrmaj
      write(*,*) 'varnrmin = ',varnrmin
      write(*,*) 'varnkap = ',varnkap
      write(*,*) 'varndelu = ',varndelu
      write(*,*) 'varndell = ',varndell
      write(*,*) 'varnvol = ',varnvol
      write(*,*) 'varnsfarea = ',varnsfa
      write(*,*) 'varncxarea = ',varncxa
      write(*,*) 'varngrho1 = ',varngr1
      write(*,*) 'varngrho2 = ',varngr2
      write(*,*) 'varnqbe = ',varnqbe
      write(*,*) 'varnqbi = ',varnqbi
      write(*,*) 'varnqrfe = ',varnqrfe
      write(*,*) 'varnqrfi = ',varnqrfi
      write(*,*) 'varnqione = ',varnqione
      write(*,*) 'varnqioni = ',varnqioni
      write(*,*) 'varnqcx = ',varnqcx
      write(*,*) 'varnqrad = ',varnqrad
      write(*,*) 'varnqohm = ',varnqohm
      write(*,*) 'varnqfuse = ',varnqfuse
      write(*,*) 'varnqfusi = ',varnqfusi
      write(*,*) 'varnwdote = ',varnwde
      write(*,*) 'varnwdoti = ',varnwdi
      write(*,*) 'varnqdelt = ',varnqdelt
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c Read in 1D quantities
c
c.. get no. grid pts
c
      call ncvinq(ncid,varnnj,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      call ncvgt(ncid,varnnj,1,1,nj_d,rcode)
      write(*,*) 'nj = ',nj_d
c
c.. get no. of ion species
c
      call ncvinq(ncid,varnnion,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      call ncvgt(ncid,varnnion,1,1,nion_d,rcode)
      write(*,*) 'nion = ',nion_d
c
c.. get no. of primary ion species
c
      call ncvinq(ncid,varnnprim,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      call ncvgt(ncid,varnnprim,1,1,nprim_d,rcode)
      write(*,*) 'nprim = ',nprim_d
c
c.. get no. of impurity ion species
c
      call ncvinq(ncid,varnnimp,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      call ncvgt(ncid,varnnimp,1,1,nimp_d,rcode)
      write(*,*) 'nimp = ',nimp_d
c
c.. get R_0
c
      call ncvinq(ncid,varnrmajor,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      call ncvgt(ncid,varnrmajor,1,1,rmajor_d,rcode)
      write(*,*) 'R_0 = ',rmajor_d
c
c.. get kappa
c
      call ncvinq(ncid,varnkappa,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      call ncvgt(ncid,varnkappa,1,1,kappa_d,rcode)
      write(*,*) 'kappa = ',kappa_d
c
c.. get delta (upper) on axis
c
      call ncvinq(ncid,varndelta,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      call ncvgt(ncid,varndelta,1,1,deltao_d,rcode)
      write(*,*) 'deltao = ',deltao_d
c
c.. get volume
c
      call ncvinq(ncid,varnvolume,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      call ncvgt(ncid,varnvolume,1,1,volo_d,rcode)
      write(*,*) 'volo = ',volo_d
c
c.. get cross sectional area
c
      call ncvinq(ncid,varnarea,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      call ncvgt(ncid,varnarea,1,1,areao_d,rcode)
      write(*,*) 'areao = ',areao_d
c
c.. get Bt, vacuum toroidal field
c
      call ncvinq(ncid,varnbt,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      call ncvgt(ncid,varnbt,1,1,btor_d,rcode)
      write(*,*) 'Btor = ',btor_d
c
c.. get total plasma current
c
      call ncvinq(ncid,varnip,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      call ncvgt(ncid,varnip,1,1,tocur_d,rcode)
      write(*,*) 'Ip = ',tocur_d
c
c.. get Te0
c
      call ncvinq(ncid,varnteo,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      call ncvgt(ncid,varnteo,1,1,te0_d,rcode)
c      write(*,*) 'Te0 = ',te0_d
c
c.. get Ti0
c
      call ncvinq(ncid,varntio,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      call ncvgt(ncid,varntio,1,1,ti0_d,rcode)
c      write(*,*) 'Ti0 = ',ti0_d
c
c... time at which data is printed          
      time_d=xp_time
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c Read in profiles on onetwo transport grid (rho_ex)
c
      call ncvinq(ncid,varnrhog,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      write(*,*) 'nvdim = ',nvdim
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
        write(*,*) 'start,count = ',start(i),count(i)
      enddo
      call ncvgt(ncid,varnrhog,start,count,rho_ex,rcode)
c
c.. get Te
c
      call ncvinq(ncid,varnte,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      write(*,*) 'nvdim = ',nvdim
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
        write(*,*) 'start,count = ',start(i),count(i)
      enddo
      call ncvgt(ncid,varnte,start,count,te_ex,rcode)
c
c.. get Ti
c
      call ncvinq(ncid,varnti,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnti,start,count,ti_ex,rcode)
c
c.. get ne
c
      call ncvinq(ncid,varnne,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnne,start,count,ne_ex,rcode)
c
c.. get nfast
c
      call ncvinq(ncid,varnnf,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnnf,start,count,nfast_ex,rcode)
c
c.. get Zeff
c
      call ncvinq(ncid,varnzeff,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnzeff,start,count,zeff_ex,rcode)
c
c.. get charges (d,t,c,he)
c
      call ncvinq(ncid,varnz,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      write(*,*) 'ndsize = ',ndsize
      call ncvgt(ncid,varnz,start,count,z_ex(1:101,:),rcode)
c
c.. get ni assuming (d,t,c,he)
c
      call ncvinq(ncid,varnni,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      write(*,*) 'nvdim = ',nvdim
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
        write(*,*) 'start,count (enion) = ',start(i),count(i)
      enddo
      write(*,*) 'ndsize = ',ndsize
      call ncvgt(ncid,varnni,start,count,en_ex(1:101,:),rcode)
c
c.. compute ni (enm_ex(i,1)) from ne subtracting off nfast, 
c   n_c (en_ex(i,3)), and n_He (en_ex(i,4))
c   and nimp (enm_ex(i,2) assuming Z_DT=1 and Z_imp=6
c
      write(*,*) 'nion = ',nion_d
      if (nion_d.eq.4) then
        do i=1,nj_d
          enm_ex(i,1)= (ne_ex(i)-nfast_ex(i)-en_ex(i,3)*z_ex(i,3)-
     >                en_ex(i,4)*z_ex(i,4))
          enm_ex(i,2)= (ne_ex(i)*zeff_ex(i)-
     >                enm_ex(i,1)-nfast_ex(i))/36.D0
          nem_ex(i)=en_ex(i,1)*z_ex(i,1)+en_ex(i,2)*z_ex(i,2)+
     >          en_ex(i,3)*z_ex(i,3)+en_ex(i,4)*z_ex(i,4)+nfast_ex(i)
        enddo
c
c        do i=1,nj_d
c          write(*,30) i, rho_ex(i), enm_ex(i,1), en_ex(i,1)+en_ex(i,2)
c     >      en_ex(i,3)*z_ex(i,3),en_ex(i,4)*z_ex(i,4),
c     >      ne_ex(i)
c           write(*,30) i, rho_ex(i),enm_ex(i,1),enm_ex(i,2),
c     >                 zeff_ex(i), (enm_ex(i,1)+
c     >                 en_ex(i,3)*z_ex(i,3)**2+
c     >                 en_ex(i,4)*z_ex(i,4)**2+nfast_ex(i))/ne_ex(i)
c        enddo
      else
        write(*,*) 'stopping cos less than 4 ions'
        stop
      endif
c
c.. get q
c
      call ncvinq(ncid,varnq,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnq,start,count,q_ex,rcode)
c
c.. get ptot, total pressure
c
      call ncvinq(ncid,varnptot,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnptot,start,count,ptot_ex,rcode)
c
c.. get pfast, fast ion pressure
c
      call ncvinq(ncid,varnpfst,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnpfst,start,count,pfast_ex,rcode)
c
c.. get angrot
c
      call ncvinq(ncid,varnang,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnang,start,count,angrot_ex,rcode)
c
c.. get total current density
c
      call ncvinq(ncid,varncur,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varncur,start,count,curden_ex,rcode)
c
c.. get torque density
c
      call ncvinq(ncid,varntoq,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varntoq,start,count,torque_ex,rcode)
c
c.. get neutral density
c
      call ncvinq(ncid,varnneu,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnneu,start,count,enn_ex,rcode)
c
c.. get neutral density from wall source
c
      call ncvinq(ncid,varnnw,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnnw,start,count,ennw_ex,rcode)
c
c.. get neutral density from volume source
c
      call ncvinq(ncid,varnnv,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnnv,start,count,ennv_ex,rcode)
c
c.. get source due to ionization
c
      call ncvinq(ncid,varnsion,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnsion,start,count,sion_ex,rcode)
c
c.. get source due to recombination
c
      call ncvinq(ncid,70,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,70,start,count,srecom_ex,rcode)
c
c.. get source due to cx thermal neutrals
c
      call ncvinq(ncid,varnscx,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnscx,start,count,scx_ex,rcode)
c
c.. get sink due to cx with beam neutrals 
c
      call ncvinq(ncid,varnsbcx,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnsbcx,start,count,sbcx_d,rcode)
c
c.. get total source rate 
c
      call ncvinq(ncid,varnstot,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnstot,start,count,s_ex,rcode)
c
c.. get sdot
c
      call ncvinq(ncid,varnsdot,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnsdot,start,count,dudtsv_ex,rcode)
c
c.. get beam thermal ion source
c
      call ncvinq(ncid,varnsbeam,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnsbeam,start,count,sbeam_ex,rcode)
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... Geometric quantities on onetwo transport grid
c    All other geo. quantities on onetwo MHD grid
c
c.. get rho grid
c
      call ncvinq(ncid,varnrhog,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      write(*,*) 'nvdim = ',nvdim
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
c        write(*,*) 'start,count = ',start(i),count(i)
      enddo
      call ncvgt(ncid,varnrhog,start,count,rho_ex,rcode)
c
c.. get fcap
c
      call ncvinq(ncid,varnfcap,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnfcap,start,count,fcap_ex,rcode)
c
c.. get gcap
c
      call ncvinq(ncid,varngcap,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varngcap,start,count,gcap_ex,rcode)
c
c.. get hcap
c
      call ncvinq(ncid,varnhcap,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnhcap,start,count,hcap_ex,rcode)
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c Geometric quantitites on onetwo MHD grid (rhomhd_ex)
c
c.. get rho_mhd, rho grid corresponding to mhd psival grid
c
      call ncvinq(ncid,varnrmhd,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnrmhd,start,count,rho_mhd,rcode)
c
c.. store no. mhd grid pts
c
      npsi_d=ndsize
      write(*,*) 'npsi = ',npsi_d
c
c.. reverse rhomhd data
c
      do i=1,npsi_d
        temp(i)=rho_mhd(npsi_d+1-i)
      enddo
      do i=1,npsi_d
        rho_mhd(i)=temp(i)
      enddo
c
c.. get rminor
c
      call ncvinq(ncid,varnrmin,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnrmin,start,count,rminavnpsi_mhd,rcode)
c
      do i=1,npsi_d
        temp(i)=rminavnpsi_mhd(npsi_d+1-i)
      enddo
      do i=1,npsi_d
        rminavnpsi_mhd(i)=temp(i)
      enddo
c
c.. get rmajor
c
      call ncvinq(ncid,varnrmaj,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnrmaj,start,count,rmajavnpsi_mhd,rcode)
c
      do i=1,npsi_d
        temp(i)=rmajavnpsi_mhd(npsi_d+1-i)
      enddo
      do i=1,npsi_d
        rmajavnpsi_mhd(i)=temp(i)
      enddo
c
c.. get kappa
c
      call ncvinq(ncid,varnkap,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnkap,start,count,elong_mhd,rcode)
c
      do i=1,npsi_d
        temp(i)=elong_mhd(npsi_d+1-i)
      enddo
      do i=1,npsi_d
        elong_mhd(i)=temp(i)
      enddo
c
c.. get delta_u, upper triangularity
c
      call ncvinq(ncid,varndelu,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varndelu,start,count,deltau_mhd,rcode)
c
      do i=1,npsi_d
        temp(i)=deltau_mhd(npsi_d+1-i)
      enddo
      do i=1,npsi_d
        deltau_mhd(i)=temp(i)
      enddo
c
c.. get delta_l, lower triangularity
c
      call ncvinq(ncid,varndell,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varndell,start,count,deltal_mhd,rcode)
c
      do i=1,npsi_d
        temp(i)=deltal_mhd(npsi_d+1-i)
      enddo
      do i=1,npsi_d
        deltal_mhd(i)=temp(i)
      enddo
c
c... avg delta
c
      do i=1,npsi_d
        delta_mhd(i)=(deltal_mhd(i)+deltau_mhd(i))/2.D0
      enddo
c
c.. get volume
c
      call ncvinq(ncid,varnvol,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnvol,start,count,psivolp_mhd,rcode)
c
      do i=1,npsi_d
        temp(i)=psivolp_mhd(npsi_d+1-i)
      enddo
      do i=1,npsi_d
        psivolp_mhd(i)=temp(i)
      enddo
c
c.. get surface area, 4*pi*pi*R0*hcap*rho*<ABS(grad rho)>
c
      call ncvinq(ncid,varnsfa,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnsfa,start,count,sfareanpsi_mhd,rcode)
c
      do i=1,npsi_d
        temp(i)=sfareanpsi_mhd(npsi_d+1-i)
      enddo
      do i=1,npsi_d
        sfareanpsi_mhd(i)=temp(i)
      enddo
c
c.. get cross sectional area (m^2)
c
      call ncvinq(ncid,varncxa,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varncxa,start,count,cxarea_mhd,rcode)
c
      do i=1,npsi_d
        temp(i)=cxarea_mhd(npsi_d+1-i)
      enddo
      do i=1,npsi_d
        cxarea_mhd(i)=temp(i)
      enddo
c
c.. get grho1
c
      call ncvinq(ncid,varngr1,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varngr1,start,count,grho1_mhd,rcode)
c
      do i=1,npsi_d
        temp(i)=grho1_mhd(npsi_d+1-i)
      enddo
      do i=1,npsi_d
        grho1_mhd(i)=temp(i)
      enddo
c
c.. get grho2
c
      call ncvinq(ncid,varngr2,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varngr2,start,count,grho2_mhd,rcode)
c
      do i=1,npsi_d
        temp(i)=grho2_mhd(npsi_d+1-i)
      enddo
      do i=1,npsi_d
        grho2_mhd(i)=temp(i)
      enddo
c
c      do i=1,npsi_d
c        write(*,30) i,rho_mhd(i),rmajavnpsi_mhd(i)
c      enddo
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... heating and loss power densities
c
c.. get qbeame
c
      call ncvinq(ncid,varnqbe,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnqbe,start,count,qbeame_ex,rcode)
c
c.. get qbeami
c
      call ncvinq(ncid,varnqbi,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnqbi,start,count,qbeami_ex,rcode)
c
c.. get qrfe
c
      call ncvinq(ncid,varnqrfe,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnqrfe,start,count,qrfe_ex,rcode)
c
c.. get qrfi
c
      call ncvinq(ncid,varnqrfi,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnqrfi,start,count,qrfi_ex,rcode)
c
c.. get qione
c
      call ncvinq(ncid,varnqione,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnqione,start,count,qione_ex,rcode)
c
c.. get qioni
c
      call ncvinq(ncid,varnqioni,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnqioni,start,count,qioni_ex,rcode)
c
c.. get qcx
c
      call ncvinq(ncid,varnqcx,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnqcx,start,count,qcx_ex,rcode)
c
c.. get qrad
c
      call ncvinq(ncid,varnqrad,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnqrad,start,count,qrad_ex,rcode)
c
c.. get qohm
c
      call ncvinq(ncid,varnqohm,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnqohm,start,count,qohm_ex,rcode)
c
c.. get qfuse,qfusi
c
      call ncvinq(ncid,varnqfuse,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnqfuse,start,count,qfuse_ex,rcode)
c
      call ncvinq(ncid,varnqfusi,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnqfusi,start,count,qfusi_ex,rcode)
c
c.. get wdot_e
c
      call ncvinq(ncid,varnwde,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnwde,start,count,dpedtc_ex,rcode)
c
c.. get wdot_i
c
      call ncvinq(ncid,varnwdi,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnwdi,start,count,dpidtc_ex,rcode)
c
c.. get qdelt
c
      call ncvinq(ncid,varnqdelt,strbuf,nctype,nvdim,vdims,
     &            nvatts,rcode)
      lenstr=1
      do i=1,nvdim
        call ncdinq(ncid,vdims(i),strbuf,ndsize,rcode)
        lenstr=lenstr*ndsize
        start(i)=1
        count(i)=ndsize
      enddo
      call ncvgt(ncid,varnqdelt,start,count,qdelt_ex,rcode)
c
c      do i=1,nj_d
c         write(*,30) i, rho_ex(i), qbeame_ex(i), qbeami_ex(i)
c      enddo
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c Map over to onetwo transport grid from MHD grid
c
      call w_lin_interp_r8(npsi_d,rho_mhd,rminavnpsi_mhd,
     &                    nj_d,rho_ex,rminavnpsi_ex,iflag,msg)
      call w_lin_interp_r8(npsi_d,rho_mhd,rmajavnpsi_mhd,
     &                    nj_d,rho_ex,rmajavnpsi_ex,iflag,msg)
      call w_lin_interp_r8(npsi_d,rho_mhd,elong_mhd,
     &                    nj_d,rho_ex,elong_ex,iflag,msg)
      call w_lin_interp_r8(npsi_d,rho_mhd,delta_mhd,
     &                    nj_d,rho_ex,delta_ex,iflag,msg)
      call w_lin_interp_r8(npsi_d,rho_mhd,psivolp_mhd,
     &                    nj_d,rho_d,psivolp_ex,iflag,msg)
      call w_lin_interp_r8(npsi_d,rho_mhd,sfareanpsi_mhd,
     &                    nj_d,rho_d,sfareanpsi_ex,iflag,msg)
      call w_lin_interp_r8(npsi_d,rho_mhd,cxarea_mhd,
     &                    nj_d,rho_d,cxarea_ex,iflag,msg)
      call w_lin_interp_r8(npsi_d,rho_mhd,grho1_mhd,
     &                    nj_d,rho_d,grho1_ex,iflag,msg)
      call w_lin_interp_r8(npsi_d,rho_mhd,grho2_mhd,
     &                    nj_d,rho_d,grho2_ex,iflag,msg)
c
c      do i=1,nj_d
c        write(*,30) i, rho_ex(i), delta_ex(i)
c      enddo
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c Map from onetwo grid to xptor grid
c
c... rho on xptor grid
c
      nex=mxgrid+1
      do j=1,nex
        aj=REAL(j-1)
        rho_d(j)=aj/REAL(jmaxm)
      enddo
c
      call w_lin_interp_r8(nj_d,rho_ex,rho_ex,
     &                    nex,rho_d,r_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,rminavnpsi_ex,
     &                    nex,rho_d,rminavnpsi_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,rmajavnpsi_ex,
     &                    nex,rho_d,rmajavnpsi_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,elong_ex,
     &                    nex,rho_d,elongx_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,delta_ex,
     &                    nex,rho_d,deltax_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,psivolp_ex,
     &                    nex,rho_d,psivolp_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,sfareanpsi_ex,
     &                    nex,rho_d,sfareanpsi_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,cxarea_ex,
     &                    nex,rho_d,cxareanpsi_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,grho1_ex,
     &                    nex,rho_d,grho1npsi_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,grho2_ex,
     &                    nex,rho_d,grho2npsi_d,iflag,msg)
c
      call w_lin_interp_r8(nj_d,rho_ex,te_ex,
     &                    nex,rho_d,te_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,ti_ex,
     &                    nex,rho_d,ti_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,ne_ex,
     &                    nex,rho_d,ene_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,enm_ex(:,1),
     &                    nex,rho_d,en_d(:,1),iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,enm_ex(:,2),
     &                    nex,rho_d,en_d(:,2),iflag,msg)
c
      call w_lin_interp_r8(nj_d,rho_ex,zeff_ex,
     &                    nex,rho_d,zeff_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,q_ex,
     &                    nex,rho_d,q_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,ptot_ex,
     &                    nex,rho_d,ptot_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,pfast_ex,
     &                    nex,rho_d,pfast_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,angrot_ex,
     &                    nex,rho_d,angrot_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,curden_ex,
     &                    nex,rho_d,curden_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,torque_ex,
     &                    nex,rho_d,torque_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,enn_ex,
     &                    nex,rho_d,enn_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,ennw_ex,
     &                    nex,rho_d,ennw_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,ennv_ex,
     &                    nex,rho_d,ennv_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,sion_ex,
     &                    nex,rho_d,sion_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,srecom_ex,
     &                    nex,rho_d,srecom_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,scx_ex,
     &                    nex,rho_d,scx_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,sbcx_ex,
     &                    nex,rho_d,sbcx_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,s_ex,
     &                    nex,rho_d,s_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,dudtsv_ex,
     &                    nex,rho_d,dudtsv_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,sbeam_ex,
     &                    nex,rho_d,sbeam_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,fcap_ex,
     &                    nex,rho_d,fcap_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,gcap_ex,
     &                    nex,rho_d,gcap_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,hcap_ex,
     &                    nex,rho_d,hcap_d,iflag,msg)
c
      call w_lin_interp_r8(nj_d,rho_ex,qbeame_ex,
     &                    nex,rho_d,qbeame_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,qbeami_ex,
     &                    nex,rho_d,qbeami_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,qrfe_ex,
     &                    nex,rho_d,qrfe_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,qrfi_ex,
     &                    nex,rho_d,qrfi_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,qione_ex,
     &                    nex,rho_d,qione_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,qioni_ex,
     &                    nex,rho_d,qioni_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,qcx_ex,
     &                    nex,rho_d,qcx_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,qrad_ex,
     &                    nex,rho_d,qrad_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,qohm_ex,
     &                    nex,rho_d,qohm_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,qfuse_ex,
     &                    nex,rho_d,qfuse_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,qfusi_ex,
     &                    nex,rho_d,qfusi_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,dpedtc_ex,
     &                    nex,rho_d,dpedtc_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,dpidtc_ex,
     &                    nex,rho_d,dpidtc_d,iflag,msg)
      call w_lin_interp_r8(nj_d,rho_ex,qdelt_ex,
     &                    nex,rho_d,qdelt_d,iflag,msg)
c
      amin_d=rminavnpsi_d(nex)
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c... rho on xptor grid
c
      nex=mxgrid+1
      nj_d=nex  ! reset nj_d to desired no. of grid pts
c      do j=1,nex
c        aj=REAL(j-1)
c        rho_d(j)=aj/REAL(jmaxm)
c      enddo
c
      do j=0,jmaxm
        aj=dfloat(j)
        rho(j)=aj/dfloat(jmaxm)
      enddo
      rho(0)=1.D-6
      rho(jmaxm)=rho(jmaxm)-1.D-6
c
c... rho grid, meters
c    from edge kappa
c
c      do j=1,nex
c        r_d(j)=rho_d(j)*amin_d*sqrt(kappa_d)
c        write(*,30) j, rho_d(j), r_d(j), qbeame_d(j), qbeami_d(j),
c     &              angrot_d(j)
c      enddo
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... Neoclassical transport (KAPISN only)
c
      if (use_xneo_m.eq.0 .or. use_xneo_m.eq.1) then
c
        if(i_proc.eq.0) write(*,30) use_xneo_m
        ng_nc=1                    !number of hyd. species
c       aplasm_nc=amassgas_exp     !array of atomic masses of hyd. species
        numzones_nc=nj_d            !number of zones
        btf_nc=DABS(btor_d)         !tf (tesla) fld at cntr of outer flux
        drshaf_nc=0.D0              !shaf shift of outer zone bdr. (cm)
        rminor_nc=rminavnpsi_d(nj_d)*1.D2
        rmajor_nc=rmajavnpsi_d(nj_d)*1.D2
c
        do j=1,ng_nc
         aplasm_nc(j)=amassgas_exp
        enddo
c
        do j=1,nj_d
          rhoel_nc(j)=dabs(ene_d(j)*1.D-6) !electron density (cm**-3)
          rhob_nc(j,1)=dabs(en_d(j,1)*1.D-6) !hyd. spec. den s (cm**-3)
          rhi_nc(j)=dabs(en_d(j,2)*1.D-6)  !z.c. array of av. impurity s density
          rhoi_nc(j)=rhi_nc(j)+rhob_nc(j,1) !z.c. array of total ion density
          te_nc(j)=dabs(te_d(j)*1.D3)      !z.c. array of Te (ev)
          ti_nc(j)=dabs(ti_d(j)*1.D3)      !z.c. array of Ti (ev)
          zeff_nc(j)=zeff_d(j)            !z.c. array of plasma zeff
          q_nc(j)=q_d(j)                  !z.c. array of safety factor
c
c...      Note that for the following arrays the rho-independent value is
c...      specified in mlt0in.
c
          aimp_nc(j)=aimp_nco
          xzimp_nc(j)=xzimp_nco
        enddo
c
        call kapisn(
     >          nkimod_nc,         !*kapai model nr desired
     >          aimp_nc,           !*atomic mass of av. impurity
     >          xzimp_nc,          !*atomic number of av. impurity
     >          aplasm_nc,         !array of atomic masses of hyd. species
     >          ng_nc,             !number of hyd. species
     >          nj_d,              !radial domain size
     >          rhoel_nc,          !zone centered electron density (cm**-3)
     >          rhob_nc,           !z.c. array of hyd. spec. den s (cm**-3)
     >          rhi_nc,            !z.c. array of av. impurity s density  
     >          rhoi_nc,           !z.c. array of total ion density       
     >          te_nc,             !z.c. array of Te (ev)
     >          ti_nc,             !z.c. array of Ti (ev)
     >          zeff_nc,           !z.c. array of plasma zeff
     >          numzones_nc,       !number of zones 
     >          q_nc,              !z.c. array of safety factor
     >          btf_nc,            !tf (tesla) fld at cntr of outer flux 
     >          drshaf_nc,         !shaf shift of outer zone bdr. (cm)
     >          rminor_nc,         !plasma minor radius (cm)
     >          rmajor_nc,         !major radius (cntr of outer flux) (cm)
     >          istringer_nc,      !Stringer correctioon
     >          xkapi_nc,          !o Neo-Class Ion thermal diff. (cm**2/sec)
     >          xnstari_nc,        !o nu-star-ions, 
     >          xnstare_nc,        !o nu-star-elecs, 
     >          ztaui_nc,          !o ion collision time
     >          zrhoi_nc,          !o ion poloidal gyro-radius
     >          zfluxlim_nc)       !o flux lim flow max temp grad length
c
        do j=1,nj_d
          xkapi_nc(j)=xkapi_nc(j)*1.D-4*en_d(j,1)/elongx_d(j) ! 1/(m*s)
          delta_nc(j)=rminavnpsi_d(j) / rmajor_d
          zrhoi_nc(j)=zrhoi_nc(j)*1.D-2 ! convert to meters
          zptineomax(j-1)=r_d(jmaxm+1)*100.D0/zfluxlim_nc(j)
        enddo
c
       if ((xkineo_d(nj_d).eq.0).or.(imyneoclass.eq.1)) then
         if (i_proc.eq.0) then
            write(6,*) 'Neoclassical xkineo calculated by kapisn'
         endif
         do j=1,nj_d
           xkineo_d(j)=xkapi_nc(j)
c fix this:
c           xkmneo_d(j)=0.1D0*delta_nc(j)**2.D0*
c     >                 (zrhoi_nc(j)**2.D0/ztaui_nc(j))
           xkmneo_d(j)=delta_nc(j)**1.5D0*xkineo_d(j)
c           write(*,*) j, rho_d(j), elongx_d(j), xkapi_nc(j)
c           write(*,60) j, rho(j), delta_nc(j), zrhoi_nc(j), ztaui_nc(j),
c     >                 xkmneo_d(j)/en_d(j,1)/grho2npsi_d(j),
c     >                 xkineo_d(j)/en_d(j,1)/grho2npsi_d(j)
         enddo
       endif
c
      endif
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      deallocate(rhob_nc, zfluxlim_nc, xnstari_nc, 
     >           xnstare_nc, xkapi_nc, zrhoi_nc, ztaui_nc)
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
 25   format(i3,2x,0p1f6.4,0p6f10.5)
 30   format(i3,2x,0p1f8.6,1p6e12.4)
 60   format(2x,i2,2x,0p1f10.6,1p6e15.7)
      end
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine handle_err(status)
      integer status
      if(status.ne.nf_noerr) then
        print *, nf_strerror(status)
        stop 'error, stopped'
      endif
      end
