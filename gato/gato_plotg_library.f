      subroutine plotg(kmesh)
c
c ----------------
c  plotting routine
c ----------------
c
      parameter (npx=200,ncx=2*npx,np1=npx+1,np2=npx+2,nc1=ncx+1)
      parameter (nxx=513,nxz=nxx)
      parameter (npp=929,ntt=2*npp-1,npb=npp)
      parameter (npq=np2)
      parameter (nw2=2*nxx,nh2=2*nxz,nvk0=nw2+nh2)
      parameter (nlx=1441)
      parameter (nvn=7)
c
      character*1   lbchpv, lbchsp, lbchrh, lbchro,
     &              lbchtf, lbchst ,lbchvl, lbchsv
      character*8   pclab,  tclab,  labpsi, 
     &              labpsiv,labpsin,labsqtp,
     &              labrhov,labrho2,labrhon,
     &              labtorf,labtorn,labsqtt,
     &              labvolm,labvoln,labsqtv,
     &              labdrh2,labdrho,labdtor,
     &              labdsqt,labdvol,labdsqv
      character*16  filout,filinp,fileql,filmap,filgrd,filrdm,filrdn,
     &              filfst,filplt
c
      common/cnst/  pi,hlfpi,twopi,amu,boltzm,echarg,pmass,roundff,bigno
      common/inpd/  ntor,ncase,idnsty,ndnxp0,ndnxp1,ndnxp2,ifastp,
     &              nmap,neqtyp,mpreset,ndoublt,ndivert,nmtype,
     &              jpsi,itht,jpsi1,jpsi2,ithtp,jsep,
     &              isym,igrid,nham1,nham2,nham3,npowr,
     &              bfieldf,rdefolt,qxin,btdes,qsurf,
     &              gamma,gamav,rmantl,betaf,zeffect,
     &              ncorr,corrfac,nccellr,peqpk0,peqpk1,peqpk2
      common/mapd/  mapmaxd,dpsisl,dpsisd,nqaxis,nwtmag,nfitmax,nfitpts,
     &              ifitrng,jfitrng,jfitchk,fitchek,cnvmag,epsaxs,
     &              maxerlp,maxerlv,delbox,delboz,delac,delav,
     &              delstsf,delstlp,delstlv,prfrac,
     &              nerstop,nerprnt,qptol,tolspln,tolbchi,tolbtor,
     &              tolsymm,tolaugm,errsep,precisn,plossmx,
     &              narcmx,ntrymx,ntdecr,ntmmin,npfit,npcmin,
     &              kuttaop,nrkmax0,nrkmax1,numstp,nwtfitp,nwtfitm,
     &              stepfac,flxstp,psispl,tolnwtp,tolnwtm,
     &              delpakf,delpakc,delpkf,delpkc,psichek,boxtnd,
     &              maptrace,norient,maxcutc,
     &              dresolv,dlclose,pntshft,endtol,
     &              narcln,nangax,nanglm,nbpmax,nwtmax,nslmax,nhfmax,
     &              bperor,sersnm,sertnm,arcmin,delgap,stepcut,
     &              novap,ioutm,iplotm
      common/size/  xdim,zdim,redge,zlowr,ipestg
      common/prof/  nprofl,npst,nfast,nspbc0,nspbc1,
     &              psimsh(npp),sf(npp),sp(npp),
     &              sffp(npp),spp(npp),sfp(npp),sdns(npp),
     &              spfst(npp),spsif(npb),sfast(npb),bcb(4),
     &              csf(3,npp),csp(3,npp),csffp(3,npp),cspp(3,npp),
     &              csfp(3,npp),csdn(3,npp),csfst(3,npp)
      common/eqd1/  nx,nz,nxd,nzd,dmx,dmz,
     &              x(nxx),z(nxz),psarray(nxx,nxz),
     &              gpx(nxx,nxz),gpy(nxx,nxz),grsq(nxx,nxz),
     &              cspl1(2,nxx,nh2),work0(nxz,2,nxx),vork0(nvk0,2)
      common/sarc/  ntmax,ntmsh,npcf,npcb,npc,xp(nlx),zp(nlx),gsq(nlx),
     &              arc(nlx),tp(nlx),arcc(nlx),tpp(nlx),bcd(4),
     &              csx(3,nlx),csz(3,nlx),cseq1(3,nlx),cseq2(3,nlx),
     &              st1(nlx),st2(nlx),st3(nlx),csveq(3,nlx),
     &              sv0(nlx),sv1(nlx),sv2(nlx),sv3(nlx),sv4(nlx),
     &              sv5(nlx)
      common/dblt/  psisep,xsep,zsep,smap(2),smax(2),smaz(2),numax,
     &              psimx(2),xax(2),zax(2),fitax(14,2),fitsep(14)
      common/ort1/  ndim1,ndim2,xma,zma,rcnt,btor,totcur,
     &              psimax0,psilim0,delpsir0,psimax1,psilim1,delpsir1,
     &              psimax,psilim,delpsir,dpsidec,dpsisv,dpedg1,dpedg2,
     &              eaxe,taxe,qaxe,fqiaxe,fqmaxe,paxe,faxe,ppaxe,ffpaxe,
     &              dnaxe,pfaxe,qppaxe,dlnaxe,dpnaxe,dhnaxe,
     &              arcaxe,dnnorm,pfnorm,
     &              qlim,fqilim,fqmlim,plim,flim,pplim,ffplim,
     &              dnlim,pflim,xlim,zlim,
     &              qmer,btmer,gkappa,arcprev,chicmax,
     &              charcmx,chpstmx,chhammx,arcnrm0,pstnrm0,hamnrm0,
     &              elongax1,elongax2,elongaxs,elongaxp,elongtn(npp)
      common/ort2/  psival(np1),pp(np1),ffp(np1),p(np1),f(np1),
     &              dnsty(np1),pfast(np1),qp(np1),fqpi(np1),fqim(np1),
     &              qdel(np1),psinitl(np1),qpinitl(np1),
     &              qprime(np1),qpint(np1),fqint(np1),
     &              psivl0(np2),rst(np2),zst(np2),rstp(np2),zstp(np2),
     &              rsinitl(np2),zsinitl(np2),
     &              arcsurf(nc1),pestchi(nc1),hamlchi(nc1),
     &              alfarc (nc1),alfpst (nc1),alfham (nc1),
     &              rs(nc1,np1),zs(nc1,np1)
      common/ort3/  rsrf(nc1),zsrf(nc1),chic(nc1),chie(nc1),dydx(nc1),
     &              xjsc(nc1)
      common/ort4/  rcell  (np1,nc1),zcell  (np1,nc1),
     &              dpsdr  (np1,nc1),dpsdz  (np1,nc1),
     &              chiarcl(np1,nc1),chipest(np1,nc1),chihaml(np1,nc1),
     &              xjcarcl(np1,nc1),xjcpest(np1,nc1),xjchaml(np1,nc1),
     &              alfarcl(np1,nc1),alfpest(np1,nc1),alfhaml(np1,nc1),
     &              arcnorm(np1),    pstnorm(np1),    hamnorm(np1),
     &              dlnlval(np1),    dpnlval(np1),    dhnlval(np1)
      common/ort5/  f3 (np1,nc1),f4 (np1,nc1),f5 (np1,nc1),
     &              f7 (np1,nc1),f8 (np1,nc1),f9 (np1,nc1),
     &              f10(np1,nc1),f11(np1,nc1),f12(np1,nc1),
     &              f13(np1,nc1),f14(np1,nc1),f15(np1,nc1),
     &              f16(np1,nc1),f17(np1,nc1),f18(np1,nc1),
     &              f19(np1,nc1),f20(np1,nc1),f21(np1,nc1),
     &              f22(np1,nc1),f23(np1,nc1),f24(np1,nc1),
     &              f25(np1,nc1),f26(np1,nc1),f27(np1,nc1),
     &              f28(np1,nc1)
      common/ort6/  savgax,tavgax,sntax0,sntax1,sntax2,sntax3,sntax4,
     &              svint(nvn),savge(np1),tavge(np1),sint0(np1),
     &              sint1(np1),sint2(np1),sint3(np1),sint4(np1)
      common/auxmsh/rh2tot,rhotot,flxtot,voltot,
     &              anltflxp,dvrtflxp,anltvolp,dvrtvolp,
     &              analtflx,divrtflx,analtvol,divrtvol,
     &              psivl1 (np2),qp1    (np2),fqpi1  (np2),
     &              qprim1 (np2),sintp0 (np2),psivmsh(np2),
     &              psivalu(np2),psinorm(np2),psisqrt(np2),
     &              psirho2(np2),psirhov(np2),psinrho(np2),
     &              psitorf(np2),psintor(np2),psisqtf(np2),
     &              psivolm(np2),psinvol(np2),psisqvl(np2),
     &              dpsirh2(np2),dpsirho(np2),dpsinrh(np2),
     &              dpsitor(np2),dpsintf(np2),dpsisqt(np2),
     &              dpsivol(np2),dpsinvl(np2),dpsisqv(np2)
       common/ratnl/jqmin, jqmax, numbqs,nq1totl,
     &              rcentr,zcentr,rminor,zminor,
     &              asprat00,asprat01,asprat10,asprat11,
     &              pminvl,qminvl,qpminv,pmaxvl,qmaxvl,qpmaxv,
     &              psivlmnq,    psivnmnq,    psisqmnq,
     &              psir2mnq,    psirhmnq,    psinrmnq,
     &              psitfmnq,    psintmnq,    psistmnq,
     &              psivmmnq,    psinvmnq,    psisvmnq,
     &              dpsr2mnq,    dpsrhmnq,    dpsnrmnq,
     &              dpstfmnq,    dpsntmnq,    dpsstmnq,
     &              dpsvmmnq,    dpsnvmnq,    dpssvmnq,
     &              psivlmxq,    psivnmxq,    psisqmxq,
     &              psir2mxq,    psirhmxq,    psinrmxq,
     &              psitfmxq,    psintmxq,    psistmxq,
     &              psivmmxq,    psinvmxq,    psisvmxq,
     &              dpsr2mxq,    dpsrhmxq,    dpsnrmxq,
     &              dpstfmxq,    dpsntmxq,    dpsstmxq,
     &              dpsvmmxq,    dpsnvmxq,    dpssvmxq,
     &              lpsiq  (npq),psivlq (npq),
     &              qprimq (npq),qvalue (npq),psimshq(npq),
     &              psivalq(npq),psinrmq(npq),psisqrq(npq),
     &              psirh2q(npq),psirhoq(npq),psinrhq(npq),
     &              psitorq(npq),psintfq(npq),psisqtq(npq),
     &              psivolq(npq),psinvlq(npq),psisqvq(npq),
     &              dpsir2q(npq),dpsirhq(npq),dpsinrq(npq),
     &              dpsitfq(npq),dpsintq(npq),dpsistq(npq),
     &              dpsivlq(npq),dpsinvq(npq),dpsisvq(npq),
     &              shearps(npq),shearrh(npq),
     &              sheartf(npq),shearvl(npq),
     &              epslrh1 (npq),shearrh1(npq),shearfrh(npq),
     &              epslvl1 (npq),shearvl1(npq),shearfvl(npq)
      common/vcal/  btnew,btave,bpave,betat,betap,betax0,betax1,volme,
     &              vhalf,bavet(3,nxx),bavep(3,nxx),pvolm(3,nxx),
     &              betav(3,nxx)
      common/volm/  pvansh,vp0,pmantl,vpm(np1),apm(np1),vcurnt(np1)
      common/labels/lbchpv, lbchsp, lbchrh, lbchro,
     &              lbchtf, lbchst, lbchvl, lbchsv,
     &              pclab,  tclab,  labpsi,
     &              labpsiv,labpsin,labsqtp,
     &              labrhov,labrho2,labrhon,
     &              labtorf,labtorn,labsqtt,
     &              labvolm,labvoln,labsqtv,
     &              labdrh2,labdrho,labdtor,
     &              labdsqt,labdvol,labdsqv
      common/pldf/  x0min,x0max,y0min,y0max
      common/inout/ kutty,kuout,kuinp,kueql,kumap,kugrd,kurdm,kurdn,
     &              kufst,krdeqlb,krdmgta,krdngta,krdfast
      common/flnm/  filout,filinp,fileql,filmap,filgrd,filrdm,filrdn,
     &              filfst,filplt
c
c
c
c 1.0 Initialize plots
c
c 1.1 Set the mapping type
c
      kmap    = nmap
      ktype   = nmtype
c
      if(ktype .lt.  0  .or.  ktype .gt. 2) then
         call abortjob
     &        ('plotg   ',  1,   'Invalid equilibrium type parameter  ' 
     &        ,'ktype   ', ktype,    'kmap    ', kmap,     -1)
         return
      endif

c
c
c 1.2 Initialize plot data
c
      x0min   = 0.0
      x0max   = 1.0
      y0min   = 0.0
      y0max   = 1.0
c
c
c 1.3 Return if no plots required
c
      if(iplotm .le. 0) return
c
c
c 1.4 Open the TV80 plot file
c
      call ncarcgm(1,filplt)
c
c
c
c 2.0 Plot the input data
c
c 2.1 Plot the namelist
c
      if(iplotm .ge. 1) call plotnml
c
c
c 2.2 Plot the mapping grid
c
      if(iplotm .ge. 2) call plotmap
c
c
c 2.3 Plot the equilibrium errors
c
      if(iplotm .ge. 3) call ploterr(ktype)
c
c
c 2.4 Plot the text data
c
      if(iplotm .ge. 4) call plottxt
c
c
c
c 3.0 Plot the packing data
c
      if(iplotm .ge. 5  .and.  kmesh .gt. 0) call plotpak
c
c
c
c 4.0 Plot the final profiles
c
c 4.1 Plot the q, pressure, and density, flux surface integrals, and qprime
c
      if(iplotm .ge. 6) then
        call plotqpr
        call plotint
        call plotqpp
      endif
c
c
c 4.2 Plot the radial flux functions
c
      if(iplotm .ge. 7) then
        call plotrho
        call plotdrh
      endif
c
c
c 4.3 Plot the profile tables
c
      if(iplotm .ge. 8) call plottab
c
c
c 4.4 Plot the equilibrium profiles
c
      if(iplotm .ge. 9) call plotfpj(ktype)
c
c
c
c 5.0 Close the TV80 plot file
c
      call plote
c
c
c
c 6.0 Return and end
c
      return
      end
      subroutine plotnml
c
c -------------------------------------------------------------
c  Plot namelist
c -------------------------------------------------------------
c
      parameter (npx=200,ncx=2*npx,np1=npx+1,np2=npx+2,nc1=ncx+1)
      parameter (nxx=513,nxz=nxx)
      parameter (npp=929,ntt=2*npp-1,npb=npp)
      parameter (nlx=1441)
      parameter (npk=np1,nrtdm=np1,nptt=2401)
      parameter (nhd1=5,nhd2=3)
c
      parameter (intfsz0=37,intfsz1=37,intfsz2=47
     &          ,intfsz3=27,intfsz4= 1,intfsz5= 1
     &          ,intfsz6=32,intfsz7=41,intfsz8=29)
      parameter (nplmx=5)
c
      character*32  intfil0(intfsz0)
      character*32  intfil1(intfsz1)
      character*32  intfil2(intfsz2)
      character*32  intfil3(intfsz3)
      character*32  intfil4(intfsz4)
      character*32  intfil5(intfsz5)
      character*32  intfil6(intfsz6)
      character*32  intfil7(intfsz7)
      character*32  intfil8(intfsz8)
      character*64  string (nhd1)
      character*64  strdat
c
      character*8   headr
      character*8   version
      character*8   verold0,verold1,verold2,verold3,verold4,verold5,
     &              verold6
      character*16  sourcnam,sourcdat
      character*16  datestamp,timestamp
      character*16  filout,filinp,fileql,filmap,filgrd,filrdm,filrdn,
     &              filfst,filplt
c
      common/vsno/  verold0,verold1,verold2,verold3,verold4,verold5,
     &              verold6,sourcnam,sourcdat,version,headr(nhd1,nhd2)
      common/cnst/  pi,hlfpi,twopi,amu,boltzm,echarg,pmass,roundff,bigno
      common/inpd/  ntor,ncase,idnsty,ndnxp0,ndnxp1,ndnxp2,ifastp,
     &              nmap,neqtyp,mpreset,ndoublt,ndivert,nmtype,
     &              jpsi,itht,jpsi1,jpsi2,ithtp,jsep,
     &              isym,igrid,nham1,nham2,nham3,npowr,
     &              bfieldf,rdefolt,qxin,btdes,qsurf,
     &              gamma,gamav,rmantl,betaf,zeffect,
     &              ncorr,corrfac,nccellr,peqpk0,peqpk1,peqpk2
      common/mapd/  mapmaxd,dpsisl,dpsisd,nqaxis,nwtmag,nfitmax,nfitpts,
     &              ifitrng,jfitrng,jfitchk,fitchek,cnvmag,epsaxs,
     &              maxerlp,maxerlv,delbox,delboz,delac,delav,
     &              delstsf,delstlp,delstlv,prfrac,
     &              nerstop,nerprnt,qptol,tolspln,tolbchi,tolbtor,
     &              tolsymm,tolaugm,errsep,precisn,plossmx,
     &              narcmx,ntrymx,ntdecr,ntmmin,npfit,npcmin,
     &              kuttaop,nrkmax0,nrkmax1,numstp,nwtfitp,nwtfitm,
     &              stepfac,flxstp,psispl,tolnwtp,tolnwtm,
     &              delpakf,delpakc,delpkf,delpkc,psichek,boxtnd,
     &              maptrace,norient,maxcutc,
     &              dresolv,dlclose,pntshft,endtol,
     &              narcln,nangax,nanglm,nbpmax,nwtmax,nslmax,nhfmax,
     &              bperor,sersnm,sertnm,arcmin,delgap,stepcut,
     &              novap,ioutm,iplotm
      common/mshc/  nmesh,npak,mpak,nedge,npkmax,nrat,mmesh,
     &              nrepeat,nppack,nqpack,nsrcheg,ncutedg,
     &              minpak,maxpak,incpak,psipak,chiwth,cspak,psincr,
     &              pkfrac,qpfrac,epsrat,sedg0,sedg1,epspak,
     &              spakmn,swidmn,swidmx,plpak(3,npk),pspak(3,npk)
      common/naml/  norm,nmod,nlt,nlim,toldrdz,
     &              iwal,iwalsym,nwall,
     &              irext,norign,nekdefn,maxitek,
     &              rext,rexmax,rcutoff,
     &              nrestrt,ndskopt,ndsktim,ndsksz,buffact,
     &              nev,neigmax,nforce,nreslv,nbrmax,
     &              nismax,ncymax,nitmax,ncyfin,mxdcomp,
     &              al0,dal0,al0bas,al0min,al0max,epschy,epscon,
     &              lineplt,lampplt,njplot,niplot,nskpj,nskpi,njedge,
     &              ntphase,ncont,ncplot,mshpsi,mshchi,nxisgn,
     &              nxiplt,nxuplt,nxrplt,nxpplt,nxdplt,ncphip,
     &              nbiplt,nbuplt,nbrplt,nbpplt,
     &              naiplt,nauplt,narplt,napplt,
     &              ndpplt,njiplt,njuplt,njrplt,njpplt,nvfft,
     &              torphase,dpltfac,dsplnrm,dspldbv,dspldav,
     &              psiscal,plsuprs,
     &              iomshp,ioeqlp,iowalp,
     &              ioeigp,iodbvp,iodavp,iodjvp,iopsip,
     &              iolinp,iolnbp,iolnap,iolnjp,
     &              iofftp,ioffbp,ioffap,ioffjp,
     &              ioconp,iodlbp,iodlap,iodljp,iodlwp,
     &              ioutw,iouta,ioute,ioutt,ioutp
      common/prof/  nprofl,npst,nfast,nspbc0,nspbc1,
     &              psimsh(npp),sf(npp),sp(npp),
     &              sffp(npp),spp(npp),sfp(npp),sdns(npp),
     &              spfst(npp),spsif(npb),sfast(npb),bcb(4),
     &              csf(3,npp),csp(3,npp),csffp(3,npp),cspp(3,npp),
     &              csfp(3,npp),csdn(3,npp),csfst(3,npp)
      common/vcal/  btnew,btave,bpave,betat,betap,betax0,betax1,volme,
     &              vhalf,bavet(3,nxx),bavep(3,nxx),pvolm(3,nxx),
     &              betav(3,nxx)
      common/volm/  pvansh,vp0,pmantl,vpm(np1),apm(np1),vcurnt(np1)
      common/geom/  aminor,rcgeom,epslon,capa,triangl,
     &              allim,bpornl,deltal,qcyl,q0surf,welln,
     &              well(np1),vdpsi(np1)
      common/psft/  rpsi,zpsi,rj0,rjp,rj3,rj5,rj7,sa0,sa2,sa3,sa4,
     &              sb2,sb3,sb4
      common/pldf/  x0min,x0max,y0min,y0max
      common/tcpu/  timcpu,timiop,timsys,timwal,datestamp,timestamp,
     &              tim0c,tim0i,tim0s,tim0w,tim1c,tim1i,tim1s,tim1w
      common/inout/ kutty,kuout,kuinp,kueql,kumap,kugrd,kurdm,kurdn,
     &              kufst,krdeqlb,krdmgta,krdngta,krdfast
      common/flnm/  filout,filinp,fileql,filmap,filgrd,filrdm,filrdn,
     &              filfst,filplt
c
c
      namelist/inputd/ntor,ncase,norm,nmod,nlt,nlim,idnsty,ifastp
     &               ,ndnxp0,ndnxp1,ndnxp2,bfieldf,rdefolt,qxin,btdes
     &               ,qsurf,gamma,rmantl,betaf,zeffect
     &               ,nmap,neqtyp,mpreset,ndoublt,ndivert,ncorr,corrfac
     &               ,jpsi,itht,isym,igrid,nham1,nham2,nham3
     &               ,nmesh,npak,mpak,nedge,npkmax,nrat,nrepeat,nppack
     &               ,nqpack,nsrcheg,ncutedg,minpak,maxpak,incpak
     &               ,psipak,chiwth,cspak,psincr,pkfrac,qpfrac,epsrat
     &               ,sedg0,sedg1,plpak,pspak,epspak,spakmn
     &               ,swidmn,swidmx
     &               ,mapmaxd,dpsisl,dpsisd
     &               ,nqaxis,nspbc0,nspbc1,nwtmag,nfitmax,nfitpts
     &               ,ifitrng,jfitrng,jfitchk,fitchek,cnvmag,epsaxs
     &               ,maxerlp,maxerlv,delbox,delboz,delac,delav
     &               ,delstsf,delstlp,delstlv,prfrac
     &               ,nerstop,nerprnt,qptol,tolspln,tolbchi
     &               ,tolbtor,tolsymm,tolaugm,errsep
     &               ,toldrdz,pvansh,precisn,plossmx,roundff,bigno
     &               ,narcmx,ntrymx,ntdecr,ntmmin
     &               ,nccellr,peqpk0,peqpk1,peqpk2,npfit,npcmin,kuttaop
     &               ,nrkmax0,nrkmax1,numstp,nwtfitp,nwtfitm
     &               ,stepfac,flxstp,psispl,tolnwtp,tolnwtm
     &               ,delpakf,delpakc,delpkf,delpkc,psichek,boxtnd
     &               ,maptrace,norient,maxcutc
     &               ,dresolv,dlclose,pntshft,endtol
     &               ,narcln,nangax,nanglm,nbpmax,nwtmax,nslmax,nhfmax
     &               ,bperor,sersnm,sertnm,arcmin,delgap,stepcut
     &               ,iwal,iwalsym,nwall,irext,norign,nekdefn,maxitek
     &               ,rext,rexmax,rcutoff
     &               ,nrestrt,ndskopt,ndsktim,ndsksz,buffact
     &               ,nev,neigmax,nforce,nreslv,nbrmax
     &               ,nismax,ncymax,nitmax,ncyfin,mxdcomp
     &               ,al0,dal0,al0bas,al0min,al0max,epschy,epscon
     &               ,lineplt,lampplt,njplot,niplot,nskpj,nskpi,njedge
     &               ,ntphase,npowr,ncont,ncplot,mshpsi,mshchi,nxisgn
     &               ,nxiplt,nxuplt,nxrplt,nxpplt,nxdplt,ncphip
     &               ,nbiplt,nbuplt,nbrplt,nbpplt
     &               ,naiplt,nauplt,narplt,napplt
     &               ,ndpplt,njiplt,njuplt,njrplt,njpplt,nvfft
     &               ,torphase,dpltfac,dsplnrm,dspldbv,dspldav
     &               ,psiscal,plsuprs
     &               ,iomshp,ioeqlp,iowalp
     &               ,ioeigp,iodbvp,iodavp,iodjvp,iopsip
     &               ,iolinp,iolnbp,iolnap,iolnjp
     &               ,iofftp,ioffbp,ioffap,ioffjp
     &               ,ioconp,iodlbp,iodlap,iodljp,iodlwp
     &               ,novap,iplotm,ioutm,ioutw,iouta,ioute,ioutt
     &               ,ioutp
c
c
c
c 1.0 Initialize the title page and print out the heading
c
c 1.1 Initialize the title page
c
c 1.1.1 Set the page limits
c
      x0mn0    = x0min
      x0mx0    = x0max
      y0mn0    = y0min
      y0mx0    = y0max
c
c 1.1.2 Initialize the page
c
      call map(x0mn0,x0mx0,y0mn0,y0mx0,x0mn0,x0mx0,y0mn0,y0mx0)
c
c
c 1.2 Set up and plot the heading
c
c 1.2.1 Set up the pointer position and print characteristics
c
      xhead    = 0.0325
      yhead    = 0.7250
      icshd    = 2
      iszhd    = 3
      iorhd    = 0
c
c
c 1.2.2 Set up the box positions
c
      hboxx0   = 0.025
      hboxy0   = 0.500
      hboxx1   = 0.975
      hboxy1   = 0.800
c
c
c 1.2.3 Plot the heading
c
      call setlch(xhead,yhead,icshd,iszhd,iorhd,-1)
c
      do 100 k2  = 1,nhd2
      write(string,1000) (headr(k1,k2),k1=1,nhd1)
      call wrtstr(string,1)
  100 continue
c
c
c 1.2.4 Plot the box around the heading
c
      call line(hboxx0,hboxy0,hboxx0,hboxy1)
      call line(hboxx0,hboxy1,hboxx1,hboxy1)
      call line(hboxx1,hboxy1,hboxx1,hboxy0)
      call line(hboxx1,hboxy0,hboxx0,hboxy0)
c
c
c
c 2.0 Plot the Version number and date and time stamps.
c
c 2.1 Plot the version number
c
c 2.1.1 Set up the pointer and print characteristics
c
      xvers0   = 0.0500
      yvers0   = 0.2000
      icsvs0   = 2
      iszvs0   = 2
      iorvs0   = 0
c
      call setlch(xvers0,yvers0,icsvs0,iszvs0,iorvs0,-1)
c
c 2.1.2 Plot the Source date and Version number
c
      write(strdat,1500) sourcnam,sourcdat,version
      call wrtstr(strdat,1)
c
c
c 2.2 Plot the date and time stamps
c
c 2.2.1 Set up the pointer and print characteristics
c
      xvers1   = 0.0525
      yvers1   = 0.1500
      icsvs1   = 2
      iszvs1   = 2
      iorvs1   = 0
c
      call setlch(xvers1,yvers1,icsvs1,iszvs1,iorvs1,-1)
c
c 2.2.2 Plot the date and time stamps
c
      write(strdat,1510) datestamp,timestamp
      call wrtstr(strdat,1)
c
c
c 2.3 Close the frame
c
      call frame(0)
c
c
c
c 3.0 Plot the first page of the namelist
c
c 3.1 Reinitialize the page
c
c 3.1.1 Set the page limits
c
      x0mn1    = x0mn0
      x0mx1    = x0mx0
      y0mn1    = y0mn0
      y0mx1    = y0mx0
c
c 3.1.2 Initialize the page
c
      call map(x0mn1,x0mx1,y0mn1,y0mx1,x0mn1,x0mx1,y0mn1,y0mx1)
c
c
c 3.2 Set up and plot the page heading
c
c 3.2.1 Set up the pointer position and print characteristics
c
      xnamhd   = 0.0750
      ynamhd   = 0.9375
      icsnmh   = 2
      isznmh   = 3
      iornmh   = 0
c
      call setlch(xnamhd,ynamhd,icsnmh,isznmh,iornmh,-1)
c
c 3.2.2 Set up the box positions
c
      tboxx0   = 0.0500
      tboxy0   = 0.9125
      tboxx1   = 0.9000
      tboxy1   = 0.9950
c
c 3.2.3 Print the page title
c
      write(strdat,2000)
      call wrtstr(strdat,1)
c
c 3.2.4 Plot the box around the heading
c
      call line(tboxx0,tboxy0,tboxx0,tboxy1)
      call line(tboxx0,tboxy1,tboxx1,tboxy1)
      call line(tboxx1,tboxy1,tboxx1,tboxy0)
      call line(tboxx1,tboxy0,tboxx0,tboxy0)
c
c
c 3.3 Define the starting points and print specifications
c
c 3.3.1 Specify the starting points
c
      xnaml0   = 0.0150
      ynaml0   = 0.8750
      xnaml1   = 0.4000
      ynaml1   = ynaml0
      xnaml2   = 0.7500
      ynaml2   = ynaml0
c
c 3.3.2 Set the print specifications
c
      icsnm    = 2
      isznm    = 0
      iornm    = 0
c
c
c 3.4 Print the mapping namelist data
c
c 3.4.1 Print the physical data, dimension data, and general mapping
c       parameters
c
      call setlch(xnaml0,ynaml0,icsnm,isznm,iornm,-1)
      write(intfil0,2100)
     &                ntor,ncase,norm,nmod,nlt,nlim,idnsty,ifastp
     &               ,ndnxp0,ndnxp1,ndnxp2,bfieldf,rdefolt,qxin,btdes
     &               ,qsurf,gamma,rmantl,betaf,zeffect
     &               ,nmap,neqtyp,mpreset,ndoublt,ndivert,ncorr,corrfac
     &               ,jpsi,itht,isym,igrid,nham1,nham2,nham3
     &               ,mapmaxd,dpsisl,dpsisd
c
      call wrtstr(intfil0,intfsz0)
c
c 3.4.2 Print the axis fitting data, equilibrium checking data and
c       tolerance data
c
      call setlch(xnaml1,ynaml1,icsnm,isznm,iornm,-1)
      write(intfil1,2200)
     &                nqaxis,nspbc0,nspbc1,nwtmag,nfitmax,nfitpts
     &               ,ifitrng,jfitrng,jfitchk,fitchek,cnvmag,epsaxs
     &               ,maxerlp,maxerlv,delbox,delboz,delac,delav
     &               ,delstsf,delstlp,delstlv,prfrac
     &               ,nerstop,nerprnt,qptol,tolspln,tolbchi
     &               ,tolbtor,tolsymm,tolaugm,errsep
     &               ,toldrdz,pvansh,precisn,plossmx,roundff,bigno
      call wrtstr(intfil1,intfsz1)
c
c 3.4.3 Print the flux surface mapping data
c
      call setlch(xnaml2,ynaml2,icsnm,isznm,iornm,-1)
      write(intfil2,2300)
     &                narcmx,ntrymx,ntdecr,ntmmin
     &               ,nccellr,peqpk0,peqpk1,peqpk2,npfit,npcmin,kuttaop
     &               ,nrkmax0,nrkmax1,numstp,nwtfitp,nwtfitm
     &               ,stepfac,flxstp,psispl,tolnwtp,tolnwtm
     &               ,delpakf,delpakc,delpkf,delpkc,psichek,boxtnd
     &               ,maptrace,norient,maxcutc
     &               ,dresolv,dlclose,pntshft,endtol
     &               ,narcln,nangax,nanglm,nbpmax,nwtmax,nslmax,nhfmax
     &               ,bperor,sersnm,sertnm,arcmin,delgap,stepcut
      call wrtstr(intfil2,intfsz2)
c
c
c 3.5 Close the frame
c
      call frame(0)
c
c
c
c 4.0 Plot the second page of the namelist
c
c 4.1 Reinitialize the page
c
c 4.1.1 Set the page limits
c
      x0mn1    = x0mn0
      x0mx1    = x0mx0
      y0mn1    = y0mn0
      y0mx1    = y0mx0
c
c 4.1.2 Initialize the page
c
      call map(x0mn1,x0mx1,y0mn1,y0mx1,x0mn1,x0mx1,y0mn1,y0mx1)
c
c
c 4.2 Set up and plot the page heading
c
c 4.2.1 Set up the pointer position and print characteristics
c
      xnamhd   = 0.0750
      ynamhd   = 0.9375
      icsnmh   = 2
      isznmh   = 3
      iornmh   = 0
c
      call setlch(xnamhd,ynamhd,icsnmh,isznmh,iornmh,-1)
c
c 4.2.2 Set up the box positions
c
      tboxx0   = 0.0250
      tboxy0   = 0.9125
      tboxx1   = 0.9625
      tboxy1   = 0.9950
c
c 4.2.3 Print the page title
c
      write(strdat,3000)
      call wrtstr(strdat,1)
c
c 4.2.4 Plot the box around the heading
c
      call line(tboxx0,tboxy0,tboxx0,tboxy1)
      call line(tboxx0,tboxy1,tboxx1,tboxy1)
      call line(tboxx1,tboxy1,tboxx1,tboxy0)
      call line(tboxx1,tboxy0,tboxx0,tboxy0)
c
c
c 4.3 Define the starting points and print specifications
c
c 4.3.1 Specify the starting points
c
      xnaml3   = 0.1500
      ynaml3   = 0.8750
      xnaml4   = 0.4000
      ynaml4   = ynaml3
      xnaml5   = 0.7500
      ynaml5   = ynaml3
c
c 4.3.2 Set the print specifications
c
      icsnm    = 2
      isznm    = 0
      iornm    = 0
c
c
c 4.4 Print the flux surface packing data
c
c 4.4.1 Scalar packing parameters
c
      call setlch(xnaml3,ynaml3,icsnm,isznm,iornm,-1)
      write(intfil3,3100)
     &                nmesh,npak,mpak,nedge,npkmax,nrat,nrepeat,nppack
     &               ,nqpack,nsrcheg,ncutedg,minpak,maxpak,incpak
     &               ,psipak,chiwth,cspak,psincr,pkfrac,qpfrac,epsrat
     &               ,sedg0,sedg1,epspak,spakmn
     &               ,swidmn,swidmx
      call wrtstr(intfil3,intfsz3)
c
c 4.4.2 Print the q packing positions, widths, and amplitudes
c
      npakp    = iabs(npak)
      npkpr    = min0(npakp,npk,nplmx)
      if(npkpr .gt. 0) then
c
c 4.4.2.1 Set starting position
        call setlch(xnaml4,ynaml4,icsnm,isznm,iornm,-1)
c
c 4.4.2.2 Loop over packing positions
        do 200 k = 1,npkpr
c
c 4.4.2.3 Print packing position
        write(intfil4,3200) k,plpak(1,k)
        call wrtstr(intfil4,intfsz4)
c
c 4.4.2.4 Print packing width
        write(intfil4,3210) k,plpak(2,k)
        call wrtstr(intfil4,intfsz4)
c
c 4.4.2.5 Print packing amplitude
        write(intfil4,3220) k,plpak(3,k)
        call wrtstr(intfil4,intfsz4)
  200   continue
      endif
c
c 4.4.3 Print the psi packing positions, widths, and amplitudes
c
      mpakp    = iabs(mpak)
      mpkpr    = min0(mpakp,npk,nplmx)
      if(mpkpr .gt. 0) then
c
c 4.4.3.1 Set starting position
        call setlch(xnaml5,ynaml5,icsnm,isznm,iornm,-1)
c
c 4.4.3.2 Loop over packing positions
        do 250 k = 1,mpkpr
c
c 4.4.3.3 Print packing position
        write(intfil5,3300) k,pspak(1,k)
        call wrtstr(intfil5,intfsz5)
c
c 4.4.3.4 Print packing width
        write(intfil5,3310) k,pspak(2,k)
        call wrtstr(intfil5,intfsz5)
c
c 4.4.3.5 Print packing amplitude
        write(intfil5,3320) k,pspak(3,k)
        call wrtstr(intfil5,intfsz5)
  250   continue
      endif
c
c
c 4.5 Close the frame
c
      call frame(0)
c
c
c
c 5.0 Plot the third page of the namelist
c
c 5.1 Reinitialize the page
c
c 5.1.1 Set the page limits
c
      x0mn1    = x0mn0
      x0mx1    = x0mx0
      y0mn1    = y0mn0
      y0mx1    = y0mx0
c
c 5.1.2 Initialize the page
c
      call map(x0mn1,x0mx1,y0mn1,y0mx1,x0mn1,x0mx1,y0mn1,y0mx1)
c
c
c 5.2 Set up and plot the page heading
c
c 5.2.1 Set up the pointer position and print characteristics
c
      xnamhd   = 0.0750
      ynamhd   = 0.9375
      icsnmh   = 2
      isznmh   = 3
      iornmh   = 0
c
      call setlch(xnamhd,ynamhd,icsnmh,isznmh,iornmh,-1)
c
c 5.2.2 Set up the box positions
c
      tboxx0   = 0.0250
      tboxy0   = 0.9125
      tboxx1   = 0.9625
      tboxy1   = 0.9950
c
c 5.2.3 Print the page title
c
      write(strdat,4000)
      call wrtstr(strdat,1)
c
c 5.2.4 Plot the box around the heading
c
      call line(tboxx0,tboxy0,tboxx0,tboxy1)
      call line(tboxx0,tboxy1,tboxx1,tboxy1)
      call line(tboxx1,tboxy1,tboxx1,tboxy0)
      call line(tboxx1,tboxy0,tboxx0,tboxy0)
c
c
c 5.3 Define the starting points and print specifications
c
c 5.3.1 Specify the starting points
c
      xnaml6   = 0.0150
      ynaml6   = 0.7875
      xnaml7   = 0.4000
      ynaml7   = ynaml6
      xnaml8   = 0.7500
      ynaml8   = ynaml6
c
c 5.3.2 Set the print specifications
c
      icsnm    = 2
      isznm    = 0
      iornm    = 0
c
c
c 5.4 Print the wall vacuum and eigenvalue search data
c
      call setlch(xnaml6,ynaml6,icsnm,isznm,iornm,-1)
      write(intfil6,4100)
     &                iwal,iwalsym,nwall
     &               ,irext,norign,nekdefn,maxitek
     &               ,rext,rexmax,rcutoff
     &               ,nrestrt,ndskopt,ndsktim,ndsksz,buffact
     &               ,nev,neigmax,nforce,nreslv,nbrmax
     &               ,nismax,ncymax,nitmax,ncyfin,mxdcomp
     &               ,al0,dal0,al0bas,al0min,al0max,epschy,epscon
      call wrtstr(intfil6,intfsz6)
c
c
c 5.5 Print the plotting specification data
c
      call setlch(xnaml7,ynaml7,icsnm,isznm,iornm,-1)
      write(intfil7,4200)
     &                lineplt,lampplt,njplot,niplot,nskpj,nskpi,njedge
     &               ,ntphase,npowr,ncont,ncplot,mshpsi,mshchi,nxisgn
     &               ,nxiplt,nxuplt,nxrplt,nxpplt,nxdplt,ncphip
     &               ,nbiplt,nbuplt,nbrplt,nbpplt
     &               ,naiplt,nauplt,narplt,napplt
     &               ,ndpplt,njiplt,njuplt,njrplt,njpplt,nvfft
     &               ,torphase,dpltfac,dsplnrm,dspldbv,dspldav
     &               ,psiscal,plsuprs
      call wrtstr(intfil7,intfsz7)
c
c
c 5.5 Print the plotting specification data
c
      call setlch(xnaml8,ynaml8,icsnm,isznm,iornm,-1)
      write(intfil8,4300)
     &                iomshp,ioeqlp,iowalp
     &               ,ioeigp,iodbvp,iodavp,iodjvp,iopsip
     &               ,iolinp,iolnbp,iolnap,iolnjp
     &               ,iofftp,ioffbp,ioffap,ioffjp
     &               ,ioconp,iodlbp,iodlap,iodljp,iodlwp
     &               ,novap,iplotm,ioutm,ioutw,iouta,ioute,ioutt
     &               ,ioutp
      call wrtstr(intfil8,intfsz8)
c
c
c
c 6.0 Close the frame, return and end
c
      call frame(0)
      return
c
 1000 format(5a8)
 1500 format(a8,2x,'Source Date: ',a16,2x,'Version:  ', a8)
 1510 format(10x,  'Run Date:    ',a10,8x,'Run Time: ',a10)
 2000 format('Physical and Mapping Namelist Input')
 2100 format(
     &       ' ntor     = ',i5,/,
     &       ' ncase    = ',i5,/,
     &       ' norm     = ',i5,/,
     &       ' nmod     = ',i5,/,
     &       ' nlt      = ',i5,/,
     &       ' nlim     = ',i5,/,
     &       ' idnsty   = ',i5,/,
     &       ' ifastp   = ',i5,/,
     &       ' ndnxp0   = ',i5,/,
     &       ' ndnxp1   = ',i5,/,
     &       ' ndnxp2   = ',i5,/,
     &       ' bfieldf  = ',e15.8,/,
     &       ' rdefolt  = ',e15.8,/,
     &       ' qxin     = ',e15.8,/,
     &       ' btdes    = ',e15.8,/,
     &       ' qsurf    = ',e15.8,/,
     &       ' gamma    = ',e15.8,/,
     &       ' rmantl   = ',e15.8,/,
     &       ' betaf    = ',e15.8,/,
     &       ' zeffect  = ',e15.8,/,
     &       ' nmap     = ',i5,/,
     &       ' neqtyp   = ',i5,/,
     &       ' mpreset  = ',i5,/,
     &       ' ndoublt  = ',i5,/,
     &       ' ndivert  = ',i5,/,
     &       ' ncorr    = ',i5,/,
     &       ' corrfac  = ',e15.8,/,
     &       ' jpsi     = ',i5,/,
     &       ' itht     = ',i5,/,
     &       ' isym     = ',i5,/,
     &       ' igrid    = ',i5,/,
     &       ' nham1    = ',i5,/,
     &       ' nham2    = ',i5,/,
     &       ' nham3    = ',i5,/,
     &       ' mapmaxd  = ',i5,/,
     &       ' dpsisl   = ',e15.8,/,
     &       ' dpsisd   = ',e15.8)
 2200 format(
     &       ' nqaxis   = ',i5,/,
     &       ' nspbc0   = ',i5,/,
     &       ' nspbc1   = ',i5,/,
     &       ' nwtmag   = ',i5,/,
     &       ' nfitmax  = ',i5,/,
     &       ' nfitpts  = ',i5,/,
     &       ' ifitrng  = ',i5,/,
     &       ' jfitrng  = ',i5,/,
     &       ' jfitchk  = ',i5,/,
     &       ' fitchek  = ',e15.8,/,
     &       ' cnvmag   = ',e15.8,/,
     &       ' epsaxs   = ',e15.8,/,
     &       ' maxerlp  = ',i5,/,
     &       ' maxerlv  = ',i5,/,
     &       ' delbox   = ',e15.8,/,
     &       ' delboz   = ',e15.8,/,
     &       ' delac    = ',e15.8,/,
     &       ' delav    = ',e15.8,/,
     &       ' delstsf  = ',e15.8,/,
     &       ' delstlp  = ',e15.8,/,
     &       ' delstlv  = ',e15.8,/,
     &       ' prfrac   = ',e15.8,/,
     &       ' nerstop  = ',i5,/,
     &       ' nerprnt  = ',i5,/,
     &       ' qptol    = ',e15.8,/,
     &       ' tolspln  = ',e15.8,/,
     &       ' tolbchi  = ',e15.8,/,
     &       ' tolbtor  = ',e15.8,/,
     &       ' tolsymm  = ',e15.8,/,
     &       ' tolaugm  = ',e15.8,/,
     &       ' errsep   = ',e15.8,/,
     &       ' toldrdz  = ',e15.8,/,
     &       ' pvansh   = ',e15.8,/,
     &       ' precisn  = ',e15.8,/,
     &       ' plossmx  = ',e15.8,/,
     &       ' roundff  = ',e15.8,/,
     &       ' bigno    = ',e15.8)
 2300 format(
     &       ' narcmx   = ',i5,/,
     &       ' ntrymx   = ',i5,/,
     &       ' ntdecr   = ',i5,/,
     &       ' ntmmin   = ',i5,/,
     &       ' nccellr  = ',i5,/,
     &       ' peqpk0   = ',e15.8,/,
     &       ' peqpk1   = ',e15.8,/,
     &       ' peqpk2   = ',e15.8,/,
     &       ' npfit    = ',i5,/,
     &       ' npcmin   = ',i5,/,
     &       ' kuttaop  = ',i5,/,
     &       ' nrkmax0  = ',i5,/,
     &       ' nrkmax1  = ',i5,/,
     &       ' numstp   = ',i5,/,
     &       ' nwtfitp  = ',i5,/,
     &       ' nwtfitm  = ',i5,/,
     &       ' stepfac  = ',e15.8,/,
     &       ' flxstp   = ',e15.8,/,
     &       ' psispl   = ',e15.8,/,
     &       ' tolnwtp  = ',e15.8,/,
     &       ' tolnwtm  = ',e15.8,/,
     &       ' delpakf  = ',e15.8,/,
     &       ' delpakc  = ',e15.8,/,
     &       ' delpkf   = ',e15.8,/,
     &       ' delpkc   = ',e15.8,/,
     &       ' psichek  = ',e15.8,/,
     &       ' boxtnd   = ',e15.8,/,
     &       ' maptrace = ',i5,/,
     &       ' norient  = ',i5,/,
     &       ' maxcutc  = ',i5,/,
     &       ' dresolv  = ',e15.8,/,
     &       ' dlclose  = ',e15.8,/,
     &       ' pntshft  = ',e15.8,/,
     &       ' endtol   = ',e15.8,/,
     &       ' narcln   = ',i5,/,
     &       ' nangax   = ',i5,/,
     &       ' nanglm   = ',i5,/,
     &       ' nbpmax   = ',i5,/,
     &       ' nwtmax   = ',i5,/,
     &       ' nslmax   = ',i5,/,
     &       ' nhfmax   = ',i5,/,
     &       ' bperor   = ',e15.8,/,
     &       ' sersnm   = ',e15.8,/,
     &       ' sertnm   = ',e15.8,/,
     &       ' arcmin   = ',e15.8,/,
     &       ' delgap   = ',e15.8,/,
     &       ' stepcut  = ',e15.8)
 3000 format('Mesh Packing data')
 3100 format(
     &       ' nmesh    = ',i5,/,
     &       ' npak     = ',i5,/,
     &       ' mpak     = ',i5,/,
     &       ' nedge    = ',i5,/,
     &       ' npkmax   = ',i5,/,
     &       ' nrat     = ',i5,/,
     &       ' nrepeat  = ',i5,/,
     &       ' nppack   = ',i5,/,
     &       ' nqpack   = ',i5,/,
     &       ' nsrcheg  = ',i5,/,
     &       ' ncutedg  = ',i5,/,
     &       ' minpak   = ',i5,/,
     &       ' maxpak   = ',i5,/,
     &       ' incpak   = ',i5,/,
     &       ' psipak   = ',e15.8,/,
     &       ' chiwth   = ',e15.8,/,
     &       ' cspak    = ',e15.8,/,
     &       ' psincr   = ',e15.8,/,
     &       ' pkfrac   = ',e15.8,/,
     &       ' qpfrac   = ',e15.8,/,
     &       ' epsrat   = ',e15.8,/,
     &       ' sedg0    = ',e15.8,/,
     &       ' sedg1    = ',e15.8,/,
     &       ' epspak   = ',e15.8,/,
     &       ' spakmn   = ',e15.8,/,
     &       ' swidmn   = ',e15.8,/,
     &       ' swidmx   = ',e15.8)
 3200 format(' plpak(1,',i3,') = ',e12.5)
 3210 format(' plpak(2,',i3,') = ',e12.5)
 3220 format(' plpak(3,',i3,') = ',e12.5)
 3300 format(' pspak(1,',i3,') = ',e12.5)
 3310 format(' pspak(2,',i3,') = ',e12.5)
 3320 format(' pspak(3,',i3,') = ',e12.5)
 4000 format('Eigenvalue and Plotting Namelist Input')
 4100 format(
     &       ' iwal     = ',i5,/,
     &       ' iwalsym  = ',i5,/,
     &       ' nwall    = ',i5,/,
     &       ' irext    = ',i5,/,
     &       ' norign   = ',i5,/,
     &       ' nekdefn  = ',i5,/,
     &       ' maxitek  = ',i5,/,
     &       ' rext     = ',e15.8,/,
     &       ' rexmax   = ',e15.8,/,
     &       ' rcutoff  = ',e15.8,/,
     &       ' nrestrt  = ',i5,/,
     &       ' ndskopt  = ',i5,/,
     &       ' ndsktim  = ',i5,/,
     &       ' ndsksz   = ',i5,/,
     &       ' buffact  = ',e15.8,/,
     &       ' nev      = ',i5,/,
     &       ' neigmax  = ',i5,/,
     &       ' nforce   = ',i5,/,
     &       ' nreslv   = ',i5,/,
     &       ' nbrmax   = ',i5,/,
     &       ' nismax   = ',i5,/,
     &       ' ncymax   = ',i5,/,
     &       ' nitmax   = ',i5,/,
     &       ' ncyfin   = ',i5,/,
     &       ' mxdcomp  = ',i5,/,
     &       ' al0      = ',e15.8,/,
     &       ' dal0     = ',e15.8,/,
     &       ' al0bas   = ',e15.8,/,
     &       ' al0min   = ',e15.8,/,
     &       ' al0max   = ',e15.8,/,
     &       ' epschy   = ',e15.8,/,
     &       ' epscon   = ',e15.8)
 4200 format(
     &       ' lineplt  = ',i5,/,
     &       ' lampplt  = ',i5,/,
     &       ' njplot   = ',i5,/,
     &       ' niplot   = ',i5,/,
     &       ' nskpj    = ',i5,/,
     &       ' nskpi    = ',i5,/,
     &       ' njedge   = ',i5,/,
     &       ' ntphase  = ',i5,/,
     &       ' npowr    = ',i5,/,
     &       ' ncont    = ',i5,/,
     &       ' ncplot   = ',i5,/,
     &       ' mshpsi   = ',i5,/,
     &       ' mshchi   = ',i5,/,
     &       ' nxisgn   = ',i5,/,
     &       ' nxiplt   = ',i5,/,
     &       ' nxuplt   = ',i5,/,
     &       ' nxrplt   = ',i5,/,
     &       ' nxpplt   = ',i5,/,
     &       ' nxdplt   = ',i5,/,
     &       ' ncphip   = ',i5,/,
     &       ' nbiplt   = ',i5,/,
     &       ' nbuplt   = ',i5,/,
     &       ' nbrplt   = ',i5,/,
     &       ' nbpplt   = ',i5,/,
     &       ' naiplt   = ',i5,/,
     &       ' nauplt   = ',i5,/,
     &       ' narplt   = ',i5,/,
     &       ' napplt   = ',i5,/,
     &       ' ndpplt   = ',i5,/,
     &       ' njiplt   = ',i5,/,
     &       ' njuplt   = ',i5,/,
     &       ' njrplt   = ',i5,/,
     &       ' njpplt   = ',i5,/,
     &       ' nvfft    = ',i5,/,
     &       ' torphase = ',e15.8,/,
     &       ' dpltfac  = ',e15.8,/,
     &       ' dsplnrm  = ',e15.8,/,
     &       ' dspldbv  = ',e15.8,/,
     &       ' dspldav  = ',e15.8,/,
     &       ' psiscal  = ',e15.8,/,
     &       ' plsuprs  = ',e15.8)
 4300 format(
     &       ' iomshp   = ',i5,/,
     &       ' ioeqlp   = ',i5,/,
     &       ' iowalp   = ',i5,/,
     &       ' ioeigp   = ',i5,/,
     &       ' iodbvp   = ',i5,/,
     &       ' iodavp   = ',i5,/,
     &       ' iodjvp   = ',i5,/,
     &       ' iopsip   = ',i5,/,
     &       ' iolinp   = ',i5,/,
     &       ' iolnbp   = ',i5,/,
     &       ' iolnap   = ',i5,/,
     &       ' iolnjp   = ',i5,/,
     &       ' iofftp   = ',i5,/,
     &       ' ioffbp   = ',i5,/,
     &       ' ioffap   = ',i5,/,
     &       ' ioffjp   = ',i5,/,
     &       ' ioconp   = ',i5,/,
     &       ' iodlbp   = ',i5,/,
     &       ' iodlap   = ',i5,/,
     &       ' iodljp   = ',i5,/,
     &       ' iodlwp   = ',i5,/,
     &       ' novap    = ',i5,/,
     &       ' iplotm   = ',i5,/,
     &       ' ioutm    = ',i5,/,
     &       ' ioutw    = ',i5,/,
     &       ' iouta    = ',i5,/,
     &       ' ioute    = ',i5,/,
     &       ' ioutt    = ',i5,/,
     &       ' ioutp    = ',i5)
      end
      subroutine plotmap
c
c ----------------------
c plot equilibrium mesh
c ----------------------
c
      parameter (npx=200,ncx=2*npx,np1=npx+1,np2=npx+2,nc1=ncx+1)
      parameter (npp=929,ntt=2*npp-1,npb=npp)
      parameter (nlx=1441)
      parameter (nvn=7)
c
c
      common/cnst/  pi,hlfpi,twopi,amu,boltzm,echarg,pmass,roundff,bigno
      common/inpd/  ntor,ncase,idnsty,ndnxp0,ndnxp1,ndnxp2,ifastp,
     &              nmap,neqtyp,mpreset,ndoublt,ndivert,nmtype,
     &              jpsi,itht,jpsi1,jpsi2,ithtp,jsep,
     &              isym,igrid,nham1,nham2,nham3,npowr,
     &              bfieldf,rdefolt,qxin,btdes,qsurf,
     &              gamma,gamav,rmantl,betaf,zeffect,
     &              ncorr,corrfac,nccellr,peqpk0,peqpk1,peqpk2
      common/mapd/  mapmaxd,dpsisl,dpsisd,nqaxis,nwtmag,nfitmax,nfitpts,
     &              ifitrng,jfitrng,jfitchk,fitchek,cnvmag,epsaxs,
     &              maxerlp,maxerlv,delbox,delboz,delac,delav,
     &              delstsf,delstlp,delstlv,prfrac,
     &              nerstop,nerprnt,qptol,tolspln,tolbchi,tolbtor,
     &              tolsymm,tolaugm,errsep,precisn,plossmx,
     &              narcmx,ntrymx,ntdecr,ntmmin,npfit,npcmin,
     &              kuttaop,nrkmax0,nrkmax1,numstp,nwtfitp,nwtfitm,
     &              stepfac,flxstp,psispl,tolnwtp,tolnwtm,
     &              delpakf,delpakc,delpkf,delpkc,psichek,boxtnd,
     &              maptrace,norient,maxcutc,
     &              dresolv,dlclose,pntshft,endtol,
     &              narcln,nangax,nanglm,nbpmax,nwtmax,nslmax,nhfmax,
     &              bperor,sersnm,sertnm,arcmin,delgap,stepcut,
     &              novap,ioutm,iplotm
      common/size/  xdim,zdim,redge,zlowr,ipestg
      common/sarc/  ntmax,ntmsh,npcf,npcb,npc,xp(nlx),zp(nlx),gsq(nlx),
     &              arc(nlx),tp(nlx),arcc(nlx),tpp(nlx),bcd(4),
     &              csx(3,nlx),csz(3,nlx),cseq1(3,nlx),cseq2(3,nlx),
     &              st1(nlx),st2(nlx),st3(nlx),csveq(3,nlx),
     &              sv0(nlx),sv1(nlx),sv2(nlx),sv3(nlx),sv4(nlx),
     &              sv5(nlx)
      common/dblt/  psisep,xsep,zsep,smap(2),smax(2),smaz(2),numax,
     &              psimx(2),xax(2),zax(2),fitax(14,2),fitsep(14)
      common/ort1/  ndim1,ndim2,xma,zma,rcnt,btor,totcur,
     &              psimax0,psilim0,delpsir0,psimax1,psilim1,delpsir1,
     &              psimax,psilim,delpsir,dpsidec,dpsisv,dpedg1,dpedg2,
     &              eaxe,taxe,qaxe,fqiaxe,fqmaxe,paxe,faxe,ppaxe,ffpaxe,
     &              dnaxe,pfaxe,qppaxe,dlnaxe,dpnaxe,dhnaxe,
     &              arcaxe,dnnorm,pfnorm,
     &              qlim,fqilim,fqmlim,plim,flim,pplim,ffplim,
     &              dnlim,pflim,xlim,zlim,
     &              qmer,btmer,gkappa,arcprev,chicmax,
     &              charcmx,chpstmx,chhammx,arcnrm0,pstnrm0,hamnrm0,
     &              elongax1,elongax2,elongaxs,elongaxp,elongtn(npp)
      common/ort2/  psival(np1),pp(np1),ffp(np1),p(np1),f(np1),
     &              dnsty(np1),pfast(np1),qp(np1),fqpi(np1),fqim(np1),
     &              qdel(np1),psinitl(np1),qpinitl(np1),
     &              qprime(np1),qpint(np1),fqint(np1),
     &              psivl0(np2),rst(np2),zst(np2),rstp(np2),zstp(np2),
     &              rsinitl(np2),zsinitl(np2),
     &              arcsurf(nc1),pestchi(nc1),hamlchi(nc1),
     &              alfarc (nc1),alfpst (nc1),alfham (nc1),
     &              rs(nc1,np1),zs(nc1,np1)
      common/ort3/  rsrf(nc1),zsrf(nc1),chic(nc1),chie(nc1),dydx(nc1),
     &              xjsc(nc1)
      common/ort4/  rcell  (np1,nc1),zcell  (np1,nc1),
     &              dpsdr  (np1,nc1),dpsdz  (np1,nc1),
     &              chiarcl(np1,nc1),chipest(np1,nc1),chihaml(np1,nc1),
     &              xjcarcl(np1,nc1),xjcpest(np1,nc1),xjchaml(np1,nc1),
     &              alfarcl(np1,nc1),alfpest(np1,nc1),alfhaml(np1,nc1),
     &              arcnorm(np1),    pstnorm(np1),    hamnorm(np1),
     &              dlnlval(np1),    dpnlval(np1),    dhnlval(np1)
      common/ort5/  f3 (np1,nc1),f4 (np1,nc1),f5 (np1,nc1),
     &              f7 (np1,nc1),f8 (np1,nc1),f9 (np1,nc1),
     &              f10(np1,nc1),f11(np1,nc1),f12(np1,nc1),
     &              f13(np1,nc1),f14(np1,nc1),f15(np1,nc1),
     &              f16(np1,nc1),f17(np1,nc1),f18(np1,nc1),
     &              f19(np1,nc1),f20(np1,nc1),f21(np1,nc1),
     &              f22(np1,nc1),f23(np1,nc1),f24(np1,nc1),
     &              f25(np1,nc1),f26(np1,nc1),f27(np1,nc1),
     &              f28(np1,nc1)
      common/ort6/  savgax,tavgax,sntax0,sntax1,sntax2,sntax3,sntax4,
     &              svint(nvn),savge(np1),tavge(np1),sint0(np1),
     &              sint1(np1),sint2(np1),sint3(np1),sint4(np1)
      common/pldf/  x0min,x0max,y0min,y0max
      common/inout/ kutty,kuout,kuinp,kueql,kumap,kugrd,kurdm,kurdn,
     &              kufst,krdeqlb,krdmgta,krdngta,krdfast
c
c
c
c 1.0 Initialization
c
c 1.1 Return if no plot requested or mesh is too large
c
      if(iplotm .lt. 2) return
c
      nlxx     = nlx
      ithtp1   = itht + 1
c
      if(jpsi   .gt. nlxx) then
         call abortjob
     &        ('plotmap ',  1,   'Skipping map plot: jpsi   > nlx     '
     &        ,'jpsi    ', jpsi,     'nlxx    ', nlxx,     -1)
         return
      endif
c
      if(ithtp1 .gt. nlxx) then
         call abortjob
     &        ('plotmap ',  2,   'Skipping map plot: ithtp1 > nlx     '
     &        ,'ithtp1  ', ithtp1,   'nlxx    ', nlxx,     -1)
         return
      endif
c
c
c
c 2.0 Set up the plot page
c
c 2.1 Set up the map
c
      xa       = 0.000
      xb       = 1.00
      ya       = 0.000
      yb       = 1.000
      x1       = 0.100
      x2       = 0.90
      y1       = 0.100
      y2       = 0.900
      xhdp     = 0.360
      yhdp     = 0.975
      xhgd     = 0.3325
      yhgd     = 0.925
c
c
c
c 2.2 Print the heading for the box
c
      call map(xa,xb,ya,yb,xa,xb,ya,yb)
      call setlch(xhdp,yhdp,2,2,0,-1)
      call gtext('EQUILIBRIUM GRID',16,-1)
c
      call setlch(xhgd,yhgd,2,2,0,-1)
      if(igrid .eq. 0) call gtext('Equal Arc coordinates',21,-1)
      if(igrid .ne. 0) call gtext('Pest Coordinates',16,-1)
c
c
c 2.3 Set up the equilibrium box with an appropriate aspect
c     ratio to fit on the page
c
      ratio    = xdim/zdim
      ymidp    = zlowr + 0.5*zdim
      if    (ratio .gt. 1.0) then
        xl       = redge
        xr       = redge +      xdim
        yb       = ymidp - 0.50*xdim
        yt       = ymidp + 0.50*xdim
        half     = 0.5*(xdim-zdim)
        xmov     = 0.0
        xmrt     = 0.0
      elseif(ratio .le. 1.0) then
        xl       = redge
        xr       = redge +      zdim
        yb       = ymidp - 0.50*zdim
        yt       = ymidp + 0.50*zdim
        half     = 0.0
        xmov     = 0.0
        xmrt     = xdim - zdim
      endif
c
      xshift   = 0.5*xmrt
      yshift   = 0.0
c
      xlmap    = xl + xshift
      xrmap    = xr + xshift
      ybmap    = yb + yshift
      ytmap    = yt + yshift
c
c
c
c 2.4.2 Reset the mapping in the box
c
      call map(xlmap,xrmap,ybmap,ytmap,x1,x2,y1,y2)
c
c
c 2.5 Draw the box
c
      xllin    = xl + xmov
      xrlin    = xr + xmrt
      yblin    = yb + half
      ytlin    = yt - half
c
      call line(xllin, yblin, xllin, ytlin)
      call line(xllin, ytlin, xrlin, ytlin)
      call line(xrlin, ytlin, xrlin, yblin)
      call line(xrlin, yblin, xllin, yblin)
c
c
c
c 3.0 Plot the poloidal angle coordinates
c
      do 100 iaxe = 1,numax
c
c 3.1 Set up loop for each magnetic axis
c
      if    (jsep .gt. jpsi) then
        istart      = 1
        iend        = itht
      elseif(jsep .le. jpsi) then
        if    (iaxe .eq. 1) then
          istart      = 1
          iend        = itht/2
        elseif(iaxe .eq. 2) then
          istart      = itht/2 + 1
          iend        = itht
        endif
      endif
c
c
c 3.2 Loop over angles and plot the chi lines
c
      do  50 i    = istart,iend
      call setcrt(smax(iaxe),smaz(iaxe))
      xp(1)       = smax(iaxe)
      zp(1)       = smaz(iaxe)
      do  20 j    = 1,jpsi
      jj          = jpsi + 1 - j
      jc          = j+1
      xp(jc)      = rs(i,jj)
      zp(jc)      = zs(i,jj)
  20  continue
      xp(jpsi2)   = rsrf(i)
      zp(jpsi2)   = zsrf(i)
c
      call trace(xp,zp,jpsi2,1,1,0.0,0.0)
  50  continue
  100 continue
c
c
c
c
c 4.0 Plot the flux coordinate
c
      do 200 j    = 1,jpsi
c
c 4.1 Account for multiple axes
c
      if    (j .lt. jsep) then
        istart      = 1
        iend        = itht
      elseif(j .ge. jsep) then
        istart      = 1
        iend        = itht/2
      endif
c
      call setcrt(rs(istart,j),zs(istart,j))
c
      do 160 i    = istart,iend
      ik          = i-istart+1
      xp(ik)      = rs(i,j)
      zp(ik)      = zs(i,j)
  160 continue
c
      itotal      = iend - istart + 2
      xp(itotal)  = xp(1)
      zp(itotal)  = zp(1)
      call trace(xp,zp,itotal,1,1,0.0,0.0)
c
c
c
c 4.2 Plot the flux coordinate for the second axis if present
c
      if(j .ge. jsep) then
        istart      = itht/2 + 1
        iend        = itht
        call setcrt(rs(istart,j),zs(istart,j))
c
        do 180 i    = istart,iend
        ik          = i-istart+1
        xp(ik)      = rs(i,j)
        zp(ik)      = zs(i,j)
  180   continue
c
        itotal      = iend - istart + 2
        xp(itotal)  = xp(1)
        zp(itotal)  = zp(1)
        call trace(xp,zp,itotal,1,1,0.0,0.0)
      endif
  200 continue
c
c
c 4.3 Plot the outer surface
c
      istart      = 1
      iend        = ithtp
      call setcrt(rsrf(istart),zsrf(istart))
c
      do 300 i    = istart,iend
      ik          = i-istart+1
      xp(ik)      = rsrf(i)
      zp(ik)      = zsrf(i)
  300 continue
c
      itotal      = iend - istart + 1
      call trace(xp,zp,itotal,1,1,0.0,0.0)
c
c
c 5.0 Return and end
c
      call frame(0)
      return
      end
      subroutine ploterr(ktype)
c
c ------------------------
c plot equilibrium errors
c ------------------------
c
      parameter (npx=200,ncx=2*npx,np1=npx+1,np2=npx+2,nc1=ncx+1)
      parameter (nxx=513,nxz=nxx)
      parameter (npp=929,ntt=2*npp-1,npb=npp)
      parameter (nw2=2*nxx,nh2=2*nxz,nvk0=nw2+nh2)
      parameter (nlx=1441)
      parameter (nck=nlx)
      parameter (nsm=13)
      parameter (nvn=7)
c
      character*1   symbl,eqsymbl,tqsymbl
      character*1   symchar
c
      common/cnst/  pi,hlfpi,twopi,amu,boltzm,echarg,pmass,roundff,bigno
      common/inpd/  ntor,ncase,idnsty,ndnxp0,ndnxp1,ndnxp2,ifastp,
     &              nmap,neqtyp,mpreset,ndoublt,ndivert,nmtype,
     &              jpsi,itht,jpsi1,jpsi2,ithtp,jsep,
     &              isym,igrid,nham1,nham2,nham3,npowr,
     &              bfieldf,rdefolt,qxin,btdes,qsurf,
     &              gamma,gamav,rmantl,betaf,zeffect,
     &              ncorr,corrfac,nccellr,peqpk0,peqpk1,peqpk2
      common/mapd/  mapmaxd,dpsisl,dpsisd,nqaxis,nwtmag,nfitmax,nfitpts,
     &              ifitrng,jfitrng,jfitchk,fitchek,cnvmag,epsaxs,
     &              maxerlp,maxerlv,delbox,delboz,delac,delav,
     &              delstsf,delstlp,delstlv,prfrac,
     &              nerstop,nerprnt,qptol,tolspln,tolbchi,tolbtor,
     &              tolsymm,tolaugm,errsep,precisn,plossmx,
     &              narcmx,ntrymx,ntdecr,ntmmin,npfit,npcmin,
     &              kuttaop,nrkmax0,nrkmax1,numstp,nwtfitp,nwtfitm,
     &              stepfac,flxstp,psispl,tolnwtp,tolnwtm,
     &              delpakf,delpakc,delpkf,delpkc,psichek,boxtnd,
     &              maptrace,norient,maxcutc,
     &              dresolv,dlclose,pntshft,endtol,
     &              narcln,nangax,nanglm,nbpmax,nwtmax,nslmax,nhfmax,
     &              bperor,sersnm,sertnm,arcmin,delgap,stepcut,
     &              novap,ioutm,iplotm
      common/size/  xdim,zdim,redge,zlowr,ipestg
      common/eqd1/  nx,nz,nxd,nzd,dmx,dmz,
     &              x(nxx),z(nxz),psarray(nxx,nxz),
     &              gpx(nxx,nxz),gpy(nxx,nxz),grsq(nxx,nxz),
     &              cspl1(2,nxx,nh2),work0(nxz,2,nxx),vork0(nvk0,2)
      common/sarc/  ntmax,ntmsh,npcf,npcb,npc,xp(nlx),zp(nlx),gsq(nlx),
     &              arc(nlx),tp(nlx),arcc(nlx),tpp(nlx),bcd(4),
     &              csx(3,nlx),csz(3,nlx),cseq1(3,nlx),cseq2(3,nlx),
     &              st1(nlx),st2(nlx),st3(nlx),csveq(3,nlx),
     &              sv0(nlx),sv1(nlx),sv2(nlx),sv3(nlx),sv4(nlx),
     &              sv5(nlx)
      common/dblt/  psisep,xsep,zsep,smap(2),smax(2),smaz(2),numax,
     &              psimx(2),xax(2),zax(2),fitax(14,2),fitsep(14)
      common/dstr/  eqderor(nxx,nxz),toqeror(npp,ntt),surfpsi(nck),
     &              surferr(nck),surfcur(nck),nsrferr(nck),ksrferr(nck),
     &              toteror,toterov
      common/anal/  kntdlf(nsm),dlfac(nsm),symbl(nsm),
     &              eqsymbl(nxx,nxz),tqsymbl(npp,ntt)
      common/ort1/  ndim1,ndim2,xma,zma,rcnt,btor,totcur,
     &              psimax0,psilim0,delpsir0,psimax1,psilim1,delpsir1,
     &              psimax,psilim,delpsir,dpsidec,dpsisv,dpedg1,dpedg2,
     &              eaxe,taxe,qaxe,fqiaxe,fqmaxe,paxe,faxe,ppaxe,ffpaxe,
     &              dnaxe,pfaxe,qppaxe,dlnaxe,dpnaxe,dhnaxe,
     &              arcaxe,dnnorm,pfnorm,
     &              qlim,fqilim,fqmlim,plim,flim,pplim,ffplim,
     &              dnlim,pflim,xlim,zlim,
     &              qmer,btmer,gkappa,arcprev,chicmax,
     &              charcmx,chpstmx,chhammx,arcnrm0,pstnrm0,hamnrm0,
     &              elongax1,elongax2,elongaxs,elongaxp,elongtn(npp)
      common/ort2/  psival(np1),pp(np1),ffp(np1),p(np1),f(np1),
     &              dnsty(np1),pfast(np1),qp(np1),fqpi(np1),fqim(np1),
     &              qdel(np1),psinitl(np1),qpinitl(np1),
     &              qprime(np1),qpint(np1),fqint(np1),
     &              psivl0(np2),rst(np2),zst(np2),rstp(np2),zstp(np2),
     &              rsinitl(np2),zsinitl(np2),
     &              arcsurf(nc1),pestchi(nc1),hamlchi(nc1),
     &              alfarc (nc1),alfpst (nc1),alfham (nc1),
     &              rs(nc1,np1),zs(nc1,np1)
      common/ort3/  rsrf(nc1),zsrf(nc1),chic(nc1),chie(nc1),dydx(nc1),
     &              xjsc(nc1)
      common/ort4/  rcell  (np1,nc1),zcell  (np1,nc1),
     &              dpsdr  (np1,nc1),dpsdz  (np1,nc1),
     &              chiarcl(np1,nc1),chipest(np1,nc1),chihaml(np1,nc1),
     &              xjcarcl(np1,nc1),xjcpest(np1,nc1),xjchaml(np1,nc1),
     &              alfarcl(np1,nc1),alfpest(np1,nc1),alfhaml(np1,nc1),
     &              arcnorm(np1),    pstnorm(np1),    hamnorm(np1),
     &              dlnlval(np1),    dpnlval(np1),    dhnlval(np1)
      common/ort5/  f3 (np1,nc1),f4 (np1,nc1),f5 (np1,nc1),
     &              f7 (np1,nc1),f8 (np1,nc1),f9 (np1,nc1),
     &              f10(np1,nc1),f11(np1,nc1),f12(np1,nc1),
     &              f13(np1,nc1),f14(np1,nc1),f15(np1,nc1),
     &              f16(np1,nc1),f17(np1,nc1),f18(np1,nc1),
     &              f19(np1,nc1),f20(np1,nc1),f21(np1,nc1),
     &              f22(np1,nc1),f23(np1,nc1),f24(np1,nc1),
     &              f25(np1,nc1),f26(np1,nc1),f27(np1,nc1),
     &              f28(np1,nc1)
      common/ort6/  savgax,tavgax,sntax0,sntax1,sntax2,sntax3,sntax4,
     &              svint(nvn),savge(np1),tavge(np1),sint0(np1),
     &              sint1(np1),sint2(np1),sint3(np1),sint4(np1)
      common/toq0/  npsi,nthet,ntht1,neqsym,axddxz,
     &              peqmsh0(npp),peqdsd0(npp),peqdss0(npp),
     &              peqmsh1(npp),peqmsh2(npp),
     &              sqvl(npp),sfqi(npp),csqvl(3,npp),csfqi(3,npp)
      common/toq1/  thchi(ntt), seqrps(npp,ntt), seqzps(npp,ntt),
     &              csrvl(npp,ntt,4), cszvl(npp,ntt,4),
     &              pdsa(6),pdsr(6),pdsz(6),pdsp(6),pdsc(6)
      common/pldf/  x0min,x0max,y0min,y0max
      common/inout/ kutty,kuout,kuinp,kueql,kumap,kugrd,kurdm,kurdn,
     &              kufst,krdeqlb,krdmgta,krdngta,krdfast
c
c
c
c 1.0 Initialization
c
c 1.1 Return if no plot requested or mesh is too large
c
      if(iplotm .lt. 3) return
c
      nlxx      = nlx
      ithtp1    = itht + 1
c
      if(jpsi   .gt. nlxx) then
         call abortjob
     &        ('ploterr ',  1,   'Skipping error plot: jpsi > nlx     '
     &        ,'jpsi    ', jpsi,     'nlxx    ', nlxx,     -1)
         return
      endif
c
      if(ithtp1 .gt. nlxx) then
         call abortjob
     &        ('ploterr ',  2,   'Skipping error plot: ithtp > nlx    '
     &        ,'ithtp1  ', ithtp1,   'nlxx    ', nlxx,     -1)
         return
      endif
c
c
c 1.2 Set the plotting constants
c
      icasch    = 1
      isizch    = 1
      itypch    = 0
      kspac     = 1
c
      numch     = 1
      incpx     = 1
      incpy     = 1
      delpx     = 0.0
      delpy     = 0.0
c
c
c
c 2.0 Set up the plot page
c
c 2.1 Set up the map
c
      xa        = 0.000
      xb        = 1.00
      ya        = 0.000
      yb        = 1.000
      x1        = 0.100
      x2        = 0.90
      y1        = 0.100
      y2        = 0.900
      xhdp      = 0.275
      yhdp      = 0.975
      xhgd      = 0.370
      yhgd      = 0.925
c
c
c 2.2 Print the heading for the box
c
      call map(xa,xb,ya,yb,xa,xb,ya,yb)
      call setlch(xhdp,yhdp,2,2,0,-1)
      call gtext('Normalized Equilibrium Errors',29,-1)
c
      call setlch(xhgd,yhgd,2,2,0,-1)
      if(ktype  .eq. 0) call gtext('Direct equilibrium',18,-1)
      if(ktype  .gt. 0) call gtext('Inverse equilibrium',19,-1)
c
c
c 2.3 Set up the equilibrium mapping
c
c 2.3.1 Set up the equilibrium box
c
      ratio     = xdim/zdim
      ymidp     = zlowr + 0.5*zdim
      if    (ratio .gt. 1.0) then
        xl        = redge
        xr        = redge +      xdim
        yb        = ymidp - 0.50*xdim
        yt        = ymidp + 0.50*xdim
        half      = 0.5*(xdim-zdim)
        xmov      = 0.0
        xmrt      = 0.0
      elseif(ratio .le. 1.0) then
        xl        = redge
        xr        = redge +      zdim
        yb        = ymidp - 0.50*zdim
        yt        = ymidp + 0.50*zdim
        half      = 0.0
        xmov      = 0.0
        xmrt      = xdim - zdim
      endif
c
      xshift    = 0.5*xmrt
      yshift    = 0.0
c
      xlmap     = xl + xshift
      xrmap     = xr + xshift
      ybmap     = yb + yshift
      ytmap     = yt + yshift
c
c 2.3.2 Reset the map for the box
c
      call map(xlmap,xrmap,ybmap,ytmap,x1,x2,y1,y2)
c
c
c 2.4 Draw the box
c
      xllin     = xl + xmov
      xrlin     = xr + xmrt
      yblin     = yb + half
      ytlin     = yt - half
c
      call line(xllin, yblin, xllin, ytlin)
      call line(xllin, ytlin, xrlin, ytlin)
      call line(xrlin, ytlin, xrlin, yblin)
      call line(xrlin, yblin, xllin, yblin)
c
c
c
c 3.0 Plot the equilibrium errors
c
c 3.1 For direct equilibria
c
      if    (ktype  .eq. 0) then
c
c 3.1.1 Loop over the grid
c
        do 120 i    = 1,nx
        iv          = i
        xch         = x(iv)
        do 100 j    = 1,nz
        jv          = j
        ych         = z(jv)
        eqlerr      = abs(eqderor(iv,jv))
c
c 3.1.2 Set the symbol
c
        symchar     = eqsymbl(iv,jv)
c
c 3.1.3 Plot the symbol
c
        call setpch(icasch,isizch,itypch,kspac)
        call point (xch,ych)
        call pointc(symchar,xch,ych,numch,incpx,incpy,delpx,delpy)
c
  100   continue
  120   continue
c
c
c 3.2 For inverse equilibria
c
      elseif(ktype  .gt. 0) then
c
c 3.2.1 Loop over the grid
c
        do 220 jj   = 1,npsi
        jjv         = jj
        do 200 ii   = 1,nthet
        iiv         = ii
c
        xch         = seqrps(jjv,iiv)
        ych         = seqzps(jjv,iiv)
        eqlerr      = abs(toqeror(jjv,iiv))
c
c
c 3.2.2 Set the symbol
c
        symchar     = tqsymbl(jjv,iiv)
c
c 3.2.3 Plot the symbol
c
        call setpch(icasch,isizch,itypch,kspac)
        call point (xch,ych)
        call pointc(symchar,xch,ych,numch,incpx,incpy,delpx,delpy)
c
  200   continue
  220   continue
      endif
c
c
c
c 4.0 Plot the plasma boundary
c
c 4.1 Set the curve
c
      do 300 i  = 1,itht
      ip        = i
      xp(ip)    = rsrf(ip)
      zp(ip)    = zsrf(ip)
 300  continue
      xp(ithtp) = rsrf( 1)
      zp(ithtp) = zsrf( 1)
c
c
c 4.2 Plot the curve
c
      call trace(xp,zp,ithtp,incpx,incpy,delpx,delpy)
c
c
c
c 5.0 Return and end
c
      call frame(0)
      return
      end
      subroutine plottxt
c
c ---------------------
c plot equilibrium data
c ----------------------
c
      parameter (npx=200,ncx=2*npx,np1=npx+1,np2=npx+2,nc1=ncx+1)
      parameter (nxx=513,nxz=nxx)
      parameter (npp=929,ntt=2*npp-1,npb=npp)
      parameter (nw2=2*nxx,nh2=2*nxz,nvk0=nw2+nh2)
      parameter (nft=5)
      parameter (nvn=7)
c
      character*8   etitl,date
      character*80  string
c
      common/cnst/  pi,hlfpi,twopi,amu,boltzm,echarg,pmass,roundff,bigno
      common/inpd/  ntor,ncase,idnsty,ndnxp0,ndnxp1,ndnxp2,ifastp,
     &              nmap,neqtyp,mpreset,ndoublt,ndivert,nmtype,
     &              jpsi,itht,jpsi1,jpsi2,ithtp,jsep,
     &              isym,igrid,nham1,nham2,nham3,npowr,
     &              bfieldf,rdefolt,qxin,btdes,qsurf,
     &              gamma,gamav,rmantl,betaf,zeffect,
     &              ncorr,corrfac,nccellr,peqpk0,peqpk1,peqpk2
      common/mapd/  mapmaxd,dpsisl,dpsisd,nqaxis,nwtmag,nfitmax,nfitpts,
     &              ifitrng,jfitrng,jfitchk,fitchek,cnvmag,epsaxs,
     &              maxerlp,maxerlv,delbox,delboz,delac,delav,
     &              delstsf,delstlp,delstlv,prfrac,
     &              nerstop,nerprnt,qptol,tolspln,tolbchi,tolbtor,
     &              tolsymm,tolaugm,errsep,precisn,plossmx,
     &              narcmx,ntrymx,ntdecr,ntmmin,npfit,npcmin,
     &              kuttaop,nrkmax0,nrkmax1,numstp,nwtfitp,nwtfitm,
     &              stepfac,flxstp,psispl,tolnwtp,tolnwtm,
     &              delpakf,delpakc,delpkf,delpkc,psichek,boxtnd,
     &              maptrace,norient,maxcutc,
     &              dresolv,dlclose,pntshft,endtol,
     &              narcln,nangax,nanglm,nbpmax,nwtmax,nslmax,nhfmax,
     &              bperor,sersnm,sertnm,arcmin,delgap,stepcut,
     &              novap,ioutm,iplotm
      common/ttle/  etitl(nft),date
      common/size/  xdim,zdim,redge,zlowr,ipestg
      common/eqd1/  nx,nz,nxd,nzd,dmx,dmz,
     &              x(nxx),z(nxz),psarray(nxx,nxz),
     &              gpx(nxx,nxz),gpy(nxx,nxz),grsq(nxx,nxz),
     &              cspl1(2,nxx,nh2),work0(nxz,2,nxx),vork0(nvk0,2)
      common/dblt/  psisep,xsep,zsep,smap(2),smax(2),smaz(2),numax,
     &              psimx(2),xax(2),zax(2),fitax(14,2),fitsep(14)
      common/ort1/  ndim1,ndim2,xma,zma,rcnt,btor,totcur,
     &              psimax0,psilim0,delpsir0,psimax1,psilim1,delpsir1,
     &              psimax,psilim,delpsir,dpsidec,dpsisv,dpedg1,dpedg2,
     &              eaxe,taxe,qaxe,fqiaxe,fqmaxe,paxe,faxe,ppaxe,ffpaxe,
     &              dnaxe,pfaxe,qppaxe,dlnaxe,dpnaxe,dhnaxe,
     &              arcaxe,dnnorm,pfnorm,
     &              qlim,fqilim,fqmlim,plim,flim,pplim,ffplim,
     &              dnlim,pflim,xlim,zlim,
     &              qmer,btmer,gkappa,arcprev,chicmax,
     &              charcmx,chpstmx,chhammx,arcnrm0,pstnrm0,hamnrm0,
     &              elongax1,elongax2,elongaxs,elongaxp,elongtn(npp)
      common/ort2/  psival(np1),pp(np1),ffp(np1),p(np1),f(np1),
     &              dnsty(np1),pfast(np1),qp(np1),fqpi(np1),fqim(np1),
     &              qdel(np1),psinitl(np1),qpinitl(np1),
     &              qprime(np1),qpint(np1),fqint(np1),
     &              psivl0(np2),rst(np2),zst(np2),rstp(np2),zstp(np2),
     &              rsinitl(np2),zsinitl(np2),
     &              arcsurf(nc1),pestchi(nc1),hamlchi(nc1),
     &              alfarc (nc1),alfpst (nc1),alfham (nc1),
     &              rs(nc1,np1),zs(nc1,np1)
      common/ort3/  rsrf(nc1),zsrf(nc1),chic(nc1),chie(nc1),dydx(nc1),
     &              xjsc(nc1)
      common/ort4/  rcell  (np1,nc1),zcell  (np1,nc1),
     &              dpsdr  (np1,nc1),dpsdz  (np1,nc1),
     &              chiarcl(np1,nc1),chipest(np1,nc1),chihaml(np1,nc1),
     &              xjcarcl(np1,nc1),xjcpest(np1,nc1),xjchaml(np1,nc1),
     &              alfarcl(np1,nc1),alfpest(np1,nc1),alfhaml(np1,nc1),
     &              arcnorm(np1),    pstnorm(np1),    hamnorm(np1),
     &              dlnlval(np1),    dpnlval(np1),    dhnlval(np1)
      common/ort5/  f3 (np1,nc1),f4 (np1,nc1),f5 (np1,nc1),
     &              f7 (np1,nc1),f8 (np1,nc1),f9 (np1,nc1),
     &              f10(np1,nc1),f11(np1,nc1),f12(np1,nc1),
     &              f13(np1,nc1),f14(np1,nc1),f15(np1,nc1),
     &              f16(np1,nc1),f17(np1,nc1),f18(np1,nc1),
     &              f19(np1,nc1),f20(np1,nc1),f21(np1,nc1),
     &              f22(np1,nc1),f23(np1,nc1),f24(np1,nc1),
     &              f25(np1,nc1),f26(np1,nc1),f27(np1,nc1),
     &              f28(np1,nc1)
      common/ort6/  savgax,tavgax,sntax0,sntax1,sntax2,sntax3,sntax4,
     &              svint(nvn),savge(np1),tavge(np1),sint0(np1),
     &              sint1(np1),sint2(np1),sint3(np1),sint4(np1)
      common/vcal/  btnew,btave,bpave,betat,betap,betax0,betax1,volme,
     &              vhalf,bavet(3,nxx),bavep(3,nxx),pvolm(3,nxx),
     &              betav(3,nxx)
      common/volm/  pvansh,vp0,pmantl,vpm(np1),apm(np1),vcurnt(np1)
      common/geom/  aminor,rcgeom,epslon,capa,triangl,
     &              allim,bpornl,deltal,qcyl,q0surf,welln,
     &              well(np1),vdpsi(np1)
      common/pldf/  x0min,x0max,y0min,y0max
      common/inout/ kutty,kuout,kuinp,kueql,kumap,kugrd,kurdm,kurdn,
     &              kufst,krdeqlb,krdmgta,krdngta,krdfast
c
c
c
c 1.0 Plot the equilibrium eqdsk title
c
      if(iplotm .lt. 4) return
      smn     = 0.00
      smx     = 1.00
      smn0    = 0.00
      smx0    = 1.00
      call map(smn,smx,smn,smx,smn0,smx0,smn0,smx0)
c
      xxtt    = 0.20000
      yytt    = 0.98500
      call setlch(xxtt,yytt,2,2,0,-1)
      write(string,1000) (etitl(it),it=1,nft)
      call wrtstr(string,1)
c
c
c
c 2.0 Plot the text string
c
      xxtx0   = 0.050
      yytx0   = 0.915
      call setlch(xxtx0,yytx0,2,2,0,-1)
c
      call gtext('   jpsi', 7,-1)
      write(string,1010) jpsi
      call wrtstr(string,1)
      call gtext('   itht', 7,-1)
      write(string,1010) itht
      call wrtstr(string,1)
      call gtext('   nmap', 7,-1)
      write(string,1010) nmap
      call wrtstr(string,1)
      call gtext('   ncase', 8,-1)
      write(string,1010) ncase
      call wrtstr(string,1)
      call gtext('   isym', 7,-1)
      write(string,1010) isym
      call wrtstr(string,1)
      call gtext('   neqtyp', 9,-1)
      write(string,1010) neqtyp
      call wrtstr(string,1)
      call gtext('  ndoublt', 9,-1)
      write(string,1010) ndoublt
      call wrtstr(string,1)
      call gtext('    redge', 9,-1)
      write(string,1020) redge
      call wrtstr(string,1)
      call gtext('    xdim', 8,-1)
      write(string,1020) xdim
      call wrtstr(string,1)
      call gtext('    zdim', 8,-1)
      write(string,1020) zdim
      call wrtstr(string,1)
      call gtext('    psilim',10,-1)
      write(string,1020) psilim
      call wrtstr(string,1)
      call gtext('    psisep',10,-1)
      write(string,1020) psisep
      call wrtstr(string,1)
      call gtext('   psimaxu',10,-1)
      write(string,1020) psimx(1)
      call wrtstr(string,1)
      call gtext('   psimaxl',10,-1)
      write(string,1020) psimx(2)
      call wrtstr(string,1)
c
c
      xxtx1   = 0.400
      yytx1   = yytx0
      call setlch(xxtx1,yytx1,2,2,0,-1)
c
      call gtext('    btor', 8,-1)
      write(string,1020) btnew
      call wrtstr(string,1)
      call gtext('    totcur',10,-1)
      write(string,1020) totcur
      call wrtstr(string,1)
      call gtext('    rcnt', 8,-1)
      write(string,1020) rcnt
      call wrtstr(string,1)
      call gtext('    rcgeom',10,-1)
      write(string,1020) rcgeom
      call wrtstr(string,1)
      call gtext('    aminor',10,-1)
      write(string,1020) aminor
      call wrtstr(string,1)
      call gtext('    epslon',10,-1)
      write(string,1020) epslon
      call wrtstr(string,1)
      call gtext('    capa',8,-1)
      write(string,1020) capa
      call wrtstr(string,1)
      call gtext('  triangl',9,-1)
      write(string,1020) triangl
      call wrtstr(string,1)
      call gtext('    xaxis', 9,-1)
      write(string,1020) xax(1)
      call wrtstr(string,1)
      call gtext('    eaxe', 8,-1)
      write(string,1020) eaxe
      call wrtstr(string,1)
      call gtext('    taxe', 8,-1)
      write(string,1020) taxe
      call wrtstr(string,1)
      call gtext('    qaxe', 8,-1)
      write(string,1020) qaxe
      call wrtstr(string,1)
      call gtext('    qmer', 8,-1)
      write(string,1020) qmer
      call wrtstr(string,1)
      call gtext('    btmer', 9,-1)
      write(string,1020) btmer
      call wrtstr(string,1)
c
      xxtx1   = 0.750
      yytx1   = yytx0
      call setlch(xxtx1,yytx1,2,2,0,-1)
c
      call gtext('    qcyl', 8,-1)
      write(string,1020) qcyl
      call wrtstr(string,1)
      call gtext('    q0surf',10,-1)
      write(string,1020) q0surf
      call wrtstr(string,1)
      call gtext('    volume',10,-1)
      write(string,1020) volme
      call wrtstr(string,1)
      call gtext('    vhalf', 9,-1)
      write(string,1020) vhalf
      call wrtstr(string,1)
      call gtext('    welln', 9,-1)
      write(string,1020) welln
      call wrtstr(string,1)
      call gtext('    bpave', 9,-1)
      write(string,1020) bpave
      call wrtstr(string,1)
      call gtext('    btave', 9,-1)
      write(string,1020) btave
      call wrtstr(string,1)
      call gtext('    betap', 9,-1)
      write(string,1020) betap
      call wrtstr(string,1)
      call gtext('    betat', 9,-1)
      write(string,1020) betat
      call wrtstr(string,1)
      call gtext('    betax0',10,-1)
      write(string,1020) betax0
      call wrtstr(string,1)
      call gtext('    betax1',10,-1)
      write(string,1020) betax1
      call wrtstr(string,1)
      call gtext('    bpornl',10,-1)
      write(string,1020) bpornl
      call wrtstr(string,1)
      call gtext('    allim',9,-1)
      write(string,1020) allim
      call wrtstr(string,1)
      call gtext('    mantle',10,-1)
      write(string,1020) pmantl
      call wrtstr(string,1)
c
      call frame(0)
      return
c
 1000 format(5a8)
 1010 format(i6)
 1020 format(1pe12.5)
      end
      subroutine plotpak
c
c ------------------------------------------------
c plot rational surfaces and mesh packing weights
c ------------------------------------------------
c
c
      parameter (npx=200,ncx=2*npx,np1=npx+1,np2=npx+2,nc1=ncx+1)
      parameter (npp=929,ntt=2*npp-1,npb=npp)
      parameter (npk=np1,nrtdm=np1,nptt=2401)
c
      character*1   chrat, chwt,  chdw
      character*8   bswlab,bcslab,brtlab
      character*8   slabel,qlabel,wlabel,dlabel
      character*80  string
c
      common/cnst/  pi,hlfpi,twopi,amu,boltzm,echarg,pmass,roundff,bigno
      common/inpd/  ntor,ncase,idnsty,ndnxp0,ndnxp1,ndnxp2,ifastp,
     &              nmap,neqtyp,mpreset,ndoublt,ndivert,nmtype,
     &              jpsi,itht,jpsi1,jpsi2,ithtp,jsep,
     &              isym,igrid,nham1,nham2,nham3,npowr,
     &              bfieldf,rdefolt,qxin,btdes,qsurf,
     &              gamma,gamav,rmantl,betaf,zeffect,
     &              ncorr,corrfac,nccellr,peqpk0,peqpk1,peqpk2
      common/mapd/  mapmaxd,dpsisl,dpsisd,nqaxis,nwtmag,nfitmax,nfitpts,
     &              ifitrng,jfitrng,jfitchk,fitchek,cnvmag,epsaxs,
     &              maxerlp,maxerlv,delbox,delboz,delac,delav,
     &              delstsf,delstlp,delstlv,prfrac,
     &              nerstop,nerprnt,qptol,tolspln,tolbchi,tolbtor,
     &              tolsymm,tolaugm,errsep,precisn,plossmx,
     &              narcmx,ntrymx,ntdecr,ntmmin,npfit,npcmin,
     &              kuttaop,nrkmax0,nrkmax1,numstp,nwtfitp,nwtfitm,
     &              stepfac,flxstp,psispl,tolnwtp,tolnwtm,
     &              delpakf,delpakc,delpkf,delpkc,psichek,boxtnd,
     &              maptrace,norient,maxcutc,
     &              dresolv,dlclose,pntshft,endtol,
     &              narcln,nangax,nanglm,nbpmax,nwtmax,nslmax,nhfmax,
     &              bperor,sersnm,sertnm,arcmin,delgap,stepcut,
     &              novap,ioutm,iplotm
      common/dblt/  psisep,xsep,zsep,smap(2),smax(2),smaz(2),numax,
     &              psimx(2),xax(2),zax(2),fitax(14,2),fitsep(14)
      common/rpak/  nratnl,nwtot,cs(np1),
     &              qnval(np1),ldone(npk),mdone(npk),
     &              swgt (nptt), wght (nptt),
     &              dwds (nptt), qwgt (nptt),
     &              srat  (nrtdm),qrat   (nrtdm),dqds (nrtdm),
     &              swidth(nrtdm),sweight(nrtdm),jqrat(nrtdm),
     &              jplpak(nrtdm),jrepeat(nrtdm)
      common/ort1/  ndim1,ndim2,xma,zma,rcnt,btor,totcur,
     &              psimax0,psilim0,delpsir0,psimax1,psilim1,delpsir1,
     &              psimax,psilim,delpsir,dpsidec,dpsisv,dpedg1,dpedg2,
     &              eaxe,taxe,qaxe,fqiaxe,fqmaxe,paxe,faxe,ppaxe,ffpaxe,
     &              dnaxe,pfaxe,qppaxe,dlnaxe,dpnaxe,dhnaxe,
     &              arcaxe,dnnorm,pfnorm,
     &              qlim,fqilim,fqmlim,plim,flim,pplim,ffplim,
     &              dnlim,pflim,xlim,zlim,
     &              qmer,btmer,gkappa,arcprev,chicmax,
     &              charcmx,chpstmx,chhammx,arcnrm0,pstnrm0,hamnrm0,
     &              elongax1,elongax2,elongaxs,elongaxp,elongtn(npp)
      common/ort2/  psival(np1),pp(np1),ffp(np1),p(np1),f(np1),
     &              dnsty(np1),pfast(np1),qp(np1),fqpi(np1),fqim(np1),
     &              qdel(np1),psinitl(np1),qpinitl(np1),
     &              qprime(np1),qpint(np1),fqint(np1),
     &              psivl0(np2),rst(np2),zst(np2),rstp(np2),zstp(np2),
     &              rsinitl(np2),zsinitl(np2),
     &              arcsurf(nc1),pestchi(nc1),hamlchi(nc1),
     &              alfarc (nc1),alfpst (nc1),alfham (nc1),
     &              rs(nc1,np1),zs(nc1,np1)
      common/ort3/  rsrf(nc1),zsrf(nc1),chic(nc1),chie(nc1),dydx(nc1),
     &              xjsc(nc1)
      common/ort4/  rcell  (np1,nc1),zcell  (np1,nc1),
     &              dpsdr  (np1,nc1),dpsdz  (np1,nc1),
     &              chiarcl(np1,nc1),chipest(np1,nc1),chihaml(np1,nc1),
     &              xjcarcl(np1,nc1),xjcpest(np1,nc1),xjchaml(np1,nc1),
     &              alfarcl(np1,nc1),alfpest(np1,nc1),alfhaml(np1,nc1),
     &              arcnorm(np1),    pstnorm(np1),    hamnorm(np1),
     &              dlnlval(np1),    dpnlval(np1),    dhnlval(np1)
      common/pldf/  x0min,x0max,y0min,y0max
      common/inout/ kutty,kuout,kuinp,kueql,kumap,kugrd,kurdm,kurdn,
     &              kufst,krdeqlb,krdmgta,krdngta,krdfast
c
c
c
c 1.0 Initialization
c
c 1.1 Return if no plot is required
c
      if(iplotm .lt. 5) return
c
c
c 1.2 Initialize miscellaneous parameters
c
c 1.2.1 Roundoff for machine precision
c
      prndff   = roundff
c
c 1.2.2 Default parameters for pointc
c
      numbpc   = 1
      incpcx   = 0
      incpcy   = 0
      delpcx   = 0.0
      delpcy   = 0.0
c
c
c
c 2.0 Set up map and initialize parameters
c
c 2.1 Normalized plot dimensions
c
      smn      = 0.00
      smx      = 1.00
      smn0     = 0.00
      smx0     = 1.00
c
c
c 2.2 Physical plot dimensions
c
c 2.2.1 Plot borders
c
c 2.2.1.1 Borders for the profiles
      x0mn     = 0.120
      x0mx     = 0.820
      y0mn     = 0.100
      y0mx     = 0.800
c
c 2.2.1.2 Borders for the mesh arrays
      x1mn     = x0mn
      x1mx     = x0mx
      y1mn     = y0mx
      y1mx     = 0.925
c
c 2.2.2 Heading position
c
      xxhd     = 0.250
      yyhd     = 0.95
c
c 2.3 Label position parameters
c
c 2.3.1 Position increments
c
      nspx0    = 16
      nspy0    = 16
      nspz0    = 16
c
c 2.3.2 Axis label position fraction
c
      msps0    = (7*nspx0)/10
      mspq0    = (7*nspy0)/8
      mspw0    = (3*nspz0)/4
      mspd0    = (1*nspz0)/2
c
c
c 2.4 Set the normalized tick parameters
c
c 2.4.1 Tick scale skipping
c
      ntikxs0  = 2
      ntikys0  = 2
      ntikzs0  = 2
c
c 2.4.2 Normalized tick lengths
c
      xtiky0   = +0.300
      ytikx0   = +0.250
      ztikx0   = -0.250
c
c 2.4.3 Tick scale positions
c
c 2.4.3.1 Label shift along relevant axis
      xsclx0   = -0.250
      yscly0   = -0.250
      zsclz0   = -0.250
c
c 2.4.3.2 Label shift away from relevant axis
      xscly0   = -1.000
      ysclx0   = -2.100
      zsclx0   = +0.250
c
c
c 2.5 Set the normalized labels and positions
c
c 2.5.1 Horizontal axis
c
      slby10   = -1.500
      slabel   = 'sw '
c
c 2.5.2 Left  Vertical axis
c
      qlbx10   = -2.800
      qlabel   = 'n*q'
c
c 2.5.3 Right Vertical axis
c
      wlbx10   = +1.050
      dlbx10   = +1.050
      wlabel   = ' wght '
      dlabel   = 'dw/ds '
c
c
c 2.6 Set the character printing parameters
c
c 2.6.1 Heading: case, size, and orientation
c
      icashd   = +2
      isizhd   = +3
      iornhd   =  0
c
c 2.6.2 Horizontal axis
c
c 2.6.2.1 Horizontal axis scale: case, size, and orientation
      icshsc   = +1
      iszhsc   = +1
      iorhsc   =  0
c
c 2.6.2.2 Horizontal axis label: case, size, and orientation
      icshlb   = +1
      iszhlb   = +3
      iorhlb   =  0
c
c 2.6.3 Left Vertical axis
c
c 2.6.3.1 Left Vertical axis scale: case, size, and orientation
      icsvsy   = +1
      iszvsy   = +1
      iorvsy   =  0
c
c 2.6.3.2 Left Vertical axis label: case, size, and orientation
      icslbq   = +1
      iszlbq   = +3
      iorlbq   =  0
c
c 2.6.4 Right Vertical axis
c
c 2.6.4.1 Right Vertical axis scale: case, size, and orientation
      icsvsz   = +1
      iszvsz   = +1
      iorvsz   =  0
c
c 2.6.4.2 Right Vertical axis label: case, size, and orientation
      icslbw   = +1
      iszlbw   = +3
      iorlbw   =  0
c
      icslbd   = +1
      iszlbd   = +3
      iorlbd   =  0
c
c
c 2.7 Line type parameters
c
c 2.7.1 Identifying characters
c
c 2.7.1.1 Rational surfaces
      icsrt    = 2
      iszrt    = 0
      itprt    = 0
      ksprt    = 1
c
      icsrat   = +1
      iszrat   = +1
      iorrat   =  0
      chrat    = '*'
c
c 2.7.1.2 Weight function
      icswt    = 2
      iszwt    = 0
      itpwt    = 0
      kspwt    = 1
c
      icschw   = +1
      iszchw   = +1
      iorchw   =  0
      chwt     = 'x'
c
c 2.7.1.3 Derivative weight function
      icsdw    = 2
      iszdw    = 0
      itpdw    = 0
      kspdw    = 1
c
      icschd   = +1
      iszchd   = +1
      iorchd   =  0
      chdw     = '+'
c
c
c 2.8 Set the parameters for the extra rational surface lines
c
c 2.8.1 Set the line type
c
      kratsp   = 2
c
c
c 2.9 Set the parameters for the box above the main plot
c
c 2.9.1 Set the number of mesh types
c
      numbox   = 3
c
c 2.9.2 Set the map parameters for the box
c
      y1min    = 0.00
      y1max    = 1.00
c
c 2.9.3 Set the normalized tick lengths for the meshes
c
      ytikln0  = 0.3
c
c 2.9.4 Set the labels and label positions
c
      bswlbx0  = -1.500
      bswlby0  = +0.500
      bswlab   = ' sw '
c
      bcslbx0  = -1.500
      bcslby0  = +0.500
      bcslab   = ' cs '
c
      brtlbx0  = -1.500
      brtlby0  = +0.500
      brtlab   = 'srat'
c
c 2.9.5 Set the character printing case, size, and orientation
c
c 2.9.5.1 Text characters
c
      ibswcs   = +1
      ibswsz   = +2
      ibswor   =  0
c
      ibcscs   = +1
      ibcssz   = +2
      ibcsor   =  0
c
      ibrtcs   = +1
      ibrtsz   = +2
      ibrtor   =  0
c
c 2.9.5.2 Curve character parameters
c
c
c
c 3.0 Find the function minima and maxima
c
c 3.1 Search the physically scaled data
c
      swmin    = swgt  (  1  )
      swmax    = swgt  (nwtot)
      qwmin    =  0.00
      qwmax    = -1.0e+10
      do 10  k = 1,nwtot
      qwmin    = amin1(qwmin,qwgt (k))
      qwmax    = amax1(qwmax,qwgt (k))
   10 continue
c
c
c 3.2 Adjust the maxima as needed
c
c 3.2.1 Physical n*q(sw)  profile maximum
c
      iqmx     = ifix (qwmax) + 1
      qwmax    = float(iqmx)
c
      prmin    = 0.00
      prmax    = 1.00
c
c
c
c 4.0 Set the physical plot parameters
c
c 4.1 Set the number of tick marks
c
c 4.1.1 Horizontal axis
      nspx     = 10
c
c 4.1.2 Left  Vertical axis
c
      nspy1    = ifix(qwmax)
      nspy1    = min0(nspy1,nspy0)
      nspy1    = max0(nspy1,1)
c
      if    (nspy1 .le. 4) then
        if    (nspy1 .le. 1) then
          nspy   = 5*nspy1
        elseif(nspy1 .le. 2) then
          nspy   = 4*nspy1
        else
          nspy   = 2*nspy1
        endif
c
      elseif(nspy1 .gt. 4) then
        nspy   =   nspy1
      endif
c
c 4.1.3 Right Vertical axis
c
      zspz     = 10.0*prmax + prndff
      nspz     = ifix(zspz)
c
c
c 4.2 Normalize the physical spacing increments to the physical range
c
c 4.2.1 Standard spacing
c
      delx0    = (swmax-swmin)/nspx0
      dely0    = (qwmax-qwmin)/nspy0
      delz0    = (prmax-prmin)/nspz0
c
c 4.2.2 Tick spacing
c
      delx     = (swmax-swmin)/nspx
      dely     = (qwmax-qwmin)/nspy
      delz     = (prmax-prmin)/nspz
c
c
c 4.3 Set the physical spacings
c
c 4.3.1 Tick scale values
c
      ntikxs   = min0(ntikxs0,nspx)
      ntikxs   = max0(ntikxs,  1  )
      if(nspx .lt. 5) ntikxs  = 1
c
      ntikys   = min0(ntikys0,nspy)
      ntikys   = max0(ntikys,  1  )
      if(nspy .lt. 5) ntikys  = 1
c
      ntikzs   = min0(ntikzs0,nspz)
      ntikzs   = max0(ntikzs,  1  )
      if(nspz .lt. 5) ntikzs  = 1
c
c 4.3.2 Axis scale positions
c
c
c 4.3.2.1 Label shift along relevant axis
      xsclx    = xsclx0*delx0
      yscly    = yscly0*dely0
      zsclz    = zsclz0*delz0
c
c 4.3.2.2 Label shift away from relevant axis
      xscly    = xscly0*dely0
      ysclx    = ysclx0*delx0
      zsclx    = zsclx0*delx0
c
c 4.3.2 Tick lengths
c
      xtiky    = xtiky0*dely0
      ytikx    = ytikx0*delx0
      ztikx    = ztikx0*delx0
c
c 4.3.4 Label positions
c
      slby1    = slby10*dely0
      qlbx1    = qlbx10*delx0
      wlbx1    = wlbx10*delx0
      dlbx1    = dlbx10*delx0
c
c 4.3.5 Label positions for the upper box
c
      bswlbx   = bswlbx0*delx0
      bswlby   = bswlby0
c
      bcslbx   = bcslbx0*delx0
      bcslby   = bcslby0
c
      brtlbx   = brtlbx0*delx0
      brtlby   = brtlby0
c
c
c
c 5.0 Print out the plot heading
c
c 5.1 Set up the plot in normalized coordinates
c
      call map(smn,smx,smn,smx,smn0,smx0,smn0,smx0)
c
c
c 5.2 Set the position and plot the text
c
      call setlch(xxhd,  yyhd,  icashd,isizhd,iornhd,-1)
      call gtext('Weight function vs sw ',22,-1)
c
c
c
c 6.0 Plot n*q = qwgt(swgt)
c
c 6.1 Renormalize the map for the n*q(sw) profile and plot the borders
c
c 6.1.1 Set the plot page in physical units
c
      call map(swmin,swmax,qwmin,qwmax,x0mn,x0mx,y0mn,y0mx)
c
c 6.1.2 Plot the four borders in physical units
c
      call line(swmin,qwmin,swmin,qwmax)
      call line(swmin,qwmax,swmax,qwmax)
      call line(swmax,qwmax,swmax,qwmin)
      call line(swmin,qwmin,swmax,qwmin)
c
c
c 6.2 Set up the Horizontal axis
c
c 6.2.1 Plot the horizontal axis tick marks and scale
c
      qpos     = qwmin
      qpos1    = qpos
      qpos2    = qpos + xtiky
      qpossc   = qpos + xscly
      do 100 i = 0,nspx
      iskpp    = (i/ntikxs)*ntikxs
      spos     = swmin + i*delx
      spossc   = spos  + xsclx
c
c 6.2.1.1 Plot the tick mark
      call line  (spos,  qpos1, spos,  qpos2 )
c
c 6.2.1.2 Plot the scale value
      if(i .eq. iskpp) then
        call setlch(spossc,qpossc,icshsc,iszhsc,iorhsc,-1)
        write(string,1000) spos
        call wrtstr(string,1)
      endif
  100 continue
c
c 6.2.2 Plot the horizontal axis label
c
      sposlb   = swmin  + msps0*delx0
      qposlb   = qwmin  + slby1
      call setlch(sposlb,qposlb,icshlb,iszhlb,iorhlb,-1)
      write(string,1010) slabel
      call wrtstr(string,1)
c
c
c 6.3 Set up the Left  vertical axis
c
c 6.3.1 Plot the  vertical  axis tick marks and scale
c
      spos     = swmin
      spos1    = spos
      spos2    = spos + ytikx
      spossc   = spos + ysclx
      do 110 j = 0,nspy
      jskpp    = (j/ntikys)*ntikys
      qpos     = qwmin + j*dely
      qpossc   = qpos  + yscly
c
c 6.3.1.1 Plot the tick mark
      call line  (spos1, qpos,  spos2, qpos  )
c
c 6.3.1.2 Plot the scale value
      if(j .eq. jskpp) then
        call setlch(spossc,qpossc,icsvsy,iszvsy,iorvsy,-1)
        write(string,1100) qpos
        call wrtstr(string,1)
      endif
  110 continue
c
c 6.3.2 Plot the  vertical  axis label
c
      sposlb   = swmin  + qlbx1
      qposlb   = qwmin  + mspq0*dely0
      call setlch(sposlb,qposlb,icslbq,iszlbq,iorlbq,-1)
      write(string,1110) qlabel
      call wrtstr(string,1)
c
c
c 6.4 Plot the n*q(sw)  profile
c
      call trace(swgt,qwgt,nwtot,-1,-1,0.0,0.0)
c
c
c 6.5 Plot the rational surfaces on the plot of n*q(sw)
c
c 6.5.1 Reset the printing parameters for pointc
c
      call setpch(icsrt,iszrt,itprt,ksprt)
      call setlch(x0mn,y0mn,icsrat,iszrat,iorrat,-1)
c
c 6.5.2 Plot the line out to n*q(srat) and down to srat
c
      do 120 lr  = 1,nratnl
      swrat      = srat (lr)
      qwrat      = qrat (lr)
      call linep (swmin,qwrat,swrat,qwrat,kratsp)
      call linep (swrat,qwrat,swrat,qwmin,kratsp)
      call pointc(chrat,swrat,qwrat,numbpc,incpcx,incpcy,delpcx,delpcy)
  120 continue
c
c
c
c 7.0 Plot the normalized profiles
c
c 7.1 Renormalize the map for the normalized profiles
c
      call map(swmin,swmax,prmin,prmax,x0mn,x0mx,y0mn,y0mx)
c
c
c 7.2 Set up the Right Vertical axis
c
c 7.2.1 Plot the  vertical  axis tick marks and scale
c
      spos     = swmax
      spos1    = spos + ztikx
      spos2    = spos
      spossc   = spos + zsclx
      do 150 j = 0,nspz
      jskpp    = (j/ntikzs)*ntikzs
      ppos     = prmin + j*delz
      ppossc   = ppos  + zsclz
c
c 7.2.1.1 Plot the tick mark
      call line  (spos1, ppos,  spos2, ppos  )
c
c 7.2.1.2 Plot the scale value
      if(j .eq. jskpp) then
        call setlch(spossc,ppossc,icsvsz,iszvsz,iorvsz,-1)
        write(string,1200) ppos
        call wrtstr(string,1)
      endif
  150 continue
c
c 7.2.2 Plot the  vertical  axis labels
c
c 7.2.2.1 Weight function label
      sposlw   = swmax  + wlbx1
      pposlw   = prmin  + mspw0*delz0
      call setlch(sposlw,pposlw,icslbw,iszlbw,iorlbw,-1)
      write(string,1210) wlabel
      call wrtstr(string,1)
c
c 7.2.2.1 Weight function derivative label
      sposld   = swmax  + dlbx1
      pposld   = prmin  + mspd0*delz0
      call setlch(sposld,pposld,icslbd,iszlbd,iorlbd,-1)
      write(string,1220) dlabel
      call wrtstr(string,1)
c
c
c
c 7.3 Plot the normalized profiles
c
c 7.3.1 Plot the weight function profile
c
c 7.3.1.1 Find the maximum value for normalizing the profile
      wtmax    = 0.0
      do 200 k = 1,nwtot
      if(abs(wght  (k)) .gt. abs(wtmax)) wtmax  = wght  (k)
  200 continue
c
      if(wtmax .eq. 0.0) wtmax  = 1.0
c
c 7.3.1.2 Reset the printing parameters for pointc
      call setpch(icswt,iszwt,itpwt,kspwt)
      call setlch(x0mn,y0mn,icschw,iszchw,iorchw,-1)
c
c 7.3.1.3 Plot the normalized values
      do 210 k = 1,nwtot
      swpt     =     swgt  (k)
      wtpt     = abs(wght  (k)/wtmax)
      call pointc(chwt,swpt,wtpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  210 continue
c
c 7.3.2 Plot the weight function derivative profile
c
c 7.3.2.1 Find the maximum value for normalizing the profile
      dwmax    = 0.0
      do 300 k = 1,nwtot
      if(abs(dwds  (k)) .gt. abs(dwmax)) dwmax  = dwds  (k)
  300 continue
c
      if(dwmax .eq. 0.0) dwmax  = 1.0
c
c 7.3.2.2 Reset the printing parameters for pointc
      call setpch(icsdw,iszdw,itpdw,kspdw)
      call setlch(x0mn,y0mn,icschd,iszchd,iorchd,-1)
c
c 7.3.2.3 Plot the normalized values
      do 310 k = 1,nwtot
      swpt     =     swgt  (k)
      dwpt     = abs(dwds  (k)/dwmax)
      call pointc(chdw,swpt,dwpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  310 continue
c
c
c
c 8.0 Plot the mesh values in the second box
c
c 8.1 Reset the map for the new box above the main plot
c
      call map(swmin,swmax,y1min,y1max,x1mn,x1mx,y1mn,y1mx)
c
c
c 8.2 Draw in the box
c
      call line(swmin,y1min,swmin,y1max)
      call line(swmin,y1max,swmax,y1max)
      call line(swmax,y1max,swmax,y1min)
      call line(swmax,y1min,swmin,y1min)
c
c
c 8.3 Set up and draw the dividing lines
c
c 8.3.1 Set up the dividing lines
c
      ytotal   = y1max - y1min
      ydivide  = ytotal/float(numbox)
      y1third  = y1min  +  1.0*ydivide
      y2third  = y1min  +  2.0*ydivide
c
      ytikln   = ytikln0*ytotal
c
      ysw0     = y1min
      ysw1     = ysw0   +  ytikln
      ycs0     = y1third
      ycs1     = ycs0   +  ytikln
      yrt0     = y2third
      yrt1     = yrt0   +  ytikln
c
c 8.3.2 Draw in the dividing lines
c
      call line  (swmin,y1third,swmax,y1third)
      call line  (swmin,y2third,swmax,y2third)
c
c
c 8.4 Plot the swgt mesh
c
c 8.4.1 Place tick marks
c
      do 400 k = 1,nwtot
      swpt     = swgt(k)
      call line  (swpt,ysw0,swpt,ysw1)
  400 continue
c
c 8.4.2 Print the label
c
      xswpos   = swmin + bswlbx
      yswpos   = ysw0  + bswlby*ydivide
c
      call setlch(xswpos,yswpos,ibswcs,ibswsz,ibswor,-1)
      write(string,2000) bswlab
      call wrtstr(string,1)
c
c
c 8.5 Plot the new cs mesh
c
c 8.5.1 Place tick marks
c
      do 420 j = 1,jpsi1
      cspt     = cs  (j)
      call line  (cspt,ycs0,cspt,ycs1)
  420 continue
c
c 8.5.2 Print the label
c
      xcspos   = swmin + bcslbx
      ycspos   = ycs0  + bcslby*ydivide
c
      call setlch(xcspos,ycspos,ibcscs,ibcssz,ibcsor,-1)
      write(string,2010) bcslab
      call wrtstr(string,1)
c
c
c
c 8.6 Plot the packing values
c
c 8.6.1 Place tick marks
c
      do 440 lr  = 1,nratnl
      rtpt       = srat(lr)
      call line  (rtpt,yrt0,rtpt,yrt1)
  440 continue
c
c 8.6.2 Print the label
c
      xrtpos     = swmin + brtlbx
      yrtpos     = yrt0  + brtlby*ydivide
c
      call setlch(xrtpos,yrtpos,ibrtcs,ibrtsz,ibrtor,-1)
      write(string,2020) brtlab
      call wrtstr(string,1)
c
c
c
c 9.0 Close frame, return, and end
c
      call frame(0)
      return
c
 1000 format(f3.1)
 1010 format(a3)
 1100 format(f7.2)
 1110 format(a3)
 1200 format(f4.2)
 1210 format(a6)
 1220 format(a6)
 2000 format(a6)
 2010 format(a6)
 2020 format(a6)
      end
      subroutine plotqpr
c
c ----------------------------------------------------
c plot q and normalized pressure and density profiles
c ----------------------------------------------------
c
      parameter (npx=200,ncx=2*npx,np1=npx+1,np2=npx+2,nc1=ncx+1)
      parameter (npp=929,ntt=2*npp-1,npb=npp)
      parameter (nlx=1441)
      parameter (npq=np2)
      parameter (nvn=7)
c
      character*1   chpp,  chdn,  chpf
      character*1   lbchpv, lbchsp, lbchrh, lbchro,
     &              lbchtf, lbchst ,lbchvl, lbchsv
      character*8   pclab,  tclab,  labpsi, 
     &              labpsiv,labpsin,labsqtp,
     &              labrhov,labrho2,labrhon,
     &              labtorf,labtorn,labsqtt,
     &              labvolm,labvoln,labsqtv,
     &              labdrh2,labdrho,labdtor,
     &              labdsqt,labdvol,labdsqv
      character*8   slabel,qlabel,plabel,dlabel,flabel
      character*80  string
c
      common/cnst/  pi,hlfpi,twopi,amu,boltzm,echarg,pmass,roundff,bigno
      common/inpd/  ntor,ncase,idnsty,ndnxp0,ndnxp1,ndnxp2,ifastp,
     &              nmap,neqtyp,mpreset,ndoublt,ndivert,nmtype,
     &              jpsi,itht,jpsi1,jpsi2,ithtp,jsep,
     &              isym,igrid,nham1,nham2,nham3,npowr,
     &              bfieldf,rdefolt,qxin,btdes,qsurf,
     &              gamma,gamav,rmantl,betaf,zeffect,
     &              ncorr,corrfac,nccellr,peqpk0,peqpk1,peqpk2
      common/mapd/  mapmaxd,dpsisl,dpsisd,nqaxis,nwtmag,nfitmax,nfitpts,
     &              ifitrng,jfitrng,jfitchk,fitchek,cnvmag,epsaxs,
     &              maxerlp,maxerlv,delbox,delboz,delac,delav,
     &              delstsf,delstlp,delstlv,prfrac,
     &              nerstop,nerprnt,qptol,tolspln,tolbchi,tolbtor,
     &              tolsymm,tolaugm,errsep,precisn,plossmx,
     &              narcmx,ntrymx,ntdecr,ntmmin,npfit,npcmin,
     &              kuttaop,nrkmax0,nrkmax1,numstp,nwtfitp,nwtfitm,
     &              stepfac,flxstp,psispl,tolnwtp,tolnwtm,
     &              delpakf,delpakc,delpkf,delpkc,psichek,boxtnd,
     &              maptrace,norient,maxcutc,
     &              dresolv,dlclose,pntshft,endtol,
     &              narcln,nangax,nanglm,nbpmax,nwtmax,nslmax,nhfmax,
     &              bperor,sersnm,sertnm,arcmin,delgap,stepcut,
     &              novap,ioutm,iplotm
      common/dblt/  psisep,xsep,zsep,smap(2),smax(2),smaz(2),numax,
     &              psimx(2),xax(2),zax(2),fitax(14,2),fitsep(14)
      common/sarc/  ntmax,ntmsh,npcf,npcb,npc,xp(nlx),zp(nlx),gsq(nlx),
     &              arc(nlx),tp(nlx),arcc(nlx),tpp(nlx),bcd(4),
     &              csx(3,nlx),csz(3,nlx),cseq1(3,nlx),cseq2(3,nlx),
     &              st1(nlx),st2(nlx),st3(nlx),csveq(3,nlx),
     &              sv0(nlx),sv1(nlx),sv2(nlx),sv3(nlx),sv4(nlx),
     &              sv5(nlx)
      common/ort1/  ndim1,ndim2,xma,zma,rcnt,btor,totcur,
     &              psimax0,psilim0,delpsir0,psimax1,psilim1,delpsir1,
     &              psimax,psilim,delpsir,dpsidec,dpsisv,dpedg1,dpedg2,
     &              eaxe,taxe,qaxe,fqiaxe,fqmaxe,paxe,faxe,ppaxe,ffpaxe,
     &              dnaxe,pfaxe,qppaxe,dlnaxe,dpnaxe,dhnaxe,
     &              arcaxe,dnnorm,pfnorm,
     &              qlim,fqilim,fqmlim,plim,flim,pplim,ffplim,
     &              dnlim,pflim,xlim,zlim,
     &              qmer,btmer,gkappa,arcprev,chicmax,
     &              charcmx,chpstmx,chhammx,arcnrm0,pstnrm0,hamnrm0,
     &              elongax1,elongax2,elongaxs,elongaxp,elongtn(npp)
      common/ort2/  psival(np1),pp(np1),ffp(np1),p(np1),f(np1),
     &              dnsty(np1),pfast(np1),qp(np1),fqpi(np1),fqim(np1),
     &              qdel(np1),psinitl(np1),qpinitl(np1),
     &              qprime(np1),qpint(np1),fqint(np1),
     &              psivl0(np2),rst(np2),zst(np2),rstp(np2),zstp(np2),
     &              rsinitl(np2),zsinitl(np2),
     &              arcsurf(nc1),pestchi(nc1),hamlchi(nc1),
     &              alfarc (nc1),alfpst (nc1),alfham (nc1),
     &              rs(nc1,np1),zs(nc1,np1)
      common/ort3/  rsrf(nc1),zsrf(nc1),chic(nc1),chie(nc1),dydx(nc1),
     &              xjsc(nc1)
      common/ort4/  rcell  (np1,nc1),zcell  (np1,nc1),
     &              dpsdr  (np1,nc1),dpsdz  (np1,nc1),
     &              chiarcl(np1,nc1),chipest(np1,nc1),chihaml(np1,nc1),
     &              xjcarcl(np1,nc1),xjcpest(np1,nc1),xjchaml(np1,nc1),
     &              alfarcl(np1,nc1),alfpest(np1,nc1),alfhaml(np1,nc1),
     &              arcnorm(np1),    pstnorm(np1),    hamnorm(np1),
     &              dlnlval(np1),    dpnlval(np1),    dhnlval(np1)
      common/ort5/  f3 (np1,nc1),f4 (np1,nc1),f5 (np1,nc1),
     &              f7 (np1,nc1),f8 (np1,nc1),f9 (np1,nc1),
     &              f10(np1,nc1),f11(np1,nc1),f12(np1,nc1),
     &              f13(np1,nc1),f14(np1,nc1),f15(np1,nc1),
     &              f16(np1,nc1),f17(np1,nc1),f18(np1,nc1),
     &              f19(np1,nc1),f20(np1,nc1),f21(np1,nc1),
     &              f22(np1,nc1),f23(np1,nc1),f24(np1,nc1),
     &              f25(np1,nc1),f26(np1,nc1),f27(np1,nc1),
     &              f28(np1,nc1)
      common/ort6/  savgax,tavgax,sntax0,sntax1,sntax2,sntax3,sntax4,
     &              svint(nvn),savge(np1),tavge(np1),sint0(np1),
     &              sint1(np1),sint2(np1),sint3(np1),sint4(np1)
      common/auxmsh/rh2tot,rhotot,flxtot,voltot,
     &              anltflxp,dvrtflxp,anltvolp,dvrtvolp,
     &              analtflx,divrtflx,analtvol,divrtvol,
     &              psivl1 (np2),qp1    (np2),fqpi1  (np2),
     &              qprim1 (np2),sintp0 (np2),psivmsh(np2),
     &              psivalu(np2),psinorm(np2),psisqrt(np2),
     &              psirho2(np2),psirhov(np2),psinrho(np2),
     &              psitorf(np2),psintor(np2),psisqtf(np2),
     &              psivolm(np2),psinvol(np2),psisqvl(np2),
     &              dpsirh2(np2),dpsirho(np2),dpsinrh(np2),
     &              dpsitor(np2),dpsintf(np2),dpsisqt(np2),
     &              dpsivol(np2),dpsinvl(np2),dpsisqv(np2)
       common/ratnl/jqmin, jqmax, numbqs,nq1totl,
     &              rcentr,zcentr,rminor,zminor,
     &              asprat00,asprat01,asprat10,asprat11,
     &              pminvl,qminvl,qpminv,pmaxvl,qmaxvl,qpmaxv,
     &              psivlmnq,    psivnmnq,    psisqmnq,
     &              psir2mnq,    psirhmnq,    psinrmnq,
     &              psitfmnq,    psintmnq,    psistmnq,
     &              psivmmnq,    psinvmnq,    psisvmnq,
     &              dpsr2mnq,    dpsrhmnq,    dpsnrmnq,
     &              dpstfmnq,    dpsntmnq,    dpsstmnq,
     &              dpsvmmnq,    dpsnvmnq,    dpssvmnq,
     &              psivlmxq,    psivnmxq,    psisqmxq,
     &              psir2mxq,    psirhmxq,    psinrmxq,
     &              psitfmxq,    psintmxq,    psistmxq,
     &              psivmmxq,    psinvmxq,    psisvmxq,
     &              dpsr2mxq,    dpsrhmxq,    dpsnrmxq,
     &              dpstfmxq,    dpsntmxq,    dpsstmxq,
     &              dpsvmmxq,    dpsnvmxq,    dpssvmxq,
     &              lpsiq  (npq),psivlq (npq),
     &              qprimq (npq),qvalue (npq),psimshq(npq),
     &              psivalq(npq),psinrmq(npq),psisqrq(npq),
     &              psirh2q(npq),psirhoq(npq),psinrhq(npq),
     &              psitorq(npq),psintfq(npq),psisqtq(npq),
     &              psivolq(npq),psinvlq(npq),psisqvq(npq),
     &              dpsir2q(npq),dpsirhq(npq),dpsinrq(npq),
     &              dpsitfq(npq),dpsintq(npq),dpsistq(npq),
     &              dpsivlq(npq),dpsinvq(npq),dpsisvq(npq),
     &              shearps(npq),shearrh(npq),
     &              sheartf(npq),shearvl(npq),
     &              epslrh1 (npq),shearrh1(npq),shearfrh(npq),
     &              epslvl1 (npq),shearvl1(npq),shearfvl(npq)
      common/labels/lbchpv, lbchsp, lbchrh, lbchro,
     &              lbchtf, lbchst, lbchvl, lbchsv,
     &              pclab,  tclab,  labpsi,
     &              labpsiv,labpsin,labsqtp,
     &              labrhov,labrho2,labrhon,
     &              labtorf,labtorn,labsqtt,
     &              labvolm,labvoln,labsqtv,
     &              labdrh2,labdrho,labdtor,
     &              labdsqt,labdvol,labdsqv
      common/pldf/  x0min,x0max,y0min,y0max
      common/inout/ kutty,kuout,kuinp,kueql,kumap,kugrd,kurdm,kurdn,
     &              kufst,krdeqlb,krdmgta,krdngta,krdfast
c
      dimension profle(nlx)
      equivalence (st1 (1),profle(1))
c
c
c
c 1.0 Initialization
c
c 1.1 Return if no plot is required
c
      if(iplotm .lt. 6) return
c
c
c 1.2 Initialize miscellaneous parameters
c
c 1.2.1 Roundoff for machine precision
c
      prndff   = roundff
c
c 1.2.2 Default parameters for pointc
c
      numbpc   = 1
      incpcx   = 0
      incpcy   = 0
      delpcx   = 0.0
      delpcy   = 0.0
c
c
c
c 2.0 Set up map and initialize parameters
c
c 2.1 Normalized plot dimensions
c
      smn      = 0.00
      smx      = 1.00
      smn0     = 0.00
      smx0     = 1.00
c
c
c 2.2 Physical plot dimensions
c
c 2.2.1 Plot borders
c
      x0mn     = 0.120
      x0mx     = 0.820
      y0mn     = 0.200
      y0mx     = 0.900
c
c 2.2.2 Heading position
c
      xxhd     = 0.125
      yyhd     = 0.95
c
c 2.3 Label position parameters
c
c 2.3.1 Position increments
c
      nspx0    = 16
      nspy0    = 16
      nspz0    = 16
c
c 2.3.2 Axis label position fraction
c
      msps0    = (2*nspx0)/3
      mspq0    = (3*nspy0)/4
      mspp0    = (3*nspz0)/4
      mspd0    = (5*nspz0)/12
      mspf0    = (1*nspz0)/5
c
c
c 2.4 Set the normalized tick parameters
c
c 2.4.1 Tick scale skipping
c
      ntikxs0  = 2
      ntikys0  = 2
      ntikzs0  = 2
c
c 2.4.2 Normalized tick lengths
c
      xtiky0   = +0.250
      ytikx0   = +0.250
      ztikx0   = -0.250
c
c 2.4.3 Tick scale positions
c
c 2.4.3.1 Label shift along relevant axis
      xsclx0   = +0.250
      yscly0   = -0.250
      zsclz0   = -0.250
c
c 2.4.3.2 Label shift away from relevant axis
      xscly0   = -2.350
      ysclx0   = -2.100
      zsclx0   = +0.250
c
c
c 2.5 Set the normalized labels and positions
c
c 2.5.1 Horizontal axis
c
      slby10   = -2.850
      slabel   = pclab
c
c 2.5.2 Left  Vertical axis
c
      qlbx10   = -2.700
      qlabel   = ' q '
c
c 2.5.3 Right Vertical axis
c
      plbx10   = +0.200
      dlbx10   = +0.200
      flbx10   = +0.200
c
      plabel   = 'pp/ppmax'
      dlabel   = 'n/nmax'
      flabel   = 'pf/pfmax'
c
c
c 2.6 Set the character printing parameters
c
c 2.6.1 Heading: case, size, and orientation
c
      icashd   = +2
      isizhd   = +3
      iornhd   =  0
c
c 2.6.2 Horizontal axis
c
c 2.6.2.1 Horizontal axis scale: case, size, and orientation
      icshsc   = +1
      iszhsc   = +1
      iorhsc   = +1
c
c 2.6.2.2 Horizontal axis label: case, size, and orientation
      icshlb   = +1
      iszhlb   = +3
      iorhlb   =  0
c
c 2.6.3 Left Vertical axis
c
c 2.6.3.1 Left Vertical axis scale: case, size, and orientation
      icsvsy   = +1
      iszvsy   = +1
      iorvsy   =  0
c
c 2.6.3.2 Left Vertical axis label: case, size, and orientation
      icslbq   = +1
      iszlbq   = +3
      iorlbq   =  0
c
c 2.6.4 Right Vertical axis
c
c 2.6.4.1 Right Vertical axis scale: case, size, and orientation
      icsvsz   = +1
      iszvsz   = +1
      iorvsz   =  0
c
c 2.6.4.2 Right Vertical axis label: case, size, and orientation
      icslbp   = +1
      iszlbp   = +2
      iorlbp   =  0
c
      icslbd   = +1
      iszlbd   = +2
      iorlbd   =  0
c
      icslbf   = +1
      iszlbf   = +2
      iorlbf   =  0
c
c
c 2.7 Line type parameters
c
c 2.7.1 Identifying characters
c
      icschp   =  2
      iszchp   = +1
      iorchp   =  0
      kspchp   =  2
      chpp     = 'p'
c
      icschd   =  2
      iszchd   = +1
      iorchd   =  0
      kspchd   =  2
      chdn     = 'n'
c
      icschf   =  2
      iszchf   = +1
      iorchf   =  0
      kspchf   =  2
      chpf     = 'f'
c
c
c
c 3.0 Find the function minima and maxima
c
c 3.1 Search the physically scaled data
c     Note that psival, qp, dnsty, and pp are all stored in
c     reverse order
c
      spmin    = psivmsh(  1  )
      spmax    = psivmsh(jpsi2)
      qpmin    =  0.00
      qpmax    = -1.0e+10
      do 10  j = 1,jpsi1
      qpmin    = amin1(qpmin,qp   (j))
      qpmax    = amax1(qpmax,qp   (j))
   10 continue
c
      qpmin    = amin1(qpmin,qaxe)
      qpmax    = amax1(qpmax,qaxe)
c
c
c 3.2 Adjust the maxima as needed
c
c 3.2.1 Physical  q(psi)  profile maximum
c
      iqmx     = ifix (qpmax) + 1
      qpmax    = float(iqmx)
c
      prmin    = 0.0
      prmax    = 1.20
c
c
c
c 4.0 Set the physical plot parameters
c
c 4.1 Set the number of tick marks
c
c 4.1.1 Horizontal axis
c
      nspx     = 10
c
c 4.1.2 Left  Vertical axis
c
      nspy1    = ifix(qpmax)
      nspy1    = min0(nspy1,nspy0)
      nspy1    = max0(nspy1,1)
c
      if    (nspy1 .le. 4) then
        if    (nspy1 .le. 1) then
          nspy   = 5*nspy1
        elseif(nspy1 .le. 2) then
          nspy   = 4*nspy1
        else
          nspy   = 2*nspy1
        endif
c
      elseif(nspy1 .gt. 4) then
        nspy   =   nspy1
      endif
c
c 4.1.3 Right Vertical axis
c
      zspz     = 10.0*prmax + prndff
      nspz     = ifix(zspz)
c
c
c 4.2 Normalize the physical spacing increments to the physical range
c
c 4.2.1 Standard spacing
c
      delx0    = (spmax-spmin)/nspx0
      dely0    = (qpmax-qpmin)/nspy0
      delz0    = (prmax-prmin)/nspz0
c
c 4.2.2 Tick spacing
c
      delx     = (spmax-spmin)/nspx
      dely     = (qpmax-qpmin)/nspy
      delz     = (prmax-prmin)/nspz
c
c
c 4.3 Set the physical spacings
c
c 4.3.1 Tick scale values
c
      ntikxs   = min0(ntikxs0,nspx)
      ntikxs   = max0(ntikxs,  1  )
      if(nspx .lt. 5) ntikxs  = 1
c
      ntikys   = min0(ntikys0,nspy)
      ntikys   = max0(ntikys,  1  )
      if(nspy .lt. 5) ntikys  = 1
c
      ntikzs   = min0(ntikzs0,nspz)
      ntikzs   = max0(ntikzs,  1  )
      if(nspz .lt. 5) ntikzs  = 1
c
c
c 4.3.2 Axis scale positions
c
c 4.3.2.1 Label shift along relevant axis
      xsclx    = xsclx0*delx0
      yscly    = yscly0*dely0
      zsclz    = zsclz0*delz0
c
c 4.3.2.2 Label shift away from relevant axis
      xscly    = xscly0*dely0
      ysclx    = ysclx0*delx0
      zsclx    = zsclx0*delx0
c
c 4.3.3 Tick lengths
c
      xtiky    = xtiky0*dely0
      ytikx    = ytikx0*delx0
      ztikx    = ztikx0*delx0
c
c 4.3.3 Label positions
c
      slby1    = slby10*dely0
      qlbx1    = qlbx10*delx0
      plbx1    = plbx10*delx0
      dlbx1    = dlbx10*delx0
      flbx1    = flbx10*delx0
c
c
c
c 5.0 Print out the plot heading
c
c 5.1 Set up the plot in normalized coordinates
c
      call map(smn,smx,smn,smx,smn0,smx0,smn0,smx0)
c
c
c 5.2 Set the position and plot the text
c
      call setlch(xxhd,  yyhd,  icashd,isizhd,iornhd,-1)
      call gtext('Plot of q pp n and pfast vs psi',31,-1)
c
c
c
c 6.0 Plot q(psi)
c
c 6.1 Renormalize the map for the  q(psi) profile and plot the borders
c
c 6.1.1 Set the plot page in physical units
c
      call map(spmin,spmax,qpmin,qpmax,x0mn,x0mx,y0mn,y0mx)
c
c 6.1.2 Plot the four borders in physical units
c
      call line(spmin,qpmin,spmin,qpmax)
      call line(spmin,qpmax,spmax,qpmax)
      call line(spmax,qpmax,spmax,qpmin)
      call line(spmin,qpmin,spmax,qpmin)
c
c
c 6.2 Set up the Horizontal axis
c
c 6.2.1 Plot the horizontal axis tick marks and scale
c
      qpos     = qpmin
      qpos1    = qpos
      qpos2    = qpos + xtiky
      qpossc   = qpos + xscly
      do 100 k = 0,nspx
      kskpp    = (k/ntikxs)*ntikxs
      spos     = spmin + k*delx
      spossc   = spos  + xsclx
c
c 6.2.1.1 Plot the tick mark
      call line  (spos,  qpos1, spos,  qpos2 )
c
c 6.2.1.2 Plot the scale value
      if(k .eq. kskpp) then
        call setlch(spossc,qpossc,icshsc,iszhsc,iorhsc,-1)
        write(string,1000) spos
        call wrtstr(string,1)
      endif
  100 continue
c
c 6.2.2 Plot the horizontal axis label
c
      sposlb   = spmin  + msps0*delx0
      qposlb   = qpmin  + slby1
      call setlch(sposlb,qposlb,icshlb,iszhlb,iorhlb,-1)
      write(string,1010) slabel
      call wrtstr(string,1)
c
c
c 6.3 Set up the Left  vertical axis
c
c 6.3.1 Plot the  vertical  axis tick marks and scale
c
      spos     = spmin
      spos1    = spos
      spos2    = spos + ytikx
      spossc   = spos + ysclx
      do 110 i = 0,nspy
      iskpp    = (i/ntikys)*ntikys
      qpos     = qpmin + i*dely
      qpossc   = qpos  + yscly
c
c 6.3.1.1 Plot the tick mark
      call line  (spos1, qpos,  spos2, qpos  )
c
c 6.3.1.2 Plot the scale value
      if(i .eq. iskpp) then
        call setlch(spossc,qpossc,icsvsy,iszvsy,iorvsy,-1)
        write(string,1100) qpos
        call wrtstr(string,1)
      endif
  110 continue
c
c 6.3.2 Plot the  vertical  axis label
c
      sposlb   = spmin  + qlbx1
      qposlb   = qpmin  + mspq0*dely0
      call setlch(sposlb,qposlb,icslbq,iszlbq,iorlbq,-1)
      write(string,1110) qlabel
      call wrtstr(string,1)
c
c
c 6.4 Plot the  q(psi)  profile
c
c 6.4.1 Set the profile
c
      do 120 j   = 1,jpsi2
      jp         = j
      qvpt       = qp1    (jp)
      profle(jp) = qvpt
  120 continue
c
c 6.4.2 Plot the curve
c
      call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
c
c
c
c 7.0 Plot the normalied profiles
c
c 7.1 Renormalize the map for the normalized profiles
c
      call map(spmin,spmax,prmin,prmax,x0mn,x0mx,y0mn,y0mx)
c
c
c 7.2 Set up the Right Vertical axis
c
c 7.2.1 Plot the  vertical  axis tick marks and scale
c
      spos     = spmax
      spos1    = spos + ztikx
      spos2    = spos
      spossc   = spos + zsclx
      do 150 i = 0,nspz
      iskpp    = (i/ntikzs)*ntikzs
      ppos     = prmin + i*delz
      ppossc   = ppos  + zsclz
c
c 7.2.1.1 Plot the tick mark
      call line  (spos1, ppos,  spos2, ppos  )
c
c 7.2.1.2 Plot the scale value
      if(i .eq. iskpp) then
        call setlch(spossc,ppossc,icsvsz,iszvsz,iorvsz,-1)
        write(string,1200) ppos
        call wrtstr(string,1)
      endif
  150 continue
c
c
c 7.3 Set up the physical plot ranges for normalizing the profiles
c
c 7.3.1 Find the maximum pressure gradient
c 
      ppmax    = ppaxe
      ppmin    = 0.0
      jppmax   = jpsi2
      do 200 j = 1,jpsi1
      if(abs(pp   (j)) .gt. abs(ppmax)) then
        ppmax    = pp   (j)
        jppmax   = j
      endif
  200 continue
c
      if(ppmax .eq. 0.0) ppmax  = 1.0
c
c 7.3.2 Find the maximum density
c
      dnmax    = dnaxe
      dnmin    = 0.0
      jdnmax   = jpsi2
      do 210 j = 1,jpsi1
      if(abs(dnsty(j)) .gt. abs(dnmax)) then
        dnmax    = dnsty(j)
        jdnmax   = j
      endif
  210 continue
c
      if(dnmax .eq. 0.0) dnmax  = 1.0
c
c 7.3.3 Find the maximum and minimum fast ion pressure
c
      pfmin    = pfaxe
      pfminj   = jpsi2
      pfmax    = pfaxe
      jpfmax   = jpsi2
      do 220 j = 1,jpsi1
      if(abs(pfast (j)) .lt. abs(pfmin)) then
        pfmin    = pfast(j)
        pfminj   = j
      endif
      if(abs(pfast (j)) .gt. abs(pfmax)) then
        pfmax    = pfast(j)
        jpfmax   = j
      endif
  220 continue
c
      if(pfmin .eq. 0.0  .and.  pfmax .eq. 0.0) nofplot = 0
      if(pfmin .ne. 0.0  .or.   pfmax .ne. 0.0) nofplot = 1
c
      if(pfmax .eq. 0.0) pfmax  = 1.0
c
c
c 7.4 Plot the  vertical  axis labels
c
c 7.4.1 Pressure gradient label
c
      sposlp   = spmax  + plbx1
      pposlp   = prmin  + mspp0*delz0
      call setlch(sposlp,pposlp,icslbp,iszlbp,iorlbp,-1)
      write(string,1210) plabel
      call wrtstr(string,1)
      write(string,1220) abs(ppmax)
      call wrtstr(string,1)
      write(string,1230) jppmax
      call wrtstr(string,1)
c
c 7.4.2 Density profile label
c
      sposld   = spmax  + dlbx1
      pposld   = prmin  + mspd0*delz0
      call setlch(sposld,pposld,icslbd,iszlbd,iorlbd,-1)
      write(string,1310) dlabel
      call wrtstr(string,1)
      write(string,1320) abs(dnmax)
      call wrtstr(string,1)
      write(string,1330) jdnmax
      call wrtstr(string,1)
c
c 7.4.3 Fast particle pressure profile label if profile is nonzero
c
      if(nofplot .ne. 0) then
        sposlf   = spmax  + flbx1
        pposlf   = prmin  + mspf0*delz0
        call setlch(sposlf,pposlf,icslbf,iszlbf,iorlbf,-1)
        write(string,1410) flabel
        call wrtstr(string,1)
        write(string,1420) abs(pfmax)
        call wrtstr(string,1)
        write(string,1430) jpfmax
        call wrtstr(string,1)
      endif
c
c
c 7.5 Plot the normalized profiles
c
c 7.5.1 Plot the pressure gradient profile
c
c 7.5.1.2 Reset the printing parameters for pointc
      call setpch(icschp,iszchp,iorchp,kspchp)
c
c 7.5.1.3 Plot the normalized values
c
c 7.5.1.1 Set the profile and plot the points
      do 250 j   = 1,jpsi2
      jp         = j
      jr         = jpsi2 - j
      sppt       =     psivmsh(jp)
      if(jr .eq. 0) pppt     = abs(paxe      /ppmax)
      if(jr .gt. 0) pppt     = abs(pp    (jr)/ppmax)
      profle(jp) = pppt
      call pointc(chpp,sppt,pppt,numbpc,incpcx,incpcy,delpcx,delpcy)
  250 continue
c
c 7.5.1.2 Plot the curve
      call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
c
c 7.5.2 Plot the density profile
c
c 7.5.2.2 Reset the printing parameters for pointc
      call setpch(icschd,iszchd,iorchd,kspchd)
c
c 7.5.2.3 Plot the normalized values
c
c 7.5.2.3.1 Set the profile and plot the points
      do 260 j   = 1,jpsi2
      jp         = j
      jr         = jpsi2 - j
      sppt       =  psivmsh(jp)
      if(jr .eq. 0) dnpt     = abs(dnaxe     /dnmax)
      if(jr .gt. 0) dnpt     = abs(dnsty (jr)/dnmax)
      profle(jp) = dnpt
      call pointc(chdn,sppt,dnpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  260 continue
c
c 7.5.2.3.2 Plot the curve
      call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
c
c 7.5.3 Plot the fast particle pressure profile if nonzero
c
c 7.5.3.2 Reset the printing parameters for pointc
      if(nofplot .ne. 0) then
        call setpch(icschf,iszchf,iorchf,kspchf)
c
c 7.5.3.3 Plot the normalized values
c
c 7.5.3.3.1 Set the profile and plot the points
        do 270 j   = 1,jpsi2
        jp         = j
        jr         = jpsi2 - j
        sppt       =  psivmsh(jp)
        if(jr .eq. 0) fppt     = abs(pfaxe     /pfmax)
        if(jr .gt. 0) fppt     = abs(pfast (jr)/pfmax)
        profle(jp) = fppt
        call pointc(chpf,sppt,fppt,numbpc,incpcx,incpcy,delpcx,delpcy)
  270   continue
c
c 7.5.3.3.2 Plot the curve
        call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
      endif
c
c
c
c 8.0 Close frame, return, and end
c
      call frame(0)
      return
c
 1000 format(f8.4)
 1010 format(a8)
 1100 format(f7.2)
 1110 format(a3)
 1200 format(f4.2)
 1210 format(a6)
 1220 format(1pe8.2)
 1230 format('j=',i3)
 1310 format(a6)
 1320 format(1pe8.2)
 1330 format('j=',i3)
 1410 format(a6)
 1420 format(1pe8.2)
 1430 format('j=',i3)
      end
      subroutine plotint
c
c ----------------------------------------------------
c plot flux surface integrals
c ----------------------------------------------------
c
      parameter (npx=200,ncx=2*npx,np1=npx+1,np2=npx+2,nc1=ncx+1)
      parameter (npp=929,ntt=2*npp-1,npb=npp)
      parameter (npq=np2)
      parameter (nvn=7)
c
      character*1   chs0,  chs1,  chs2,  chs3,  chs4,
     &              chsv,  chtv,  chdq,  chdl
      character*8   slabel,s0labl,s1labl,s2labl,s3labl,s4labl,
     &                     svlabl,tvlabl,dqlabl,dllabl
      character*80  string
c
      character*1   lbchpv, lbchsp, lbchrh, lbchro,
     &              lbchtf, lbchst ,lbchvl, lbchsv
      character*8   pclab,  tclab,  labpsi, 
     &              labpsiv,labpsin,labsqtp,
     &              labrhov,labrho2,labrhon,
     &              labtorf,labtorn,labsqtt,
     &              labvolm,labvoln,labsqtv,
     &              labdrh2,labdrho,labdtor,
     &              labdsqt,labdvol,labdsqv
c
      common/cnst/  pi,hlfpi,twopi,amu,boltzm,echarg,pmass,roundff,bigno
      common/inpd/  ntor,ncase,idnsty,ndnxp0,ndnxp1,ndnxp2,ifastp,
     &              nmap,neqtyp,mpreset,ndoublt,ndivert,nmtype,
     &              jpsi,itht,jpsi1,jpsi2,ithtp,jsep,
     &              isym,igrid,nham1,nham2,nham3,npowr,
     &              bfieldf,rdefolt,qxin,btdes,qsurf,
     &              gamma,gamav,rmantl,betaf,zeffect,
     &              ncorr,corrfac,nccellr,peqpk0,peqpk1,peqpk2
      common/mapd/  mapmaxd,dpsisl,dpsisd,nqaxis,nwtmag,nfitmax,nfitpts,
     &              ifitrng,jfitrng,jfitchk,fitchek,cnvmag,epsaxs,
     &              maxerlp,maxerlv,delbox,delboz,delac,delav,
     &              delstsf,delstlp,delstlv,prfrac,
     &              nerstop,nerprnt,qptol,tolspln,tolbchi,tolbtor,
     &              tolsymm,tolaugm,errsep,precisn,plossmx,
     &              narcmx,ntrymx,ntdecr,ntmmin,npfit,npcmin,
     &              kuttaop,nrkmax0,nrkmax1,numstp,nwtfitp,nwtfitm,
     &              stepfac,flxstp,psispl,tolnwtp,tolnwtm,
     &              delpakf,delpakc,delpkf,delpkc,psichek,boxtnd,
     &              maptrace,norient,maxcutc,
     &              dresolv,dlclose,pntshft,endtol,
     &              narcln,nangax,nanglm,nbpmax,nwtmax,nslmax,nhfmax,
     &              bperor,sersnm,sertnm,arcmin,delgap,stepcut,
     &              novap,ioutm,iplotm
      common/dblt/  psisep,xsep,zsep,smap(2),smax(2),smaz(2),numax,
     &              psimx(2),xax(2),zax(2),fitax(14,2),fitsep(14)
      common/ort1/  ndim1,ndim2,xma,zma,rcnt,btor,totcur,
     &              psimax0,psilim0,delpsir0,psimax1,psilim1,delpsir1,
     &              psimax,psilim,delpsir,dpsidec,dpsisv,dpedg1,dpedg2,
     &              eaxe,taxe,qaxe,fqiaxe,fqmaxe,paxe,faxe,ppaxe,ffpaxe,
     &              dnaxe,pfaxe,qppaxe,dlnaxe,dpnaxe,dhnaxe,
     &              arcaxe,dnnorm,pfnorm,
     &              qlim,fqilim,fqmlim,plim,flim,pplim,ffplim,
     &              dnlim,pflim,xlim,zlim,
     &              qmer,btmer,gkappa,arcprev,chicmax,
     &              charcmx,chpstmx,chhammx,arcnrm0,pstnrm0,hamnrm0,
     &              elongax1,elongax2,elongaxs,elongaxp,elongtn(npp)
      common/ort2/  psival(np1),pp(np1),ffp(np1),p(np1),f(np1),
     &              dnsty(np1),pfast(np1),qp(np1),fqpi(np1),fqim(np1),
     &              qdel(np1),psinitl(np1),qpinitl(np1),
     &              qprime(np1),qpint(np1),fqint(np1),
     &              psivl0(np2),rst(np2),zst(np2),rstp(np2),zstp(np2),
     &              rsinitl(np2),zsinitl(np2),
     &              arcsurf(nc1),pestchi(nc1),hamlchi(nc1),
     &              alfarc (nc1),alfpst (nc1),alfham (nc1),
     &              rs(nc1,np1),zs(nc1,np1)
      common/ort3/  rsrf(nc1),zsrf(nc1),chic(nc1),chie(nc1),dydx(nc1),
     &              xjsc(nc1)
      common/ort4/  rcell  (np1,nc1),zcell  (np1,nc1),
     &              dpsdr  (np1,nc1),dpsdz  (np1,nc1),
     &              chiarcl(np1,nc1),chipest(np1,nc1),chihaml(np1,nc1),
     &              xjcarcl(np1,nc1),xjcpest(np1,nc1),xjchaml(np1,nc1),
     &              alfarcl(np1,nc1),alfpest(np1,nc1),alfhaml(np1,nc1),
     &              arcnorm(np1),    pstnorm(np1),    hamnorm(np1),
     &              dlnlval(np1),    dpnlval(np1),    dhnlval(np1)
      common/ort5/  f3 (np1,nc1),f4 (np1,nc1),f5 (np1,nc1),
     &              f7 (np1,nc1),f8 (np1,nc1),f9 (np1,nc1),
     &              f10(np1,nc1),f11(np1,nc1),f12(np1,nc1),
     &              f13(np1,nc1),f14(np1,nc1),f15(np1,nc1),
     &              f16(np1,nc1),f17(np1,nc1),f18(np1,nc1),
     &              f19(np1,nc1),f20(np1,nc1),f21(np1,nc1),
     &              f22(np1,nc1),f23(np1,nc1),f24(np1,nc1),
     &              f25(np1,nc1),f26(np1,nc1),f27(np1,nc1),
     &              f28(np1,nc1)
      common/ort6/  savgax,tavgax,sntax0,sntax1,sntax2,sntax3,sntax4,
     &              svint(nvn),savge(np1),tavge(np1),sint0(np1),
     &              sint1(np1),sint2(np1),sint3(np1),sint4(np1)
      common/auxmsh/rh2tot,rhotot,flxtot,voltot,
     &              anltflxp,dvrtflxp,anltvolp,dvrtvolp,
     &              analtflx,divrtflx,analtvol,divrtvol,
     &              psivl1 (np2),qp1    (np2),fqpi1  (np2),
     &              qprim1 (np2),sintp0 (np2),psivmsh(np2),
     &              psivalu(np2),psinorm(np2),psisqrt(np2),
     &              psirho2(np2),psirhov(np2),psinrho(np2),
     &              psitorf(np2),psintor(np2),psisqtf(np2),
     &              psivolm(np2),psinvol(np2),psisqvl(np2),
     &              dpsirh2(np2),dpsirho(np2),dpsinrh(np2),
     &              dpsitor(np2),dpsintf(np2),dpsisqt(np2),
     &              dpsivol(np2),dpsinvl(np2),dpsisqv(np2)
       common/ratnl/jqmin, jqmax, numbqs,nq1totl,
     &              rcentr,zcentr,rminor,zminor,
     &              asprat00,asprat01,asprat10,asprat11,
     &              pminvl,qminvl,qpminv,pmaxvl,qmaxvl,qpmaxv,
     &              psivlmnq,    psivnmnq,    psisqmnq,
     &              psir2mnq,    psirhmnq,    psinrmnq,
     &              psitfmnq,    psintmnq,    psistmnq,
     &              psivmmnq,    psinvmnq,    psisvmnq,
     &              dpsr2mnq,    dpsrhmnq,    dpsnrmnq,
     &              dpstfmnq,    dpsntmnq,    dpsstmnq,
     &              dpsvmmnq,    dpsnvmnq,    dpssvmnq,
     &              psivlmxq,    psivnmxq,    psisqmxq,
     &              psir2mxq,    psirhmxq,    psinrmxq,
     &              psitfmxq,    psintmxq,    psistmxq,
     &              psivmmxq,    psinvmxq,    psisvmxq,
     &              dpsr2mxq,    dpsrhmxq,    dpsnrmxq,
     &              dpstfmxq,    dpsntmxq,    dpsstmxq,
     &              dpsvmmxq,    dpsnvmxq,    dpssvmxq,
     &              lpsiq  (npq),psivlq (npq),
     &              qprimq (npq),qvalue (npq),psimshq(npq),
     &              psivalq(npq),psinrmq(npq),psisqrq(npq),
     &              psirh2q(npq),psirhoq(npq),psinrhq(npq),
     &              psitorq(npq),psintfq(npq),psisqtq(npq),
     &              psivolq(npq),psinvlq(npq),psisqvq(npq),
     &              dpsir2q(npq),dpsirhq(npq),dpsinrq(npq),
     &              dpsitfq(npq),dpsintq(npq),dpsistq(npq),
     &              dpsivlq(npq),dpsinvq(npq),dpsisvq(npq),
     &              shearps(npq),shearrh(npq),
     &              sheartf(npq),shearvl(npq),
     &              epslrh1 (npq),shearrh1(npq),shearfrh(npq),
     &              epslvl1 (npq),shearvl1(npq),shearfvl(npq)
      common/labels/lbchpv, lbchsp, lbchrh, lbchro,
     &              lbchtf, lbchst, lbchvl, lbchsv,
     &              pclab,  tclab,  labpsi,
     &              labpsiv,labpsin,labsqtp,
     &              labrhov,labrho2,labrhon,
     &              labtorf,labtorn,labsqtt,
     &              labvolm,labvoln,labsqtv,
     &              labdrh2,labdrho,labdtor,
     &              labdsqt,labdvol,labdsqv
      common/pldf/  x0min,x0max,y0min,y0max
      common/inout/ kutty,kuout,kuinp,kueql,kumap,kugrd,kurdm,kurdn,
     &              kufst,krdeqlb,krdmgta,krdngta,krdfast
c
c
c
c 1.0 Initialization
c
c 1.1 Return if no plot is required
c
      if(iplotm .lt. 6) return
c
c
c 1.2 Initialize miscellaneous parameters
c
c 1.2.1 Roundoff for machine precision
c
      prndff   = roundff
c
c 1.2.2 Default parameters for lines
c
      kskp0    = 2
      kskp1    = 3
      kskp2    = 2
      kskp3    = 3
      zero     = 0.0
c
c 1.2.3 Default parameters for pointc
c
      numbpc   = 1
      incpcx   = 0
      incpcy   = 0
      delpcx   = 0.0
      delpcy   = 0.0
c
c
c 1.3 Vertical axis position increments
c
      nspx0    = 16
      nspy0    = 16
      nspz0    = 16
c
c
c
c 2.0 Set up map and initialize parameters
c
c 2.1 Normalized plot dimensions
c
      smn      = 0.00
      smx      = 1.00
      smn0     = 0.00
      smx0     = 1.00
c
c
c 2.2 Set physical plot dimensions
c
c 2.2.1 Set plot borders
c
      x0mn     = 0.120
      x0mx     = 0.820
      y0mn     = 0.200
      y0mx     = 0.900
c
c 2.2.2 Heading position
c
      xxhd     = 0.225
      yyhd     = 0.95
c
c
c 2.3 Set the character printing parameters
c
c 2.3.1 Heading: case, size, and orientation
c
      icashd   = +2
      isizhd   = +3
      iornhd   =  0
c
c
c
c 3.0 Set up data for the horizontal axis
c
c 3.1 Set the normalized tick parameters
c
c 3.1.1 Set the number of horizontal axis tick marks
c
      nspx     = 10
c
c 3.1.2 Tick scale skipping
c
      ntikxsm  = 2
      ntikysm  = 2
      ntikzsm  = 2
c
c 3.1.3 Tick scale values
c
      ntikxs   = min0(ntikxsm,nspx)
      ntikxs   = max0(ntikxs,  1  )
      if(nspx .lt. 5) ntikxs  = 1
c
c 3.1.4 Normalized tick lengths
c
      xtiky0   = +0.250
      ytikx0   = +0.250
      ztikx0   = -0.250
c
c 3.1.5 Tick scale positions
c
c 3.1.5.1 Label shift along relevant axis
      xsclx0   = +0.250
      yscly0   = -0.250
      zsclz0   = -0.250
c
c 3.1.5.2 Label shift away from relevant axis
      xscly0   = -2.350
      ysclx0   = -2.250
      zsclx0   = +0.225
c
c
c 3.2 Set the horizontal axis label parameters
c
c 3.2.1 Horizontal axis
c
      slby10   = -2.850
      slabel   = pclab
c
c 3.2.2 Horizontal axis: case, size, and orientation
c
c 3.2.2.1 Scale: case, size, and orientation
      icshsc   = +1
      iszhsc   = +1
      iorhsc   = +1
c
c 3.2.2.2 Label: case, size, and orientation
      icshlb   = +1
      iszhlb   = +3
      iorhlb   =  0
c
c 3.2.3 Label position fraction
c
      mspps    = (2*nspx0)/3
c
c
c 3.3 Set the minimum and maximum
c
      spmin    = psivmsh(  1  )
      spmax    = psivmsh(jpsi1)
c
c
c 3.4 Set the vertical axis scale character parameters
c
c 3.4.1 Left vertical axis: case, size, and orientation
c
      icsvsy   = +1
      iszvsy   = +1
      iorvsy   =  0
c
c 3.4.2 Right vertical axis: case, size, and orientation
c
      icsvsz   = +1
      iszvsz   = +1
      iorvsz   =  0
c
c
c
c 4.0 Plot the surface integrals
c
c 4.1 Set the left vertical axis parameters
c
c 4.1.1 Left axis label position fraction
c
      mss00    = (5*nspy0)/6
c
c 4.1.2 Set the normalized labels and positions
c
      s0lbx10  = -2.70
      s0labl   = 's0'
c
c 4.1.3 Left vertical axis label: case, size, and orientation
c
      icslbs0  = +1
      iszlbs0  = +3
      iorlbs0  =  0
c
c 4.1.4 Left vertical axis: specific line type parameters 
c
      icschs0  =  0
      iszchs0  = +1
      iorchs0  =  0
      chs0     = '0'
c
c
c 4.2 Set the right vertical axis parameters
c
c 4.2.1 Right axis label position fraction
c
      mss10    = (7*nspz0)/8
      mss20    = (4*nspz0)/5
      mss30    = (2*nspz0)/3
      mss40    = (1*nspz0)/3
c
c 4.2.2 Set the normalized labels and positions
c
      s1lbx10  = +1.050
      s2lbx10  = +1.050
      s3lbx10  = +1.050
      s4lbx10  = +1.050
c
      s1labl   = 's1/s0'
      s2labl   = 's2/s0'
      s3labl   = 's3/s0'
      s4labl   = 's4/s0'
c
c 4.2.3 Right vertical axis label: case, size, and orientation
c
      icslbs1  = +2
      iszlbs1  = +3
      iorlbs1  =  0
c
      icslbs2  = +2
      iszlbs2  = +3
      iorlbs2  =  0
c
      icslbs3  = +2
      iszlbs3  = +3
      iorlbs3  =  0
c
      icslbs4  = +2
      iszlbs4  = +3
      iorlbs4  =  0
c
c 4.2.4 Right vertical axis: specific line type parameters
c
      icschs1  = +2
      iszchs1  = +1
      iorchs1  =  0
      kspchs1  =  2
      chs1     = '1'
c
      icschs2  = +2
      iszchs2  = +1
      iorchs2  =  0
      kspchs2  =  2
      chs2     = '2'
c
      icschs3  = +2
      iszchs3  = +1
      iorchs3  =  0
      kspchs3  =  2
      chs3     = '3'
c
      icschs4  = +2
      iszchs4  = +1
      iorchs4  =  0
      kspchs4  =  2
      chs4     = '4'
c
c
c 4.3 Set the physical plot parameters
c
c 4.3.1 Find the function minima and maxima
c
c 4.3.1.1 Search the physically scaled data
c         Note that psivmsh and the integrals are stored in forward order
      s0min    =  0.00
      s0max    = -bigno
      do 10  j = 1,jpsi1
      s0min    = amin1(s0min,sint0(j))
      s0max    = amax1(s0max,sint0(j))
   10 continue
c
c 4.3.1.2 Adjust the maxima as needed
      is0mx    = ifix (s0max) + 1
      s0max    = float(is0mx)
c
c 4.3.1.3 Set the left and right axis extrema
c       Set the left axis quantities normalized to 1 with 20% headroom
      slftmin  = s0min
      slftmax  = s0max
      if(slftmin .gt. zero) slftmin = zero
      if(slftmax .lt. zero) slftmax = zero
c
      srgtmin  = 0.0
      srgtmax  = 1.20
c
c 4.3.2 Set the number of tick marks
c
c 4.3.2.1 Left  Vertical axis
      nspy1    = ifix(slftmax)
      nspy1    = min0(nspy1,nspy0)
      nspy1    = max0(nspy1,1)
      if    (nspy1 .le. 4) then
        if    (nspy1 .le. 1) then
          nspy   = 5*nspy1
        elseif(nspy1 .le. 2) then
          nspy   = 4*nspy1
        else
          nspy   = 2*nspy1
        endif
      elseif(nspy1 .gt. 4) then
        nspy   =   nspy1
      endif
c
c 4.3.2.2 Right Vertical axis
      nspz1    = ifix(srgtmax)
      nspz1    = min0(nspz1,nspz0)
      nspz1    = max0(nspz1,1)
      if    (nspz1 .le. 4) then
        if    (nspz1 .le. 1) then
          nspz   = 5*nspz1
        elseif(nspz1 .le. 2) then
          nspz   = 4*nspz1
        else
          nspz   = 2*nspz1
        endif
      elseif(nspz1 .gt. 4) then
        nspz   =   nspz1
      endif
c
c 4.3.3 Normalize the physical spacing increments to the physical range
c
c 4.3.3.1 Standard spacing
      delx0    = (spmax-spmin)/nspx0
      dely0    = (slftmax-slftmin)/nspy0
      delz0    = (srgtmax-srgtmin)/nspz0
c
c 4.3.3.2 Tick spacing
      delx     = (spmax-spmin)/nspx
      dely     = (slftmax-slftmin)/nspy
      delz     = (srgtmax-srgtmin)/nspz
c
c 4.3.3.3 Tick scale values for vertical axes
      ntikys   = min0(ntikysm,nspy)
      ntikys   = max0(ntikys,  1  )
      if(nspy .lt. 5) ntikys  = 1
c
      ntikzs   = min0(ntikzsm,nspz)
      ntikzs   = max0(ntikzs,  1  )
      if(nspz .lt. 5) ntikzs  = 1
c
c 4.3.3.4 Tick lengths
      xtiky    = xtiky0*dely0
      ytikx    = ytikx0*delx0
      ztikx    = ztikx0*delx0
c
c 4.3.4 Axis scale positions
c
c 4.3.4.1 Label shift along relevant axis
      xsclx    = xsclx0*delx0
      yscly    = yscly0*dely0
      zsclz    = zsclz0*delz0
c
c 4.3.4.2 Label shift away from relevant axis
      xscly    = xscly0*dely0
      ysclx    = ysclx0*delx0
      zsclx    = zsclx0*delx0
c
c 4.3.5 Axis label positions
c
c 4.3.5.1 Horizontal axis label
      slby1    = slby10*dely0
c
c 4.3.5.2 Left vertical axis label
      s0lbx1   = s0lbx10*delx0
c
c 4.3.5.3 Right vertical axis labels
      s1lbx1   = s1lbx10*delx0
      s2lbx1   = s2lbx10*delx0
      s3lbx1   = s3lbx10*delx0
      s4lbx1   = s4lbx10*delx0
c
c
c 4.4 Print out the plot heading
c
c 4.4.1 Set up the plot in normalized coordinates
c
      call map(smn,smx,smn,smx,smn0,smx0,smn0,smx0)
c
c 4.4.2 Set the position and plot the text
c
      call setlch(xxhd,  yyhd,  icashd,isizhd,iornhd,-1)
      call gtext('Plot of surface integrals',25,-1)
c
c
c 4.5 Plot the left vertical axis functions
c
c 4.5.1 Renormalize the map in physical units
c
      call map(spmin,spmax,slftmin,slftmax,x0mn,x0mx,y0mn,y0mx)
c
c 4.5.2 Plot the four borders in physical units
c
      call line(spmin,slftmin,spmin,slftmax)
      call line(spmin,slftmax,spmax,slftmax)
      call line(spmax,slftmax,spmax,slftmin)
      call line(spmin,slftmin,spmax,slftmin)
c
c 4.5.3 Set up the Horizontal axis
c
c 4.5.3.1 Plot the horizontal axis tick marks and scale
c
c 4.5.3.1.1 Find the tick positions
      slftpos  = slftmin
      slftpos1 = slftpos
      slftpos2 = slftpos + xtiky
      slftpsc  = slftpos + xscly
      do 20 k  = 0,nspx
      kskpp    = (k/ntikxs)*ntikxs
      spos     = spmin + k*delx
      spossc   = spos  + xsclx
c
c 4.5.3.1.2 Plot the tick mark
      call line  (spos,  slftpos1, spos,  slftpos2)
c
c 4.5.3.1.3 Plot the scale value
      if(k .eq. kskpp) then
        call setlch(spossc,slftpsc,icshsc,iszhsc,iorhsc,-1)
        write(string,1000) spos
        call wrtstr(string,1)
      endif
  20  continue
c
c 4.5.3.2 Plot the horizontal axis label
      sposlb   = spmin    + mspps*delx0
      slftpslb = slftmin  + slby1
      call setlch(sposlb,slftpslb,icshlb,iszhlb,iorhlb,-1)
      write(string,1010) slabel
      call wrtstr(string,1)
c
c 4.5.4 Set up the Left  vertical axis
c
c 4.5.4.1 Plot the  vertical  axis tick marks and scale
c
c 4.5.4.1.1 Find the tick positions
      spos     = spmin
      spos1    = spos
      spos2    = spos + ytikx
      spossc   = spos + ysclx
      do 30 i  = 0,nspy
      iskpp    = (i/ntikys)*ntikys
      slftpos  = slftmin + i*dely
      slftpssc = slftpos + yscly
c
c 4.5.4.1.2 Plot the tick mark
      call line  (spos1, slftpos,  spos2, slftpos )
c
c 4.5.4.1.3 Plot the scale value
      if(i .eq. iskpp) then
        call setlch(spossc,slftpssc,icsvsy,iszvsy,iorvsy,-1)
        write(string,1100) slftpos
        call wrtstr(string,1)
      endif
  30  continue
c
c 4.5.4.2 Plot the zero line if present
      if(slftmin .lt. zero  .and.  slftmax .gt. zero) then
        call linep(spmin,zero,spmax,zero,kskp0)
      endif
c
c 4.5.5 Plot the  vertical  axis label
c
      sposlb   = spmin    + s0lbx1
      s0pslb   = slftmin  + mss00*dely0
      call setlch(sposlb,s0pslb,icslbs0,iszlbs0,iorlbs0,-1)
      write(string,1110) s0labl
      call wrtstr(string,1)
c
c 4.5.6 Plot the  profile for the left axis: sint0(psi)
c
      call trace(psivmsh,sint0,jpsi1,-1,-1,0.0,0.0)
c
c
c 4.6 Plot the right axis functions 
c
c 4.6.1 Renormalize the map for the right vertical axis function
c
      call map(spmin,spmax,srgtmin,srgtmax,x0mn,x0mx,y0mn,y0mx)
c
c 4.6.2 Set up the Right Vertical axis
c
c 4.6.2.1 Plot the  vertical  axis tick marks and scale
c
c 4.6.2.1.1 Find the tick positions
      spos     = spmax
      spos1    = spos + ztikx
      spos2    = spos
      spossc   = spos + zsclx
      do 50 i  = 0,nspz
      iskpp    = (i/ntikzs)*ntikzs
      srgtpos  = srgtmin + i*delz
      srgtpssc = srgtpos + zsclz
c
c 4.6.2.1.2 Plot the tick mark
      call line  (spos1, srgtpos,  spos2, srgtpos)
c
c 4.6.2.1.3 Plot the scale value
      if(i .eq. iskpp) then
        call setlch(spossc,srgtpssc,icsvsz,iszvsz,iorvsz,-1)
        write(string,1200) srgtpos
        call wrtstr(string,1)
      endif
  50  continue
c
c 4.6.2.2 Plot the zero line if present
      if(srgtmin .lt. zero  .and.  srgtmax .gt. zero) then
        call linep(spmin,zero,spmax,zero,kskp1)
      endif
c
c 4.6.3 Plot the  vertical  axis labels
c
c 4.6.3.1 First surface integral label
      sposlp   = spmax    + s1lbx1
      s1oslp   = srgtmin  + mss10*delz0
      call setlch(sposlp,s1oslp,icslbs1,iszlbs1,iorlbs1,-1)
      write(string,1210) s1labl
      call wrtstr(string,1)
c
c 4.6.3.2 Second surface integral label
      sposlp   = spmax    + s2lbx1
      s2oslp   = srgtmin  + mss20*delz0
      call setlch(sposlp,s2oslp,icslbs2,iszlbs2,iorlbs2,-1)
      write(string,1220) s2labl
      call wrtstr(string,1)
c
c 4.6.3.3 Third surface integral label
      sposlp   = spmax    + s3lbx1
      s3oslp   = srgtmin  + mss30*delz0
      call setlch(sposlp,s3oslp,icslbs3,iszlbs3,iorlbs3,-1)
      write(string,1230) s3labl
      call wrtstr(string,1)
c
c 4.6.3.4 Fourth surface integral label
      sposlp   = spmax    + s4lbx1
      s4oslp   = srgtmin  + mss40*delz0
      call setlch(sposlp,s4oslp,icslbs4,iszlbs4,iorlbs4,-1)
      write(string,1240) s4labl
      call wrtstr(string,1)
c
c
c 4.7 Plot the right vertical axis functions
c
c 4.7.1 Plot the first surface integral
c
c 4.7.1.1 Find the maximum value for normalizing the profile
      s1max    = 0.0
      do 100 j = 1,jpsi1
      if(abs(sint1(j)/sint0(j)) .gt. abs(s1max)) 
     &                s1max  = abs(sint1(j)/sint0(j))
  100 continue
c
      if(s1max .le. 0.0) s1max  = 1.0
c
c 4.7.1.2 Reset the printing parameters for pointc
      call setpch(icschs1,iszchs1,iorchs1,kspchs1)
c
c 4.7.1.3 Plot the normalized values
      do 110 j = 1,jpsi1
      sspt     =     psivmsh(j)
      s1pt     = abs((sint1(j)/sint0(j)) / s1max)
      call pointc(chs1,sspt,s1pt,numbpc,incpcx,incpcy,delpcx,delpcy)
  110 continue
c
c 4.7.2 Plot the second surface integral
c
c 4.7.2.1 Find the maximum value for normalizing the profile
      s2max    = 0.0
      do 150 j = 1,jpsi1
      if(abs(sint2(j)/sint0(j)) .gt. abs(s2max))
     &                s2max  = abs(sint2(j)/sint0(j))
  150 continue
c
      if(s2max .le. 0.0) s2max  = 1.0
c
c 4.7.2.2 Reset the printing parameters for pointc
      call setpch(icschs2,iszchs2,iorchs2,kspchs2)
c
c 4.7.2.3 Plot the normalized values
      do 160 j = 1,jpsi1
      sspt     =     psivmsh(j)
      s2pt     = abs((sint2(j)/sint0(j)) / s2max)
      call pointc(chs2,sspt,s2pt,numbpc,incpcx,incpcy,delpcx,delpcy)
  160 continue
c
c 4.7.3 Plot the third surface integral
c
c 4.7.3.1 Find the maximum value for normalizing the profile
      s3max    = 0.0
      do 200 j = 1,jpsi1
      if(abs(sint3(j)/sint0(j)) .gt. abs(s3max)) 
     &                s3max  = abs(sint3(j)/sint0(j))
  200 continue
c
      if(s3max .le. 0.0) s3max  = 1.0
c
c 4.7.3.2 Reset the printing parameters for pointc
      call setpch(icschs3,iszchs3,iorchs3,kspchs3)
c
c 4.7.3.3 Plot the normalized values
      do 210 j = 1,jpsi1
      sspt     =     psivmsh(j)
      s3pt     = abs((sint3(j)/sint0(j)) / s3max)
      call pointc(chs3,sspt,s3pt,numbpc,incpcx,incpcy,delpcx,delpcy)
  210 continue
c
c 4.7.4 Plot the fourth surface integral
c
c 4.7.4.1 Find the maximum value for normalizing the profile
      s4max    = 0.0
      do 250 j = 1,jpsi1
      if(abs(sint4(j)/sint0(j)) .gt. abs(s4max))
     &                s4max  = abs(sint4(j)/sint0(j))
  250 continue
c
      if(s4max .le. 0.0) s4max  = 1.0
c
c 4.7.4.2 Reset the printing parameters for pointc
      call setpch(icschs4,iszchs4,iorchs4,kspchs4)
c
c 4.7.4.3 Plot the normalized values
      do 260 j = 1,jpsi1
      sspt     =     psivmsh(j)
      s4pt     = abs((sint4(j)/sint0(j)) / s4max)
      call pointc(chs4,sspt,s4pt,numbpc,incpcx,incpcy,delpcx,delpcy)
  260 continue
c
c
c 4.8 End the plot page
c
      call frame(0)
c
c
c
c 5.0 Plot the surface averaged local shear quantities
c
c 5.1 Set the left vertical axis parameters
c
c 5.1.1 Left axis label position fraction
c
      mssv0    = (5*nspy0)/6
      mstv0    = (1*nspy0)/3
c
c 5.1.2 Set the normalized labels and positions
c
      svlbx10  = -2.85
      svlabl   = 'save'
c
      tvlbx10  = -2.85
      tvlabl   = 'tave'
c
c 5.1.3 Left Vertical axis label: case, size, and orientation
c
      icslbsv  = +2
      iszlbsv  = +3
      iorlbsv  =  0
c
      icslbtv  = +2
      iszlbtv  = +3
      iorlbtv  =  0
c
c 5.1.4 Left vertical axis: specific line type parameters
c
      icschsv  = +2
      iszchsv  = +1
      iorchsv  =  0
      kspchsv  =  2
      chsv     = 'S'
c
      icschtv  = +2
      iszchtv  = +1
      iorchtv  =  0
      kspchtv  =  2
      chtv     = 'T'
c
c
c 5.2 Set the right vertical axis parameters
c
c 5.2.1 Right axis label position fraction
c
      msdq0    = (5*nspz0)/6
      msdl0    = (1*nspz0)/3
c
c 5.2.2 Set the normalized labels and positions
c
      dqlbx10  = +1.050
      dllbx10  = +1.050
      dqlabl   = 'qprime'
      dllabl   = 'dlnl'
c
c 5.2.3 Right vertical axis label: case, size, and orientation
c
      icslbdq  = +2
      iszlbdq  = +3
      iorlbdq  =  0
c
      icslbdl  = +2
      iszlbdl  = +3
      iorlbdl  =  0
c
c 5.2.4 Right vertical axis: specific line type parameters
c
      icschdq  = +2
      iszchdq  = +1
      iorchdq  =  0
      kspchdq  =  2
      chdq     = 'q'
c
      icschdl  = +2
      iszchdl  = +1
      iorchdl  =  0
      kspchdl  =  2
      chdl     = 'L'
c
c
c 5.3 Set the physical plot parameters
c
c 5.3.1 Find the function minima and maxima
c
c 5.3.1.1 Search the physically scaled data
c         Note that psivmsh and the integrals are stored in forward order
      svmin    =  0.00
      svmax    = -bigno
      tvmin    =  0.00
      tvmax    = -bigno
      do 300 j = 1,jpsi1
      svmin    = amin1(svmin,savge(j))
      svmax    = amax1(svmax,savge(j))
      tvmin    = amin1(tvmin,tavge(j))
      tvmax    = amax1(tvmax,tavge(j))
  300 continue
c
      dqmin    =  0.00
      dqmax    = -bigno
      do 310 j = 1,jpsi1
      dqmin    = amin1(dqmin,qprime(j))
      dqmax    = amax1(dqmax,qprime(j))
  310 continue
c
      dlmin    =  0.00
      dlmax    = -bigno
      do 320 j = 1,jpsi1
      dlmin    = amin1(dlmin,dlnlval(j))
      dlmax    = amax1(dlmax,dlnlval(j))
  320 continue
c
c 5.3.1.2 Adjust the maxima as needed
      isvmx    = ifix (svmax) + 1
      svmax    = float(isvmx)
c
      itvmx    = ifix (tvmax) + 1
      stmax    = float(itvmx)
c
      idqmx    = ifix (dqmax) + 1
      dqmax    = float(idqmx)
c
      idlmx    = ifix (dlmax) + 1
      dlmax    = float(idlmx)
c
c 5.3.1.3 Set the left and right axis extrema
      slftmin  = amin1(svmin,tvmin)
      slftmax  = amax1(svmax,tvmax)
      if(slftmin .gt. zero) slftmin = zero
      if(slftmax .lt. zero) slftmax = zero
c
      srgtmin  = amin1(dqmin,dlmin)
      srgtmax  = amax1(dqmax,dlmax)
      if(srgtmin  .gt. zero) srgtmin = zero
      if(srgtmax  .lt. zero) srgtmax = zero
c
c 5.3.2 Set the number of tick marks
c
c 5.3.2.1 Left  Vertical axis
      nspy1    = ifix(slftmax)
      nspy1    = min0(nspy1,nspy0)
      nspy1    = max0(nspy1,1)
      if    (nspy1 .le. 4) then
        if    (nspy1 .le. 1) then
          nspy   = 5*nspy1
        elseif(nspy1 .le. 2) then
          nspy   = 4*nspy1
        else
          nspy   = 2*nspy1
        endif
      elseif(nspy1 .gt. 4) then
        nspy   =   nspy1
      endif
c
c 5.3.2.2 Right Vertical axis
      nspz1    = ifix(srgtmax)
      nspz1    = min0(nspz1,nspz0)
      nspz1    = max0(nspz1,1)
      if    (nspz1 .le. 4) then
        if    (nspz1 .le. 1) then
          nspz   = 5*nspz1
        elseif(nspz1 .le. 2) then
          nspz   = 4*nspz1
        else
          nspz   = 2*nspz1
        endif
      elseif(nspz1 .gt. 4) then
        nspz   =   nspz1
      endif
c
c 5.3.3 Normalize the physical spacing increments to the physical range
c
c 5.3.3.1 Standard spacing
      delx0    = (spmax-spmin)/nspx0
      dely0    = (slftmax-slftmin)/nspy0
      delz0    = (srgtmax-srgtmin)/nspz0
c
c 5.3.3.2 Tick spacing
      delx     = (spmax-spmin)/nspx
      dely     = (slftmax-slftmin)/nspy
      delz     = (srgtmax-srgtmin)/nspz
c
c 5.3.3.3 Tick scale values for vertical axes
      ntikys   = min0(ntikysm,nspy)
      ntikys   = max0(ntikys,  1  )
      if(nspy .lt. 5) ntikys  = 1
c
      ntikzs   = min0(ntikzsm,nspz)
      ntikzs   = max0(ntikzs,  1  )
      if(nspz .lt. 5) ntikzs  = 1
c
c 5.3.3.4 Tick lengths
      xtiky    = xtiky0*dely0
      ytikx    = ytikx0*delx0
      ztikx    = ztikx0*delx0
c
c 5.3.4 Axis scale positions
c
c 5.3.4.1 Label shift along relevant axis
      xsclx    = xsclx0*delx0
      yscly    = yscly0*dely0
      zsclz    = zsclz0*delz0
c
c 5.3.4.2 Label shift away from relevant axis
      xscly    = xscly0*dely0
      ysclx    = ysclx0*delx0
      zsclx    = zsclx0*delx0
c
c 5.3.5 Axis label positions
c
c 5.3.5.1 Horizontal axis label
      slby1    = slby10*dely0
c
c 5.3.5.2 Left vertical axis labels
      svlbx1   = svlbx10*delx0
      tvlbx1   = tvlbx10*delx0
c
c 5.3.5.3 Right vertical axis labels
      dqlbx1   = dqlbx10*delx0
      dllbx1   = dllbx10*delx0
c
c 5.3.5.4 Scale format selection
c
c 5.3.5.4.1 Format padding
      iformlf0 = 3
      iformrt0 = 3
c
c 5.3.5.4.2 Left vertical axis
      if(abs(slftmax) .eq. 0.0) aloglfmx = 0.0 
      if(abs(slftmax) .ne. 0.0) aloglfmx = alog10(abs(slftmax))
      if(abs(slftmin) .eq. 0.0) aloglfmn = 0.0
      if(abs(slftmin) .ne. 0.0) aloglfmn = alog10(abs(slftmin))
      loglfmn  = ifix(aloglfmn)
      loglfmx  = ifix(aloglfmx)
      lformlf  = max0(loglfmn,loglfmx)
      if(lformlf .le. iformlf0) iformlf = 0
      if(lformlf .gt. iformlf0) iformlf = lformlf - iformlf0
c
c 5.3.5.4.3 Right vertical axis
      if(abs(srgtmax) .eq. 0.0) alogrtmx = 0.0 
      if(abs(srgtmax) .ne. 0.0) alogrtmx = alog10(abs(srgtmax))
      if(abs(srgtmin) .eq. 0.0) alogrtmn = 0.0
      if(abs(srgtmin) .ne. 0.0) alogrtmn = alog10(abs(srgtmin))
      logrtmn  = ifix(alogrtmn)
      logrtmx  = ifix(alogrtmx)
      lformrt  = max0(logrtmn,logrtmx)
      if(lformrt .le. iformrt0) iformrt = 0
      if(lformrt .gt. iformrt0) iformrt = lformrt - iformrt0
c
c
c 5.4 Print out the plot heading
c
c 5.4.1 Set up the plot in normalized coordinates
c
      call map(smn,smx,smn,smx,smn0,smx0,smn0,smx0)
c
c 5.4.2 Set the position and plot the text
c
      call setlch(xxhd,  yyhd,  icashd,isizhd,iornhd,-1)
      call gtext('Plot of averaged local shear',28,-1)
c
c
c 5.5 Plot the left vertical axis functions
c
c 5.5.1 Renormalize the map in physical units
c
      call map(spmin,spmax,slftmin,slftmax,x0mn,x0mx,y0mn,y0mx)
c
c 5.5.2 Plot the four borders in physical units
c
      call line(spmin,slftmin,spmin,slftmax)
      call line(spmin,slftmax,spmax,slftmax)
      call line(spmax,slftmax,spmax,slftmin)
      call line(spmin,slftmin,spmax,slftmin)
c
c 5.5.3 Set up the Horizontal axis
c
c 5.5.3.1 Plot the horizontal axis tick marks and scale
c
c 5.5.3.1.1 Find the tick positions
      slftpos  = slftmin
      slftpos1 = slftpos
      slftpos2 = slftpos + xtiky
      slftpsc  = slftpos + xscly
      do 330 k = 0,nspx
      kskpp    = (k/ntikxs)*ntikxs
      spos     = spmin + k*delx
      spossc   = spos  + xsclx
c
c 5.5.3.1.2 Plot the tick mark
      call line  (spos,  slftpos1, spos,  slftpos2)
c
c 5.5.3.1.3 Plot the scale value
      if(k .eq. kskpp) then
        call setlch(spossc,slftpsc,icshsc,iszhsc,iorhsc,-1)
        write(string,2000) spos
        call wrtstr(string,1)
      endif
  330 continue
c
c 5.5.3.2 Plot the horizontal axis label
      sposlb   = spmin    + mspps*delx0
      slftpslb = slftmin  + slby1
      call setlch(sposlb,slftpslb,icshlb,iszhlb,iorhlb,-1)
      write(string,2010) slabel
      call wrtstr(string,1)
c
c 5.5.4 Set up the Left  vertical axis
c
c 5.5.4.1 Plot the  vertical  axis tick marks and scale
c
c 5.5.4.1.1 Find the tick positions
      spos     = spmin
      spos1    = spos
      spos2    = spos + ytikx
      spossc   = spos + ysclx
      do 340 i = 0,nspy
      iskpp    = (i/ntikys)*ntikys
      slftpos  = slftmin + i*dely
      slftpssc = slftpos + yscly
c
c 5.5.4.1.2 Plot the tick mark
      call line  (spos1, slftpos,  spos2, slftpos )
c
c 5.5.4.1.3 Plot the scale value
      if(i .eq. iskpp) then
        call setlch(spossc,slftpssc,icsvsy,iszvsy,iorvsy,-1)
        if(iformlf .eq. 0) write(string,2100) slftpos
        if(iformlf .eq. 1) write(string,2102) slftpos
        if(iformlf .eq. 2) write(string,2105) slftpos
        if(iformlf .gt. 2) write(string,2107) slftpos
        call wrtstr(string,1)
      endif
 340  continue
c
c 5.5.4.2 Plot the zero line if present
      if(slftmin .lt. zero  .and.  slftmax .gt. zero) then
        call linep(spmin,zero,spmax,zero,kskp2)
      endif
c
c 5.5.5 Plot the  vertical  axis labels
c
c 5.5.5.1 First surface average label
      sposlp   = spmin    + svlbx1
      svoslp   = slftmin  + mssv0*dely0
      call setlch(sposlp,svoslp,icslbsv,iszlbsv,iorlbsv,-1)
      write(string,2110) svlabl
      call wrtstr(string,1)
c
c 5.5.5.2 Second surface average label
      sposlp   = spmin    + tvlbx1
      tvoslp   = slftmin  + mstv0*dely0
      call setlch(sposlp,tvoslp,icslbtv,iszlbtv,iorlbtv,-1)
      write(string,2120) tvlabl
      call wrtstr(string,1)
c
c 5.5.6 Plot the  profile for the left axis: savge(psi) and tavge(psi)
c
      call trace(psivmsh,savge,jpsi1,-1,-1,0.0,0.0)
      call trace(psivmsh,tavge,jpsi1,-1,-1,0.0,0.0)
c
c 5.5.7 Plot the profile characters for savge
c
c 5.5.7.1 Reset the printing parameters for pointc
      call setpch(icschsv,iszchsv,iorchsv,kspchsv)
c
c 5.5.7.2 Plot the characters
      do 410 j = 1,jpsi1
      sspt     =     psivmsh(j)
      svpt     =     savge (j)
      call pointc(chsv,sspt,svpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  410 continue
c
c 5.5.8 Plot the profile characters for tavge
c
c 5.5.8.2 Reset the printing parameters for pointc
      call setpch(icschtv,iszchtv,iorchtv,kspchtv)
c
c 5.5.8.3 Plot the characters
      do 420 j = 1,jpsi1
      sspt     =     psivmsh(j)
      tvpt     =     tavge (j)
      call pointc(chtv,sspt,tvpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  420 continue
c
c
c 5.6 Plot the right axis functions
c
c 5.6.1 Renormalize the map for the right vertical axis function
c
      call map(spmin,spmax,srgtmin,srgtmax,x0mn,x0mx,y0mn,y0mx)
c
c 5.6.2 Set up the Right Vertical axis
c
c 5.6.2.1 Plot the  vertical  axis tick marks and scale
c
c 5.6.2.1.1 Find the tick marks
      spos     = spmax
      spos1    = spos + ztikx
      spos2    = spos
      spossc   = spos + zsclx
      do 350 i = 0,nspz
      iskpp    = (i/ntikzs)*ntikzs
      srgtpos  = srgtmin + i*delz
      srgtpssc = srgtpos + zsclz
c
c 5.6.2.1.2 Plot the tick mark
      call line  (spos1, srgtpos,  spos2, srgtpos)
c
c 5.6.2.1.3 Plot the scale value
      if(i .eq. iskpp) then
        call setlch(spossc,srgtpssc,icsvsz,iszvsz,iorvsz,-1)
        if(iformrt .eq. 0) write(string,2200) srgtpos
        if(iformrt .eq. 1) write(string,2202) srgtpos
        if(iformrt .eq. 2) write(string,2205) srgtpos
        if(iformrt .gt. 2) write(string,2207) srgtpos
        call wrtstr(string,1)
      endif
 350  continue
c
c 5.6.2.2 Plot the zero line if present
      if(srgtmin .lt. zero  .and.  srgtmax .gt. zero) then
        call linep(spmin,zero,spmax,zero,kskp3)
      endif
c
c 5.6.3 Plot the  vertical  axis labels
c
      sposlp   = spmax    + dqlbx1
      dgoslp   = srgtmin  + msdq0*delz0
      call setlch(sposlp,dgoslp,icslbdq,iszlbdq,iorlbdq,-1)
      write(string,2210) dqlabl
      call wrtstr(string,1)
c
      sposlp   = spmax    + dllbx1
      dgoslp   = srgtmin  + msdl0*delz0
      call setlch(sposlp,dgoslp,icslbdl,iszlbdl,iorlbdl,-1)
      write(string,2220) dllabl
      call wrtstr(string,1)
c
c
c 5.7 Plot the right vertical axis functions: qprime(psi) and dlnlval(psi)
c
c 5.7.1 Plot the curves for the right vertical axis
c
      call trace(psivmsh,qprime, jpsi1,-1,-1,0.0,0.0)
      call trace(psivmsh,dlnlval,jpsi1,-1,-1,0.0,0.0)
c
c 5.7.2 Plot the profile characters for qprime
c
c 5.7.2.1 Reset the printing parameters for pointc
      call setpch(icschdq,iszchdq,iorchdq,kspchdq)
c
c 5.7.2.2 Plot the characters
      do 450 j = 1,jpsi1
      sspt     =     psivmsh(j)
      dqpt     =     qprime(j)
      call pointc(chdq,sspt,dqpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  450 continue
c
c 5.7.3 Plot the profile characters for dlnlval
c
c 5.7.3.1 Reset the printing parameters for pointc
      call setpch(icschdl,iszchdl,iorchdl,kspchdl)
c
c 5.7.3.2 Plot the characters
      do 460 j = 1,jpsi1
      sspt     =     psivmsh(j)
      dlpt     =     dlnlval(j)
      call pointc(chdl,sspt,dlpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  460 continue
c
c
c 5.8  End the plot page
c
      call frame(0)
c
c
c
c 6.0 Return, and end
c
      return
c
 1000 format(f8.4)
 1010 format(a8)
 1100 format(f7.2)
 1110 format(a3)
 1200 format(f4.2)
 1210 format(a6)
 1220 format(a6)
 1230 format(a6)
 1240 format(a6)
 2000 format(f8.4)
 2010 format(a8)
 2100 format(f8.2)
 2102 format(f9.2)
 2105 format(f10.2)
 2107 format(f16.2)
 2110 format(a6)
 2120 format(a6)
 2200 format(f8.2)
 2202 format(f9.2)
 2205 format(f10.2)
 2207 format(f16.2)
 2210 format(a6)
 2220 format(a6)
      end
      subroutine plotqpp
c
c -------------------------------------------------------------
c
c plot q fqpi and qprime profile values
c
c -------------------------------------------------------------
c
      parameter (npx=200,ncx=2*npx,np1=npx+1,np2=npx+2,nc1=ncx+1)
      parameter (npp=929,ntt=2*npp-1,npb=npp)
      parameter (npq=np2)
      parameter (nlx=1441)
      parameter (nvn=7)
c
      character*132 string
c
      character*1   chqv,   chfq,   chqp,   chql
      character*8   qvlabel,fqlabel,qplabel,qllabel
      character*8   qvalue
      character*8   splabel
c
      character*1   lbchpv, lbchsp, lbchrh, lbchro,
     &              lbchtf, lbchst ,lbchvl, lbchsv
      character*8   pclab,  tclab,  labpsi, 
     &              labpsiv,labpsin,labsqtp,
     &              labrhov,labrho2,labrhon,
     &              labtorf,labtorn,labsqtt,
     &              labvolm,labvoln,labsqtv,
     &              labdrh2,labdrho,labdtor,
     &              labdsqt,labdvol,labdsqv
      character*16  filout,filinp,fileql,filmap,filgrd,filrdm,filrdn,
     &              filfst,filplt
c
      common/cnst/  pi,hlfpi,twopi,amu,boltzm,echarg,pmass,roundff,bigno
      common/inpd/  ntor,ncase,idnsty,ndnxp0,ndnxp1,ndnxp2,ifastp,
     &              nmap,neqtyp,mpreset,ndoublt,ndivert,nmtype,
     &              jpsi,itht,jpsi1,jpsi2,ithtp,jsep,
     &              isym,igrid,nham1,nham2,nham3,npowr,
     &              bfieldf,rdefolt,qxin,btdes,qsurf,
     &              gamma,gamav,rmantl,betaf,zeffect,
     &              ncorr,corrfac,nccellr,peqpk0,peqpk1,peqpk2
      common/mapd/  mapmaxd,dpsisl,dpsisd,nqaxis,nwtmag,nfitmax,nfitpts,
     &              ifitrng,jfitrng,jfitchk,fitchek,cnvmag,epsaxs,
     &              maxerlp,maxerlv,delbox,delboz,delac,delav,
     &              delstsf,delstlp,delstlv,prfrac,
     &              nerstop,nerprnt,qptol,tolspln,tolbchi,tolbtor,
     &              tolsymm,tolaugm,errsep,precisn,plossmx,
     &              narcmx,ntrymx,ntdecr,ntmmin,npfit,npcmin,
     &              kuttaop,nrkmax0,nrkmax1,numstp,nwtfitp,nwtfitm,
     &              stepfac,flxstp,psispl,tolnwtp,tolnwtm,
     &              delpakf,delpakc,delpkf,delpkc,psichek,boxtnd,
     &              maptrace,norient,maxcutc,
     &              dresolv,dlclose,pntshft,endtol,
     &              narcln,nangax,nanglm,nbpmax,nwtmax,nslmax,nhfmax,
     &              bperor,sersnm,sertnm,arcmin,delgap,stepcut,
     &              novap,ioutm,iplotm
      common/sarc/  ntmax,ntmsh,npcf,npcb,npc,xp(nlx),zp(nlx),gsq(nlx),
     &              arc(nlx),tp(nlx),arcc(nlx),tpp(nlx),bcd(4),
     &              csx(3,nlx),csz(3,nlx),cseq1(3,nlx),cseq2(3,nlx),
     &              st1(nlx),st2(nlx),st3(nlx),csveq(3,nlx),
     &              sv0(nlx),sv1(nlx),sv2(nlx),sv3(nlx),sv4(nlx),
     &              sv5(nlx)
      common/ort1/  ndim1,ndim2,xma,zma,rcnt,btor,totcur,
     &              psimax0,psilim0,delpsir0,psimax1,psilim1,delpsir1,
     &              psimax,psilim,delpsir,dpsidec,dpsisv,dpedg1,dpedg2,
     &              eaxe,taxe,qaxe,fqiaxe,fqmaxe,paxe,faxe,ppaxe,ffpaxe,
     &              dnaxe,pfaxe,qppaxe,dlnaxe,dpnaxe,dhnaxe,
     &              arcaxe,dnnorm,pfnorm,
     &              qlim,fqilim,fqmlim,plim,flim,pplim,ffplim,
     &              dnlim,pflim,xlim,zlim,
     &              qmer,btmer,gkappa,arcprev,chicmax,
     &              charcmx,chpstmx,chhammx,arcnrm0,pstnrm0,hamnrm0,
     &              elongax1,elongax2,elongaxs,elongaxp,elongtn(npp)
      common/ort2/  psival(np1),pp(np1),ffp(np1),p(np1),f(np1),
     &              dnsty(np1),pfast(np1),qp(np1),fqpi(np1),fqim(np1),
     &              qdel(np1),psinitl(np1),qpinitl(np1),
     &              qprime(np1),qpint(np1),fqint(np1),
     &              psivl0(np2),rst(np2),zst(np2),rstp(np2),zstp(np2),
     &              rsinitl(np2),zsinitl(np2),
     &              arcsurf(nc1),pestchi(nc1),hamlchi(nc1),
     &              alfarc (nc1),alfpst (nc1),alfham (nc1),
     &              rs(nc1,np1),zs(nc1,np1)
      common/ort3/  rsrf(nc1),zsrf(nc1),chic(nc1),chie(nc1),dydx(nc1),
     &              xjsc(nc1)
      common/ort4/  rcell  (np1,nc1),zcell  (np1,nc1),
     &              dpsdr  (np1,nc1),dpsdz  (np1,nc1),
     &              chiarcl(np1,nc1),chipest(np1,nc1),chihaml(np1,nc1),
     &              xjcarcl(np1,nc1),xjcpest(np1,nc1),xjchaml(np1,nc1),
     &              alfarcl(np1,nc1),alfpest(np1,nc1),alfhaml(np1,nc1),
     &              arcnorm(np1),    pstnorm(np1),    hamnorm(np1),
     &              dlnlval(np1),    dpnlval(np1),    dhnlval(np1)
      common/ort5/  f3 (np1,nc1),f4 (np1,nc1),f5 (np1,nc1),
     &              f7 (np1,nc1),f8 (np1,nc1),f9 (np1,nc1),
     &              f10(np1,nc1),f11(np1,nc1),f12(np1,nc1),
     &              f13(np1,nc1),f14(np1,nc1),f15(np1,nc1),
     &              f16(np1,nc1),f17(np1,nc1),f18(np1,nc1),
     &              f19(np1,nc1),f20(np1,nc1),f21(np1,nc1),
     &              f22(np1,nc1),f23(np1,nc1),f24(np1,nc1),
     &              f25(np1,nc1),f26(np1,nc1),f27(np1,nc1),
     &              f28(np1,nc1)
      common/ort6/  savgax,tavgax,sntax0,sntax1,sntax2,sntax3,sntax4,
     &              svint(nvn),savge(np1),tavge(np1),sint0(np1),
     &              sint1(np1),sint2(np1),sint3(np1),sint4(np1)
      common/auxmsh/rh2tot,rhotot,flxtot,voltot,
     &              anltflxp,dvrtflxp,anltvolp,dvrtvolp,
     &              analtflx,divrtflx,analtvol,divrtvol,
     &              psivl1 (np2),qp1    (np2),fqpi1  (np2),
     &              qprim1 (np2),sintp0 (np2),psivmsh(np2),
     &              psivalu(np2),psinorm(np2),psisqrt(np2),
     &              psirho2(np2),psirhov(np2),psinrho(np2),
     &              psitorf(np2),psintor(np2),psisqtf(np2),
     &              psivolm(np2),psinvol(np2),psisqvl(np2),
     &              dpsirh2(np2),dpsirho(np2),dpsinrh(np2),
     &              dpsitor(np2),dpsintf(np2),dpsisqt(np2),
     &              dpsivol(np2),dpsinvl(np2),dpsisqv(np2)
       common/ratnl/jqmin, jqmax, numbqs,nq1totl,
     &              rcentr,zcentr,rminor,zminor,
     &              asprat00,asprat01,asprat10,asprat11,
     &              pminvl,qminvl,qpminv,pmaxvl,qmaxvl,qpmaxv,
     &              psivlmnq,    psivnmnq,    psisqmnq,
     &              psir2mnq,    psirhmnq,    psinrmnq,
     &              psitfmnq,    psintmnq,    psistmnq,
     &              psivmmnq,    psinvmnq,    psisvmnq,
     &              dpsr2mnq,    dpsrhmnq,    dpsnrmnq,
     &              dpstfmnq,    dpsntmnq,    dpsstmnq,
     &              dpsvmmnq,    dpsnvmnq,    dpssvmnq,
     &              psivlmxq,    psivnmxq,    psisqmxq,
     &              psir2mxq,    psirhmxq,    psinrmxq,
     &              psitfmxq,    psintmxq,    psistmxq,
     &              psivmmxq,    psinvmxq,    psisvmxq,
     &              dpsr2mxq,    dpsrhmxq,    dpsnrmxq,
     &              dpstfmxq,    dpsntmxq,    dpsstmxq,
     &              dpsvmmxq,    dpsnvmxq,    dpssvmxq,
     &              lpsiq  (npq),psivlq (npq),
     &              qprimq (npq),qvalue (npq),psimshq(npq),
     &              psivalq(npq),psinrmq(npq),psisqrq(npq),
     &              psirh2q(npq),psirhoq(npq),psinrhq(npq),
     &              psitorq(npq),psintfq(npq),psisqtq(npq),
     &              psivolq(npq),psinvlq(npq),psisqvq(npq),
     &              dpsir2q(npq),dpsirhq(npq),dpsinrq(npq),
     &              dpsitfq(npq),dpsintq(npq),dpsistq(npq),
     &              dpsivlq(npq),dpsinvq(npq),dpsisvq(npq),
     &              shearps(npq),shearrh(npq),
     &              sheartf(npq),shearvl(npq),
     &              epslrh1 (npq),shearrh1(npq),shearfrh(npq),
     &              epslvl1 (npq),shearvl1(npq),shearfvl(npq)
      common/labels/lbchpv, lbchsp, lbchrh, lbchro,
     &              lbchtf, lbchst, lbchvl, lbchsv,
     &              pclab,  tclab,  labpsi,
     &              labpsiv,labpsin,labsqtp,
     &              labrhov,labrho2,labrhon,
     &              labtorf,labtorn,labsqtt,
     &              labvolm,labvoln,labsqtv,
     &              labdrh2,labdrho,labdtor,
     &              labdsqt,labdvol,labdsqv
      common/pldf/  x0min,x0max,y0min,y0max 
      common/inout/ kutty,kuout,kuinp,kueql,kumap,kugrd,kurdm,kurdn,
     &              kufst,krdeqlb,krdmgta,krdngta,krdfast
      common/flnm/  filout,filinp,fileql,filmap,filgrd,filrdm,filrdn,
     &              filfst,filplt
c
      dimension profle(nlx)
      equivalence (st1 (1),profle(1))
c
c -------------------------------------------------------------
c
c
c
c 1.0 Initialization
c
c 1.1 Return if no plot is required
c
      if(iplotm .lt. 6) return
c
c
c 1.2 Constants
c
      bignum    = bigno
      smallno   = abs(roundff)
c
c
c 1.3 Physical plot dimensions
c
c 1.3.1 Plot borders
c
      x0mn      = 0.120
      x0mx      = 0.820
      y0mn      = 0.200
      y0mx      = 0.900
c
c 1.3.2 Heading position
c
      xxhd      = 0.0875
      yyhd      = 0.9500
c
c
c 1.4 Heading: case, size, and orientation
c
      icashd    = +2
      isizhd    = +3
      iornhd    =  0
c
c
c 1.5 Label position parameters
c
c 1.5.1 Position increments
c
      nspx0     = 16
      nspy0     = 16
      nspz0     = 30
c
c 1.5.2 Axis label position fraction
c
      mspps0    = (2*nspx0)/3
      mspqv0    = (3*nspy0)/4
      mspfq0    = (2*nspy0)/3
      mspqp0    = (3*nspz0)/4
      mspql0    = (1*nspz0)/2
c
c
c 1.6 Default parameters for pointc
c
      numbpc    = 1
      incpcx    = 0
      incpcy    = 0
      delpcx    = 0.0
      delpcy    = 0.0
c
c
c
c 2.0 Axis settings
c
c 2.1 Set the normalized tick parameters
c
c 2.1.1 Tick scale skipping
c
      ntikxs0   = 2
      ntikys0   = 2
      ntikzs0   = 4
c
c 2.1.2 Normalized tick lengths
c
      xtiky0    = +0.250
      ytikx0    = +0.250
      ztikx0    = -0.250
c
c 2.1.3 Axis extensions
c
      xtndl     = 1.000
      xtndr     = 1.000
c
c 2.1.4 Tick scale positions
c
c 2.1.4.1 Label shift along relevant axis
      xsclx0    = +0.250
      yscly0    = -0.250
      zsclz0    = -0.250
c
c 2.1.4.2 Label shift away from relevant axis
      xscly0    = -2.350
      ysclx0    = -2.100
      zsclx0    = -0.200
c
c
c 2.2 Set the axis scales: case, size, and orientation
c
c 2.2.1 Horizontal axis
c
      icshsc    = +1
      iszhsc    = +1
      iorhsc    = +1
c
c 2.2.2 Left Vertical axis
c
      icsvsy    = +1
      iszvsy    = +1
      iorvsy    =  0
c
c 2.2.3 Right Vertical axis
c
      icsvsz    = +1
      iszvsz    = +1
      iorvsz    =  0
c
c
c 2.3 Axis labels specific to quantities plotted
c
c 2.3.1 Horizontal axis
c
c 2.3.1.1 Label and label position
      splby10   = -2.850
      splabel   = pclab
c
c 2.3.1.2 Horizontal axis label: case, size, and orientation
      icshlb    = +1
      iszhlb    = +3
      iorhlb    =  0
c
c 2.3.2 Left  Vertical axis
c
c 2.3.2.1 Label and label positions
c
c 2.3.2.1.1 Label and label position for q(psi)
      qvlbx10   = -2.900
      qvlabel   = ' q '
c
c 2.3.2.1.2 Label and label and scale position for fqpi(psi)
      fqlbx10   = -2.900
      fqscx10   = -0.5000
      fqscy10   = -0.6250
      fqlabel   = 'fqpi'
c
c 2.3.2.2 Left Vertical axis label: case, size, and orientation
c
c 2.3.2.2.1 Axis label: case, size, and orientation for q(psi)
      icslbqv   = +1
      iszlbqv   = +3
      iorlbqv   =  0
c
c 2.3.2.2.2 Axis label: case, size, and orientation for fqpi(psi)
      icslbfq   = +1
      iszlbfq   = +3
      iorlbfq   =  0
c
c 2.3.2.2.3 Axis scale: case, size, and orientation for fqpi(psi)
      icsscfq   = +1
      iszscfq   = +1
      iorscfq   =  0
c
c 2.3.3 Right Vertical axis
c
c 2.3.3.1 Label and label position 
c
c 2.3.3.1.1 Label and label position for qprime(psi)
      qplbx10   = +1.0000
      qplabel   = 'dq/dpsi'
c
c 2.3.3.1.2 Label and label and scale position for modified qprime
c           profile
      qllbx10   = +0.2500
      qlscx10   = -0.1000
      qlscy10   = -1.2500
      qllabel   = 'qp (mod)'
c
c 2.3.3.2.2 Right Vertical axis label: case, size, and orientation
c
c 2.3.3.2.2.1 Axis label: case, size, and orientation for qp(psi)
      icslbqp   = +1
      iszlbqp   = +3
      iorlbqp   =  0
c
c 2.3.3.2.2.2 Axis label: case, size, and orientation for modified qprime
c           profile
      icslbql   = +1
      iszlbql   = +3
      iorlbql   =  0
c
c 2.3.3.2.3 Axis scale: case, size, and orientation for modified qprime
c           profile
      icsscql   = +1
      iszscql   = +1
      iorscql   =  0
c
c
c 2.4 Line type parameters
c
c 2.4.1 q profile
c
      icschqv   = +2
      iszchqv   = +1
      iorchqv   =  0
      kspchqv   = +2
      chqv      = 'q'
c
c 2.4.2 fqpi profile
c
      icschfq   = +2
      iszchfq   = +1
      iorchfq   =  0
      kspchfq   = +2
      chfq      = 'f'
c
c 2.4.3 qprime profile
c
      icschqp   = +2
      iszchqp   = +1
      iorchqp   =  0
      kspchqp   = +2
      chqp      = '+'
c
c 2.4.4 Modified qprime profile
c
      icschql   = +2
      iszchql   = +1
      iorchql   =  0
      kspchql   = +2
      chql      = 'm'
c
c
c
c 3.0 Find the function minima and maxima
c
c 3.1 Search the physically scaled data
c
      spmin      = psivmsh(  1  )
      spmax      = psivmsh(jpsi2)
      qvmin      = +bignum
      qvmax      = -bignum
      fqmin      = +bignum
      fqmax      = -bignum
      qpmin      = +bignum
      qpmax      = -bignum
      qlmin      = +bignum
      qlmax      = -bignum
      sprange    =  spmax-spmin
      do 10  j   = 1,jpsi2
      jp0        = j
      psvlu      = psivmsh(jp0)
      qvvlu      = qp1    (jp0)
      fqvlu      = fqpi1  (jp0)
c
      qpvlu      = qprim1 (jp0)
      ppvl       = (spmax-psvlu)/sprange
      ppvlu      = amax1(ppvl,smallno)
      if(npowr .lt. -1) qlvlu   = alog10(amax1(qpvlu*qpvlu,smallno))
      if(npowr .eq. -1) qlvlu   = qpvlu/(1.0 +  abs(alog10(ppvlu)))
      if(npowr .eq.  0) qlvlu   = qpvlu*ppvl
      if(npowr .gt.  0) qlvlu   = qpvlu/(1.0 + (abs(1.0/ppvlu))**npowr)
c
      qvmin      = amin1(qvmin,qvvlu)
      qvmax      = amax1(qvmax,qvvlu)
      fqmin      = amin1(fqmin,fqvlu)
      fqmax      = amax1(fqmax,fqvlu)
      qpmin      = amin1(qpmin,qpvlu)
      qpmax      = amax1(qpmax,qpvlu)
      qlmin      = amin1(qlmin,qlvlu)
      qlmax      = amax1(qlmax,qlvlu)
   10 continue
c
c
c 3.2 Adjust the maxima as needed
c
c 3.2.1 Horizontal axis extrema
c
      spmin0     = spmin
      spmax0     = spmax
      sprang0    = spmax0 - spmin0
c
c 3.2.2 q(psi)  profile extrema
c
c 3.2.2.1 Push the extrema to the next integer value
      iqvmx      = ifix (qvmax) + 1
      qvmax0     = float(iqvmx)
c
      if(qvmin .eq. 0.0) iqvmn  = ifix (qvmin)
      if(qvmin .ne. 0.0) iqvmn  = ifix (qvmin) - 1
      qvmin0     = float(iqvmn)
c
c 3.2.2.2 Ensure the origin is included if not in the range
      if(qvmax0 .gt. 0.0) qvmin0  = amin1(qvmin0,0.0)
      if(qvmin0 .lt. 0.0) qvmax0  = amax1(qvmax0,0.0)
c
c 3.2.3 fqpi(psi)  profile extrema
c
c 3.2.3.1 Push the extrema to the next integer value
      ifqmx      = ifix (fqmax) + 1
      fqmax0     = float(ifqmx)
c
      if(fqmin .eq. 0.0) ifqmn  = ifix (fqmin)
      if(fqmin .ne. 0.0) ifqmn  = ifix (fqmin) - 1
      fqmin0     = float(ifqmn)
c
c 3.2.3.2 Ensure the origin is included if not in the range
      if(fqmax0 .gt. 0.0) fqmin0  = amin1(fqmin0,0.0)
      if(fqmin0 .lt. 0.0) fqmax0  = amax1(fqmax0,0.0)
c
c 3.2.4 qprime(psi)  profile extrema
c
c 3.2.4.1 Push the extrema to the next nice number value
      qplog1     = alog10(amax1(abs(qpmax),smallno))
      iqplog1    = ifix (qplog1)
      if    (qplog1 .gt. 0.0) then
         incqp1   = 10**iqplog1
      elseif(qplog1 .lt. 0.0) then
         incqp1  = 10**iqplog1
      elseif(qplog1 .eq. 0) then
         incqp1  = 0
      endif
      iqpmx      = ifix (qpmax) + incqp1
      qpmax0     = float(iqpmx)
c
      qplog0     = alog10(amax1(abs(qpmin),smallno))
      iqplog0    = ifix (qplog0)
      if    (qplog0 .gt. 0.0) then
         incqp0   = 10**iqplog0
      elseif(qplog0 .lt. 0.0) then
         incqp0  = 10**iqplog0
      elseif(qplog0 .eq. 0) then
         incqp0  = 0
      endif
      iqpmn      = ifix (qpmin) - incqp0
      qpmin0     = float(iqpmn)
c
c 3.2.4.2 Ensure the origin is included if not in the range
      if(qpmax0 .gt. 0.0) qpmin0  = amin1(qpmin0,0.0)
      if(qpmin0 .lt. 0.0) qpmax0  = amax1(qpmax0,0.0)
c
c 3.2.5 Log modified qprime(psi)  profile extrema
c
c 3.2.5.1 Push the extrema to the next integer value
      iqlmx      = ifix (qlmax) + 1
      qlmax0     = float(iqlmx)
c
      if(qlmin .eq. 0.0) iqlmn  = ifix (qlmin)
      if(qlmin .ne. 0.0) iqlmn  = ifix (qlmin) - 1
      qlmin0     = float(iqlmn)
c
c 3.2.5.2 Ensure the origin is included if not in the range
      if(qlmax0 .gt. 0.0) qlmin0  = amin1(qlmin0,0.0)
      if(qlmin0 .lt. 0.0) qlmax0  = amax1(qlmax0,0.0)
c
c
c 3.3 Set the range parameters
c
c 3.3.1 Range for qp
c
      if(qvmin0 .ge. 0.0  .and.  qvmax0 .gt. 0.0) lqpval  = +1
      if(qvmin0 .lt. 0.0  .and.  qvmax0 .gt. 0.0) lqpval  =  0
      if(qvmin0 .lt. 0.0  .and.  qvmax0 .le. 0.0) lqpval  = -1
c
c 3.3.2 Range for fqpi
c
      if(fqmin0 .ge. 0.0  .and.  fqmax0 .gt. 0.0) lfqval  = +1
      if(fqmin0 .lt. 0.0  .and.  fqmax0 .gt. 0.0) lfqval  =  0
      if(fqmin0 .lt. 0.0  .and.  fqmax0 .le. 0.0) lfqval  = -1
c
c 3.3.3 Range for qprime
c
      if(qpmin0 .ge. 0.0  .and.  qpmax0 .gt. 0.0) lqprim  = +1
      if(qpmin0 .lt. 0.0  .and.  qpmax0 .gt. 0.0) lqprim  =  0
      if(qpmin0 .lt. 0.0  .and.  qpmax0 .le. 0.0) lqprim  = -1
c
c 3.3.4 Range for modified qprime
c
      if(qlmin0 .ge. 0.0  .and.  qlmax0 .gt. 0.0) lqplog  = +1
      if(qlmin0 .lt. 0.0  .and.  qlmax0 .gt. 0.0) lqplog  =  0
      if(qlmin0 .lt. 0.0  .and.  qlmax0 .le. 0.0) lqplog  = -1
c
c
c 3.4 Set the scale factors between plots
c
c 3.4.1 Scale factor between qp and fqpi
c
      if    (lqpval .gt. 0  .and.  lfqval .gt. 0) then
        lscalfq  = +1
        ratiofq  = qvmax0/fqmax0
c
      elseif(lqpval .gt. 0  .and.  lfqval .eq. 0) then
        lscalfq  = +1
        ratiofq  = qvmax0/fqmax0
c
      elseif(lqpval .gt. 0  .and.  lfqval .lt. 0) then
        lscalfq  = +2
        ratiofq  = abs(qvmax0/fqmin0)
c
      elseif(lqpval .eq. 0  .and.  lfqval .gt. 0) then
        lscalfq  = +1
        ratiofq  = qvmax0/fqmax0
c
      elseif(lqpval .eq. 0  .and.  lfqval .eq. 0) then
        if    (abs(fqmax0) .ge. abs(fqmin0)) then
          lscalfq  = +1
          ratiofq  = qvmax0/fqmax0
        elseif(abs(fqmax0) .lt. abs(fqmin0)) then
          lscalfq  = -1
          ratiofq  = qvmin0/fqmin0
        endif
c
      elseif(lqpval .eq. 0  .and.  lfqval .lt. 0) then
        lscalfq  = -1
        ratiofq  = qvmin0/fqmin0
c
      elseif(lqpval .lt. 0  .and.  lfqval .gt. 0) then
        lscalfq  = -2
        ratiofq  = abs(qvmin0/fqmax0)
c
      elseif(lqpval .lt. 0  .and.  lfqval .eq. 0) then
        lscalfq  = -1
        ratiofq  = qvmin0/fqmin0
c
      elseif(lqpval .lt. 0  .and.  lfqval .lt. 0) then
        lscalfq  = -1
        ratiofq  = qvmin0/fqmin0
      else
        lscalfq  = -1
        ratiofq  = qvmin0/fqmin0
        call abortjob
     &        ('plotqpp ',  1,   'Undetermined scale factor for fq    '
     &        ,'lqpval  ', lqpval,   'lfqval  ', lfqval,   -1)
      endif
c
      if(ratiofq .le. 0.0) then
        call abortjob
     &        ('plotqpp ',  2,   'Nonpositive scale factor for fq     '
     &        ,'lqpval  ', lqpval,   'lfqval  ', lfqval,   -1)
      endif
c
c 3.4.2 Scale factor between qprime and the modified qprime
c
      if    (lqprim .gt. 0  .and.  lqplog .gt. 0) then
        lscalqq  = +1
        ratioqq  = qpmax0/qlmax0
c
      elseif(lqprim .gt. 0  .and.  lqplog .eq. 0) then
        lscalqq  = +1
        ratioqq  = qpmax0/qlmax0
c
      elseif(lqprim .gt. 0  .and.  lqplog .lt. 0) then
        lscalqq  = +2
        ratioqq  = abs(qpmax0/qlmin0)
c
      elseif(lqprim .eq. 0  .and.  lqplog .gt. 0) then
        lscalqq  = +1
        ratioqq  = qpmax0/qlmax0
c
      elseif(lqprim .eq. 0  .and.  lqplog .eq. 0) then
        if    (abs(qlmax0) .ge. abs(qlmin0)) then
          lscalqq  = +1
          ratioqq  = qpmax0/qlmax0
        elseif(abs(qlmax0) .lt. abs(qlmin0)) then
          lscalqq  = -1
          ratioqq  = qpmin0/qlmin0
        endif
c
      elseif(lqprim .eq. 0  .and.  lqplog .lt. 0) then
        lscalqq  = -1
        ratioqq  = qpmin0/qlmin0
c
      elseif(lqprim .lt. 0  .and.  lqplog .gt. 0) then
        lscalqq  = -2
        ratioqq  = abs(qpmin0/qlmax0)
c
      elseif(lqprim .lt. 0  .and.  lqplog .eq. 0) then
        lscalqq  = -1
        ratioqq  = qpmin0/qlmin0
c
      elseif(lqprim .lt. 0  .and.  lqplog .lt. 0) then
        lscalqq  = -1
        ratioqq  = qpmin0/qlmin0
      else
        lscalqq  = -1
        ratioqq  = qpmin0/qlmin0
        call abortjob
     &        ('plotqpp ',  3,   'Undetermined scale for mod qprim    '
     &        ,'lqprim  ', lqprim,   'lqplog  ', lqplog,   -1)
      endif
c
      if(ratioqq .le. 0.0) then
        call abortjob
     &        ('plotqpp ',  4,   'Nonpositive scale for mod qprime    '
     &        ,'lqprim  ', lqprim,   'lqplog  ', lqplog,   -1)
      endif
c
c
c
c 4.0 Set the physical plot parameters
c
c 4.1 Set the map extrema from the extrema
c
c 4.1.1 Left scale extrema for qp1 and fqpi1
c
      alfmin    =       qvmin0
      alfmax    = xtndl*qvmax0
c
      alfmin    = amin1(alfmin,ratiofq*fqmin0)
      alfmax    = amax1(alfmax,ratiofq*fqmax0)
c
c 4.1.2 Right scale extrema for qprim1 and modified qprim1
c
      artmin    =       qpmin0
      artmax    = xtndr*qpmax0
c
      artmin    = amin1(artmin,ratioqq*qlmin0)
      artmax    = amax1(artmax,ratioqq*qlmax0)
c
c
c 4.2 Set the number of tick marks
c
c 4.2.1 Horizontal axis
c
      nspx      = 10
c
c 4.2.2 Left  Vertical axis
c
      nspy1     = ifix(alfmax - alfmin)
      nspy1     = min0(nspy1,nspy0)
      nspy1     = max0(nspy1,1)
c
      if    (nspy1 .le. 4) then
        if    (nspy1 .le. 1) then
          nspy    = 5*nspy1
        elseif(nspy1 .le. 2) then
          nspy    = 4*nspy1
        else
          nspy    = 2*nspy1
        endif
c
      elseif(nspy1 .gt. 4) then
        nspy    =   nspy1
      endif
c
c 4.2.3 Right Vertical axis
c
      nspz1     = ifix(artmax - artmin)
      nspz1     = min0(nspz1,nspz0)
      nspz1     = max0(nspz1,1)
c
      if    (nspz1 .le. 4) then
        if    (nspz1 .le. 1) then
          nspz    = 5*nspz1
        elseif(nspz1 .le. 2) then
          nspz    = 4*nspz1
        else
          nspz    = 2*nspz1
        endif
c
      elseif(nspz1 .gt. 4) then
        nspz    =   nspz1
      endif
c
c
c 4.3 Normalize the physical spacing increments to the physical range
c
c 4.3.1 Standard spacing
c
      delx0     = (spmax0 - spmin0)/nspx0
      dely0     = (alfmax - alfmin)/nspy0
      delz0     = (artmax - artmin)/nspz0
c
c 4.3.2 Tick spacing
c
      delx      = (spmax0 - spmin0)/nspx
      dely      = (alfmax - alfmin)/nspy
      delz      = (artmax - artmin)/nspz
c
c
c 4.4 Set the physical spacings
c
c 4.4.1 Tick scale values
c
c 4.4.1.1 Horizontal axis
      ntikxs    = min0(ntikxs0,nspx)
      ntikxs    = max0(ntikxs,  1  )
      if(nspx .lt. 5) ntikxs  = 1
c
c 4.4.1.2 Left Vertical axis
      ntikys    = min0(ntikys0,nspy)
      ntikys    = max0(ntikys,  1  )
      if(nspy .lt. 5) ntikys  = 1
c
c 4.4.1.3 Right Vertical axis
      ntikzs    = min0(ntikzs0,nspz)
      ntikzs    = max0(ntikzs,  1  )
      if(nspz .lt. 5) ntikzs  = 1
c
c 4.4.2 Axis scale positions
c
c 4.4.2.1 Label shift along relevant axis
      xsclx     = xsclx0*delx0
      yscly     = yscly0*dely0
      zsclz     = zsclz0*delz0
c
c 4.4.2.2 Label shift away from relevant axis
      xscly     = xscly0*dely0
      ysclx     = ysclx0*delx0
      zsclx     = zsclx0*delx0
c
c 4.4.3 Tick lengths
c
      xtiky     = xtiky0*dely0
      ytikx     = ytikx0*delx0
      ztikx     = ztikx0*delx0
c
c 4.4.4 Label positions specific to profiles
c
c 4.4.4.1 Horizontal axis
      splby1    = splby10*dely0
c
c 4.4.4.2 Left Vertical Axis
      qvlbx1    = qvlbx10*delx0
      fqlbx1    = fqlbx10*delx0
      fqscx1    = fqscx10*delx0
      fqscy1    = fqscy10*dely0
c
c 4.4.4.3 Right Vertical Axis
      qplbx1    = qplbx10*delx0
      qllbx1    = qllbx10*delx0
      qlscx1    = qlscx10*delx0
      qlscy1    = qlscy10*delz0
c
c 4.4.5 Format selection
c
c 4.4.5.1 Format padding
      iformlf0  = 2
      iformrt0  = 3
      bigscale  = 1.0e+04
c
c 4.4.5.2 Left Vertical Axis
      if(abs(alfmax) .eq. 0.0) aloglfmx = 0.0 
      if(abs(alfmax) .ne. 0.0) aloglfmx = alog10(abs(alfmax))
      if(abs(alfmin) .eq. 0.0) aloglfmn = 0.0
      if(abs(alfmin) .ne. 0.0) aloglfmn = alog10(abs(alfmin))
      loglfmn   = ifix(aloglfmn)
      loglfmx   = ifix(aloglfmx)
      lformlf   = max0(loglfmn,loglfmx)
      if(lformlf .le. iformlf0) iformlf = 0
      if(lformlf .gt. iformlf0) iformlf = lformlf - iformlf0
c
c 4.4.5.3 Right Vertical Axis
      if(abs(artmax) .eq. 0.0) alogrtmx = 0.0 
      if(abs(artmax) .ne. 0.0) alogrtmx = alog10(abs(artmax))
      if(abs(artmin) .eq. 0.0) alogrtmn = 0.0
      if(abs(artmin) .ne. 0.0) alogrtmn = alog10(abs(artmin))
      logrtmn   = ifix(alogrtmn)
      logrtmx   = ifix(alogrtmx)
      lformrt   = max0(logrtmn,logrtmx)
      if(lformrt .le. iformrt0) iformrt = 0
      if(lformrt .gt. iformrt0) iformrt = lformrt - iformrt0
c
c
c 4.5 Additional lines
c
      zero      = 0.0
      kspzero   =  2
c
c
c
c 5.0 Print out the plot heading
c
c 5.1 Set up the plot in normalized coordinates
c
      call map(x0min,x0max,y0min,y0max,x0min,x0max,y0min,y0max)
c
c
c 5.2 Plot the four borders
c
      call line(x0mn,y0mn,x0mn,y0mx)
      call line(x0mn,y0mx,x0mx,y0mx)
      call line(x0mx,y0mx,x0mx,y0mn)
      call line(x0mn,y0mn,x0mx,y0mn)
c
c
c 5.3 Set the position and plot the text
c
      call setlch(xxhd,  yyhd,  icashd,isizhd,iornhd,-1)
      write(string,1000) splabel 
      call wrtstr(string,1)
c
c
c
c 6.0 Plot q(psi) and fqpi(psi)
c
c 6.1 Set the plot page in physical units (left axis)
c
      call map(spmin0,spmax0,alfmin,alfmax,x0mn,x0mx,y0mn,y0mx)
c
c
c 6.2 Set up the Horizontal axis
c
c 6.2.1 Plot the horizontal axis tick marks and scale
c
c 6.2.1.1 Set the base tick parameters
      vrtpos     = alfmin
      vrtpos1    = vrtpos
      vrtpos2    = vrtpos + xtiky
      vrtpossc   = vrtpos + xscly
      do 100 k   = 0,nspx
      kskpp      = (k/ntikxs)*ntikxs
      spos       = spmin0 + k*delx
      spossc     = spos   + xsclx
c
c 6.2.1.2 Plot the tick mark
      call line  (spos, vrtpos1, spos, vrtpos2 )
c
c 6.2.1.3 Plot the scale value
      if(k .eq. kskpp) then
        call setlch(spossc,vrtpossc,icshsc,iszhsc,iorhsc,-1)
        write(string,2000) spos
        call wrtstr(string,1)
      endif
  100 continue
c
c 6.2.2 Plot the horizontal axis label
c
      sposlb    = spmin0  + mspps0*delx0
      vrtposlb  = vrtpos  + splby1
      call setlch(sposlb,vrtposlb,icshlb,iszhlb,iorhlb,-1)
      write(string,2010) splabel
      call wrtstr(string,1)
c
c
c 6.3 Set up the Left  vertical axis
c
c 6.3.1 Plot the  vertical  axis tick marks and scale
c
c 6.3.1.1 Set the base tick parameters
      spos       = spmin0
      spos1      = spos
      spos2      = spos + ytikx
      spossc     = spos + ysclx
      do 110 i   = 0,nspy
      iskpp      = (i/ntikys)*ntikys
      alfpos     = alfmin  + i*dely
      alfpossc   = alfpos  + yscly
c
c 6.3.1.2 Plot the tick mark
      call line  (spos1, alfpos,  spos2, alfpos  )
c
c 6.3.1.3 Plot the scale value
      if(i .eq. iskpp) then
        call setlch(spossc,alfpossc,icsvsy,iszvsy,iorvsy,-1)
        if(iformlf .eq. 0) write(string,2100) alfpos
        if(iformlf .eq. 1) write(string,2102) alfpos
        if(iformlf .eq. 2) write(string,2105) alfpos
        if(iformlf .gt. 2) write(string,2107) alfpos
        call wrtstr(string,1)
      endif
  110 continue
c
c
c 6.4 Plot the q(psi)  profile
c
c 6.4.1 Plot the  vertical  axis label
c
      sposlb   = spmin0  + qvlbx1
      qvposlb  = alfmin  + mspqv0*dely0
      call setlch(sposlb,qvposlb,icslbqv,iszlbqv,iorlbqv,-1)
      write(string,2110) qvlabel
      call wrtstr(string,1)
c
c 6.4.2 Reset the printing parameters for pointc
c
      call setpch(icschqv,iszchqv,iorchqv,kspchqv)
c
c 6.4.3 Plot the q(psi)  profile
c
c 6.4.3.1 Plot the points and store in profle
      do 200 j   = 1,jpsi2
      jp         = j
      sppt       =         psivmsh(jp)
      qvpt       =         qp1    (jp)
      profle(jp) =         qvpt
      call pointc(chqv,sppt,qvpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  200 continue
c
c 6.4.3.2 Plot the curve
      call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
c
c
c 6.5 Plot the fqpi profile
c
c 6.5.1 Plot the  vertical  axis label
c
c 6.5.1.1 Write the label itself
      sposlb   = spmin0  + fqlbx1
      fqposlb  = alfmin  + mspfq0*dely0
      call setlch(sposlb,fqposlb,icslbfq,iszlbfq,iorlbfq,-1)
      write(string,2120) fqlabel
      call wrtstr(string,1)
c
c 6.5.1.2 Plot the scale factor
      sposscf  = sposlb    + fqscx1
      fqpossc  = fqposlb   + fqscy1
      call setlch(sposscf,fqpossc,icsscfq,iszscfq,iorscfq,-1)
c
      if(ratiofq .gt. 1.0) write(string,2130) ratiofq
      if(ratiofq .lt. 1.0) write(string,2140) ratiofq
      if(ratiofq .ne. 1.0) call wrtstr(string,1)
c
c 6.5.2 Reset the printing parameters for pointc
c
      call setpch(icschfq,iszchfq,iorchfq,kspchfq)
c
c 6.5.3 Plot the renormalized values
c
c 6.5.3.1 Plot the points and store in profle
      do 250 j   = 1,jpsi2
      jp         = j
      sppt       =         psivmsh(jp)
      fqpt       = ratiofq*fqpi1  (jp)
      profle(jp) =         fqpt
      call pointc(chfq,sppt,fqpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  250 continue
c
c 6.5.3.2 Plot the curve
      call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
c
c
c 6.6 Plot a line if the origin is included in the range
c
      if(lqpval .eq. 0  .or.  lfqval .eq. 0) then
        call linep(spmin0,zero,spmax0,zero,kspzero)
      endif
c
c
c
c 7.0 Plot the derivative of q and the modified derivative
c
c 7.1 Renormalize the map for the derivative
c
      call map(spmin0,spmax0,artmin,artmax,x0mn,x0mx,y0mn,y0mx)
c
c
c 7.2 Set up the Right Vertical axis
c
c 7.2.1 Plot the  vertical  axis tick marks and scale
c
c 7.2.1.1 Set the base tick parameters
      spos       = spmax0
      spos1      = spos + ztikx
      spos2      = spos
      spossc     = spos + zsclx
      do 300 i   = 0,nspz
      iskpp      = (i/ntikzs)*ntikzs
      artpos     = artmin  + i*delz
      artpossc   = artpos  + zsclz
c
c 7.2.1.2 Plot the tick mark
      call line  (spos1, artpos,  spos2, artpos  )
c
c 7.2.1.3 Plot the scale value
      if(i .eq. iskpp) then
        call setlch(spossc,artpossc,icsvsz,iszvsz,iorvsz,-1)
        if(iformrt .eq. 0) write(string,2200) artpos
        if(iformrt .eq. 1) write(string,2202) artpos
        if(iformrt .eq. 2) write(string,2205) artpos
        if(iformrt .gt. 2) write(string,2207) artpos
        call wrtstr(string,1)
      endif
  300 continue
c
c
c 7.3 Plot the profile for qprime
c
c 7.3.1 Plot the  vertical  axis labels for the q derivative
      sposlb    = spmax0  + qplbx1
      qpposlb   = artmin  + mspqp0*delz0
      call setlch(sposlb,qpposlb,icslbqp,iszlbqp,iorlbqp,-1)
      write(string,2210) qplabel
      call wrtstr(string,1)
c
c 7.3.2 Reset the printing parameters for pointc
c
      call setpch(icschqp,iszchqp,iorchqp,kspchqp)
c
c 7.3.3 Plot the qprime profile
c
c 7.3.3.1 Plot the points and store in profle
      do 310 j   = 1,jpsi2
      jp         = j
      sppt       =         psivmsh(jp)
      qppt       =         qprim1 (jp)
      profle(jp) =         qppt
      call pointc(chqp,sppt,qppt,numbpc,incpcx,incpcy,delpcx,delpcy)
  310 continue
c
c 7.3.3.2 Plot the curve
      call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
c
c
c 7.4 Plot the profile for the modified qprime
c
c 7.4.1 Plot the  vertical  axis labels for the modified derivative
c
c 7.4.1.1 Write the label itself
      sposlb    = spmax0  + qllbx1
      qlposlb   = artmin  + mspql0*delz0
      call setlch(sposlb,qlposlb,icslbql,iszlbql,iorlbql,-1)
      write(string,2220) qllabel
      call wrtstr(string,1)
c
c 7.4.1.2 Plot the scale factor
      sposscl  = sposlb    + qlscx1
      qlpossc  = qlposlb   + qlscy1
      call setlch(sposscl,qlpossc,icsscql,iszscql,iorscql,-1)
c
      if(ratioqq .gt. bigscale)      write(string,2230) ratioqq
      if(ratioqq .gt.   1.0    .and.  ratioqq .le. bigscale)
     &                               write(string,2235) ratioqq
      if(ratioqq .lt.   1.0   )      write(string,2240) ratioqq
c
      if(ratioqq .ne.   1.0   )   call wrtstr(string,1)
c
c 7.4.2 Reset the printing parameters for pointc
c
      call setpch(icschql,iszchql,iorchql,kspchql)
c
c 7.4.3 Plot the modified qprime profile
c
c 7.4.3.1 Plot the points and store in profle
      do 320 j   = 1,jpsi2
      jp         = j
      sppt       =         psivmsh(jp)
      qppt       =         qprim1 (jp)
c
      ppps       = (spmax0-sppt)/sprang0
      pppt       = amax1(ppps,smallno)
      if(npowr .lt. -1) qlvl   = alog10(amax1(qppt*qppt,smallno))
      if(npowr .eq. -1) qlvl   = qppt/(1.0 +  abs(alog10(pppt)))
      if(npowr .eq.  0) qlvl   = qppt*ppps
      if(npowr .gt.  0) qlvl   = qppt/(1.0 + (abs(1.0/pppt))**npowr)
      qlpt       = ratioqq*qlvl
      profle(jp) =         qlpt
      call pointc(chql,sppt,qlpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  320 continue
c
c 7.4.3.2 Plot the curve
      call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
c
c
c 7.5 Plot a line if the origin is included in the range
c
      if(lqprim .eq. 0  .or.  lqplog .eq. 0) then
        call linep(spmin0,zero,spmax0,zero,kspzero)
      endif
c
c
c
c 8.0 Close frame, return, and end
c
      call frame(0)
      return
 1000 format(1x,'Plot of q fqpi and qprime vs',1x,a8)
 2000 format(f8.4)
 2010 format(a8)
 2100 format(f7.2)
 2102 format(f8.2)
 2105 format(f9.2)
 2107 format(f16.2)
 2110 format(a3)
 2120 format(a3)
 2130 format('(x',f6.2,')')
 2140 format('(x',f6.4,')')
 2200 format(f8.2)
 2202 format(f9.2)
 2205 format(f10.2)
 2207 format(f16.2)
 2210 format(a7)
 2220 format(a7)
 2230 format('(x',e10.3,')')
 2235 format('(x',f10.2,')')
 2240 format('(x',f10.6,')')
      end
      subroutine plotrho
c
c -------------------------------------------------------------
c
c plot radial mesh profiles
c
c -------------------------------------------------------------
c
      parameter (npx=200,ncx=2*npx,np1=npx+1,np2=npx+2,nc1=ncx+1)
      parameter (npp=929,ntt=2*npp-1,npb=npp)
      parameter (npq=np2)
      parameter (nlx=1441)
      parameter (nvn=7)
c
      character*1   chpn, chsq, chrh, chro,
     &              chtf, chst, chvl, chsv
      character*8   qvalue
      character*8   splabel
      character*8   pnlabel, sqlabel, rhlabel, rolabel,
     &              tflabel, stlabel, vllabel, svlabel
c
      character*132 string
c
      character*1   lbchpv, lbchsp, lbchrh, lbchro,
     &              lbchtf, lbchst ,lbchvl, lbchsv
      character*8   pclab,  tclab,  labpsi, 
     &              labpsiv,labpsin,labsqtp,
     &              labrhov,labrho2,labrhon,
     &              labtorf,labtorn,labsqtt,
     &              labvolm,labvoln,labsqtv,
     &              labdrh2,labdrho,labdtor,
     &              labdsqt,labdvol,labdsqv
      character*16  filout,filinp,fileql,filmap,filgrd,filrdm,filrdn,
     &              filfst,filplt
c
      common/cnst/  pi,hlfpi,twopi,amu,boltzm,echarg,pmass,roundff,bigno
      common/inpd/  ntor,ncase,idnsty,ndnxp0,ndnxp1,ndnxp2,ifastp,
     &              nmap,neqtyp,mpreset,ndoublt,ndivert,nmtype,
     &              jpsi,itht,jpsi1,jpsi2,ithtp,jsep,
     &              isym,igrid,nham1,nham2,nham3,npowr,
     &              bfieldf,rdefolt,qxin,btdes,qsurf,
     &              gamma,gamav,rmantl,betaf,zeffect,
     &              ncorr,corrfac,nccellr,peqpk0,peqpk1,peqpk2
      common/mapd/  mapmaxd,dpsisl,dpsisd,nqaxis,nwtmag,nfitmax,nfitpts,
     &              ifitrng,jfitrng,jfitchk,fitchek,cnvmag,epsaxs,
     &              maxerlp,maxerlv,delbox,delboz,delac,delav,
     &              delstsf,delstlp,delstlv,prfrac,
     &              nerstop,nerprnt,qptol,tolspln,tolbchi,tolbtor,
     &              tolsymm,tolaugm,errsep,precisn,plossmx,
     &              narcmx,ntrymx,ntdecr,ntmmin,npfit,npcmin,
     &              kuttaop,nrkmax0,nrkmax1,numstp,nwtfitp,nwtfitm,
     &              stepfac,flxstp,psispl,tolnwtp,tolnwtm,
     &              delpakf,delpakc,delpkf,delpkc,psichek,boxtnd,
     &              maptrace,norient,maxcutc,
     &              dresolv,dlclose,pntshft,endtol,
     &              narcln,nangax,nanglm,nbpmax,nwtmax,nslmax,nhfmax,
     &              bperor,sersnm,sertnm,arcmin,delgap,stepcut,
     &              novap,ioutm,iplotm
      common/sarc/  ntmax,ntmsh,npcf,npcb,npc,xp(nlx),zp(nlx),gsq(nlx),
     &              arc(nlx),tp(nlx),arcc(nlx),tpp(nlx),bcd(4),
     &              csx(3,nlx),csz(3,nlx),cseq1(3,nlx),cseq2(3,nlx),
     &              st1(nlx),st2(nlx),st3(nlx),csveq(3,nlx),
     &              sv0(nlx),sv1(nlx),sv2(nlx),sv3(nlx),sv4(nlx),
     &              sv5(nlx)
      common/ort1/  ndim1,ndim2,xma,zma,rcnt,btor,totcur,
     &              psimax0,psilim0,delpsir0,psimax1,psilim1,delpsir1,
     &              psimax,psilim,delpsir,dpsidec,dpsisv,dpedg1,dpedg2,
     &              eaxe,taxe,qaxe,fqiaxe,fqmaxe,paxe,faxe,ppaxe,ffpaxe,
     &              dnaxe,pfaxe,qppaxe,dlnaxe,dpnaxe,dhnaxe,
     &              arcaxe,dnnorm,pfnorm,
     &              qlim,fqilim,fqmlim,plim,flim,pplim,ffplim,
     &              dnlim,pflim,xlim,zlim,
     &              qmer,btmer,gkappa,arcprev,chicmax,
     &              charcmx,chpstmx,chhammx,arcnrm0,pstnrm0,hamnrm0,
     &              elongax1,elongax2,elongaxs,elongaxp,elongtn(npp)
      common/ort2/  psival(np1),pp(np1),ffp(np1),p(np1),f(np1),
     &              dnsty(np1),pfast(np1),qp(np1),fqpi(np1),fqim(np1),
     &              qdel(np1),psinitl(np1),qpinitl(np1),
     &              qprime(np1),qpint(np1),fqint(np1),
     &              psivl0(np2),rst(np2),zst(np2),rstp(np2),zstp(np2),
     &              rsinitl(np2),zsinitl(np2),
     &              arcsurf(nc1),pestchi(nc1),hamlchi(nc1),
     &              alfarc (nc1),alfpst (nc1),alfham (nc1),
     &              rs(nc1,np1),zs(nc1,np1)
      common/ort3/  rsrf(nc1),zsrf(nc1),chic(nc1),chie(nc1),dydx(nc1),
     &              xjsc(nc1)
      common/ort4/  rcell  (np1,nc1),zcell  (np1,nc1),
     &              dpsdr  (np1,nc1),dpsdz  (np1,nc1),
     &              chiarcl(np1,nc1),chipest(np1,nc1),chihaml(np1,nc1),
     &              xjcarcl(np1,nc1),xjcpest(np1,nc1),xjchaml(np1,nc1),
     &              alfarcl(np1,nc1),alfpest(np1,nc1),alfhaml(np1,nc1),
     &              arcnorm(np1),    pstnorm(np1),    hamnorm(np1),
     &              dlnlval(np1),    dpnlval(np1),    dhnlval(np1)
      common/ort5/  f3 (np1,nc1),f4 (np1,nc1),f5 (np1,nc1),
     &              f7 (np1,nc1),f8 (np1,nc1),f9 (np1,nc1),
     &              f10(np1,nc1),f11(np1,nc1),f12(np1,nc1),
     &              f13(np1,nc1),f14(np1,nc1),f15(np1,nc1),
     &              f16(np1,nc1),f17(np1,nc1),f18(np1,nc1),
     &              f19(np1,nc1),f20(np1,nc1),f21(np1,nc1),
     &              f22(np1,nc1),f23(np1,nc1),f24(np1,nc1),
     &              f25(np1,nc1),f26(np1,nc1),f27(np1,nc1),
     &              f28(np1,nc1)
      common/ort6/  savgax,tavgax,sntax0,sntax1,sntax2,sntax3,sntax4,
     &              svint(nvn),savge(np1),tavge(np1),sint0(np1),
     &              sint1(np1),sint2(np1),sint3(np1),sint4(np1)
      common/auxmsh/rh2tot,rhotot,flxtot,voltot,
     &              anltflxp,dvrtflxp,anltvolp,dvrtvolp,
     &              analtflx,divrtflx,analtvol,divrtvol,
     &              psivl1 (np2),qp1    (np2),fqpi1  (np2),
     &              qprim1 (np2),sintp0 (np2),psivmsh(np2),
     &              psivalu(np2),psinorm(np2),psisqrt(np2),
     &              psirho2(np2),psirhov(np2),psinrho(np2),
     &              psitorf(np2),psintor(np2),psisqtf(np2),
     &              psivolm(np2),psinvol(np2),psisqvl(np2),
     &              dpsirh2(np2),dpsirho(np2),dpsinrh(np2),
     &              dpsitor(np2),dpsintf(np2),dpsisqt(np2),
     &              dpsivol(np2),dpsinvl(np2),dpsisqv(np2)
       common/ratnl/jqmin, jqmax, numbqs,nq1totl,
     &              rcentr,zcentr,rminor,zminor,
     &              asprat00,asprat01,asprat10,asprat11,
     &              pminvl,qminvl,qpminv,pmaxvl,qmaxvl,qpmaxv,
     &              psivlmnq,    psivnmnq,    psisqmnq,
     &              psir2mnq,    psirhmnq,    psinrmnq,
     &              psitfmnq,    psintmnq,    psistmnq,
     &              psivmmnq,    psinvmnq,    psisvmnq,
     &              dpsr2mnq,    dpsrhmnq,    dpsnrmnq,
     &              dpstfmnq,    dpsntmnq,    dpsstmnq,
     &              dpsvmmnq,    dpsnvmnq,    dpssvmnq,
     &              psivlmxq,    psivnmxq,    psisqmxq,
     &              psir2mxq,    psirhmxq,    psinrmxq,
     &              psitfmxq,    psintmxq,    psistmxq,
     &              psivmmxq,    psinvmxq,    psisvmxq,
     &              dpsr2mxq,    dpsrhmxq,    dpsnrmxq,
     &              dpstfmxq,    dpsntmxq,    dpsstmxq,
     &              dpsvmmxq,    dpsnvmxq,    dpssvmxq,
     &              lpsiq  (npq),psivlq (npq),
     &              qprimq (npq),qvalue (npq),psimshq(npq),
     &              psivalq(npq),psinrmq(npq),psisqrq(npq),
     &              psirh2q(npq),psirhoq(npq),psinrhq(npq),
     &              psitorq(npq),psintfq(npq),psisqtq(npq),
     &              psivolq(npq),psinvlq(npq),psisqvq(npq),
     &              dpsir2q(npq),dpsirhq(npq),dpsinrq(npq),
     &              dpsitfq(npq),dpsintq(npq),dpsistq(npq),
     &              dpsivlq(npq),dpsinvq(npq),dpsisvq(npq),
     &              shearps(npq),shearrh(npq),
     &              sheartf(npq),shearvl(npq),
     &              epslrh1 (npq),shearrh1(npq),shearfrh(npq),
     &              epslvl1 (npq),shearvl1(npq),shearfvl(npq)
      common/labels/lbchpv, lbchsp, lbchrh, lbchro,
     &              lbchtf, lbchst, lbchvl, lbchsv,
     &              pclab,  tclab,  labpsi,
     &              labpsiv,labpsin,labsqtp,
     &              labrhov,labrho2,labrhon,
     &              labtorf,labtorn,labsqtt,
     &              labvolm,labvoln,labsqtv,
     &              labdrh2,labdrho,labdtor,
     &              labdsqt,labdvol,labdsqv
      common/pldf/  x0min,x0max,y0min,y0max 
      common/inout/ kutty,kuout,kuinp,kueql,kumap,kugrd,kurdm,kurdn,
     &              kufst,krdeqlb,krdmgta,krdngta,krdfast
      common/flnm/  filout,filinp,fileql,filmap,filgrd,filrdm,filrdn,
     &              filfst,filplt
c
      dimension profle(nlx)
      equivalence (st1 (1),profle(1))
c
c -------------------------------------------------------------
c
c
c
c 1.0 Initialization
c
c 1.1 Return if no plot is required
c
      if(iplotm .lt. 7) return
c
c
c 1.2 Constants
c
      bignum    = bigno
      smallno   = abs(roundff)
c
c
c 1.3 Plot dimensions
c
      x0mn      = 0.120
      x0mx      = 0.820
      y0mn      = 0.200
      y0mx      = 0.900
c
c
c 1.4 Heading and Legend
c
c 1.4.1 Heading position
c
      xxhd      = 0.0625
      yyhd      = 0.9500
c
c 1.4.2 Heading: case, size, and orientation
c
      icashd    = +2
      isizhd    = +3
      iornhd    =  0
c
c 1.4.3 Legend position
c
      xxlg      = 0.1500
      yylg      = 0.8625
c
c 1.4.4 Legend: case, size, and orientation
c
      icaslg    = +2
      isizlg    = +1
      iornlg    =  0
c
c
c 1.5 Label position parameters
c
c 1.5.1 Position increments
c
      nspx0     = 16
      nspy0     = 30
      nspz0     = 30
c
c 1.5.2 Axis label position fractions
c
      mspps0    = (2*nspx0)/3
      msppn0    = (5*nspy0)/6
      mspsq0    = (3*nspy0)/4
      msprh0    = (2*nspy0)/3
      mspro0    = (1*nspy0)/2
      msptf0    = (5*nspz0)/6
      mspst0    = (3*nspz0)/4
      mspvl0    = (2*nspz0)/3
      mspsv0    = (1*nspz0)/2
c
c
c 1.6 Default parameters for pointc
c
      numbpc    = 1
      incpcx    = 0
      incpcy    = 0
      delpcx    = 0.0
      delpcy    = 0.0
c
c
c
c 2.0 Axis settings
c
c 2.1 Set the normalized tick parameters
c
c 2.1.1 Tick scale skipping
c
      ntikxs0   = 2
      ntikys0   = 1
      ntikzs0   = 1
c
c 2.1.2 Axis extensions
c
      xtndl     = 1.200
      xtndr     = 1.200
c
c 2.1.3 Normalized tick lengths
c
      xtiky0    = +0.250
      ytikx0    = +0.250
      ztikx0    = -0.250
c
c 2.1.4 Tick scale positions
c
c 2.1.4.1 Label shift along relevant axis
      xsclx0    = +0.300
      yscly0    = -0.250
      zsclz0    = -0.250
c
c 2.1.4.2 Label shift away from relevant axis
      xscly0    = -4.500
      ysclx0    = -1.800
      zsclx0    = +0.020
c
c
c 2.2 Set the axis scales: case, size, and orientation
c
c 2.2.1 Horizontal axis
c
      icshsc    = +1
      iszhsc    = +1
      iorhsc    = +1
c
c 2.2.2 Left Vertical axis
c
      icsvsy    = +1
      iszvsy    = +1
      iorvsy    =  0
c
c 2.2.3 Right Vertical axis
c
      icsvsz    = +1
      iszvsz    = +1
      iorvsz    =  0
c
c
c 2.3 Axis labels specific to quantities plotted
c
c 2.3.1 Horizontal axis
c
c 2.3.1.1 Label and label position
      splby10   = -4.80
      splabel   = pclab
c
c 2.3.1.2 Horizontal axis label: case, size, and orientation
      icshlb    = +1
      iszhlb    = +3
      iorhlb    =  0
c
c 2.3.2 Left  Vertical axis
c
c 2.3.2.1 Label and label positions
c
c 2.3.2.1.1 Label and label position for psinorm
      pnlbx10   = -2.900
      pnlabel   = labpsin
c
c 2.3.2.1.2 Label and label and scale position for psisqrt
      sqlbx10   = -2.900
      sqscx10   = -0.5000
      sqscy10   = -0.7500
      sqlabel   = labsqtp
c
c 2.3.2.1.3 Label and label and scale position for psirho2
      rhlbx10   = -2.900
      rhscx10   = -0.5000
      rhscy10   = -0.7500
      rhlabel   = labrho2
c
c 2.3.2.1.4 Label and label and scale position for psirhov
      rolbx10   = -2.900
      roscx10   = -0.5000
      roscy10   = -0.7500
      rolabel   = labrhov
c
c 2.3.2.2 Left Vertical axis label and scale characteristics
c
c 2.3.2.2.1 Axis label: case, size, and orientation for psinorm
      icslbpn   = +1
      iszlbpn   = +2
      iorlbpn   =  0
c
c 2.3.2.2.2 Case, size, and orientation for psisqrt
c
c 2.3.2.2.2.1 Axis label
      icslbsq   = +1
      iszlbsq   = +2
      iorlbsq   =  0
c
c 2.3.2.2.2.2 Axis scale
      icsscsq   = +1
      iszscsq   = +1
      iorscsq   =  0
c
c 2.3.2.2.3 Case, size, and orientation for psirho2
c
c 2.3.2.2.3.1 Axis label
      icslbrh   = +1
      iszlbrh   = +2
      iorlbrh   =  0
c
c 2.3.2.2.3.2 Axis scale
      icsscrh   = +1
      iszscrh   = +1
      iorscrh   =  0
c
c 2.3.2.2.4 Case, size, and orientation for psirhov
c
c 2.3.2.2.4.1 Axis label
      icslbro   = +1
      iszlbro   = +2
      iorlbro   =  0
c
c 2.3.2.2.4.2 Axis scale
      icsscro   = +1
      iszscro   = +1
      iorscro   =  0
c
c 2.3.3 Right Vertical axis
c
c 2.3.3.1 Label and label position 
c
c 2.3.3.1.1 Label and label position for psintor
      tflbx10   = +1.0000
      tflabel   = labtorf
c
c 2.3.3.1.2 Label and label and scale position for psisqtf
      stlbx10   = +1.0000
      stscx10   = -0.5000
      stscy10   = -0.7500
      stlabel   = labsqtv
c
c 2.3.3.1.3 Label and label and scale position for psinvol
      vllbx10   = +1.0000
      vlscx10   = -0.5000
      vlscy10   = -0.7500
      vllabel   = labvolm
c
c 2.3.3.1.4 Label and label and scale position for psisqvl
      svlbx10   = +1.0000
      svscx10   = -0.5000
      svscy10   = -0.7500
      svlabel   = labsqtv
c
c 2.3.3.2 Right Vertical axis label and scale characteristics
c
c 2.3.3.2.1 Axis label: case, size, and orientation for psintor
      icslbtf   = +1
      iszlbtf   = +2
      iorlbtf   =  0
c
c 2.3.3.2.2 Case, size, and orientation for psisqtf
c
c 2.3.3.2.2.1 Axis label
      icslbst   = +1
      iszlbst   = +2
      iorlbst   =  0
c
c 2.3.3.2.2.2 Axis scale
      icsscst   = +1
      iszscst   = +1
      iorscst   =  0
c
c 2.3.3.2.3 Case, size, and orientation for psinvol
c
c 2.3.3.2.3.1 Axis label
      icslbvl   = +1
      iszlbvl   = +2
      iorlbvl   =  0
c
c 2.3.3.2.3.2 Axis scale
      icsscvl   = +1
      iszscvl   = +1
      iorscvl   =  0
c
c 2.3.3.2.4 Case, size, and orientation for psisqvl
c
c 2.3.3.2.4.1 Axis label
      icslbsv   = +1
      iszlbsv   = +2
      iorlbsv   =  0
c
c 2.3.3.2.4.2 Axis scale
      icsscsv   = +1
      iszscsv   = +1
      iorscsv   =  0
c
c
c 2.4 Line type parameters
c
c 2.4.1 psinorm profile
c
      icschpn   = +2
      iszchpn   = +1
      iorchpn   =  0
      kspchpn   = +2
      chpn      = lbchpv
c
c 2.4.2 psisqrt profile
c
      icschsq   = +2
      iszchsq   = +1
      iorchsq   =  0
      kspchsq   = +2
      chsq      = lbchsp
c
c 2.4.3 psirho2 profile
c
      icschrh   = +2
      iszchrh   = +1
      iorchrh   =  0
      kspchrh   = +2
      chrh      = lbchrh
c
c 2.4.4 psirhov profile
c
      icschro   = +2
      iszchro   = +1
      iorchro   =  0
      kspchrh   = +2
      chro      = lbchro
c
c 2.4.5 psintor profile
c
      icschtf   = +2
      iszchtf   = +1
      iorchtf   =  0
      kspchtf   = +2
      chtf      = lbchtf
c
c 2.4.6 psisqtf profile
c
      icschst   = +2
      iszchst   = +1
      iorchst   =  0
      kspchst   = +2
      chst      = lbchst
c
c 2.4.7 psinvol profile
c
      icschvl   = +2
      iszchvl   = +1
      iorchvl   =  0
      kspchvl   = +2
      chvl      = lbchvl
c
c 2.4.8 psisqvl profile
c
      icschsv   = +2
      iszchsv   = +1
      iorchsv   =  0
      kspchsv   = +2
      chsv      = lbchsv
c
c
c
c 3.0 Find the function minima and maxima
c
c 3.1 Search the physically scaled data
c
      spmin      = psivmsh(  1  )
      spmax      = psivmsh(jpsi2)
      pnmin      = +bignum
      pnmax      = -bignum
      sqmin      = +bignum
      sqmax      = -bignum
      rhmin      = +bignum
      rhmax      = -bignum
      romin      = +bignum
      romax      = -bignum
      tfmin      = +bignum
      tfmax      = -bignum
      stmin      = +bignum
      stmax      = -bignum
      vlmin      = +bignum
      vlmax      = -bignum
      svmin      = +bignum
      svmax      = -bignum
      sprange    =  spmax-spmin
      do 10  j   = 1,jpsi2
      jp0        = j
      psvlu      = psivmsh(jp0)
      pnvlu      = psinorm(jp0)
      sqvlu      = psisqrt(jp0)
      rhvlu      = psirho2(jp0)
      rovlu      = psirhov(jp0)
      tfvlu      = psintor(jp0)
      stvlu      = psisqtf(jp0)
      vlvlu      = psinvol(jp0)
      svvlu      = psisqvl(jp0)
c
      pnmin      = amin1(pnmin,pnvlu)
      pnmax      = amax1(pnmax,pnvlu)
      sqmin      = amin1(sqmin,sqvlu)
      sqmax      = amax1(sqmax,sqvlu)
      rhmin      = amin1(rhmin,rhvlu)
      rhmax      = amax1(rhmax,rhvlu)
      romin      = amin1(romin,rovlu)
      romax      = amax1(romax,rovlu)
      tfmin      = amin1(tfmin,tfvlu)
      tfmax      = amax1(tfmax,tfvlu)
      stmin      = amin1(stmin,stvlu)
      stmax      = amax1(stmax,stvlu)
      vlmin      = amin1(vlmin,vlvlu)
      vlmax      = amax1(vlmax,vlvlu)
      svmin      = amin1(svmin,svvlu)
      svmax      = amax1(svmax,svvlu)
   10 continue
c
c
c 3.2 Adjust the maxima as needed
c
c 3.2.1 Horizontal axis extrema
c
      spmin0     = spmin
      spmax0     = spmax
      sprang0    = spmax0 - spmin0
c
c 3.2.2 Vertical axis extrema
c       All profiles but psirho2 and psirhov should be normalized between zero and 1.0
c
      pnmin0  = amin1(pnmin,0.0)
      pnmax0  = amax1(pnmax,1.0)
c
      sqmin0  = amin1(sqmin,0.0)
      sqmax0  = amax1(sqmax,1.0)
c
      rhmin0  = amin1(rhmin,0.0)
      rhmax0  = amax1(rhmax,rh2tot)
c
      romin0  = amin1(romin,0.0)
      romax0  = amax1(romax,rhotot)
c
      tfmin0  = amin1(tfmin,0.0)
      tfmax0  = amax1(tfmax,1.0)
c
      stmin0  = amin1(stmin,0.0)
      stmax0  = amax1(stmax,1.0)
c
      vlmin0  = amin1(vlmin,0.0)
      vlmax0  = amax1(vlmax,1.0)
c
      svmin0  = amin1(svmin,0.0)
      svmax0  = amax1(svmax,1.0)
c
c
c 3.3 Set the range parameters
c
c 3.3.1 Range for psinorm
c
      if(pnmin0 .ge. 0.0  .and.  pnmax0 .gt. 0.0) lpnval  = +1
      if(pnmin0 .lt. 0.0  .and.  pnmax0 .gt. 0.0) lpnval  =  0
      if(pnmin0 .lt. 0.0  .and.  pnmax0 .le. 0.0) lpnval  = -1
c
c 3.3.2 Range for psisqrt
c
      if(sqmin0 .ge. 0.0  .and.  sqmax0 .gt. 0.0) lsqval  = +1
      if(sqmin0 .lt. 0.0  .and.  sqmax0 .gt. 0.0) lsqval  =  0
      if(sqmin0 .lt. 0.0  .and.  sqmax0 .le. 0.0) lsqval  = -1
c
c 3.3.3 Range for psirho2
c
      if(rhmin0 .ge. 0.0  .and.  rhmax0 .gt. 0.0) lrhval  = +1
      if(rhmin0 .lt. 0.0  .and.  rhmax0 .gt. 0.0) lrhval  =  0
      if(rhmin0 .lt. 0.0  .and.  rhmax0 .le. 0.0) lrhval  = -1
c
c 3.3.4 Range for psirhov
c
      if(romin0 .ge. 0.0  .and.  romax0 .gt. 0.0) lroval  = +1
      if(romin0 .lt. 0.0  .and.  romax0 .gt. 0.0) lroval  =  0
      if(romin0 .lt. 0.0  .and.  romax0 .le. 0.0) lroval  = -1
c
c 3.3.5 Range for psintor
c
      if(tfmin0 .ge. 0.0  .and.  tfmax0 .gt. 0.0) ltfval  = +1
      if(tfmin0 .lt. 0.0  .and.  tfmax0 .gt. 0.0) ltfval  =  0
      if(tfmin0 .lt. 0.0  .and.  tfmax0 .le. 0.0) ltfval  = -1
c
c 3.3.6 Range for psisqtf
c
      if(stmin0 .ge. 0.0  .and.  stmax0 .gt. 0.0) lstval  = +1
      if(stmin0 .lt. 0.0  .and.  stmax0 .gt. 0.0) lstval  =  0
      if(stmin0 .lt. 0.0  .and.  stmax0 .le. 0.0) lstval  = -1
c
c 3.3.7 Range for psinvol
c
      if(vlmin0 .ge. 0.0  .and.  vlmax0 .gt. 0.0) lvlval  = +1
      if(vlmin0 .lt. 0.0  .and.  vlmax0 .gt. 0.0) lvlval  =  0
      if(vlmin0 .lt. 0.0  .and.  vlmax0 .le. 0.0) lvlval  = -1
c
c 3.3.8 Range for psisqvl
c
      if(svmin0 .ge. 0.0  .and.  svmax0 .gt. 0.0) lsvval  = +1
      if(svmin0 .lt. 0.0  .and.  svmax0 .gt. 0.0) lsvval  =  0
      if(svmin0 .lt. 0.0  .and.  svmax0 .le. 0.0) lsvval  = -1
c
c 3.3.9 Check the ranges are valid for normalized profiles
c
      if(lpnval .le. 0  .or.  lsqval .le. 0) then
        call abortjob
     &        ('plotrho ',  1,   'Poloidal flux scale non positive    '
     &        ,'lpnval  ', lpnval,   'lsqval  ', lsqval,   -1)
      endif
c
      if(lrhval .le. 0  .or.  lroval .le. 0) then
        call abortjob
     &        ('plotrho ',  2,   'Poloidal rho  scale non positive    '
     &        ,'lrhval  ', lrhval,   'lroval  ', lroval,   -1)
      endif
c
      if(ltfval .le. 0  .or.  lstval .le. 0) then
        call abortjob
     &        ('plotrho ',  3,   'Toroidal flux scale non positive    '
     &        ,'ltfval  ', ltfval,   'lstval  ', lstval,   -1)
      endif
c
      if(lvlval .le. 0  .or.  lsvval .le. 0) then
        call abortjob
     &        ('plotrho ',  4,   'Flux volume  scale non positive     '
     &        ,'lvlval  ', lvlval,   'lsvval  ', lsvval,   -1)
      endif
c
c
c 3.4 Vertical axis scale factors
c
      ratiosq  = sqmax0/pnmax0
      ratiorh  = rhmax0/pnmax0
      ratioro  = romax0/pnmax0
      ratiotf  = tfmax0/pnmax0
      ratiost  = stmax0/pnmax0
      ratiovl  = vlmax0/pnmax0
      ratiosv  = svmax0/pnmax0
c
      if(ratiosq .le. 0.0) then
        call abortjob
     &        ('plotrho ',  5,   'Nonpositive scale for psisqrt       '
     &        ,'lpnval  ', lpnval,   'lsqval  ', lsqval,   -1)
        if(abs(ratiosq) .lt. smallno) ratiosq  = -1.0
      endif
c
      if(ratiorh .le. 0.0) then
        call abortjob
     &        ('plotrho ',  6,   'Nonpositive scale for psirho2       '
     &        ,'lpnval  ', lpnval,   'lrhval  ', lrhval,   -1)
        if(abs(ratiorh) .lt. smallno) ratiorh  = -1.0
      endif
c
      if(ratioro .le. 0.0) then
        call abortjob
     &        ('plotrho ',  7,   'Nonpositive scale for psirhov       '
     &        ,'lpnval  ', lpnval,   'lroval  ', lroval,   -1)
        if(abs(ratioro) .lt. smallno) ratioro  = -1.0
      endif
c
      if(ratiotf .le. 0.0) then
        call abortjob
     &        ('plotrho ',  8,   'Nonpositive scale for psintor       '
     &        ,'lpnval  ', lpnval,   'ltfval  ', ltfval,   -1)
        if(abs(ratiotf) .lt. smallno) ratiotf  = -1.0
      endif
c
      if(ratiost .le. 0.0) then
        call abortjob
     &        ('plotrho ',  9,   'Nonpositive scale for psisqtf       '
     &        ,'lpnval  ', lpnval,   'lstval  ', lstval,   -1)
        if(abs(ratiost) .lt. smallno) ratiost  = -1.0
      endif
c
      if(ratiovl .le. 0.0) then
        call abortjob
     &        ('plotrho ', 10,   'Nonpositive scale for psinvol       '
     &        ,'lpnval  ', lpnval,   'lvlval  ', lvlval,   -1)
        if(abs(ratiovl) .lt. smallno) ratiovl  = -1.0
      endif
c
      if(ratiosv .le. 0.0) then
        call abortjob
     &        ('plotrho ', 11,   'Nonpositive scale for psisqvl       '
     &        ,'lpnval  ', lpnval,   'lsvval  ', lsvval,   -1)
        if(abs(ratiosv) .lt. smallno) ratiosv  = -1.0
      endif
c
      ratsqpn  =   1.0  /ratiosq
      ratrhpn  =   1.0  /ratiorh
      ratropn  =   1.0  /ratioro
c
      ratsttf  = ratiotf/ratiost
      ratvltf  = ratiotf/ratiovl
      ratsvtf  = ratiotf/ratiosv
c
c
c
c 4.0 Set the physical plot parameters
c
c 4.1 Set the map extrema from the extrema
c
c 4.1.1 Left scale extrema
c
      alfmin    =       pnmin0
      alfmax    = xtndl*pnmax0
c
      alfmin    = amin1(alfmin,ratsqpn*sqmin0)
      alfmin    = amin1(alfmin,ratrhpn*rhmin0)
      alfmin    = amin1(alfmin,ratropn*romin0)
      alfmax    = amax1(alfmax,ratsqpn*sqmax0)
      alfmax    = amax1(alfmax,ratrhpn*rhmax0)
      alfmax    = amax1(alfmax,ratropn*romax0)
c
c 4.1.2 Right scale extrema
c
      artmin    =       tfmin0
      artmax    = xtndr*tfmax0
c
      artmin    = amin1(artmin,ratsttf*stmin0)
      artmin    = amin1(artmin,ratvltf*vlmin0)
      artmin    = amin1(artmin,ratsvtf*stmin0)
      artmax    = amax1(artmax,ratsttf*stmax0)
      artmax    = amax1(artmax,ratvltf*vlmax0)
      artmax    = amax1(artmax,ratsvtf*svmax0)
c
c
c 4.2 Set the number of tick marks
c
c 4.2.1 Horizontal axis
c
      nspx      = 10
c
c 4.2.2 Left  Vertical axis
c
      nspy1     = ifix(alfmax - alfmin)
      nspy1     = min0(nspy1,nspy0)
      nspy1     = max0(nspy1,1)
c
      if    (nspy1 .le. 4) then
        if    (nspy1 .le. 1) then
          nspy    = 5*nspy1
        elseif(nspy1 .le. 2) then
          nspy    = 4*nspy1
        else
          nspy    = 2*nspy1
        endif
c
      elseif(nspy1 .gt. 4) then
        nspy    =   nspy1
      endif
c
c 4.2.3 Right Vertical axis
c
      nspz1     = ifix(artmax - artmin)
      nspz1     = min0(nspz1,nspz0)
      nspz1     = max0(nspz1,1)
c
      if    (nspz1 .le. 4) then
        if    (nspz1 .le. 1) then
          nspz    = 5*nspz1
        elseif(nspz1 .le. 2) then
          nspz    = 4*nspz1
        else
          nspz    = 2*nspz1
        endif
c
      elseif(nspz1 .gt. 4) then
        nspz    =   nspz1
      endif
c
c
c 4.3 Normalize the physical spacing increments to the physical range
c
c 4.3.1 Standard spacing
c
      delx0     = (spmax0 - spmin0)/nspx0
      dely0     = (alfmax - alfmin)/nspy0
      delz0     = (artmax - artmin)/nspz0
c
c 4.3.2 Tick spacing
c
      delx      = (spmax0 - spmin0)/nspx
      dely      = (alfmax - alfmin)/nspy
      delz      = (artmax - artmin)/nspz
c
c
c 4.4 Set the physical spacings
c
c 4.4.1 Tick scale values
c
c 4.4.1.1 Horizontal axis
      ntikxs    = min0(ntikxs0,nspx)
      ntikxs    = max0(ntikxs,  1  )
      if(nspx .lt. 5) ntikxs  = 1
c
c 4.4.1.2 Left Vertical axis
      ntikys    = min0(ntikys0,nspy)
      ntikys    = max0(ntikys,  1  )
      if(nspy .lt. 5) ntikys  = 1
c
c 4.4.1.3 Right Vertical axis
      ntikzs    = min0(ntikzs0,nspz)
      ntikzs    = max0(ntikzs,  1  )
      if(nspz .lt. 5) ntikzs  = 1
c
c 4.4.2 Axis scale positions
c
c 4.4.2.1 Label shift along relevant axis
      xsclx     = xsclx0*delx0
      yscly     = yscly0*dely0
      zsclz     = zsclz0*delz0
c
c 4.4.2.2 Label shift away from relevant axis
      xscly     = xscly0*dely0
      ysclx     = ysclx0*delx0
      zsclx     = zsclx0*delx0
c
c 4.4.3 Tick lengths
c
      xtiky     = xtiky0*dely0
      ytikx     = ytikx0*delx0
      ztikx     = ztikx0*delx0
c
c 4.4.4 Label positions specific to profiles
c
c 4.4.4.1 Horizontal axis
      splby1    = splby10*dely0
c
c 4.4.4.2 Left Vertical Axis
      pnlbx1    = pnlbx10*delx0
      sqlbx1    = sqlbx10*delx0
      sqscx1    = sqscx10*delx0
      sqscy1    = sqscy10*dely0
      rhlbx1    = rhlbx10*delx0
      rhscx1    = rhscx10*delx0
      rhscy1    = rhscy10*dely0
      rolbx1    = rolbx10*delx0
      roscx1    = roscx10*delx0
      roscy1    = roscy10*dely0
c
c 4.4.4.3 Right Vertical Axis
      tflbx1    = tflbx10*delx0
      stlbx1    = stlbx10*delx0
      stscx1    = stscx10*delx0
      stscy1    = stscy10*delz0
      vllbx1    = vllbx10*delx0
      vlscx1    = vlscx10*delx0
      vlscy1    = vlscy10*delz0
      svlbx1    = svlbx10*delx0
      svscx1    = svscx10*delx0
      svscy1    = svscy10*delz0
c
c
c 4.5 Additional lines
c
      zero      = 0.0
      kspzero   =  2
c
c
c 5.0 Print out the plot heading
c
c 5.1 Set up the plot in normalized coordinates
c
      call map(x0min,x0max,y0min,y0max,x0min,x0max,y0min,y0max)
c
c
c 5.2 Plot the four borders
c
      call line(x0mn,y0mn,x0mn,y0mx)
      call line(x0mn,y0mx,x0mx,y0mx)
      call line(x0mx,y0mx,x0mx,y0mn)
      call line(x0mn,y0mn,x0mx,y0mn)
c
c
c 5.3 Set the position and plot the heading text
c
      call setlch(xxhd,  yyhd,  icashd,isizhd,iornhd,-1)
      write(string,1000) splabel 
      call wrtstr(string,1)
c
c
c 5.4 Set the position and plot the legend
c
      call setlch(xxlg,  yylg,  icaslg,isizlg,iornlg,-1)
      write(string,1100) pnlabel,chpn
      call wrtstr(string,1)
      write(string,1110) sqlabel,chsq
      call wrtstr(string,1)
      write(string,1120) rhlabel,chrh
      call wrtstr(string,1)
      write(string,1130) rolabel,chro
      call wrtstr(string,1)
c
      write(string,1140) tflabel,chtf
      call wrtstr(string,1)
      write(string,1150) stlabel,chst
      call wrtstr(string,1)
      write(string,1160) vllabel,chvl
      call wrtstr(string,1)
      write(string,1170) svlabel,chsv
      call wrtstr(string,1)
c
c
c
c 6.0 Plot psinorm, psisqrt, psirho2, and psirhov
c
c 6.1 Set the plot page in physical units (left axis)
c
      call map(spmin0,spmax0,alfmin,alfmax,x0mn,x0mx,y0mn,y0mx)
c
c
c 6.2 Set up the Horizontal axis
c
c 6.2.1 Plot the horizontal axis tick marks and scale
c
c 6.2.1.1 Set the base tick parameters
      vrtpos     = alfmin
      vrtpos1    = vrtpos
      vrtpos2    = vrtpos + xtiky
      vrtpossc   = vrtpos + xscly
      do 100 k   = 0,nspx
      kskpp      = (k/ntikxs)*ntikxs
      spos       = spmin0 + k*delx
      spossc     = spos   + xsclx
c
c 6.2.1.2 Plot the tick mark
      call line  (spos, vrtpos1, spos, vrtpos2 )
c
c 6.2.1.3 Plot the scale value
      if(k .eq. kskpp) then
        call setlch(spossc,vrtpossc,icshsc,iszhsc,iorhsc,-1)
        write(string,2000) spos
        call wrtstr(string,1)
      endif
  100 continue
c
c 6.2.2 Plot the horizontal axis label
c
      sposlb    = spmin0  + mspps0*delx0
      vrtposlb  = vrtpos  + splby1
      call setlch(sposlb,vrtposlb,icshlb,iszhlb,iorhlb,-1)
      write(string,2010) splabel
      call wrtstr(string,1)
c
c
c 6.3 Set up the Left  vertical axis
c
c 6.3.1 Plot the  vertical  axis tick marks and scale
c
c 6.3.1.1 Set the base tick parameters
      spos       = spmin0
      spos1      = spos
      spos2      = spos + ytikx
      spossc     = spos + ysclx
      do 110 i   = 0,nspy
      iskpp      = (i/ntikys)*ntikys
      alfpos     = alfmin  + i*dely
      alfpossc   = alfpos  + yscly
c
c 6.3.1.2 Plot the tick mark
      call line  (spos1, alfpos,  spos2, alfpos  )
c
c 6.3.1.3 Plot the scale value
      if(i .eq. iskpp) then
        call setlch(spossc,alfpossc,icsvsy,iszvsy,iorvsy,-1)
        write(string,2020) alfpos
        call wrtstr(string,1)
      endif
  110 continue
c
c
c 6.4 Plot the psinorm profile
c
c 6.4.1 Plot the  vertical  axis label
c
      sposlb   = spmin0  + pnlbx1
      pnposlb  = alfmin  + msppn0*dely0
      call setlch(sposlb,pnposlb,icslbpn,iszlbpn,iorlbpn,-1)
      write(string,2100) pnlabel
      call wrtstr(string,1)
c
c 6.4.2 Reset the printing parameters for pointc
c
      call setpch(icschpn,iszchpn,iorchpn,kspchpn)
c
c 6.4.3 Plot the profile
c
c 6.4.3.1 Plot the points and store in profle
      do 200 j   = 1,jpsi2
      jp         = j
      sppt       =         psivmsh(jp)
      pnpt       =         psinorm(jp)
      profle(jp) =         pnpt
      call pointc(chpn,sppt,pnpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  200 continue
c
c 6.4.3.2 Plot the curve
      call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
c
c
c 6.5 Plot the psisqrt profile
c
c 6.5.1 Plot the  vertical  axis label
c
c 6.5.1.1 Write the label itself
      sposlb   = spmin0  + sqlbx1
      sqposlb  = alfmin  + mspsq0*dely0
      call setlch(sposlb,sqposlb,icslbsq,iszlbsq,iorlbsq,-1)
      write(string,2200) sqlabel
      call wrtstr(string,1)
c
c 6.5.1.2 Plot the scale factor
      sposscs  = sposlb    + sqscx1
      sqpossc  = sqposlb   + sqscy1
      call setlch(sposscs,sqpossc,icsscsq,iszscsq,iorscsq,-1)
c
      if(ratsqpn .gt. 1.0) write(string,2210) ratsqpn
      if(ratsqpn .lt. 1.0) write(string,2220) ratsqpn
      if(ratsqpn .ne. 1.0) call wrtstr(string,1)
c
c 6.5.2 Reset the printing parameters for pointc
c
      call setpch(icschsq,iszchsq,iorchsq,kspchsq)
c
c 6.5.3 Plot the profile
c
c 6.5.3.1 Plot the points and store in profle
      do 210 j   = 1,jpsi2
      jp         = j
      sppt       =         psivmsh(jp)
      sqpt       =         psisqrt(jp)
      profle(jp) =         sqpt
      call pointc(chsq,sppt,sqpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  210 continue
c
c 6.5.3.2 Plot the curve
      call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
c
c
c 6.6 Plot the psirho2 profile
c
c 6.6.1 Plot the  vertical  axis label
c
c 6.6.1.1 Write the label itself
      sposlb   = spmin0  + rhlbx1
      rhposlb  = alfmin  + msprh0*dely0
      call setlch(sposlb,rhposlb,icslbrh,iszlbrh,iorlbrh,-1)
      write(string,2300) rhlabel
      call wrtstr(string,1)
c
c 6.6.1.2 Plot the scale factor
      sposscr  = sposlb    + rhscx1
      rhpossc  = rhposlb   + rhscy1
      call setlch(sposscr,rhpossc,icsscrh,iszscrh,iorscrh,-1)
c
      if(ratrhpn .gt. 1.0) write(string,2310) ratrhpn
      if(ratrhpn .lt. 1.0) write(string,2320) ratrhpn
      if(ratrhpn .ne. 1.0) call wrtstr(string,1)
c
c 6.6.2 Reset the printing parameters for pointc
c
      call setpch(icschrh,iszchrh,iorchrh,kspchrh)
c
c 6.6.3 Plot the renormalized values
c
c 6.6.3.1 Plot the points and store in profle
      do 220 j   = 1,jpsi2
      jp         = j
      sppt       =         psivmsh(jp)
      rhpt       = ratrhpn*psirho2(jp)
      profle(jp) =         rhpt
      call pointc(chrh,sppt,rhpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  220 continue
c
c 6.6.3.2 Plot the curve
      call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
c
c
c 6.7 Plot the psirhov profile
c
c 6.7.1 Plot the  vertical  axis label
c
c 6.7.1.1 Write the label itself
      sposlb   = spmin0  + rolbx1
      roposlb  = alfmin  + mspro0*dely0
      call setlch(sposlb,roposlb,icslbro,iszlbro,iorlbro,-1)
      write(string,2400) rolabel
      call wrtstr(string,1)
c
c 6.7.1.2 Plot the scale factor
      sposscu  = sposlb    + roscx1
      ropossc  = roposlb   + roscy1
      call setlch(sposscu,ropossc,icsscro,iszscro,iorscro,-1)
c
      if(ratropn .gt. 1.0) write(string,2410) ratropn
      if(ratropn .lt. 1.0) write(string,2420) ratropn
      if(ratropn .ne. 1.0) call wrtstr(string,1)
c
c 6.7.2 Reset the printing parameters for pointc
c
      call setpch(icschro,iszchro,iorchro,kspchro)
c
c 6.7.3 Plot the renormalized values
c
c 6.7.3.1 Plot the points and store in profle
      do 230 j   = 1,jpsi2
      jp         = j
      sppt       =         psivmsh(jp)
      ropt       = ratropn*psirhov(jp)
      profle(jp) =         ropt
      call pointc(chro,sppt,ropt,numbpc,incpcx,incpcy,delpcx,delpcy)
  230 continue
c
c 6.7.3.2 Plot the curve
      call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
c
c
c 6.8 Plot a line if the origin is included in the range
c
      if(lpnval .eq. 0  .or.  lsqval .eq. 0  .or.
     &   lrhval .eq. 0  .or.  lroval .eq. 0) then
        call linep(spmin0,zero,spmax0,zero,kspzero)
      endif
c
c
c
c 7.0 Plot psintor, psisqtf, psivolm, and psisqvl
c
c 7.1 Renormalize the map for the derivative
c
      call map(spmin0,spmax0,artmin,artmax,x0mn,x0mx,y0mn,y0mx)
c
c
c 7.2 Set up the Right Vertical axis
c
c 7.2.1 Plot the  vertical  axis tick marks and scale
c
c 7.2.1.1 Set the base tick parameters
      spos       = spmax0
      spos1      = spos + ztikx
      spos2      = spos
      spossc     = spos + zsclx
      do 300 i   = 0,nspz
      iskpp      = (i/ntikzs)*ntikzs
      artpos     = artmin  + i*delz
      artpossc   = artpos  + zsclz
c
c 7.2.1.2 Plot the tick mark
      call line  (spos1, artpos,  spos2, artpos  )
c
c 7.2.1.3 Plot the scale value
      if(i .eq. iskpp) then
        call setlch(spossc,artpossc,icsvsz,iszvsz,iorvsz,-1)
        write(string,2500) artpos
        call wrtstr(string,1)
      endif
  300 continue
c
c
c 7.3 Plot the profile for psintor
c
c 7.3.1 Plot the  vertical  axis labels for psintor
      sposlb    = spmax0  + tflbx1
      tfposlb   = artmin  + msptf0*delz0
      call setlch(sposlb,tfposlb,icslbtf,iszlbtf,iorlbtf,-1)
      write(string,2600) tflabel
      call wrtstr(string,1)
c
c 7.3.2 Reset the printing parameters for pointc
c
      call setpch(icschtf,iszchtf,iorchtf,kspchtf)
c
c 7.3.3 Plot the profile
c
c 7.3.3.1 Plot the points and store in profle
      do 400 j   = 1,jpsi2
      jp         = j
      sppt       =         psivmsh(jp)
      tfpt       =         psintor(jp)
      profle(jp) =         tfpt
      call pointc(chtf,sppt,tfpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  400 continue
c
c 7.3.3.2 Plot the curve
      call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
c
c
c 7.4 Plot the profile for psisqtf
c
c 7.4.1 Plot the  vertical  axis labels for psisqtf
c
c 7.4.1.1 Write the label itself
      sposlb    = spmax0  + stlbx1
      stposlb   = artmin  + mspst0*delz0
      call setlch(sposlb,stposlb,icslbst,iszlbst,iorlbst,-1)
      write(string,2700) stlabel
      call wrtstr(string,1)
c
c 7.4.1.2 Plot the scale factor
      sposscw  = sposlb    + stscx1
      stpossc  = stposlb   + stscy1
      call setlch(sposscw,stpossc,icsscst,iszscst,iorscst,-1)
c
      if(ratsttf .gt. 1.0) write(string,2710) ratsttf
      if(ratsttf .lt. 1.0) write(string,2720) ratsttf
      if(ratsttf .ne. 1.0) call wrtstr(string,1)
c
c 7.4.2 Reset the printing parameters for pointc
c
      call setpch(icschst,iszchst,iorchst,kspchst)
c
c 7.4.3 Plot the profile
c
c 7.4.3.1 Plot the points and store in profle
      do 410 j   = 1,jpsi2
      jp         = j
      sppt       =         psivmsh(jp)
      stpt       = ratsttf*psisqtf(jp)
      profle(jp) =         stpt
      call pointc(chst,sppt,stpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  410 continue
c
c 7.4.3.2 Plot the curve
      call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
c
c
c 7.5 Plot the profile for psinvol
c
c 7.5.1 Plot the  vertical  axis labels for psinvol
c
c 7.5.1.1 Write the label itself
      sposlb    = spmax0  + vllbx1
      vlposlb   = artmin  + mspvl0*delz0
      call setlch(sposlb,vlposlb,icslbvl,iszlbvl,iorlbvl,-1)
      write(string,2800) vllabel
      call wrtstr(string,1)
c
c 7.5.1.2 Plot the scale factor
      sposscv  = sposlb    + vlscx1
      vlpossc  = vlposlb   + vlscy1
      call setlch(sposscv,vlpossc,icsscvl,iszscvl,iorscvl,-1)
c
      if(ratvltf .gt. 1.0) write(string,2810) ratvltf
      if(ratvltf .lt. 1.0) write(string,2820) ratvltf
      if(ratvltf .ne. 1.0) call wrtstr(string,1)
c
c 7.5.2 Reset the printing parameters for pointc
c
      call setpch(icschvl,iszchvl,iorchvl,kspchvl)
c
c 7.5.3 Plot the profile
c
c 7.5.3.1 Plot the points and store in profle
      do 420 j   = 1,jpsi2
      jp         = j
      sppt       =         psivmsh(jp)
      vlpt       = ratvltf*psinvol(jp)
      profle(jp) =         vlpt
      call pointc(chvl,sppt,vlpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  420 continue
c
c 7.5.3.2 Plot the curve
      call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
c
c
c 7.6 Plot the profile for psisqvl
c
c 7.6.1 Plot the  vertical  axis labels for psisqvl
c
c 7.6.1.1 Write the label itself
      sposlb    = spmax0  + svlbx1
      svposlb   = artmin  + mspsv0*delz0
      call setlch(sposlb,svposlb,icslbsv,iszlbsv,iorlbsv,-1)
      write(string,2900) svlabel
      call wrtstr(string,1)
c
c 7.6.1.2 Plot the scale factor
      sposscx  = sposlb    + svscx1
      svpossc  = svposlb   + svscy1
      call setlch(sposscx,svpossc,icsscsv,iszscsv,iorscsv,-1)
c
      if(ratsvtf .gt. 1.0) write(string,2910) ratsvtf
      if(ratsvtf .lt. 1.0) write(string,2920) ratsvtf
      if(ratsvtf .ne. 1.0) call wrtstr(string,1)
c
c 7.6.2 Reset the printing parameters for pointc
c
      call setpch(icschsv,iszchsv,iorchsv,kspchsv)
c
c 7.6.3 Plot the profile
c
c 7.6.3.1 Plot the points and store in profle
      do 430 j   = 1,jpsi2
      jp         = j
      sppt       =         psivmsh(jp)
      svpt       = ratsvtf*psisqvl(jp)
      profle(jp) =         svpt
      call pointc(chsv,sppt,svpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  430 continue
c
c 7.6.3.2 Plot the curve
      call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
c
c
c 7.7 Plot a line if the origin is included in the range
c
      if(ltfval .eq. 0  .or.  lstval .eq. 0  .or.
     &   lvlval .eq. 0  .or.  lsvval .eq. 0) then
        call linep(spmin0,zero,spmax0,zero,kspzero)
      endif
c
c
c
c 8.0 Close frame, return, and end
c
      call frame(0)
      return
 1000 format(1x,'Normalized radial profiles vs',1x,a8)
 1100 format('psinorm=',1x,a8,':',2x,a1)
 1110 format('psisqrt=',1x,a8,':',2x,a1)
 1120 format('psirho2=',1x,a8,':',2x,a1)
 1130 format('psirhov=',1x,a8,':',2x,a1)
 1140 format('psintor=',1x,a8,':',2x,a1)
 1150 format('psisqtf=',1x,a8,':',2x,a1)
 1160 format('psinvol=',1x,a8,':',2x,a1)
 1170 format('psisqvl=',1x,a8,':',2x,a1)
 2000 format(f8.4)
 2010 format(a8)
 2020 format(f6.2)
 2100 format(a8)
 2200 format(a8)
 2210 format('(x',f6.2,')')
 2220 format('(x',f6.4,')')
 2300 format(a6)
 2310 format('(x',f6.2,')')
 2320 format('(x',f6.4,')')
 2400 format(a3)
 2410 format('(x',f6.2,')')
 2420 format('(x',f6.4,')')
 2500 format(f6.2)
 2600 format(a8)
 2700 format(a8)
 2710 format('(x',f8.2,')')
 2720 format('(x',f6.4,')')
 2800 format(a8)
 2810 format('(x',f8.2,')')
 2820 format('(x',f6.4,')')
 2900 format(a8)
 2910 format('(x',f8.2,')')
 2920 format('(x',f6.4,')')
      end
      subroutine plotdrh
c
c -------------------------------------------------------------
c
c plot derivatives of psi with respect to radial mesh profiles
c
c -------------------------------------------------------------
c
      parameter (npx=200,ncx=2*npx,np1=npx+1,np2=npx+2,nc1=ncx+1)
      parameter (npp=929,ntt=2*npp-1,npb=npp)
      parameter (npq=np2)
      parameter (nlx=1441)
      parameter (nvn=7)
c
      character*1   chdrh, chdro, chdtf, chdst, chdvl, chdsv
      character*8   qvalue
      character*8   splabel
      character*8   drhlabel,drolabel
      character*8   dtflabel,dstlabel,dvllabel,dsvlabel
c
      character*132 string
c
      character*1   lbchpv, lbchsp, lbchrh, lbchro,
     &              lbchtf, lbchst ,lbchvl, lbchsv
      character*8   pclab,  tclab,  labpsi, 
     &              labpsiv,labpsin,labsqtp,
     &              labrhov,labrho2,labrhon,
     &              labtorf,labtorn,labsqtt,
     &              labvolm,labvoln,labsqtv,
     &              labdrh2,labdrho,labdtor,
     &              labdsqt,labdvol,labdsqv
      character*16  filout,filinp,fileql,filmap,filgrd,filrdm,filrdn,
     &              filfst,filplt
c
      common/cnst/  pi,hlfpi,twopi,amu,boltzm,echarg,pmass,roundff,bigno
      common/inpd/  ntor,ncase,idnsty,ndnxp0,ndnxp1,ndnxp2,ifastp,
     &              nmap,neqtyp,mpreset,ndoublt,ndivert,nmtype,
     &              jpsi,itht,jpsi1,jpsi2,ithtp,jsep,
     &              isym,igrid,nham1,nham2,nham3,npowr,
     &              bfieldf,rdefolt,qxin,btdes,qsurf,
     &              gamma,gamav,rmantl,betaf,zeffect,
     &              ncorr,corrfac,nccellr,peqpk0,peqpk1,peqpk2
      common/mapd/  mapmaxd,dpsisl,dpsisd,nqaxis,nwtmag,nfitmax,nfitpts,
     &              ifitrng,jfitrng,jfitchk,fitchek,cnvmag,epsaxs,
     &              maxerlp,maxerlv,delbox,delboz,delac,delav,
     &              delstsf,delstlp,delstlv,prfrac,
     &              nerstop,nerprnt,qptol,tolspln,tolbchi,tolbtor,
     &              tolsymm,tolaugm,errsep,precisn,plossmx,
     &              narcmx,ntrymx,ntdecr,ntmmin,npfit,npcmin,
     &              kuttaop,nrkmax0,nrkmax1,numstp,nwtfitp,nwtfitm,
     &              stepfac,flxstp,psispl,tolnwtp,tolnwtm,
     &              delpakf,delpakc,delpkf,delpkc,psichek,boxtnd,
     &              maptrace,norient,maxcutc,
     &              dresolv,dlclose,pntshft,endtol,
     &              narcln,nangax,nanglm,nbpmax,nwtmax,nslmax,nhfmax,
     &              bperor,sersnm,sertnm,arcmin,delgap,stepcut,
     &              novap,ioutm,iplotm
      common/sarc/  ntmax,ntmsh,npcf,npcb,npc,xp(nlx),zp(nlx),gsq(nlx),
     &              arc(nlx),tp(nlx),arcc(nlx),tpp(nlx),bcd(4),
     &              csx(3,nlx),csz(3,nlx),cseq1(3,nlx),cseq2(3,nlx),
     &              st1(nlx),st2(nlx),st3(nlx),csveq(3,nlx),
     &              sv0(nlx),sv1(nlx),sv2(nlx),sv3(nlx),sv4(nlx),
     &              sv5(nlx)
      common/ort1/  ndim1,ndim2,xma,zma,rcnt,btor,totcur,
     &              psimax0,psilim0,delpsir0,psimax1,psilim1,delpsir1,
     &              psimax,psilim,delpsir,dpsidec,dpsisv,dpedg1,dpedg2,
     &              eaxe,taxe,qaxe,fqiaxe,fqmaxe,paxe,faxe,ppaxe,ffpaxe,
     &              dnaxe,pfaxe,qppaxe,dlnaxe,dpnaxe,dhnaxe,
     &              arcaxe,dnnorm,pfnorm,
     &              qlim,fqilim,fqmlim,plim,flim,pplim,ffplim,
     &              dnlim,pflim,xlim,zlim,
     &              qmer,btmer,gkappa,arcprev,chicmax,
     &              charcmx,chpstmx,chhammx,arcnrm0,pstnrm0,hamnrm0,
     &              elongax1,elongax2,elongaxs,elongaxp,elongtn(npp)
      common/ort2/  psival(np1),pp(np1),ffp(np1),p(np1),f(np1),
     &              dnsty(np1),pfast(np1),qp(np1),fqpi(np1),fqim(np1),
     &              qdel(np1),psinitl(np1),qpinitl(np1),
     &              qprime(np1),qpint(np1),fqint(np1),
     &              psivl0(np2),rst(np2),zst(np2),rstp(np2),zstp(np2),
     &              rsinitl(np2),zsinitl(np2),
     &              arcsurf(nc1),pestchi(nc1),hamlchi(nc1),
     &              alfarc (nc1),alfpst (nc1),alfham (nc1),
     &              rs(nc1,np1),zs(nc1,np1)
      common/ort3/  rsrf(nc1),zsrf(nc1),chic(nc1),chie(nc1),dydx(nc1),
     &              xjsc(nc1)
      common/ort4/  rcell  (np1,nc1),zcell  (np1,nc1),
     &              dpsdr  (np1,nc1),dpsdz  (np1,nc1),
     &              chiarcl(np1,nc1),chipest(np1,nc1),chihaml(np1,nc1),
     &              xjcarcl(np1,nc1),xjcpest(np1,nc1),xjchaml(np1,nc1),
     &              alfarcl(np1,nc1),alfpest(np1,nc1),alfhaml(np1,nc1),
     &              arcnorm(np1),    pstnorm(np1),    hamnorm(np1),
     &              dlnlval(np1),    dpnlval(np1),    dhnlval(np1)
      common/ort5/  f3 (np1,nc1),f4 (np1,nc1),f5 (np1,nc1),
     &              f7 (np1,nc1),f8 (np1,nc1),f9 (np1,nc1),
     &              f10(np1,nc1),f11(np1,nc1),f12(np1,nc1),
     &              f13(np1,nc1),f14(np1,nc1),f15(np1,nc1),
     &              f16(np1,nc1),f17(np1,nc1),f18(np1,nc1),
     &              f19(np1,nc1),f20(np1,nc1),f21(np1,nc1),
     &              f22(np1,nc1),f23(np1,nc1),f24(np1,nc1),
     &              f25(np1,nc1),f26(np1,nc1),f27(np1,nc1),
     &              f28(np1,nc1)
      common/ort6/  savgax,tavgax,sntax0,sntax1,sntax2,sntax3,sntax4,
     &              svint(nvn),savge(np1),tavge(np1),sint0(np1),
     &              sint1(np1),sint2(np1),sint3(np1),sint4(np1)
      common/auxmsh/rh2tot,rhotot,flxtot,voltot,
     &              anltflxp,dvrtflxp,anltvolp,dvrtvolp,
     &              analtflx,divrtflx,analtvol,divrtvol,
     &              psivl1 (np2),qp1    (np2),fqpi1  (np2),
     &              qprim1 (np2),sintp0 (np2),psivmsh(np2),
     &              psivalu(np2),psinorm(np2),psisqrt(np2),
     &              psirho2(np2),psirhov(np2),psinrho(np2),
     &              psitorf(np2),psintor(np2),psisqtf(np2),
     &              psivolm(np2),psinvol(np2),psisqvl(np2),
     &              dpsirh2(np2),dpsirho(np2),dpsinrh(np2),
     &              dpsitor(np2),dpsintf(np2),dpsisqt(np2),
     &              dpsivol(np2),dpsinvl(np2),dpsisqv(np2)
       common/ratnl/jqmin, jqmax, numbqs,nq1totl,
     &              rcentr,zcentr,rminor,zminor,
     &              asprat00,asprat01,asprat10,asprat11,
     &              pminvl,qminvl,qpminv,pmaxvl,qmaxvl,qpmaxv,
     &              psivlmnq,    psivnmnq,    psisqmnq,
     &              psir2mnq,    psirhmnq,    psinrmnq,
     &              psitfmnq,    psintmnq,    psistmnq,
     &              psivmmnq,    psinvmnq,    psisvmnq,
     &              dpsr2mnq,    dpsrhmnq,    dpsnrmnq,
     &              dpstfmnq,    dpsntmnq,    dpsstmnq,
     &              dpsvmmnq,    dpsnvmnq,    dpssvmnq,
     &              psivlmxq,    psivnmxq,    psisqmxq,
     &              psir2mxq,    psirhmxq,    psinrmxq,
     &              psitfmxq,    psintmxq,    psistmxq,
     &              psivmmxq,    psinvmxq,    psisvmxq,
     &              dpsr2mxq,    dpsrhmxq,    dpsnrmxq,
     &              dpstfmxq,    dpsntmxq,    dpsstmxq,
     &              dpsvmmxq,    dpsnvmxq,    dpssvmxq,
     &              lpsiq  (npq),psivlq (npq),
     &              qprimq (npq),qvalue (npq),psimshq(npq),
     &              psivalq(npq),psinrmq(npq),psisqrq(npq),
     &              psirh2q(npq),psirhoq(npq),psinrhq(npq),
     &              psitorq(npq),psintfq(npq),psisqtq(npq),
     &              psivolq(npq),psinvlq(npq),psisqvq(npq),
     &              dpsir2q(npq),dpsirhq(npq),dpsinrq(npq),
     &              dpsitfq(npq),dpsintq(npq),dpsistq(npq),
     &              dpsivlq(npq),dpsinvq(npq),dpsisvq(npq),
     &              shearps(npq),shearrh(npq),
     &              sheartf(npq),shearvl(npq),
     &              epslrh1 (npq),shearrh1(npq),shearfrh(npq),
     &              epslvl1 (npq),shearvl1(npq),shearfvl(npq)
      common/labels/lbchpv, lbchsp, lbchrh, lbchro,
     &              lbchtf, lbchst, lbchvl, lbchsv,
     &              pclab,  tclab,  labpsi,
     &              labpsiv,labpsin,labsqtp,
     &              labrhov,labrho2,labrhon,
     &              labtorf,labtorn,labsqtt,
     &              labvolm,labvoln,labsqtv,
     &              labdrh2,labdrho,labdtor,
     &              labdsqt,labdvol,labdsqv
      common/pldf/  x0min,x0max,y0min,y0max 
      common/inout/ kutty,kuout,kuinp,kueql,kumap,kugrd,kurdm,kurdn,
     &              kufst,krdeqlb,krdmgta,krdngta,krdfast
      common/flnm/  filout,filinp,fileql,filmap,filgrd,filrdm,filrdn,
     &              filfst,filplt
c
      dimension profle(nlx)
      equivalence (st1 (1),profle(1))
c
c -------------------------------------------------------------
c
c
c
c 1.0 Initialization
c
c 1.1 Return if no plot is required
c
      if(iplotm .lt. 7) return
c
c
c 1.2 Constants
c
      bignum    = bigno
      smallno   = abs(roundff)
c
c
c 1.3 Plot dimensions
c
      x0mn      = 0.120
      x0mx      = 0.820
      y0mn      = 0.200
      y0mx      = 0.900
c
c
c 1.4 Heading and Legend
c
c 1.4.1 Heading position
c
      xxhd      = 0.0600
      yyhd      = 0.9500
c
c 1.4.2 Heading: case, size, and orientation
c
      icashd    = +2
      isizhd    = +3
      iornhd    =  0
c
c 1.4.3 Legend position
c
      xxlgl     = 0.2000
      yylgl     = 0.2550
      xxlgr     = 0.4750
      yylgr     = 0.3000
c
c 1.4.4 Legend: case, size, and orientation
c
      icaslgl   = +2
      isizlgl   = +1
      iornlgl   =  0
c
      icaslgr   = +2
      isizlgr   = +1
      iornlgr   =  0
c
c
c 1.5 Label position parameters
c
c 1.5.1 Position increments
c
      nspx0     = 16
      nspy0     = 30
      nspz0     = 30
c
c 1.5.2 Axis label position fractions
c
      mspps0    = (2*nspx0)/3
      mspdrh0   = (3*nspy0)/4
      mspdro0   = (2*nspy0)/3
      mspdtf0   = (5*nspz0)/6
      mspdst0   = (3*nspz0)/4
      mspdvl0   = (2*nspz0)/3
      mspdsv0   = (1*nspz0)/2
c
c
c 1.6 Default parameters for pointc
c
      numbpc    = 1
      incpcx    = 0
      incpcy    = 0
      delpcx    = 0.0
      delpcy    = 0.0
c
c
c
c 2.0 Axis settings
c
c 2.1 Set the normalized tick parameters
c
c 2.1.1 Tick scale skipping
c
      ntikxs0   = 2
      ntikys0   = 1
      ntikzs0   = 1
c
c 2.1.2 Axis extensions
c
      xtndl     = 1.200
      xtndr     = 1.200
c
c 2.1.3 Normalized tick lengths
c
      xtiky0    = +0.250
      ytikx0    = +0.250
      ztikx0    = -0.250
c
c 2.1.4 Tick scale positions
c
c 2.1.4.1 Label shift along relevant axis
      xsclx0    = +0.300
      yscly0    = -0.250
      zsclz0    = -0.250
c
c 2.1.4.2 Label shift away from relevant axis
      xscly0    = -4.500
      ysclx0    = -1.800
      zsclx0    = +0.020
c
c
c 2.2 Set the axis scales: case, size, and orientation
c
c 2.2.1 Horizontal axis
c
      icshsc    = +1
      iszhsc    = +1
      iorhsc    = +1
c
c 2.2.2 Left Vertical axis
c
      icsvsy    = +1
      iszvsy    = +1
      iorvsy    =  0
c
c 2.2.3 Right Vertical axis
c
      icsvsz    = +1
      iszvsz    = +1
      iorvsz    =  0
c
c
c 2.3 Axis labels specific to quantities plotted
c
c 2.3.1 Horizontal axis
c
c 2.3.1.1 Label and label position
      splby10   = -4.80
      splabel   = pclab
c
c 2.3.1.2 Horizontal axis label: case, size, and orientation
      icshlb    = +1
      iszhlb    = +3
      iorhlb    =  0
c
c 2.3.2 Left  Vertical axis
c
c 2.3.2.1 Label and label positions
c
c 2.3.2.1.1 Label and label and scale position for dpsirh2
      drhlbx10  = -2.900
      drhlabel  = labdrh2
c
c 2.3.2.1.2 Label and label and scale position for dpsirho
      drolbx10  = -2.900
      droscx10  = -0.5000
      droscy10  = -0.7500
      drolabel  = labdrho
c
c 2.3.2.2 Left Vertical axis label and scale characteristics
c
c 2.3.2.2.1 Axis label: case, size, and orientation for dpsirh2
      icslbdrh  = +1
      iszlbdrh  = +2
      iorlbdrh  =  0
c
c 2.3.2.2.2 Case, size, and orientation for dpsirho
c
c 2.3.2.2.2.1 Axis label
      icslbdro  = +1
      iszlbdro  = +2
      iorlbdro  =  0
c
c 2.3.2.2.2.2 Axis scale
      icsscdro  = +1
      iszscdro  = +1
      iorscdro  =  0
c
c 2.3.3 Right Vertical axis
c
c 2.3.3.1 Label and label position 
c
c 2.3.3.1.1 Label and label position for dpsitor
      dtflbx10  = +1.0000
      dtflabel  = labdtor
c
c 2.3.3.1.2 Label and label and scale position for dpsisqt
      dstlbx10  = +1.0000
      dstscx10  = -0.0000
      dstscy10  = -0.7500
      dstlabel  = labdsqt
c
c 2.3.3.1.3 Label and label and scale position for dpsivol
      dvllbx10  = +1.0000
      dvlscx10  = -0.0000
      dvlscy10  = -0.7500
      dvllabel  = labdvol
c
c 2.3.3.1.4 Label and label and scale position for dpsisqv
      dsvlbx10  = +1.0000
      dsvscx10  = -0.0000
      dsvscy10  = -0.7500
      dsvlabel  = labdsqv
c
c 2.3.3.2 Right Vertical axis label and scale characteristics
c
c 2.3.3.2.1 Axis label: case, size, and orientation for dpsitor
      icslbdtf  = +1
      iszlbdtf  = +2
      iorlbdtf  =  0
c
c 2.3.3.2.2 Case, size, and orientation for dpsisqt
c
c 2.3.3.2.2.1 Axis label
      icslbdst  = +1
      iszlbdst  = +2
      iorlbdst  =  0
c
c 2.3.3.2.2.2 Axis scale
      icsscdst  = +1
      iszscdst  = +1
      iorscdst  =  0
c
c 2.3.3.2.3 Case, size, and orientation for dpsivol
c
c 2.3.3.2.3.1 Axis label
      icslbdvl  = +1
      iszlbdvl  = +2
      iorlbdvl  =  0
c
c 2.3.3.2.3.2 Axis scale
      icsscdvl  = +1
      iszscdvl  = +1
      iorscdvl  =  0
c
c 2.3.3.2.4 Case, size, and orientation for dpsisqv
c
c 2.3.3.2.4.1 Axis label
      icslbdsv  = +1
      iszlbdsv  = +2
      iorlbdsv  =  0
c
c 2.3.3.2.4.2 Axis scale
      icsscdsv  = +1
      iszscdsv  = +1
      iorscdsv  =  0
c
c
c 2.4 Line type parameters
c
c 2.4.1 dpsirh2 profile
c
      icschdrh  = +2
      iszchdrh  = +1
      iorchdrh  =  0
      kspchdrh  = +2
      chdrh     = lbchrh
c
c 2.4.2 dpsirho profile
c
      icschdro  = +2
      iszchdro  = +1
      iorchdro  =  0
      kspchdro  = +2
      chdro     = lbchro
c
c 2.4.3 dpsitor profile
c
      icschdtf  = +2
      iszchdtf  = +1
      iorchdtf  =  0
      kspchdtf  = +2
      chdtf     = lbchtf
c
c 2.4.4 dpsisqt profile
c
      icschdst  = +2
      iszchdst  = +1
      iorchdst  =  0
      kspchdst  = +2
      chdst     = lbchst
c
c 2.4.5 dpsivol profile
c
      icschdvl  = +2
      iszchdvl  = +1
      iorchdvl  =  0
      kspchdvl  = +2
      chdvl     = lbchvl
c
c 2.4.6 dpsisqv profile
c
      icschdsv  = +2
      iszchdsv  = +1
      iorchdsv  =  0
      kspchdsv  = +2
      chdsv     = lbchsv
c
c
c
c 3.0 Find the function minima and maxima
c
c 3.1 Search the physically scaled data
c
      spmin      = psivmsh(  1  )
      spmax      = psivmsh(jpsi2)
      drhmin     = +bignum
      drhmax     = -bignum
      dromin     = +bignum
      dromax     = -bignum
      dtfmin     = +bignum
      dtfmax     = -bignum
      dstmin     = +bignum
      dstmax     = -bignum
      dvlmin     = +bignum
      dvlmax     = -bignum
      dsvmin     = +bignum
      dsvmax     = -bignum
      sprange    =  spmax-spmin
      do 10  j   = 1,jpsi2
      jp0        = j
      psvlu      = psivmsh(jp0)
      drhvlu     = dpsirh2(jp0)
      drovlu     = dpsirho(jp0)
      dtfvlu     = dpsitor(jp0)
      dstvlu     = dpsisqt(jp0)
      dvlvlu     = dpsivol(jp0)
      dsvvlu     = dpsisqv(jp0)
c
      drhmin     = amin1(drhmin,drhvlu)
      drhmax     = amax1(drhmax,drhvlu)
      dromin     = amin1(dromin,drovlu)
      dromax     = amax1(dromax,drovlu)
      dtfmin     = amin1(dtfmin,dtfvlu)
      dtfmax     = amax1(dtfmax,dtfvlu)
      dstmin     = amin1(dstmin,dstvlu)
      dstmax     = amax1(dstmax,dstvlu)
      dvlmin     = amin1(dvlmin,dvlvlu)
      dvlmax     = amax1(dvlmax,dvlvlu)
      dsvmin     = amin1(dsvmin,dsvvlu)
      dsvmax     = amax1(dsvmax,dsvvlu)
   10 continue
c
c
c 3.2 Adjust the maxima as needed
c
c 3.2.1 Horizontal axis extrema
c
      spmin0     = spmin
      spmax0     = spmax
      sprang0    = spmax0 - spmin0
c
c 3.2.2 Vertical axis extrema: Include origin in the range
c
      drhmin0    = amin1(drhmin,0.0)
      drhmax0    = amax1(drhmax,0.0)
c
      dromin0    = amin1(dromin,0.0)
      dromax0    = amax1(dromax,0.0)
c
      dtfmin0    = amin1(dtfmin,0.0)
      dtfmax0    = amax1(dtfmax,0.0)
c
      dstmin0    = amin1(dstmin,0.0)
      dstmax0    = amax1(dstmax,0.0)
c
      dvlmin0    = amin1(dvlmin,0.0)
      dvlmax0    = amax1(dvlmax,0.0)
c
      dsvmin0    = amin1(dsvmin,0.0)
      dsvmax0    = amax1(dsvmax,0.0)
c
c
c 3.3 Set the range parameters
c
c 3.3.1 Range for dpsirh2
c
      if(drhmin0 .ge. 0.0  .and.  drhmax0 .gt. 0.0) ldrhvl  = +1
      if(drhmin0 .lt. 0.0  .and.  drhmax0 .gt. 0.0) ldrhvl  =  0
      if(drhmin0 .lt. 0.0  .and.  drhmax0 .le. 0.0) ldrhvl  = -1
c
c 3.3.2 Range for dpsirho
c
      if(dromin0 .ge. 0.0  .and.  dromax0 .gt. 0.0) ldrovl  = +1
      if(dromin0 .lt. 0.0  .and.  dromax0 .gt. 0.0) ldrovl  =  0
      if(dromin0 .lt. 0.0  .and.  dromax0 .le. 0.0) ldrovl  = -1
c
c 3.3.3 Range for dpsitor
c
      if(dtfmin0 .ge. 0.0  .and.  dtfmax0 .gt. 0.0) ldtfvl  = +1
      if(dtfmin0 .lt. 0.0  .and.  dtfmax0 .gt. 0.0) ldtfvl  =  0
      if(dtfmin0 .lt. 0.0  .and.  dtfmax0 .le. 0.0) ldtfvl  = -1
c
c 3.3.4 Range for dpsisqt
c
      if(dstmin0 .ge. 0.0  .and.  dstmax0 .gt. 0.0) ldstvl  = +1
      if(dstmin0 .lt. 0.0  .and.  dstmax0 .gt. 0.0) ldstvl  =  0
      if(dstmin0 .lt. 0.0  .and.  dstmax0 .le. 0.0) ldstvl  = -1
c
c 3.3.5 Range for dpsivol
c
      if(dvlmin0 .ge. 0.0  .and.  dvlmax0 .gt. 0.0) ldvlvl  = +1
      if(dvlmin0 .lt. 0.0  .and.  dvlmax0 .gt. 0.0) ldvlvl  =  0
      if(dvlmin0 .lt. 0.0  .and.  dvlmax0 .le. 0.0) ldvlvl  = -1
c
c 3.3.6 Range for dpsisqv
c
      if(dsvmin0 .ge. 0.0  .and.  dsvmax0 .gt. 0.0) ldsvvl  = +1
      if(dsvmin0 .lt. 0.0  .and.  dsvmax0 .gt. 0.0) ldsvvl  =  0
      if(dsvmin0 .lt. 0.0  .and.  dsvmax0 .le. 0.0) ldsvvl  = -1
c
c 3.3.7 Check the ranges are valid
c
      if(ldrhvl .le. 0) then
        call abortjob
     &        ('plotdrh ',  1,   'dpsi/drho2 scale non positive       '
     &        ,'ldrhvl  ', ldrhvl,   'ldrovl  ', ldrovl,   -1)
      endif
c
      if(ldrovl .le. 0) then
        call abortjob
     &        ('plotdrh ',  2,   'dpsi/drhov scale non positive       '
     &        ,'ldrhvl  ', ldrhvl,   'ldrovl  ', ldrovl,   -1)
      endif
c
      if(ldtfvl .le. 0) then
        call abortjob
     &        ('plotdrh ',  3,   'dpsi/dtorf scale non positive       '
     &        ,'ldtfvl  ', ldtfvl,   'ldstvl  ', ldstvl,   -1)
      endif
c
      if(ldstvl .le. 0) then
        call abortjob
     &        ('plotdrh ',  4,   'dpsi/dsqtf scale non positive       '
     &        ,'ldtfvl  ', ldtfvl,   'ldstvl  ', ldstvl,   -1)
      endif
c
      if(ldvlvl .le. 0) then
        call abortjob
     &        ('plotdrh ',  5,   'dpsi/dvolm scale non positive       '
     &        ,'ldvlvl  ', ldvlvl,   'ldsvvl  ', ldsvvl,   -1)
      endif
c
      if(ldsvvl .le. 0) then
        call abortjob
     &        ('plotdrh ',  6,   'dpsi/dsqvl scale non positive       '
     &        ,'ldvlvl  ', ldvlvl,   'ldsvvl  ', ldsvvl,   -1)
      endif

c
c
c 3.4 Vertical axis scale factors
c
      ratiodro = dromax0/drhmax0
      ratiodtf = dtfmax0/drhmax0
      ratiodst = dstmax0/drhmax0
      ratiodvl = dvlmax0/drhmax0
      ratiodsv = dsvmax0/drhmax0
c
      if(abs(ratiodro) .le. smallno) then
        call abortjob
     &        ('plotdrh ',  7,   'Vanishing scale for dpsirho         '
     &        ,'ldrhvl  ', ldrhvl,   'ldrovl  ', ldrovl,   -1)
        ratiodro   = 1.0
      endif
c
      if(abs(ratiodtf) .le. smallno) then
        call abortjob
     &        ('plotdrh ',  8,   'Vanishing scale for dpsitor         '
     &        ,'ldrhvl  ', ldrhvl,   'ldtfvl  ', ldtfvl,   -1)
        ratiodtf   = 1.0
      endif
c
      if(abs(ratiodst) .le. smallno) then
        call abortjob
     &        ('plotdrh ',  9,   'Vanishing scale for dpsisqt         '
     &        ,'ldrhvl  ', ldrhvl,   'ldstvl  ', ldstvl,   -1)
        ratiodst   = 1.0
      endif
c
      if(abs(ratiodvl) .le. smallno) then
        call abortjob
     &        ('plotdrh ', 10,   'Vanishing scale for dpsivol         '
     &        ,'ldrhvl  ', ldrhvl,   'ldvlvl  ', ldvlvl,   -1)
        ratiodvl   = 1.0
      endif
c
      if(abs(ratiodsv) .le. smallno) then
        call abortjob
     &        ('plotdrh ', 11,   'Vanishing scale for dpsisqv         '
     &        ,'ldrhvl  ', ldrhvl,   'ldsvvl  ', ldsvvl,   -1)
        ratiodsv   = 1.0
      endif
c
      ratrorh  =   1.0   /ratiodro
      ratsttf  = ratiodtf/ratiodst
      ratvltf  = ratiodtf/ratiodvl
      ratsvtf  = ratiodtf/ratiodsv
c
c
c
c 4.0 Set the physical plot parameters
c
c 4.1 Set the map extrema from the extrema
c
c 4.1.1 Left scale extrema
c
      alfmin    =       drhmin0
      alfmax    = xtndl*drhmax0
c
      alfmin    = amin1(alfmin,ratrorh*dromin0)
      alfmax    = amax1(alfmax,ratrorh*dromax0)
c
c 4.1.2 Right scale extrema
c
      artmin    =       dtfmin0
      artmax    = xtndr*dtfmax0
c
      artmin    = amin1(artmin,ratsttf*dstmin0)
      artmin    = amin1(artmin,ratvltf*dvlmin0)
      artmin    = amin1(artmin,ratsvtf*dstmin0)
      artmax    = amax1(artmax,ratsttf*dstmax0)
      artmax    = amax1(artmax,ratvltf*dvlmax0)
      artmax    = amax1(artmax,ratsvtf*dsvmax0)
c
c
c 4.2 Set the number of tick marks
c
c 4.2.1 Horizontal axis
c
      nspx      = 10
c
c 4.2.2 Left  Vertical axis
c
      nspy1     = ifix(alfmax - alfmin)
      nspy1     = min0(nspy1,nspy0)
      nspy1     = max0(nspy1,1)
c
      if    (nspy1 .le. 4) then
        if    (nspy1 .le. 1) then
          nspy    = 5*nspy1
        elseif(nspy1 .le. 2) then
          nspy    = 4*nspy1
        else
          nspy    = 2*nspy1
        endif
c
      elseif(nspy1 .gt. 4) then
        nspy    =   nspy1
      endif
c
c 4.2.3 Right Vertical axis
c
      nspz1     = ifix(artmax - artmin)
      nspz1     = min0(nspz1,nspz0)
      nspz1     = max0(nspz1,1)
c
      if    (nspz1 .le. 4) then
        if    (nspz1 .le. 1) then
          nspz    = 5*nspz1
        elseif(nspz1 .le. 2) then
          nspz    = 4*nspz1
        else
          nspz    = 2*nspz1
        endif
c
      elseif(nspz1 .gt. 4) then
        nspz    =   nspz1
      endif
c
c
c 4.3 Normalize the physical spacing increments to the physical range
c
c 4.3.1 Standard spacing
c
      delx0     = (spmax0 - spmin0)/nspx0
      dely0     = (alfmax - alfmin)/nspy0
      delz0     = (artmax - artmin)/nspz0
c
c 4.3.2 Tick spacing
c
      delx      = (spmax0 - spmin0)/nspx
      dely      = (alfmax - alfmin)/nspy
      delz      = (artmax - artmin)/nspz
c
c
c 4.4 Set the physical spacings
c
c 4.4.1 Tick scale values
c
c 4.4.1.1 Horizontal axis
      ntikxs    = min0(ntikxs0,nspx)
      ntikxs    = max0(ntikxs,  1  )
      if(nspx .lt. 5) ntikxs  = 1
c
c 4.4.1.2 Left Vertical axis
      ntikys    = min0(ntikys0,nspy)
      ntikys    = max0(ntikys,  1  )
      if(nspy .lt. 5) ntikys  = 1
c
c 4.4.1.3 Right Vertical axis
      ntikzs    = min0(ntikzs0,nspz)
      ntikzs    = max0(ntikzs,  1  )
      if(nspz .lt. 5) ntikzs  = 1
c
c 4.4.2 Axis scale positions
c
c 4.4.2.1 Label shift along relevant axis
      xsclx     = xsclx0*delx0
      yscly     = yscly0*dely0
      zsclz     = zsclz0*delz0
c
c 4.4.2.2 Label shift away from relevant axis
      xscly     = xscly0*dely0
      ysclx     = ysclx0*delx0
      zsclx     = zsclx0*delx0
c
c 4.4.3 Tick lengths
c
      xtiky     = xtiky0*dely0
      ytikx     = ytikx0*delx0
      ztikx     = ztikx0*delx0
c
c 4.4.4 Label positions specific to profiles
c
c 4.4.4.1 Horizontal axis
      splby1    = splby10*dely0
c
c 4.4.4.2 Left Vertical Axis
      drhlbx1   = drhlbx10*delx0
      drolbx1   = drolbx10*delx0
      droscx1   = droscx10*delx0
      droscy1   = droscy10*dely0
c
c 4.4.4.3 Right Vertical Axis
      dtflbx1   = dtflbx10*delx0
      dstlbx1   = dstlbx10*delx0
      dstscx1   = dstscx10*delx0
      dstscy1   = dstscy10*delz0
      dvllbx1   = dvllbx10*delx0
      dvlscx1   = dvlscx10*delx0
      dvlscy1   = dvlscy10*delz0
      dsvlbx1   = dsvlbx10*delx0
      dsvscx1   = dsvscx10*delx0
      dsvscy1   = dsvscy10*delz0
c
c
c 4.5 Additional lines
c
      zero      = 0.0
      kspzero   =  2
c
c
c 5.0 Print out the plot heading
c
c 5.1 Set up the plot in normalized coordinates
c
      call map(x0min,x0max,y0min,y0max,x0min,x0max,y0min,y0max)
c
c
c 5.2 Plot the four borders
c
      call line(x0mn,y0mn,x0mn,y0mx)
      call line(x0mn,y0mx,x0mx,y0mx)
      call line(x0mx,y0mx,x0mx,y0mn)
      call line(x0mn,y0mn,x0mx,y0mn)
c
c
c 5.3 Set the position and plot the heading text
c
      call setlch(xxhd,  yyhd,  icashd,isizhd,iornhd,-1)
      write(string,1000) splabel 
      call wrtstr(string,1)
c
c
c 5.4 Set the position and plot the legend
c
c 5.4.1 Left axis
c
      call setlch(xxlgl, yylgl, icaslgl,isizlgl,iornlgl,-1)
      write(string,1100) drhlabel,chdrh
      call wrtstr(string,1)
      write(string,1110) drolabel,chdro
      call wrtstr(string,1)
c
c 5.4.2 Right axis
c
      call setlch(xxlgr, yylgr, icaslgr,isizlgr,iornlgr,-1)
      write(string,1120) dtflabel,chdtf
      call wrtstr(string,1)
      write(string,1130) dstlabel,chdst
      call wrtstr(string,1)
      write(string,1140) dvllabel,chdvl
      call wrtstr(string,1)
      write(string,1150) dsvlabel,chdsv
      call wrtstr(string,1)
c
c
c
c 6.0 Plot dpsirh2, and dpsirho
c
c 6.1 Set the plot page in physical units (left axis)
c
      call map(spmin0,spmax0,alfmin,alfmax,x0mn,x0mx,y0mn,y0mx)
c
c
c 6.2 Set up the Horizontal axis
c
c 6.2.1 Plot the horizontal axis tick marks and scale
c
c 6.2.1.1 Set the base tick parameters
      vrtpos     = alfmin
      vrtpos1    = vrtpos
      vrtpos2    = vrtpos + xtiky
      vrtpossc   = vrtpos + xscly
      do 100 k   = 0,nspx
      kskpp      = (k/ntikxs)*ntikxs
      spos       = spmin0 + k*delx
      spossc     = spos   + xsclx
c
c 6.2.1.2 Plot the tick mark
      call line  (spos, vrtpos1, spos, vrtpos2 )
c
c 6.2.1.3 Plot the scale value
      if(k .eq. kskpp) then
        call setlch(spossc,vrtpossc,icshsc,iszhsc,iorhsc,-1)
        write(string,2000) spos
        call wrtstr(string,1)
      endif
  100 continue
c
c 6.2.2 Plot the horizontal axis label
c
      sposlb    = spmin0  + mspps0*delx0
      vrtposlb  = vrtpos  + splby1
      call setlch(sposlb,vrtposlb,icshlb,iszhlb,iorhlb,-1)
      write(string,2010) splabel
      call wrtstr(string,1)
c
c
c 6.3 Set up the Left  vertical axis
c
c 6.3.1 Plot the  vertical  axis tick marks and scale
c
c 6.3.1.1 Set the base tick parameters
      spos       = spmin0
      spos1      = spos
      spos2      = spos + ytikx
      spossc     = spos + ysclx
      do 110 i   = 0,nspy
      iskpp      = (i/ntikys)*ntikys
      alfpos     = alfmin  + i*dely
      alfpossc   = alfpos  + yscly
c
c 6.3.1.2 Plot the tick mark
      call line  (spos1, alfpos,  spos2, alfpos  )
c
c 6.3.1.3 Plot the scale value
      if(i .eq. iskpp) then
        call setlch(spossc,alfpossc,icsvsy,iszvsy,iorvsy,-1)
        write(string,2020) alfpos
        call wrtstr(string,1)
      endif
  110 continue
c
c
c 6.4 Plot the dpsirh2 profile
c
c 6.4.1 Plot the  vertical  axis label
c
      sposlb   = spmin0  + drhlbx1
      drhposlb = alfmin  + mspdrh0*dely0
      call setlch(sposlb,drhposlb,icslbdrh,iszlbdrh,iorlbdrh,-1)
      write(string,2100) drhlabel
      call wrtstr(string,1)
c
c 6.4.2 Reset the printing parameters for pointc
c
      call setpch(icschdrh,iszchdrh,iorchdrh,kspchdrh)
c
c 6.4.3 Plot the renormalized values
c
c 6.4.3.1 Plot the points and store in profle
      do 200 j   = 1,jpsi2
      jp         = j
      sppt       =         psivmsh(jp)
      drhpt      =         dpsirh2(jp)
      profle(jp) =         drhpt
      call pointc(chdrh,sppt,drhpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  200 continue
c
c 6.4.3.2 Plot the curve
      call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
c
c
c 6.5 Plot the dpsirho profile
c
c 6.5.1 Plot the  vertical  axis label
c
c 6.5.1.1 Write the label itself
      sposlb   = spmin0  + drolbx1
      droposlb = alfmin  + mspdro0*dely0
      call setlch(sposlb,droposlb,icslbdro,iszlbdro,iorlbdro,-1)
      write(string,2200) drolabel
      call wrtstr(string,1)
c
c 6.5.1.2 Plot the scale factor
      sposscu  = sposlb    + droscx1
      dropossc = droposlb  + droscy1
      call setlch(sposscu,dropossc,icsscdro,iszscdro,iorscdro,-1)
c
      if(ratrorh .gt. 1.0) write(string,2210) ratrorh
      if(ratrorh .lt. 1.0) write(string,2220) ratrorh
      if(ratrorh .ne. 1.0) call wrtstr(string,1)
c
c 6.5.2 Reset the printing parameters for pointc
c
      call setpch(icschdro,iszchdro,iorchdro,kspchdro)
c
c 6.5.3 Plot the renormalized values
c
c 6.5.3.1 Plot the points and store in profle
      do 210 j   = 1,jpsi2
      jp         = j
      sppt       =         psivmsh(jp)
      dropt      = ratrorh*dpsirho(jp)
      profle(jp) =         dropt
      call pointc(chdro,sppt,dropt,numbpc,incpcx,incpcy,delpcx,delpcy)
  210 continue
c
c 6.5.3.2 Plot the curve
      call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
c
c
c 6.6 Plot a line if the origin is included in the range
c
      if(ldrhvl .eq. 0  .or.  ldrovl .eq. 0) then
        call linep(spmin0,zero,spmax0,zero,kspzero)
      endif
c
c
c
c 7.0 Plot dpsitor, dpsisqt, dpsivol, and dpsisqv
c
c 7.1 Renormalize the map for the derivative
c
      call map(spmin0,spmax0,artmin,artmax,x0mn,x0mx,y0mn,y0mx)
c
c
c 7.2 Set up the Right Vertical axis
c
c 7.2.1 Plot the  vertical  axis tick marks and scale
c
c 7.2.1.1 Set the base tick parameters
      spos       = spmax0
      spos1      = spos + ztikx
      spos2      = spos
      spossc     = spos + zsclx
      do 300 i   = 0,nspz
      iskpp      = (i/ntikzs)*ntikzs
      artpos     = artmin  + i*delz
      artpossc   = artpos  + zsclz
c
c 7.2.1.2 Plot the tick mark
      call line  (spos1, artpos,  spos2, artpos  )
c
c 7.2.1.3 Plot the scale value
      if(i .eq. iskpp) then
        call setlch(spossc,artpossc,icsvsz,iszvsz,iorvsz,-1)
        write(string,2500) artpos
        call wrtstr(string,1)
      endif
  300 continue
c
c
c 7.3 Plot the profile for dpsitor
c
c 7.3.1 Plot the  vertical  axis labels for dpsitor
      sposlb    = spmax0  + dtflbx1
      dtfposlb  = artmin  + mspdtf0*delz0
      call setlch(sposlb,dtfposlb,icslbdtf,iszlbdtf,iorlbdtf,-1)
      write(string,2600) dtflabel
      call wrtstr(string,1)
c
c 7.3.2 Reset the printing parameters for pointc
c
      call setpch(icschdtf,iszchdtf,iorchdtf,kspchdtf)
c
c 7.3.3 Plot the profile
c
c 7.3.3.1 Plot the points and store in profle
      do 400 j   = 1,jpsi2
      jp         = j
      sppt       =         psivmsh(jp)
      dtfpt      =         dpsitor(jp)
      profle(jp) =         dtfpt
      call pointc(chdtf,sppt,dtfpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  400 continue
c
c 7.3.3.2 Plot the curve
      call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
c
c
c 7.4 Plot the profile for dpsisqt
c
c 7.4.1 Plot the  vertical  axis labels for dpsisqt
c
c 7.4.1.1 Write the label itself
      sposlb    = spmax0  + dstlbx1
      dstposlb  = artmin  + mspdst0*delz0
      call setlch(sposlb,dstposlb,icslbdst,iszlbdst,iorlbdst,-1)
      write(string,2700) dstlabel
      call wrtstr(string,1)
c
c 7.4.1.2 Plot the scale factor
      sposscw  = sposlb    + dstscx1
      dstpossc = dstposlb  + dstscy1
      call setlch(sposscw,dstpossc,icsscdst,iszscdst,iorscdst,-1)
c
      if(ratsttf .gt. 1.0) write(string,2710) ratsttf
      if(ratsttf .lt. 1.0) write(string,2720) ratsttf
      if(ratsttf .ne. 1.0) call wrtstr(string,1)
c
c 7.4.2 Reset the printing parameters for pointc
c
      call setpch(icschdst,iszchdst,iorchdst,kspchdst)
c
c 7.4.3 Plot the profile
c
c 7.4.3.1 Plot the points and store in profle
      do 410 j   = 1,jpsi2
      jp         = j
      sppt       =         psivmsh(jp)
      dstpt      = ratsttf*dpsisqt(jp)
      profle(jp) =         dstpt
      call pointc(chdst,sppt,dstpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  410 continue
c
c 7.4.3.2 Plot the curve
      call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
c
c
c 7.5 Plot the profile for dpsivol
c
c 7.5.1 Plot the  vertical  axis labels for dpsivol
c
c 7.5.1.1 Write the label itself
      sposlb    = spmax0  + dvllbx1
      dvlposlb  = artmin  + mspdvl0*delz0
      call setlch(sposlb,dvlposlb,icslbdvl,iszlbdvl,iorlbdvl,-1)
      write(string,2800) dvllabel
      call wrtstr(string,1)
c
c 7.5.1.2 Plot the scale factor
      sposscv  = sposlb    + dvlscx1
      dvlpossc = dvlposlb  + dvlscy1
      call setlch(sposscv,dvlpossc,icsscdvl,iszscdvl,iorscdvl,-1)
c
      if(ratvltf .gt. 1.0) write(string,2810) ratvltf
      if(ratvltf .lt. 1.0) write(string,2820) ratvltf
      if(ratvltf .ne. 1.0) call wrtstr(string,1)
c
c 7.5.2 Reset the printing parameters for pointc
c
      call setpch(icschdvl,iszchdvl,iorchdvl,kspchdvl)
c
c 7.5.3 Plot the profile
c
c 7.5.3.1 Plot the points and store in profle
      do 420 j   = 1,jpsi2
      jp         = j
      sppt       =         psivmsh(jp)
      dvlpt      = ratvltf*dpsivol(jp)
      profle(jp) =         dvlpt
      call pointc(chdvl,sppt,dvlpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  420 continue
c
c 7.5.3.2 Plot the curve
      call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
c
c
c 7.6 Plot the profile for dpsisqv
c
c 7.6.1 Plot the  vertical  axis labels for dpsisqv
c
c 7.6.1.1 Write the label itself
      sposlb    = spmax0  + dsvlbx1
      dsvposlb  = artmin  + mspdsv0*delz0
      call setlch(sposlb,dsvposlb,icslbdsv,iszlbdsv,iorlbdsv,-1)
      write(string,2900) dsvlabel
      call wrtstr(string,1)
c
c 7.6.1.2 Plot the scale factor
      sposscx  = sposlb    + dsvscx1
      dsvpossc = dsvposlb  + dsvscy1
      call setlch(sposscx,dsvpossc,icsscdsv,iszscdsv,iorscdsv,-1)
c
      if(ratsvtf .gt. 1.0) write(string,2910) ratsvtf
      if(ratsvtf .lt. 1.0) write(string,2920) ratsvtf
      if(ratsvtf .ne. 1.0) call wrtstr(string,1)
c
c 7.6.2 Reset the printing parameters for pointc
c
      call setpch(icschdsv,iszchdsv,iorchdsv,kspchdsv)
c
c 7.6.3 Plot the profile
c
c 7.6.3.1 Plot the points and store in profle
      do 430 j   = 1,jpsi2
      jp         = j
      sppt       =         psivmsh(jp)
      dsvpt      = ratsvtf*dpsisqv(jp)
      profle(jp) =         dsvpt
      call pointc(chdsv,sppt,dsvpt,numbpc,incpcx,incpcy,delpcx,delpcy)
  430 continue
c
c 7.6.3.2 Plot the curve
      call trace(psivmsh,profle,jpsi2,-1,-1,0.0,0.0)
c
c
c 7.7 Plot a line if the origin is included in the range
c
      if(ldtfvl .eq. 0  .or.  ldstvl .eq. 0  .or.
     &   ldvlvl .eq. 0  .or.  ldsvvl .eq. 0) then
        call linep(spmin0,zero,spmax0,zero,kspzero)
      endif
c
c
c
c 8.0 Close frame, return, and end
c
      call frame(0)
      return
 1000 format(1x,'Profiles of psi derivatives vs',1x,a8)
 1100 format('dpsirh2=',1x,a8,':',2x,a1)
 1110 format('dpsirho=',1x,a8,':',2x,a1)
 1120 format('dpsitor=',1x,a8,':',2x,a1)
 1130 format('dpsisqt=',1x,a8,':',2x,a1)
 1140 format('dpsivol=',1x,a8,':',2x,a1)
 1150 format('dpsisqv=',1x,a8,':',2x,a1)
 2000 format(f8.4)
 2010 format(a8)
 2020 format(f6.2)
 2100 format(a8)
 2200 format(a8)
 2210 format('(x',f6.2,')')
 2220 format('(x',f6.4,')')
 2500 format(f6.2)
 2600 format(a8)
 2700 format(a8)
 2710 format('(x',f8.2,')')
 2720 format('(x',f6.4,')')
 2800 format(a8)
 2810 format('(x',f8.2,')')
 2820 format('(x',f6.4,')')
 2900 format(a8)
 2910 format('(x',f8.2,')')
 2920 format('(x',f6.4,')')
      end
      subroutine plottab
c
c ------------------------
c plot the profile tables
c ------------------------
c
      parameter (npx=200,ncx=2*npx,np1=npx+1,np2=npx+2,nc1=ncx+1)
      parameter (nxx=513,nxz=nxx)
      parameter (npp=929,ntt=2*npp-1,npb=npp)
      parameter (nvn=7)
c
      character*136 string
c
      common/cnst/  pi,hlfpi,twopi,amu,boltzm,echarg,pmass,roundff,bigno
      common/inpd/  ntor,ncase,idnsty,ndnxp0,ndnxp1,ndnxp2,ifastp,
     &              nmap,neqtyp,mpreset,ndoublt,ndivert,nmtype,
     &              jpsi,itht,jpsi1,jpsi2,ithtp,jsep,
     &              isym,igrid,nham1,nham2,nham3,npowr,
     &              bfieldf,rdefolt,qxin,btdes,qsurf,
     &              gamma,gamav,rmantl,betaf,zeffect,
     &              ncorr,corrfac,nccellr,peqpk0,peqpk1,peqpk2
      common/mapd/  mapmaxd,dpsisl,dpsisd,nqaxis,nwtmag,nfitmax,nfitpts,
     &              ifitrng,jfitrng,jfitchk,fitchek,cnvmag,epsaxs,
     &              maxerlp,maxerlv,delbox,delboz,delac,delav,
     &              delstsf,delstlp,delstlv,prfrac,
     &              nerstop,nerprnt,qptol,tolspln,tolbchi,tolbtor,
     &              tolsymm,tolaugm,errsep,precisn,plossmx,
     &              narcmx,ntrymx,ntdecr,ntmmin,npfit,npcmin,
     &              kuttaop,nrkmax0,nrkmax1,numstp,nwtfitp,nwtfitm,
     &              stepfac,flxstp,psispl,tolnwtp,tolnwtm,
     &              delpakf,delpakc,delpkf,delpkc,psichek,boxtnd,
     &              maptrace,norient,maxcutc,
     &              dresolv,dlclose,pntshft,endtol,
     &              narcln,nangax,nanglm,nbpmax,nwtmax,nslmax,nhfmax,
     &              bperor,sersnm,sertnm,arcmin,delgap,stepcut,
     &              novap,ioutm,iplotm
      common/prof/  nprofl,npst,nfast,nspbc0,nspbc1,
     &              psimsh(npp),sf(npp),sp(npp),
     &              sffp(npp),spp(npp),sfp(npp),sdns(npp),
     &              spfst(npp),spsif(npb),sfast(npb),bcb(4),
     &              csf(3,npp),csp(3,npp),csffp(3,npp),cspp(3,npp),
     &              csfp(3,npp),csdn(3,npp),csfst(3,npp)
      common/dblt/  psisep,xsep,zsep,smap(2),smax(2),smaz(2),numax,
     &              psimx(2),xax(2),zax(2),fitax(14,2),fitsep(14)
      common/ort1/  ndim1,ndim2,xma,zma,rcnt,btor,totcur,
     &              psimax0,psilim0,delpsir0,psimax1,psilim1,delpsir1,
     &              psimax,psilim,delpsir,dpsidec,dpsisv,dpedg1,dpedg2,
     &              eaxe,taxe,qaxe,fqiaxe,fqmaxe,paxe,faxe,ppaxe,ffpaxe,
     &              dnaxe,pfaxe,qppaxe,dlnaxe,dpnaxe,dhnaxe,
     &              arcaxe,dnnorm,pfnorm,
     &              qlim,fqilim,fqmlim,plim,flim,pplim,ffplim,
     &              dnlim,pflim,xlim,zlim,
     &              qmer,btmer,gkappa,arcprev,chicmax,
     &              charcmx,chpstmx,chhammx,arcnrm0,pstnrm0,hamnrm0,
     &              elongax1,elongax2,elongaxs,elongaxp,elongtn(npp)
      common/ort2/  psival(np1),pp(np1),ffp(np1),p(np1),f(np1),
     &              dnsty(np1),pfast(np1),qp(np1),fqpi(np1),fqim(np1),
     &              qdel(np1),psinitl(np1),qpinitl(np1),
     &              qprime(np1),qpint(np1),fqint(np1),
     &              psivl0(np2),rst(np2),zst(np2),rstp(np2),zstp(np2),
     &              rsinitl(np2),zsinitl(np2),
     &              arcsurf(nc1),pestchi(nc1),hamlchi(nc1),
     &              alfarc (nc1),alfpst (nc1),alfham (nc1),
     &              rs(nc1,np1),zs(nc1,np1)
      common/ort3/  rsrf(nc1),zsrf(nc1),chic(nc1),chie(nc1),dydx(nc1),
     &              xjsc(nc1)
      common/ort4/  rcell  (np1,nc1),zcell  (np1,nc1),
     &              dpsdr  (np1,nc1),dpsdz  (np1,nc1),
     &              chiarcl(np1,nc1),chipest(np1,nc1),chihaml(np1,nc1),
     &              xjcarcl(np1,nc1),xjcpest(np1,nc1),xjchaml(np1,nc1),
     &              alfarcl(np1,nc1),alfpest(np1,nc1),alfhaml(np1,nc1),
     &              arcnorm(np1),    pstnorm(np1),    hamnorm(np1),
     &              dlnlval(np1),    dpnlval(np1),    dhnlval(np1)
      common/ort5/  f3 (np1,nc1),f4 (np1,nc1),f5 (np1,nc1),
     &              f7 (np1,nc1),f8 (np1,nc1),f9 (np1,nc1),
     &              f10(np1,nc1),f11(np1,nc1),f12(np1,nc1),
     &              f13(np1,nc1),f14(np1,nc1),f15(np1,nc1),
     &              f16(np1,nc1),f17(np1,nc1),f18(np1,nc1),
     &              f19(np1,nc1),f20(np1,nc1),f21(np1,nc1),
     &              f22(np1,nc1),f23(np1,nc1),f24(np1,nc1),
     &              f25(np1,nc1),f26(np1,nc1),f27(np1,nc1),
     &              f28(np1,nc1)
      common/ort6/  savgax,tavgax,sntax0,sntax1,sntax2,sntax3,sntax4,
     &              svint(nvn),savge(np1),tavge(np1),sint0(np1),
     &              sint1(np1),sint2(np1),sint3(np1),sint4(np1)
      common/vcal/  btnew,btave,bpave,betat,betap,betax0,betax1,volme,
     &              vhalf,bavet(3,nxx),bavep(3,nxx),pvolm(3,nxx),
     &              betav(3,nxx)
      common/pldf/  x0min,x0max,y0min,y0max
      common/inout/ kutty,kuout,kuinp,kueql,kumap,kugrd,kurdm,kurdn,
     &              kufst,krdeqlb,krdmgta,krdngta,krdfast
c
c
c
c 1.0 Print out the profiles
c
      if(iplotm .lt. 7) return
      smn      = 0.00
      smx      = 1.00
      smn0     = 0.00
      smx0     = 1.00
      call map(smn,smx,smn,smx,smn0,smx0,smn0,smx0)
c
      xxt1     = 0.0035
      yyt1     = 0.985
      call setlch(xxt1,yyt1,2,2,0,-1)
      call gtext('  j   psi      q       n',24,-1)
c
      xxt2     = 0.500
      yyt2     = yyt1
      call setlch(xxt2,yyt2,2,2,0,-1)
      call gtext(' p     pprime     f      ffprime',32,-1)
c
c
      xxpr     = 0.0010
      yypr     = 0.935
      call setlch(xxpr,yypr,2,1,0,-1)
c
      jw1      = 0
      psilimw1 = psilim
      qlimw1   = qlim
      dnlimw1  = dnlim
      plimw1   = plim
      pplimw1  = pplim
      flimw1   = flim
      ffplimw1 = ffplim
      write(string,1000) jw1,psilimw1,qlimw1,dnlimw1,plimw1,pplimw1
     &                      ,flimw1,ffplimw1
      call wrtstr(string,1)
c
      do 10 jj = 1,jpsi
      write(string,1000) jj,psival(jj),qp(jj),dnsty(jj),p(jj),pp(jj)
     &                     ,f(jj),ffp(jj)
      call wrtstr(string,1)
   10 continue
c
      jw0      = jpsi1
      psimaxw0 = psimax
      qaxew0   = qaxe
      dnaxew0  = dnaxe
      paxew0   = paxe
      ppaxew0  = ppaxe
      faxew0   = faxe
      ffpaxew0 = ffpaxe
      write(string,1000) jw0,psimaxw0,qaxew0,dnaxew0,paxew0,ppaxew0
     &                      ,faxew0,ffpaxew0
      call wrtstr(string,1)
c
c
c
c 2.0 Call frame, return, and end
c
      call frame(0)
      return
c
 1000 format(i3,1x,e11.4,1x,f9.5,5(1x,e11.4))
      end
      subroutine plotfpj(ktype)
c
c --------------------------
c plot equilibrium profiles
c --------------------------
c
      parameter (npx=200,ncx=2*npx,np1=npx+1,np2=npx+2,nc1=ncx+1)
      parameter (nxx=513,nxz=nxx)
      parameter (npp=929,ntt=2*npp-1,npb=npp)
      parameter (nw2=2*nxx,nh2=2*nxz,nvk0=nw2+nh2)
      parameter (nlx=1441)
c
c

      common/cnst/  pi,hlfpi,twopi,amu,boltzm,echarg,pmass,roundff,bigno
      common/inpd/  ntor,ncase,idnsty,ndnxp0,ndnxp1,ndnxp2,ifastp,
     &              nmap,neqtyp,mpreset,ndoublt,ndivert,nmtype,
     &              jpsi,itht,jpsi1,jpsi2,ithtp,jsep,
     &              isym,igrid,nham1,nham2,nham3,npowr,
     &              bfieldf,rdefolt,qxin,btdes,qsurf,
     &              gamma,gamav,rmantl,betaf,zeffect,
     &              ncorr,corrfac,nccellr,peqpk0,peqpk1,peqpk2
      common/mapd/  mapmaxd,dpsisl,dpsisd,nqaxis,nwtmag,nfitmax,nfitpts,
     &              ifitrng,jfitrng,jfitchk,fitchek,cnvmag,epsaxs,
     &              maxerlp,maxerlv,delbox,delboz,delac,delav,
     &              delstsf,delstlp,delstlv,prfrac,
     &              nerstop,nerprnt,qptol,tolspln,tolbchi,tolbtor,
     &              tolsymm,tolaugm,errsep,precisn,plossmx,
     &              narcmx,ntrymx,ntdecr,ntmmin,npfit,npcmin,
     &              kuttaop,nrkmax0,nrkmax1,numstp,nwtfitp,nwtfitm,
     &              stepfac,flxstp,psispl,tolnwtp,tolnwtm,
     &              delpakf,delpakc,delpkf,delpkc,psichek,boxtnd,
     &              maptrace,norient,maxcutc,
     &              dresolv,dlclose,pntshft,endtol,
     &              narcln,nangax,nanglm,nbpmax,nwtmax,nslmax,nhfmax,
     &              bperor,sersnm,sertnm,arcmin,delgap,stepcut,
     &              novap,ioutm,iplotm
      common/size/  xdim,zdim,redge,zlowr,ipestg
      common/prof/  nprofl,npst,nfast,nspbc0,nspbc1,
     &              psimsh(npp),sf(npp),sp(npp),
     &              sffp(npp),spp(npp),sfp(npp),sdns(npp),
     &              spfst(npp),spsif(npb),sfast(npb),bcb(4),
     &              csf(3,npp),csp(3,npp),csffp(3,npp),cspp(3,npp),
     &              csfp(3,npp),csdn(3,npp),csfst(3,npp)
      common/eqd1/  nx,nz,nxd,nzd,dmx,dmz,
     &              x(nxx),z(nxz),psarray(nxx,nxz),
     &              gpx(nxx,nxz),gpy(nxx,nxz),grsq(nxx,nxz),
     &              cspl1(2,nxx,nh2),work0(nxz,2,nxx),vork0(nvk0,2)
      common/sarc/  ntmax,ntmsh,npcf,npcb,npc,xp(nlx),zp(nlx),gsq(nlx),
     &              arc(nlx),tp(nlx),arcc(nlx),tpp(nlx),bcd(4),
     &              csx(3,nlx),csz(3,nlx),cseq1(3,nlx),cseq2(3,nlx),
     &              st1(nlx),st2(nlx),st3(nlx),csveq(3,nlx),
     &              sv0(nlx),sv1(nlx),sv2(nlx),sv3(nlx),sv4(nlx),
     &              sv5(nlx)
      common/dblt/  psisep,xsep,zsep,smap(2),smax(2),smaz(2),numax,
     &              psimx(2),xax(2),zax(2),fitax(14,2),fitsep(14)
      common/ort1/  ndim1,ndim2,xma,zma,rcnt,btor,totcur,
     &              psimax0,psilim0,delpsir0,psimax1,psilim1,delpsir1,
     &              psimax,psilim,delpsir,dpsidec,dpsisv,dpedg1,dpedg2,
     &              eaxe,taxe,qaxe,fqiaxe,fqmaxe,paxe,faxe,ppaxe,ffpaxe,
     &              dnaxe,pfaxe,qppaxe,dlnaxe,dpnaxe,dhnaxe,
     &              arcaxe,dnnorm,pfnorm,
     &              qlim,fqilim,fqmlim,plim,flim,pplim,ffplim,
     &              dnlim,pflim,xlim,zlim,
     &              qmer,btmer,gkappa,arcprev,chicmax,
     &              charcmx,chpstmx,chhammx,arcnrm0,pstnrm0,hamnrm0,
     &              elongax1,elongax2,elongaxs,elongaxp,elongtn(npp)
      common/ort2/  psival(np1),pp(np1),ffp(np1),p(np1),f(np1),
     &              dnsty(np1),pfast(np1),qp(np1),fqpi(np1),fqim(np1),
     &              qdel(np1),psinitl(np1),qpinitl(np1),
     &              qprime(np1),qpint(np1),fqint(np1),
     &              psivl0(np2),rst(np2),zst(np2),rstp(np2),zstp(np2),
     &              rsinitl(np2),zsinitl(np2),
     &              arcsurf(nc1),pestchi(nc1),hamlchi(nc1),
     &              alfarc (nc1),alfpst (nc1),alfham (nc1),
     &              rs(nc1,np1),zs(nc1,np1)
      common/pldf/  x0min,x0max,y0min,y0max
      common/inout/ kutty,kuout,kuinp,kueql,kumap,kugrd,kurdm,kurdn,
     &              kufst,krdeqlb,krdmgta,krdngta,krdfast
c
c
c
c
c 1.0 Initialization
c
      if(iplotm .lt. 8) return
c
      ntiks   = 10
      ksp     = 2
      yz0     = 0.00
c
      rndff   = roundff
c
c
c
c 2.0 Plot of ffprime and pprime vs. r
c
      if(ktype .eq. 0) then
c
c 2.1 Set up plots
c
        xmap0   = 0.275
        xmap1   = 0.945
        ymap0   = 0.200
        ymap1   = 0.870
c
        ratio   = xdim/zdim
        xl      = redge
        xgap0   = 0.0*xdim
        xgap1   = 0.0*zdim
        if    (ratio .gt. 1.0) then
          xr      =  xl + xdim
          xmov    = -xgap0
        elseif(ratio .le. 1.0) then
          xr      =  xl    + zdim
          xmov    =  xdim  - zdim - xgap1
        endif
c
        xrp     = xr + xmov
        dyy     = zdim/(nz-1)
        nzax    = ifix(zax(1)/dyy)
        jax     = nzax + (nz/2+1)
c
c
c 2.2 ffprime vs. r
c
        do 15 i = 1,nx
        ivl     = i
        wff     = psarray(ivl,jax)
        if    (wff .gt. psilim) then
          st1(ivl) = 0.0
        elseif(wff .le. psilim) then
          if    (npst .eq. 0) then
            ffp1     = sterpl(wff,psimsh,sffp,csffp,nprofl,npp,0,
     &                            rndff, ierf)
            if(ierf .ne. 0) call abortjob
     &        ('plotfpj ',  1,   'Spline evaluation error:ffp(psi)    '
     &        ,'ierf    ', ierf,     'ivl     ', ivl,      -1)
          elseif(npst .eq. 1) then
            f1       = sterpl(wff,psimsh,sf,  csf,  nprofl,npp,0,
     &                            rndff, ierf)
            fp1      = sterpp(wff,psimsh,sf,  csf,  nprofl,npp,0,
     &                            rndff, ierfp)
            if(ierf .ne. 0) call abortjob
     &        ('plotfpj ',  2,   'Interpolation error for f(i,j)      '
     &        ,'ierf    ', ierf,     'ivl     ', ivl,      -1)
            if(ierfp .ne. 0) call abortjob
     &        ('plotfpj ',  3,   'Differentiation error for f(i,j)    '
     &        ,'ierfp   ', ierfp,    'ivl     ', ivl,      -1)
            ffp1     = f1*fp1
          endif
          st1(ivl) = ffp1
        endif
   15   continue
c
        sffmx   = st1(1)
        sffmn   = st1(1)
        do 20 i = 1,nx
        sffmx   = amax1(sffmx,st1(i))
        sffmn   = amin1(sffmn,st1(i))
   20   continue
c
        sffmn   = amin1(0.0,sffmn)
        sffmx   = amax1(0.0,sffmx)
c
        if(sffmn .eq. sffmx) then
          if    (sffmn .eq. 0.0) then
            sffmn  = -1.0
            sffmx  = +1.0
          elseif(sffmn .gt. 0.0) then
            sffmn  =  0.0
            sffmx  = +1.5*sffmx
          elseif(sffmn .lt. 0.0) then
            sffmn  =  1.5*sffmn
            sffmx  =  0.0
          endif
        endif
c
        yl      = sffmn
        yu      = sffmx
        yll     = 2*yl - yu
        call map(xl,xr,yll,yu,xmap0,xmap1,ymap0,ymap1)
c
        call gridm(xl,yl,xrp,yu,'r-axis',  ntiks,-1)
        call gridm(xl,yl,xr, yu,'ff-prime',ntiks,+2)
        call line(xl, yu,xrp,yu)
        call line(xrp,yl,xrp,yu)
        if((yl-yz0)*(yu-yz0) .lt. 0.0) call linep(xl,yz0,xrp,yz0,ksp)
c
        call trace(x,st1,nx,1,1,0.0,0.0)
c
c
c 2.3 pprime vs. r
c
        do 30 i = 1,nx
        ivl     = i
        wpp     = psarray(ivl,jax)
        if    (wpp .gt. psilim) then
          st1(ivl) = 0.0
        elseif(wpp .le. psilim) then
          if    (npst .eq. 0) then
            pp1      = sterpl(wpp,psimsh,spp, cspp,  nprofl,npp,0,
     &                            rndff, ierp)
            if(ierp  .ne. 0) call abortjob
     &        ('plotfpj ',  4,   'Spline evaluation error: pp(psi)    '
     &        ,'ierp    ', ierp,     'ivl     ', ivl,      -1)
          elseif(npst .eq. 1) then
            pp1      = sterpp(wpp,psimsh,sp,  csp,   nprofl,npp,0,
     &                            rndff, ierpp)
            if(ierpp .ne. 0) call abortjob
     &        ('plotfpj ',  5,   'Differentiation error for p(i,j)    '
     &        ,'ierpp   ', ierpp,    'ivl     ', ivl,      -1)
          endif
          st1(ivl) = pp1
        endif
   30   continue
c
        sppmx   = st1(1)
        sppmn   = st1(1)
        do 35 i = 1,nx
        sppmx   = amax1(sppmx,st1(i))
        sppmn   = amin1(sppmn,st1(i))
   35   continue
c
        sppmn   = amin1(0.0,sppmn)
        sppmx   = amax1(0.0,sppmx)
c
        if(sppmn .eq. sppmx) then
          if    (sppmn .eq. 0.0) then
            sppmn  = -1.0
            sppmx  = +1.0
          elseif(sppmn .gt. 0.0) then
            sppmn  =  0.0
            sppmx  = +1.5*sppmx
          elseif(sppmn .lt. 0.0) then
            sppmn  =  1.5*sppmn
            sppmx  =  0.0
          endif
        endif
c
        yu      = sppmx
        yl      = sppmn
        yuu     = 2*yu-yl
        call map(xl,xr,yl,yuu,xmap0,xmap1,ymap0,ymap1)
c
        call gridm(xl,yl,xrp,yu,'r-axis', ntiks,+1)
        call gridm(xl,yl,xr, yu,'p-prime',ntiks,+2)
        call line(xl,yu, xrp,yu)
        call line(xrp,yl,xrp,yu)
        if((yl-yz0)*(yu-yz0) .lt. 0.0) call linep(xl,yz0,xrp,yz0,ksp)
c
        call trace(x,st1,nx,1,1,0.0,0.0)
        call frame(0)
c
c
c
c 3.0 Plot of pressure and current vs. r
c
c 3.1 Set up plots
c
        xmap0   = 0.275
        xmap1   = 0.945
        ymap0   = 0.200
        ymap1   = 0.870
c
        ratio   = xdim/zdim
        xl      = redge
        xgap0   = 0.0*xdim
        xgap1   = 0.0*zdim
        if    (ratio .gt. 1.0) then
          xr      =  xl + xdim
          xmov    = -xgap0
        elseif(ratio .le. 1.0) then
          xr      =  xl    + zdim
          xmov    =  xdim  - zdim - xgap1
        endif
c
        xrp     = xr + xmov
        dyy     = zdim/(nz-1)
        nzax    = ifix(zax(1)/dyy)
        jax     = nzax + (nz/2+1)
c
c
c 3.2 pressure vs. r
c
        do 40 i = 1,nx
        ivl     = i
        wpr     = psarray(ivl,jax)
        if    (wpr .gt. psilim) then
          st1(ivl) = 0.0
        elseif(wpr .le. psilim) then
          st1(ivl) = sterpl(wpr,psimsh,sp,  csp,  nprofl,npp,0,
     &                          rndff, ierp)
          if(ierp .ne. 0) call abortjob
     &        ('plotfpj ',  6,   'Spline evaluation error:  p(psi)    '
     &        ,'ierp    ', ierp,     'ivl     ', ivl,      -1)
        endif
   40   continue
c
        sprmx   = st1(1)
        sprmn   = st1(1)
        do 45 i = 1,nx
        sprmx   = amax1(sprmx,st1(i))
        sprmn   = amin1(sprmn,st1(i))
   45   continue
c
        sprmn   = amin1(0.0,sprmn)
        sprmx   = amax1(0.0,sprmx)
c
        if(sprmn .eq. sprmx) then
          if    (sprmn .eq. 0.0) then
            sprmn  = -1.0
            sprmx  = +1.0
          elseif(sprmn .gt. 0.0) then
            sprmn  =  0.0
            sprmx  = +1.5*sprmx
          elseif(sprmn .lt. 0.0) then
            sprmn  =  1.5*sprmn
            sprmx  =  0.0
          endif
        endif
c
        yl      = sprmn
        yu      = sprmx
        yll     = 2*yl - yu
        call map(xl,xr,yll,yu,xmap0,xmap1,ymap0,ymap1)
c
        call gridm(xl,yl,xrp,yu,'r-axis',  ntiks,-1)
        call gridm(xl,yl,xr, yu,'pressure',ntiks,+2)
        call line(xl, yu,xrp,yu)
        call line(xrp,yl,xrp,yu)
        if((yl-yz0)*(yu-yz0) .lt. 0.0) call linep(xl,yz0,xrp,yz0,ksp)
c
        call trace(x,st1,nx,1,1,0.0,0.0)
c
c
c 3.3 Current density vs. r
c
        do 50 i = 1,nx
        ivl     = i
        r1      = x(ivl)
        wjr     = psarray(ivl,jax)
c
        if    (wjr .gt. psilim) then
          st1(ivl)  = 0.0
c
        elseif(wjr .le. psilim) then
          if    (npst .eq. 0) then
            pp1       = sterpl(wjr,psimsh,spp, cspp, nprofl,npp,0,
     &                             rndff, ierp)
            ffp1      = sterpl(wjr,psimsh,sffp,csffp,nprofl,npp,0,
     &                             rndff, ierf)
            if(ierp .ne. 0) call abortjob
     &        ('plotfpj ',  7,   'Spline evaluation error: pp(psi)    '
     &        ,'ierp    ', ierp,     'ivl     ', ivl,      -1)
            if(ierf .ne. 0) call abortjob
     &        ('plotfpj ',  8,   'Spline evaluation error:ffp(psi)    '
     &        ,'ierf    ', ierf,     'ivl     ', ivl,      -1)
          elseif(npst .eq. 1) then
            pp1      = sterpp(wjr,psimsh,sp,  csp,  nprofl,npp,0,
     &                            rndff, ierpp)
            f1       = sterpl(wjr,psimsh,sf,  csf,  nprofl,npp,0,
     &                            rndff, ierf )
            fp1      = sterpp(wjr,psimsh,sf,  csf,  nprofl,npp,0,
     &                            rndff, ierfp)
            if(ierpp .ne. 0) call abortjob
     &        ('plotfpj ',  9,   'Differentiation error for p(i,j)    '
     &        ,'ierpp   ', ierpp,    'ivl     ', ivl,      -1)
            if(ierf  .ne. 0) call abortjob
     &        ('plotfpj ', 10,   'Interpolation error for f(i,j)      '
     &        ,'ierf    ', ierf,     'ivl     ', ivl,      -1)
            if(ierfp .ne. 0) call abortjob
     &        ('plotfpj ', 11,   'Differentiation error for f(i,j)    '
     &        ,'ierfp   ', ierfp,    'ivl     ', ivl,      -1)
            ffp1     = f1*fp1
          endif
c
          st1(ivl)  = -(r1*pp1 +ffp1/(r1*amu))
        endif
   50   continue
c
        srjmx   = st1(1)
        srjmn   = st1(1)
        do 55 i = 1,nx
        srjmx   = amax1(srjmx,st1(i))
        srjmn   = amin1(srjmn,st1(i))
   55   continue
c
        srjmn   = amin1(0.0,srjmn)
        srjmx   = amax1(0.0,srjmx)
c
        if(srjmn .eq. srjmx) then
          if    (srjmn .eq. 0.0) then
            srjmn  = -1.0
            srjmx  = +1.0
          elseif(srjmn .gt. 0.0) then
            srjmn  =  0.0
            srjmx  = +1.5*srjmx
          elseif(srjmn .lt. 0.0) then
            srjmn  =  1.5*srjmn
            srjmx  =  0.0
          endif
        endif
c
        yu      = srjmx
        yl      = srjmn
        yuu     = 2*yu-yl
        call map(xl,xr,yl,yuu,xmap0,xmap1,ymap0,ymap1)
c
        call gridm(xl,yl,xrp,yu,'r-axis', ntiks,+1)
        call gridm(xl,yl,xr, yu,'current',ntiks,+2)
        call line(xl, yu,xrp,yu)
        call line(xrp,yl,xrp,yu)
        if((yl-yz0)*(yu-yz0) .lt. 0.0) call linep(xl,yz0,xrp,yz0,ksp)
c
        call trace(x,st1,nx,1,1,0.0,0.0)
c
        call frame(0)
      endif
c
c
c
c 4.0 Plot ffprime and pprime vs. psi
c
c 4.1 Set up plots
c
      if(iplotm .lt. 9) return
c
      xmap0   = 0.20
      xmap1   = 0.87
      ymap0   = 0.20
      ymap1   = 0.87
c
      xl      = psimsh(   1  )
      xr      = psimsh(nprofl)
      xmov    = 0.0
      xrp     = xr + xmov
c
c
c 4.2 ffprime vs. psi (npst = 0) or f vs. psi (npst = 1)
c
      if    (npst .eq. 0) then
        sfpmx   = sffp(1)
        sfpmn   = sffp(1)
      elseif(npst .eq. 1) then
        sfpmx   = sf  (1)
        sfpmn   = sf  (1)
      endif

      do 60 i = 1,nprofl
      if    (npst .eq. 0) then
        sfpmx   = amax1(sfpmx,sffp(i))
        sfpmn   = amin1(sfpmn,sffp(i))
      elseif(npst .eq. 1) then
        sfpmx   = amax1(sfpmx,sf  (i))
        sfpmn   = amin1(sfpmn,sf  (i))
      endif
   60 continue
c
      sfpmn   = amin1(0.0,sfpmn)
      sfpmx   = amax1(0.0,sfpmx)
c
      if(sfpmn .eq. sfpmx) then
        if    (sfpmn .eq. 0.0) then
          sfpmn  = -1.0
          sfpmx  = +1.0
        elseif(sfpmn .gt. 0.0) then
          sfpmn  =  0.0
          sfpmx  = +1.5*sfpmx
        elseif(sfpmn .lt. 0.0) then
          sfpmn  =  1.5*sfpmn
          sfpmx  =  0.0
        endif
      endif
c
      yl      = sfpmn
      yu      = sfpmx
      yll     = 2.*yl-yu
      call map(xl,xr,yll,yu,xmap0,xmap1,ymap0,ymap1)
c
      call gridm(xl,yl,xrp,yu,'psi',     ntiks,-1)
      call gridm(xl,yl,xr, yu,'ff-prime',ntiks,+2)
      call line(xl,yu, xrp,yu)
      call line(xrp,yl,xrp,yu)
      if((yl-yz0)*(yu-yz0) .lt. 0.0) call linep(xl,yz0,xrp,yz0,ksp)
c
      if(npst .eq. 0) call trace(psimsh,sffp,nprofl,1,1,0.0,0.0)
      if(npst .eq. 1) call trace(psimsh,sf,  nprofl,1,1,0.0,0.0)
c
c
c 4.3 pprime vs. psi (npst = 0) or p vs. psi (npst = 1)
c
      if(npst .eq. 0) then
        sdpmx   = spp (1)
        sdpmn   = spp (1)
      elseif(npst .eq. 1) then
        sdpmx   = sp  (1)
        sdpmn   = sp  (1)
      endif
      do 65 i = 1,nprofl
      if    (npst .eq. 0) then
        sdpmx   = amax1(sdpmx ,spp (i))
        sdpmn   = amin1(sdpmn ,spp (i))
      elseif(npst .eq. 1) then
        sdpmx   = amax1(sdpmx ,sp  (i))
        sdpmn   = amin1(sdpmn ,sp  (i))
      endif
   65 continue
c
      sdpmn   = amin1(0.0,sdpmn)
      sdpmx   = amax1(0.0,sdpmx)
c
      if(sdpmn .eq. sdpmx) then
        if    (sdpmn .eq. 0.0) then
          sdpmn  = -1.0
          sdpmx  = +1.0
        elseif(sdpmn .gt. 0.0) then
          sdpmn  =  0.0
          sdpmx  = +1.5*sdpmx
        elseif(sdpmn .lt. 0.0) then
          sdpmn  =  1.5*sdpmn
          sdpmx  =  0.0
        endif
      endif
c
      yu      = sdpmx
      yl      = sdpmn
      yuu     = 2.*yu-yl
      call map(xl,xr,yl,yuu,xmap0,xmap1,ymap0,ymap1)
c
      call gridm(xl,yl,xrp,yu,'psi',    ntiks,+1)
      call gridm(xl,yl,xr, yu,'p-prime',ntiks,+2)
      call line(xl,yu, xrp,yu)
      call line(xrp,yl,xrp,yu)
      if((yl-yz0)*(yu-yz0) .lt. 0.0) call linep(xl,yz0,xrp,yz0,ksp)
c
      if(npst .eq. 0) call trace(psimsh,spp,  nprofl,1,1,0.0,0.0)
      if(npst .eq. 1) call trace(psimsh,sp,   nprofl,1,1,0.0,0.0)
c
c
c
c 5.0 Return and end
c
      return
      end
      subroutine gridm(x1,y1,x2,y2,axlabl,ntiks,nv)
c
c     nv   = +1: draw horizontal axis with labels
c     nv   = -1: draw horizontal axis without labels
c     nv   = +2: draw vertical axis with labels
c     nv   = -2: draw vertical axis without labels
c
      character*(*) axlabl
      character*80  string
c
      nptv    =  ntiks + 1
      nvabs   =  iabs(nv)
      defy0   =  0.040
      defx0   =  0.020
      dshyysc = -5.000
      dshyylb = -6.250
      dshxxsc = -6.600
      dshxxlb = -7.000
      xlbpos  = +0.60
      ylbpos  = +0.60
c
c
c     Horizontal axis
      if    (nvabs .eq. 1) then
c
        call line(x1,y1,x2,y1)
c
        delx    = (x2 - x1)/(nptv-1.0)
        defy    = defy0*(y2 - y1)
        dyh     = 0.50*defy
        yy      = y1
        yyu     = yy + dyh
        yysc    = y1 + dshyysc*defy
c
        do 10 i = 1,nptv
        xx      = x1 + (i-1)*delx
        if(i .gt.  1) call line(xx,yy,xx,yyu)
        if(nv .gt. 0) then
          call setlch(xx,yysc,1,0,1,-1)
          write(string,100) xx
          call wrtstr(string,1)
        endif
   10   continue
c
        if(nv .gt. 0) then
          xxlb    = x1 + xlbpos*(x2-x1)
          yylb    = y1 + dshyylb*defy
          call setlch(xxlb,yylb,1,1,0,-1)
          write(string,110) axlabl
          call wrtstr(string,1)
        endif
c
c
c     Vertical axis
      elseif(nvabs .eq. 2) then
c
        call line(x1,y1,x1,y2)
c
        dely    = (y2 - y1)/(nptv-1.0)
        defx    = defx0*(x2 - x1)
        dxh     = 0.50*defx
        xx      = x1
        xxr     = xx + dxh
        xxsc    = x1 + dshxxsc*defx
c
        do 20 i = 1,nptv
        yy      = y1 + (i-1)*dely
        if(i .gt.  1) call line(xx,yy,xxr,yy)
        if(i .lt. nptv  .and.  nv .gt. 0) then
          call setlch(xxsc,yy,1,0,0,-1)
          write(string,200) yy
          call wrtstr(string,1)
        endif
   20   continue
c
        if(nv .gt. 0) then
          yylb    = y1 + ylbpos*(y2-y1)
          xxlb    = x1 + dshxxlb*defx
          call setlch(xxlb,yylb,1,1,1,-1)
          write(string,210) axlabl
          call wrtstr(string,1)
        endif
      endif
c
c
      return
c
  100 format(f7.4)
  110 format(a8)
  200 format(e11.4)
  210 format(a8)
      end
