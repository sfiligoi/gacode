C------------------------------------------
      subroutine nbi_rbloat

      use nbi_com
      implicit NONE

      integer iznbmri

      call mcgrid_getnumr(id_mcgrid,NZNBMR,iznbmri)
      XBMBND=XMINBM
      NZNBLO=IFIX(XMINBM/DXIZC(LCENTR) + 0.5)

      end
C******************** END FILE BMZINI.FOR ; GROUP BMZINI ***************
      subroutine nbi_shrink

      use nbi_com
      implicit NONE

      NZNBMR=NZONES/NXSKPB
      NZNBLO=NZONES
      XBMBND=1.0

      end
C******************** START FILE BMZINI.FOR ; GROUP BMZINI *************
C..........................................
C
      SUBROUTINE BMZINI(ierr)
C
      use nbi_com
      use nbi_dimensions
      use nbi_idhdr
C
      implicit NONE
C
      integer, intent(out) :: ierr  ! completion code, 0=OK
C
C  beam code version of old TRANSP routine
C  get 2d MC grid data from XPLASMA...
C  dmc 7 Dec 2001
C
C--------------------------------------------------
      integer iznbmr,iznbmri
      integer iudsym
C
      integer izones,jznbmr,iabort,i,ith,iwarn
C
      integer, parameter :: mjth = 40
C
C--------------------------------------------------
C
      ierr = 0
C
      print *,'in bmzini, iudsym = ',iudsym
      call mcgrid_getsym(id_mcgrid,iudsym)
      print *,'in bmzini, iudsym ,after getsym= ',iudsym
      if(iudsym.eq.1) then
         nlsym2b=.TRUE.
      else
         nlsym2b=.FALSE.
      endif
C
      call mcgrid_getnumr(id_mcgrid,iznbmr,iznbmri)
      nxskpb=nzones/iznbmri
      nxtbzn=iznbmr-iznbmri
C
C  test...
C
      IZONES = (NZONES/NXSKPB) * NXSKPB
      IF(IZONES.NE.NZONES) THEN
         write(nonlin,9901) NZONES,NXSKPB
         ierr=1
         return
      ENDIF
 9901 format(' ?BMZINI -- NZONES=',i3,' and NXSKPB=',i3/
     >'   NXSKPB must divide NZONES without remainder for TRANSP'/
     >'   beam code indexing to work properly.')
C
      DXITK = 1.0/iznbmri
      NZNBMR = iznbmr
      NZNBMRI = iznbmri
      jznbmr = nznbmr+lcm1
C
C  count the zones; double if not updown symmetric.
C
      if(allocated(nthzsm)) deallocate(nthzsm)
      allocate(nthzsm(jznbmr))
C
      NTHZSM(1:LCM1) = 0
      call mcgrid_getnumzns(id_mcgrid,NTHZSM(lcentr:jznbmr))
C
C  convert to sum
C
      do i=lcentr,jznbmr
         nthzsm(i)=nthzsm(i-1)+nthzsm(i)
      enddo
C
      NFBZNS=NTHZSM(JZNBMR)
      NFBZNSI=NTHZSM(JZNBMR-NXTBZN)
      NTHZSM0=NTHZSM(lcentr)
C
C--------------------------------------------------------
C  *** if this is the first call ***
C
      if(.NOT.NBI_INIT_FLAG) then

         ! set some array dimensions

         MIMXBZ=NFBZNS
         MIMXBZF=NFBZNSI
         MINB=max(nptcl_max,nptcls,nptclf,nptclh) + 1000

         call nbi_alloc2_init   ! allocate rest of NBI_COM ...!

         !  nbi_alloc2_init also sets NBI_INIT_FLAG

         call NBINIT0           ! setup NBNDEX & other misc. initializations

         call nbi_lun_msgs(nonlin)

         call boxn0_setup(iwarn) ! beam-in-box calculation setup:  ignore error

         call nbi_random_init(nseed,0,1)  ! init C. Karney RNG
         
         if(ierr.ne.0) return

      else if(nbi_idhdr_skip) then

         MIMXBZ=NFBZNS
         MIMXBZF=NFBZNSI

      endif
C
C--------------------------------------------------------
C
      IF(NFBZNS.NE.MIMXBZ) THEN
        write(nonlin,*) '?BMZINI -- NFBZNS/MIMXBZ mismatch: ',
     >      nfbzns,mimxbe
        IERR=1
      ENDIF
      IF(NFBZNSI.NE.MIMXBZF) THEN
        write(nonlin,*) '?BMZINI -- NFBZNSI/MIMXBZF mismatch: ',
     >      nfbznsi, mimxbzf
        IERR=1
      ENDIF
      IF(IERR.EQ.1) return
C
      XBMBND = NZNBMR * DXITK
      XMINBM = XBMBND
C
      NZNBLO=IFIX(XMINBM/DXIZC(LCENTR) + 0.5)
      if(allocated(xiblo)) deallocate(xiblo)
      allocate(xiblo(mtbl))
      xiblo(1:lcm1)=0.0
      do i=1,mtbl-lcm1
         xiblo(lcm1+i)=(i-1)*DXIZC(LCENTR)
      enddo
      xiblo(lep1)=xi(lep1,1)
C
C  a poloidal angle grid -- added dmc 6 June 1994
C
      if(NLSYM2B) then
         thbdy0=0.0
         thbdy1=xpi
         nthzons=mjth/2
      else
         thbdy0=-xpi
         thbdy1=xpi
         nthzons=mjth
      endif
C
      if(allocated(thzons)) deallocate(thzons)
      if(allocated(gxtabl)) deallocate(gxtabl)
      allocate(thzons(nthzons),gxtabl(nthzons,mtbl))
      thzons=0.0
C
      do ith=1,nthzons
         thzons(ith)=thbdy0 +
     >        (thbdy1-thbdy0)*(float(ith)-0.5)/float(nthzons)
      enddo
C
      xbjaco=0.0                        ! clear XBJACO
C
      RETURN
      END
