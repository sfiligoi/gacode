! This is another toq12 interface.

      SUBROUTINE set_gfile_to_toq
!---------------------------------------------------------------------
! 1)reads existing dskfixb
! 2)uses psir (in normalized form from Onetwo) to interpolate pprim
!   and ffprim onto toq grid
! 3)creates new dskfixb file
! 4)spwans script fidb129xy.sc (which copies dskfixb into toqfile
!   that is read by fixb code)
! 5)fixb code creates standard eqdsk file and
!   eqdskfilename is set to this eqdsk name
!-------------------------------------------------------------HSJ-04/03/05/

!save this
      USE param
      USE fusion
      USE numbrs
      USE ions
      USE solcon
      USE soln
      USE mhdpar
      USE nub
      USE nub2
      USE extra
      USE mesh
      USE geom
      USE psig
      USE rhog
      USE soln2d
      USE ename, ONLY: eqdskfilename
      USE bicube
      USE etc
      USE io, ONLY : ncrt,nout,nitre
      USE  toq_12_interface,  ONLY : toq_base_out
      USE ext_prog_info, ONLY : get_toq_base
      IMPLICIT  INTEGER (i-n), REAL*8 (a-h, o-z)

      INCLUDE 'imsl.i'
!      INCLUDE 'etc.i'
      CHARACTER *50 err_string 
      CHARACTER*50 command
      REAL*8 xnorm_toq(259,514),znorm_toq(259,514)
      REAL*8 fnorm_toq(259),fprimenorm_toq(259)
      REAL*8 pponly_toq(259), psir_norm(kj)
      INTEGER i,j,nfixb, npsi_toq, nbdry_toq,len_str
      REAL*8 amu,pprim_toq, curtot_toq, xval_toq(259)
      REAL*8 ffpnorm_toq(259),btorfixb_toq,rcfixb_toq

      err_string = 'ERROR,file dskfixb required but not found'
      nfixb=73
      CALL getioun(nfixb,nfixb)
      OPEN(unit=nfixb,file="dskfixb",status="unknown",ERR =300 )
      READ(nfixb,'(a)')
      READ(nfixb,'(a)')
      READ(nfixb,100) npsi_toq,nbdry_toq
      READ(nfixb,'(a)')
      READ(nfixb,200) rcfixb_toq,xnorm_toq(1,1),znorm_toq(1,1), &
            btorfixb_toq
      READ(nfixb,200) curtot_toq
      READ(nfixb,'(a)')
      READ(nfixb,200) (xval_toq(i),i=1,npsi_toq) ! psi normalized
      READ(nfixb,'(a)')
      READ(nfixb,200) (ffpnorm_toq(i),i=1,npsi_toq) ! ff'
      READ(nfixb,'(a)')
      READ(nfixb,200) (pponly_toq(i),i=1,npsi_toq) ! p'
      READ(nfixb,'(a)')
      READ(nfixb,200) (xnorm_toq(npsi_toq,j),j=1,nbdry_toq)
      READ(nfixb,'(a)')
      READ(nfixb,200) (znorm_toq(npsi_toq,j),j=1,nbdry_toq)


      CLOSE(nfixb)
      CALL giveupus(nfixb)


      DO j=1, nj
           psir_norm(j) = (psir(j)-psir(1))/(psir(nj)-psir(1))
      ENDDO

      CALL intrp (1, 1, psir_norm,  pprim, nj, &
             xval_toq, pponly_toq, npsi_toq)
      CALL intrp (1, 1, psir_norm,  ffprim, nj, &
                  xval_toq, ffpnorm_toq, npsi_toq)
      CALL getioun(nfixb,nfixb)
      OPEN(unit=nfixb,file="dskfixb",status="unknown")

! write dskfixb

      WRITE(nfixb,'("toq equilibria info for fixb and refeqd")')
      WRITE(nfixb,'("npsi,ngeom")')
      WRITE(nfixb,100) npsi_toq, nbdry_toq ! number of mesh points
      WRITE(nfixb,'("rc,xaxis,zaxis,bvac at rc,totcur")')
      WRITE(nfixb,200) rcfixb_toq,xnorm_toq(1,1), &
                        znorm_toq(1,1),btorfixb_toq
      WRITE(nfixb,200) curtot_toq
      WRITE(nfixb,'("normalized psi from center to edge")')
      WRITE(nfixb,200) (xval_toq(i),i=1,npsi_toq) ! psi normalized
      WRITE(nfixb,'("ffprime from center to edge")')
      WRITE(nfixb,200) (ffpnorm_toq(i),i=1,npsi_toq) ! ff'
      WRITE(nfixb,'("pprime from center to edge")')
      WRITE(nfixb,200) (pponly_toq(i),i=1,npsi_toq) ! the new p'
      WRITE(nfixb,'("R of boundary points--equal arc length")')
      WRITE(nfixb,200) (xnorm_toq(npsi_toq,j),j=1,nbdry_toq)
      WRITE(nfixb,'("Z of boundary points--equal arc length")')
      WRITE(nfixb,200) (znorm_toq(npsi_toq,j),j=1,nbdry_toq)

      CLOSE(nfixb)
      CALL giveupus(nfixb)

      CALL get_toq_base(ncrt,nout,toq_base_out,len_str)
      command = ADJUSTL(toq_base_out(1:LEN_TRIM( toq_base_out)))
      command = command(1:LEN_TRIM(command))//'/fixb129xy.sc'
      command =  command(1:LEN_TRIM(command))//' 12345 '//'1234'
      IF(ISHELL(command).LT.0) &
        CALL STOP("failure to invoke fixbdry to create gfile", 1)
!      eqdskfilename = 'g812345.01234'    

 100  FORMAT(2i5)
 200  FORMAT(4e19.12)

      RETURN
300   WRITE(nitre,FMT ='(a)')err_string
      PRINT *,err_string 
      CALL STOP('set_gfile_to_toq',1)
      END





            SUBROUTINE read_toq_output
!-------------------------------------------------------HSJ/8/05/03----
!--- read the results  from toq equilibrium calculation (file dskonetwo)
!----------------------------------------------------------------------
      USE param
      USE io, ONLY : ncrt,nout
      USE solcon
      USE soln
      USE contour
      USE limiter
      USE mhdpar
      USE mhdgrid
      USE rf
      USE ename
      USE extra
      USE yoka
      USE numbrs
      USE mesh
      USE machin
      USE geom
      USE psig
      USE mhdcom
      USE rhog
      USE soln2d
      USE tordlrot
      USE constnts
      USE shapctr
      USE etc
      USE flxav
      USE metrics
      USE neo2d
      USE neo2dp
      USE gpsi
      USE  toq_12_interface,  ONLY : toq_base_out
      USE  ext_prog_info,     ONLY : get_toq_base
      IMPLICIT  INTEGER (i-n), REAL*8 (a-h, o-z)

!      INCLUDE 'etc.i'
!      INCLUDE 'small.i'

      INTEGER  i, j, npsi_toq, nthet_toq, ndskonetwo
      INTEGER  nbndry_toq, nw_fixb, nh_fixb, nfixb,len_str
      REAL*8 rax_toq, zax_toq, rmin_toq, rmax_toq
      REAL*8 ali_toq, flim_toq, zmin_toq, zmax_toq
      REAL*8  widp_toq, hitp_toq, circum_toq, kappa_toq
      REAL*8,DIMENSION(:),ALLOCATABLE :: cxarea_toq,sfarea_toq,          & 
               torflux_toq,curden_toq,xbndry_toq,zbndry_toq, rho_toq,    & 
               drhodpsi_toq,vprime_toq, fpsi_toq, psivolp_toq,           & 
               invr2av_toq, rav_toq, r2av_toq, invb2av_toq, invbav_toq,  & 
               bav_toq, b2av_toq,psic_toq, psiv_toq,                     & 
               gradrho2r2_toq, fnorm_toq, xhm2_toq, xi11_toq,dummy,      & 
               xi33_toq, xips_toq,epsp_toq,gradpsi2_toq, elong_toq,      & 
               gradpsi_toq,grth_toq, bsq_toq, bmsq_toq,  grbmsq_toq        
      REAL*8   btgeom_toq,rdim_fixb, zdim_fixb, rmin_fixb, zmid_fixb,    &
               psiax_fixb, dr_fixb, dz_fixb,rgeom_toq,                   &
               a1, a2, a3, a4, psib_fixb
      REAL*8,  DIMENSION(:,:),ALLOCATABLE :: gfm_toq,psi_fixb
      CHARACTER*8 ntitle(6), dat
      CHARACTER(len =80) line, command

      ndskonetwo=45
      CALL getioun(ndskonetwo,ndskonetwo)
      OPEN(unit=ndskonetwo,file="dskonetwo",status="unknown", &
                   err=1000)

      READ(ndskonetwo,600) line

      READ(ndskonetwo, 100) npsi_toq, ali_toq, flim_toq, &
                          rax_toq, zax_toq
      READ(ndskonetwo,600) line
      READ(ndskonetwo, 100) nbndry_toq, rmin_toq, rmax_toq, &
                              rgeom_toq, btgeom_toq

      IF(.NOT. ALLOCATED(cxarea_toq))ALLOCATE(cxarea_toq(npsi_toq))     !!1
      IF(.NOT. ALLOCATED(sfarea_toq))ALLOCATE(sfarea_toq(npsi_toq))     !!2
      IF(.NOT. ALLOCATED(elong_toq))ALLOCATE(elong_toq(npsi_toq))       !!3
      IF(.NOT. ALLOCATED(epsp_toq))ALLOCATE(epsp_toq(npsi_toq))         !!4
      IF(.NOT. ALLOCATED(torflux_toq))ALLOCATE(torflux_toq(npsi_toq))   !!5
      IF(.NOT. ALLOCATED(curden_toq))ALLOCATE(curden_toq(npsi_toq))     !!6
      IF(.NOT. ALLOCATED(grth_toq))ALLOCATE(grth_toq(npsi_toq))         !!7
      IF(.NOT. ALLOCATED(bmsq_toq))ALLOCATE(bmsq_toq(npsi_toq))         !!8
      IF(.NOT. ALLOCATED(grbmsq_toq))ALLOCATE(grbmsq_toq(npsi_toq))     !!9
      IF(.NOT. ALLOCATED(psiv_toq))ALLOCATE(psiv_toq(npsi_toq))         !10
      IF(.NOT. ALLOCATED(invr2av_toq))ALLOCATE(invr2av_toq(npsi_toq))   !11
      IF(.NOT. ALLOCATED(rav_toq))ALLOCATE(rav_toq(npsi_toq))           !12
      IF(.NOT. ALLOCATED(r2av_toq))ALLOCATE(r2av_toq(npsi_toq))         !13
      IF(.NOT. ALLOCATED(invb2av_toq))ALLOCATE(invb2av_toq(npsi_toq))   !14
      IF(.NOT. ALLOCATED(b2av_toq))ALLOCATE(b2av_toq(npsi_toq))         !15
      IF(.NOT. ALLOCATED(psivolp_toq))ALLOCATE(psivolp_toq(npsi_toq))   !16
      IF(.NOT. ALLOCATED(gradpsi_toq))ALLOCATE(gradpsi_toq(npsi_toq))   !17
      IF(.NOT. ALLOCATED(gradpsi2_toq))ALLOCATE(gradpsi2_toq(npsi_toq)) !18
      IF(.NOT. ALLOCATED(vprime_toq))ALLOCATE(vprime_toq(npsi_toq))     !19
      IF(.NOT. ALLOCATED(fpsi_toq))ALLOCATE(fpsi_toq(npsi_toq))         !20
      IF(.NOT. ALLOCATED(rho_toq))ALLOCATE(rho_toq(npsi_toq))           !21
      IF(.NOT. ALLOCATED(xhm2_toq))ALLOCATE(xhm2_toq(npsi_toq))         !22
      IF(.NOT. ALLOCATED(xi11_toq))ALLOCATE(xi11_toq(npsi_toq))         !23
      IF(.NOT. ALLOCATED(xi33_toq))ALLOCATE(xi33_toq(npsi_toq))         !24
      IF(.NOT. ALLOCATED(xips_toq))ALLOCATE(xips_toq(npsi_toq))         !25
      IF(.NOT. ALLOCATED(bsq_toq)) ALLOCATE(bsq_toq(npsi_toq))          !26
      IF(.NOT. ALLOCATED(xbndry_toq))ALLOCATE(xbndry_toq(nbndry_toq))   !27
      IF(.NOT. ALLOCATED(zbndry_toq))ALLOCATE(zbndry_toq(nbndry_toq))   !28
      IF(.NOT. ALLOCATED(psic_toq))ALLOCATE(psic_toq(npsi_toq))        !29 not used ??
      IF(.NOT. ALLOCATED(drhodpsi_toq))ALLOCATE(drhodpsi_toq(npsi_toq)) !30 not used ??
      IF(.NOT. ALLOCATED(invbav_toq))ALLOCATE(invbav_toq(npsi_toq))     !31 not used ??
      IF(.NOT. ALLOCATED(gradrho2r2_toq))ALLOCATE(gradrho2r2_toq(npsi_toq)) !32 not used ??
      IF(.NOT. ALLOCATED(fnorm_toq))ALLOCATE(fnorm_toq(npsi_toq))       !33 not used ??
      IF(.NOT. ALLOCATED(bav_toq))ALLOCATE(bav_toq(npsi_toq))           !34 not used ??
      IF(.NOT. ALLOCATED(gfm_toq))ALLOCATE(gfm_toq(3,npsi_toq))         !35



      READ(ndskonetwo,600) line
      READ(ndskonetwo, 100) nthet_toq, widp_toq, hitp_toq, &
                                 circum_toq, kappa_toq
      READ(ndskonetwo,600) line
      READ(ndskonetwo,610) (j,psic_toq(i),psiv_toq(i),gradrho2r2_toq(i), &
                               fnorm_toq(i),i=1,npsi_toq)
      READ(ndskonetwo,600) line
      READ(ndskonetwo,610) (j,invr2av_toq(i),psivolp_toq(i),rav_toq(i), &
                              r2av_toq(i),i=1,npsi_toq)
      READ(ndskonetwo,600) line
      READ(ndskonetwo,610) (j,invb2av_toq(i),invbav_toq(i),bav_toq(i), &
                              b2av_toq(i),i=1,npsi_toq)
      READ(ndskonetwo,600) line

      READ (ndskonetwo,610) (j, rho_toq(i),drhodpsi_toq(i), &
                              vprime_toq(i),fpsi_toq(i), i=1,npsi_toq)
      READ (ndskonetwo,600 )  line
      READ (ndskonetwo,610) (j, xhm2_toq(i), xi11_toq(i), xi33_toq(i), &
                              xips_toq(i), i=1, npsi_toq)
      READ (ndskonetwo,600 ) line
      READ (ndskonetwo,610) (j, epsp_toq(i), elong_toq(i), &
             gradpsi_toq(i),  gradpsi2_toq(i), i = 1, npsi_toq)
      READ (ndskonetwo,600 ) line
      READ(ndskonetwo, 610) (i, cxarea_toq(j), sfarea_toq(j), &
                          torflux_toq(j), curden_toq(j),j=1, npsi_toq)
      READ (ndskonetwo,600 ) line
      READ(ndskonetwo, 610) (i, grth_toq(j), bsq_toq(j), &
                           bmsq_toq(j), grbmsq_toq(j), j=1, npsi_toq)
      READ (ndskonetwo,600 ) line
      READ (ndskonetwo, 611) (i, (gfm_toq(j, jj), j=1,3), &
                                 jj=1, npsi_toq)
      READ (ndskonetwo,600 ) line
      READ(ndskonetwo, 620) (xbndry_toq(i), zbndry_toq(i), &
          i=1, nbndry_toq)
!    npsi = npsi_toq         ! in flxav.i
     CLOSE(ndskonetwo)



      rma = rax_toq
      zma = zax_toq
      xax(1) = rax_toq
      yax(1) = zax_toq
        xmagn1 = rma
        ymagn1 = zma
      flim = flim_toq
      psiaxis = psiv_toq(1)
      psibdry = psiv_toq(npsi_toq)
      elongax = elong_toq(1)
      ali = ali_toq
        rgeom = rgeom_toq
        btgeom = btgeom_toq
        hitep = hitp_toq
        widep = widp_toq
        kappa = kappa_toq
        circum = circum_toq
      DO  i = npsi_toq, 1, -1
      psival(i) = psiv_toq(npsi_toq+1-i)
      ratave(i) = invr2av_toq(npsi_toq+1-i)
      ravg(i) = rav_toq(npsi_toq+1-i)
      ratavei(i) = r2av_toq(npsi_toq+1-i)
       bsqinvavg(i) = invb2av_toq(npsi_toq+1-i)
      bsq_avg(i) = b2av_toq(npsi_toq+1-i)
      psivolp(i) = psivolp_toq(npsi_toq+1-i)
      grho1npsi(i) = gradpsi_toq(npsi_toq+1-i)
      grho2npsi(i) = gradpsi2_toq(npsi_toq+1-i)
      vprime(i) =  vprime_toq(npsi_toq+1-i)
      fpsi(i) = fpsi_toq(npsi_toq+1-i)
      rho(i) = rho_toq(npsi_toq+1-i)
      xhm2p(i) = xhm2_toq(npsi_toq+1-i)
      xi11p(i) = xi11_toq(npsi_toq+1-i)
      xi33p(i) = xi33_toq(npsi_toq+1-i)
      xipsp(i) = xips_toq(npsi_toq+1-i)
      epsp(i) = epsp_toq(npsi_toq+1-i)
      elongx(i) =      elong_toq(npsi_toq+1-i)
        sfareanpsi(i) = sfarea_toq(npsi_toq+1-i)
        cxareanpsi(i) = cxarea_toq(npsi_toq+1-i)
        torfluxnpsi(i) = torflux_toq(npsi_toq+1-i)
        grth(i) = grth_toq(npsi_toq+1-i)
        bsq(i) = bsq_toq(npsi_toq+1-i)
        bmsq(i) = bmsq_toq(npsi_toq+1-i)
        grbmsq(i) = grbmsq_toq(npsi_toq+1-i)
        DO jj=1,3
         gfm(jj,i) = gfm_toq(jj,npsi_toq+1-i)
        ENDDO
      ENDDO




        DO i=1,nbndry_toq
           rplasbdry(i) = xbndry_toq(i)
           zplasbdry(i) = zbndry_toq(i)
        ENDDO





!     To get psi(i,j) and psi1d from fixb.
      CALL get_toq_base(ncrt,nout,toq_base_out,len_str)
      command = ADJUSTL(toq_base_out(1:LEN_TRIM( toq_base_out)))
      command = command(1:LEN_TRIM(command))//'/fixb129xy.sc'
      command =  command(1:LEN_TRIM(command))//' 12345 '//'1234'
      PRINT *,'running : '//command
      IF(ISHELL(command).LT.0) &
        CALL STOP("failure to invoke fixbdry to create gfile", 1)


      nfixb = 44
      CALL getioun(nfixb,nfixb)
      OPEN(unit=nfixb, file='g812345.01234', status='unknown', &
           err=1000)
      READ(nfixb,8190) (ntitle(i), i=1,6), j, nw_fixb, nh_fixb, dat
      READ(nfixb,8200) rdim_fixb, zdim_fixb,a1,rmin_fixb, zmid_fixb
      READ(nfixb, 8200) a1, a2, psiax_fixb,psib_fixb, a4
      IF(.NOT. ALLOCATED(dummy)) ALLOCATE(dummy(MAX(nw_fixb,5)))
      IF(.NOT. ALLOCATED(psi_fixb)) ALLOCATE(psi_fixb(nw_fixb,nh_fixb))
      READ(nfixb, 8200) (dummy(i), i=1,5)
      READ(nfixb, 8200) (dummy(i), i=1,5)
      READ(nfixb, 8200) (dummy(i), i=1,nw_fixb)
      READ(nfixb, 8200) (dummy(i), i=1,nw_fixb)
      READ(nfixb, 8200) (dummy(i), i=1,nw_fixb)
      READ(nfixb, 8200) (dummy(i), i=1,nw_fixb)
      READ(nfixb, 8200) ((psi_fixb(i, j),i=1, nw_fixb), j=1,nh_fixb)
      CLOSE(nfixb)
      CALL giveupus(nfixb)


      dr_fixb = rdim_fixb/(nw_fixb-1)
      dz_fixb = zdim_fixb/(nh_fixb-1)
      rmhdgrid(1)=rmin_fixb
      zmhdgrid(1)=zmid_fixb - 0.5*zdim_fixb
      DO i =2,nw_fixb
      rmhdgrid(i) = rmhdgrid(1)+(i-1)*dr_fixb
      ENDDO
      DO i=1, nh_fixb
         zmhdgrid(i) = zmhdgrid(1)+(i-1)*dz_fixb
      ENDDO
      DO j=1, nh_fixb
         DO i=1,nw_fixb
            k=(i-1)*nh_fixb+j
            psi(i,j) = psi_fixb(i,j)+(psiaxis-psiax_fixb)
            psi1d(k) = psi(i,j)
            p(i,j) = psi(i,j)
         ENDDO
       ENDDO
       rplasmax = rmhdgrid(nw_fixb)
        rplasmin = rmin_fixb
        zplasmax = zmhdgrid(nh_fixb)
        zplasmin = zmhdgrid(1)
        rmax = 100.0*rplasmax
        rmin = 100.0*rplasmin
        zmax = 100.0*zplasmax
        zmin =100.0*zplasmin
      PRINT*, "rmax=", rmax,"rmin=", rmin
      PRINT*, "zmax=", zmax,"zmin=", zmin
      isignpsi = -1
      CALL psiscale (psiaxis, psibdry, psir, nj, isignpsi)
!            call intrp (1, 1, psir, curden, nj,
!     &      psiv_toq, curden_toq, npsi_toq)

      
      DEALLOCATE(cxarea_toq)    ; DEALLOCATE(sfarea_toq)   ; DEALLOCATE(elong_toq )
      DEALLOCATE(epsp_toq )     ; DEALLOCATE(torflux_toq ) ; DEALLOCATE(grth_toq )
      DEALLOCATE( bmsq_toq )    ; DEALLOCATE(grbmsq_toq )  ; DEALLOCATE(bsq_toq )
      DEALLOCATE(xbndry_toq)    ; DEALLOCATE(zbndry_toq )  ; DEALLOCATE(psiv_toq )
      DEALLOCATE(invr2av_toq)   ; DEALLOCATE( rav_toq)     ; DEALLOCATE(r2av_toq )
      DEALLOCATE(invb2av_toq)   ; DEALLOCATE(bav_toq)      ; DEALLOCATE(psivolp_toq)
      DEALLOCATE(gfm_toq )      ; DEALLOCATE(b2av_toq)     ; DEALLOCATE(gradpsi_toq)
      DEALLOCATE(gradpsi2_toq)  ; DEALLOCATE(vprime_toq)   ; DEALLOCATE(xhm2_toq)
      DEALLOCATE(xi11_toq)      ; DEALLOCATE(fpsi_toq)     ; DEALLOCATE(rho_toq)
      DEALLOCATE(xi33_toq)      ; DEALLOCATE(xips_toq)     ; DEALLOCATE(dummy)
      DEALLOCATE(psic_toq)      ; DEALLOCATE(drhodpsi_toq) ; DEALLOCATE(invbav_toq)
      DEALLOCATE(gradrho2r2_toq); DEALLOCATE(psi_fixb)     ; DEALLOCATE(curden_toq)    
      DEALLOCATE(fnorm_toq)
 100  FORMAT(i5, 4e16.8)
 600  FORMAT(80a)
 610  FORMAT(i3, 4e16.8)
 611  FORMAT(i3, 3e16.8)
 620  FORMAT(2e16.8)
 8190 FORMAT(6a8,3i4,t73,a)
 8200 FORMAT(5e16.9)

      RETURN

1000  CALL STOP("error in open dskonetwo or fixb, check it", 2)
      END






      SUBROUTINE use_eqdsk
!------------------------------------------------------------------------
! 
! We first read dskeqdata, then substitute pprime and ffprime with new
! values from ONETWO run (after a interp from grid nj to npsi_dsk),
! finally write dskeqdata for following usage(ieqtoq=0) in toq.
! Dzhou 3/12/04
!-------------------------------------------------------------------------

      USE param
      USE fusion
      USE ions
      USE solcon
      USE soln
      USE mhdpar
      USE nub
      USE ename, only  :  eqdskfilename
      USE numbrs
      USE nub2
      USE etc
      USE extra
      USE mesh
      USE geom
      USE psig
      USE rhog
      USE bicube
      USE io,       only : nitre
      IMPLICIT  INTEGER (i-n), REAL*8 (a-h, o-z)

      INCLUDE 'imsl.i'
!      INCLUDE 'etc.i'


      INTEGER neqdata
      INTEGER npsi_dsk,nbndry_dsksk
      REAL*8,DIMENSION(:),ALLOCATABLE ::  psiv_dsk, p_dsk, pp_dsk,        &
             f_dsk, ffp_dsk, q_dsk,grot_dsk, grotp_dsk,                   &
             cd_dsk, psiv_norm,dummy,xbndry, zbndry
      REAL*8  psir_norm(kj)
      CHARACTER*129 line
      CHARACTER *48 error_string
      COMMON /dskeqdat/ modelcd, modelbnd




      error_string = 'ERROR: dskeqdata not found'
      neqdata=90
      CALL getioun(neqdata,neqdata)
      OPEN(unit=neqdata,file='dskeqdata',status='old',iostat=ios)
      IF(ios.NE.0)THEN
         write(nitre,FMT='(a)')error_string
         PRINT *,error_string
         CALL STOP('use_eqdsk ',1)
      ENDIF

      READ(neqdata,'(a)') line
      READ(neqdata,'(a)') line
      READ(neqdata,'(i3)') npsi_dsk
      READ(neqdata,'(i3)') nbndry_dsk
      IF(.NOT. ALLOCATED(psiv_dsk))ALLOCATE(psiv_dsk(npsi_dsk))
      IF(.NOT. ALLOCATED(p_dsk))ALLOCATE(p_dsk(npsi_dsk))
      IF(.NOT. ALLOCATED(pp_dsk))ALLOCATE(pp_dsk(npsi_dsk))
      IF(.NOT. ALLOCATED(f_dsk))ALLOCATE(f_dsk(npsi_dsk))
      IF(.NOT. ALLOCATED(ffp_dsk))ALLOCATE(ffp_dsk(npsi_dsk))
      IF(.NOT. ALLOCATED(q_dsk))ALLOCATE(q_dsk(npsi_dsk))
      IF(.NOT. ALLOCATED(grot_dsk))ALLOCATE(grot_dsk(npsi_dsk))
      IF(.NOT. ALLOCATED(grotp_dsk))ALLOCATE(grotp_dsk(npsi_dsk))
      IF(.NOT. ALLOCATED(cd_dsk))ALLOCATE(cd_dsk(npsi_dsk))
      IF(.NOT. ALLOCATED(psiv_norm))ALLOCATE(psiv_norm(npsi_dsk))
      IF(.NOT. ALLOCATED(dummy))ALLOCATE(dummy(nbndry_dsk))
      IF(.NOT. ALLOCATED(xbndry))ALLOCATE(xbndry(nbndry_dsk))
      IF(.NOT. ALLOCATED(zbndry))ALLOCATE(zbndry(nbndry_dsk))


      READ (neqdata,100)(psiv_dsk(j),j=1,npsi_dsk)
      READ (neqdata,100)(p_dsk(j),j=1,npsi_dsk)
      READ (neqdata,100)(pp_dsk(j),j=1,npsi_dsk) ! actual dp/dchi
      READ (neqdata,100)(f_dsk(j),j=1,npsi_dsk)
      READ (neqdata,100)(ffp_dsk(j),j=1,npsi_dsk)
      READ (neqdata,100)(q_dsk(j),j=1,npsi_dsk)
      READ (neqdata,100)(grot_dsk(j),j=1,npsi_dsk)
      READ (neqdata,100)(grotp_dsk(j),j=1,npsi_dsk) ! actual dg/dchi
      PRINT*, 'ndbry=', nbndry_dsk, modelbnd
      IF (modelbnd.EQ.2) THEN
         nbndry = nbndry_dsk
         READ (neqdata,100) (xbndry (i), i=1,nbndry_dsk)
         READ (neqdata,100) (zbndry (i), i=1,nbndry_dsk)

      ELSE
         READ (neqdata,100) (dummy (i), i=1,nbndry_dsk)
         READ (neqdata,100) (dummy (i), i=1,nbndry_dsk)
      ENDIF

      IF (modelcd .EQ. 2) THEN
         READ (neqdata,100) (cd_dsk(j),j=1,npsi_dsk)
      ENDIF
      CLOSE(neqdata)
      CALL giveupus(neqdata)

! now we interpolate pprim and ffprim from psir grid to
! uniform psi grid(129)

       DO j=1, nj
           psir_norm(j) = (psir(j)-psir(1))/(psir(nj)-psir(1))
        ENDDO
      DO j= 1, npsi_dsk
         psiv_norm(j)= (psiv_dsk(j)-psiv_dsk(1))/ &
                               (psiv_dsk(npsi_dsk)-psiv_dsk(1))
      ENDDO
!
! --- interpolate pprim, ffprim, pressb and q onto eqdsk psi grid


      CALL intrp (1, 1, psir_norm,  pprim, nj, &
                    psiv_norm, pp_dsk, npsi_dsk)

      CALL intrp (1, 1, psir_norm, ffprim, nj, &
                     psiv_norm, ffp_dsk  , npsi_dsk)

      CALL intrp (1, 1, psir_norm,  q, nj, psiv_norm, q_dsk, npsi_dsk)

      CALL intrp (1, 1, psir_norm, fpsi, nj, &
                     psiv_norm, f_dsk  , npsi_dsk)

        CALL intrp (1, 1, psir_norm, press, nj, &
                     psiv_norm, p_dsk  , npsi_dsk)


      CALL getioun(neqdata,neqdata)
      OPEN(unit=neqdata,file='dskeqdata',status='unknown',iostat=ios)

      WRITE(neqdata,'("eqdsk conversion by eqdtoq")')
      WRITE(neqdata,'(" &end")')

      WRITE (neqdata,'(i3)') npsi_dsk
      WRITE (neqdata,'(i3)') nbndry_dsk
      WRITE (neqdata,100)(psiv_dsk(j),j=1,npsi_dsk)
      WRITE (neqdata,100)(10.*p_dsk(j),j=1,npsi_dsk)
      WRITE (neqdata,100)(pp_dsk(j)*1.e-7,j=1,npsi_dsk) ! actual dp/dchi
      WRITE (neqdata,100)(f_dsk(j)*1.e6,j=1,npsi_dsk)
      WRITE (neqdata,100)(ffp_dsk(j)*1.e4,j=1,npsi_dsk)
      WRITE (neqdata,100)(q_dsk(j),j=1,npsi_dsk)
      WRITE (neqdata,100)(grot_dsk(j),j=1,npsi_dsk)
      WRITE (neqdata,100)(grotp_dsk(j),j=1,npsi_dsk) ! actual dg/dchi

      IF (modelbnd.EQ.2) THEN
         nbndry = nbndry_dsk
         WRITE (neqdata,100) (xbndry (i), i=1,nbndry_dsk)
         WRITE (neqdata,100) (zbndry (i), i=1,nbndry_dsk)

      ELSE
         WRITE (neqdata,100) (dummy (i), i=1,nbndry_dsk)
         WRITE (neqdata,100) (dummy (i), i=1,nbndry_dsk)
      ENDIF
      IF (modelcd .EQ. 2) THEN
         WRITE (neqdata,100) (cd_dsk(j),j=1,npsi_dsk)
      ENDIF
      CLOSE(neqdata)
      CALL giveupus(neqdata)


 100  FORMAT (4e24.16)


      DEALLOCATE(psiv_dsk)  ;   DEALLOCATE(p_dsk)     ; DEALLOCATE(pp_dsk)
      DEALLOCATE(f_dsk)     ;   DEALLOCATE(ffp_dsk)   ; DEALLOCATE(q_dsk)
      DEALLOCATE(grot_dsk)  ;   DEALLOCATE(grotp_dsk) ; DEALLOCATE(cd_dsk)
      DEALLOCATE(dummy)     ;   DEALLOCATE(xbndry)    ; DEALLOCATE(zbndry)
      DEALLOCATE(psiv_norm)



      RETURN
       END




      SUBROUTINE using_toq(ic)
        USE param
        USE io, ONLY:  eqdskin,nitre
        USE ename, ONLY: eqdskfilename,eqdsk_name
        USE toq_12_interface, ONLY:toq_drive, ieqdtoq, &
        ieqdsk_toqp,ieqdsk_toq,fneqdsk
        USE soln2d
        USE solcon, ONLY : time 
        USE etc 
        IMPLICIT REAL *8 (a-h, o-z)
!        INCLUDE 'etc.i'

        CHARACTER(len = 132)command

        WRITE(nitre,FMT='("using_toq,time,dteq  =",2(2x,1pe12.6))')time,dteq
        WRITE(nitre,FMT='("ic,ieq =",3(2x,i5))')ic,ieq

!        IF(ieq .EQ. 0 .AND. ic .EQ. 0 .AND. it .EQ. 0 .AND. ieqdtoq .EQ. 1) THEN 
        IF(ieq .EQ. 0 .AND. ic .EQ. 0  .AND. ieqdtoq .EQ. 1) THEN  !04/18/05 HSJ
           WRITE(nitre,FMT='("TOQ case1,using input g eqdsk")')
            fneqdsk = '"'//eqdskin(1:LEN_TRIM(eqdskin))//'"'
!        ELSE IF(ieq .EQ. 0 .AND. ic .EQ. 0  .AND. it .EQ. 0 ) THEN 
        ELSE IF(ieq .EQ. 0 .AND. ic .EQ. 0 ) THEN  !04/18/05 HSJ
           !here ieqdtoq .ne. 1 , meaning that the startup equilibirum
           !is taken from dskeqi file (which must exist)
           WRITE(nitre,FMT='("TOQ case2,using dskeqi")')
           ieqdsk_toq = 2              ! possibly overwrites inone file value
            fneqdsk = "'not relevant'" ! Toq will read dskeqi
        ELSE IF (ieqdsk_toqp.EQ.2) THEN
           !requires that dskeqdata exists
           WRITE(nitre,FMT='("TOQ case3 using dskeqdata")')
           CALL use_eqdsk       !rewrite dskeqdata with current pprime,ffprime
                                !stop if dskeqdata is  not found
           ieqdtoq=0            !Toq will not read g file
           ieqdsk_toq = 0       !Toq will not read dskeqi
            fneqdsk = "'not relevant'"
        ELSE
           !here we require  that dskfixb exists
           WRITE(nitre,FMT='("TOQ case4,using gfile made from dskfixb")')
           CALL set_gfile_to_toq       !create g file from current 
                                       !pprim,ffprim(uses fixb code)
           ieqdtoq = 1                 !Toq will read g file
           ieqdsk_toq = 0              !Toq will not read dskeqi
           fneqdsk ="'"//'g812345.01234'//"'" 
        ENDIF


        !create Toq  namelist file with eqdsk name given by 
        !fneqdsk = eqdskfilename and spawn Toq:
        !fneqdsk may or may not be used depending on setting of
        !ieqdtoq,ieqdsk_toq,ieqdsk_toqp:
          CALL toq_drive(tocur)


        !read file dskonetwo,create input to fixb and
        !run fixb (through script fixb129xy.sc)
        !to generate eqdsk called g812345.01234
          CALL read_toq_output
          CALL eqdsk_name(eqdskfilename)
        !move the fixed file g812345.01234 created by toq to 
        !something more appropriate
       
       command = 'cp  g812345.01234 '//eqdskfilename(1:LEN_TRIM(eqdskfilename))
        IF (ISHELL(command) .LT. 0) &
              CALL STOP ('SUBROUTINE USING_TOQ', 67)



      RETURN
      END





      SUBROUTINE using_toq_in_prepar
      USE param
      USE fusion
      USE ions
      USE solcon
      USE soln
      USE mhdpar
      USE nub
      USE ename
      USE nub2
      USE     extra
      USE mesh
      USE geom
      USE psig
      USE rhog
      USE soln2d
      USE bicube
      USE etc
      USE toq_12_interface, ONLY:toq_drive, ieqdtoq, &
         ieqdsk_toqp
      IMPLICIT REAL *8 (a-h, o-z)
      INCLUDE 'imsl.i'
!      INCLUDE 'etc.i'

       CHARACTER(len =256)       :: command
       CHARACTER(len =256)       :: directory
       CHARACTER(len =256)       :: filesuf
       CHARACTER itchar*6
      itms = ABS (1000 * time)
      IF (itms .LT. 10) THEN
            WRITE (intfl, 9000) itms
            READ  (intfl, 8025) itchar
          ELSE IF (itms .LT. 100   ) THEN
            WRITE (intfl, 9020) itms
            READ  (intfl, 8025) itchar
          ELSE IF (itms .LT. 1000  ) THEN
            WRITE (intfl, 9030) itms
            READ  (intfl, 8025) itchar
          ELSE IF (itms .LT. 10000 ) THEN
            WRITE (intfl, 9040) itms
            READ  (intfl, 8025) itchar
          ELSE IF (itms .LT. 100000) THEN
            WRITE (intfl, 9050) itms
            READ  (intfl, 8025) itchar
          ELSE
            itms = itms / 10
  200       IF (itms .GE. 100000) THEN
              itms = itms / 10
              go to 200
            END IF
            WRITE (intfl, 9050) itms
            READ  (intfl, 8025) itchar
      END IF
       directory = 'dds'//itchar
      command = 'rm -rf '//directory(1:9)
        IF (ISHELL(command) .LT. 0) &
            PRINT*, " The dir is not existing, make it"
        command = 'mkdir '//directory(1:9)
        IF (ISHELL(command) .LT. 0) &
        CALL STOP ('SUBROUTINE TOQ_DRIVE: failure of mkdir dsk', 67)
!      if (ieqdsk_toqp.eq.2) then
!      call use_eqdsk
!      ieqdtoq=0
!      call toq_drive(tocur)
!      call read_toq_output
!          else
!            call set_gfile_to_toq !removed eqdskfilename from this routine
                                   !HSJ 05/03/05
!            ieqdtoq = 1
!            call toq_drive(tocur)
!      call read_toq_output
!      endif
          command = 'cp dsk* '//directory
        IF (ISHELL(command) .LT. 0) &
         CALL STOP ('SUBROUTINE TOQ_DRIVE: failure of cp dsk*', 67)


 8025 FORMAT (a)
 9000 FORMAT ('k0000' , i1)
 9020 FORMAT ('k000'  , i2)
 9030 FORMAT ('k00'   , i3)
 9040 FORMAT ('k0'    , i4)
 9050 FORMAT ('k'     , i5)


      RETURN
      END








