      subroutine write_ascii(lprint,lprint_pulse)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c     This routine writes out results in ASCII format
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      implicit none
c
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/ptor.m'
      include '../inc/glf.m'
c
      integer j, lprint, lprint_pulse
c
      if (lprint.eq. 3) then
        open (34,file='iterdbprofs.dat',status='unknown')
        write(34,'(a12)')' Te - iterdb'
        write(34,500) (te_exp(j),j=0,nj_d-1)
        write(34,'(a12)')' Ti - iterdb'
        write(34,500) (ti_exp(j),j=0,nj_d-1)
        write(34,'(a12)')' ne - iterdb'
        write(34,500) (ne_exp(j),j=0,nj_d-1)
        write(34,'(a12)')' ni - iterdb'
        write(34,500) (ni_exp(j),j=0,nj_d-1)
c
        write(34,'(a12)')' Te - xptor'
        write(34,500) (te_m(j),j=0,nj_d-1)
        write(34,'(a12)')' Ti - xptor'
        write(34,500) (ti_m(j),j=0,nj_d-1)
        write(34,'(a12)')' ne - xptor'
        write(34,500) (ne_m(j),j=0,nj_d-1)
        write(34,'(a12))')' ni - xptor'
        write(34,500) (ni_m(j),j=0,nj_d-1)
      endif
c
      if (lprint_pulse .eq. 2) then
        open (9,file='tepulse.dat',status='unknown')
        open (12,file='tipulse.dat',status='unknown')
      else
        open (7,file='out',status='unknown')
        open (8,file='t0.dat',status='unknown')
        open (9,file='nepulse.dat',status='unknown')
        open (10,file='tepulse.dat',status='unknown')
        open (12,file='tipulse.dat',status='unknown')
        open (15,file='vexbpulse.dat',status='unknown')
        open (31,file='vpolpulse.dat',status='unknown')
        open (16,file='exb.dat',status='unknown')
        open (17,file='vexb.dat',status='unknown')
        open (18,file='vetor.dat',status='unknown')
        open (19,file='vepol.dat',status='unknown')
        open (20,file='vstar.dat',status='unknown')
        open (21,file='egamma.dat',status='unknown')
        open (22,file='anrate.dat',status='unknown')
        open (23,file='etaphi.dat',status='unknown')
        open (24,file='gammap.dat',status='unknown')
        open (27,file='egamma-vphi.dat',status='unknown')
        open (28,file='egamma-vpol.dat',status='unknown')
        open (29,file='egamma-vstar.dat',status='unknown')
        open (32,file='chie.dat',status='unknown')
        open (33,file='chii.dat',status='unknown')
c
        write(7,150) mxgrid
        do j=1,mxgrid
          write(7,200) j, rho(j), te_exp_sav(j), 
     &                 ti_exp_sav(j), te_m(j), ti_m(j),ne_exp(j),ne_m(j)
        enddo
        do j=1,mxgrid
          write(7,210) j, rho(j), vphi_exp(j), 
     &                 angrotp_exp_sav(j)*rmajor_exp, 
     &               angrot_exp(j), angrotp_exp(j), vphi_m(j)
        enddo
        do j=1,mxgrid
          write(7,210) j, rho(j), zpti_exp(j), rhosda_exp(j),
     &         csda_exp(j), drhodr(j)
        enddo
        do j=1,mxgrid
          write(7,205) j, rho(j), chiegb_m(j)*cgyrobohm_m(j),
     &       chiigb_m(j)*cgyrobohm_m(j), 
     &       etagb_phi_m(j)*cgyrobohm_m(j),
     &       chieneogb_m(j)*cgyrobohm_m(j),
     &       chiineogb_m(j)*cgyrobohm_m(j)
        enddo
        write(16,150) mxgrid
        do j=1,mxgrid
c          write(16,215) j, rho(j), vexb_m(j), vstar_m(j),
c     &       vepol_m(j), vetor_m(j), anrate_m(j)
          write(16,215) j, rho(j), egamma_m(j), egamma_vphi(j),
     &       egamma_vpol(j), egamma_vstar(j), anrate_m(j)
        enddo
c
        write(30,150) jout_m
        do j=0,jout_m
          write(30,202) j, rho(j), vphi_exp(j),vphi_m(j)
        enddo
c
        write(8,150) nsteps_v+1
        do j=0,nsteps_v
          write(8,250) j, te_t(0,j), ti_t(0,j)
        enddo
c
      endif
c
      if (lprint.eq. 2 .and. lprint_pulse.ne.2) then
        open (25,file='out2',status='unknown')
        write(25,150) mxgrid
        do j=1,mxgrid
          write(25,200) j, rho(j), ne_exp(j), ni_exp(j),
     &       nz_exp(j), nfst_exp(j), zeff_exp(j), q_exp(j)
        enddo
        close(25)
      endif
c
      if (lprint_pulse .eq. 1 .and. dvflag .gt. 0) then
        write(9,150) ntime_t, ngrid
        write(10,150) ntime_t, ngrid
        write(12,150) ntime_t, ngrid
        write(15,150) ntime_t, ngrid
        write(31,150) ntime_t, ngrid
        write(32,150) ntime_t, ngrid
        write(33,150) ntime_t, ngrid
        do j=0,ntime_t
          write(9,310) time_t(j), ne_t(0:ngrid,j)
          write(10,310) time_t(j), te_t(0:ngrid,j)
          write(12,310) time_t(j), ti_t(0:ngrid,j)
          write(15,310) time_t(j), vexb_t(0:ngrid,j)
          write(31,310) time_t(j), vpol_t(0:ngrid,j)
          write(32,310) time_t(j), chie_t(0:ngrid,j)
          write(33,310) time_t(j), chii_t(0:ngrid,j)
        enddo
      elseif (lprint_pulse .ge. 2) then
        write(10,150) nsteps_v+1
        write(12,150) nsteps_v+1
        write(10,320) rho(0:ngrid)
        write(12,320) rho(0:ngrid)
        do j=0,nsteps_v
          write(10,315) time_t(j), te_t(0:ngrid,j)
          write(12,315) time_t(j), ti_t(0:ngrid,j)
        enddo
      else
        write(9,150) nsteps_v+1
        write(10,150) nsteps_v+1
        write(12,150) nsteps_v+1
        write(15,150) nsteps_v+1
        write(17,150) nsteps_v+1
        write(18,150) nsteps_v+1
        write(19,150) nsteps_v+1
        write(20,150) nsteps_v+1
        write(21,150) nsteps_v+1
        write(22,150) nsteps_v+1
        write(23,150) nsteps_v+1
        write(24,150) nsteps_v+1
        write(27,150) nsteps_v+1
        write(28,150) nsteps_v+1
        write(29,150) nsteps_v+1
        write(32,150) nsteps_v+1
        write(33,150) nsteps_v+1
        do j=0,nsteps_v
          write(9,300) j, ne_t(0:ngrid,j)
          write(10,300) j, te_t(0:ngrid,j)
          write(12,300) j, ti_t(0:ngrid,j)
          write(15,350) j, vphi_t(0:ngrid,j)
          write(17,350) j, vexb_t(0:ngrid,j)
          write(18,350) j, vetor_t(0:ngrid,j)
          write(19,350) j, vepol_t(0:ngrid,j)
          write(20,350) j, vstar_t(0:ngrid,j)
          write(21,350) j, egamma_t(0:ngrid,j)
          write(22,350) j, anratem_t(0:ngrid,j)
          write(23,350) j, etaphi_t(0:ngrid,j)
          write(24,350) j, gammap_t(0:ngrid,j)
          write(27,350) j, egamma_vphi_t(0:ngrid,j)
          write(28,350) j, egamma_vpol_t(0:ngrid,j)
          write(29,350) j, egamma_vstar_t(0:ngrid,j)
          write(32,350) j, chie_t(0:ngrid,j)
          write(33,350) j, chii_t(0:ngrid,j)
        enddo
      endif
c
 150  format(2i4)
 200  format(i2,2x,0p1f4.2,0p6f10.5)
 202  format(i2,2x,0p1f4.2,0p5f11.2)
 205  format(i2,2x,0p1f4.2,0p6f12.5)
 210  format(i2,2x,0p1f4.2,0p5f14.5)
 215  format(i2,2x,0p1f4.2,0p4f14.5,2x,1pf10.5)
 250  format(i4,2x,0p4f10.5)
 300  format(i4,2x,0p301f10.5)
 310  format(0pe16.7,2x,0p301e20.11)
 315  format(0pf10.5,2x,0p301f12.5)
 320  format(3x,'time',3x,0p301f12.5)
 350  format(i4,2x,0p301f14.5)
 500  format(5(2x,1pe14.4))     ! iterdb format
c
      close(7)
      close(8)
      close(9)
      close(10)
      close(12)
      close(15)
      close(16)
      close(17)
      close(18)
      close(19)
      close(20)
      close(21)
      close(22)
      close(23)
      close(24)
      close(27)
      close(28)
      close(29)
      close(31)
      close(32)
      close(33)
      close(34)
c
      return
      end
