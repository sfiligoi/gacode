c@datavg.f
c jek 19-Jan-11 version 2.0
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c... Smoothes data via 7-pt averaging if(ismooth_all.ne.0)
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c  
      subroutine datavg
c
      implicit none
      include '../inc/data.m'
      include '../inc/input.m'
c
      integer j, jj, jn
c
       call average7_1d(fcap_d,nj_d)
       call average7_1d(fcap_d,nj_d)
       call average7_1d(gcap_d,nj_d)
       call average7_1d(gcap_d,nj_d)
       call average7_1d(hcap_d,nj_d)
       call average7_1d(hcap_d,nj_d)
c
       call average7_1d(te_d,nj_d)
       call average7_1d(ti_d,nj_d)
c
       call average7_1d(q_d,nj_d)
       call average7_1d(q_d,nj_d)
c
       call average7_1d(ene_d,nj_d)
       do jj=1,nion_d
         call average7_1d(en_d(1,jj),nj_d)
       enddo
       do jj=1,nprim_d
          call average7_1d(sion_d(1,jj),nj_d)
          call average7_1d(srecom_d(1,jj),nj_d)
          call average7_1d(scx_d(1,jj),nj_d)
          call average7_1d(sbcx_d(1,jj),nj_d)
          call average7_1d(s_d(1,jj),nj_d)
          call average7_1d(dudtsv_d(1,jj),nj_d)
       enddo
       do jn=1,nneu_d
          call average7_1d(enn_d(1,jn),nj_d)
          call average7_1d(ennw_d(1,jn),nj_d)
          call average7_1d(ennv_d(1,jn),nj_d)
c          call average7_1d(volsn_d(1,jn),nj_d)
       enddo
       call average7_1d(enbeam_d,nj_d)
       call average7_1d(enbeam_d,nj_d)
c
       call average7_1d(sbion_d,nj_d)
       call average7_1d(sbeam_d,nj_d)
c
       call average7_1d(curden_d,nj_d)
       call average7_1d(curohm_d,nj_d)
       call average7_1d(curboot_d,nj_d)
       call average7_1d(curbeam_d,nj_d)
       call average7_1d(currf_d,nj_d)
c
       call average7_1d(rbp_d,nj_d)
c
       call average7_1d(zeff_d,nj_d)
       call average7_1d(zeff_d,nj_d)
c
       call average7_1d(angrot_d,nj_d)
c
       if(ismooth_all.ge.2) call average7_1d(dpedtc_d,nj_d)
       if(ismooth_all.ge.2) call average7_1d(dpidtc_d,nj_d)
c
       call average7_1d(qdelt_d,nj_d)
       call average7_1d(qbeame_d,nj_d)
       call average7_1d(qbeami_d,nj_d)
       call average7_1d(qrfe_d,nj_d)
       call average7_1d(qrfi_d,nj_d)
       call average7_1d(qfuse_d,nj_d)
       call average7_1d(qfusi_d,nj_d)
       if(ismooth_all.ge.2) call average7_1d(qrad_d,nj_d)
       if(ismooth_all.ge.2) call average7_1d(qohm_d,nj_d)
       if(ismooth_all.ge.2) call average7_1d(qione_d,nj_d)
       if(ismooth_all.ge.2) call average7_1d(qioni_d,nj_d)
       if(ismooth_all.ge.2) call average7_1d(qcx_d,nj_d)
c
       call average7_1d(rmajavnpsi_d,nj_d)
       call average7_1d(rmajavnpsi_d,nj_d)
c
       call average7_1d(rminavnpsi_d,nj_d)
       call average7_1d(rminavnpsi_d,nj_d)
c
       call average7_1d(psivolp_d,nj_d)
       call average7_1d(psivolp_d,nj_d)
c
       call average7_1d(elongx_d,nj_d)
       call average7_1d(elongx_d,nj_d)
c
       call average7_1d(deltax_d,nj_d)
       call average7_1d(deltax_d,nj_d)
c
       call average7_1d(sfareanpsi_d,nj_d)
       call average7_1d(sfareanpsi_d,nj_d)
c
       call average7_1d(cxareanpsi_d,nj_d)
       call average7_1d(cxareanpsi_d,nj_d)
c      
       call average7_1d(grho1npsi_d,nj_d)
       call average7_1d(grho1npsi_d,nj_d)
c
       call average7_1d(grho2npsi_d,nj_d)
       call average7_1d(grho2npsi_d,nj_d)
c
       if(itorque.ne.0 .or. iptotr.eq.2)then
         call average7_1d(torque_d,nj_d)
       endif
c
       call average7_1d(ptot_d,nj_d)
       if(iptotr.eq.2) then
         call average7_1d(pfast_d,nj_d)
       endif
c
       if(ncl_flag.eq.1) then
         call average7_1d(xb2_d,nj_d)
         call average7_1d(xbm2_d,nj_d)
         call average7_1d(xb2_d,nj_d)
         call average7_1d(xbm2_d,nj_d)
         call average7_1d(xngrth_d,nj_d)
         call average7_1d(xgrbm2_d,nj_d)
         call average7_1d(fm1_d,nj_d)
         call average7_1d(fm2_d,nj_d)
         call average7_1d(fm3_d,nj_d)
         call average7_1d(fhat_d,nj_d)
       endif
c
      end
