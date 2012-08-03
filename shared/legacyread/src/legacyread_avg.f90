! jek 19-Jan-11 version 2.0
! jmc 20-Jul-12 mapped to f90
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!
!... Smoothes data via 7-pt averaging if(ismooth.ne.0)
!
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!  
subroutine datavg(ismooth,ncl_flag)

  use legacyread_interface

  implicit none

  integer :: ismooth,ncl_flag
  integer :: j,jj,jn

  call average7_1d(fcap_d,nj_d)
  call average7_1d(fcap_d,nj_d)
  call average7_1d(gcap_d,nj_d)
  call average7_1d(gcap_d,nj_d)
  call average7_1d(hcap_d,nj_d)
  call average7_1d(hcap_d,nj_d)
  
  call average7_1d(te_d,nj_d)
  call average7_1d(ti_d,nj_d)
  
  call average7_1d(q_d,nj_d)
  call average7_1d(q_d,nj_d)
  
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
  enddo
  call average7_1d(enbeam_d,nj_d)
  call average7_1d(enbeam_d,nj_d)
  
  call average7_1d(sbion_d,nj_d)
  call average7_1d(sbeam_d,nj_d)
  
  call average7_1d(curden_d,nj_d)
  call average7_1d(curohm_d,nj_d)
  call average7_1d(curboot_d,nj_d)
  call average7_1d(curbeam_d,nj_d)
  call average7_1d(currf_d,nj_d)
  
  call average7_1d(rbp_d,nj_d)
  
  call average7_1d(zeff_d,nj_d)
  call average7_1d(zeff_d,nj_d)
  
  call average7_1d(angrot_d,nj_d)
  
  if(ismooth.ge.2) call average7_1d(dpedtc_d,nj_d)
  if(ismooth.ge.2) call average7_1d(dpidtc_d,nj_d)
  
  call average7_1d(qdelt_d,nj_d)
  call average7_1d(qbeame_d,nj_d)
  call average7_1d(qbeami_d,nj_d)
  call average7_1d(qrfe_d,nj_d)
  call average7_1d(qrfi_d,nj_d)
  call average7_1d(qfuse_d,nj_d)
  call average7_1d(qfusi_d,nj_d)
  if(ismooth.ge.2) call average7_1d(qrad_d,nj_d)
  if(ismooth.ge.2) call average7_1d(qohm_d,nj_d)
  if(ismooth.ge.2) call average7_1d(qione_d,nj_d)
  if(ismooth.ge.2) call average7_1d(qioni_d,nj_d)
  if(ismooth.ge.2) call average7_1d(qcx_d,nj_d)
  
  call average7_1d(rmajavnpsi_d,nj_d)
  call average7_1d(rmajavnpsi_d,nj_d)
  
  call average7_1d(rminavnpsi_d,nj_d)
  call average7_1d(rminavnpsi_d,nj_d)
  
  call average7_1d(psivolp_d,nj_d)
  call average7_1d(psivolp_d,nj_d)
  
  call average7_1d(elongx_d,nj_d)
  call average7_1d(elongx_d,nj_d)
  
  call average7_1d(deltax_d,nj_d)
  call average7_1d(deltax_d,nj_d)
  
  call average7_1d(sfareanpsi_d,nj_d)
  call average7_1d(sfareanpsi_d,nj_d)
  
  call average7_1d(cxareanpsi_d,nj_d)
  call average7_1d(cxareanpsi_d,nj_d)
       
  call average7_1d(grho1npsi_d,nj_d)
  call average7_1d(grho1npsi_d,nj_d)
  
  call average7_1d(grho2npsi_d,nj_d)
  call average7_1d(grho2npsi_d,nj_d)
  
  call average7_1d(torque_d,nj_d)
  
  call average7_1d(ptot_d,nj_d)
  
  call average7_1d(pfast_d,nj_d)
  
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
  
end subroutine datavg

subroutine average7_1d(f,nk)

  !*****************************************************
  !
  ! performs a two-pass seven-point average 
  !
  !*****************************************************

  implicit none

  integer :: nk,k,i,m
  real :: f(nk),g(nk)

  ! check if grid is too small
  if (nk < 7) return

  g(1)=f(1)
  g(nk)=f(nk)
  k=2
  g(k)=(f(k)+f(k+1)+f(k-1))/3.0
  k=3
  g(k)=(f(k)+f(k+1)+f(k-1)+f(k+2)+f(k-2))/5.0
  k=nk-1
  g(k)=(f(k)+f(k+1)+f(k-1))/3.0
  k=nk-2
  g(k)=(f(k)+f(k+1)+f(k-1)+f(k+2)+f(k-2))/5.0
  do k=4,nk-3
     g(k)=(f(k)+f(k-1)+f(k+1)+f(k-2)+f(k+2)+f(k+3)+f(k-3))/7.0
  enddo
  k=2
  f(k)=(g(k)+g(k+1)+g(k-1))/3.0
  k=3
  f(k)=(g(k)+g(k+1)+g(k-1)+g(k+2)+g(k-2))/5.0
  k=nk-1
  f(k)=(g(k)+g(k+1)+g(k-1))/3.0
  k=nk-2
  f(k)=(g(k)+g(k+1)+g(k-1)+g(k+2)+g(k-2))/5.0
  do k=4,nk-3
     f(k)= (g(k)+g(k+1)+g(k-1)+g(k+2)+g(k-2)+g(k+3)+g(k-3))/7.0
  enddo

end subroutine average7_1d
