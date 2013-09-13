      MODULE TDEM
!
!
      USE param, only : kjjjj => kj
      USE mhdpar, only : mxtbcmhd, nw, nh 
      USE contour ,only : nconmax
      USE limiter ,only: maxlimpt
      implicit none
     
!
!
      save
!
!     stores quantities related to time-dependent eqdsk mode
!     of operation, has to be used in conjunction with param.i
!     and mhdpar.i ... HSJ 6/13/96
!
!
!
      integer, parameter   ::  miscelan = 50
      logical,public       ::  ncdfile_open
      integer *4 ,  public,dimension(mxtbcmhd)  :: shot_tdem, time_tdem
      integer *4 ,  public,dimension(mxtbcmhd)  :: ncontr_tdem 
      integer *4 ,  public,dimension(miscelan )  :: intg_tdem
      integer *4,   public :: nlim_tdem, ntime_tdem, &
                              il_tdem, ncid, tdemvb,iostat
      data         ncdfile_open, iostat /.false., 0/
      real*8 ,  public, dimension(:) ::                           &
                   xlim_tdem(maxlimpt), ylim_tdem(maxlimpt),      &
                   real_tdem(miscelan),   vol_tdem(mxtbcmhd),     &
                   circum_tdem(mxtbcmhd), rma_tdem(mxtbcmhd),     &
                   zma_tdem(mxtbcmhd),    psilim_tdem(mxtbcmhd),  &
                   psimag_tdem(mxtbcmhd), rsep_tdem(mxtbcmhd),    &
                   zsep_tdem(mxtbcmhd),   toteqd_tdem(mxtbcmhd),  &
                   beqd_tdem(mxtbcmhd),   btorax_tdem(mxtbcmhd),  &
                   psisep_tdem(mxtbcmhd), area_tdem(mxtbcmhd),    &
                   fpsi_tdem(nw),                                 &
                   qpsi_tdem(nw),         ffprimpsi_tdem(nw),     &
                   pprimpsi_tdem(nw),     presspsi_tdem(nw),      &
                   rcontr_tdem(nconmax),  zcontr_tdem(nconmax),   &
                   rtime_tdem(mxtbcmhd),  curden_tdem(kjjjj),        &
                   hcap_tdem(kjjjj),         r2cap_tdem(kjjjj),         &
                   rhoa_tdem(mxtbcmhd),   q_tdem(kjjjj),             &
                   bp0_tdem(kjjjj),          volpsi_tdem(nw),        &
                   psigrid_tdem(nw),      dbdt_tdem(kjjjj),          &
                   drbpdt_tdem(kjjjj),                               &
                   bpo_save_tdem(kjjjj),  dpsidt_const_rho_tdem(kjjjj), &
                   rhomax_smooth_tdem(mxtbcmhd),                  &
                   dpsidrho_const_t_tdem(kjjjj),                     &
                   dpsidt_const_zeta_tdem(kjjjj),                    &
                   const_zeta_psigrid_smooth_tdem(kjjjj)    
      
      real*8 ,  public, dimension(:,:) ::  psi_tdem(nw,nh),       &
                                     csp_rho_fit_tdem(mxtbcmhd,3)
      real*8 ,  public :: dpsidt_tdem,psimag_t,psilim_t,rma_t,zma_t
   
!
 
   END MODULE TDEM
