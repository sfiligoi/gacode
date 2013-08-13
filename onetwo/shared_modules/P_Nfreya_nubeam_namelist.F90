 
  MODULE P_Nfreya_nubeam_namelist
 !--------------------------------------------------------------------------
 ! THIS NAMELIST MODULE IS USED FOR Nubeam  input.
 ! IT is based on the new nubeam version 201112 namelits input
 ! and can not be used with the older nubeam (eg nubeam prior to
 ! 201112)
 !------------------------------------------------------------HSJ-2/14//2012
       USE nrtype,                                 ONLY : Dp,I4B,SP
       USE param,                                  ONLY : kj,kprim,kimp,kb
       USE dnubeam_mod



!       USE P_nfreya_interface,   ONLY :                           &
!                           beam_data => beam_datao  
!                           nubeam_nclass,nlfhe3,nlfhe4,   &
!                         nlfst,nlfsp, nlusf3,xdatsfa, xdatsf3,    &
!                         xdatsft,xdatsfp,nlusfa,nlusfp,nlusft,    &
!                         plfhe3,plfhe4,plfst,plfsp,dxbsmoo,wghta, &
!                         nptcls,nptclf,ndep0,goocon,dtn_orbit,    &
!                         nsigexc,nlminsv,nlbbcx,nlebei,dn0out,    &
!                         nmsigx,xcfanbi,xdfanbi,xefanbi,xcfafus,  &
!                         xdfafus,xefafus,nlbcde,nlbcoh,nlbcpa,    &
!                         nlorbo,nlbflr,nlfbmflr,nlbgflr,gflr_min, &
!                         nkdifb,fdifbe,edifbe,plhgt,nlcprb,       &
!                         nmcurb,xdepmod,xp_nt1,inzri,             &
!                         ndifbep,nubeam_restart,                  &
!                         nubeam_dt,nubeam0_dt,d_fast_ion

! 
       USE numbrs,    ONLY   : nj,nprim,nimpc=>nimp
!
!
       IMPLICIT NONE
       SAVE
!        INTEGER kj_check
!       INTEGER,PARAMETER   :: nbeamx     = kb ! HSJ 8/10/11
        INTEGER,PARAMETER   :: nbeamx     = nbmax ! use nbmax 2/15//2012 HSJ
!       INTEGER,PARAMETER   :: nrhixm     = 4   !set to beam_data%nrhix ??
!       INTEGER,PARAMETER   :: nfracionsx = 1 ! dummy
!       INTEGER,PARAMETER   :: nfracimpsx = 1 !dummy 
!
!      nfracionsx,nfracimpsx and frac_ions, frac_imp are not
!      used to generate profiles in the Onetwo version.
!      thes e quantities are retained only so that the
!      namelist reads will work without having to modify the
!      namelists. See ...../nbdrive/nbdrive_prof.f90
!      on how these quantities are eliminated 
!
!       INTEGER,PARAMETER :: nneutx = 10
!       INTEGER,PARAMETER :: ngmaxx = 10
        REAL(DP) tbona(nbeamx),tboffa(nbeamx) !required for sub nblist_read only
!
!
  ! the PYTHON generated namelist values have a one-to-one correspondence
  ! with elements of NUBEAM input data structures -- units and meaning are
  ! documented in nubeam/nbspec.dat (** not repeated here **).
!
  ! only a subset of the NUBEAM controls are exposed in this NUBEAM
  ! test driver namelist.  This is an effort at simplification.  The
  ! omitted controls correspond to specialized "code developer"
  ! adjustments, or, to "advanced" features which (for the time
  ! being at least) have been omitted from the test driver.
!
  ! inputs needed specifically by the driver code (but not directly by
  ! NUBEAM itself), without this one-to-one correspondence, are declared
  ! and described in the "hand maintained section", below...


    REAL(DP) :: ebdmax
  !****NAMELIST/nbdrive_naml/ nzones,nznbma,nznbme,ebdmax
 
  ! --> nubeam/nbspec.dat block: beams
  !  INTEGER  :: ngmax
  !  REAL(DP) :: backz(ngmaxx)  ! (1:ngmax)
  !   REAL(DP) :: aplasm(ngmaxx)  ! (1:ngmax)
 



  ! added to this code by hsj 5/9/11 nlbdat 
  ! this appears to be a new variable in namelist which is written by
  ! C. Greenfield's IDL code )
   LOGICAL nlbdat
   REAL(DP)  :: nubeam_dt

! use the namelist from the new nubeam
! dnubeam_driver.f90 contains this namelist as well
  
  include './shared_modules/nbnamelist_post_201112.inc'

  END MODULE P_Nfreya_nubeam_namelist
