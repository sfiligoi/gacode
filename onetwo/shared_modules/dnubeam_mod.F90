
  MODULE dnubeam_mod
!-----------------------------------------------------------------------------------------
! -- ORIGNAL FROM JM PARK, REWORKED FOR COMPATIBILITY WITH ONETWO/P_Nfreya, HSJ 2/15/2012
!-----------------------------------------------------------------------------------------
  USE nrtype,                 ONLY : DP,I4B
  USE param,                  ONLY : kb
#ifdef ONETWO
  USE io,                     ONLY : lun_nubeam
#endif
#ifdef P_Nfreya
  USE io_gcnmp,               ONLY : lun_nubeam                
#endif
!  INTEGER(I4B), PARAMETER :: nbmax = 32  ! HSJ
  INTEGER(I4B),  PARAMETER :: nbmax = kb  ! HSJ
  INTEGER(I4B),  PARAMETER :: nionmax = 10
  INTEGER(I4B),  PARAMETER :: nrhomax = 201
  INTEGER(I4B),  PARAMETER :: ntimemax = 30 
!  INTEGER(I4B),  PARAMETER :: lun = 99 HSJ

!---nbi configuration 

  INTEGER(I4B) nbeam
  LOGICAL, DIMENSION(nbmax)      :: nlco
  REAL(DP),  DIMENSION(nbmax)    :: abeama,xzbeama
  INTEGER(I4B), DIMENSION(nbmax) :: nbshapa
  REAL(DP),  DIMENSION(nbmax)    :: rtcena,xlbtna,xybsca
  REAL(DP),  DIMENSION(nbmax)    :: bmwidra,bmwidza
  REAL(DP),  DIMENSION(nbmax)    :: divra,divza,foclza,foclra
  INTEGER(I4B), DIMENSION(nbmax) :: nbapsha
  REAL(DP),  DIMENSION(nbmax)    :: rapedga,xzpedga,xlbapa,xybapa
  REAL(DP),  DIMENSION(nbmax)    :: xrapoffa,xzapoffa
  INTEGER(I4B), DIMENSION(nbmax) :: nbapsh2
  REAL(DP),  DIMENSION(nbmax)    :: rapedg2,xzpedg2,xlbapa2
  REAL(DP),  DIMENSION(nbmax)    :: xrapoff2,xzapoff2
  REAL(DP),  DIMENSION(nbmax)    :: xbzeta
  INTEGER(I4B), DIMENSION(nbmax) :: ntrace

!---nbi power, energy, species mix

  REAL(DP),  DIMENSION(nbmax) :: pinja,einja,ffulla,fhalfa

!---nubeam input profile

  REAL(DP) :: t0_nubeam, t1_nubeam          ! sec
 
  INTEGER :: nrho_in                      ! number of input grid

  INTEGER(I4B) :: nion_in                      ! number of thermal ions (including inmpurities)
  REAL(DP), DIMENSION(nionmax) :: aion_in   ! A
  REAL(DP), DIMENSION(nionmax) :: zion_in   ! Z
 
  REAL(DP), DIMENSION(nrhomax) :: rho_in    ! input grid [normalized]
  REAL(DP), DIMENSION(nrhomax) :: te_in     ! [keV]
  REAL(DP), DIMENSION(nrhomax) :: ti_in     ! [keV]
  REAL(DP), DIMENSION(nrhomax) :: ne_in     ! [m^-3]
  REAL(DP), DIMENSION(nrhomax) :: zeff_in   ! []
  REAL(DP), DIMENSION(nrhomax) :: omega_in  ! [rad/sec]
 
  REAL(DP), DIMENSION(nrhomax) :: vloop_in  ! [V]
  REAL(DP), DIMENSION(nrhomax) :: n0w_in    ! [m^-3]
  REAL(DP), DIMENSION(nrhomax) :: n0v_in    ! [m^-3]
  REAL(DP), DIMENSION(nrhomax) :: t0_in     ! [kev]
  REAL(DP)                     :: f0w_in    ! [#/sec]
 
!---nubeam output profile

  INTEGER(I4B) :: nrho_out

  REAL(DP), DIMENSION(nrhomax) :: rho_out      ! output grid [normalized]
  REAL(DP), DIMENSION(nrhomax) :: pbe_out      ! W/m^3
  REAL(DP), DIMENSION(nrhomax) :: pbi_out      ! W/m^3
  REAL(DP), DIMENSION(nrhomax) :: curb_out     ! A/m^2
  REAL(DP), DIMENSION(nrhomax) :: ucurb_out    ! A/m^2
  REAL(DP), DIMENSION(nrhomax) :: curbdotb_out ! A/m^2/T
  REAL(DP), DIMENSION(nrhomax) :: bdep_out     ! #/m^3/sec
  REAL(DP), DIMENSION(nrhomax) :: bdenss_out   ! #/m^3
  REAL(DP), DIMENSION(nrhomax) :: pbth_out     ! #/m^3
  REAL(DP), DIMENSION(nrhomax) :: tqbsum_out   ! NT-M/M**3 
  REAL(DP), DIMENSION(nrhomax) :: tqbi_out     ! NT-M/M**3 
  REAL(DP), DIMENSION(nrhomax) :: tqbe_out     ! NT-M/M**3 
  REAL(DP), DIMENSION(nrhomax) :: udenspl_out  ! J/M**3
  REAL(DP), DIMENSION(nrhomax) :: udenspp_out  ! J/M**3

!---nubeam output scalar

  REAL(DP)  pbe_tot
  REAL(DP)  pbi_tot
  REAL(DP)  pbth_tot
  REAL(DP)  curb_tot
  REAL(DP)  ucurb_tot
  REAL(DP)  bdep_tot
  REAL(DP)  bdenss_tot
  REAL(DP)  tqbsum_tot,tqbi_tot,tqbe_tot
  REAL(DP)  udenspl_tot,udenspp_tot,wbeam_tot
  REAL(DP)  bbntot_out
  REAL(DP)  btneut_out

!---internal data

  REAL(DP), DIMENSION(nrhomax,2) :: rho_save 
  REAL(DP), DIMENSION(nrhomax,2) :: pbe_save
  REAL(DP), DIMENSION(nrhomax,2) :: pbi_save
  REAL(DP), DIMENSION(nrhomax,2) :: curb_save
  REAL(DP), DIMENSION(nrhomax,2) :: ucurb_save
  REAL(DP), DIMENSION(nrhomax,2) :: curbdotb_save
  REAL(DP), DIMENSION(nrhomax,2) :: bdep_save
  REAL(DP), DIMENSION(nrhomax,2) :: bdenss_save
  REAL(DP), DIMENSION(nrhomax,2) :: pbth_save
  REAL(DP), DIMENSION(nrhomax,2) :: tqbsum_save
  REAL(DP), DIMENSION(nrhomax,2) :: tqbi_save
  REAL(DP), DIMENSION(nrhomax,2) :: tqbe_save
  REAL(DP), DIMENSION(nrhomax,2) :: udenspl_save
  REAL(DP), DIMENSION(nrhomax,2) :: udenspp_save

  REAL(DP), DIMENSION(2) :: pbe_tot_save
  REAL(DP), DIMENSION(2) :: pbi_tot_save
  REAL(DP), DIMENSION(2) :: pbth_tot_save
  REAL(DP), DIMENSION(2) :: curb_tot_save
  REAL(DP), DIMENSION(2) :: ucurb_tot_save
  REAL(DP), DIMENSION(2) :: bdep_tot_save
  REAL(DP), DIMENSION(2) :: bdenss_tot_save
  REAL(DP), DIMENSION(2) :: tqbsum_tot_save,tqbi_tot_save,tqbe_tot_save
  REAL(DP), DIMENSION(2) :: udenspl_tot_save,udenspp_tot_save
  REAL(DP), DIMENSION(2) :: bbntot_save
  REAL(DP), DIMENSION(2) :: btneut_save

!--control

  REAL(DP) :: dn0out  ! now state-file input
  INTEGER:: nkdifb  ! now staae-file input

  INTEGER(I4B) :: adiff_ntime
  REAL(DP), DIMENSION(ntimemax)  :: adiff_a,adiff_0
  REAL(DP), DIMENSION(ntimemax)  :: adiff_xpin,adiff_xpout,adiff_time
  REAL(DP) :: difb_0,difb_a,difb_in,difb_out

  !character*120 :: mission,workpath
  REAL(DP) :: orbrzv_toric_prftot0(15),orbrzv_toric_frnm
  REAL(DP) :: xdatsf3,xdatsfa,xdatsft,xdatsfp
  REAL(DP) :: plfhe3,plfhe4,plfst,plfsp,cxsplt,dxbsmoo,dxbsm_nc,xpand_nptcl
  REAL(DP) :: frac_depmax,frac_dep_lim,frac_depmin,frac_orbrr
  REAL(DP) :: aimp,xzimp
  REAL(DP) :: wghta,goocon,cxpcon,fppcon,fdtnxy,dtn,dt_acc,orbrzv_zzerr_con
  REAL(DP) :: xdepmod,xcfanbi,xdfanbi,xefanbi,xcfafus,xdfafus,xefafus
  REAL(DP) :: gflr_min,gflr_rl,gflr_ll,gflr_xv(2),gflr_op(15)
  REAL(DP) :: xswfrac_allfast,xswfrac_beam,xswfrac_fusn,fporcelli
  REAL(DP) :: taurip,asrd,bsrd
  REAL(DP) :: erngfi(10),fbemin,fbemax,fvpvmn,fvpvmx,fbltim,fshper,fshwid
  REAL(DP) :: tfshon,tfshof,xfishmin,xfishmax,fbtrap_depth
  REAL(DP) :: xboxhw,yboxhw,xlbox1,xlbox2
  REAL(DP) :: rtube(200),ytube(200),xzetatube(200)
  REAL(DP) :: phitube(200),thetatube(200),xl1tube(200),xl2tube(200)
  REAL(DP) :: rhotube,edbfac

  INTEGER(I4B) :: nonlin,nseed,nlbout,lunnbx,lunres,mrstrt,nclass
  INTEGER(I4B) :: only_io,ref_namelist,quasi_check,nltest_output
  INTEGER(I4B) :: blk_mpi_bcast_int,blk_mpi_bcast_r8
  INTEGER(I4B) :: nzones,nzone_fb,nlsym2b,nth0,xplasma_in_memory
  INTEGER(I4B) :: nznbma,nznbme,ngyro
 ! INTEGER(I4B) :: nlusf3,nlusfa,nlusft,nlusfp,nlfatom,nlbfpp
  LOGICAL    nlusf3,nlusfa,nlusft,nlusfp,nlfatom,nlbfpp
  INTEGER(I4B) :: nptcls,nptclf,nptclh,nptcl_max
  INTEGER(I4B) :: nltrk_dep0,nldep0_gather,ndep0,ndep0_max,ndep_set_beam,ndep_set_ien
  INTEGER(I4B) :: nbbcal,ncx0,nper_cx,lev_nbidep,levmod_halo,nmimp
  INTEGER(I4B) :: ngoocon_vpvbin,ngoocon_ebin,ngoocon_rbin,ndtorb,nchdvp
  INTEGER(I4B) :: orbrzv_option,orbrzv_rf_option
  INTEGER(I4B) :: nlminsv,nlebei,nsigexc,nlbbcx,nbbcx_avg,nbbcx_bb,nbfallgr
!  INTEGER(I4B) :: nlbeamcx,nlhvion,nlbcde,nlbcoh,nlbcpa,nlorbo,nlbflr,nlfbmflr,nlbgflr
  LOGICAL nlbflr,nlbgflr,nlfbmflr,nlorbo,nlbcpa,nlbcoh  ! HSJ 3/1/2013
  INTEGER(I4B) ::  nlbeamcx,nlhvion,nlbcde
  INTEGER(I4B) :: sawflag,nlsawb,nlsawf,nmix_kdsaw,ngradd_opt
  INTEGER(I4B) :: nrip,nerngfi,nsdbgb,nlfbon,nfbon_vpvopt,nfbon_species
  INTEGER(I4B) :: nbbox,nbsbox(100),nbebox(100)
  INTEGER(I4B) :: nxbox,nybox,nlbox,ndepbox
  INTEGER(I4B) :: nbtube,nbstube(200),nbetube(200),lmidtube(200),nsegtube(200)
!  INTEGER(I4B) :: ndeptube,nlcprb,nlpsirz,nmcurb,nlfdep
    INTEGER(I4B) :: ndeptube,nlpsirz,nmcurb,nlfdep
  LOGICAL nlcprb
END MODULE dnubeam_mod
