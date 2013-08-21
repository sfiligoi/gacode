

  subroutine set_nbprofs
!_________________________________________________
! store profiles as MKS; temperatures in KeV
! convert from TRANSP units
!.................................................

  USE transp
  use nbi_types
  USE ions
  USE neut
  USE soln
  use nbi_dimensions ! get mj from here
  implicit NONE



  type (nbitype_profiles) :: zprf

  ! the xplasma ids returned for each profile fit are stored
  ! as elements of "zprf" which can then be passed conveniently as
  ! a block to the NUBEAM module...

  integer, parameter :: ispline   = 2  ! spline fit
  integer, parameter :: ihermite  = 1  ! hermite fit
  integer, parameter :: ig00      = 1  ! f'=0 at axis, not-a-knot @edge
  integer, parameter :: iknot     = 0  ! not-a-knot BC, axis & edge
  integer, parameter :: izcen     = 2  ! zone centered quantity
  integer, parameter :: izbdy     = 1  ! bdy centered quanity
  real zwork(mj)

  integer            :: ii, ig, is, itype, isb, str_length
  integer            :: iz, ia, inst, ierr
  character(len=3)   :: chid
  character(len=10)  :: zid

  integer index_nbfi     ! indexing function

  integer j
  real *4 one,onemill

  !-------------------------

  ierr=0
  one = 1.0
  !note:  xpload_12  will load boundary oriented data using
  !id_rhonj
  !zone centered data will be loaded using id_rhozc
  !both are passed by xpload_module
  

   call nbi_init_profiles(zprf)  !zero the ids

!-------------------------------------------------------------------
!temperatures:  Te, Ti, Tmajority (=Ti if defaulted)
!  nubeam  wants these in KEV in r4 form:
!--------------------------------------------------------------------
     print *,'in set_nbprofs, calling xpload_12'
     !register te profile,zprf%id_te is determined by xpload_12:
     call xpload_12(te,one,ispline,ig00,izbdy,'TE',zprf%id_te)
     !register ti profile,zprf%id_ti is determined by xpload_12:
    call xpload_12(ti,one,ispline,ig00,izbdy,'TI',zprf%id_ti)
    print *,'zprf%id_te,zprf%id_ti =',zprf%id_te,zprf%id_ti
   if(zprf%id_ti .eq. 0 .or. zprf%id_te .eq. 0)then
       print *,'in set_nbprofs,zprf%id_ti or te  = 0'
       call STOP('set_nbprofs ',0)
   endif


  if(nlntmj) then
     call xpload_12(tmjsm,one,ispline,ig00,izbdy,'TMJSM',zprf%id_tifusn)
  else
     zprf%id_tifusn=zprf%id_ti
  endif
   if(zprf%id_tifusn .eq. 0)then
       print *,'in set_nbprofs,zprf%id_tifusion = 0'
       call STOP("sub set_nbprofs,Tmajority not defined",0)
   endif
  print *,'processed te,ti,tmjsm'

!--------------------------------------------------------------------
!  electron density
!--------------------------------------------------------------------
  onemill = 1.e6
  call xpload_12(ene,onemill,ispline,ig00,izbdy,'RHOEL',zprf%id_ne)
   if(zprf%id_ne .eq. 0)then
       print *,'in set_nbprofs,zprf%id_ne = 0'
       call STOP("sub set_nbprofs,ne  not defined",0)
   endif
  print *,'processed electron density'

!---------------------------------------------------------------------
!thermal primary ion  densities
!----------------------------------------------------------------------
  call set_thermal_den(zprf,ispline,ig00,izbdy)
  print *,'processed thermal ion densities'



!-----------------------------------------------------------------------
!  impurity density 
!  spatially varying Z & A allowed
!  ONETWO note: z varies due to te dep. but A is constant:
!------------------------------------------------------------------------
!   call set_impurity_den2(zprf,ispline,iknot,ig00,izcen) 
   call set_impurity_den(zprf,ispline,iknot,ig00,izbdy) 
   print *,'processed impurity densities'


!--------------------------------------------------------------------------
! neutral densities
! in Onetwo the # neutral species  is given by nneu and need not
! tbe the same as the number  of thermal ion species
! note that if there are 2 ions and only 1 neutral then
! ineut(1) = 1,ineut(2) =0 if neutral is species 1 type ions
! and ineut(1) =0,ineut(2) = 1 if neutral is species 2 type ions
!--------------------------------------------------------------------------
  call set_neutral_density(zprf,ispline,iknot,ig00,izbdy)

  print *,'processed neutral densities'



!---------------------------------------------------------------------------
!  minority densities
!  (ordering should match that defined in set_nbinputs)
!---------------------------------------------------------------------------
  zprf%id_nrf=0
  print *,'minority density ,nmini =',nmini
  do is=1,nmini
     zwork=rhmin(1:mj,2)*fracmini(is)
     iz = xzmini(is)+0.5
     ia = amini(is)+0.5
     call namels(iz,ia,zid)
     call STOP("sub set_nbprofs minority  not implemented",0)
     call xpload_12(zwork,onemill,ispline,ig00,izcen,'RHMIN_'//zid,zprf%id_nrf(is))
  enddo
  print *,'processed minority  densities'




!---------------------------------------------------------------------------
!  FPP fast ion densities
!  --this applies of the FPP code is used for the slowing down calculation
!---------------------------------------------------------------------------
  zprf%id_nfast = 0
  print *,'FPP, nsbeam =',nsbeam
  do is=1,nsbeam
     if(nlfbmfpp(is)) then
        iz = beam_data%xzbeams(is)+0.5   !note this is not set at presen, we have FP code
        ia = beam_data%abeams(is)+0.5
        if(beam_data%nlfprod(is) .gt. 0) then
           itype=4
        else
           itype=2
        endif

        isb = index_nbfi(iz,ia,itype)  ! NUBEAM index

        call namels(iz,ia,zid)
        call STOP("sub set_nbprofs FPP  not implemented",0)
        call xpload_12(rhbs(1,2,is),onemill,ispline,ig00,izcen,'RHBS_'//zid, &
             zprf%id_nfast(isb))
     endif
  enddo



!--------------------------------------------------------------------------
!  angular velocity (toroidal rotation)
!---------------------------------------------------------------------------
  call set_omegag
  call xpload_12(omegag,one,ispline,ig00,izbdy,'OMEGAG',zprf%id_omega)
  print *,'set_nbprof,id_omega =',zprf%id_omega



!---------------------------------------------------------------------------
!  radial electrostatic potential
!---------------------------------------------------------------------------
  call xpload_12(phiprg,one,ihermite,iknot,izbdy,'PHIPRG',zprf%id_epot)
  print *,'set_nbprof,id_epot =',zprf%id_epot


!---------------------------------------------------------------------------
!  area integral toroidal current
!---------------------------------------------------------------------------
  call set_curt
  call xpload_12(curt,one,ispline,ig00,izbdy,'CURT',zprf%id_curt)
  print *,'set_nbprof,id_curt =',zprf%id_curt



!---------------------------------------------------------------------------
!  toroidal loop voltage profile
!---------------------------------------------------------------------------
  call set_vpoh
  call xpload_12(vpoh,one,ispline,ig00,izbdy,'VPOH',zprf%id_vpoh)
  print *,'set_nbprof,id_vpoh =',zprf%id_vpoh




!---------------------------------------------------------------------------
!  anomolous diffusivity profile
!---------------------------------------------------------------------------
  call eq_gfnum("DIFB",zprf%id_difb)
  print *,'fast ion diffusivity profile',zprf%id_difb


!----------------------------------------------------------------------------
!  transfer ids
!----------------------------------------------------------------------------
!  call peek_nbi_com('a')
  call nbi_update_profiles(zprf)
     print *,'updated profiles'


!---------------------------------------------------------------------------
!  and now, interpolate profiles to the MC grid inside 
!---------------------------------------------------------------------------
!   call peek_nbi_com('b')
 call nbi_interp_profiles(ierr)  !in nbi_getprofiles.f90
!   call peek_nbi_com('c')
  if(ierr.ne.0) then
     call errmsg_exit(' ?set_nbprofs:  nbi_interp_profiles error.')
  endif


     print *,'leaving set_nbprofs'
end subroutine set_nbprofs


 subroutine update_nbprofs
!----------------------------------------------------------------------

  use nbi_types



  type (nbitype_profiles) :: zprf


!  call nbi_init_profiles(zprf)  !zero the ids

  call set_nbprofs              !re define the profiles ??

  return
  end subroutine update_nbprofs
