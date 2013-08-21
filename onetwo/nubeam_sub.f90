subroutine r8_nbfusn_emax(iZ,iA,emax)

  implicit NONE

  integer, intent(in) :: iZ  ! Z of fusion product
  integer, intent(in) :: iA  ! A of fusion product

  real*8, intent(out) :: emax  ! a number safely above the birth energy

  !-------------------------------------

  emax = 0.0d0

  if(iz.eq.1) then
     if(ia.eq.1) then
        emax = 5.00d6
     else if(ia.eq.3) then
        emax = 1.25d6
     endif
  else if(iz.eq.2) then
     if(ia.eq.3) then
        emax = 1.00d6
     else if(ia.eq.4) then
        emax = 5.00d6
     endif
  endif

end subroutine r8_nbfusn_emax

subroutine r8_nbi_efbm(ilfprod,inum,zemax,ilin,zelin,efbmb,efbm)
  !
  !  compute the energy grid for a MC fast ion specie dist. fcn.
  !  (code taken from old trcore/datckb.for -- dmc 10/04/01)
 
  implicit NONE
 
  logical,intent(in) :: ilfprod  ! .TRUE. if this is a fusion product specie
  integer,intent(in) ::  inum    ! total no. of zones
  real*8,intent(in) ::  zemax      ! max energy
 
  integer,intent(out) ::  ilin   ! no. of linearly spaced zones; ilin.ge.1
  real*8,intent(out) ::  zelin     ! max. linearly spaced energy (min is zero)
 
  real*8,intent(out) ::  efbmb(inum+1)  ! energy grid bdys
  real*8,intent(out) ::  efbm(inum)     ! energy grid zone ctrs

  !  beam species -- linear grid
  !  fusion products -- linear at low energies, otherwise logarithmic
 
  !$r8real_input: zelin,zemax
  !$r8real_output: efbmb,efbm
 
  ! note if ilin.eq.inum, zelin.eq.zemax is expected; no checks here!
 
 
  ! local...
 
  integer il,ilp1
  real*8 zfac
 
  !---------------------------------
  if(ilfprod) then
     ilin=inum/10
     zelin=zemax/100
  else
     ilin=inum
     zelin=zemax
  endif
 
  EFBMB(1)=0.0
 
!  LINEAR PART
  DO IL=1,ILIN
     ILP1=IL+1
     EFBMB(ILP1)=IL*ZELIN/ILIN
  enddo
 
!  LOG PART
  IF(INUM.GT.ILIN) THEN
     ZFAC=EXP(LOG(ZEMAX/ZELIN)/(INUM-ILIN))
     DO IL=ILIN+1,INUM
        ILP1=IL+1
        EFBMB(ILP1)=ZFAC*EFBMB(IL)
     enddo
  ENDIF
 
!  FILL IN ZONE CTR VALUES
  DO IL=1,INUM
     EFBM(IL)=(EFBMB(IL)+EFBMB(IL+1))/2
  enddo
 
end subroutine r8_nbi_efbm
