program locpargen

  use EXPRO_interface

  implicit none

  real :: r0
  real :: a
  real, dimension(1) :: x
  real, dimension(1) :: y
  real, dimension(3) :: z


  open(unit=1,file='input.locpargen',status='old')
  read(1,*) r0
  read(1,*) z(1)
  read(1,*) z(2)
  read(1,*) z(3)
  close(1)

  EXPRO_ctrl_density_method = 1
  EXPRO_ctrl_z(1:3) = z(1:3)
  EXPRO_ctrl_numeq_flag = 0 
  EXPRO_ctrl_signq = (-1)*(-1)
  EXPRO_ctrl_signb = -(-1)
  EXPRO_ctrl_rotation_method = 1

  call EXPRO_alloc('./',1) 
  call EXPRO_read

  ! Minor radius
  a = EXPRO_rmin(EXPRO_n_exp)
  x(1) = r0

  print *
  print '(a)','Local input parameters'
  print *

  ! RADIUS
  print 10,'RADIUS=',r0

  ! ASPECT_RATIO
  call cub_spline(EXPRO_rmin/a,EXPRO_rmaj/a,EXPRO_n_exp,x,y,1)
  print 10,'ASPECT_RATIO=',y(1)

  ! SHIFT
  call cub_spline(EXPRO_rmin/a,EXPRO_drmaj,EXPRO_n_exp,x,y,1)
  print 10,'SHIFT=',y(1)

  ! ZMAG
  call cub_spline(EXPRO_rmin/a,EXPRO_zmag/a,EXPRO_n_exp,x,y,1)
  print 10,'ZMAG=',y(1)

  ! DZMAG
  call cub_spline(EXPRO_rmin/a,EXPRO_dzmag,EXPRO_n_exp,x,y,1)
  print 10,'DZMAG=',y(1)

  ! SHEAR
  call cub_spline(EXPRO_rmin/a,EXPRO_s,EXPRO_n_exp,x,y,1)
  print 10,'SHEAR=',y(1)

  ! SAFETY_FACTOR
  call cub_spline(EXPRO_rmin/a,EXPRO_q,EXPRO_n_exp,x,y,1)
  print 10,'SAFETY_FACTOR=',y(1)

  ! KAPPA
  call cub_spline(EXPRO_rmin/a,EXPRO_kappa,EXPRO_n_exp,x,y,1)
  print 10,'KAPPA=',y(1)

  ! S_KAPPA
  call cub_spline(EXPRO_rmin/a,EXPRO_skappa,EXPRO_n_exp,x,y,1)
  print 10,'S_KAPPA=',y(1)

  ! DELTA
  call cub_spline(EXPRO_rmin/a,EXPRO_delta,EXPRO_n_exp,x,y,1)
  print 10,'DELTA=',y(1)

  ! S_DELTA
  call cub_spline(EXPRO_rmin/a,EXPRO_sdelta,EXPRO_n_exp,x,y,1)
  print 10,'S_DELTA=',y(1)

  ! ZETA
  call cub_spline(EXPRO_rmin/a,EXPRO_zeta,EXPRO_n_exp,x,y,1)
  print 10,'ZETA=',y(1)

  ! S_ZETA
  call cub_spline(EXPRO_rmin/a,EXPRO_szeta,EXPRO_n_exp,x,y,1)
  print 10,'S_ZETA=',y(1)

  ! TI_OVER_TE
  call cub_spline(EXPRO_rmin/a,EXPRO_ti(1,:)/EXPRO_te,EXPRO_n_exp,x,y,1)
  print 10,'TI_OVER_TE=',y(1)

  ! TI_OVER_TE_2
  call cub_spline(EXPRO_rmin/a,EXPRO_ti(2,:)/EXPRO_te,EXPRO_n_exp,x,y,1)
  print 10,'TI_OVER_TE_2=',y(1)

  ! NI_OVER_NE
  call cub_spline(EXPRO_rmin/a,EXPRO_ni(1,:)/EXPRO_ne,EXPRO_n_exp,x,y,1)
  print 10,'NI_OVER_NE=',y(1)

  ! NI_OVER_NE_2
  call cub_spline(EXPRO_rmin/a,EXPRO_ni(2,:)/EXPRO_ne,EXPRO_n_exp,x,y,1)
  print 10,'NI_OVER_NE_2=',y(1)

  ! DLNNDR
  call cub_spline(EXPRO_rmin/a,a*EXPRO_dlnnidr(1,:),EXPRO_n_exp,x,y,1)
  print 10,'DLNNDR=',y(1)

  ! DLNNDR_2
  call cub_spline(EXPRO_rmin/a,a*EXPRO_dlnnidr(2,:),EXPRO_n_exp,x,y,1)
  print 10,'DLNNDR_2=',y(1)

  ! DLNNDR_ELECTRON
  call cub_spline(EXPRO_rmin/a,a*EXPRO_dlnnedr,EXPRO_n_exp,x,y,1)
  print 10,'DLNNDR_ELECTRON=',y(1)

  ! DLNTDR
  call cub_spline(EXPRO_rmin/a,a*EXPRO_dlntidr(1,:),EXPRO_n_exp,x,y,1)
  print 10,'DLNTDR=',y(1)

  ! DLNTDR_2
  call cub_spline(EXPRO_rmin/a,a*EXPRO_dlntidr(2,:),EXPRO_n_exp,x,y,1)
  print 10,'DLNTDR_2=',y(1)

  ! DLNTDR_ELECTRON
  call cub_spline(EXPRO_rmin/a,a*EXPRO_dlntedr,EXPRO_n_exp,x,y,1)
  print 10,'DLNTDR_ELECTRON=',y(1)

  ! GAMMA_E
  call cub_spline(EXPRO_rmin/a,EXPRO_gamma_e*a/EXPRO_cs,EXPRO_n_exp,x,y,1)
  print 10,'GAMMA_E=',y(1)

  ! RHO_STAR
  call cub_spline(EXPRO_rmin/a,EXPRO_rhos/a,EXPRO_n_exp,x,y,1)
  print 10,'RHO_STAR=',y(1)

  print *
  print '(a)','Additional dimensional quantities'
  print *

  !---------------------------------------
  ! Some added physical quantities
  call cub_spline(EXPRO_rmin/a,EXPRO_bunit,EXPRO_n_exp,x,y,1)
  print 10,'B_unit [T]  : ',y(1)
  call cub_spline(EXPRO_rmin/a,EXPRO_cs,EXPRO_n_exp,x,y,1)
  print 10,'c_s   [m/s] : ',y(1)
  call cub_spline(EXPRO_rmin/a,EXPRO_rhos,EXPRO_n_exp,x,y,1)
  print 10,'rhos   [m]  : ',y(1)
  !---------------------------------------

  call EXPRO_alloc('./',0) 

10 format(a,sp1pe12.5)

end program locpargen
