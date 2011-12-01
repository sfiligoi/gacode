program locpargen

  use locpargen_globals
  use EXPRO_interface

  implicit none
  real :: r0
  real :: a
  real, dimension(1) :: x
  real, dimension(1) :: y


  open(unit=1,file='input.locpargen',status='old')
  read(1,*) r0
  close(1)

  EXPRO_ctrl_density_method = 1
  EXPRO_ctrl_z = 1.0
  EXPRO_ctrl_numeq_flag = 0 
  EXPRO_ctrl_signq = (-1)*(-1)
  EXPRO_ctrl_signb = -(-1)
  EXPRO_ctrl_rotation_method = 1
 
  call EXPRO_alloc('./',1) 
  call EXPRO_read

  ! Minor radius
  a = EXPRO_rmin(EXPRO_n_exp)
  x(1) = r0

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

  !---------------------------------------

  ! B_unit
  call cub_spline(EXPRO_rmin/a,EXPRO_bunit,EXPRO_n_exp,x,y,1)
  print 10,'B_unit : ',y(1)

  call EXPRO_alloc('./',0) 

10 format(a,sp1pe12.5)

end program locpargen
