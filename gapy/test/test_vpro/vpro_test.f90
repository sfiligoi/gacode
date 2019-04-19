program test

  use vpro
  use expro_locsim_interface

  implicit none

  real :: ipccw,btccw,a_meters

  call vpro_read_legacy
  call vpro_write

  call vpro_read('./')

  print '(a,i0)'            ,'nexp    ',expro_n_exp
  print '(a,1pe12.5)'       ,'bt_exp  ',expro_b_ref
  print '(a,1pe12.5)'       ,'arho_exp',expro_arho
  print '(a,60(1pe12.5,1x))','ni      ',expro_ni(1,:)
  print '(a,60(1pe12.5,1x))','bunit   ',expro_bunit

  call expro_locsim_profiles('./',&
       -1,&
       0,&
       1,&
       0,&
       expro_n_ion+1,&
       0.5,&
       btccw,&
       ipccw,&
       a_meters)

  print '(a,1pe12.5)'       ,'kappa_loc',kappa_loc

end program test
