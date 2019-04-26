program test

  use expro
  use expro_locsim_interface

  implicit none

  real :: ipccw,btccw,a_meters

  !call expro_read_legacy
  !call expro_write

  call expro_read('input.gacode')

  print '(2a)'               ,'name1    ',expro_name(1)
  print '(2a)'               ,'name2    ',expro_name(2)
  print '(a,i0)'            ,'nexp     ',expro_n_exp
  print '(a,1pe12.5)'       ,'torfluxa ',expro_torfluxa
  print '(a,60(1pe12.5,1x))','ni       ',expro_ni(1,:)
  print '(a,60(1pe12.5,1x))','bunit    ',expro_bunit

  stop
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
