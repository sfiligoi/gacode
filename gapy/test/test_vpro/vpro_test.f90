program test

  use vpro

  implicit none

  call vpro_read_legacy
  call vpro_write
  call vpro_init(0)

  call vpro_read()

  print '(a,i0)'            ,'nexp    ',EXPRO_n_exp
  print '(a,1pe12.5)'       ,'bt_exp  ',EXPRO_b_ref
  print '(a,1pe12.5)'       ,'arho_exp',EXPRO_arho
  print '(a,60(1pe12.5,1x))','ni      ',EXPRO_ni(1,:)

  call vpro_compute_derived()

  print '(a,60(1pe12.5,1x))','bunit   ',EXPRO_bunit

end program test
