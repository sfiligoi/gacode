subroutine expromake_read_input

  use expromake_globals

  implicit none

  open(unit=1,file='input.expromake.gen',status='old')
  read(1,*) exm_z(1)
  read(1,*) exm_z(2)
  read(1,*) exm_z(3)
  read(1,*) exm_b_ref
  read(1,*) exm_arho
  read(1,*) exm_kappa
  read(1,*) exm_delta
  read(1,*) exm_te_axis
  read(1,*) exm_alte
  read(1,*) exm_ti_axis
  read(1,*) exm_alti
  read(1,*) exm_set_b_ref
  read(1,*) exm_set_arho
  read(1,*) exm_set_kappa
  read(1,*) exm_set_delta
  read(1,*) exm_set_te_axis
  read(1,*) exm_set_alte
  read(1,*) exm_set_ti_axis
  read(1,*) exm_set_alti
  close(1)

end subroutine expromake_read_input
