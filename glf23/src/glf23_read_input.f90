!--------------------------------------------------------------
! glf23_read_input.f90
!
! PURPOSE:
!  Complete read of glf23 input parameters as specified by 
!  glf23_interface 
!--------------------------------------------------------------

subroutine glf23_read_input

  use glf23_interface

  implicit none


  open(unit=1,file=trim(glf23_path_in)//'input.glf23.gen',status='old')

  read(1,*) glf23_use_transport_model_in

  ! Data passed to: put_model_parameters
  read(1,*) glf23_alpha_p_in
  read(1,*) glf23_alpha_quench_in
  read(1,*) glf23_version_in

  ! Data passed to: put_species
  read(1,*) glf23_ns_in
  read(1,*) glf23_mass_in(1)
  read(1,*) glf23_mass_in(2)
  read(1,*) glf23_mass_in(3)
  read(1,*) glf23_zs_in(1)
  read(1,*) glf23_zs_in(2)
  read(1,*) glf23_zs_in(3)

  ! Data passed to: put_kys
  read(1,*) glf23_ky_in

  ! Data passed to: put_gradients
  read(1,*) glf23_rlns_in(1)
  read(1,*) glf23_rlns_in(2)
  read(1,*) glf23_rlns_in(3)
  read(1,*) glf23_rlts_in(1)
  read(1,*) glf23_rlts_in(2)
  read(1,*) glf23_rlts_in(3)
  read(1,*) glf23_vpar_shear_in(1)
  read(1,*) glf23_vpar_shear_in(2)
  read(1,*) glf23_vpar_shear_in(3)
  read(1,*) glf23_vexb_shear_in

  ! Data passed to: put_averages
  read(1,*) glf23_taus_in(1)
  read(1,*) glf23_taus_in(2)
  read(1,*) glf23_taus_in(3)
  read(1,*) glf23_as_in(1)
  read(1,*) glf23_as_in(2)
  read(1,*) glf23_as_in(3)
  read(1,*) glf23_betae_in
  read(1,*) glf23_xnue_in

  ! Data passed to: put_s_alpha_geometry
  read(1,*) glf23_rmin_sa_in
  read(1,*) glf23_rmaj_sa_in
  read(1,*) glf23_q_sa_in
  read(1,*) glf23_shat_sa_in
  read(1,*) glf23_alpha_sa_in
  read(1,*) glf23_xwell_sa_in
  read(1,*) glf23_theta0_sa_in

  close(1)

end subroutine glf23_read_input
