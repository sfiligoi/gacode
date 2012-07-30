subroutine allocate_plasmastate_vars

  use prgen_globals

  implicit none

  allocate(plst_all_name(plst_dp1_nspec_all))
  allocate(plst_alla_name(plst_dp1_nspec_alla))
  allocate(plst_ts(nx,plst_dp1_nspec_th+1))
  allocate(plst_ns(nx,plst_dp1_nspec_th+1))
  allocate(plst_ptowb(nx))
  allocate(plst_nb(nx))
  allocate(plst_eparb(nx))
  allocate(plst_eperpb(nx))
  allocate(plst_vol(nx))
  allocate(plst_rho(nx))
  allocate(plst_grho1(nx))
  allocate(plst_phit(nx))
  allocate(plst_psipol(nx))
  allocate(plst_elong(nx))
  allocate(plst_triang(nx))
  allocate(plst_iota(nx))
  allocate(plst_r_midp_in(nx))
  allocate(plst_r_midp_out(nx))
  allocate(plst_z_midp(nx))
  allocate(plst_zeff(nx))
  allocate(plst_epot(nx))
  allocate(plst_omegat(nx))
  allocate(plst_pbe(nx))
  allocate(plst_pbi(nx))
  allocate(plst_pe_trans(nx))
  allocate(plst_pi_trans(nx))
  allocate(plst_pei_trans(nx))
  allocate(plst_tq_trans(nx))
  allocate(plst_sn_trans(nx))

end subroutine allocate_plasmastate_vars
