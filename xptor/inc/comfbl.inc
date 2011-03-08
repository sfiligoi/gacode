!COMFBL-radial variables for FORCEBAL
!Current density and electrical resistivity
!  eta_par_r-parallel electrical resistivity (Ohm*m)
!  xj_bs_r-bootstrap current density=<Jbc.B>/Bt0 (A/m**2)
!  xj_ex_r-<F.B> driven current=<Jex.B>/bt0 (A/m**2)
!Particle transport - single charge state
!  d_eff_r-eff part diffusivity=-gam_r/(dn/(dr_t)) (m**2/s)
!  d_eff_g_r-eff part diffusivity=d_eff_r/<grad(r_t)^2> (m**2/s)
!  d_n_r-diagonal diffusivity (m**2/s)
!  d_n_g_r-d_n_r/<grad(r_t)^2> (m**2/s)
!  g_n_r-coefficient for density driven flux (-)
!  g_te_r-coefficient for electron temperature driven flux (-)
!  g_ti_r-coefficient for ion temperature driven flux (-)
!  gam_r-net radial flux=<Gam.grad(r_t)> (/m**2/s)
!  v_eb_r-<E.B> convective particle velocity (m/s)
!  v_eb_g_r-v_eb_r/(<grad(r_t)^2>)**0.5 (m/s)
!  v_ex_r-<F.B> convective particle velocity (m/s)
!  v_ex_g_r-v_ex_r/(<grad(r_t)^2>)**0.5 (m/s)
!  v_nt_r-nondiffusive (non-WP) velocity (m/s)
!  v_nt_g_r-v_n_r/(<grad(r_t)^2>)**0.5 (m/s)
!Particle transport - multiple charge state (MCS)
!  d_n_g_im2_r-d_n_r/<grad(r_t)^2> for MCS imp (m**2/s)
!  gam_im2_r-net nc radial flux=<Gam.grad(r_t)> for MCS imp (/m**2/s)
!  v_eb_g_im2_r-<E.B> convective particle velocity for MCS imp (m/s)
!  v_ex_g_im2_r-<F.B> convective particle velocity for MCS imp (m/s)
!  v_nt_g_im2_r-v_n_r/(<grad(r_t)^2>)**0.5 for MCS imp (m/s)
!Heat transport
!  chi_eff_r-eff cond=-q_con_r/(ne*dkTe/(drt)) (m**2/s)
!  chi_eff_g_r-chi_eff_r/<grad(r_t)^2> (m**2/s)
!  q_con_r-net nc cond flux=<q_con.grad(r_t)> (w/m**2)
!Rotation velocities
!  u_par_r-<u.B> (T*m/s)
!  u_p_r-<u.grad(theta)>/<B.grad(theta)> (m/s/T)
!  v_p_o_r-poloidal velocity on outside of torus (m/s)
!  v_t_o_r-toroidal velocity on outside of torus (m/s)
!  xm_p_o_r-poloidal Mach number on outside of torus
!  xm_p_ex_im_o_r-expt poloidal Mach number of im on outside of torus
!  xm_t_o_r-toroidal Mach number on outside of torus
!  xm_t_ex_im_o_r-expt toroidal Mach number of im on outside of torus
!  e_rad_r-E_r and components=-dPhi/(dr_t) (V/m)
!  e_rad_ex_r-E_r using expt p', vp and vt (V/m)
!  e_rad_o_r-E_r and components on outside of torus (V/m)
!  e_rad_ex_o_r-E_r on outside using expt p', vp and vt (V/m)
!  edotb_r-<E.B> (V*T/m)
!  omexb_o_r-ExB shear damping rate on outside of torus (rad/s)
!  omexb_ex_o_r-ExB shear damping rate on outside of torus using expt Er  (rad/s)
!  sqz_r-orbit squeezing factor
!  omega_r-toroidal rotation frequency (rad/s)
      REAL           eta_par_r,               xj_bs_r,
     &               xj_ex_r,
     &               d_eff_r,                 d_eff_g_r,
     &               d_n_r,                   d_n_g_r,
     &               g_n_r,
     &               g_te_r,                  g_ti_r,
     &               gam_r,
     &               v_eb_r,                  v_eb_g_r,
     &               v_ex_r,                  v_ex_g_r,
     &               v_nt_r,                  v_nt_g_r,
     &               d_n_g_im2_r,             gam_im2_r,
     &               v_eb_g_im2_r,            v_ex_g_im2_r,
     &               v_nt_g_im2_r,
     &               chi_eff_r,               chi_eff_g_r,
     &               q_con_r,
     &               e_rad_r,                 e_rad_ex_r,
     &               e_rad_o_r,               e_rad_ex_o_r,
     &               edotb_r,                 e_par_r,
     &               omega_r,                 sqz_r,
     &               u_par_r,                 u_p_r,
     &               v_p_o_r,                 v_t_o_r,
     &               xm_p_o_r,                xm_p_im_ex_o_r,
     &               xm_t_o_r,                xm_t_im_ex_o_r,
     &               omexb_o_r,               omexb_ex_o_r
      COMMON/COMFBL/ eta_par_r(mxnr_r),       xj_bs_r(mxnr_r),
     &               xj_ex_r(mxnr_r),
     &               d_eff_r(mxnr_r,mx_ms),   d_eff_g_r(mxnr_r,mx_ms),
     &               d_n_r(mxnr_r,mx_ms),     d_n_g_r(mxnr_r,mx_ms),
     &               g_n_r(mxnr_r,mx_ms,mx_ms),
     &               g_te_r(mxnr_r,mx_ms),    g_ti_r(mxnr_r,mx_ms),
     &               gam_r(mxnr_r,mx_ms),
     &               v_eb_r(mxnr_r,mx_ms),    v_eb_g_r(mxnr_r,mx_ms),
     &               v_ex_r(mxnr_r,mx_ms),    v_ex_g_r(mxnr_r,mx_ms),
     &               v_nt_r(mxnr_r,mx_ms),    v_nt_g_r(mxnr_r,mx_ms),
     &               d_n_g_im2_r(mxnr_r,mx_mz),
     &               gam_im2_r(mxnr_r,mx_mz),
     &               v_eb_g_im2_r(mxnr_r,mx_mz),
     &               v_ex_g_im2_r(mxnr_r,mx_mz),
     &               v_nt_g_im2_r(mxnr_r,mx_mz),
     &               chi_eff_r(mxnr_r,2),     chi_eff_g_r(mxnr_r,2),
     &               q_con_r(mxnr_r,2),
     &               e_rad_r(mxnr_r,4),       e_rad_ex_r(mxnr_r),
     &               e_rad_o_r(mxnr_r,4),     e_rad_ex_o_r(mxnr_r),
     &               e_par_r(mxnr_r),         edotb_r(mxnr_r),
     &               omega_r(mxnr_r,mx_ms),   sqz_r(mxnr_r,mx_ms),
     &               u_par_r(mxnr_r,mx_ms),   u_p_r(mxnr_r,mx_ms),
     &               v_p_o_r(mxnr_r,mx_ms),   v_t_o_r(mxnr_r,mx_ms),
     &               xm_p_o_r(mxnr_r,mx_ms),  xm_p_im_ex_o_r(mxnr_r),
     &               xm_t_o_r(mxnr_r,mx_ms),  xm_t_im_ex_o_r(mxnr_r),
     &               omexb_o_r(mxnr_r),       omexb_ex_o_r(mxnr_r)
