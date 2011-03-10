! variables for read_eqdsk
      character*52   etitle
      integer        nx_xy, ny_xy
      integer        n_bdry, n_lim
      real           rdim, zdim, a0, r0, rmin, zmid,
     &               rmag, zmag, psimag, psilim, 
     &               bt0, q0_efit, el0, current, dum1
c
c... Limiter surface-lim
      integer        mxn_lim
      parameter     (mxn_lim=200)
      real           x_lim(mxn_lim), y_lim(mxn_lim)
c
c... Plasma boundary-bdry
      integer        mxn_bdry
      parameter     (mxn_bdry=1500)
      real           x_bdry(mxn_bdry), y_bdry(mxn_bdry)
c
c... Radial profiles
c    nr_r=number of radial pts
      integer        mxnr_r, nr_r
      parameter     (mxnr_r=300)
      real           fm_r(3,mxnr_r)
      real           psi_r(mxnr_r), psip_r(mxnr_r),
     &               rhop_r(mxnr_r), frb_r(mxnr_r),
     &               p_r(mxnr_r), ffp_r(mxnr_r), pp_r(mxnr_r),
     &               q_r(mxnr_r), rhot_r(mxnr_r), rin_r(mxnr_r), 
     &               rout_r(mxnr_r), elong_r(mxnr_r), vol_r(mxnr_r),
     &               phit_r(mxnr_r), grho1_r(mxnr_r), grho2_r(mxnr_r),
     &               bpout_r(mxnr_r), btout_r(mxnr_r), 
     &               rhor_r(mxnr_r), f_r(mxnr_r), vp_r(mxnr_r),
     &               grth_r(mxnr_r), gph_r(mxnr_r), gth_r(mxnr_r),
     &               rm1_r(mxnr_r), rm2_r(mxnr_r),
     &               b2_r(mxnr_r), bm2_r(mxnr_r), fhat_r(mxnr_r),
     &               ftrap_r(mxnr_r), gr2bm2_r(mxnr_r), 
     &               xj_r(mxnr_r), xj_nb_ex_r(mxnr_r)
c
c... EFIT Profiles
      integer        mxn_e
      parameter      (mxn_e=130)
      real           psi_e(mxn_e), rhop_e(mxn_e),
     &               f_e(mxn_e), ffp_e(mxn_e),      
     &               p_e(mxn_e), pp_e(mxn_e),
     &               q_e(mxn_e), bp_e(mxn_e)
c
      integer        mx_ni, izim1, izim2
      parameter      (mx_ni=5)
      integer        k_edotb, k_potato, izi0(mx_ni)
      real           k_sqz, c_den,
     &               amuim1, amuim2, amui0(mx_ni)
      parameter      (c_den=1.e10)
c
      common /efitdat/ x_bdry, y_bdry, fm_r,
     &               psi_r, psip_r, rhop_r,
     &               frb_r, p_r, ffp_r, pp_r, q_r, rhot_r,
     &               rin_r, rout_r, elong_r, vol_r, phit_r,
     &               grho1_r, grho2_r, 
     &               bpout_r, btout_r, rhor_r, f_r, vp_r,
     &               grth_r,gph_r,gth_r, rm1_r, rm2_r, b2_r, bm2_r,
     &               fhat_r, ftrap_r, gr2bm2_r, xj_r, xj_nb_ex_r,
     &               x_lim, y_lim,
     &               psi_e, rhop_e, f_e, ffp_e,    
     &               p_e, pp_e, q_e, bp_e,
     &               rdim, zdim, a0, r0, rmin, zmid,
     &               rmag, zmag, psimag, psilim, 
     &               bt0, q0_efit, el0, current, dum1,
     &               nx_xy, ny_xy, n_bdry, n_lim,
     &               nr_r,etitle
      common /ncldat/ amui0, k_edotb, k_potato, k_sqz, 
     &               amuim1, amuim2, izi0,izim1, izim2
