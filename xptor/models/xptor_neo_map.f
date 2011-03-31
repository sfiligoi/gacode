      subroutine xptor_neo_map

      use neo_interface
      IMPLICIT NONE
      include 'mpif.h'
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/share.m'
      include '../inc/sharegk.m'
      include '../inc/ptor.m'
      include '../inc/glf.m'
c
      real :: mu1, mu3, xnue, dr_loc, drhodr_loc
      real :: taue,rminm,rmajm,qm,lnlamda
      real :: n0,a0,T0,v0,w0

      mu1 = sqrt(amassgas_exp)*42.851
      mu3 = sqrt(amassimp_exp)*42.851

! Initialize NEO
      call neo_init

! Simulation mode (dke solve + analytic)
      neo_sim_model_in = 1

!  units: 
      n0 = nem                !10^19/m^3
      a0 = rmin_exp(mxgrid)   !m
      T0 = tim                !kev
      v0 = 9.79D3*DSQRT(T0*1.D3)/DSQRT(amassgas_exp) !m/s
      w0 = v0/a0              !1/s
! electron collision frequency
       lnlamda=15.94D0-0.5*LOG(nem)+LOG(tem)
       taue = 1.088D-3*(tem**1.5)/(nem*lnlamda)
!  xnue = 3/4 (Pi**0.5)/taue
       xnue = 1.329/taue

! neo resolution
      neo_energy_max_in = 16.0

! Geometry
      rminm=(rmin_exp(jm+1)+rmin_exp(jm))/2.D0
      rmajm=(rmaj_exp(jm+1)+rmaj_exp(jm))/2.0
      qm = (q_exp(jm+1)+q_exp(jm))/2.0
      dr_loc = rmin_exp(jm+1)-rmin_exp(jm)
      drhodr_loc = arho_exp*(rho(jm+1)-rho(jm))/dr_loc
      neo_write_out_mode_in    = 0
      neo_equilibrium_model_in = 2
      neo_rmin_over_a_in = rminm/a0
      neo_rmaj_over_a_in = rmajm/a0
      neo_q_in           = qm*sign_Bt_exp
      neo_shear_in       = (rminm/qm)*(q_exp(jm+1)-q_exp(jm))/dr_loc
      neo_shift_in       = (rmaj_exp(jm+1)-rmaj_exp(jm))/dr_loc
      neo_kappa_in       =  0.5*(elong_exp(jm+1)+elong_exp(jm))
      neo_s_kappa_in     =  (rminm/neo_kappa_in)*
     >                (elong_exp(jm+1)-elong_exp(jm))/dr_loc
      neo_delta_in       =  0.5*(delta_exp(jm+1)+delta_exp(jm))
      neo_s_delta_in     = rminm*(delta_exp(jm+1)-delta_exp(jm))/dr_loc

      neo_ipccw_in = 1.0   !xptor chooses signs relative to It
      neo_btccw_in = sign_Bt_exp

      neo_n_species_in = 3
      neo_rho_star_in  = 0.001

  ! Electrons
      neo_z_2_in      = -1
      neo_mass_2_in   = 1.0/mu1**2
      neo_dens_2_in   = nem/n0
      neo_temp_2_in   = tem/tim
      neo_dlnndr_2_in = -a0*gradnem/nem
      neo_dlntdr_2_in = -a0*gradtem/tem
      neo_nu_2_in     = xnue/w0

  ! Main ions
      neo_z_1_in      = 1.0
      neo_mass_1_in   = 1.0
      neo_dens_1_in   = nim/n0
      neo_temp_1_in   = 1.0
      neo_dlnndr_1_in = -a0*gradnim/nim
      neo_dlntdr_1_in = -a0*gradtim/tim
      neo_nu_1_in     = (xnue/w0)/(mu1*(tim/tem)**1.5)

  ! Impurity
      neo_z_3_in      = zimp_exp
      neo_mass_3_in   = amassimp_exp/amassgas_exp
      neo_dens_3_in   = nzm/n0
      neo_dlnndr_3_in = -a0*gradnzm/nzm
      neo_dlntdr_3_in = -a0*gradtim/tim
      neo_nu_3_in = (xnue/w0)/(mu3*(tim/tem)**1.5)

  ! Rotation is always active
      neo_rotation_model_in = 1
      neo_omega_rot_in = vexbm*cv/w0
      neo_omega_rot_deriv_in = gradvexbm*cv/w0

  ! Parameter only used for global runs.
      neo_rmin_over_a_2_in = neo_rmin_over_a_in

  ! General geometry Fourier coefficients
c not yet availible for xptor
c      if (loc_num_equil_flag == 1) then
c     ! Resetting if numerical equilibrium wanted.
c       neo_equilibrium_model_in = 3
c      endif
c
c      neo_geo_ny_in = n_fourier_geo
c      neo_geo_yin_in(:,:) = a_fourier_geo(:,:,i_r)

      end subroutine xptor_neo_map
