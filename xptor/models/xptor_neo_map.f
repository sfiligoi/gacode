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
      real :: mu2, mu3, xnue, dr_loc, drhodr_loc
      real :: taue,rminm,rhom,rmajm,qm,lnlamda
      real :: b_unit
      real :: m0,n0,a0,T0,v0,w0
      real :: nbm,gradnbm,tbm,gradtbm
   
      mu2 = sqrt(amassgas_exp)*42.851
      mu3 = sqrt(amassimp_exp)*42.851

! Initialize NEO
      call neo_init("path")

! Simulation mode (dke solve + analytic)
      neo_sim_model_in = 1

!  units: 
      m0 = amassgas_exp       !proton mass
      n0 = nem                !10^19/m^3
      a0 = rmin_exp(mxgrid)   !m
      T0 = tim                !kev
      v0 = 9.79D3*DSQRT(T0*1.D3)/DSQRT(m0) !m/s
      w0 = v0/a0              !1/s
! electron collision frequency
       lnlamda=15.94D0-0.5*LOG(nem)+LOG(tem)
       taue = 1.088D-3*(tem**1.5)/(nem*lnlamda)
!  xnue = 3/4 (Pi**0.5)/taue
       xnue = 1.329/taue
! Geometry
      rhom = arho_exp*(rho(jm+1)+rho(jm))/2.0
      rminm=(rmin_exp(jm+1)+rmin_exp(jm))/2.0
      rmajm=(rmaj_exp(jm+1)+rmaj_exp(jm))/2.0
      qm = (q_exp(jm+1)+q_exp(jm))/2.0
      dr_loc = rmin_exp(jm+1)-rmin_exp(jm)
      drhodr_loc = arho_exp*(rho(jm+1)-rho(jm))/dr_loc
!      neo_write_out_mode_in    = 0
      neo_equilibrium_model_in = 2
      neo_rmin_over_a_in = rminm/a0
      neo_rmaj_over_a_in = rmajm/a0
      neo_q_in           = ABS(qm)  ! absolute value of q
      neo_shear_in       = (rminm/qm)*(q_exp(jm+1)-q_exp(jm))/dr_loc
      neo_shift_in       = (rmaj_exp(jm+1)-rmaj_exp(jm))/dr_loc
      neo_kappa_in       =  0.5*(elong_exp(jm+1)+elong_exp(jm))
      neo_s_kappa_in     =  (rminm/neo_kappa_in)*
     >                (elong_exp(jm+1)-elong_exp(jm))/dr_loc
      neo_delta_in       =  0.5*(delta_exp(jm+1)+delta_exp(jm))
      neo_s_delta_in     = rminm*(delta_exp(jm+1)-delta_exp(jm))/dr_loc

      neo_ipccw_in = sign_It_exp  !current direction counter clockwise from above
      neo_btccw_in = sign_Bt_exp  !magnetic field direction "

!      neo_n_species_in = 4
      neo_n_species_in = 3

      b_unit = ABS(bt_exp)*(rhom/rminm)*drhodr(jm)
      neo_rho_star_in  = (1.02D2*DSQRT(m0*T0*1.D3)/
     >  (b_unit*1.D4))/(a0*100.D0)

  ! Electrons
      neo_z_1_in      = -1
      neo_mass_1_in   = 1.0/mu2**2
      neo_dens_1_in   = nem/n0
      neo_temp_1_in   = tem/tim
      neo_dlnndr_1_in = -drhodr_loc*a0*gradnem/nem
      neo_dlntdr_1_in = -drhodr_loc*a0*gradtem/tem
      neo_nu_1_in     = xnue/w0

  ! Main ions
      neo_z_2_in      = 1.0
      neo_mass_2_in   = 1.0
      neo_dens_2_in   = nim/n0
      neo_temp_2_in   = 1.0
      neo_dlnndr_2_in = -drhodr_loc*a0*gradnim/nim
      neo_dlntdr_2_in = -drhodr_loc*a0*gradtim/tim
      neo_nu_2_in     = (xnue/w0)*(nim/nem)/(mu2*(tim/tem)**1.5)

  ! Impurity
      neo_z_3_in      = zimp_exp
      neo_mass_3_in   = amassimp_exp/amassgas_exp
      neo_dens_3_in   = nzm/n0
      neo_temp_3_in   = 1.0
      neo_dlnndr_3_in = -drhodr_loc*a0*gradnzm/nzm
      neo_dlntdr_3_in = -drhodr_loc*a0*gradtim/tim
      neo_nu_3_in = (xnue/w0)*(nzm/nem)*(zimp_exp**4)
     >  /(mu3*(tim/tem)**1.5)

  ! fast ions
      nbm = 0.5*(nfast_exp(jm+1)+nfast_exp(jm))
      tbm = 0.5*(tfast_exp(jm+1)+tfast_exp(jm))
      gradnbm = (nfast_exp(jm+1)-nfast_exp(jm))/dr_loc
      gradtbm = (tfast_exp(jm+1)-tfast_exp(jm))/dr_loc
      neo_z_4_in      = 1.0
      neo_mass_4_in   = 1.0
      neo_dens_4_in   = nfast_exp(jm)/n0
      neo_temp_4_in   = tbm/t0
      neo_dlnndr_4_in = -a0*gradnbm/nbm
      neo_dlntdr_4_in = -a0*gradtbm/tbm
      neo_nu_4_in = (xnue/w0)*(nbm/nem)/(mu2*(tbm/tem)**1.5)

  ! Rotation is always active
      neo_rotation_model_in = 2
!      neo_rotation_model_in = 1
      neo_omega_rot_in = vexbm*cv/(w0*rmajor_exp)
      neo_omega_rot_in = 0.0  ! eliminate mach number corrections
      neo_omega_rot_deriv_in = drhodr_loc*a0*gradvexbm*cv
     > /(w0*rmajor_exp)
!      neo_omega_rot_deriv_in = 0.0

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
