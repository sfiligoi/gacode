!-----------------------------------------------------------------
!  input routines
!-----------------------------------------------------------------
      SUBROUTINE put_species(nsp,zsp,msp)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN):: nsp
      REAL,INTENT(IN) :: zsp(nsm),msp(nsm)
      INTEGER :: is
!
      ns_in = nsp
      if(ns_in.lt.2.or.ns_in.gt.nsm)then
        write(*,*)"number of species must be >=2 or <=",nsm
      else
! transfer values
        do is=1,nsp
          zs_in(is)=zsp(is)
          mass_in(is)=msp(is)
        enddo
        new_matrix = .TRUE.
      endif
!
      END SUBROUTINE put_species
!
!-----------------------------------------------------------------
!
      SUBROUTINE put_kys(kys)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      REAL,INTENT(IN) :: kys
!
! check for changes and update flow controls
! 
      if(kys.ne.ky_in)new_matrix = .TRUE.
!
! transfer values
!
      ky_in = kys
!
      END SUBROUTINE put_kys
!
!-----------------------------------------------------------------
!
      SUBROUTINE put_gaussian_width(width,width_min,nwidth,find_width)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: find_width
      INTEGER,INTENT(IN) :: nwidth
      REAL,INTENT(IN) :: width,width_min
!
! check for changes and update flow controls
! 
      if(width.ne.width_in)new_width = .TRUE.
!
! transfer values
!
      width_in = width
      width_min_in=width_min
      find_width_in = find_width
      nwidth_in=MIN(nwidth,nt0)
!
      END SUBROUTINE put_gaussian_width
!
!
      SUBROUTINE put_eikonal(new_eikonal)
!*********************************************
!
!*********************************************
      USE tglf_internal_interface
!
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: new_eikonal
!
!  set flow control switch
!
      new_eikonal_in = new_eikonal
! check consistency
      if(new_eikonal_in)then
        eikonal_unsaved=.TRUE.
      else
        if(eikonal_unsaved)then
          write(*,*)"warning put_eikonal:"
          write(*,*)"new_eikonal = .FALSE.attempted before call with new_eikonal=.TRUE."
          new_eikonal_in = .TRUE.
        endif 
      endif    
!
      END SUBROUTINE put_eikonal
!
!-----------------------------------------------------------------
!
      SUBROUTINE put_gradients(rln,rlt,vexb_shear,vpar_shear)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      REAL,INTENT(IN) :: rln(nsm),rlt(nsm),vexb_shear,vpar_shear
      INTEGER :: is
!
! transfer values
!
      do is=1,nsm
        rlns_in(is) = rln(is)
        rlts_in(is) = rlt(is)
      enddo
      vexb_shear_in = vexb_shear
      vpar_shear_in = vpar_shear
!    
      END SUBROUTINE put_gradients
!
!-----------------------------------------------------------------
!
      SUBROUTINE put_averages(tsp,asp,betae,xnuei,zeff,debye,vexb,vpar)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      REAL,INTENT(IN) :: tsp(nsm),asp(nsm),betae,xnuei,zeff,debye,vexb,vpar
      INTEGER :: is
!
! set flow control switch
      new_matrix = .TRUE.
! transfer values
      do is=1,nsm
        taus_in(is) = tsp(is)
        as_in(is) = asp(is)
      enddo
      betae_in = betae
      xnuei_in = xnuei
      zeff_in = zeff
      debye_in = debye
      vexb_in = vexb
      vpar_in = vpar
!      
      END SUBROUTINE put_averages
!
!-----------------------------------------------------------------
!
      SUBROUTINE put_switches(iflux,use_bper,use_bpar,use_bisection, &
      ibranch,nmodes,nb_max,nb_min,nxgrid,nky,kygrid_model,xnu_model)
!
      USE tglf_internal_interface
      USE tglf_dimensions
!
      IMPLICIT NONE
      LOGICAL :: iflux,use_bper,use_bpar,use_bisection
      INTEGER :: ibranch,nmodes,nb_max,nb_min
      INTEGER :: nxgrid,nky
      INTEGER :: kygrid_model,xnu_model
!
! validaty checks
! reset to defaults if invlaid
!
      if(nb_max.lt.0.or.nb_max.gt.nb)nb_max=nbasis_max_in
      if(nb_min.lt.0.or.nb_min.gt.nb)nb_min=nbasis_min_in
      if(nb_max.lt.nb_min)nb_max=nb_min
      if(ibranch.lt.-1.or.ibranch.gt.2)ibranch=ibranch_in
      if(nxgrid.lt.1.or.2*nxgrid-1.gt.nxm)nxgrid=nxgrid_in
      if(nmodes.gt.maxmodes)nmodes=nmodes_in
      if(nky.lt.2.or.nky.gt.nkym)nky=nky_in
      if(kygrid_model.lt.0.or.kygrid_model.gt.3)kygrid_model = kygrid_model_in
      if(xnu_model.lt.0.or.xnu_model.gt.3)xnu_model = xnu_model_in

!      write(*,*)nb_max,nb_min,ibranch,nxgrid,nmodes,nky
!
! check for changes and update flow controls
!
!      if(nxgrid.ne.nxgrid_in)new_start = .TRUE.
!      if(nb_max.ne.nbasis_max_in)new_start = .TRUE.
!
! transfer values
!
      iflux_in = iflux
      use_bper_in = use_bper
      use_bpar_in = use_bpar
      use_bisection_in = use_bisection
      ibranch_in = ibranch
      nmodes_in = nmodes
      nbasis_max_in = nb_max
      nbasis_min_in = nb_min
      nxgrid_in = nxgrid
      nky_in = nky
      xnu_model_in = xnu_model
      kygrid_model_in = kygrid_model
!
      END SUBROUTINE put_switches
!
!-----------------------------------------------------------------
!
      SUBROUTINE put_rare_switches(rquench,rpark,rghat,rgchat,rwd_zero,   &
                 rLinsker,rgradB,rx_psi,rsat_rule,ralpha_kx0)
!
      USE tglf_internal_interface
      USE tglf_dimensions
!
      IMPLICIT NONE
      REAl,INTENT(IN) :: rquench,rpark,rghat,rgchat,rwd_zero
      REAl,INTENT(IN) :: rLinsker,rgradB,rx_psi,ralpha_kx0
      INTEGER,INTENT(IN) :: rsat_rule
!
! transfer values
!
      alpha_quench_in = rquench
      park_in = rpark
      ghat_in = rghat
      gchat_in = rgchat
      wd_zero_in = rwd_zero
      Linsker_factor_in = rLinsker
      gradB_factor_in = rgradB
      x_psi_in = rx_psi
      sat_rule_in = rsat_rule
      alpha_kx0_in = ralpha_kx0
!
      END SUBROUTINE put_rare_switches
!
!-----------------------------------------------------------------
!
      SUBROUTINE put_model_parameters(adi_elec,alpha_p,alpha_e,  &
         theta_trap,xnu_fac,debye_fac,etg_fac,filter)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: adi_elec
      REAL,INTENT(IN) :: alpha_p,alpha_e,etg_fac,filter
      REAL,INTENT(IN) :: theta_trap,xnu_fac,debye_fac
!
! check for changes and update flow controls
!
      if(adi_elec .NEQV. adiabatic_elec_in)new_matrix = .TRUE.
!
! transfer values
!
      adiabatic_elec_in = adi_elec
      alpha_p_in = alpha_p
      alpha_e_in = alpha_e
      theta_trapped_in = theta_trap
      xnu_factor_in = xnu_fac
      debye_factor_in = debye_fac
      etg_factor_in = etg_fac
      filter_in = filter
!
      END SUBROUTINE put_model_parameters
!
!-----------------------------------------------------------------
!
      SUBROUTINE put_s_alpha_geometry(rmin,rmaj,q,shat,alpha,xwell, &
        theta0,b_model,ft_model)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN):: b_model,ft_model
      REAL,INTENT(IN) :: rmin,rmaj,q,shat,alpha,theta0,xwell
!
! set geometry type flag for shifted circle
      igeo = 0
! set flow control switch
      new_geometry = .TRUE.
!
! transfer values
!
      rmin_sa = rmin
      rmaj_sa = rmaj
      q_sa = q
      q_unit = q   ! needed for kygrid_model_in = 3
      shat_sa = shat
      alpha_sa = alpha
      xwell_sa = xwell
      theta0_sa = theta0
      b_model_sa = b_model
      ft_model_sa = ft_model
!
! validatiy checks
!
      if(rmin_sa.ge.rmaj_sa)rmin_sa=0.999*rmaj_sa   
       if(ft_model_sa.lt.0.or.ft_model_sa.gt.3)then
         write(*,*)"******* ERROR ft_model_sa invalid *******"  
         ft_model_sa = 1
       endif
!
      END SUBROUTINE put_s_alpha_geometry
!
!-----------------------------------------------------------------
!

      SUBROUTINE put_Miller_geometry(rmin,rmaj,q,q_prime,p_prime,shift, &
       kappa,s_kappa,delta,s_delta)
!
! This routine eliminates the need for subroutine miller_init 
! and the miller.dat input file.
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      REAL,INTENT(IN) :: rmin,rmaj,q,q_prime,p_prime,shift
      REAL,INTENT(IN) :: kappa,s_kappa,delta,s_delta
!
! set geometry type flag for Miller
!
      igeo = 1
!
! set flow control switch
!
      new_geometry = .TRUE.
!
! transfer values
!
      rmin_loc = rmin
      rmaj_loc = rmaj
      q_loc = q
      q_unit = q   ! needed for kygrid_model_in = 3
      p_prime_loc = p_prime
      q_prime_loc = q_prime
      shift_loc = shift
      kappa_loc = kappa
      s_kappa_loc = s_kappa
      delta_loc = delta
      s_delta_loc = s_delta
!
! validatiy checks
!
      if(rmin_loc.ge.rmaj_loc)rmin_loc = 0.999*rmaj_loc     
!
      END SUBROUTINE put_Miller_geometry
!
!-----------------------------------------------------------------
!
      SUBROUTINE put_ELITE_geometry(nc,q,q_prime,p_prime,r_c,z_c,bp_c)
!
! This routine requires having read the data from an ELITE geometry file
! giving R,Z,Bp on a flux surface contour with nc points.
!
      USE tglf_internal_interface
      USE tglf_sgrid
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: nc
      REAL,INTENT(IN) :: q,q_prime,p_prime
      REAL,INTENT(IN) :: r_c(0:max_ELITE),z_c(0:max_ELITE),bp_c(0:max_ELITE)
!
      INTEGER i
!
! set geometry type flag for ELITE
!
      igeo = 2
!
! set flow control switch
!
      new_geometry = .TRUE.
!
! validatiy checks
!
      if(nc.gt.max_ELITE)then
        write(*,*)"error in put_ELITE: number of points exceeds ",max_ELITE
        STOP
      endif    
      if(nc.lt.ms)then
        write(*,*)"error in put_ELITE: number of points less than ",ms
        STOP
      endif    
!
! transfer values
!
      n_ELITE = nc
      q_ELITE = q
      q_prime_ELITE = q_prime
      p_prime_ELITE = p_prime
! note direction change for contour angle
      do i=0,n_ELITE
        R_ELITE(n_ELITE-i) = r_c(i)
        Z_ELITE(n_ELITE-i) = z_c(i)
        Bp_ELITE(n_ELITE-i) = bp_c(i)
      enddo
!
      END SUBROUTINE put_ELITE_geometry
!
!-----------------------------------------------------------------
!  output routines
!-----------------------------------------------------------------
!
      REAL FUNCTION get_growthrate(index)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: index
      INTEGER :: i3
!
      i3 = SIZE(gamma_out)
      if(index.gt.i3)then
        write(*,*)"requested growthrate index out of bounds",i3
        get_growthrate = 0.0
      else
        get_growthrate = gamma_out(index)
      endif
!
      END FUNCTION get_growthrate
!
!-----------------------------------------------------------------
!
      REAL FUNCTION get_frequency(index)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: index
      INTEGER :: i3
!
      i3 = SIZE(freq_out)
      if(index.gt.i3)then
        write(*,*)"requested frequency index is of bounds",i3
        get_frequency = 0.0
      else
        get_frequency = freq_out(index)
      endif
!
      END FUNCTION get_frequency
!
!-----------------------------------------------------------------
!
      REAL FUNCTION get_QL_particle_flux(i1,i2,i3)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i1,i2,i3
      INTEGER :: i4
!
      i4=SIZE(particle_QL_out)
      if(i1*i2*i3.gt.i4)then
        write(*,*)"requested QL particle flux index is of bounds",i4
        get_QL_particle_flux = 0.0
      else
        get_QL_particle_flux = particle_QL_out(i1,i2,i3)
      endif
!
      END FUNCTION get_QL_particle_flux
!
!-----------------------------------------------------------------
!
      REAL FUNCTION get_QL_energy_flux(i1,i2,i3)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i1,i2,i3
      INTEGER :: i4
!
      i4=SIZE(energy_QL_out)
      if(i1*i2*i3.gt.i4)then
        write(*,*)"requested QL energy flux index is of bounds",i3
        get_QL_energy_flux = 0.0
      else
        get_QL_energy_flux = energy_QL_out(i1,i2,i3)
      endif
!
      END FUNCTION get_QL_energy_flux
!
!-----------------------------------------------------------------
!
      REAL FUNCTION get_QL_stress_par(i1,i2,i3)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i1,i2,i3
      INTEGER :: i4
!
      i4=SIZE(stress_par_QL_out)
      if(i1*i2*i3.gt.i4)then
        write(*,*)"requested QL stress_par index is of bounds",i3
        get_QL_stress_par = 0.0
      else
        get_QL_stress_par = stress_par_QL_out(i1,i2,i3)
      endif
!
      END FUNCTION get_QL_stress_par
!
!
!-----------------------------------------------------------------
!
      REAL FUNCTION get_QL_stress_tor(i1,i2,i3)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i1,i2,i3
      INTEGER :: i4
!
      i4=SIZE(stress_tor_QL_out)
      if(i1*i2*i3.gt.i4)then
        write(*,*)"requested QL stress_tor index is of bounds",i3
        get_QL_stress_tor = 0.0
      else
        get_QL_stress_tor = stress_tor_QL_out(i1,i2,i3)
      endif
!
      END FUNCTION get_QL_stress_tor
!
!-----------------------------------------------------------------
!
      REAL FUNCTION get_QL_exchange(i1,i2,i3)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i1,i2,i3
      INTEGER :: i4
!
      i4=SIZE(exchange_QL_out)
      if(i1*i2*i3.gt.i4)then
        write(*,*)"requested QL exchange index is of bounds",i3
        get_QL_exchange = 0.0
      else
        get_QL_exchange = exchange_QL_out(i1,i2,i3)
      endif
!
      END FUNCTION get_QL_exchange
!-----------------------------------------------------------------
!
      REAL FUNCTION get_QL_phi(i1)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i1
      INTEGER :: i3
!
      i3 = SIZE(phi_QL_out)
      if(i1.gt.i3)then
        write(*,*)"requested QL phi index is of bounds",i3
        get_QL_phi = 0.0
      else
        get_QL_phi = phi_QL_out(i1)
      endif
!
      END FUNCTION get_QL_phi
!
!-----------------------------------------------------------------
!
      REAL FUNCTION get_QL_density(i1,i2)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i1,i2
      INTEGER :: i3
!
      i3=SIZE(N_QL_out)
      if(i1*i2.gt.i3)then
        write(*,*)"requested QL_density index is of bounds",i3
        get_QL_density = 0.0
      else
        get_QL_density = N_QL_out(i1,i2)
      endif
!
      END FUNCTION get_QL_density
!-----------------------------------------------------------------
!
      REAL FUNCTION get_QL_temperature(i1,i2)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i1,i2
      INTEGER :: i3
!
      i3=SIZE(T_QL_out)
      if(i1*i2.gt.i3)then
        write(*,*)"requested QL_temperature index is of bounds",i3
        get_QL_temperature = 0.0
      else
        get_QL_temperature = T_QL_out(i1,i2)
      endif
!
      END FUNCTION get_QL_temperature
!-----------------------------------------------------------------
!
      REAL FUNCTION get_phi_bar(i1)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i1
      INTEGER :: i3
!
      i3=SIZE(phi_bar_out)
      if(i1.gt.i3)then
        write(*,*)"requested phi_bar index is of bounds",i3
        get_phi_bar = 0.0
      else
        get_phi_bar = phi_bar_out(i1)
      endif
!
      END FUNCTION get_phi_bar
!-----------------------------------------------------------------
!
      REAL FUNCTION get_v_bar(i1)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i1
      INTEGER :: i3
!
      i3=SIZE(v_bar_out)
      if(i1.gt.i3)then
        write(*,*)"requested v_bar index is of bounds",i3
        get_v_bar = 0.0
      else
        get_v_bar = v_bar_out(i1)
      endif
!
      END FUNCTION get_v_bar
!-----------------------------------------------------------------
!
      REAL FUNCTION get_n_bar(i1,i2)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i1,i2
      INTEGER :: i3
!
      i3=SIZE(n_bar_out)
      if(i1*i2.gt.i3)then
        write(*,*)"requested n_bar index is of bounds",i3
        get_n_bar = 0.0
      else
        get_n_bar = n_bar_out(i1,i2)
      endif
!
      END FUNCTION get_n_bar
!-----------------------------------------------------------------
!
      REAL FUNCTION get_t_bar(i1,i2)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i1,i2
      INTEGER :: i3
!
      i3=SIZE(t_bar_out)
      if(i1*i2.gt.i3)then
        write(*,*)"requested t_bar index is of bounds",i3
        get_t_bar = 0.0
      else
        get_t_bar = t_bar_out(i1,i2)
      endif
!
      END FUNCTION get_t_bar
!-----------------------------------------------------------------
!
      REAL FUNCTION get_gaussian_width()
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
!
      get_gaussian_width = width_in
!
      END FUNCTION get_gaussian_width
!-----------------------------------------------------------------
!
      REAL FUNCTION get_R_unit()
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
!
      get_R_unit = R_unit
!
      END FUNCTION get_R_unit
!-----------------------------------------------------------------
!
      REAL FUNCTION get_q_unit()
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
!
      get_q_unit = q_unit
!
      END FUNCTION get_q_unit
!-----------------------------------------------------------------
!
      REAL FUNCTION get_B2_ave()
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
!
      get_B2_ave = B2_ave_out
!
      END FUNCTION get_B2_ave
!-----------------------------------------------------------------
!
      REAL FUNCTION get_R2_ave()
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
!
      get_R2_ave = R2_ave_out
!
      END FUNCTION get_R2_ave
!-----------------------------------------------------------------
!
      REAL FUNCTION get_RBt_ave()
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
!
      get_RBt_ave = RBt_ave_out
!
      END FUNCTION get_RBt_ave
!-----------------------------------------------------------------
!
      REAL FUNCTION get_ave_wd(n1,n2)
!
      USE tglf_internal_interface
      USE tglf_coeff
!
      IMPLICIT NONE
! 
      INTEGER,INTENT(IN):: n1,n2
!
      get_ave_wd = ave_wd(n1,n2)
!
      END FUNCTION get_ave_wd
!-----------------------------------------------------------------
!
      REAL FUNCTION get_wd_bar(n1)
!
      USE tglf_internal_interface
      USE tglf_coeff
!
      IMPLICIT NONE
! 
      INTEGER,INTENT(IN):: n1
!
      get_wd_bar = wd_bar_out(n1)
!
      END FUNCTION get_wd_bar
!-----------------------------------------------------------------
!
      REAL FUNCTION get_ave_b0(n1,n2)
!
      USE tglf_internal_interface
      USE tglf_coeff
!
      IMPLICIT NONE
! 
      INTEGER,INTENT(IN):: n1,n2
!
      get_ave_b0 = ave_b0(n1,n2)
!
      END FUNCTION get_ave_b0
!-----------------------------------------------------------------
!
      REAL FUNCTION get_particle_flux(i1,i2)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i1,i2
      INTEGER :: i3
!
      i3=SIZE(particle_flux_out)
      if(i1*i2.gt.i3)then
        write(*,*)"requested particle flux index is of bounds",i3
        get_particle_flux = 0.0
      else
        get_particle_flux = particle_flux_out(i1,i2)
      endif
!
      END FUNCTION get_particle_flux
!
!-----------------------------------------------------------------
!
      REAL FUNCTION get_energy_flux(i1,i2)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i1,i2
      INTEGER :: i3
!
      i3=SIZE(energy_flux_out)
      if(i1*i2.gt.i3)then
        write(*,*)"requested energy flux index is of bounds",i3
        get_energy_flux = 0.0
      else
        get_energy_flux = energy_flux_out(i1,i2)
      endif
!
      END FUNCTION get_energy_flux
!-----------------------------------------------------------------
!
      REAL FUNCTION get_stress_par(i1,i2)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i1,i2
      INTEGER :: i3
!
      i3=SIZE(stress_par_out)
      if(i1*i2.gt.i3)then
        write(*,*)"requested stress_par index is of bounds",i3
        get_stress_par = 0.0
      else
        get_stress_par = stress_par_out(i1,i2)
      endif
!
      END FUNCTION get_stress_par
!!-----------------------------------------------------------------
!
      REAL FUNCTION get_stress_tor(i1,i2)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i1,i2
      INTEGER :: i3
!
      i3=SIZE(stress_tor_out)
      if(i1*i2.gt.i3)then
        write(*,*)"requested stress_tor index is of bounds",i3
        get_stress_tor = 0.0
      else
        get_stress_tor = stress_tor_out(i1,i2)
      endif
!
      END FUNCTION get_stress_tor
!-----------------------------------------------------------------
!
      REAL FUNCTION get_exchange(i1,i2)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i1,i2
      INTEGER :: i3
!
      i3=SIZE(exchange_out)
      if(i1*i2.gt.i3)then
        write(*,*)"requested exchange index is of bounds",i3
        get_exchange = 0.0
      else
        get_exchange = exchange_out(i1,i2)
      endif
!
      END FUNCTION get_exchange
!
!-----------------------------------------------------------------
!
      REAL FUNCTION get_n_bar_sum(i1)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i1
      INTEGER :: i3
!
      i3=SIZE(n_bar_sum_out)
      if(i1.gt.i3)then
        write(*,*)"requested n_bar_sum index is of bounds",i3
        get_n_bar_sum = 0.0
      else
        get_n_bar_sum = n_bar_sum_out(i1)
      endif
!
      END FUNCTION get_n_bar_sum
!-----------------------------------------------------------------
!
      REAL FUNCTION get_q_low(i1)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i1
      INTEGER :: i3
!
      i3=SIZE(q_low_out)
      if(i1.gt.i3)then
        write(*,*)"requested energy flux index is of bounds",i3
        get_q_low = 0.0
      else
        get_q_low = q_low_out(i1)
      endif
!
      END FUNCTION get_q_low
!-----------------------------------------------------------------
!
      REAL FUNCTION get_t_bar_sum(i1)
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i1
      INTEGER :: i3
!
      i3=SIZE(t_bar_sum_out)
      if(i1.gt.i3)then
        write(*,*)"requested t_bar_sum index is of bounds",i3
        get_t_bar_sum = 0.0
      else
        get_t_bar_sum = t_bar_sum_out(i1)
      endif
!
      END FUNCTION get_t_bar_sum
!-----------------------------------------------------------------
!
      REAL FUNCTION get_phi_bar_sum()
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
!
      get_phi_bar_sum = phi_bar_sum_out
!
      END FUNCTION get_phi_bar_sum
!-----------------------------------------------------------------
!
      REAL FUNCTION get_v_bar_sum()
!
      USE tglf_internal_interface
!
      IMPLICIT NONE
!
      get_v_bar_sum = v_bar_sum_out
!
      END FUNCTION get_v_bar_sum
!-----------------------------------------------------------------
!

