module prgen_read_globals

  ! Parameters
  real, parameter :: kevdsecpmw=1.6022e-19*1e3*1e-6
  real, parameter :: pi=3.14159265358979323846
  integer, parameter :: nfourier=8
  
  ! List of possible data files
  character (len=70) :: date
  character (len=70) :: raw_data_file
  character (len=70) :: cer_file

  !----------------------------------------------------------
  ! Internal control variables
  !
  integer :: nx
  integer :: n_indx
  integer :: n_indx2
  integer :: extra_powers_flag = 0
  integer :: format_type
  integer :: gato_flag
  integer :: nogatoq_flag
  integer :: verbose_flag
  character (len=16), dimension(:), allocatable :: tag
  character (len=16), dimension(:), allocatable :: tag2
  character (len=70) :: efit_header
  integer, dimension(5) :: reorder_vec
  real :: dpsi_gato
  real :: dpsi_data
  real :: signpsi
  !----------------------------------------------------------

  ! Internal physics variables
  real, dimension(:), allocatable :: rho
  real, dimension(:), allocatable :: dpsi
  real, dimension(:), allocatable :: rmin
  real, dimension(:), allocatable :: rmaj
  real, dimension(:), allocatable :: q
  real, dimension(:), allocatable :: q_gato
  real, dimension(:), allocatable :: vphi_exp
  real, dimension(:), allocatable :: vphi_carbon
  real, dimension(:), allocatable :: vpolc_exp
  real, dimension(:), allocatable :: vtorc_exp
  real, dimension(:), allocatable :: omega0
  real, dimension(:), allocatable :: kappa
  real, dimension(:), allocatable :: delta
  real, dimension(:), allocatable :: zmag
  real, dimension(:), allocatable :: zeta
  real, dimension(:), allocatable :: sion_d
  real, dimension(:), allocatable :: sbcx_d
  real, dimension(:), allocatable :: powe_beam_exp
  real, dimension(:), allocatable :: powi_beam_exp
  real, dimension(:), allocatable :: powe_rf_exp
  real, dimension(:), allocatable :: powi_rf_exp
  real, dimension(:), allocatable :: powe_oh_exp
  real, dimension(:), allocatable :: powe_rad_exp
  real, dimension(:), allocatable :: powe_ion_exp
  real, dimension(:), allocatable :: powi_ion_exp
  real, dimension(:), allocatable :: powe_wdot_exp
  real, dimension(:), allocatable :: powi_wdot_exp
  real, dimension(:), allocatable :: powe_fus_exp
  real, dimension(:), allocatable :: powi_fus_exp
  real, dimension(:), allocatable :: pow_ei_exp
  real, dimension(:), allocatable :: pow_e
  real, dimension(:), allocatable :: pow_i
  real, dimension(:), allocatable :: powi_cx_exp
  real, dimension(:), allocatable :: flow_wall_exp
  real, dimension(:), allocatable :: flow_beam
  real, dimension(:), allocatable :: flow_mom
  real, dimension(:,:), allocatable :: vec
  real, dimension(:,:), allocatable :: vec2

  character (len=70), dimension(5) :: ion_name

  !---------------------------------------------------------
  ! ONETWO variables
  !
  integer :: onetwo_ishot
  integer :: onetwo_nj
  integer :: onetwo_nb
  integer :: onetwo_npsi
  integer :: onetwo_nion
  integer :: onetwo_nprim
  integer :: onetwo_nimp
  integer :: onetwo_nneu
  integer :: onetwo_nbion
  character (len=2), dimension(5) :: onetwo_namep
  character (len=2), dimension(5) :: onetwo_namei
  character (len=2), dimension(5) :: onetwo_nameb
  real :: onetwo_time
  real :: onetwo_Rgeom
  real :: onetwo_Rmag
  real :: onetwo_R0
  real :: onetwo_kappa
  real :: onetwo_delta
  real :: onetwo_volo
  real :: onetwo_cxareao
  real :: onetwo_Btor
  !
  real, dimension(:), allocatable :: onetwo_rho_grid
  real, dimension(:), allocatable :: onetwo_volume
  real, dimension(:), allocatable :: onetwo_hcap
  real, dimension(:), allocatable :: onetwo_qbeame
  real, dimension(:), allocatable :: onetwo_qrfe
  real, dimension(:), allocatable :: onetwo_qohm
  real, dimension(:), allocatable :: onetwo_qrad
  real, dimension(:), allocatable :: onetwo_qione
  real, dimension(:), allocatable :: onetwo_qioni
  real, dimension(:), allocatable :: onetwo_dpedt
  real, dimension(:), allocatable :: onetwo_qfuse
  real, dimension(:), allocatable :: onetwo_qbeami
  real, dimension(:), allocatable :: onetwo_qrfi
  real, dimension(:), allocatable :: onetwo_qcx
  real, dimension(:), allocatable :: onetwo_dpidt
  real, dimension(:), allocatable :: onetwo_qfusi
  real, dimension(:), allocatable :: onetwo_qdelt
  real, dimension(:), allocatable :: onetwo_sbeam
  real, dimension(:), allocatable :: onetwo_te
  real, dimension(:), allocatable :: onetwo_ti
  real, dimension(:), allocatable :: onetwo_q
  real, dimension(:), allocatable :: onetwo_angrot
  real, dimension(:), allocatable :: onetwo_zeff
  real, dimension(:), allocatable :: onetwo_ene
  real, dimension(:), allocatable :: onetwo_psi
  real, dimension(:), allocatable :: onetwo_storqueb
  real, dimension(:), allocatable :: onetwo_press
  real, dimension(:,:), allocatable :: onetwo_enion
  real, dimension(:,:), allocatable :: onetwo_enbeam
  real, dimension(:,:), allocatable :: onetwo_pressb
  real, dimension(:,:), allocatable :: onetwo_sion
  real, dimension(:,:), allocatable :: onetwo_srecom
  real, dimension(:,:), allocatable :: onetwo_scx
  real, dimension(:,:), allocatable :: onetwo_sbcx
  real, dimension(:,:), allocatable :: onetwo_s
  real, dimension(:,:), allocatable :: onetwo_dudt
  !
  ! psi-grid vectors
  !
  real, dimension(:), allocatable :: onetwo_rho_mhd_gridnpsi
  real, dimension(:), allocatable :: onetwo_rmajavnpsi
  real, dimension(:), allocatable :: onetwo_rminavnpsi
  real, dimension(:), allocatable :: onetwo_psivolpnpsi
  real, dimension(:), allocatable :: onetwo_elongxnpsi
  real, dimension(:), allocatable :: onetwo_triangnpsi_u
  real, dimension(:), allocatable :: onetwo_triangnpsi_l
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! PLASMA STATE variables
  !
  integer :: plst_shot_number
  character (len=32) :: plst_tokamak_id
  integer :: plst_dim_nrho
  integer :: plst_dim_nrho_eq_geo
  integer :: plst_dp1_nspec_th
  integer :: plst_dp1_nspec_tha
  integer :: plst_dp1_nspec_all
  integer :: plst_dp1_nspec_alla
  character (len=32) :: plst_tag
  character (len=32), dimension(:), allocatable :: plst_all_name
  character (len=32), dimension(:), allocatable :: plst_alla_name
  real :: plst_b_axis_vac
  integer :: plst_kccw_bphi
  integer :: plst_kccw_jphi
  real, dimension(:,:), allocatable :: plst_ts
  real, dimension(:,:), allocatable :: plst_ns
  real, dimension(:), allocatable :: plst_ptowb
  real, dimension(:), allocatable :: plst_nb
  real, dimension(:), allocatable :: plst_eperpb
  real, dimension(:), allocatable :: plst_eparb
  real, dimension(:), allocatable :: plst_vol
  real, dimension(:), allocatable :: plst_rho
  real, dimension(:), allocatable :: plst_grho1
  real, dimension(:), allocatable :: plst_phit
  real, dimension(:), allocatable :: plst_psipol
  real, dimension(:), allocatable :: plst_elong
  real, dimension(:), allocatable :: plst_triang
  real, dimension(:), allocatable :: plst_iota
  real, dimension(:), allocatable :: plst_r_midp_in
  real, dimension(:), allocatable :: plst_r_midp_out
  real, dimension(:), allocatable :: plst_z_midp
  real, dimension(:), allocatable :: plst_zeff
  real, dimension(:), allocatable :: plst_epot
  real, dimension(:), allocatable :: plst_omegat
  real, dimension(:), allocatable :: plst_pbe
  real, dimension(:), allocatable :: plst_pbi
  real, dimension(:), allocatable :: plst_pe_trans
  real, dimension(:), allocatable :: plst_pi_trans
  real, dimension(:), allocatable :: plst_pei_trans
  real, dimension(:), allocatable :: plst_tq_trans
  real, dimension(:), allocatable :: plst_sn_trans
  !---------------------------------------------------------
 
  !---------------------------------------------------------
  ! PEQDSK variables
  !
  integer :: peqdsk_nj
  real :: peqdsk_bref
  real :: peqdsk_arho
  real, dimension(:), allocatable :: peqdsk_psi
  real, dimension(:), allocatable :: peqdsk_ne
  real, dimension(:), allocatable :: peqdsk_te
  real, dimension(:), allocatable :: peqdsk_ni
  real, dimension(:), allocatable :: peqdsk_ti
  real, dimension(:), allocatable :: peqdsk_omegat
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! null input variables
  !
  real :: null_bref
  real :: null_arho
  !---------------------------------------------------------

  !--------------------------------------------------------
  ! Some GATO variables
  !
  integer :: gato_npsi
  integer :: gato_ntheta
  !--------------------------------------------------------

  !---------------------------------------------------------
  ! CORSICA variables
  !
  real :: corsica_time
  real :: corsica_current
  real :: corsica_loop_voltage
  real :: corsica_fusion_power
  real :: corsica_bootstrap_fraction
  real :: corsica_betat
  real :: corsica_betan
  real :: corsica_taue
  real :: corsica_q95
  real :: corsica_r
  real :: corsica_a
  real :: corsica_kappa
  real :: corsica_delta
  real :: corsica_fusion_gain
  real :: corsica_betap
  real :: corsica_li3

  integer, parameter :: corsica_nvals = 81
  real :: corsica_bref
  real :: corsica_arho
  real, allocatable :: corsica_rho(:)
  real, allocatable :: corsica_r_a(:)
  real, allocatable :: corsica_psin(:)
  real, allocatable :: corsica_vl(:)
  real, allocatable :: corsica_te(:)
  real, allocatable :: corsica_ti(:)
  real, allocatable :: corsica_ne(:)
  real, allocatable :: corsica_ndt(:)
  real, allocatable :: corsica_nz(:)
  real, allocatable :: corsica_nalpha(:)
  real, allocatable :: corsica_zeff(:)
  real, allocatable :: corsica_q(:)
  real, allocatable :: corsica_j(:)
  real, allocatable :: corsica_jbs(:)
end module prgen_read_globals
