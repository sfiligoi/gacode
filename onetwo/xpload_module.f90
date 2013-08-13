module xpload_module
!
!  module to store some xplasma id information
!
  integer,save :: id_rho         ! standard flux surface grid id
  integer,save :: id_chi         ! poloidal angle grid id
  integer,save :: id_rhozc       ! grid for zone centered functions
!
  integer,save :: id_Rpol,id_Zpol ! ids for R(rho,chi), Z(rho,chi)
!
  integer,save :: id_R,id_Z      ! ids for R,Z grids (e.g. for fast inv. map)
!
  integer,save :: id_br, id_bz, id_bphi ! ids for B field components
!
  integer,save :: id_g           ! id for 1D profile of 'G' = R*Bt (T*m)
  integer,save :: id_psi         ! id for 1D profile of poloidal flux 'PSI_ORB' (Webers)
  integer,save :: id_dpsi        ! id for 1D profile of derivative of poloidal flux
                                 ! 'DPSI_ORB' (Webers)
  integer,save :: id_q           ! id for 1D profile of 'Q_ORB'
  integer,save :: id_phi         ! id for 1D profile of electric potential 'PHI_ORB' (volts)
  integer,save :: id_bmod        ! id for internal BMOD 2D profile -- 'BMOD' (T)
  integer,save :: id_invbmod     ! id for 1./BMOD 2D profile -- 'INVBMOD_ORB' (1/T)
  integer,save :: id_detjacob    ! id for the jacobian determinant 2D profile
                                 ! 'DET_JACOB_ORB' (m**3)
  integer,save :: id_detjacobpsi ! id for the jacobian/d(psi)/dxi 2D profile
                                 ! 'DET_JACOB_P_ORB' (m**3/Weber=m/tesla)
  integer,save :: id_mcgrid      !HSJ added
!
end module xpload_module
