MODULE common_constants     ! Define Universal Constants
  USE nrtype,                      ONLY : DP,I4B
  IMPLICIT NONE

  ! ========== Mathematical Constants ==========
 
  REAL(DP), PARAMETER:: Deg_Per_Rad  = 57.295779513082320876798155_DP
  REAL(DP), PARAMETER:: Rad_Per_Deg  = 0.017453292519943295769237_DP

  REAL(DP), PARAMETER:: e_Value      =  2.71828182845904523560287_DP
  REAL(DP), PARAMETER:: e_Recip      =  0.3678794411714423215955238_DP
  REAL(DP), PARAMETER:: e_Squared    =  7.389056098930650227230427_DP
  REAL(DP), PARAMETER:: Log10_of_e   =  0.4342944819032518276511289_DP

  REAL(DP), PARAMETER:: Euler        =  0.5772156649015328606_DP
  REAL(DP), PARAMETER:: Euler_Log    = -0.5495393129816448223_DP
  REAL(DP), PARAMETER:: Gamma        =  0.577215664901532860606512_DP
  REAL(DP), PARAMETER:: Gamma_Log    = -0.549539312981644822337662_DP
  REAL(DP), PARAMETER:: Golden_Ratio =  1.618033988749894848_DP

  REAL(DP), PARAMETER:: Ln_2         =  0.6931471805599453094172321_DP
  REAL(DP), PARAMETER:: Ln_10        =  2.3025850929940456840179915_DP
  REAL(DP), PARAMETER:: Log10_of_2   =  0.3010299956639811952137389_DP


  REAL(DP), PARAMETER :: SQRT2    = 1.41421356237309504880168872420969807856967_DP

  REAL(DP), PARAMETER :: PI       = 3.141592653589793238462643383279502884197_DP
  REAL(DP), PARAMETER :: PIO2     = 1.57079632679489661923132169163975144209858_DP
  REAL(DP), PARAMETER :: PISQ     = 9.869604401089358618834491_DP
  REAL(DP), PARAMETER :: TWOPI    = 6.283185307179586476925286766559005768394_DP
  REAL(DP), PARAMETER :: ROOTPI   = 1.772453850905516027298167_DP
  REAL(DP), PARAMETER :: ROOTPIO2 = 1.253314137315500251207883_DP
  REAL(DP), PARAMETER :: ROOT2PI  = 2.506628274631000502415765_DP
  REAL(DP), PARAMETER :: ZEROC    = 0.0_DP
  INTEGER(I4B),PARAMETER :: IZERO = 0_I4B

  ! ========== Physics Constants and units ==========

  REAL(DP), PARAMETER :: AMU_Value          = 1.6605402E-27_DP   ! kg
  REAL(DP), PARAMETER :: Atmosphere_Pres    = 9.80665E+04_DP     ! Pa
  REAL(DP), PARAMETER :: Avogadro           = 6.0221367E+23_DP   ! 1/mol
  REAL(DP), PARAMETER :: Boltzmann          = 1.380657E-23_DP    ! J/K
  REAL(DP), PARAMETER :: c_Light            = 2.997924580E+8_DP  ! m/s
  REAL(DP), PARAMETER :: Electron_Charge    =-1.60217738E-19_DP  ! coul
  REAL(DP), PARAMETER :: Electron_Rest_Mass = 9.1093897E-31_DP   ! kg
  REAL(DP), PARAMETER :: Faraday            = 9.6485309E+04_DP   ! C/mo
  REAL(DP), PARAMETER :: Neutron_Mass       = 1.6749286E-27_DP   ! kg
  REAL(DP), PARAMETER :: Permeability       = 1.25663706143E-06_DP ! H/m
  REAL(DP), PARAMETER :: Permittivity       = 8.85418781762E-12_DP ! F/m
  REAL(DP), PARAMETER :: Planck_Const       = 6.6260754E-34_DP   ! J*s
  REAL(DP), PARAMETER :: Proton_Mass        = 1.6726230E-27_DP   ! kg
  REAL(DP), PARAMETER :: Stefan_Boltzmann   = 5.67050E-08_DP     ! W/(m^2*K^4)
  REAL(DP), PARAMETER :: Thomson_cross_sect = 6.6516E-29_DP      ! m^2
  REAL(DP), PARAMETER :: Universal_Gas_C    = 8.314510_DP        ! J/mol*K

! conversion constants
  REAL(DP), PARAMETER :: g2kg               = 0.001_DP           ! kg/g
  REAL(DP), PARAMETER :: M2cm               = 100._DP            ! cm/m
  REAL(DP), PARAMETER :: M32cm3             = 1.e6_DP            ! m^3 to cm^3
  REAL(DP), PARAMETER :: cm2M               = 0.01_DP            ! m/cm
  REAL(DP), PARAMETER :: M22cm2             = 1.e4_DP            ! cm2/m2
  REAL(DP), PARAMETER :: im32icm3           = 1.e-6_DP           ! m3/cm3
  REAL(DP), PARAMETER :: im22icm2           = 1.e-4_DP           ! m2/cm2
  REAL(DP), PARAMETER :: T2gauss            = 1.e4_DP            ! Gauss/Tesla
  REAL(DP), PARAMETER :: Tm2gcm             = 1.e-6_DP           ! Tesla-m/ Gauss cm
  REAL(DP), PARAMETER :: vs2kgcm2           = 1.e5_DP            ! vols sec to kgauss cm^2
  REAL(DP), PARAMETER :: convert            = 1.60217733e-10_DP  ! keV/sec/cm*3 > watts/m**3
  REAL(DP), PARAMETER :: kgauss2t           = 0.1_DP             ! Tesla/kilogauss
  REAL(DP), PARAMETER :: pconvert            = 0.1_DP             ! pascal/(dyne.cm^2)
  REAL(DP), PARAMETER :: ppmks              = 1.e-7              ! pprime convesion to mks
  REAL(DP), PARAMETER :: ffpmks             = 1.e-4              ! ff prime convesion to mks                                                     !
  REAL(DP), PARAMETER :: kevperg            = 6.2415097e+08      ! Kev/erg
  REAL(DP), PARAMETER :: joupkev            = 1.6021765e-16      ! jou/Kev
  REAL(DP), PARAMETER :: kevpjou            = 6.2415097e+15      ! Kev/jou

END MODULE common_constants
 
