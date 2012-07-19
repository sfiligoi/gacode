!-------------------------------------------------------------------
! prgen.f90
!
! PURPOSE:
!  Generate input.profiles from prexisting data.  This executable is
!  managed and called by gacode/shared/bin/profiles_gen.
!
!  Acceptable input formats:
!
!  1. ONETWO iterdb text (no extenstion)  
!  2. ONETWO iterdb NetCDF (*.nc)
!  3. PEQDSK text (*.peq)
!  4. PLASMA STATE NetCDF (*.cdf)
!  5. CORSICA text (*.corsica)
!
! Note that ASTRA (*.astra) format is handled by a separate python 
! routine.
!--------------------------------------------------------------------

program prgen

  use prgen_read_globals

  !--------------------------------------------------
  implicit none
  !
  logical :: back=0
  !--------------------------------------------------

  !--------------------------------------------------
  ! ** Define number of tags here: 
  !
  n_indx  = 40
  n_indx2 = 15  
  !--------------------------------------------------

  !--------------------------------------------------
  ! Parse the config file
  !
  open(unit=1,file='.config',status='old')
  read(1,'(a)') date
  read(1,'(a)') raw_data_file
  read(1,'(a)') cer_file
  read(1,*) gato_flag
  read(1,*) nogatoq_flag
  read(1,*) verbose_flag
  read(1,*) reorder_vec(:)

  close(1)
  !--------------------------------------------------

  !--------------------------------------------------
  ! Define tags
  !
  allocate(tag(n_indx))
  !
  tag(1)  = 'rho(-)'
  tag(2)  = 'rmin(m)'
  tag(3)  = 'rmaj(m)'
  tag(4)  = 'q(-)'
  tag(5)  = 'kappa(-)'
  tag(6)  = 'delta(-)'
  tag(7)  = 'Te(keV)'
  tag(8)  = 'ne(10^19/m^3)'
  tag(9)  = 'z_eff(-)'
  tag(10) = 'omega0(1/s)'
  tag(11) = 'flow_mom(N-m)'
  tag(12) = 'pow_e(MW)'
  tag(13) = 'pow_i(MW)'
  tag(14) = 'pow_ei(MW)'
  tag(15) = 'zeta(-)'
  tag(16) = 'flow_beam(kW/eV)'
  tag(17) = 'flow_wall(kW/eV)'
  tag(18) = 'zmag(m)'
  tag(19) = 'ptot(Pa)'
  tag(20) = 'polflux(Wb/rad)'
  tag(21) = 'ni_1(10^19/m^3)'
  tag(22) = 'ni_2(10^19/m^3)'
  tag(23) = 'ni_3(10^19/m^3)'
  tag(24) = 'ni_4(10^19/m^3)'
  tag(25) = 'ni_5(10^19/m^3)'
  tag(26) = 'Ti_1 (keV)'
  tag(27) = 'Ti_2 (keV)'
  tag(28) = 'Ti_3 (keV)'
  tag(29) = 'Ti_4 (keV)'
  tag(30) = 'Ti_5 (keV)'
  tag(31) = 'vtor_1 (m/s)'
  tag(32) = 'vtor_2 (m/s)'
  tag(33) = 'vtor_3 (m/s)'
  tag(34) = 'vtor_4 (m/s)'
  tag(35) = 'vtor_5 (m/s)'
  tag(36) = 'vpol_1 (m/s)'
  tag(37) = 'vpol_2 (m/s)'
  tag(38) = 'vpol_3 (m/s)'
  tag(39) = 'vpol_4 (m/s)'
  tag(40) = 'vpol_5 (m/s)'
  !
  ! Transport power components
  !  
  allocate(tag2(n_indx2))
  !
  tag2(1) = 'powe_beam(MW)' 
  tag2(2) = 'powe_RF(MW)'
  tag2(3) = 'powe_oh_RF(MW)'
  tag2(4) = 'powe_rad_RF(MW)'
  tag2(5) = 'powe_ion(MW)'
  tag2(6) = 'powe_wdot(MW)'
  tag2(7) = 'powe_fus(MW)'
  tag2(8) = '[tr-pow_e]'
  tag2(9) = '[ex-pow_ei_exp]' 
  tag2(10) = '[tr-pow_i]'
  tag2(11) = 'powi_beam(MW)'
  tag2(12) = 'powi_ion(MW)'
  tag2(13) = 'powi_wdot(MW)'
  tag2(14) = 'powi_fus(MW)'
  tag2(15) = 'powi_cx(MW)'
  !--------------------------------------------------


  !------------------------------------------------------------------
  ! Read the iterdb file and define standard variables.
  !
  ! Note that nx will be the experimental vector length in ALL cases:
  !
  if (trim(raw_data_file) == 'null') then
 
    ! Pure gfile parsing

    format_type = 0

    ! Minimal processing required to merge gfile data into otherwise
    ! empty input.profiles output file.
    call prgen_read_null

  else if (index(raw_data_file,'.nc',back) /= 0) then

     ! New NetCDF format
     print '(a)','INFO: (prgen) Assuming iterdb NetCDF format.'

     format_type = 1

     call prgen_read_iterdb_nc

  else if (index(raw_data_file,'.cdf',back) /= 0) then

     ! Plasmastate format
     print '(a)','INFO: (prgen) Assuming plasma_state format.'

     format_type = 2

     call prgen_read_plasmastate

  else if (index(raw_data_file,'.peq',back) /= 0) then

     ! peqdsk format
     print '(a)','INFO: Assuming peqdsk format.'

     format_type = 3

     if (gato_flag /= 1) then
        print '(a)','ERROR: (prgen) geqdsk must be provided for peqdsk format'
        stop
     endif

     call prgen_read_peqdsk

  else if (index(raw_data_file,'.corsica',back) /= 0) then

     ! corsica format
     print '(a)','INFO: Assuming corsica format.'

     format_type = 5

     if (gato_flag /= 1) then
        print '(a)','WARNING: (prgen) geqdsk must be provided for corsica format'
     endif

     call prgen_read_corsica

  else if (index(raw_data_file,'.ufile',back) /= 0) then

     ! UFILE format
     print '(a)','INFO: Assuming UFILE format.'

     format_type = 6

     call prgen_read_ufile

  else

     ! Old text format
     print '(a)','INFO: (prgen) Assuming old iterdb text format.'

     format_type = 1

     call prgen_read_iterdb

  endif
  !------------------------------------------------------------------

  !---------------------------------------------------
  ! Read the GATO file for "better" geometry.  At this
  ! point, GATO has already run and we are just reading 
  ! the output.
  !
  if (gato_flag == 1) call prgen_read_gato
  !---------------------------------------------------

  if (format_type == 0) then
     call prgen_map_null
  else if (index(raw_data_file,'.cdf',back) /= 0) then
     call prgen_map_plasmastate
  else if (index(raw_data_file,'.peq',back) /= 0) then
     call prgen_map_peqdsk
  else if (index(raw_data_file,'.corsica',back) /= 0) then
     call prgen_map_corsica
  else
     call prgen_map_iterdb
  endif

  call prgen_write

  ! Successful completion
  open(unit=1,file='success',status='replace')
  write(1,*) 1
  close(1)

end program prgen
