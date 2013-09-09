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
!  3. PEQDSK text (*.peq, *.peq2)
!  4. PLASMA STATE NetCDF (*.cdf, *.CDF)
!  5. CORSICA text (*.corsica)
!  6. UFILE text (UFILE)
!
! Note that ASTRA (*.astra) format is handled by a separate python 
! routine.
!--------------------------------------------------------------------

program prgen

  use prgen_globals
  use EXPRO_interface

  !--------------------------------------------------
  implicit none
  !--------------------------------------------------

  !--------------------------------------------------
  ! Parse the config file
  !
  open(unit=1,file='.config',status='old')
  read(1,'(a)') date
  read(1,'(a)') raw_data_file
  read(1,'(a)') cer_file
  read(1,*) efit_method
  read(1,*) nogatoq_flag
  read(1,*) verbose_flag
  read(1,*) gmerge_flag
  read(1,*) ipccw
  read(1,*) btccw
  read(1,*) nfourier
  read(1,*) reorder_vec(:)
  close(1)
  !--------------------------------------------------

  n_indx  = size(EXPRO_tag) 

  !------------------------------------------------------------------
  ! Read the iterdb file and define standard variables.
  !
  ! Note that nx will be the experimental vector length in ALL cases:
  !
  if (gmerge_flag == 1) then

     call prgen_read_inputprofiles

     format_type = 7

  else if (trim(raw_data_file) == 'null') then

     ! Pure gfile parsing

     format_type = 0

     ! Minimal processing required to merge gfile data into otherwise
     ! empty input.profiles output file.
     call prgen_read_null

  else if (index(raw_data_file,'.nc') /= 0) then

     ! New NetCDF format
     print '(a)','INFO: (prgen) Assuming iterdb NetCDF format.'

     format_type = 1

     call prgen_read_iterdb_nc

  else if (index(raw_data_file,'.cdf') /= 0 .or. index(raw_data_file,'.CDF') /= 0) then

     ! Plasmastate format
     print '(a)','INFO: (prgen) Assuming plasma_state format.'

     format_type = 2

     call prgen_read_plasmastate

  else if (index(raw_data_file,'.peq') /= 0) then

     ! peqdsk format
     print '(a)','INFO: (prgen) Assuming peqdsk format.'

     format_type = 3

     if (efit_method == 0) then
        print '(a)','ERROR: (prgen) geqdsk must be provided for peqdsk format'
        stop
     endif

     call prgen_read_peqdsk

  else if (index(raw_data_file,'.corsica') /= 0) then

     ! corsica format
     print '(a)','INFO: (prgen) Assuming corsica format.'

     format_type = 5

     if (efit_method == 0) then
        print '(a)','WARNING: (prgen) geqdsk must be provided for corsica format'
     endif

     call prgen_read_corsica

  else if (index(raw_data_file,'UFILE') /= 0) then

     ! UFILE format
     print '(a)','INFO: (prgen) Assuming UFILE format.'

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
  select case (efit_method)
  case (1)
     ! Use geometry data contained in profile data 
  case (2)
     ! Use GATO-EFIT mapper
     call prgen_read_gato
  case (3)
     ! Use OMFIT-EFIT mapper
     call prgen_read_omfit
  end select
  !---------------------------------------------------

  select case (format_type)

  case (0)
     call prgen_map_null
  case (1)
     call prgen_map_iterdb
  case (2) 
     call prgen_map_plasmastate
  case (3) 
     call prgen_map_peqdsk
  case (5) 
     call prgen_map_corsica
  case (6)
     call prgen_map_ufile
  case (7)
     call prgen_map_inputprofiles

  end select

  call prgen_write

  ! Successful completion
  open(unit=1,file='success',status='replace')
  write(1,*) 1
  close(1)

end program prgen
