  MODULE io_gcnmp
    USE nrtype,   ONLY : I2B,I4B,DP
    INTEGER(I2B), PARAMETER :: ncrt      = 6    ! output to screen
    INTEGER(I2B), PARAMETER :: io_repeat = 9    ! plotfile handling
    INTEGER(I2B), PARAMETER :: ioplot    = 700  ! plotfile handling 88888899999
    !INTEGER(I2B), PARAMETER :: maxchr   = 132
    INTEGER(I2B), PARAMETER :: maxchr    = 256  ! HSJ 2/6/13
    INTEGER(I2B)  niterdb,nout,nlog,switch_iterdb_output,iot,io_fd
    INTEGER(I4B)  rpc
    INTEGER(I4B)  lun_nubeam

    CHARACTER(maxchr)                                              &
                  output_filename,namelist_filename,               &
                  runlog_filename,iterdb_12_input_file,            &
                  statefile_output_name,fixt_filename

    LOGICAL switch_statefile_output
 
    ! NOTE switch_statefile_output is new logical version of 
    ! switch_iterdb_output. (P_nfreya uses switch_statefile_output)

  END MODULE io_gcnmp
