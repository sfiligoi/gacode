    MODULE gcnmp_input

      USE nrtype,                   ONLY : I4B,Dp
      USE param,                    ONLY : kk
      USE solcon_gcnmp,             ONLY : std_file_size
      IMPLICIT NONE
      INTEGER(I4B)                         eq_split_gcnmp,                &
                                           save_incr_restart_file_gcnmp
      INTEGER(I4B)                         itran_gcnmp(kk),gcnmp_nprocs,  &
                                           switch_iterdb_output
      REAL(DP)                             gcnmp_macro_dt,dtmin_gcnmp
      CHARACTER(LEN = std_file_size)       iterdbf,gcnmp_host,            &
                                           gcnmp_nml_filename,            &
                                           gcnmp_remote_dir,              &
                                           gcnmp_iterdb_filename

      LOGICAL                              write_iterdb_txt

    END MODULE gcnmp_input
