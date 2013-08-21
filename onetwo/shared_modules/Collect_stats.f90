  MODULE statistics
    !------------------------------------------------------------------------------------
    ! collect information about code runs.
    ! At present only certain execution times are monitored
    !--------------------------------------------------------HSJ--------------------------

    USE nrtype,                             ONLY : DP,I4B,SP

    IMPLICIT NONE
    INTEGER, PARAMETER :: nstats = 5        ! no of different types of stats
                                            ! to be collected. EG nfreya run time,
                                            ! nubeam run times, time step calcualtions,etc


    LOGICAL   p_nf_mon_set,nubeam_mon_set,ech_mon_set,o12_mon_set,tdem_mon_set
    LOGICAL,  DIMENSION(nstats) :: start_timer,stop_timer,initialize
    REAL(SP), DIMENSION(nstats) :: elapsed_time
    REAL(SP)  elsv,clock_rate_sv
    INTEGER(I4B), DIMENSION(nstats) :: time_start,time_end,imon
    INTEGER(I4B)   :: clock_start,clock_end,clock_rate,last_mon_index,max_ct,cycle
    INTEGER(I4B)   :: p_nf_index,nubeam_index,ech_index,o12_index,tdem_index
    CHARACTER(len = 64), DIMENSION(nstats) :: descrip_mon    ! description of data to be collected
                                                             ! defined in main code as necessary

    DATA initialize(1:nstats) /nstats*.TRUE./

    DATA last_mon_index /0/
    DATA cycle /0/

    DATA p_nf_index /0/
    DATA nubeam_index /0/
    DATA ech_index /0/
    DATA o12_index /0/
    DATA tdem_index /0/

    DATA p_nf_mon_set /.FALSE./
    DATA nubeam_mon_set /.FALSE./
    DATA ech_mon_set /.FALSE./
    DATA o12_mon_set /.FALSE./
    DATA tdem_mon_set /.FALSE./


    CONTAINS

      SUBROUTINE collect_stats(stat_index)
      !-------------------------------------------------------------------------------
      ! stat_index points to particular event times defined by the calling code
      !-------------------------------------------------------------------------------
        USE error_handler,                          ONLY : terminate,lerrno,iomaxerr
        USE io_gcnmp,                               ONLY : nlog,ncrt


        INTEGER(I4B) stat_index


        IF(start_timer(stat_index)  .AND. stop_timer(stat_index))THEN
            WRITE(ncrt, FMT='("Error: sub collect_stats, start and stop on simultaneously for ",a)')descrip_mon(stat_index)
            lerrno = 178 + iomaxerr
            CALL terminate(lerrno,nlog)
        ENDIF
        IF(.NOT. start_timer(stat_index)  .AND. .NOT. stop_timer(stat_index))THEN
            WRITE(ncrt, FMT='("Error: sub collect_stats, start and stop off simultaneously for ",a)')descrip_mon(stat_index)
            lerrno = 178 + iomaxerr
            CALL terminate(lerrno,nlog)
        ENDIF

        IF(stat_index .LT. 1 .OR. stat_index .GT. nstats)THEN
            WRITE(ncrt, FMT='("Error: sub collect_stats,stat_index = ",i4)')stat_index
            WRITE(ncrt, FMT='("Error: sub collect_stats,nstats parameter  = ",i4)')nstats
            lerrno = 179 + iomaxerr
            CALL terminate(lerrno,nlog)
        ENDIF

  if(o12_index .Gt. 1) THEN
!     write(888,fmt='("o12 index =",i5)')o12_index
     call stop("o12 index prob",1)
  ENDIF

        IF(start_timer(stat_index) .AND. .NOT. initialize(stat_index) )THEN
           CALL SYSTEM_CLOCK(COUNT=clock_start) ! Start timing
           time_start(stat_index) = clock_start
!           write(888,FMT='("time_start (non init) = ",x,i10)') time_start(stat_index)
!           write(888,FMT='(" event = ",a)')descrip_mon(stat_index)
        ELSEIF(start_timer(stat_index) .AND.  initialize(stat_index) )THEN
          ! Get the clock rate on the machine we are running on.
          ! We will assume that the rate is constant, not changing
          ! due to load,temperature,etc.
           CALL SYSTEM_CLOCK(COUNT=clock_start,COUNT_RATE=clock_rate,COUNT_MAX = max_ct) 
           clock_rate_sv = REAL(clock_rate,sp)
           time_start(stat_index) = clock_start
!           write(888,FMT='("time_start (init) = ",x,i10)') time_start(stat_index)
!           write(888,FMT='(" event = ",a)')descrip_mon(stat_index)
           elapsed_time(stat_index) = 0.0_SP 
           imon(stat_index) = 0
           initialize(stat_index) = .FALSE.
        ENDIF

        IF(stop_timer(stat_index) )THEN
           elsv = elapsed_time(stat_index)
           CALL SYSTEM_CLOCK(COUNT=clock_end) ! end timing
           time_end(stat_index) = clock_end   ! will be 0 if max_ct is exceeded
           IF(clock_end == 0)CALL reset_timer
           elapsed_time(stat_index) =  elapsed_time(stat_index) + &
           REAL(time_end(stat_index)-time_start(stat_index),sp)/clock_rate_sv
!           write(888,FMT='("elt_sv,elt_cur =",3(x,1pe14.6))')elsv, elapsed_time(stat_index),clock_rate_sv
!              write(888,FMT='("tim_e tim_s, =",2(x,i10))')time_end(stat_index),time_start(stat_index)
!              write(888,FMT='(" event = ",a)')descrip_mon(stat_index) ! 888888999
        ENDIF
        
        
        imon(stat_index) = 1  ! keeps track of which indecies contain useful info

   
        RETURN
      END SUBROUTINE collect_stats

      SUBROUTINE reset_timer
      !---------------------------------------------------------------------------------
      ! integer time has reached max integer value, need to reset
      ! will loose some fraction of  timing cycle because of this 
      !---------------------------------------------------------------------------------
        INTEGER(I4B) j
        REAL(SP) tmax_cy
            time_end(:)     = 0 ! reset all because clock count is exceeded
                                ! for all cases simultaneously 
            time_start(:)   = 0
            cycle = cycle+1
            ! account for o12_index which is overall run time
            ! NOTE: this assumes that the o12_index case is called only at
            ! begining and end of run. if run time is less than tmax_cy then
            ! cycle =0 and this subrotuine is not active. 
            ! otherwise cycle .GE. 1 but elapsed time counts only the elapsed time in the 
            ! current cycle for this case. For other cases elapsed time is
            ! accumulated and this correction is not necessary.
            tmax_cy = REAL(max_ct)/REAL(clock_rate)
            elapsed_time(o12_index) = elapsed_time(o12_index) + tmax_cy
!            write(888,FMT='("cycle ,tmax_cy,o12_index = ",i5,x,1pe12.5,x,i5)')cycle,tmax_cy,o12_index
!            DO j=1,nstats 
!               write(888,FMT='(a,x,1pe12.5)')descrip_mon(j),elapsed_time(j) ! 888889999
!            ENDDO
        RETURN 
      END SUBROUTINE reset_timer

      SUBROUTINE print_stats
      !-------------------------------------------------------------------------------
      ! Print out collected statistics
      !-------------------------------------------------------------------------------
        USE io_gcnmp,                               ONLY : nlog,ncrt


        INTEGER(I4B) j
        REAL(SP) model_time,run_time,other_time
        model_time =0.0_SP
        DO j=1,nstats
           IF(imon(j) .GT. 0)THEN
              WRITE(ncrt,FMT='(/,a,x,1pe12.3)')descrip_mon(j),elapsed_time(j)
!              WRITE(nlog,FMT='(/,a,x,1pe12.3)')descrip_mon(j),elapsed_time(j)
              IF( descrip_mon(j)(1:6) .NE. 'Onetwo')THEN
                 model_time = model_time + elapsed_time(j)
              ELSE
                   run_time = elapsed_time(j) 
              ENDIF
           ENDIF
        ENDDO

        other_time =  run_time - model_time
        WRITE(ncrt,FMT='(/,"Remaining compute time = ",1pe12.2)')other_time
!        WRITE(nlog,FMT='(/,"Remaining compute time = ",1pe12.2)')other_time

        RETURN
      END SUBROUTINE print_stats

  END MODULE statistics
