 MODULE events

   USE nrtype,                             ONLY : DP,I4B

   USE common_constants,                   ONLY : zeroc,izero

   IMPLICIT NONE
   INTEGER(I4B),PARAMETER               :: nevents = 10
   INTEGER(I4b) tot_no_events            ! no of valid entries in event_name
   CHARACTER(len=24),DIMENSION(nevents) :: event_name,cur_events
   ! cur_event is list  of events particular to this run
   ! event_name is list of all possible events, independent of 
   ! any run.


   CONTAINS
    
     SUBROUTINE event_list
     !-------------------------------------------------------------------------
     ! -- list of all known events of interest in transport runs.
     ! -- event_name must  have all events up front,the event NONE
     ! -- may not be interspersed with legitimate events.
     ! -- once NONE is encountered the event scan stops
     ! -- changes in this routine must be coordianted with changes in
     ! -- subrotuinecheck_events
     !-----------------------------------------------------------------HSJ------
       tot_no_events       = izero ! counts # defined events
       event_name(1:nevents) = 'NONE'

       event_name(1)       = 'neutral_beam_injection'
       tot_no_events       = tot_no_events +1

       event_name(2)       = 'rf_heating_cd'
       tot_no_events       = tot_no_events +1

       event_name(3)       = 'mhd_equilibrium'
       tot_no_events       = tot_no_events +1

       event_name(4)       = 'sawtooth_crash'
       tot_no_events       = tot_no_events +1

       RETURN

     END SUBROUTINE event_list

     SUBROUTINE set_cur_events
     !-------------------------------------------------------------------------------
     ! -- setup array cur_events, which is lsit of events for a particualr run
     ! -- (event_name is list of all events that this code knows about)
     ! -- give up to nevent event names
     ! -- cur_events  must have all events up front,the event NONE
     ! -- may not be interspersed with legitimate events.
     ! -- once NONE is encountered the event scan stops
     !-----------------------------------------------------------------HSJ-----------
       INTEGER(I4B) k

       DO k =1,nevents
          cur_events(k)  ='NONE'   ! flesh out as needed
       ENDDO
       RETURN
     END SUBROUTINE set_cur_events



     SUBROUTINE check_events(dt)
     !-------------------------------------------------------------------------------
     ! -- all events that can modify the time step dt are governed by this subroutine
     ! -- Note that there is an implicit association between events and 
     ! -- subroutines  **_event
     ! --
     ! -- INPUT/OUTPUT
     !       dt     current time step which may be modified by this  routine
     ! -- INPUT
     !--------------------------------------------------------------------HSJ--------



        REAL(DP) dt,dt_event,dt_eventl
        INTEGER(I4B) k

        dt_event = HUGE(1._DP)

        DO k = 1,nevents
           IF(cur_events(k)  == 'NONE')EXIT  ! exit  the k loop on first encounter of this
           IF(cur_events(k) == event_name(1))THEN
                  CALL neutral_beam_event(dt_eventl)
                  dt_event = MIN(dt_eventl,dt_event)

           ELSE IF(cur_events(k) == event_name(2))THEN
                  CALL rf_heating_cd_event(dt_eventl)
                  dt_event = MIN(dt_eventl,dt_event)

           ELSE IF(cur_events(k) == event_name(3))THEN
                  CALL mhd_equilibrium_event(dt_eventl)
                  dt_event = MIN(dt_eventl,dt_event)

           ELSE IF(cur_events(k) == event_name(4))THEN
                  CALL sawtooth_crash_event(dt_eventl)
                  dt_event = MIN(dt_eventl,dt_event)

           ELSE
                  CALL set_event_error(k)
           ENDIF

        ENDDO

        ! dt_event is now the smallest dt for all active events
        ! reset the time step dt to this value:
        dt = MIN(dt,dt_event) ! ? set time+dt  here??

        RETURN

        END SUBROUTINE check_events



        SUBROUTINE neutral_beam_event(dt_eventl)
        !-------------------------------------------------------------------------------
        ! -- 
        !--------------------------------------------------------------------HSJ--------

            REAL(DP) dt_eventl
               dt_eventl = HUGE(1._DP)
            RETURN
        END SUBROUTINE neutral_beam_event

        SUBROUTINE rf_heating_cd_event(dt_eventl)
        !-------------------------------------------------------------------------------
        ! -- 
        !--------------------------------------------------------------------HSJ--------
            REAL(DP) dt_eventl
               dt_eventl = HUGE(1._DP)
            RETURN
        END SUBROUTINE rf_heating_cd_event
  
        SUBROUTINE mhd_equilibrium_event(dt_eventl)
        !-------------------------------------------------------------------------------
        ! -- 
        !--------------------------------------------------------------------HSJ--------
            REAL(DP) dt_eventl
               dt_eventl = HUGE(1._DP)
            RETURN
        END SUBROUTINE mhd_equilibrium_event
 
        SUBROUTINE sawtooth_crash_event(dt_eventl)
        !-------------------------------------------------------------------------------
        ! -- 
        !--------------------------------------------------------------------HSJ--------
            REAL(DP) dt_eventl
               dt_eventl = HUGE(1._DP)
            RETURN
        END SUBROUTINE sawtooth_crash_event

        SUBROUTINE set_event_error(index)
        !-------------------------------------------------------------------------------
        ! -- call error handler and exit code
        ! -- INPUT
        !       index ==>   cur_event(index) caused error
        !--------------------------------------------------------------------HSJ--------

         USE error_handler,              ONLY : iomaxerr, lerrno,terminate

         USE io_gcnmp,                   ONLY : nlog,ncrt

         USE MPI_data,                   ONLY : myid,master

         INTEGER(I4B) index,k

         IF(myid == master)THEN
            WRITE(ncrt,FMT='(" error, event not recognized:",a)')cur_events(index)
            WRITE(ncrt,FMT='(" The list of known events is")')
            DO k=1,tot_no_events
               WRITE(ncrt,FMT='(a)')event_name(k)
            ENDDO
         ENDIF

         lerrno = iomaxerr+183
         CALL terminate(lerrno,nlog)
         RETURN

     END SUBROUTINE set_event_error

 END MODULE events
