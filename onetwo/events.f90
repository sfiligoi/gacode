
  MODULE events
      USE param,                                       ONLY : kevents
      USE mhdpar,                                      ONLY : mxtbcmhd
	integer                                                        &
                index_bc,index_pbon,index_pboff,index_beam_ul,         &
                index_beam_ll, nevents, ievent(kevents),               &
                index_p_nf_start,index_p_nf_stop,isave_event

        real *8                                                        &
                timevent(kevents+mxtbcmhd),toray_run_times(kevents+mxtbcmhd),            &
                statefile_write_times(kevents+mxtbcmhd),tdem_run_times(kevents+mxtbcmhd)

        CHARACTER*64,event_array(kevents)

      data index_p_nf_start,index_p_nf_stop /78,141/

      data event_array /"timmax","3D plot point (timplt)",               &
             "print point (prtlst,timprt)","beam plot point (timbplt)",  &
             "mhd equilibrium time step (dteq)","plot point (pltlst)",   &
             "rf plot point (timrfp)","pellet injection time",           &
             "beam on","beam off","rf(1) on","rf(1) off",                &
             "rf(2) on" ,"rf(2) off","rf(3) on" ,"rf(3) off",            &
             "rf(4) on" ,"rf(4) off","rf(5) on" ,"rf(5) off",            &
             "rf(6) on" ,"rf(6) off","rf(7) on" ,"rf(7) off",            &
            "rf(8) on" ,"rf(8) off","rf(9) on" ,"rf(9) off",             &
             "rf(10) on" ,"rf(10) off","rf(11) on" ,"rf(11) off",        &
             "rf(12) on" ,"rf(12) off","rf(13) on" ,"rf(13) off",        &
             "rf(14) on" ,"rf(14) off","rf(15) on" ,"rf(15) off",        &
             "rf(16) on" ,"rf(16) off","rf(17) on" ,"rf(17) off",        &
             "rf(18) on" ,"rf(18) off","rf(19) on" ,"rf(19) off",        &
             "rf(20) on" ,"rf(20) off","rf(21) on" ,"rf(21) off",        &
             "rf(22) on" ,"rf(22) off","rf(23) on" ,"rf(23) off",        &
             "rf(24) on" ,"rf(24) off","rf(25) on" ,"rf(25) off",        &
             "rf(26) on" ,"rf(26) off","rf(27) on" ,"rf(27) off",        &
             "rf(28) on" ,"rf(28) off","rf(29) on" ,"rf(29) off",        &
             "rf(30) on" ,"rf(30) off",                                  &
             "time dependent boundary conditions",                       &
             "pulsed beam on","pulsed beam off",                         &
             "beam integral upper limit", "beam integral lower limit",   &
             "nubeam switching  time", " " ,                             &
             "P_Nf_beam(1) on","P_Nf_beam(1) off",                       &
             "P_Nf_beam(2) on","P_Nf_beam(2) off",                       &
             "P_Nf_beam(3) on","P_Nf_beam(3) off",                       &
             "P_Nf_beam(4) on","P_Nf_beam(4) off",                       & 
             "P_Nf_beam(5) on","P_Nf_beam(5) off",                       &
             "P_Nf_beam(6) on","P_Nf_beam(6) off",                       &
             "P_Nf_beam(7) on","P_Nf_beam(7) off",                       &
             "P_Nf_beam(8) on","P_Nf_beam(8) off",                       & 
             "P_Nf_beam(9) on","P_Nf_beam(9) off",                       &
             "P_Nf_beam(10) on","P_Nf_beam(10) off",                     &
             "P_Nf_beam(11) on","P_Nf_beam(11) off",                     &
             "P_Nf_beam(12) on","P_Nf_beam(12) off",                     & 
             "P_Nf_beam(13) on","P_Nf_beam(13) off",                     &
             "P_Nf_beam(14) on","P_Nf_beam(14) off",                     &
             "P_Nf_beam(15) on","P_Nf_beam(15) off",                     &
             "P_Nf_beam(16) on","P_Nf_beam(16) off",                     & 
             "P_Nf_beam(17) on","P_Nf_beam(17) off",                     &
             "P_Nf_beam(18) on","P_Nf_beam(18) off",                     &
             "P_Nf_beam(19) on","P_Nf_beam(19) off",                     & 
             "P_Nf_beam(20) on","P_Nf_beam(20) off",                     &
             "P_Nf_beam(21) on","P_Nf_beam(21) off",                     &
             "P_Nf_beam(22) on","P_Nf_beam(22) off",                     &
             "P_Nf_beam(23) on","P_Nf_beam(23) off",                     &
             "P_Nf_beam(24) on","P_Nf_beam(24) off",                     &
             "P_Nf_beam(25) on","P_Nf_beam(25) off",                     &
             "P_Nf_beam(26) on","P_Nf_beam(26) off",                     &
             "P_Nf_beam(27) on","P_Nf_beam(27) off",                     &
             "P_Nf_beam(28) on","P_Nf_beam(28) off",                     &
             "P_Nf_beam(29) on","P_Nf_beam(29) off",                     &
             "P_Nf_beam(30) on","P_Nf_beam(30) off",                     &
             "P_Nf_beam(31) on","P_Nf_beam(31) off",                     &
             "P_Nf_beam(32) on","P_Nf_beam(32) off"    /

!
! index_bc,index_pbon,index_pboff,index_beam_ul,index_beam_ll
! are indecies into the timevent
! array see cray102.f




  CONTAINS

      SUBROUTINE assign_events
!---------------------------------------------------------------------------
!-- Given prtlst,rtime_tdem,rfon,rftime
!-- create a new master print list that reflects all these times
!-- for printout.
!-- create a new prtlist time that contains all of the times
!---------------------------------------------------------------------------
      USE nrtype,                                      ONLY : DP,I4B
      USE param,                                       ONLY : krf,kevents
      USE rf,                                          ONLY : rfon
      USE mhdcom,                                      ONLY : mhdmethd
      USE tdem,                                        ONLY : rtime_tdem,ntime_tdem
      USE io,                                          ONLY : nout,ncrt,nitre,prtlst, &
                                                              prtlst_used
      USE solcon,                                      ONLY : timmax

      IMPLICIT NONE

      INTEGER j,nprt,ntotl,nrf

      REAL(DP),ALLOCATABLE,DIMENSION(:) :: temp_list
      REAL(DP),ALLOCATABLE,DIMENSION(:) :: rfon_cpy

      INTERFACE
       SUBROUTINE  sort_unique(a,km)
         REAL *8, intent(INOUT),dimension(:) :: a
         INTEGER, INTENT(OUT) :: km
       END SUBROUTINE  sort_unique
      END INTERFACE


      IF(SIZE(tdem_run_times) .lt. ntime_tdem)THEN
         WRITE(nout,10)kevents,ntime_tdem
         WRITE(nitre,10)kevents,ntime_tdem
         WRITE(ncrt,10)kevents,ntime_tdem
10       FORMAT(2x,'Error ,kevents  arrays size is too small',/, &
                2x,'Need at least ',2(x,i5))
         call STOP('recompile code, kevents  size problem',1)
      ENDIF
      ntotl = 0
      IF( mhdmethd == 'tdem')THEN
         IF( ntime_tdem > 0 )THEN
            tdem_run_times(1:ntime_tdem)          = rtime_tdem(1:ntime_tdem)
            toray_run_times(1:ntime_tdem)         = rtime_tdem(1:ntime_tdem)
            statefile_write_times(1:ntime_tdem)   = rtime_tdem(1:ntime_tdem) 
         ELSE
            WRITE(ncrt,1)
            WRITE(nout,1)
            WRITE(nitre,1) ! runlog file
1           Format(2x,'ERROR in  sub assign_events: array size problem')
            CALL STOP('sub assign_events',1)
         ENDIF
      ENDIF

      CALL sort_unique(prtlst,nprt)
      prtlst_used(:,:) = 0
!      CALL sort_unique(prtlst,nrf)
      ALLOCATE(rfon_cpy(krf)) 
      CALL sort_unique(rfon_cpy,nrf) ! get nrf
      DEALLOCATE(rfon_cpy)
      ntotl = ntotl + nprt + nrf + ntime_tdem 

      ALLOCATE(temp_list(ntotl))

      temp_list(1:nprt) = prtlst(1:nprt)
      temp_list(nprt+1:nprt+ntime_tdem)= rtime_tdem(1:ntime_tdem)
      temp_list(nprt+ntime_tdem+1:nprt+ntime_tdem +nrf) = rfon(1:nrf)
 
!      print *,'Size temp_list =',SIZE(temp_list)
      ! sort temp_list in ascending order:
      CALL sort_unique(temp_list,nprt)
!      print *,'new_size =',nprt
      IF(SIZE(prtlst) .GE. nprt)THEN
         prtlst(1:nprt) = temp_list(1:nprt)
      ELSE
         Write(ncrt,12)SIZE(prtlst),nprt
  12     FORMAT(2x,'this case can not be run without increasing',/, &
                2x,' the size of the prtlst array:',/,              &
                2x,' current size ',i5,/,                           &
                2x,'required size after merger with rf,tdem,times:',&
                i5)
         WRITE(nitre,10)SIZE(prtlst),nprt
         call STOP('recompile code, prtlist size problem',1)
      ENDIF
      DEALLOCATE(temp_list)
      RETURN
      END SUBROUTINE assign_events



      SUBROUTINE check_prtlst(flag,model)
!---------------------------------------------------------------------------
!-- Given prtlst and current time
!-- return 1 if ABS(time-prtlist(i)) < tol)
!-- otherwise return 0
!--------------------------------------------------------------------HSJ----
      USE nrtype,                                      ONLY : DP,I4B
      USE io,                                          ONLY : prtlst,prtlst_used
      USE solcon,                                      ONLY : time,timmax ,time_tol
      USE file_proc,                                   ONLY : temp_file_test

      IMPLICIT NONE

      INTEGER j,flag,js,model
      REAL(DP) diff,diffmin
      flag = 0
      DO j=1,SIZE(prtlst)
         IF( ABS(time - prtlst(j)) < 5._DP*time_tol )THEN
            flag = 1
            EXIT
         ENDIF
      ENDDO

      diff = 1.e30 ; flag =0
      do j=1,SIZE(prtlst)
          if(prtlst_used(j,model) == 0)THEN
             diffmin =  ABS(time - prtlst(j))
             diff = MIN(diff,diffmin)
             if(diff == diffmin)js = j
          ENDIF
      enddo

      IF(diff .lt. 1.e-6_DP)THEN
           prtlst_used(js,model) = 1
           flag = 1
      ENDIF
!        write(940,FMT='("check_prt,time,tol,diff,lst=",4(x,1pe12.4))') &
!                time,time_tol,diff,prtlst(js)
!       write(940,FMT='("in check_prtlst,fla,abs",i3,x,1pe12.4)') &
!                flag,ABS(time - prtlst(js))
!      write(940,FMT='("in check_prtlst,diff=",1pe12.4,x,1pe12.4)') &
!               time,time_tol
    END SUBROUTINE check_prtlst


  END MODULE events
