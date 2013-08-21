
     SUBROUTINE set_eqdsk_statefile_values
!----------------------------------------------------------------
! -- create output for state file from eqdsk info
! -- can be supplied by reading single eqdsk or from tdem mode
! -- put eqdsk quantites on
! -- mhd_dat%psivalnpsi%data(1:mhd_dat%npsi) edge to axis grid
! --
!---------------------------------------------------------HSJ----

     USE nrtype,                           ONLY : DP,I4B

     USE gpsi,                             ONLY : pppsi_eqdsk,presspsi_eqdsk, &
                                                  fpsi_eqdsk,ffppsi_eqdsk,    &
                                                  qpsi_eqdsk,nxeqd_eqdsk,     &
                                                  psimag_eqdsk,psilim_eqdsk,  &
                                                  psival_eqdsk

     USE plasma_properties,                ONLY : mhd_dat,dischg

     USE cubic_spline,                     ONLY : spl1d_inst,eval_1d_cubic_spline, &
                                                  create_1d_spline_object,         &
                                                  done_1dspline

     USE vector_class,                     ONLY : zero_vector ,delete_vector_nf

     USE grid_class,                       ONLY : reverse

     USE mhdcom,                           ONLY :  psiaxis,psibdry,rma,zma
    
     USE constnts,                         ONLY :  psikgaus,psimks 
 
     USE flxav,                            ONLY :  npsi

     USE psig,                             ONLY : presspsi,fpsi,ffppsi,qpsi,pppsi,psival

     IMPLICIT NONE

     REAL(DP),DIMENSION(:),ALLOCATABLE :: work_idep,work_prof,work_idep2

     REAL(DP) dpsi,ax,bx
  
     INTEGER(I4B) j,ib,ok,npsi_lcl

     TYPE(spl1d_inst) gen_spline 
             
          CALL done_1dspline(gen_spline)  ! clear any existing instance of gen_spline

          mhd_dat%npsi  = npsi  ;   nxeqd_eqdsk       = mhd_dat%npsi

          ALLOCATE(work_idep(nxeqd_eqdsk),work_prof(nxeqd_eqdsk),work_idep2(mhd_dat%npsi))
          IF(ALLOCATED(psival_eqdsk))DEALLOCATE(psival_eqdsk)
          ALLOCATE(psival_eqdsk(nxeqd_eqdsk))
          psival_eqdsk(:) = psival(1:npsi)*psimks

          ok = delete_vector_nf(mhd_dat%psivalnpsi)
          mhd_dat%psivalnpsi = zero_Vector(mhd_dat%npsi)
          mhd_dat%psiaxis    = psiaxis*psimks
          mhd_dat%psibdry    = psibdry*psimks

          mhd_dat%psivalnpsi%data(1:npsi) = psival(1:npsi)*psimks
          ! state file convention is that mhd_dat%psivalnpsi%data(1) = edge,
          ! mhd_dat%psivalnpsi%data(npsi) = axis

 write(777,FMT='("output to statefile")')
 write(777,FMT='("npsi,mhd_dat%npsi,nxeqd_eqdsk =",3(2x,i5))')npsi,mhd_dat%npsi,nxeqd_eqdsk
 write(777,FMT='("psiaxis,psibdry =",2(1pe14.8,2x))')psiaxis,psibdry
 write(777,FMT='("mhd_dat%psivalnpsi =")')
 write(777,FMT='(5(2x,1pe14.6))')(mhd_dat%psivalnpsi%data(j),j = 1,mhd_dat%psivalnpsi%size)
 write(777,FMT='("psival =")')
 write(777,FMT='(5(2x,1pe14.6))')(psival(j),j = 1,npsi)

          work_idep2(:)           = mhd_dat%psivalnpsi%data(:)
          CALL reverse(mhd_dat%npsi,work_idep2)
          work_idep(1)            =  mhd_dat%psivalnpsi%data(1)
          work_idep(nxeqd_eqdsk)  =  mhd_dat%psivalnpsi%data(mhd_dat%npsi)

          dpsi = (work_idep(nxeqd_eqdsk) - work_idep(1))/(nxeqd_eqdsk -1)
          DO j=2,nxeqd_eqdsk -1
             work_idep(j) = work_idep(j-1) + dpsi
          ENDDO
          CALL reverse(nxeqd_eqdsk,work_idep)

! ---   press
          IF(ALLOCATED(presspsi_eqdsk))DEALLOCATE(presspsi_eqdsk)
          ALLOCATE(presspsi_eqdsk(nxeqd_eqdsk))
          presspsi_eqdsk(1:npsi) = presspsi(1:npsi)/10._DP  ! erg/cm^3 ==> nt/m^3


          !work_prof(1:nxeqd_eqdsk) = presspsi_eqdsk(1:nxeqd_eqdsk)
          !CALL reverse(nxeqd_eqdsk,work_prof)
 
          !ax = (work_prof(2)-work_prof(1))/(work_idep(2)-work_idep(1)) 
          !bx = (work_prof(nxeqd_eqdsk)-work_prof(nxeqd_eqdsk-1))/ &
           !                  (work_idep(nxeqd_eqdsk)-work_idep(nxeqd_eqdsk -1))
          !CALL create_1d_spline_object(gen_spline,nxeqd_eqdsk,work_idep, &
          !work_prof,ax,bx,ib)
          ib =1 
          ax = (presspsi_eqdsk(2)-  presspsi_eqdsk(1))/(psival(2)-psival(1))
          bx = (presspsi_eqdsk(npsi)-  presspsi_eqdsk(npsi-1))/(psival(npsi)-psival(npsi-1))

          CALL create_1d_spline_object(gen_spline,npsi,psival_eqdsk, &
                                        presspsi_eqdsk,ax,bx,ib)

          ok = delete_Vector_nf (mhd_dat%pressnpsi)
          mhd_dat%pressnpsi      = zero_Vector(mhd_dat%npsi)
          DO j=2, mhd_dat%npsi-1
               gen_spline%xx = mhd_dat%psivalnpsi%data(j)
               CALL eval_1d_cubic_spline(gen_spline)
               mhd_dat%pressnpsi%data(j) = gen_spline%sp
          ENDDO
           mhd_dat%pressnpsi%data(1)             = presspsi_eqdsk(nxeqd_eqdsk)
           mhd_dat%pressnpsi%data(mhd_dat%npsi)  = presspsi_eqdsk(1)
          !CALL reverse(mhd_dat%npsi,mhd_dat%pressnpsi%data)

 write(777,FMT='("mhd_dat%pressnpsi =",5(2x,1pe14.6))')(mhd_dat%pressnpsi%data(j),j = 1,mhd_dat%pressnpsi%size)
 write(777,FMT='("presspsi =",5(2x,1pe14.6))')(presspsi(j),j = 1,npsi)
        call stop('in set_eqdsk statfile',1)


! ---   pppsi
          IF(.NOT. ALLOCATED(pppsi_eqdsk))THEN
             ALLOCATE(pppsi_eqdsk(mhd_dat%npsi))
           !  pppsi_eqdsk(:) = ??
          ELSE
            ! maxj = SIZE(pppsi_eqdsk)
          ENDIF
          work_prof(1:nxeqd_eqdsk) = pppsi_eqdsk(1:nxeqd_eqdsk)
          CALL reverse(nxeqd_eqdsk,work_prof)
          ib =1 ; ax = (work_prof(2)-work_prof(1))/(work_idep(2)-work_idep(1)) 
          bx = (work_prof(nxeqd_eqdsk)-work_prof(nxeqd_eqdsk-1))/ &
                             (work_idep(nxeqd_eqdsk)-work_idep(nxeqd_eqdsk -1))
          CALL create_1d_spline_object(gen_spline,nxeqd_eqdsk,work_idep, &
                                       work_prof,ax,bx,ib)
          ok = delete_Vector_nf (mhd_dat%pprimnpsi)
          mhd_dat%pprimnpsi      = zero_Vector(mhd_dat%npsi)
          DO j=2, mhd_dat%npsi-1
               gen_spline%xx = work_idep2(j)
               CALL eval_1d_cubic_spline(gen_spline)
               mhd_dat%pprimnpsi%data(j) = gen_spline%sp
          ENDDO
          mhd_dat%pprimnpsi%data(1)             = work_prof(1)
          mhd_dat%pprimnpsi%data(mhd_dat%npsi) = work_prof(nxeqd_eqdsk)
          CALL reverse(mhd_dat%npsi,mhd_dat%pprimnpsi%data)



! ---   fpsi
          work_prof(1:nxeqd_eqdsk) = fpsi_eqdsk(1:nxeqd_eqdsk)
          CALL reverse(nxeqd_eqdsk,work_prof)
          ib =1 ; ax = (work_prof(2)-work_prof(1))/(work_idep(2)-work_idep(1)) 
          bx = (work_prof(nxeqd_eqdsk)-work_prof(nxeqd_eqdsk-1))/ &
                             (work_idep(nxeqd_eqdsk)-work_idep(nxeqd_eqdsk -1))
          CALL create_1d_spline_object(gen_spline,nxeqd_eqdsk,work_idep, &
                                       work_prof,ax,bx,ib)
          ok = delete_Vector_nf (mhd_dat%fpsinpsi)
          mhd_dat%fpsinpsi      = zero_Vector(mhd_dat%npsi)
          DO j=2, mhd_dat%npsi-1
               gen_spline%xx = work_idep2(j)
               CALL eval_1d_cubic_spline(gen_spline)
               mhd_dat%fpsinpsi%data(j) = gen_spline%sp
          ENDDO
          mhd_dat%fpsinpsi%data(1)             = work_prof(1)
          mhd_dat%fpsinpsi%data(mhd_dat%npsi)  = work_prof(nxeqd_eqdsk)
          CALL reverse(mhd_dat%npsi,mhd_dat%fpsinpsi%data)

! ---   ffppsi
          work_prof(1:nxeqd_eqdsk) = ffppsi_eqdsk(1:nxeqd_eqdsk)
          CALL reverse(nxeqd_eqdsk,work_prof)
          ib =1 ; ax = (work_prof(2)-work_prof(1))/(work_idep(2)-work_idep(1)) 
          bx = (work_prof(nxeqd_eqdsk)-work_prof(nxeqd_eqdsk-1))/ &
                             (work_idep(nxeqd_eqdsk)-work_idep(nxeqd_eqdsk -1))
          CALL create_1d_spline_object(gen_spline,nxeqd_eqdsk,work_idep, &
                                       work_prof,ax,bx,ib)
          ok = delete_Vector_nf (mhd_dat%ffprimnpsi)
          mhd_dat%ffprimnpsi      = zero_Vector(mhd_dat%npsi)
          DO j=2, mhd_dat%npsi-1
               gen_spline%xx = work_idep2(j)
               CALL eval_1d_cubic_spline(gen_spline)
               mhd_dat%ffprimnpsi%data(j) = gen_spline%sp
          ENDDO
          mhd_dat%ffprimnpsi%data(1)             = work_prof(1)
          mhd_dat%ffprimnpsi%data(mhd_dat%npsi) = work_prof(nxeqd_eqdsk)
          CALL reverse(mhd_dat%npsi,mhd_dat%ffprimnpsi%data)


! ---   qpsi
          work_prof(1:nxeqd_eqdsk) = qpsi_eqdsk(1:nxeqd_eqdsk)
          CALL reverse(nxeqd_eqdsk,work_prof)
          ib =1 ; ax = (work_prof(2)-work_prof(1))/(work_idep(2)-work_idep(1)) 
          bx = (work_prof(nxeqd_eqdsk)-work_prof(nxeqd_eqdsk-1))/ &
                             (work_idep(nxeqd_eqdsk)-work_idep(nxeqd_eqdsk -1))
          CALL create_1d_spline_object(gen_spline,nxeqd_eqdsk,work_idep, &
                                       work_prof,ax,bx,ib)
          ok = delete_Vector_nf (mhd_dat%qpsinpsi)
          mhd_dat%qpsinpsi      = zero_Vector(mhd_dat%npsi)
          DO j=2, mhd_dat%npsi-1
               gen_spline%xx = work_idep2(j)
               CALL eval_1d_cubic_spline(gen_spline)
               mhd_dat%qpsinpsi%data(j) = gen_spline%sp
          ENDDO
          mhd_dat%qpsinpsi%data(1)             = work_prof(1)
          mhd_dat%qpsinpsi%data(mhd_dat%npsi) = work_prof(nxeqd_eqdsk)
          CALL reverse(mhd_dat%npsi,mhd_dat%qpsinpsi%data)


          DEALLOCATE(work_idep,work_prof,work_idep2)
          CALL done_1dspline(gen_spline)  

      RETURN

     END SUBROUTINE set_eqdsk_statefile_values
