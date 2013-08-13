
      SUBROUTINE line_avg_den
!-------------------------------------------------------------------------------
! -- get electron line avg density by integrating
! -- from inboard edge to  axis and axis  to outboard edge 
! -- (both at elevation of magnetic axis)
! -- routine is also used to put out densities and tmeperatures on rmajor grid
! --  OUTPUT
!         profile%ene_bar_inbd 
!         profile%ene_bar_outbd 
!         R_chd Major radius chord array at elevation of zelev (= mag axis z)
!               Extends from inboard to outboard side
!         psi_chord psi     on chord
!         rho_chord rho     on chord
!         ene_chord density on chord
!         te_chord  elect temp on chord
!         ti_chord  ion temp on chord
!         
!         IF an error is detected (because there are inconsistencies
!         in the eqdsk psi values) hen ierr = 1 is used to set the
!         above output to all zeros. The code does not terminate 
!         with this error however.
!--------------------------------------------------------------HSJ------

        USE nrtype,                      ONLY : Dp,I4B

        USE bicubic_spline,              ONLY : set_2dspline_dim,done_2dspline, &
                                                sets2d_spline,eval_2d_spline

        USE common_constants,            ONLY : zeroc,izero

        USE cubic_spline,                ONLY : eval_1d_cubic_spline,      &
                                                create_1d_spline_object,   &
                                                done_1dspline, spl1d_inst
        USE Plasma_properties,           ONLY : dischg,mhd_dat,profile
 
        USE grid_class,                  ONLY : rho_grid,psir_grid,nj, &
                                                R_chord,psi_chord,     &
                                                ene_chord,rho_chord,   &
                                                te_chord,ti_chord

        USE io_gcnmp,                    ONLY : ncrt,nlog

        USE Vector_class,                ONLY : assign_Vector,zero_Vector

        USE MPI_data,                    ONLY : myid,master 
  
        USE solcon_gcnmp,                ONLY : time

        IMPLICIT NONE

        INTEGER(I4B) ib,ier,icalc,i,j,k,kax,ierr
        REAL(DP) psi1d(dischg%nr_mhd*dischg%nz_mhd) ! temp array
        REAL(DP) fs(6),Rl,Rr,dr,drhodpsi,drho2dpsi,ax,bx,zelev,xx,  &
                 psisgn,rax,raxmin,ene_inbd,ene_outbd,psibdry_lcl
        REAL(DP) r_chd(dischg%nr_mhd),psi_rchd(dischg%nr_mhd),       &
                 rho_rchd(dischg%nr_mhd),ene_rchd(dischg%nr_mhd),   &
                 workx(nj),worky(nj),te_rchd(dischg%nr_mhd),        &
                 ti_rchd(dischg%nr_mhd)
        TYPE(spl1d_inst)spline_objct

         profile%ene_bar_inbd = zeroc  ; profile%ene_bar_outbd = zeroc
        !----------------------------------------------------------------
        ! 2d spline processing of psi(R,Z), no error return through ier ?
        !----------------------------------------------------------------

        ! set up required arrays:
          CALL set_2dspline_dim(dischg%nr_mhd,dischg%nz_mhd)

        ! get spline coefficients, store  in back end
          DO i=1,dischg%nr_mhd
             DO j=1,dischg%nz_mhd
                k=(i-1)*dischg%nz_mhd+j
                psi1d(k) = mhd_dat%psi(i,j)
             ENDDO
          ENDDO
 
          CALL sets2d_spline(psi1d,dischg%rmhdgrid%data,dischg%nr_mhd, &
                             dischg%zmhdgrid%data,dischg%nz_mhd)

          ! 2d function is now ready for evaluation

          ! find inboard,Rl, and outboard,Rr, rmajor of plasma boundary
          ! at elevation zelev"
          zelev = dischg%zma
          psisgn =  mhd_dat%psiaxis-mhd_dat%psibdry
          psibdry_lcl = mhd_dat%psibdry
          CALL find_Rb(Rl,Rr,zelev,psisgn,psibdry_lcl,ierr)
 
          IF(ierr == 0)THEN
             
             !define a grid for the chord from Rl to Rr at zelev
             !  and get psi along this chord:
             dr = (Rr -Rl)/(dischg%nr_mhd-1)
             r_chd(1) = Rl
             icalc =1
             raxmin = HUGE(1._DP)
             DO j=1,dischg%nr_mhd
                r_chd(j) =  Rl  + (j-1)*dr
                rax = abs(r_chd(j) - dischg%rma)
                raxmin = MIN(rax,raxmin)
                IF(rax .EQ. raxmin) kax = j
                CALL eval_2d_spline(r_chd(j),zelev,fs,icalc)
                psi_rchd(j) = fs(1)
             ENDDO
             ! reset kax value to magnetic axis exactly:
             r_chd(kax) =  dischg%rma
             r_chd(dischg%nr_mhd) = Rr
             !CALL eval_2d_spline(r_chd(kax),zelev,fs,icalc)
             psi_rchd(kax) = mhd_dat%psiaxis
             !--------------------------------------------------------------
             ! get rho on this chord. Recall that there is a 1 to 1 map
             ! between psir_grid,rho_grid,and all profile%** defined 
             ! on the rho grid
             !---------------------------------------------------------------

             CALL done_1dspline(spline_objct) 
             ib = 4      ! Third deriv con at left anf right end poiints
             ib = 1      ! use deriv at rho = 0,1.0
             ax = 1.e5_DP ! not used if ib =3,4
             ax =  rho_grid%data(2)/(psir_grid%data(2)- psir_grid%data(1) )
             bx = (rho_grid%data(nj)- rho_grid%data(nj-1))/(psir_grid%data(nj)- psir_grid%data(nj-1) ) ! not used if ib =3,4
             IF(psisgn .LT. zeroc)THEN ! psi is min on axis
                workx(:)  = psir_grid%data(:)   !  rmhd_dat%psivalnpsi%data
                worky(:)  = rho_grid%data(:)    !rho_mhd_gridnpsi%data

             ELSE     ! psi is max on axis
                ! spline requires monotonic increasing abscissa
                ! reverse psi and rho
                DO j=nj,1,-1
                   workx(nj-j+1)  = psir_grid%data(j)
                   worky(nj-j+1)  = rho_grid%data(j)

                ENDDO

             ENDIF

             CALL create_1d_spline_object(spline_objct,nj,workx,worky,ax,bx,ib)

             rho_rchd(1) = rho_grid%data(nj)
             rho_rchd(dischg%nr_mhd)    = rho_grid%data(nj)
             ene_rchd(1)                = profile%ene%data(nj)
             ene_rchd(dischg%nr_mhd)    = profile%ene%data(nj)   
             DO j=2,dischg%nr_mhd-1
                ! given rho_grid as function of psir_grid get rho_rchd as function
                ! of psi_rchd
                 spline_objct%xx = psi_rchd(j)
                 CALL eval_1d_cubic_spline(spline_objct)
                 rho_rchd(j) = spline_objct%sp
             ENDDO
             rho_rchd(kax) = zeroc ! avoids roundoff problems

             CALL done_1dspline(spline_objct) 

             !--------------------------------------------------------
             ! given rho on r_chd we can now get ene on rchd:
             !--------------------------------------------------------
             ib = 4      ! Third deriv con at left anf right end poiints
             ax = 0.0_DP ! not used if ib =3,4
             bx = 0.0_DP ! not used if ib =3,4
             !spline_objct%x(:) = rho_grid%data(:)
             !spline_objct%y(:) = profile%ene%data(:)

             CALL create_1d_spline_object(spline_objct,nj,rho_grid%data,profile%ene%data,ax,bx,ib)
             ene_rchd(1)               = profile%ene%data(nj)
             ene_rchd(dischg%nr_mhd)   = profile%ene%data(nj)  

             DO j=2 , dischg%nr_mhd-1
                 spline_objct%xx = rho_rchd(j) 
                 CALL eval_1d_cubic_spline(spline_objct)
                 ene_rchd(j) = spline_objct%sp
             ENDDO
             ene_rchd(kax) = profile%ene%data(1)


             !-----------------------------------------------------------
             !-- finally we can do the line average integrals
             !-----------------------------------------------------------

             ene_inbd = zeroc   ; ene_outbd = zeroc
             Do j = 1,kax-1
                ene_inbd = ene_inbd + (r_chd(j+1)-r_chd(j))*(ene_rchd(j)+ene_rchd(j+1))
             ENDDO
             Do j = kax,dischg%nr_mhd-1
                ene_outbd = ene_outbd + (r_chd(j+1)-r_chd(j))*(ene_rchd(j)+ene_rchd(j+1))
             ENDDO
             ene_inbd =0.5_DP*ene_inbd/(r_chd(kax)-r_chd(1)) 
             ene_outbd = 0.5_DP*ene_outbd/(r_chd(dischg%nr_mhd)-r_chd(kax))
             profile%ene_bar_inbd   = ene_inbd
             profile%ene_bar_outbd  = ene_outbd
          ELSE ! ierr =1 calcs not done
            profile%ene_bar_inbd  = zeroc
            profile%ene_bar_outbd = zeroc
            R_chord    = zero_Vector(SIZE(r_chd))
            rho_chord  = zero_Vector(SIZE(rho_rchd))
            ene_chord  = zero_Vector(SIZE(ene_rchd))
            psi_chord  = zero_Vector(SIZE(psi_rchd))
          ENDIF

          IF(myid == master)THEN
              WRITE(nlog,FMT='(" line avg ene inboard  = ",1pe12.3," at time ",1pe14.6)')profile%ene_bar_inbd,time
              WRITE(nlog,FMT='(" line avg ene outboard = ",1pe12.3," at time ", 1pe14.6)')profile%ene_bar_outbd,time
          ENDIF
 

          ! WRITE(999,FMT='("j ,r psi,rho,ene")')
          ! DO j=1,dischg%nr_mhd
          !    WRITE(999,FMT='(i5,4(2x,1pe14.6))')j,r_chd(j),psi_rchd(j),rho_rchd(j),ene_rchd(j)
          ! ENDDO
 
          ! release spline arrays:
          CALL done_2dspline
          CALL done_1dspline(spline_objct)



             !--------------------------------------------------------
             ! given rho on r_chd we can now get Te on r_chd:
             !--------------------------------------------------------
             ib = 4      ! Third deriv con at left anf right end poiints
             ax = 0.0_DP ! not used if ib =3,4
             bx = 0.0_DP ! not used if ib =3,4


             CALL create_1d_spline_object(spline_objct,nj,rho_grid%data,profile%te%data,ax,bx,ib)
             te_rchd(1)               = profile%te%data(nj)
             te_rchd(dischg%nr_mhd)   = profile%te%data(nj)  

             DO j=2 , dischg%nr_mhd-1
                 spline_objct%xx = rho_rchd(j) 
                 CALL eval_1d_cubic_spline(spline_objct)
                 te_rchd(j) = spline_objct%sp
             ENDDO
             te_rchd(kax) = profile%te%data(1)
             CALL done_1dspline(spline_objct)



             !--------------------------------------------------------
             ! given rho on r_chd we can now get Ti on r_chd:
             !--------------------------------------------------------
             ib = 4      ! Third deriv con at left anf right end poiints
             ax = 0.0_DP ! not used if ib =3,4
             bx = 0.0_DP ! not used if ib =3,4


             CALL create_1d_spline_object(spline_objct,nj,rho_grid%data,profile%ti%data,ax,bx,ib)
             ti_rchd(1)               = profile%ti%data(nj)
             ti_rchd(dischg%nr_mhd)   = profile%ti%data(nj)  

             DO j=2 , dischg%nr_mhd-1
                 spline_objct%xx = rho_rchd(j) 
                 CALL eval_1d_cubic_spline(spline_objct)
                 ti_rchd(j) = spline_objct%sp
             ENDDO
             ti_rchd(kax) = profile%ti%data(1)

             CALL done_1dspline(spline_objct)


             !-----------------------------------------------------------
             ! -- Permanent storage (for output to summary and statefile)
             !-----------------------------------------------------------
             profile%ene_bar_inbd   = 0.5_DP*ene_inbd /(dischg%rma - Rl)
             profile%ene_bar_outbd  = 0.5_DP*ene_outbd/(Rr - dischg%rma)
             R_chord      = assign_Vector(r_chd)
             rho_chord    = assign_Vector(rho_rchd)
             ene_chord    = assign_Vector(ene_rchd)
             psi_chord    = assign_Vector(psi_rchd)
             ti_chord     = assign_Vector(ti_rchd)
             te_chord     = assign_Vector(te_rchd)


 
 ! dump the data to file:
       IF(myid == master)THEN
         WRITE(ncrt,FMT='("            R major      psi          rho       nrho        ene          te          ti")')
         DO j=1,SIZE(r_chd)
            WRITE(ncrt,FMT='(2x,i3,2x,7(1pe12.2))') j,R_chord%data(j),psi_chord%data(j),   &
                 rho_chord%data(j),rho_chord%data(j)/rho_chord%data(dischg%nr_mhd),        &
                 ene_chord%data(j),te_chord%data(j),ti_chord%data(j) 
         ENDDO
       ENDIF

        RETURN
     END SUBROUTINE line_avg_den


     SUBROUTINE find_Rb(Rl,Rr,zelev,psisgn,psibdry_lcl,ierr)
     !--------------------------------------------------------------------------
     ! -- INPUT 
     !     zelev      elevation at whcih horizontal chod across plasma is desired
     ! -- OUTPUT
     !     Rl,Rr      inboard/outboard major radius of plasma edge at zelev
     !-------------------------------------------------------HSJ------------------
        USE nrtype,                      ONLY : Dp,I4B

        USE bicubic_spline,              ONLY : eval_2d_spline

        USE Plasma_properties,           ONLY : dischg,profile

        USE common_constants,            ONLY : zeroc,izero
 
        IMPLICIT NONE

        REAL(DP) fs(6),zelev,rguess_l,rguess_r,dr,psiguess_l,psiguess_r
        REAL(DP) Rl,Rr,Rs,psisgn,psibdry_lcl
        INTEGER(I4B) icalc,jstep,jstep_max,ierr

        icalc = 1 ; jstep_max = 2*dischg%nr_mhd ; ierr = izero



        dr = dischg%rmhdgrid%data(2)-dischg%rmhdgrid%data(1)
        rguess_l = dischg%rmhdgrid%data(1)
        CALL eval_2d_spline(rguess_l,zelev,fs,icalc)
        psiguess_l = fs(1)
        rguess_r = dischg%rmhdgrid%data(1)+dr
        jstep = izero
        DO WHILE( 1 .GT. 0) ! do until explicit exit
            jstep = jstep+1
           CALL eval_2d_spline(rguess_r,zelev,fs,icalc)
           psiguess_r   = fs(1)

           IF(psisgn .GT. 0.0)THEN
              !psi is increasing toward axis
               IF( psiguess_l .LE.  psibdry_lcl .AND.  fs(1) .GE.  psibdry_lcl)EXIT
              
           ELSE
               IF( psiguess_l .GE.  psibdry_lcl .AND.  fs(1)  .LE.  psibdry_lcl)EXIT
           ENDIF

           

           rguess_l     = rguess_r
           psiguess_l   = psiguess_r
           rguess_r     = rguess_r + dr
           IF( jstep .GT. jstep_max)ierr =1 
           IF(ierr == 1)EXIT
        ENDDO


        IF(ierr == 1)RETURN ! error in eqdsk psi proicessing skip calculation

        ! we have psi boundary between rguess_l and rguess_r
        ! refine estimate using Newton:
          IF(ABS( psiguess_r - psibdry_lcl) .LE. ABS( psiguess_l - psibdry_lcl))THEN
            Rs = rguess_r
          ELSE
            Rs = rguess_l
          ENDIF
         CALL refine_Rs(Rs,zelev,psibdry_lcl)
         Rl = Rs

! Repeat for outboard side:

        psiguess_l = fs(1)
        rguess_r = dischg%rma+dr
        jstep = izero
        DO WHILE( 1 .GT. 0) ! do until explicit exit
            jstep = jstep+1
           CALL eval_2d_spline(rguess_r,zelev,fs,icalc)
           psiguess_r   = fs(1)

           IF(psisgn .GT. 0.0)THEN
              !psi is increasing toward axis
               IF( psiguess_l .GE.  psibdry_lcl .AND.  fs(1) .LE.  psibdry_lcl)EXIT
              
           ELSE
               IF( psiguess_l .LE.  psibdry_lcl .AND.  fs(1)  .GE.  psibdry_lcl)EXIT
           ENDIF

           

           rguess_l     = rguess_r
           psiguess_l   = psiguess_r
           rguess_r     = rguess_r + dr
           IF( jstep .GT. jstep_max)ierr =1 
           IF(ierr == 1)EXIT
        ENDDO
        IF(ierr == 1)RETURN ! error in eqdsk psi proicessing skip calculation


        ! we have psi boundary between rguess_l and rguess_r
        ! refine estimate using Newton:
          IF(ABS( psiguess_r - psibdry_lcl) .LE. ABS( psiguess_l - psibdry_lcl))THEN
            Rs = rguess_r
          ELSE
            Rs = rguess_l
          ENDIF



         CALL refine_Rs(Rs,zelev,psibdry_lcl)
         Rr = Rs





     END SUBROUTINE find_Rb



     SUBROUTINE refine_Rs(Rs,zelev,psibdry_lcl)
!-------------------------------------------------------------------------------------
! -- INPUT/OUTPUT
! -- Rs current estmate on input
! --    converged solution on output
!-------------------------------------------------------------------------------------

        USE nrtype,                      ONLY : Dp,I4B

        USE bicubic_spline,              ONLY : eval_2d_spline

        USE Plasma_properties,           ONLY : dischg,profile

        USE common_constants,            ONLY : zeroc,izero

        IMPLICIT NONE
        
        INTEGER, PARAMETER :: itermx = 20
        REAL(DP) fs(6),zelev,Rs,deltar,tol,psibdry_lcl
        INTEGER(I4B) icalc,iter


           icalc =2 ; iter = izero ; tol = 1.e-8
           DO WHILE(iter .le. itermx)
               CALL eval_2d_spline(Rs,zelev,fs,icalc)
               deltar = -(fs(1)-psibdry_lcl)/fs(2)
               Rs =  Rs + deltar
               IF(ABS(deltar/Rs)  .LT. tol)EXIT
               iter =iter +1
           ENDDO




           IF(iter .GE. itermx)THEN

           ENDIF
     END SUBROUTINE refine_Rs
