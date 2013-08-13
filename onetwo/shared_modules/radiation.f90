     MODULE radiation
      

        USE nrtype,                               ONLY : I4B,DP

        USE grid_class,                           ONLY : nj,r
     
        USE ions_gcnmp,                           ONLY : namep,namei,nion,    &
                                                         nprim,nimp,nneu,z,zsq
        USE dep_var,                              ONLY : en,te,ene

        USE source_terms_gcnmp,                   ONLY : qrad,prad_imp,brems_nions, &
                                                         qrad_tot,brems_tot

        USE common_constants,                     ONLY : zeroc,kevpjou,joupkev,Pi 

        USE Plasma_properties,                    ONLY : dischg,mhd_dat
    
        USE solcon_gcnmp,                         ONLY : cyclotron_reflection

        IMPLICIT NONE
        INTEGER(I4b) j,k
  

        CONTAINS
       
        SUBROUTINE load_rad_loss
!----------------------------------------------------------------------------
!-- determine qrad,watts/m^3 energy lost due to radiation (electron channel)
! -- includes bremssthralung,cyclotron,and impurity radiation
! ---------------------------------------------------------------------------
      IMPLICIT NONE


      INTEGER(I4B)j,k

      REAL(DP) cyclo,teerg,wpsq,chi,phi,vdotsq,wbx, &
               rm,cee,xmasse,charge,btor,dene



!---------------------------------------------------------------------------
!--   calculate the radiative energy loss due to bremsstrahlung from
!--   primary ions
!- ---------------------------------------------------------------------------
         qrad(:) = zeroc ; prad_imp(:) = zeroc ; brems_nions(:,:) = zeroc
         DO j=1,nj
             DO k=1,nprim
                IF (namep(k) .EQ. 'he') THEN
                   brems_nions(j,k) = fradhe(te(j))*en(j,k)    
                   qrad(j) = qrad(j) + fradhe(te(j))*en(j,k)
                ELSE ! ok if namep ='dt'
                   brems_nions(j,k) =  Pwr_brem_factr(te(j),z(j,k))*en(j,k)
                   qrad(j) = qrad(j) + brems_nions(j,k) 
                END IF
             END DO
          END DO
         !NOTE: NOW qrad,brems_nions  are in watts per electron   (not watts/m^3 ,
         !      because still need to multiply by ene)



!----------------------------------------------------------------------
!  calculate the radiative energy loss due to bremsstrahlung,
!  radiative recombination, and line radiation from impurities
!  radfit returns prad_imp  in watts 
!  NOTE radfit sums over all impurities (nprim+1 to nion) 
!----------------------------------------------------------------------
          IF (nimp .NE. 0)                        &
               CALL radfit (prad_imp, en(1,nprim+1), te, nj)
   


! ----------------------------------------------------------------------
!  add radiative energy loss terms and convert units from
!  W/electron  to  keV/m**3-s
! ----------------------------------------------------------------------
          DO  j=1,nj
             DO k=1,nion
                brems_nions(j,k) = brems_nions(j,k)*kevpjou  * ene(j) ! kev/m^3 sec
             ENDDO
             prad_imp(j) = kevpjou  * ene(j) * prad_imp(j)   ! kev/M^3 sec
             qrad(j) = kevpjou  * ene(j) * qrad(j) + prad_imp(j)
          ENDDO
          
! ----------------------------------------------------------------------
!  add radiative energy loss due to cyclotron (synchrotron) radiation.
!  a realistic treatment of this term would require knowledge of the
!  emission and absoption characteristics of the plasma.  the present
!  calcultion is a simplified model which will probably only give
!  good global results. calculation is in CGS (Gaussian)
!  then converted to keV/cm**3-s.
!  ref: trubnikov: jetp letters 16, 25 (1972)
!
!  cyclotron_reflection == wall reflection coeffiecient, If cyclotron_reflection =1 then
!  all the synchrotraon radiation is assumed reabsorbed so
!  no addition to qrad is made:
! ----------------------------------------------------------------------

          cyclo = 0.0_DP
          IF (cyclotron_reflection .LT. 1.0_DP)THEN
             charge = 4.8032067e-10_DP        ! electron charge (esu)
             xmasse = 9.1093897e-28_DP        ! electron mass (g)
             cee    = 2.99792458e10_DP        ! speed of light (cm/s)
             rm     = r(nj)/dischg%rmajor     ! 
             btor   = mhd_dat%btor*1.e4       ! gauss
             wbx     = charge * ABS (btor) / (xmasse * cee)
             DO j=1,nj
                dene    = ene(j)*1.e-6_DP
                teerg   = te(j)*1.6e-9
                wpsq    = 4.0_DP * pi*dene*charge**2/xmasse
 !               chi     = r(nj)/rmajor * SQRT (xmasse*cee**2 /teerg)
                chi     = rm  * SQRT (xmasse*cee**2 /teerg)
                phi     = 60.0_DP*(teerg/(xmasse*cee**2))**1.5 &
                     * SQRT (cee*wbx/(r(nj)*100._DP*wpsq)*(1.0_DP-cyclotron_reflection)*(1.0_DP+chi))
                vdotsq  = wbx**2*2.0_DP*teerg/xmasse
                cyclo   = dene*1.5_DP*charge**2/cee**3*vdotsq*phi/1.6e-9_DP
                cyclo   = cyclo *1.e6_DP  ! kev/(m^3 sec)
                qrad(j) = qrad(j) + cyclo 
             END DO
          ENDIF


       RETURN

   END  SUBROUTINE load_rad_loss






      SUBROUTINE radfit (prad_imp, imp_den, te, nj)
! ----------------------------------------------------------------------
!     this subroutine computes prad_imp, the impurity radiation per
!     electron in watts summed over all imurities.
!     Individual impurity species are returned in  brems_nions
!     coeffficents for the curve fits used were
!     obtained from a Princeton PPL paper, "steady state radiative
!     cooling rates for low-density high-temperature plasmas"
!     (pppl-1352)
!
!     changed logic to use largest energy range in the table
!     instead of stopping code with error message.  HSJ 09/11/06
!
!     name(i)    = chemical abbreviation for element
!     coeff(i,1) = lower temperature range
!     coeff(i,2) = upper temperature range
!     coeff(i,3) ---> coeff(i,8) are the coefficients of a 5th degree
! -----------------------------------------------------------------------

          USE error_handler,               ONLY : terminate,lerrno,iomaxerr
          USE io_gcnmp,                    ONLY : ncrt,nlog

          IMPLICIT  NONE



          INTEGER, PARAMETER   ::     krows = 49
          INTEGER, PARAMETER   ::     nrows = krows
          CHARACTER*2 name(krows)
          REAL(DP)    coeff(krows,8)
          REAL(DP)    prad_imp(*), imp_den(nj,*),te(*)
          REAL(DP) x,rad,y
          INTEGER(I4B) imp,i,ki,j,nj,i0,irow
          REAL(DP) tkev
         

! ----------------------------------------------------------------------
! *********  argon             **********
          DATA name( 1) /'ar'/, (coeff( 1,i),i = 1,8)/ &
               3.000e-02,2.000e-01,-2.05304e+01,-2.83429e+00, &
               1.50690e+01, 3.51718e+01, 2.40012e+01, 5.07272e+00/
          DATA name( 2) /'ar'/, (coeff( 2,i),i = 1,8)/ &
               2.000e-01,2.000e+00,-1.96520e+01,-1.17276e-01, &
               7.83322e+00,-6.35158e+00,-3.05885e+01,-1.52853e+01/
          DATA name( 3) /'ar'/, (coeff( 3,i),i = 1,8)/ &
               2.000e+00,2.000e+01,-1.97488e+01, 2.96484e+00, &
               -8.82939e+00, 9.79100e+00,-4.96002e+00, 9.82003e-01/
          DATA name( 4) /'ar'/, (coeff( 4,i),i = 1,8)/ &
               2.000e+01,1.000e+02,-2.11794e+01, 5.19148e+00, &
               -7.43972e+00, 4.96902e+00,-1.55318e+00, 1.87705e-01/
! *********  carbon            **********
          DATA name( 5) /'c '/, (coeff( 5,i),i = 1,8)/ &
               2.000e-03,2.000e-02, 2.26196e+03, 5.27071e+03, &
               4.80768e+03, 2.16700e+03, 4.83472e+02, 4.27841e+01/
          DATA name( 6) /'c '/, (coeff( 6,i),i = 1,8)/ &
               2.000e-02,2.000e-01, 5.01266e+01, 3.32680e+02, &
               6.00316e+02, 5.17230e+02, 2.13035e+02, 3.36380e+01/
          DATA name( 7) /'c '/, (coeff( 7,i),i = 1,8)/ &
               2.000e-01,2.000e+00,-2.12371e+01,-3.13175e-01, &
               7.60470e-01,-3.01657e-01, 7.63153e-02, 1.40114e-01/
          DATA name( 8) /'c '/, (coeff( 8,i),i = 1,8)/ &
               2.000e+00,2.000e+01,-2.12367e+01,-3.44789e-01, &
               1.03546e+00,-1.01249e+00, 5.92047e-01,-1.43559e-01/
          DATA name( 9) /'c '/, (coeff( 9,i),i = 1,8)/ &
               2.00000e+01, 1.00000e+04, -2.47680e+01, 9.40818e+00, &
               -9.65745e+00, 4.99916e+00, -1.23738e+00, 1.16061e-01/
! *********  chromium          **********
          DATA name(10) /'cr'/, (coeff(10,i),i = 1,8)/ &
               2.000e-02,2.000e-01,-1.04662e+01, 4.84777e+01, &
               1.03531e+02, 9.96556e+01, 4.53392e+01, 7.96357e+00/
          DATA name(11) /'cr'/, (coeff(11,i),i = 1,8)/ &
               2.000e-01,2.000e+00,-1.85628e+01,-2.86143e+00, &
               -7.07898e+00, 1.04095e+01, 4.27995e+01, 2.95598e+01/
          DATA name(12) /'cr'/, (coeff(12,i),i = 1,8)/ &
               2.000e+00,2.000e+01,-1.70389e+01,-1.81938e+01, &
               5.07803e+01,-6.45257e+01, 3.81862e+01,-8.58752e+00/
          DATA name(13) /'cr'/, (coeff(13,i),i = 1,8)/ &
               2.000e+01,1.000e+02,-1.41239e+01,-1.39958e+01, &
               1.49161e+01,-8.16669e+00, 2.30136e+00,-2.62443e-01/
! *********  iron              **********
          DATA name(14) /'fe'/, (coeff(14,i),i = 1,8)/ &
               2.000e-03,2.000e-02, 5.40647e+02, 1.31358e+03, &
               1.22549e+03, 5.67632e+02, 1.30545e+02, 1.19286e+01/
          DATA name(15) /'fe'/, (coeff(15,i),i = 1,8)/ &
               2.000e-02,2.000e-01,-1.19062e+01, 2.59325e+01, &
               4.12580e+01, 3.03780e+01, 1.04477e+01, 1.39669e+00/
          DATA name(16) /'fe'/, (coeff(16,i),i = 1,8)/ &
               2.000e-01,2.000e+00,-1.82996e+01,-1.55942e+00, &
               -6.00441e+00,-1.19285e-01, 1.70088e+01, 1.17635e+01/
          DATA name(17) /'fe'/, (coeff(17,i),i = 1,8)/ &
               2.000e+00,2.000e+01,-1.67982e+01,-1.59773e+01, &
               3.82500e+01,-4.25281e+01, 2.22674e+01,-4.46320e+00/
          DATA name(18) /'fe'/, (coeff(18,i),i = 1,8)/ &
               2.000e+01,1.000e+02,-2.45396e+01, 1.79522e+01, &
               -2.35636e+01, 1.48450e+01,-4.54232e+00, 5.47746e-01/
! *********  helium            **********
          DATA name(19) /'he'/, (coeff(19,i),i = 1,8)/ &
               2.000e-03,1.000e-02, 3.84322e+03, 8.93072e+03, &
               8.17947e+03, 3.71287e+03, 8.35739e+02, 7.46792e+01/
          DATA name(20) /'he'/, (coeff(20,i),i = 1,8)/ &
               1.000e-02,2.000e-01,-2.25831e+01, 1.15730e-02, &
               -8.32355e-01,-1.17916e+00,-4.74033e-01,-8.64483e-02/
          DATA name(21) /'he'/, (coeff(21,i),i = 1,8)/ &
               2.000e-01,2.000e+00,-2.25560e+01, 3.28276e-01, &
               1.39226e-01,-1.22085e-01,-2.76602e-01,-2.90494e-01/
          DATA name(22) /'he'/, (coeff(22,i),i = 1,8)/ &
               2.000e+00,2.000e+01,-2.25798e+01, 4.57017e-01, &
               -1.59274e-01, 2.71952e-01,-1.71810e-01, 3.86609e-02/
          DATA name(23) /'he'/, (coeff(23,i),i = 1,8)/ &
               2.000e+01,1.000e+02,-1.73046e+01,-1.62462e+01, &
               2.10079e+01,-1.31207e+01, 4.06935e+00,-5.00944e-01/
! *********  krypton           **********
          DATA name(24) /'kr'/, (coeff(24,i),i = 1,8)/ &
               5.000e-02,2.000e-01, 5.67564e+01, 4.00852e+02, &
               8.51654e+02, 8.89400e+02, 4.54768e+02, 9.12158e+01/
          DATA name(25) /'kr'/, (coeff(25,i),i = 1,8)/ &
               2.000e-01,2.000e+00,-1.81364e+01,-1.51898e+00, &
               2.83187e+00, 6.12136e+00,-1.07998e+01,-1.57415e+01/
          DATA name(26) /'kr'/, (coeff(26,i),i = 1,8)/ &
               2.000e+00,2.000e+01,-2.33830e+01, 3.90879e+01, &
               -1.05886e+02, 1.26461e+02,-7.04126e+01, 1.50163e+01/
          DATA name(27) /'kr'/, (coeff(27,i),i = 1,8)/ &
               2.000e+01,1.000e+02,-2.63927e+01, 1.95442e+01, &
               -2.06960e+01, 1.09324e+01,-2.88191e+00, 3.06312e-01/
! *********  molybdenum        **********
          DATA name(28) /'mo'/, (coeff(28,i),i = 1,8)/ &
               6.000e-02,2.000e-01,-1.39105e+02,-6.49334e+02, &
               -1.36584e+03,-1.40646e+03,-7.08621e+02,-1.40057e+02/
          DATA name(29) /'mo'/, (coeff(29,i),i = 1,8)/ &
               2.000e-01,2.000e+00,-1.77259e+01,-1.05822e+00, &
               -3.58317e+00, 1.66009e+00, 8.56537e+00, 4.53291e+00/
          DATA name(30) /'mo'/, (coeff(30,i),i = 1,8)/ &
               2.000e+00,2.000e+01,-1.38510e+01,-3.67845e+01, &
               1.14059e+02,-1.63563e+02, 1.07626e+02,-2.64249e+01/
          DATA name(31) /'mo'/, (coeff(31,i),i = 1,8)/ &
               2.000e+01,1.000e+02, 3.99268e+01,-1.75709e+02, &
               2.07493e+02,-1.21459e+02, 3.53180e+01,-4.08383e+00/
! *********  nickel            **********
          DATA name(32) /'ni'/, (coeff(32,i),i = 1,8)/ &
               3.000e-02,2.000e-01,-1.20325e+01, 3.25391e+01, &
               6.79077e+01, 6.52992e+01, 2.97346e+01, 5.27128e+00/
          DATA name(33) /'ni'/, (coeff(33,i),i = 1,8)/ &
               2.000e-01,2.000e+00,-1.83048e+01,-3.31924e-03, &
               -3.33231e+00,-1.11280e+01, 1.05307e-01, 9.44891e+00/
          DATA name(34) /'ni'/, (coeff(34,i),i = 1,8)/ &
               2.000e+00,2.000e+01,-1.69768e+01,-9.49547e+00, &
               1.10936e+01, 4.04590e-02,-6.52193e+00, 2.65491e+00/
          DATA name(35) /'ni'/, (coeff(35,i),i = 1,8)/ &
               2.000e+01,1.000e+02,-2.86408e+01, 2.99929e+01, &
               -3.72608e+01, 2.25806e+01,-6.71660e+00, 7.91169e-01/
! *********  oxygen            **********
          DATA name(36) /'o '/, (coeff(36,i),i = 1,8)/ &
               2.000e-03,2.000e-02,-4.52975e+02,-9.70819e+02, &
               -8.54356e+02,-3.70311e+02,-7.89495e+01,-6.59517e+00/
          DATA name(37) /'o '/, (coeff(37,i),i = 1,8)/ &
               2.000e-02,2.000e-01,-5.85047e+01,-1.59971e+02, &
               -2.41395e+02,-1.60654e+02,-4.43559e+01,-3.33047e+00/
          DATA name(38) /'o '/, (coeff(38,i),i = 1,8)/ &
               2.000e-01,2.000e+00,-2.07261e+01,-7.67836e-01, &
               8.75731e-01,-1.08357e+00, 2.41126e+00, 5.75616e+00/
          DATA name(39) /'o '/, (coeff(39,i),i = 1,8)/ &
               2.000e+00,2.000e+01,-2.06891e+01,-1.09982e+00, &
               2.00613e+00,-1.58614e+00, 6.69091e-01,-1.11969e-01/
          DATA name(40) /'o '/, (coeff(40,i),i = 1,8)/ &
               2.000e+01,1.000e+02,-2.78060e+01, 2.14606e+01, &
               -2.66591e+01, 1.67083e+01,-5.19194e+00, 6.41029e-01/
! *********  silicon           **********
          DATA name(41) /'si'/, (coeff(41,i),i = 1,8)/ &
               2.000e-03,1.800e-02, 1.55215e+03, 3.32531e+03, &
               2.75818e+03, 1.12074e+03, 2.23142e+02, 1.74141e+01/
          DATA name(42) /'si'/, (coeff(42,i),i = 1,8)/ &
               1.800e-02,2.000e-01,-3.33775e+01,-4.15621e+01, &
               -3.35750e+01,-3.86101e-01, 9.55951e+00, 2.92688e+00/
          DATA name(43) /'si'/, (coeff(43,i),i = 1,8)/ &
               2.000e-01,2.000e+00,-1.94215e+01,-1.56426e-01, &
               -5.29353e+00,-1.24578e+00, 2.42390e+01, 1.86013e+01/
          DATA name(44) /'si'/, (coeff(44,i),i = 1,8)/ &
               2.000e+00,2.000e+01,-1.92475e+01,-1.98358e+00, &
               8.42427e-01, 7.54957e-01,-6.74371e-01, 1.46365e-01/
          DATA name(45) /'si'/, (coeff(45,i),i = 1,8)/ &
               2.000e+01,1.000e+02,-1.91107e+01,-2.66860e+00, &
               2.43841e+00,-1.00956e+00, 2.21647e-01,-2.07805e-02/
! *********  tungsten          **********
          DATA name(46) /'w '/, (coeff(46,i),i = 1,8)/ &
               1.000e-01,2.000e-01, 5.34083e+00, 1.56088e+02, &
               4.17170e+02, 5.50258e+02, 3.56758e+02, 9.04279e+01/
          DATA name(47) /'w '/, (coeff(47,i),i = 1,8)/ &
               2.000e-01,2.000e+00,-1.72389e+01, 5.42375e-02, &
               -1.22107e+00, 4.41181e-01,-4.48582e+00,-7.83614e+00/
          DATA name(48) /'w '/, (coeff(48,i),i = 1,8)/ &
               2.000e+00,2.000e+01,-1.47488e+01,-1.43954e+01, &
               2.10585e+01,-4.39475e+00,-1.10601e+01, 5.61699e+00/
          DATA name(49) /'w '/, (coeff(49,i),i = 1,8)/ &
               2.000e+01,1.000e+02,-2.62426e+02, 7.12559e+02, &
               -8.25017e+02, 4.74241e+02,-1.35517e+02, 1.54189e+01/

          
! ----------------------------------------------------------------------

          DO  j=1,nj ! loop over grid points
             tkev    = te(j)
             prad_imp(j) = 0.0_DP
             DO  imp=1,nimp
                !
                DO i=1,nrows
                   ki = i
                   IF (name(i) .EQ. namei(imp))  go to 2010
                END DO
                lerrno = iomaxerr + 183
                WRITE  (ncrt, 10) namei(imp)
                WRITE  (nlog, 10) namei(imp)
 10             FORMAT (' FATAL ERROR in RADFIT' / &
                     ' impurity ', a8, ' has not been added to the code' /)
                CALL terminate(lerrno,nlog)
!-----------------------------------------------------------------------
!     label 2010 takes care of lower limit.
!     Implicit assumption here is that table is
!     appropriately ordered(eg first entry for each impurity
!     is the lowest energy one )!!!!! HSJ
!     Use this lower limit even if tkev is less than c(i,1)
!-----------------------------------------------------------------------

2010            i = ki
                IF (tkev .GE. coeff(i,1))  go to 2050
                irow = i
                y    = coeff(irow,8)
                x    = LOG10 (coeff(irow,1))
                DO i=8,4,-1
                   y  = y*x + coeff(irow,i-1)
                END DO
                rad  = 10.0**y
                rad  = tkev*(rad/coeff(irow,1))
                go to 5000
                !
2050            i0 = i
                DO i=i0,nrows
                   ki = i
                   IF (tkev .GE. coeff(i,1) .AND. tkev .LE. coeff(i,2)) go to 2100
                   !        if (name(i) .ne. namei(imp))  go to 9000

                   IF (name(i) .NE. namei(imp))THEN
                      ! back up to last valid range for impurity i, HSJ,09/11/06:
                      ki = ki-1
                      WRITE(nlog,8020)name(ki), namei(imp), tkev, coeff(ki,1), &
                           coeff(ki,2)
                      WRITE(ncrt,8020)name(ki), namei(imp), tkev, coeff(ki,1), &
                           coeff(ki,2)

                      GO TO 2100
                   ENDIF
                END DO
                go to 9000
                !
2100            irow = ki
                y    = coeff(irow,8)
                x    = LOG10 (tkev)
                DO i=8,4,-1
                   y  = y*x + coeff(irow,i-1)
                END DO
                rad  = 10.0**y   ! think this is in erg/(sec  cm^3)  HSJ
                !
!5000                brems_nions(j,nprim+imp) =  rad*imp_den(j,imp)*1.0e-7
                    ! 1.e-13 = 1e-7 joules/erg X 1e-6 M^3/cm^3:
5000                brems_nions(j,nprim+imp) =  rad*imp_den(j,imp)*1.0e-13_DP
                !prad_imp (j) = prad_imp(j) + rad*imp_den(j,imp)*1.0e-7
                prad_imp(j) = prad_imp(j) + brems_nions(j,nprim+imp) 
              ENDDO ! over impurities
            ENDDO   ! over grid points
            RETURN
                
! ----------------------------------------------------------------------
! fatal errors
! ----------------------------------------------------------------------
                !
9000            WRITE  (nlog, 8010)  name(i), namei(imp), tkev, coeff(i,1), &
                     coeff(i,2)
                WRITE  (ncrt, 8010)  name(i), namei(imp), tkev, coeff(i,1), &
                     coeff(i,2)
8010            FORMAT (' FATAL ERROR in RADFIT'                              / &
                     ' you have exceeded temperature range for curve fits' / &
                     ' name(i), namei(imp)', 2x, a, 2x, a                  / &
                     ' tkev, lower limit, upper limit =', 3(2x, 1pe12.4))
                lerrno = iomaxerr + 184
                CALL terminate(lerrno,nlog)

                !
8020            FORMAT ('WARNING ERROR in RADFIT'                              / &
                     ' you have exceeded temperature range for curve fits' / &
                     ' name(i), namei(imp)', 2x, a, 2x, a                  / &
                     ' tkev, lower limit, upper limit =', 3(2x, 1pe12.4),  / &
                     ' code will continue by using largest energy  value ', / &
                     ' in the tables')
      END  SUBROUTINE radfit




      REAL(DP)  FUNCTION fradhe (temkev)
!------------------------------------------------------------------------
!
!     this function calculates the radiation losses from helium as a
!     function of elec. temp. (keV) using a formula derived by
!     George Hopkins and John M. Rawls (ga-a14141 impurity radiations)
!
!     temkev = temperature of electrons in keV
!     value returned is in units of  (watts-m**3)
!-------------------------------------------------------------------------


  
      IMPLICIT NONE
      REAL(DP) temkev,x,y,sqrtem

      x = temkev
      x = MIN(100._DP,temkev)
      x= MAX(0.01_DP,x)
      IF ( x .GT.   0.4_DP )  go to 10
      x = LOG10 (x)
      y = 10.0_DP**(-35.632_DP+x*(0.473114_DP+x*0.521255_DP))
      go to 20
   10 sqrtem = SQRT (x)
      y      = 2.136e-36_DP*sqrtem+2.08e-37_DP/sqrtem+4.8e-38_DP/(x**1.34_DP)
   20 fradhe  = y
      RETURN
!
      END FUNCTION fradhe



      REAL(DP)  FUNCTION Pwr_brem_factr(Te, Zi)
!----------------------------------------------------------------------------------
!-- this is power radiated in bremsstrhalung divided by electron and ion density
!-- eg bremssthralung power(watts/m^3)  = Pwr_brem_factr *ne*ni
!-- peculiar way of doing this is due to onetwo compatibility
!-- See Wesson,Tokamaks,pg 100 eq 4.9.2, (1987 )
!--------------------------------------------------------------------------HJ------ 
      IMPLICIT NONE
      REAL(DP) Te,Zi  ! Te is in kev in this formula

        pwr_brem_factr = 5.35E-37_DP*SQRT(Te)*Zi**2

       RETURN

      END FUNCTION Pwr_brem_factr

      SUBROUTINE tot_radiation
!----------------------------------------------------------------------------------
!     get the volume integrated values of qrad and brems_nions(1:nion)
!     INPUT   assumed units are kev/(m^3 sec)
!       qrad
!       brems_nions
!       r (r is rho grid)
!       hcap
!       volfac
! 
!     OUTPUT
!       qrad_tot watts
!       brems_tot(i:nion) watts 
!----------------------------------------------------------------------------------
      USE grid_class,                               ONLY : nj,r,hcap,volfac

      USE ions_gcnmp,                               ONLY : nion


      INTEGER(I4B) k

        CALL trapv(r,qrad,hcap,nj,qrad_tot)
        qrad_tot = qrad_tot*volfac *joupkev
        IF(.NOT. ALLOCATED(brems_tot))ALLOCATE(brems_tot(nion))
        DO k=1,nion
           CALL trapv(r,brems_nions(1,k),hcap,nj,brems_tot(k))
           brems_tot(k) = brems_tot(k)*volfac *joupkev
        ENDDO

        RETURN
        END SUBROUTINE tot_radiation


     END MODULE radiation
