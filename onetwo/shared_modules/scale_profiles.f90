    SUBROUTINE scale_den
!----------------------------------------------------------------------------
!-- Scale input density profiles: profile%ene, profile%en(1:nion)
!-- use erf functions to get edge shape
!----------------------------------------------------------------------------
        USE  nrtype,                                   ONLY : DP,I4B,I2B
        USE  Plasma_properties ,                       ONLY : profile
        USE  ions_gcnmp,                               ONLY : nion,ni_sc
        USE  common_constants,                         ONLY : SQRT2 ,zeroc 
        USE  grid_class,                               ONLY : nj,rho_gridn
    
   USE MPI_data,               ONLY : myid,numprocs,master,mpiierr !temp for debug
        IMPLICIT none
        REAL(DP) nifactr,one,arg,myerfc,derf,derfc
        INTEGER(I4B) j,k
        

        If(  ni_sc(1) .LT. zeroc .OR. ni_sc(2) .LT. zeroc  .OR.   &
             ni_sc(3) .LT. zeroc .OR. ni_sc(4) .LT. zeroc) RETURN

        DO j=1,nj
           arg = (rho_gridn%data(j)-ni_sc(2))/(SQRT2*ni_sc(3))
           myerfc  = 2._DP - derfc(arg)
           myerfc = MAX(zeroc,myerfc)
           myerfc = MIN(2._DP,myerfc)
           nifactr = 1._DP + (ni_sc(1)-1._DP)*(1._DP-0.5_DP*myerfc)
           !nifactr = ni_sc(1)*(1._DP - 0.5_DP*(1._DP+derf((rho_gridn%data(j)-ni_sc(2))/(SQRT2*ni_sc(3))))+1._DP/ni_sc(1))
           IF(rho_gridn%data(j) .ge. ni_sc(4)) EXIT
           profile%ene%data(j) = nifactr*profile%ene%data(j)
           DO k=1,nion
             profile%en(k)%data(j) = nifactr*profile%en(k)%data(j)
           ENDDO
        ENDDO

    RETURN 
    END SUBROUTINE scale_den


    SUBROUTINE scale_betan_ped
!----------------------------------------------------------------------------
!-- Scale input density profiles: profile%ene, profile%en(1:nion)
!-- to get desired betan at pedestal point ==> rho_gridn,rhon_betan_ped
!----------------------------------------------------------------------------
        USE  nrtype,                                   ONLY : DP,I4B,I2B
        USE  Plasma_properties ,                       ONLY : profile,mhd_dat,dischg
        USE  ions_gcnmp,                               ONLY : nion
        USE  common_constants,                         ONLY : izero,zeroc,Permeability 
        USE  grid_class,                               ONLY : nj,rho_gridn,dr,r,j_ped_index
        USE  vector_class,                             ONLY : zero_vector
        USE  bc_values_gcnmp,                          ONLY : fix_edge_te_bc,fix_edge_ti_bc,  &
                                                              fix_edge_rot_bc,fix_edge_ni_bc, &
                                                              te_var_edge,ti_var_edge,        &
                                                              rot_var_edge,ni_var_edge,       &
                                                              te_index,ti_index,rot_index,    & 
                                                              ni_index
        USE error_handler,                             ONLY : iomaxerr, lerrno,terminate

        USE io_gcnmp,                                  ONLY : nlog

        USE MPI_data,                                  ONLY : myid,master

        IMPLICIT none
        REAL(DP) nifactr,one,arg,myerfc,derf,derfc,den_factor,old_den,dra
        REAL(DP) delta_ni(nion),delta_ne
        INTEGER(I4B) j,k,j_ped,fi_inc,p_type,error
        LOGICAL conserve_zeff,optiona


        optiona = .FALSE.                   ! false means scale primary  ions only
                                            ! true tries to get self consistent set
                                            ! for electrons and al ions

        fi_inc = 0                          ! =1 ==> include fast ions in tot_pressure
        p_type = 1                          ! use profile%... instead of ene,etc.

        CALL get_betan(fi_inc,p_type)       ! get local betan across entire  grid


       !find r that gets closest to rhon_betan_ped (r(j_ped):
       j_ped = izero ; 
       DO j =1,nj
          j_ped=j
          dra = 0.5_DP*dr(j)/r(nj)
          if(ABS(rho_gridn%data(j)- mhd_dat%rhon_betan_ped) .LT. dra)EXIT
       ENDDO
       IF(j_ped == nj)THEN
         IF(myid == master)THEN
            PRINT *,'ERROR: You have selected betan_ped density scaling '
            PRINT *,'       But rhon_betan_ped is out of range (must be less than 1.0)'
            PRINT *,'       rhon_betan_ped input =',mhd_dat%rhon_betan_ped
         ENDIF
         lerrno = iomaxerr + 157
         CALL terminate(lerrno,nlog)
       ENDIF

       ! we will scale the density out to grid point j_ped to make betan_ped = mhd_dat%betan_ped
       ! at that point.
       den_factor =  mhd_dat%betan_ped/mhd_dat%betan%data(j_ped)
 
       old_den = profile%en(1)%data(j_ped)
       conserve_zeff = .TRUE.
       CALL get_den_ptrb(den_factor,j_ped,conserve_zeff,error,delta_ni,delta_ne)
       CALL mod_profiles(delta_ni,delta_ne,j_ped,error)      ! change ne,ni nimp  at grid pt j_ped to ge betanped

       IF(error .GT. izero .AND. optiona )THEN
          conserve_zeff = .FALSE.                            ! this option not currently used (optiona = false)
          CALL get_den_ptrb(den_factor,j_ped,conserve_zeff,error,delta_ni,delta_ne)
          CALL mod_profiles(delta_ni,delta_ne,j_ped,error) 

       ELSEIF(error .GT. izero .AND.  .NOT. optiona )THEN
          conserve_zeff = .FALSE.                            ! perturb only first two primary  ion species
                                                             ! leave impurity alone, recal zeff ene below
          CALL get_den_ptrb_prim(den_factor,j_ped,error,delta_ni,delta_ne)
          CALL mod_profiles(delta_ni,delta_ne,j_ped,error)   ! note that only grid point j_ped is done 
       ENDIF
 
       IF(error .GT. izero)THEN       ! modified densities are negative
           IF(myid == master)THEN
                PRINT *,'ERROR in sub scale_betan_ped'
                PRINT *,'Even if zeff is not conserved the solution'
                PRINT *,'for some of the  scaled densities was negative'
                WRITE(*,FMT='("Input Densities, required perturbed density")')
                DO j=1,nion
                   WRITE(*,FMT='(1pe12.4 ,8x,1pe12.4)')profile%en(j)%data(j_ped),profile%en(j)%data(j_ped)+delta_ni(j)
                ENDDO
                WRITE(*,FMT='(1pe12.4 ,8x,1pe12.4)')profile%ene%data(j_ped),profile%ene%data(j_ped)+ delta_ne
                WRITE(*,FMT='("At  rho =",1pe12.6," current betan is ",1pe12.4)')rho_gridn%data(j_ped),mhd_dat%betan%data(j_ped)
               WRITE(*,FMT='("requested rho = ",1pe12.6," requested betan is ",1pe12.4)')mhd_dat%rhon_betan_ped,mhd_dat%betan_ped
                WRITE(*,FMT='("Suggested solution is to change location of pedestal",/ &
                              "and/or change Te,Ti values on the pedestal")')
           ENDIF
           lerrno = iomaxerr + 158
           CALL terminate(lerrno,nlog)
       ENDIF

       IF(myid == master)THEN
 

          WRITE(nlog,FMT="('Rescaled the densities to give desired betan_ped:',/,       &
                         '    Nearest rho grid value           = ',1pe12.6,/,           &
                         '    Old value of betan_ped           = ',1pe12.6,/,           &
                         '    Old value of primary ion density = ',1pe12.6,/,           &
                         '    New value betan_ped              = ',1pe12.6,/,           &
                         '    New value of primary ion density = ',1pe12.6,/,           &
                         '    Scale factor                     = ',1pe12.6,/,           &
                         '    Conserved zeff                   = ',l6 )")               &
                    rho_gridn%data(j_ped),mhd_dat%betan%data(j_ped),old_den,            &
                    mhd_dat%betan_ped,profile%en(1)%data(j_ped),den_factor,conserve_zeff

           IF(conserve_zeff)THEN
               WRITE(nlog,FMT="('  primary , impurity and electron densities were scaled')") 
           ELSE
              WRITE(nlog,FMT="('  only primary  ions were scaled')") 
           ENDIF
  
       ENDIF
        CALL get_betan(fi_inc,p_type)
        j_ped_index = j_ped                      ! store in grid_class

    RETURN 
    END SUBROUTINE scale_betan_ped



    SUBROUTINE get_den_ptrb(den_factor,j_ped,conserve_zeff,error,delta_ni,delta_ne)
!---------------------------------------------------------------------------------
! -- GET Change in density profiles at constant  Zeff
! -- (conserve_zeff = true) to get betan_ped
! -- if conserve_zeff = false then let zeff float
!----------------------------------------------------------------------------------
        USE  nrtype,                                    ONLY : DP,I4B,I2B

        USE  Plasma_properties ,                        ONLY : profile,mhd_dat,dischg

        USE  ions_gcnmp,                                ONLY : nion,zeff,z,zsq

        USE  common_constants,                          ONLY : zeroc,izero,joupkev

        USE error_handler,                              ONLY : iomaxerr, lerrno,terminate

        USE io_gcnmp,                                   ONLY : nlog



        IMPLICIT NONE
        LOGICAL conserve_zeff
        INTEGER(I4B) j,neq,j_ped,nrhs,k,info,error,ul
        INTEGER(I4B) ipiv(nion)
        REAL(DP) sum_den,sum_zeff,zeff_ped,delta_ne,delta_ptot,den_factor, &
                 Ti_ped,Te_ped,ratio
        REAL(DP) a(nion,nion),rhs(nion),delta_ni(nion)



         error = izero
         delta_ni(:) = zeroc  ; delta_ne = zeroc
        ! Set up equations to find the change in ion densities for each ion
        ! species. retaining  the given zeff. 

        delta_ptot = (den_factor-1._DP)*profile%tot_press%data(j_ped) ! change in pressure required

        rhs(1) = delta_ptot
        neq = nion        ! # equations to solve for change in ion density at pedestal
        sum_den = zeroc  ; sum_zeff = zeroc 
        Ti_ped  = profile%Ti%data(j_ped)*joupkev
        Te_ped  = profile%Te%data(j_ped)*joupkev
        zeff_ped  = zeff(j_ped)
        DO j=1,nion
           sum_den  = sum_den + profile%en(j)%data(j_ped)
           a(1,j) = Ti_ped + Te_ped*z(j_ped,j)                ! first equation
           sum_zeff =sum_zeff + profile%en(j)%data(j_ped)*(zsq(j_ped,j)-zeff_ped*z(j_ped,j))
        ENDDO

        ! last equation conserves zeff
        IF(conserve_zeff)THEN
           rhs(neq) = sum_zeff/zeff(j_ped)
           DO j=1,nion
              a(neq,j) = z(j_ped,j)
           ENDDO
           ul = nion -2
        ELSE
           ul = nion -1      ! zeff eq not used so use ratio of prim to impurity instead below
        ENDIF

        ! several intermediary equations
        
        IF(nion .GT. 2)THEN  
           DO j=1,ul                                      
              ratio  = profile%en(j)%data(j_ped)/sum_den
              rhs(j+1) = zeroc
              DO k=1,neq
                 a(j+1,k) = 1._Dp
                 IF(k == j) a(j+1,k) = ratio - a(j+1,k)
              ENDDO
           ENDDO
        ENDIF


        IF( nion == 1)THEN  
           ! use single equation generated above
           ! Note that this means zeff  =1 since single species primary ion
           ! assumes  fully stripped charge state
           delta_ni(1) = rhs(1)/(Ti_ped + Te_ped*z(j_ped,1))
           delta_ne    = delta_ni(1)*z(j_ped,1)

        ELSE
           ! solve  a*delta_ni = rhs for delta_ni:
             nrhs =1 
             CALL DGETRF( neq, neq, a,neq,ipiv,info )
             CALL DGETRS('N', neq,nrhs, a,neq, ipiv, rhs,neq,info)
             IF(info .NE. izero)THEN
                lerrno = iomaxerr + 159
                CALL terminate(lerrno,nlog)
             ENDIF
             delta_ni(:) = rhs(:)
             delta_ne = zeroc
             Do j=1,nion
                delta_ne = delta_ne + delta_ni(j)*z(j_ped,j)
             ENDDO
        ENDIF



        RETURN

    END SUBROUTINE get_den_ptrb



    SUBROUTINE get_den_ptrb_prim(den_factor,j_ped,error,delta_ni,delta_ne)
!----------------------------------------------------------------------------
! -- perturb only first two primary ions species
! -- Set up equations to find the change in ion densities for up to 2 primary ions
!----------------------------------------------------------------------------
        USE  nrtype,                                    ONLY : DP,I4B,I2B

        USE  Plasma_properties ,                        ONLY : profile,mhd_dat,dischg

        USE  ions_gcnmp,                                ONLY : nion,nprim,zeff,z,zsq

        USE  common_constants,                          ONLY : zeroc,izero,joupkev

        USE error_handler,                              ONLY : iomaxerr, lerrno,terminate

        USE io_gcnmp,                                   ONLY : nlog

        IMPLICIT NONE

        INTEGER(I4B) j,neq,j_ped,nrhs,k,info,error,ul
        INTEGER(I4B) ipiv(nion)
        REAL(DP) zeff_ped,delta_ne,delta_ptot,den_factor, &
                 Ti_ped,Te_ped,ratio
        REAL(DP) delta_ni(nion)

         error = izero
         delta_ni(:) = zeroc  ; delta_ne = zeroc
         delta_ptot = (den_factor-1._DP)*profile%tot_press%data(j_ped) ! change in pressure required
         Ti_ped  = profile%Ti%data(j_ped)*joupkev
         Te_ped  = profile%Te%data(j_ped)*joupkev
         IF(nprim .GE. 2)THEN
            ratio = profile%en(2)%data(j_ped)/profile%en(1)%data(j_ped)
         ELSE
            ratio = zeroc
         ENDIF

         delta_ni(1) = delta_ptot/(Ti_ped*(1._DP+ratio) +       &
                                   Te_ped*(z(j_ped,1)+ratio*z(j_ped,2)))
         IF(nprim .GE. 2)delta_ni(2) = ratio*delta_ni(1)
         delta_ne    = delta_ni(1)*z(j_ped,1)+delta_ni(2)*z(j_ped,2)

         RETURN

    END SUBROUTINE get_den_ptrb_prim


    SUBROUTINE mod_profiles(delta_ni,delta_ne,j_ped,error)
!----------------------------------------------------------------------------
! -- return perturbed profiles  (at grid point j_ped) if allpositive.otherwise
! -- return with input input profiles unchanged
!----------------------------------------------------------------------------
        USE  nrtype,                                   ONLY : DP,I4B,I2B
        USE  Plasma_properties ,                       ONLY : profile,mhd_dat,dischg
        USE  ions_gcnmp,                               ONLY : nion
        USE  common_constants,                         ONLY : izero


        IMPLICIT NONE
        REAL(DP) save_pt,crit_den,delta_ni(nion),delta_ne,den_save(nion+1)
        INTEGER(I4B) j_ped,j,error

        crit_den = 1.e15                              ! arbitrary #/m^3
        error = izero
        DO j=1,nion       
           den_save(j)   = profile%en(j)%data(j_ped)
        ENDDO
        den_save(nion+1)   = profile%ene%data(j_ped)
        DO j=1,nion
           profile%en(j)%data(j_ped) = profile%en(j)%data(j_ped)+ delta_ni(j)
           IF(profile%en(j)%data(j_ped) .LT. crit_den)error = j
        ENDDO
        profile%ene%data(j_ped) = profile%ene%data(j_ped) + delta_ne
        IF(profile%ene%data(j_ped) .lt. crit_den)error = nion+1

        IF(error .NE. izero)THEN ! restore input densities
              DO j=1,nion
                 profile%en(j)%data(j_ped) = den_save(j) 
              ENDDO
              profile%ene%data(j_ped)        = den_save(nion+1)
        ENDIF

        RETURN

    END SUBROUTINE mod_profiles



    SUBROUTINE scale_te
!----------------------------------------------------------------------------
!-- Scale input te profile
!----------------------------------------------------------------------------
        USE nrtype,                                   ONLY : DP,I4B,I2B
        USE Plasma_properties ,                      ONLY : profile
        USE ions_gcnmp,                               ONLY : nion

        IMPLICIT none

    RETURN 
    END SUBROUTINE scale_te



    SUBROUTINE scale_ti
!----------------------------------------------------------------------------
!-- Scale input ti profile
!----------------------------------------------------------------------------
        USE nrtype,                                   ONLY : DP,I4B,I2B
        USE Plasma_properties ,                      ONLY : profile
        USE ions_gcnmp,                               ONLY : nion

        IMPLICIT none

    RETURN 
    END SUBROUTINE scale_ti





    SUBROUTINE scale_trot
!----------------------------------------------------------------------------
!-- Scale input toroidal rotation profile
!----------------------------------------------------------------------------
        USE nrtype,                                   ONLY : DP,I4B,I2B
        USE Plasma_properties ,                       ONLY : profile
        USE ions_gcnmp,                               ONLY : nion

        IMPLICIT none

    RETURN 
    END SUBROUTINE scale_trot



    SUBROUTINE get_betan(fi_inc,p_type)
!-----------------------------------------------------------------------------------
! -- GET LOCAL NORMALIZED BETA, mhd_dat%betan%data
!-----------------------------------------------------------------------------------
        USE nrtype,                                   ONLY : DP,I4B,I2B

        USE Plasma_properties ,                       ONLY : profile,dischg,mhd_dat

        USE  vector_class,                            ONLY : zero_vector

        USE  grid_class,                              ONLY : nj,rho_gridn

        USE  common_constants,                        ONLY : Permeability 
 
        IMPLICIT NONE
        INTEGER(I4B) fi_inc,p_type,j

        CALL get_tot_pressure(fi_inc,p_type)

        mhd_dat%betan = zero_Vector(nj)
        DO j =1,nj
           mhd_dat%betan%data(j) = profile%tot_press%data(j)*dischg%rminor       &
                                   *2._DP*Permeability/dischg%btgeom/mhd_dat%tot_cur*1.e8_DP 
                                   ! NOTE 1.e8 = 1.e6(A/MA) x 1.e2(percent)
        ENDDO

        RETURN
 
    END SUBROUTINE get_betan


    SUBROUTINE calc_densities
!--------------------------------------------------------------------------
! -- scale ion densities and get ene,zeff
! -- scaling parameters are from gcnmp namelist:   
! --    nictr = ni_sc(1) central density multiplier
! --    nispline no of knots in interval[0, rho_gridn%data(j_ped_index)]
! --    nirho = ni_sc(4) poiint of inflection
!--------------------------------------------------------------------------
        USE nrtype,                                   ONLY : DP,I4B,I2B

       USE MPI_data,                                  ONLY : myid,master

        USE Plasma_properties ,                       ONLY : profile,dischg,mhd_dat

        USE  vector_class,                            ONLY : zero_vector

        USE  grid_class,                              ONLY : nj,rho_gridn,j_ped_index
   
        USE  ions_gcnmp,                              ONLY : ni_sc,nion,nprim,zeff,z,zsq

        USE  common_constants,                        ONLY : izero,zeroc

        USE  tension_spline,                          ONLY : Smspline

        USE error_handler,                            ONLY : iomaxerr, lerrno,terminate

        USE io_gcnmp,                                 ONLY : nlog

        IMPLICIT NONE

        INTEGER(I4B) j,k,ind,ib,n,fi_inc,p_type
 
        REAL(DP) spline_grid(mhd_dat%ni_spline),drho_lcl,derivr,derivl, &
                 xx,sp,dsp,d2sp
        REAL(DP) ,DIMENSION(mhd_dat%ni_spline)    ::  a,b,c,d,e,g,zm,aa,dd,x,y,weights
        REAL(DP) ,DIMENSION(2*mhd_dat%ni_spline)  :: u,v,w,p, q,r,s,t
        REAL(DP), DIMENSION(nj,nion)              :: ratio

        !save old ratios of ion densities
        DO k=2,nion
           DO j= 1,nj
              ratio(j,k) = profile%en(k)%data(j)/profile%en(1)%data(j)
           ENDDO
        ENDDO

        IF(mhd_dat%ni_spline .LT. 3)THEN
           IF(myid == master)THEN
              PRINT *, 'ERROR:  spline is implementation in calc_densitiesrequires minimum 3 knots'
           ENDIF
           lerrno = iomaxerr + 162
           CALL terminate(lerrno,nlog)
        ENDIF


       ! get index of mhd_dat%rhon_spline that goes with j_ped_index
       ! force first and last knot in density spline to be at  axis(j=1)  and j_ped_index.
       ! intermediate values are taken from input namelist:
         mhd_dat%rhon_spline%data(mhd_dat%ni_spline) = rho_gridn%data(j_ped_index)
         mhd_dat%rhon_spline%data(1)   = zeroc
         y(mhd_dat%ni_spline)          = profile%en(1)%data(j_ped_index)                           ! set pedestal density
         y(1:mhd_dat%ni_spline-1)      = ni_sc(1)*mhd_dat%ni_knot%data(1:mhd_dat%ni_spline-1)           ! set interior densities 
         x(1:mhd_dat%ni_spline)        = mhd_dat%rhon_spline%data(1:mhd_dat%ni_spline)


       ib  = 1                                   ! bc cond (zero derivate at ends)
       derivl = zeroc  ; derivr = zeroc
       ind = 0                                   ! get spline coefficients
       n   = mhd_dat%ni_spline
       weights(:) = zeroc                        ! all zero means interpolating spline instead of smoothing
       CALL Smspline(n,ib,ind,x,y,a,b,c,d,e,g,weights, &
                     derivl,derivr,u,v,w,p,q,r,s,t,aa,dd,zm,xx,sp,dsp,d2sp)



       ind = 1                                   ! evaluate spline
       DO j=1,j_ped_index-1
          xx = rho_gridn%data(j)
          CALL Smspline(n,ib,ind,x,y,a,b,c,d,e,g,weights, &
                        derivl,derivr,u,v,w,p,q,r,s,t,aa,dd,zm,xx,sp,dsp,d2sp)
          profile%en(1)%data(j) = sp
       ENDDO


       DO k=2,nprim                              ! scale remaining primary ion densities
          DO j=1,j_ped_index-1
             profile%en(k)%data(j) = profile%en(1)%data(j)*ratio(j,k)
          ENDDO
       ENDDO
       
! get electron density:
       DO j=1,j_ped_index-1
          profile%ene%data(j) = zeroc
          Do k=1,nion
            profile%ene%data(j) = profile%ene%data(j) + profile%en(k)%data(j)*z(j,k)
          ENDDO
       ENDDO

! get zeff :
       DO j=1,j_ped_index-1
          profile%zeff%data(j) = zeroc
          Do k=1,nion
            profile%zeff%data(j) = profile%ene%data(j) + profile%en(k)%data(j)*zsq(j,k)
          ENDDO
            profile%zeff%data(j) = profile%zeff%data(j)/profile%ene%data(j)
       ENDDO

       fi_inc = 0                          ! =1 ==> include fast ions in tot_pressure
       p_type = 1                          ! use profile%... instead of ene,etc.
       CALL get_betan(fi_inc,p_type)

    END SUBROUTINE calc_densities

