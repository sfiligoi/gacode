#ifdef ONETWO
    SUBROUTINE   set_forcebal_namelist_input(time,e0_b,id_shot,             &
                                              cid_mcsion,cid_scsion,          &
                                              cid_diag,cid_run,cid_mhdeq,     &
                                              cid_device,c_den,k_edotb,mx_ni)
       ! dummy link
    END SUBROUTINE   set_forcebal_namelist_input
    SUBROUTINE set_forcebal_default_switches
       ! dummy link
    END SUBROUTINE set_forcebal_default_switches
#endif

#ifdef GCNMP
     SUBROUTINE   set_forcebal_namelist_input(time,e0_b,id_shot,             &
                                              cid_mcsion,cid_scsion,          &
                                              cid_diag,cid_run,cid_mhdeq,     &
                                              cid_device,c_den,k_edotb,mx_ni)
!-----------------------------------------------------------------------------
! --  time, id_shot,and e0_b are passed back to forcebal main calling routine 
! --  (SUBROUTINE forcebal, forcebal.f90)
!
! takes the place of namelist input for forcebal code.
! ALL input is derived from GCMP/statefile quantities so no namelist read is involved.
!
! Declaration of namelist input variables
!See FORCEBAL_DATA_MOD for declarations of:
! l_banana,            & !option to include banana viscosity [logical]
! l_pfirsch,           & !option to include Pfirsch-Schluter transport [logical]
! l_classical,         & !option to include classical transport [logical]
! l_potato,            & !option to include potato orbit effects [logical]
! l_squeeze,           & !option to include orbit squeezing [logical]
! l_reduce_out           !reduce output for mult charge impurities [logical]
!                        !=.TRUE. only highest charge state
!                        !=.FALSE. all charge states
! k_order                !order of v moments to be solved [integer]
!                        !=2 u and q (default)
!                        !=3 u, q, and u2
!                        !=ebalse error
! nr_r                   !number of radial nodes analysis [integer]
!Following are JET-specific for rerieving data from PPFs:
! cid_uid_efit           !UID name for EFIT solution [character]
! cid_dda_efit           !DDA name for EFIT solution [character]
! cid_uid(25)            !UID names for profiles [character]
! cid_dda(25)            !DDA names for profiles [character]
! cid_dtype(25)          !dtype names for profiles [character]
! id_seqp(25)            !sequence numbers for profiles [integer]
!Following are for fast ions from NBI:
! amu_b(3),            & !atomic mass number of beam ions [-]
! pmw_b(3)               !beam injection power [MW] 
! z_b(3)                 !beam ion charge [-]
!--------------------------------------------------------------------


     USE nrtype,                               ONLY : I4B,DP
     USE FORCEBAL_DATA_MOD,                    ONLY :  l_banana,       & ! forcebal_data_mod.f90
                                                       l_pfirsch,      &
                                                       l_classical,    &
                                                       l_potato,       &
                                                       l_squeeze,      &
                                                       l_reduce_out,   &
                                                       k_order,        &
                                                       nr_r,           &
                                                       cid_uid_efit,   &
                                                       cid_dda_efit,   &
                                                       cid_dtype,      &
                                                       id_seqp,        &
                                                       amu_b,          &
                                                       pmw_b,          &
                                                       z_b
    
     USE grid_class,                           ONLY :  nj

     USE Plasma_properties,                    ONLY : neut_beam

     USE neutral_beams,                        ONLY : nbion,nameb,fdbeam,nbeams

     USE common_constants,                     ONLY : izero,zeroc

     USE solcon_gcnmp,                         ONLY : state_time => time

     USE ions_gcnmp,                           ONLY : namep,namei,nion,     &
                                                      namen,nprim,nimp,nneu

     USE error_handler,                       ONLY : lerrno,terminate,dbg_print

     USE MPI_data,                            ONLY : master,myid

     USE io_gcnmp,                            ONLY : nlog,ncrt

     USE shot_info,                           ONLY : shot_id

     IMPLICIT NONE

     INTEGER(i4B) j,id_shot,mx_ni,k_edotb,k
     REAL(DP) e0_b,e0_avg,pbeam_tot,pbeam_tot2,pbeam_tot3,time,c_den

     CHARACTER(len=3) ::                      &
                        cid_mcsion(mx_ni),    & !names of multiple charge state ions [character] 
                        cid_scsion(mx_ni),    & !names of single charge state ions [character]
                        cid_diag,             & !name of diagnositc ion [character]
                        cid_run                 !run identification [character]

     CHARACTER(len=15) ::                     &
                        cid_mhdeq,            & !name of equilibrium code data used [character]
                                                !='efit' for EFIT equilibria
                                                !=else failure
                        cid_device              !name of device being analyzed [character]


     INTERFACE
        SUBROUTINE to_upper_case(string)
          USE nrtype,            ONLY : I4B
          IMPLICIT NONE
          INTEGER(I4B) l
          CHARACTER*(*), INTENT (INOUT) :: string
        END SUBROUTINE to_upper_case
     END INTERFACE


     CALL set_forcebal_default_namelist(e0_b,id_shot,cid_mcsion,    &
                                        cid_scsion,cid_diag,cid_run,&
                                        cid_mhdeq,cid_device,       &
                                        c_den,k_edotb,mx_ni)


     id_shot = shot_id%shot_nmbr
     nr_r     = nj
     time     = state_time
     k =0
     DO j=1,nion
        IF(j .GT. nprim)THEN
           k=k+1
           IF(namei(k) == 'he')THEN
             cid_mcsion(k) = 'He4'
           ELSEIF(namei(k) == 'c')THEN
             cid_mcsion(k) = 'C  '
           ELSEIF(namei(k) == 'o')THEN
             cid_mcsion(k) = 'O  '
           ELSEIF(namei(k) == 'si')THEN
             cid_mcsion(k) = 'Si '
           ELSEIF(namei(k) == 'ar')THEN
             cid_mcsion(k) = 'Ar '
           ELSEIF(namei(k) == 'ti')THEN
             cid_mcsion(k) = 'Ti '
           ELSEIF(namei(k) == 'cr')THEN
             cid_mcsion(k) = 'Cr '
           ELSEIF(namei(k) == 'fe')THEN
             cid_mcsion(k) = 'Fe '
           ELSEIF(namei(k) == 'ni')THEN
             cid_mcsion(k) = 'Ni  '
           ELSEIF(namei(k) == 'kr')THEN
             cid_mcsion(k) = 'Kr '
           ELSEIF(namei(k) == 'mo')THEN
             cid_mcsion(k) = 'Mo '
           ELSEIF(namei(k) == 'w')THEN
             cid_mcsion(k) = 'W '
           ELSE
              cid_mcsion(k) = 'crp'   ! forces error exit
           ENDIF
!           CALL to_upper_case(cid_mcsion(k))
        ELSE
           IF(namep(j) == 'he')THEN
              cid_scsion(j) = 'He4'
           ELSE
              cid_scsion(j) = namep(j)
              CALL to_upper_case(cid_scsion(j))
           ENDIF
        ENDIF
     ENDDO

     cid_mhdeq ='efit'
     cid_device='scrt'
     pbeam_tot = zeroc ;pbeam_tot2 = zeroc ; pbeam_tot3 = zeroc ; e0_avg = zeroc
     ! NOTE: neut_beam%pbeam(1:3,1:nbeams), full half third power for each beamlet

     IF(nbion .NE. izero .AND. neut_beam%nbeams .NE. izero)THEN
        DO j = 1, neut_beam%nbeams  ! typically 14 beamlests for DIII-D
           ! density weighted beam energy since we need single beam: 
           e0_avg     = e0_avg + neut_beam%pbeam(1,j)*neut_beam%ebeam(1,j)
           pbeam_tot  = pbeam_tot + neut_beam%pbeam(1,j)
           pbeam_tot2 = pbeam_tot2 + neut_beam%pbeam(2,j)
           pbeam_tot3 = pbeam_tot3 + neut_beam%pbeam(3,j)
        ENDDO

        e0_b   = e0_avg/pbeam_tot

        DO j=1,nbion
           IF( nameb(j) == 'd' ) THEN
              amu_b(j) = 2._DP
           ELSE
              amu_b(j) = 1._DP
           ENDIF
        ENDDO

        pmw_b(1) =   pbeam_tot/1.e6_DP   ! MW
        pmw_b(2) =   pbeam_tot2/1.e6_DP
        pmw_b(3) =   pbeam_tot3/1.e6_DP
        z_b(:)   =   1._DP

     ELSE
        amu_b(:)  =   zeroc
        e0_b     =   zeroc            ! full beam energy (half third taken in forcebal(
        pmw_b(:) =   zeroc
        z_b(:)   =   zeroc
     ENDIF
        IF(dbg_print .AND. myid == master)THEN
          WRITE(nlog,10)nbion,nbeams,neut_beam%nbeams,e0_b
          WRITE(ncrt,10)nbion,nbeams,neut_beam%nbeams,e0_b

10        FORMAT("nbion,nbeams,neut_beam%nbeams,e0_b in set_forcebal_namelist_input",i3,x,i3,x,i3,x,1pe12.4) 
!          lerrno = 100
!          call terminate(lerrno,nlog)
        ENDIF
      RETURN
 
     END SUBROUTINE  set_forcebal_namelist_input


     SUBROUTINE  set_forcebal_default_namelist(e0_b,id_shot,cid_mcsion,     &
                                               cid_scsion,cid_diag,cid_run, &
                                               cid_mhdeq,cid_device,c_den,  &
                                               k_edotb,mx_ni)
!-------------------------------------------------------------------------------
! -- set the defaults for forcebal namelist
! -- Note that set_forcebal_switches is called from gcnmp_namelist defaults
! -- hence frocebal l_* switches are set at this point.
!---------------------------------------------------------------------HSJ-------  

     USE nrtype,                               ONLY : I4B,DP
     USE SPEC_KIND_MOD
     USE FORCEBAL_DATA_MOD,                    ONLY :  l_banana,       & ! forcebal_data_mod.f90
                                                       l_pfirsch,      &
                                                       l_classical,    &
                                                       l_potato,       &
                                                       l_squeeze,      &
                                                       l_reduce_out,   &
                                                       k_order,        &
                                                       nr_r,           &
                                                       cid_uid_efit,   &
                                                       cid_dda_efit,   &
                                                       cid_uid,        &
                                                       cid_dda ,       &
                                                       cid_dtype ,     &
                                                       id_seqp,        &
                                                       amu_b,          &
                                                       pmw_b,          &
                                                       z_b    

      USE common_constants,                     ONLY : izero,zeroc
 
      IMPLICIT NONE

      REAL(DP) e0_b,time,c_den

      INTEGER(I4B) id_shot,mx_ni,k_edotb

      CHARACTER(len=3) :: &
                          cid_mcsion(mx_ni),     & !names of multiple charge state ions 
                          cid_scsion(mx_ni),     & !names of single charge state ions 
                          cid_diag,              & !name of diagnositc ion 

                          cid_run                  !run identification 

      CHARACTER(len=15) ::                       &
                          cid_mhdeq,             &  !name of equilibrium code data used 
                                                    !='efit' for EFIT equilibria
                                                    !=else failure
                          cid_device                !name of device being analyzed 


     !Namelist input - character
     cid_device='diiid'
     cid_mhdeq='efit'
     cid_run='A01'
     cid_diag=''
     cid_scsion(:)=''
     cid_mcsion(:)=''
     cid_uid_efit=''
     cid_dda_efit=''
     cid_uid(:)=''
     cid_dda(:)=''
     cid_dtype(:)=''

     !Namelist input - integer
     k_edotb=0                      ! = 1 means  use calc parallel field
                                    ! read_pro ( file read_pro_mod.f90)
                                    ! sub forcebal_nclass will use 
                                    ! total ,beam and bootstrap and resistivity
                                    ! to get parallel electric field.
                                    ! =0 means e_par_ex_r(1:nj)  will be 
                                    ! used directly ( not set up for this)

!     k_order=2        set in forcebal_switches
     id_shot=10000
     nr_r=41
     id_seqp(:)=0

     !Namelist input - real
     time     = zeroc
     amu_b(:) = zeroc
     e0_b     = zeroc
     pmw_b(:) = zeroc
     z_b(:)   = zeroc

     !Other
     c_den=1.0e10_rspec       !density cutoff below which species is ignored [/m**3]




     RETURN

     END SUBROUTINE  set_forcebal_default_namelist


     SUBROUTINE set_forcebal_default_switches
     !----------------------------------------------------------------
     ! -- this subroutine  sets those items in forcebal namelist that
     ! -- are actually set in gcnmp namelist
     !-----------------------------------------------------------HSJ--
       USE FORCEBAL_DATA_MOD,                ONLY :  l_banana,       & ! forcebal_data_mod.f90
                                                     l_pfirsch,      &
                                                     l_classical,    &
                                                     l_potato,       &
                                                     l_squeeze,      &
                                                     l_reduce_out,   &
                                                     k_order
   
       USE solcon_gcnmp,                     ONLY :  use_forcebal   

       !Namelist input 
       l_banana=.TRUE.
       l_pfirsch=.TRUE.
       l_classical=.TRUE.
       l_potato=.FALSE.
       l_squeeze=.FALSE.
       l_reduce_out=.TRUE.
       k_order = 2
       use_forcebal = 0

     RETURN

     END SUBROUTINE set_forcebal_default_switches



     SUBROUTINE dump_forcebal_namelist(time,n_sum,e0_b,id_shot,       &
                                       cid_mcsion,                    &
                                       cid_scsion,cid_diag,cid_run,   &
                                       cid_mhdeq,cid_device,          &
                                       k_edotb,c_den,mx_ni)

     USE nrtype,                               ONLY : I4B,DP,I2B
     USE FORCEBAL_DATA_MOD,                    ONLY :  l_banana,       & !forcebal_data_mod.f90
                                                       l_pfirsch,      &
                                                       l_classical,    &
                                                       l_potato,       &
                                                       l_squeeze,      &
                                                       l_reduce_out,   &
                                                       k_order,        &
                                                       nr_r,           &
                                                       cid_uid_efit,   &
                                                       cid_uid,        &
                                                       cid_dda,        &
                                                       cid_dtype,      &
                                                       cid_dda_efit,   &
                                                       id_seqp,        &
                                                       amu_b,          &
                                                       pmw_b,          &
                                                       z_b 
    
     USE MPI_data,                             ONLY :  myid,master
 
     IMPLICIT NONE

     REAL(DP) e0_b,c_den,time

     INTEGER(I4B)id_shot,mx_ni,k_edotb,mx_min
     INTEGER(I2B)n_sum

     CHARACTER(len=3) :: &
          cid_run,             & !run identification 
          cid_mcsion(mx_ni),   & !names of multiple charge state ions 
          cid_scsion(mx_ni),   & !names of single charge state ions 
          cid_diag               !name of diagnositc ion 

     CHARACTER(len=3) :: &
          cidc_mcsion(6),  & !names of multiple charge state ions 
          cidc_scsion(6)

     CHARACTER(len=15) ::      &
          cid_mhdeq,           &  !name of equilibrium code data used 
                                  !='efit' for EFIT equilibria
                                  !=else failure
          cid_device              !name of device being analyzed 



     NAMELIST /inforce/l_banana,l_pfirsch,l_classical,l_potato,l_squeeze, &
                      l_reduce_out, &
                      cid_device,cid_mhdeq,cid_run,cid_diag,& 
                      cidc_scsion,cidc_mcsion, & 
                      k_edotb,k_order, &
                      id_shot,nr_r, &
                      time, &
                      cid_uid_efit,cid_dda_efit,         &
                      cid_uid,cid_dda,cid_dtype,id_seqp, &
                      amu_b,z_b,e0_b,pmw_b,c_den
 

          mx_min = MIN(6,mx_ni)
         cidc_scsion(1:mx_min) = cid_scsion(1:mx_min)    ! cant dump variable length arrays in namelist
         cidc_mcsion(1:mx_min) = cid_mcsion(1:mx_min)

          IF(myid == master) WRITE(n_sum,inforce)

     
     RETURN

   END SUBROUTINE dump_forcebal_namelist
#endif

