 SUBROUTINE iter_dbase_nc_compat
  ! ----------------------------------------------------------------------
  ! this is a special version of  the netcdf file hat is identical to the
  ! iter_dbase_nc subrotien except for some name chnges that DIII-d uses.
  ! This routine  is called if ifmt = 1
  !CCA  NOTE: this is a common component architecture routine  shared by 
  !(currently) Onetwo and gcnmp. Changes must be consistent with all
  !associated codes .  !!!!!!!!!!!!!!!!!!!!! 
  !This almost certainly menas that in order to use this component you will 
  !have to create and interface . (See set_12_gcnmp_vars.f90 for an
  !example).
  !
  !
  !
  ! --- collects and writes out (or reads in) data for ITER database
  ! --- netcdf version. The netcdf input has different units,names, etc.
  ! --- from the Onetwo variables. Hence we use the modules from  gcnmp, 
  ! --- read in and isolate these quantities and then move them
  ! --- appropriately to Onetwo variables (in set_onetwo_vars)
  ! --- This routine is considerably extended over the older iterdb files
  ! --- and consequently may not work with codes that expect the older versions.
  ! --- On the other hand if the netcdf implementation for reading the 
  ! --- iterdb file is done  correctly in these other codes  then the
  ! ---  additions presented here will 
  ! --- not cause problems since netcdf will ignore variables that are not
  ! --- defined. 
  ! -------------------------------------------------------HSJ-8/15/07------
  ! 
  ! --- INPUT/OUTPUT
  !

  !     irwflag          read/write flag, included to make the reading
  !                      of the file created by this subroutine
  !                      consistent with the writing of the file.
  !                      irwflag = 0  ==>  WRITE the data
  !                      irwflag = 1  ==>  READ  the data
  !
  !     iterdb_file_name   name of file which will be read 
  !                      (if called in read mode)
  !     iterdb_outpt     name of output file (created in this routine,
  !                      iterdb_outpt = iterdb_file_name (could add //time)
  !     iterdsc          switch used to output description of parameters
  !
  !
  !     nplasbdry
  !     rplasbdry(j)
  !     zplasbdry(j)     j=1,2..nplasbdry describes the plasma boundary
  !
  !     q(j)             j=1,2,...nj, the safety factor on rho grid
  !     betap
  !     beta
  !
  !     tocur            total plasma current, amps
  !

  !     rho(j)           j=1,2,..npsi the rho values corresponding
  !                      to the psi grid in meters
  !     qpsi(j)          the corresponding q values
  !     psivolpnpsi(j)       volume, meter**3, enclosed by flux surface
  !     grho1(j)         < ABS (grad rho)>
  !     grho2(j)         <(grad rho)**2>
  !     rmajavnpsi(j)    j=1,2,..npsi average major radius
  !                      at elevation of magnetic axis
  !     rminavnpsi(j)    same for minor radius
  !     psivalnpsi(j)        j=1,2..npsi
  !
  !     triangnpsi_u(j)    (upper)triangularity
  !     triangnpsi_l(j)  (lower)  triangularity
  !     triangnpsi_lcl   array created here for compatibility with
  !                      older versions. it is the aveerage of upper
  !                      and lower tirangularity on nj rho grid.
  !
  !     pindentnpsi(j)   indentation
  !     sfareanpsi(j)    surface area of flux surfaces
  !     cxareanpsi(j)    cross-sectional area
  !
  !     btor             mag field at R0, in vacuum, tesla
  !     elong(j)         elongation
  !
  !
  !     fcap(j)
  !     gcap(j)
  !     hcap(j)          the geometric factors
  !     ali              plasma inductance
  !

  !     namep(i)         i=1,..nprim
  !     fd_thermal       if primary species is 'dt' then fd_thermal gives
  !                      fraction of d in 'dt' mixture. If primary species 
  !                      is/is not 'd'  then set fd_thermal 1/0
  !     namei(i)         i=1,..nimp
  !     namen(i)         i=1,..nneu 
  !     zeff(j)          j=1,2...nj
  !     nameb(i)         i = 1,nbion name of beam ions ('dt', 'd','t','h','he')
  !     nbion            number of fast beam ion  species. For a 'dt'
  !                      beam mixture use nbion =1. The fraction of
  !     fd_beam          d in the 'dt' mixture is fd_beam. For single species
  !                      beam fd_beam set fdbeam = 1.0
  !     sbeame            electron source form beam
  !     sbeam(j,i)       thermal ion source from beam, grid point j, species nameb(i)
  !
  !     enbeam(j,i)      j=1,2..nj,i=1..nbion  fast ion density due to beam,
  !                      grid point j, species nameb(i).
  !     wbeam(j)         beam stored energy density
  !     enalp(J)         alpha particel den sity
  !     walp(j)          alpsh particle stored energy density
  !
  !     enn(j,k)         j=1,2..nj(njs),k=1,2 neutral density
  !     ennw             neutral density due to wall source
  !     ennv                                    volume
  !     volsn            volume source of neutrals
  ! 
  !     rmajor 
  !     kappa
  !     btor
  !
  !     r(j)             j=1,2...nj,rho grid
  !
  !     time            current time, sec
  !     tGCNMf        the transport equations are evolved from time
  !                     to tGCNMf (sec).
  !                     The interval (time,tGCNMf] will be broken into
  !                     smaller intervals by the solver as appropriate.
  !                     results will be returned in the form of a new
  !                     iterdb file at tGCNMf. It is  assumed that
  !                     all particle and  energy sources are constant
  !                     in this time interval. Metrics introduced by
  !                     flux surface averaging (eq *cap quantites)
  !                     are also assumed constant in time over 
  !                     this interval
  ! 
  !     time_bc         boundary condition time. This time must be 
  !                     >= tGCNMf!!
  !                     It is assumed that the
  !                     interval (time,time_bc) is small enough that
  !                     linear inerpolation of boundary conditions is
  !                     reasonble. 
  !                     Initial conditions required
  !                     at t = time are taken from en(j,k),ene,zeff,te,ti,
  !                     curden,angrot (see below)  as appropriate. 
  !                     Boundary conditions
  !                     required to evolve from time to time_bc are
  !                     obtained by linear interpolation  between these 
  !                     profiles and the values given in en_bc(j,k)
  !                     ene_bc,zeff_bc,te_bc,ti_bc,angrot_bc
  !                     (see below).  Although only
  !                     information for the plasma edge is required for
  !                     the boundary conditions it is assumed that
  !                     the *_bc profiles are given over the entire rho grid
  !                     This allows us to apply boundary conditions 
  !                     for each dependent variable at a different rho
  !                     value, rather than limiting the boundary to just 
  !                     rho= 1.0. The switches  
  !                         fix-edge_te_bc  (= any value in (0.0,1.0] )     
  !                         fix_edge_ti_bc  (= any value in (0.0,1.0] )
  !                         fix_edge_rot_bc (= any value in (0.0,1.0] )
  !                         fix_edge_ni_bc  (= any value in (0.0,1.0] )
  !                         interact with the *_bc profiles in that
  !                         for rho > fix_edge_** the value of the *_bc
  !                         profile is used. To be specifie consider
  !                         for example the te profile: 
  !                         Suppose  fix_eddge_te =0.91 . Lets assume 
  !                         that the closest
  !                         normalized rho grid point is at rho =0.925.
  !                         then the boundary condtion for te at t= time is
  !                         takens as te(rho=0.925) . At t= time_bc the
  !                         bc is taken as te_bc(rho = 0.925). For arbitrary
  !                         times in the interval [time,time_bc] the bc for
  !                         te is obtained by lineat interpolation in time 
  !                         at the grid point rho = 0.925. (Because the rho grid
  !                         may change in time normalized rho is used.)
  !                         For all rho grid points > 0.925 the te profile
  !                         is also obtained by this linear inerpolation process
  !                         But note that the temperature at these points
  !                         is given a priori by the profiles te and te_bc.
  !                         The te profile for rho < 0.925 is determined 
  !                         by solving the PDE's. Note that we always solve
  !                         Faraday's law out to the plasma edge( rho =1.0) and
  !                         rhis requires te, etc all the way out to the edge.
  !                         Hence zpecification of te_bc in the edge zone will 
  !                         still influence the solution. Finally if te is run
  !                         in analysis mode then the linear inerpolation in 
  !                         time discussed above is done for all rho values from
  !                         0.0 to 1.0  inclusive. (The above holds for all
  !                         profiles, not just te).
  !                      



  !                      
  !     eqtime          time of last equilibirum 
  !     nj              transport grid size
  !     nprim
  !     nneu
  !     nimp
  !     nion
  !
  !     psir_grid(j)         j=1,2...nj psi on r (i.e., rho) grid
  !
  !     profiles at t = time (used as inital conditions)
  !     te(j)           electron and ion temperatures,keV
  !     ti(j)
  !     ene(j)          electron density
  !     en(j,k)         j=1,2..nj,k=1..nion,where nion=nprim+nimp is the
  !                     number of primary plus impurity ion species
  !     rbp(j)          rho*bp field
  !     profiles at t=time_bc (used to determine boundary condtions)
  !         te_bc
  !         ti_bc
  !         ene_bc
  !         en(j,k)_bc
  !         zeff_bc!
  !         angrot_bc

  !     
  !     fixe_edge_**  switches allow  solving the problem from
  !     rho =0 to rho = fixedge_**  .le. 1.0
  !     for fix_edge < rho <=1.0 the code takes the input profile values.  
  !     these values are cosntant during a call to gcnm but could
  !     change from one call to the next,as could rho,etcc so this
  !     is a fully dynamicallay controlled model.     
  !     fix-edge_te_bc      	    
  !     fix_edge_ti_bc           
  !     fix_edge_rot_bc 
  !     fix_edge_ni_bc        





  !     currf(j)
  !     curbe(i)
  !     curbi(j)
  !     qrad(j)         radiated power
  !     totohm
  !     totbeam
  !     totboot
  !     totrf           total currents, amps
  !     qconde,qcondi
  !     qconve,qconvi
  !     qdelt
  !     qexch    not currently used
  !     qione,qioni
  !     qbeame,qbeami
  !     qrfe, qrfi
  !     qe2d, qi2d
  !     qtfuse,qtfusi
  !     qmag
  !     qsawe,qsawi
  !     qbfusi
  !     sbeam
  !     dnedt(j)         change in electron density at this time 1/(M**3sec)
  !     dangrotdt(J)     change in toroidal rotation
  !     cheinv(j)
  !     chiinv(j)        electron and ion thermal cond.
  !     xkineo(j)        ion neoclassical thermal conductivity
  !     xkeneo(j)        electron neoclassical thermal conductivity
  !
  !     angrot(j)       angular rotation speed, radian/second
  !     storqueb(j)     beam torque density nt-m/m**3 
  !
  !     eps(j)          horizontal inverse aspect ratio = (rmax-rmin)/(rmax+rmin)
  !                     on each flux surface. Add to iterdb file
  !     xhm2
  !     xi11
  !     xi33
  !     xips 
  ! -------------------------------------------------------- HSJ
  !
  USE typeSizes                  !part of netcdf 
  USE netcdf                     !ditto

  !
  USE nrtype,                    ONLY : DP,I4B,I2B


  USE Plasma_properties ,        ONLY : dischg,profile,mhd_dat,get_values,   &
                                        diffuse,pwrden, wpdot,prtcl_src,     &
                                        pellet

  USE shot_info,                 ONLY : shot_id
  USE vector_class,              ONLY : new_Vector,get_element,              &
       real_mult_Vector,list
  USE grid_class,                ONLY : nj,psir_grid,rho_grid,rho_gridn,     &
                                        eps,xhm2,xi11,xi33,xips,xhm20,xi110,   &
                                        xi330,xips0,rho_mhd_gridnpsi

  USE error_handler,             ONLY : lerrno,terminate

#ifdef GCNMP
  USE io_gcnmp,                  ONLY : ncrt,nlog
  USE gcnmp_version,             ONLY : gcnmp_ver
#else
  USE io,                        ONLY : ncrt,versid,nlog => nitre
#endif


  USE source_terms_gcnmp,        ONLY : stsource,scx,sion,srecom,sbcx,sbeame, &
                                        dudtsv, sbeam

  USE iterdbmd_gcnmp,            ONLY : iterdsc,irwflag,iterdb_file_name,      &
                                        iterdb_outpt,ncid

  USE ions_gcnmp,                ONLY : namep,namei,nion,                     &
                                        namen,nprim,nimp,nneu,z,zsq,zeff,     &
                                        name_size,fd_thermal

  USE MPI_data,                  ONLY : mpiierr,myid



  USE solcon_gcnmp,              ONLY : time,eqtime,tGCNMf,tGCNMs



  USE neutral_data,              ONLY : enn,ennw,volsn,ennv

  USE fast_ion_data_gcnmp,       ONLY : walp,enalp

  USE neutral_beams,             ONLY : enbeam,storqueb,enbeam_tot,wbeam,   &
                                        fd_beam,nbion,nameb


  USE curden_terms,              ONLY : currf,curbeam,ibcur,irfc

  USE common_constants,          ONLY : convert,zeroc,izero

  USE bc_values_gcnmp,           ONLY : totcur_bc,time_bc,                 &
                                        fix_edge_te_bc,fix_edge_ti_bc,     &
                                        fix_edge_rot_bc,fix_edge_ni_bc,    &
                                        zeff_bc,ene_bc,te_bc,ti_bc,        &
                                        angrot_bc,en_bc,vloop_bc,flux_bc,  &
                                        convert_edge_indecies

  USE dep_var,                   ONLY : dp4

  USE fdyimp,                    ONLY : dfdt,dgdt,dhdt

  USE tglfin,                    ONLY : tglf_p_output,tglf_e_output

  USE cer,                       ONLY : kpol_c,kpol_d,kpol_exp,angrot_c,    &
                                        angrot_d,ave_vpar_d,udia_d,udia_c,  &
                                        ugrt,sqz_d,epsi,epsi_exp,cer_btdr,  &
                                        cer_bp

  IMPLICIT  NONE


  REAL(DP), ALLOCATABLE,DIMENSION(:)   :: work_nj,               &
       work_npsi,             &
       work_nplasbdry,        &
       work_nlimiter,         &
       work_nr,               &
       work_nz
  REAL(DP), ALLOCATABLE,DIMENSION(:,:) :: work_nj_nion,work_nj_ntot
!
!     F90 automatic (temporary) arrays:
  REAL(DP)             aspline(mhd_dat%npsi), bspline(mhd_dat%npsi),    &
                       cspline(mhd_dat%npsi),cs2spline(mhd_dat%npsi,3), &
                       dspline(mhd_dat%npsi), espline(mhd_dat%npsi),    &
                       fspline(mhd_dat%npsi), bpar(4), work(nj),        &
                       eni(mhd_dat%npsi,nimp),triangnpsi_lcl(mhd_dat%npsi)

  REAL(DP) elmt,tmax,tension
  INTEGER(I4B)  j,jp,jj,ji,jn,ntot,titlen
  INTEGER(I2B) i
  INTEGER(I4B), PARAMETER :: one = 1_I4B
  INTEGER(I4B), PARAMETER :: nf_max_name = 128 !should be nf90_max_name but cant find
  !any definition of this in netcdf modules
  !(v3.6)
  CHARACTER(len = nf_max_name) :: gen_name=' '
  
  INTEGER(I4B)           &         !netcdf required values
       rcode,            &
       id_title,         &
       id_shot,          &
       id_nj,            &
       id_time,          &
       id_rgeom,         &
       id_btgeom,        &
       id_rmag,          &
       id_zma,           &
       id_rsep,          &
       id_zsep,          &
       id_rmajor,        &
       id_rplasmin,      &
       id_rplasmax,      &
       id_zplasmin,      &
       id_zplasmax,      &
       id_kappa,         &
       id_psiaxis,       &
       id_psibdry,       &
       id_deltao,        &
       id_pindento,      &
       id_volume,        &
       id_circum,        &
       id_areao,         &
       id_nion,          &
       id_nprim,         &
       id_nneu,          &
       id_enn,           &
       id_ennw,          &
       id_ennv,          &
       id_nbion ,        &
       id_namep,         &
       id_namen,         &
       id_namepel,       &
       id_tgcnmf,        &
       id_time_bc,       &
       id_fd_thermal,    &
       id_nimp,          &
       id_fd_beam,       &
       id_namei,         &
       id_nameb,         &
       id_fpsinpsi,      &
       id_psivalnpsi,    &
       id_ravgnpsi,      &
       id_ravginpsi,     &
       id_ffprim,        &
       id_pprim,         &
       id_btor,          &
       id_bp,            &
       id_bprmaj,        &
       id_btotrmaj,      &
       id_tot_cur,       &
       id_totohm_cur,    &
       id_totboot_cur,   &
       id_totbeam_cur,   &
       id_totrf_cur,     &
       id_betap,         &
       id_beta,          &
       id_ali,           &
       id_te0,           &
       id_ti0,           &
       id_psi,           &
       id_psir_grid,     &
       id_rho_grid,      &
       id_rho_mhd_gridnpsi,  &
       id_rmhdgrid,      &
       id_zmhdgrid,      &
       id_fcap,          &
       id_gcap,          &
       id_hcap,          &
       id_te,            &
       id_ti,            &
       id_q_value,       & 
       id_ene,           &
       id_enp,           &
       id_eni,           &
       id_pflux,         &
       id_efluxe,        &
       id_efluxi,        &
       id_fdyflux,       &
       id_rotflux,       &
       id_tglf_p_fluxe,  &
       id_tglf_e_fluxe,  &
       id_tglf_p_fluxp,  &
       id_tglf_e_fluxp,  &
       id_tglf_p_fluxi,  &
       id_tglf_e_fluxi,  &
       id_sion,          &
       id_srecom,        &
       id_scx,           &
       id_sbcx,          &
       id_stsource,      &
       id_dudtsv,        &
       id_volsn,         &
       id_stfuse,        &
       id_sbfuse,        &
       id_sbeame,        &
       id_sbeam,         &
       id_spellet,       &
       id_curden,        &
       id_curohm,        &
       id_curboot,       &
       id_curbeam,       &
       id_currf,         &
       id_etor,          &
       id_rbp,           &
       id_zeff,          &
       id_angrot,        &
       id_d,             &
       id_chieinv,       &
       id_chiinve,       &
       id_xkeneo,        &
       id_dpedt,         &
       id_dpidt,         &
       id_qconde,        &
       id_qcondi,        &
       id_qconve,        &
       id_qconvi,        &
       id_qbeame,        &
       id_qbeami,        &
       id_qdelt,         &
       id_qexch,         &
       id_qrfe,          &
       id_qrfi,          &
       id_qione,         &
       id_qioni,         &
       id_qcx,           &
       id_qe2d,          &
       id_qi2d,          &
       id_qfuse,         &
       id_qfusi,         &
       id_qbfue,         &
       id_qbfusi,        &
       id_qmag,          &
       id_qsawe,         &
       id_qsawi,         &
       id_qrad,          &
       id_qohm,          &
       id_rmajavnpsi,    &
       id_chiinv,        &
       id_xkineo,        &
       id_qbfuse,        &
       id_qbfusei,       &
       id_rminavnpsi,    &
       id_psivolpnpsi,   &
       id_psivolp,       &
       id_elongxnpsi,    &
       id_elongx,        &
       id_triangnpsi_u,  &
       id_triangnpsi_l,  &
       id_triangnpsi_lcl,&
       id_pindentnpsi,   &
       id_sfareanpsi,    &
       id_cxareanpsi,    &
       id_grho1npsi,     &
       id_grho2npsi,     &
       id_nplasbdry,     &
       id_rplasbdry,     &
       id_zplasbdry,     &
       id_rlimiter,      &
       id_zlimiter,      &
       id_storqueb,      &
       id_totcur_bc,     &
       id_vloop_bc,      &
       id_fix_edge_te_bc,   &
       id_fix_edge_ti_bc,   &
       id_fix_edge_rot_bc,  &
       id_fix_edge_ni_bc,   &
       id_te_bc,            &
       id_ti_bc,            &
       id_press,            &
       id_pressb,           &
       id_ene_bc,           &
       id_zeff_bc,          &
       id_angrot_bc,        &
       id_en_bc,            &
       id_flux_bc,          &
       id_z,                &
       id_zsq,              &
       id_wbeam,            &
       id_walp,             &
       id_enalp,            &
       id_dnedt,            &
       id_eps,              &
       id_rcap,             &
       id_rcapi,            &
       id_r2capi,           &
       id_r2cap,            &
       id_xi11,             &
       id_xi33,             &
       id_xips,             &
       id_xhm2,             &
       id_enbeam,           &
       id_dfdt,             &
       id_dgdt,             &
       id_dhdt,             &
       id_kpol_c,           &
       id_kpol_d,           &
       id_kpol_exp,         &
       id_angrot_c,         &
       id_angrot_d,         &
       id_ave_vpar_d,       &
       id_udia_d,           &
       id_udia_c,           &
       id_ugrt,             &
       id_sqz_d,            &
       id_epsi,             &
       id_epsi_exp,         &
       id_cer_btdr,         &
       id_cer_bp,           &
       idim_char,           &
       idim_rho,            &
       idim_nr_mhd,         &
       idim_nz_mhd,         &
       idim_nion,           &
       idim_nprim,          &
       idim_nimp,           &
       idim_nneu,           &
       idim_nbion,          &
       idim_nplasbdry,      &
       idim_ntot,           &
       idim_npsi,           &
       idim_nlimiter,       &
       nDimensions,         &
       nVariables,          &
       nAttributes,         &
       unlimitedDimId,      &
       formatNum,           &
       idlen,               &
       index,               &
       kpsi,                &
       ier

  INTEGER(I4B),PARAMETER :: max_dims = 13
  INTEGER(I4B),PARAMETER :: dim_size = 16
  CHARACTER(LEN=dim_size) dim_names(max_dims)
  LOGICAL    monotonic
  INTEGER(I4B),PARAMETER :: strlen = 132
  CHARACTER  starflag*2, headerline*132,line*132,time_str*36
  CHARACTER  (LEN= strlen) label ,tlabel,base_label
  CHARACTER   bc_asc_time*24,st_asc_time*24 
  CHARACTER(len = name_size) :: tname
  CHARACTER(LEN = strlen) :: title

  DATA dim_names /"dim_rho","dim_ion","dim_nprim","dim_char",                &
       "dim_nimp", "dim_neu","dim_fi","dim_nplasbdry",            &
       "dim_ntot","dim_nr_mhd","dim_nz_mhd","dim_npsi",           &
       "dim_nlimiter"/
  INTERFACE 
     SUBROUTINE check_monotonic (array, n, monotonic, incr) 
       USE nrtype,   ONLY : I4B,DP
       INTEGER(I4B), INTENT(IN) :: n,incr
       REAL(DP),INTENT(IN),DIMENSION(:) :: array(n)
       LOGICAL, INTENT(OUT) ::    monotonic
     END SUBROUTINE check_monotonic

     SUBROUTINE netcdf_err(status,flag)
       USE nrtype,   ONLY : I4B,DP
       IMPLICIT NONE
       INTEGER, INTENT ( in) :: status
       INTEGER(I4B), OPTIONAL, INTENT(in)     :: flag
     END SUBROUTINE netcdf_err


  END INTERFACE

  iterdsc = 1                   ! always write descriptor,reads will fail otherwise

7 FORMAT (         a   )
8 FORMAT (5(2x,    a  ))        ! common character write/read format
9 FORMAT (5(2x,   i6  ))        ! common integer   write/read format
10 FORMAT (5(2x,1pe14.4))        ! common floating  write/read format
11 FORMAT (2x,i6,2x,1pe14.4)

  !     IF(ALLOCATED(work))DEALLOCATE(work)
  !     ALLOCATE (work(MAX(nj,mhd_dat%npsi,dischg%nplasbdry,dischg%nlimiter))) 
  !     universal array wont work, netcdf checks sizes internally:
  IF(ALLOCATED(work_nj))DEALLOCATE(work_nj)
  IF(ALLOCATED(work_npsi))DEALLOCATE(work_npsi)
  IF(ALLOCATED(work_nplasbdry))DEALLOCATE(work_nplasbdry)
  IF(ALLOCATED(work_nlimiter))DEALLOCATE(work_nlimiter)
  IF(ALLOCATED(work_nr))DEALLOCATE(work_nr)
  IF(ALLOCATED(work_nz))DEALLOCATE(work_nz)
  IF( ALLOCATED(work_nj_nion))DEALLOCATE(work_nj_nion)






  !
  IF (irwflag .EQ. 0) THEN      ! WRITE to netcdf file 
     iterdb_outpt = iterdb_file_name(1:LEN_TRIM(iterdb_file_name))
     rcode  = NF90_CREATE ( iterdb_outpt,NF90_CLOBBER, ncid) !enters define mode
     IF (rcode .NE.  nf90_noerr )THEN
        lerrno = 3_I4B
        CALL  terminate(lerrno,nlog)
     ENDIF
     WRITE(label,FMT='(1pe16.8)')time_bc          !get  bc time in string form for output
     READ(label,FMT='(a)')bc_asc_time
     bc_asc_time = ADJUSTL(bc_asc_time)
  ELSE                     ! READ from existing file 

     rcode =  NF90_OPEN (iterdb_file_name,NF90_NOWRITE, ncid)
     IF (rcode == nf90_noerr ) go to 3
2    WRITE  (ncrt, 4)  iterdb_file_name(1:LEN_TRIM(iterdb_file_name))
     WRITE  (nlog, 4)  iterdb_file_name(1:LEN_TRIM(iterdb_file_name))
4    FORMAT (                                                     / &
          ' ERROR: subroutine ITER_DBASE_NC has encountered an error' / &
          '        the ITER database file "', a, '" cannot be opened')
     lerrno = 3_I4B
     CALL  terminate(lerrno,nlog)
3    CONTINUE

  END IF



  !   ------------------------------------------------------------------

  !     Netcdf file is ready for reading or writing:

  !   ------------------------------------------------------------------

  rwncd:  IF (irwflag .EQ. 0) THEN       ! write te data

     ALLOCATE(work_nj(nj))
     ALLOCATE(work_npsi(mhd_dat%npsi))
     ALLOCATE(work_nplasbdry(dischg%nplasbdry))
     ALLOCATE(work_nlimiter(dischg%nlimiter))
     ALLOCATE(work_nr(dischg%nr_mhd))
     ALLOCATE(work_nz(dischg%nz_mhd))
     ALLOCATE(work_nj_nion(nj,nion))
     !     in the future we may introduce time dependent metrics.
     !     then dfdt (= d/dt FCAP ,etc) will be calculated internally in this code
     !     by reading in fcap_bc, etc.
     !     in demo code just assume time independent for each
     !     slice that the solver is called.
     IF(.NOT. ALLOCATED(dfdt))THEN 
        ALLOCATE(dfdt(nj))
        dfdt(:) = zeroc
     ENDIF
     IF(.NOT. ALLOCATED(dgdt))THEN 
        ALLOCATE(dgdt(nj))
        dgdt(:) = zeroc
     ENDIF
     IF(.NOT. ALLOCATED(dhdt))THEN
        ALLOCATE(dhdt(nj))
        dhdt(:) = zeroc
     ENDIF
     ! title of dataset:
#ifdef GCNMP
     title = 'GCNMP  netCDF output file '//gcnmp_ver(1:LEN_TRIM(gcnmp_ver))
#else
     title = 'ONETWO  netCDF state file '//versid
#endif
     CALL netcdf_err( nf90_put_att(ncid,  NF90_GLOBAL,"title",title))

     ! define dimensions:

     CALL netcdf_err( nf90_def_dim(ncid, dim_names(1), nj, idim_rho),idim_rho)                      !"dim_rho" 

     CALL netcdf_err( nf90_def_dim(ncid, dim_names(2), nion, idim_nion),idim_nion)                   !"dim_ion"

     CALL netcdf_err( nf90_def_dim(ncid, dim_names(3), nprim, idim_nprim),idim_nprim)                 !"dim_nprim"

     CALL netcdf_err( nf90_def_dim(ncid, dim_names(4), name_size, idim_char),idim_char)              !"dim_char"

     CALL netcdf_err( nf90_def_dim(ncid, dim_names(5), nimp, idim_nimp),idim_nimp)                   !"dim_nimp"

     CALL netcdf_err( nf90_def_dim(ncid, dim_names(6), nneu, idim_nneu),idim_nneu)                   !"dim_neu"

     CALL netcdf_err( nf90_def_dim(ncid, dim_names(7), nbion, idim_nbion),idim_nbion)                 !"dim_fi"

     CALL netcdf_err( nf90_def_dim(ncid, dim_names(8), dischg%nplasbdry, &
          idim_nplasbdry),idim_nplasbdry)                   !"dim_nplasbdry"

     ntot = nion + dp4
     CALL netcdf_err( nf90_def_dim(ncid, dim_names(9), ntot, idim_ntot),idim_ntot) !"dim_ntot"

     CALL netcdf_err( nf90_def_dim(ncid, dim_names(10),dischg%nr_mhd, idim_nr_mhd),idim_nr_mhd) !"dim_nr_mhd"
     CALL netcdf_err( nf90_def_dim(ncid, dim_names(11),dischg%nz_mhd, idim_nz_mhd),idim_nz_mhd) !"dim_nz_mhd"

     CALL netcdf_err( nf90_def_dim(ncid, dim_names(12),mhd_dat%npsi, idim_npsi),idim_npsi) !"dim_npsi"

     CALL netcdf_err( nf90_def_dim(ncid, dim_names(13),dischg%nlimiter, idim_nlimiter),idim_nlimiter) !"dim_nlimiter"



     ! define variables:

     CALL netcdf_err( nf90_def_var(ncid, "shot", nf90_int,id_shot),id_shot)
     CALL netcdf_err(  nf90_put_att(ncid,id_shot,'long_name','shot number'),id_shot )

     CALL netcdf_err( nf90_def_var(ncid, "nj", nf90_int, id_nj))
     CALL netcdf_err(  nf90_put_att(ncid,id_nj,'long_name',           &
          'the size of quantities defined on rho grid ') )

     CALL netcdf_err( nf90_def_var(ncid, "time", nf90_double,id_time))
     CALL netcdf_err(  nf90_put_att(ncid,id_time,'long_name',           &
          'time at which data is printed ') )

     CALL netcdf_err( nf90_def_var(ncid, "tGCNMf", nf90_double,id_tGCNMf))
     CALL netcdf_err(  nf90_put_att(ncid,id_tGCNMf,'long_name',           &
          'GCNMP will evolve solution from time to fGCNMf (unless changed in namelist) ') )

     CALL netcdf_err( nf90_def_var(ncid, "time_bc", nf90_double,id_time_bc))
     CALL netcdf_err(  nf90_put_att(ncid,id_time_bc,'long_name',           &
          'Boundary condition time') )

 

     CALL netcdf_err( nf90_def_var(ncid, "psiaxis", nf90_double,id_psiaxis))
     CALL netcdf_err(  nf90_put_att(ncid,id_psiaxis,'long_name',           &
          'psi value on magnetic axis,volt sec/rad') )
     CALL netcdf_err(  nf90_put_att(ncid,id_psiaxis,'units', 'Volt sec/rad'  ))


     CALL netcdf_err( nf90_def_var(ncid, "psibdry", nf90_double,id_psibdry))
     CALL netcdf_err(  nf90_put_att(ncid,id_psibdry,'long_name',           &
          'psi value on plasma edge (separatrix),volt sec/rad') )
     CALL netcdf_err(  nf90_put_att(ncid,id_psibdry,'units', 'Volt sec/rad'  ))


     CALL netcdf_err( nf90_def_var(ncid, "rgeom", nf90_double,id_rgeom))
     CALL netcdf_err(  nf90_put_att(ncid,id_rgeom,'long_name',           &
          'major radius of geometric center at elevation of magnetic axis') )
     CALL netcdf_err(  nf90_put_att(ncid,id_rgeom,'units', 'meters'  ))


     CALL netcdf_err( nf90_def_var(ncid, "btgeom", nf90_double,id_btgeom))
     CALL netcdf_err(  nf90_put_att(ncid,id_btgeom,'long_name',           &
          'toroidal b field at geometric center rgeom,tesla') )
     CALL netcdf_err(  nf90_put_att(ncid,id_btgeom,'units', 'Tesla'  ))


     CALL netcdf_err( nf90_def_var(ncid, "rmag", nf90_double,id_rmag))


     CALL netcdf_err(  nf90_put_att(ncid,id_rmag,'long_name',           &
          'major radius of magnetic axis') )
     CALL netcdf_err(  nf90_put_att(ncid,id_rmag,'units', 'meters'  ))

     CALL netcdf_err( nf90_def_var(ncid, "zma", nf90_double,id_zma))
     CALL netcdf_err(  nf90_put_att(ncid,id_zma,'long_name',           &
          'major radius of magnetic axis') )
     CALL netcdf_err(  nf90_put_att(ncid,id_zma,'units', 'meters'  ))

     CALL netcdf_err( nf90_def_var(ncid, "rsep", nf90_double,id_rsep))
     CALL netcdf_err(  nf90_put_att(ncid,id_rsep,'long_name',           &
          'r of separatrix x point,m') )
     CALL netcdf_err(  nf90_put_att(ncid,id_rsep,'units', 'meters'  ))


     CALL netcdf_err( nf90_def_var(ncid, "zsep", nf90_double,id_zsep))
     CALL netcdf_err(  nf90_put_att(ncid,id_zsep,'long_name',           &
          'z of separatrix x point,m') )
     CALL netcdf_err(  nf90_put_att(ncid,id_zsep,'units', 'meters'  ))


     CALL netcdf_err( nf90_def_var(ncid, "rmajor", nf90_double,id_rmajor))
     CALL netcdf_err(  nf90_put_att(ncid,id_rmajor,'long_name',           &
          'major radius of vacuum BT0 reference == R0') )
     CALL netcdf_err(  nf90_put_att(ncid,id_rmajor,'units', 'meters'  ))

     CALL netcdf_err( nf90_def_var(ncid, "rplasmin", nf90_double,id_rplasmin))
     CALL netcdf_err(  nf90_put_att(ncid,id_rplasmin,'long_name',        &
          'min R of plasma,m') )
     CALL netcdf_err(  nf90_put_att(ncid,id_rplasmin,'units', 'meters'  ))

     CALL netcdf_err( nf90_def_var(ncid, "rplasmax", nf90_double,id_rplasmax))
     CALL netcdf_err(  nf90_put_att(ncid,id_rplasmax,'long_name',        &
          'max R of plasma,m') )
     CALL netcdf_err(  nf90_put_att(ncid,id_rplasmax,'units', 'meters'  ))

     CALL netcdf_err( nf90_def_var(ncid, "zplasmin", nf90_double,id_zplasmin))
     CALL netcdf_err(  nf90_put_att(ncid,id_zplasmin,'long_name',        &
          'min Z of plasma,m') )
     CALL netcdf_err(  nf90_put_att(ncid,id_zplasmin,'units', 'meters'  ))

     CALL netcdf_err( nf90_def_var(ncid, "zplasmax", nf90_double,id_zplasmax))
     CALL netcdf_err(  nf90_put_att(ncid,id_zplasmax,'long_name',        &
          'max Z of plasma,m') )
     CALL netcdf_err(  nf90_put_att(ncid,id_zplasmax,'units', 'meters'  ))



     CALL netcdf_err( nf90_def_var(ncid, "kappa", nf90_double,id_kappa))
     CALL netcdf_err(  nf90_put_att(ncid,id_kappa,'long_name',           &
          'plasma elongation') )

     CALL netcdf_err( nf90_def_var(ncid, "deltao", nf90_double,id_deltao))
     CALL netcdf_err(  nf90_put_att(ncid,id_deltao,'long_name',           &
          '*  deltao  : plasma(upper ) triangularity on axis ') )


     CALL netcdf_err( nf90_def_var(ncid, "pindento", nf90_double,id_pindento))
     CALL netcdf_err(  nf90_put_att(ncid,id_pindento,'long_name',           &
          '*  pindento  : on axis plasma indentation  ') )


     CALL netcdf_err( nf90_def_var(ncid, "volo", nf90_double,id_volume))
     CALL netcdf_err(  nf90_put_att(ncid,id_volume,'long_name',           &
          'plasma volume') )
     CALL netcdf_err(  nf90_put_att(ncid,id_volume,'units', 'meters^3'  ))

     CALL netcdf_err( nf90_def_var(ncid, "circum", nf90_double,id_circum))
     CALL netcdf_err(  nf90_put_att(ncid,id_circum,'long_name',           &
          'plasma crcumference,meters') )
     CALL netcdf_err(  nf90_put_att(ncid,id_circum,'units', 'meters'  ))


     CALL netcdf_err( nf90_def_var(ncid, "areao", nf90_double,id_areao))
     CALL netcdf_err(  nf90_put_att(ncid,id_areao,'long_name',           &
          'plasma cross sectional area') )
     CALL netcdf_err(  nf90_put_att(ncid,id_areao,'units', 'meters^2'  ))

     CALL netcdf_err( nf90_def_var(ncid, "nion", nf90_int, id_nion))
     CALL netcdf_err(  nf90_put_att(ncid,id_nion,'long_name',      &
          'the number of ion  species') )

     CALL netcdf_err( nf90_def_var(ncid, "nprim", nf90_int, id_nprim))
     CALL netcdf_err(  nf90_put_att(ncid,id_nprim,'long_name',      &
          'the number of primary ion species') )

     CALL netcdf_err( nf90_def_var(ncid, "fd_thermal", nf90_double, id_fd_thermal))
     CALL netcdf_err(  nf90_put_att(ncid,id_fd_thermal,'long_name',      &
          'fraction of d in thermal dt mixture if input as one species') )

     CALL netcdf_err( nf90_def_var(ncid, "nimp", nf90_int, id_nimp))
     CALL netcdf_err(  nf90_put_att(ncid,id_nimp,'long_name',      &
          'the number of impurity ion species') )

     CALL netcdf_err( nf90_def_var(ncid, "nneu", nf90_int, id_nneu))
     CALL netcdf_err(  nf90_put_att(ncid,id_nneu,'long_name',      &
          'the number of neutral species') )
     CALL netcdf_err( nf90_def_var(ncid, "ibion", nf90_int, id_nbion))
     CALL netcdf_err(  nf90_put_att(ncid,id_nbion,'long_name',      &
          'number of fast beam ion  species ') )

     CALL netcdf_err( nf90_def_var(ncid, "fd_beam", nf90_double, id_fd_beam))
     CALL netcdf_err(  nf90_put_att(ncid,id_fd_beam,'long_name',      &
          'fraction of d in dt fast ion mixture') )

     CALL netcdf_err( nf90_def_var(ncid, "namep", nf90_char,        &
          DIMIDS = (/idim_char,idim_nprim/),VARID = id_namep))
     CALL netcdf_err(  nf90_put_att(ncid,id_namep,'long_name',      &
          'name of primary ion species') )

     CALL netcdf_err( nf90_def_var(ncid, "namei", nf90_char,        &
          DIMIDS = (/idim_char,idim_nimp/),VARID= id_namei))
     CALL netcdf_err(  nf90_put_att(ncid,id_namei,'long_name',      &
          'name of imputity ion species') )

     CALL netcdf_err( nf90_def_var(ncid, "namen", nf90_char,        &
          DIMIDS = (/idim_char,idim_nneu/),VARID=id_namen))
     CALL netcdf_err(  nf90_put_att(ncid,id_namen,'long_name',      &
          'name of neutral species') )


     CALL netcdf_err( nf90_def_var(ncid, "nameb", nf90_char,        &
          DIMIDS = (/idim_char,idim_nbion/),VARID=   id_nameb))
     CALL netcdf_err(  nf90_put_att(ncid,id_nameb,'long_name',      &
          'name of beam  species') )


     CALL netcdf_err( nf90_def_var(ncid, "namepel", nf90_char,        &
          DIMIDS = (/idim_char/),VARID=   id_namepel))
     CALL netcdf_err(  nf90_put_att(ncid,id_namepel,'long_name',      &
          'name of pellet  species') )



     CALL netcdf_err( nf90_def_var(ncid, "btor", nf90_double,id_btor))
     CALL netcdf_err(  nf90_put_att(ncid,id_btor,'long_name',           &
          '*  Btor : vacuum toroidal field at R0, tesla') )
     CALL netcdf_err(  nf90_put_att(ncid,id_btor,'units', 'tesla'  ))


     CALL netcdf_err( nf90_def_var(ncid, "totcur", nf90_double,id_tot_cur))
     CALL netcdf_err(  nf90_put_att(ncid,id_tot_cur,'long_name',           &
          '*  total plasma current,amps') )
     CALL netcdf_err(  nf90_put_att(ncid,id_tot_cur,'units', 'amps'  ))

     CALL netcdf_err( nf90_def_var(ncid, "totohm", nf90_double,id_totohm_cur))
     CALL netcdf_err(  nf90_put_att(ncid,id_totohm_cur,'long_name',           &
          '* total ohmic plasma current,amps') )
     CALL netcdf_err(  nf90_put_att(ncid,id_totohm_cur,'units', 'amps'  ))


     CALL netcdf_err( nf90_def_var(ncid, "totboot", nf90_double,id_totboot_cur))
     CALL netcdf_err(  nf90_put_att(ncid,id_totboot_cur,'long_name',           &
          '* total bootstrap  current,amps') )
     CALL netcdf_err(  nf90_put_att(ncid,id_totboot_cur,'units', 'amps'  ))

     CALL netcdf_err( nf90_def_var(ncid, "totbeam", nf90_double,id_totbeam_cur))
     CALL netcdf_err(  nf90_put_att(ncid,id_totbeam_cur,'long_name',           &
          '* beam driven  current,amps') )
     CALL netcdf_err(  nf90_put_att(ncid,id_totbeam_cur,'units', 'amps'  ))

     CALL netcdf_err( nf90_def_var(ncid, "totrf", nf90_double,id_totrf_cur))
     CALL netcdf_err(  nf90_put_att(ncid,id_totrf_cur,'long_name',           &
          '* rf  driven  current,amps') )
     CALL netcdf_err(  nf90_put_att(ncid,id_totrf_cur,'units', 'amps'  ))

     CALL netcdf_err( nf90_def_var(ncid, "betap", nf90_double,id_betap))
     CALL netcdf_err(  nf90_put_att(ncid,id_betap,'long_name',           &
          '*  betap : poloidal beta') )


     CALL netcdf_err( nf90_def_var(ncid, "beta", nf90_double,id_beta))
     CALL netcdf_err(  nf90_put_att(ncid,id_beta,'long_name',           &
          '*  beta : toroidal beta') )


     CALL netcdf_err( nf90_def_var(ncid, "ali", nf90_double,id_ali))
     CALL netcdf_err(  nf90_put_att(ncid,id_ali,'long_name',           &
          '*  ali : plasma inductance') )


     CALL netcdf_err( nf90_def_var(ncid, "te0", nf90_double,id_te0))
     CALL netcdf_err(  nf90_put_att(ncid,id_te0,'long_name',           &
          '*  te0 : central electron temperature') )
     CALL netcdf_err(  nf90_put_att(ncid,id_te0,'units', 'Kev'  ))


     CALL netcdf_err( nf90_def_var(ncid, "ti0", nf90_double,id_ti0))
     CALL netcdf_err(  nf90_put_att(ncid,id_ti0,'long_name',           &
          '*  ti0 : central ion  temperature') )
     CALL netcdf_err(  nf90_put_att(ncid,id_ti0,'units', 'Kev'  ))


     CALL netcdf_err( nf90_def_var(ncid, "psi", nf90_double,    &
          DIMIDS = (/idim_nr_mhd,idim_nz_mhd/),VARID=id_psi))
     CALL netcdf_err( nf90_put_att(ncid,id_psi,'long_name',    &
          '*  psi on (R,Z) grid, volt*second/radian'))
     CALL netcdf_err(  nf90_put_att(ncid,id_psi,'units',       &
          'volt*second/radian'))


     CALL netcdf_err( nf90_def_var(ncid, "psir", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_psir_grid))
     CALL netcdf_err( nf90_put_att(ncid,id_psir_grid,'long_name',    &
          '*  psi on rho grid, volt*second/radian'))
     CALL netcdf_err(  nf90_put_att(ncid,id_psir_grid,'units',       &
          'volt*second/radian'))

     CALL netcdf_err( nf90_def_var(ncid, "rho_mhd_gridnpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_rho_mhd_gridnpsi))
     CALL netcdf_err( nf90_put_att(ncid,id_rho_mhd_gridnpsi,'long_name',    &
          '* rho grid corresponding to mhd psival grid, meters '))
     CALL netcdf_err(  nf90_put_att(ncid,id_rho_mhd_gridnpsi,'units',       &
          'meters'))


     CALL netcdf_err( nf90_def_var(ncid, "rmhdgrid", nf90_double,    &
          DIMIDS = (/idim_nr_mhd/),VARID=id_rmhdgrid))
     CALL netcdf_err( nf90_put_att(ncid,id_rmhdgrid,'long_name',    &
          '* R grid, meters '))
     CALL netcdf_err(  nf90_put_att(ncid,id_rmhdgrid,'units',       &
          '* meters'))

     CALL netcdf_err( nf90_def_var(ncid, "zmhdgrid", nf90_double,    &
          DIMIDS = (/idim_nz_mhd/),VARID=id_zmhdgrid))
     CALL netcdf_err( nf90_put_att(ncid,id_zmhdgrid,'long_name',    &
          '* Z grid, meters '))
     CALL netcdf_err(  nf90_put_att(ncid,id_zmhdgrid,'units',       &
          'meters'))

     CALL netcdf_err( nf90_def_var(ncid, "r", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_rho_grid))
     CALL netcdf_err( nf90_put_att(ncid,id_rho_grid,'long_name',    &
          '* rho grid, meters '))
     CALL netcdf_err(  nf90_put_att(ncid,id_rho_grid,'units',       &
          'meters'))


     CALL netcdf_err( nf90_def_var(ncid, "fcap", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_fcap))
     CALL netcdf_err( nf90_put_att(ncid,id_fcap,'long_name',    &
          '*  fcap, (i.e., f(psilim)/f(psi))'))

     CALL netcdf_err( nf90_def_var(ncid, "gcap", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_gcap))
     CALL netcdf_err( nf90_put_att(ncid,id_gcap,'long_name',    &
          '*  gcap, (i.e., <(grad rho)**2*(R0/R)**2>)'))



     CALL netcdf_err( nf90_def_var(ncid, "hcap", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_hcap))
     CALL netcdf_err( nf90_put_att(ncid,id_hcap,'long_name',    &
          '*  hcap, (i.e., (dvolume/drho)/(4*pi*pi*R0*rho'))

     CALL netcdf_err( nf90_def_var(ncid, "Te", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_te))
     CALL netcdf_err( nf90_put_att(ncid,id_te,'long_name',    &
          '*  electron temperature, keV'  ))
     CALL netcdf_err(  nf90_put_att(ncid,id_te,'units',       &
          'Kev'))

     CALL netcdf_err( nf90_def_var(ncid, "Ti", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_ti))
     CALL netcdf_err( nf90_put_att(ncid,id_ti,'long_name',    &
          '*  ion temperature, keV'  ))
     CALL netcdf_err(  nf90_put_att(ncid,id_ti,'units',       &
          'Kev'))

     CALL netcdf_err( nf90_def_var(ncid, "press", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_press))
     CALL netcdf_err( nf90_put_att(ncid,id_press,'long_name',    &
          '*  total pressure on transport rho grid nt/m^2'  ))
     CALL netcdf_err(  nf90_put_att(ncid,id_press,'units',       &
          'nt/m^2'))

     CALL netcdf_err( nf90_def_var(ncid, "pressb", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_pressb))
     CALL netcdf_err( nf90_put_att(ncid,id_pressb,'long_name',    &
          '* beam  pressure on transport rho grid nt/m^2'  ))
     CALL netcdf_err(  nf90_put_att(ncid,id_pressb,'units',       &
          'nt/m^2'))

     CALL netcdf_err( nf90_def_var(ncid, "q", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_q_value))
     CALL netcdf_err( nf90_put_att(ncid,id_q_value,'long_name',    &
          '*  safety factor '  ))

     CALL netcdf_err( nf90_def_var(ncid, "ene", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_ene))
     CALL netcdf_err( nf90_put_att(ncid,id_ene,'long_name',    &
          '*  electron density, #/meter**3'  ))
     CALL netcdf_err(  nf90_put_att(ncid,id_ene,'units',       &
          '1./meter^3'))



     !thermal primary ion densities:
     base_label = '* thermal primary ion  densities, species: '
     CALL set_label(label,base_label,strlen,'nprim')
     CALL netcdf_err( nf90_def_var(ncid, "enp", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nprim/),VARID=id_enp))
     CALL netcdf_err( nf90_put_att(ncid,id_enp,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_enp,'units',       &
          '1./meter^3'))

     !thermal impurity ion densities:
     base_label = '* thermal impurity ion  densities, species: '
     CALL set_label(label,base_label,strlen,'nimp')
     CALL netcdf_err( nf90_def_var(ncid, "eni", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nimp/),VARID=id_eni))
     CALL netcdf_err( nf90_put_att(ncid,id_eni,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_eni,'units',       &
          '1./meter^3'))

     !thermal ion particle flux:
     base_label = '* thermal ion  particle flux, species: '
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "p_flux", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_pflux))
     CALL netcdf_err( nf90_put_att(ncid,id_pflux,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_pflux,'units',       &
          '1./(meter^2 sec)'))

     !electron energy  flux:
     base_label = '* electron energy flux :  '
     CALL netcdf_err( nf90_def_var(ncid, "e_fluxe", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_efluxe))
     CALL netcdf_err( nf90_put_att(ncid,id_efluxe,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_efluxe,'units',       &
          'Kev/(meter^2 sec)'))

     !total thermal ion  energy  flux:
     base_label = '* total thermal ion energy flux :  '
     CALL netcdf_err( nf90_def_var(ncid, "e_fluxi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_efluxi))
     CALL netcdf_err( nf90_put_att(ncid,id_efluxi,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_efluxi,'units',       &
          'Kev/(meter^2 sec)'))

     !Flux associated with fdays law:
     base_label = '* Faradays law associated  flux :  '
     CALL netcdf_err( nf90_def_var(ncid, "fday_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_fdyflux))
     CALL netcdf_err( nf90_put_att(ncid,id_fdyflux,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_fdyflux,'units',       &
          'Tesla/(sec)'))

     !toroidal rotation flux:
     base_label = '* flux associated with toroidal rotation :  '
     CALL netcdf_err( nf90_def_var(ncid, "rot_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_rotflux))
     CALL netcdf_err( nf90_put_att(ncid,id_rotflux,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_rotflux,'units',       &
          'Kg/(sec^2)'))


     !tglf turbulent particle  flux:
     base_label = '* tglf turbulent  particle flux, species: '
     CALL set_label(label,base_label,strlen,'electron')
     CALL netcdf_err( nf90_def_var(ncid,"tglf_elct_p_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_tglf_p_fluxe))
     CALL netcdf_err( nf90_put_att(ncid,id_tglf_p_fluxe,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_tglf_p_fluxe,'units',       &
          '1./(meter^2 sec)'))

    !tglf turbulent particle  flux:
     base_label = '* tglf turbulent  particle flux, species: '
     CALL set_label(label,base_label,strlen,'effective prim ion')
     CALL netcdf_err( nf90_def_var(ncid, "tglf_ion_p_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_tglf_p_fluxp))
     CALL netcdf_err( nf90_put_att(ncid,id_tglf_p_fluxp,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_tglf_p_fluxp,'units',       &
          '1./(meter^2 sec)'))

    !tglf turbulent particle  flux:
     base_label = '* tglf turbulent  particle flux, species: '
     CALL set_label(label,base_label,strlen,'effective impurity')
     CALL netcdf_err( nf90_def_var(ncid,"tglf_imp_p_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_tglf_p_fluxi))
     CALL netcdf_err( nf90_put_att(ncid,id_tglf_p_fluxi,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_tglf_p_fluxi,'units',       &
          '1./(meter^2 sec)'))

    !tglf turbulent energy  flux:
     base_label = '* tglf turbulent  energy  flux, species: '
     CALL set_label(label,base_label,strlen,'electron')
     CALL netcdf_err( nf90_def_var(ncid, "tglf_elc_e_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_tglf_e_fluxe))
     CALL netcdf_err( nf90_put_att(ncid,id_tglf_e_fluxe,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_tglf_e_fluxe,'units',       &
          'KEV/(meter^2 sec)'))

    !tglf turbulent energy  flux:
     base_label = '* tglf turbulent  energy  flux, species: '
     CALL set_label(label,base_label,strlen,'effective prim ion')
     CALL netcdf_err( nf90_def_var(ncid, "tglf_ion_e_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_tglf_e_fluxp))
     CALL netcdf_err( nf90_put_att(ncid,id_tglf_e_fluxp,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_tglf_e_fluxp,'units',       &
          'KEV/(meter^2 sec)'))

    !tglf turbulent energy  flux:
     base_label = '* tglf turbulent  energy  flux, species: '
     CALL set_label(label,base_label,strlen,'effective impurity')
     CALL netcdf_err( nf90_def_var(ncid,"tglf_imp_e_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_tglf_e_fluxi))
     CALL netcdf_err( nf90_put_att(ncid,id_tglf_e_fluxi,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_tglf_e_fluxi,'units',       &
          'KEV/(meter^2 sec)'))

     !electron density time deriv:
     CALL netcdf_err( nf90_def_var(ncid,"dnedt", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_dnedt))
     WRITE(label,FMT='(1pe16.8)')tGCNMs
     READ(label,FMT='(a)')st_asc_time
     st_asc_time = ADJUSTL(st_asc_time)
     label = '* rate of change of ene  1/(m**3 sec) at time ' &
          //st_asc_time(1:LEN_TRIM(st_asc_time))
     CALL netcdf_err(  nf90_put_att(ncid,id_dnedt,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_dnedt ,'units',        &
          '1/meter^3 sec'))

!!$     !thermal ion densities time derivatives:
!!$     base_label = '* thermal ion dni/dt by  species: '
!!$     CALL set_label(label,base_label,strlen,'nion')
!!$     CALL netcdf_err( nf90_def_var(ncid, "dnidt", nf90_double,    &
!!$          DIMIDS = (/idim_rho,idim_nion/),VARID=id_dnidt))
!!$     CALL netcdf_err( nf90_put_att(ncid,id_dnidt,'long_name',label ))
!!$     CALL netcdf_err(  nf90_put_att(ncid,id_dnidt,'units',       &
!!$          '1./(meter^3 sec)'))

     !thermal ion densities boundary condititons:
     base_label = " bc profles at time:"//bc_asc_time(1:LEN_TRIM(bc_asc_time)) &
          //',species: '
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "en_bc", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_en_bc))

     CALL netcdf_err( nf90_put_att(ncid,id_en_bc,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_en_bc,'units','1/meter^3'))

     !thermal ion flux(or mixed)  boundary condititons:
     base_label = " bc profles at time:"//bc_asc_time(1:LEN_TRIM(bc_asc_time)) &
          //',flux for species: '
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "flux_bc", nf90_double,    &
          DIMIDS = (/idim_rho,idim_ntot/),VARID=id_flux_bc))

     CALL netcdf_err( nf90_put_att(ncid,id_flux_bc,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_flux_bc,'units','?/(meter^2sec)'))

     !thermal ion source:
     base_label = '  sion : source due to ionization,'
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "sion", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_sion))
     CALL netcdf_err( nf90_put_att(ncid,id_sion,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_sion,'units','1/(meter^3 sec)'))

     !thermal ion recombination sink:
     base_label = '  srecom : (-)source due to recombinatioin,'
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "srecom", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_srecom))
     CALL netcdf_err( nf90_put_att(ncid,id_srecom,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_srecom,'units','1/(meter^3 sec)'))

     !source due to cx thermal neut
     base_label = '*  scx : source due to cx thermal neut.,'
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "scx", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_scx))
     CALL netcdf_err( nf90_put_att(ncid,id_scx,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_scx,'units','1/(meter^3 sec)'))

     base_label = '*  sbcx : sink due to cx with beam neut.,'
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "sbcx", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_sbcx))
     CALL netcdf_err( nf90_put_att(ncid,id_sbcx,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_sbcx,'units','1/(meter^3 sec)'))


     base_label = '*  stsource  : total source rate,'
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "s", nf90_double,    &
          DIMIDS = (/idim_nion,idim_rho/),VARID=id_stsource))
     CALL netcdf_err( nf90_put_att(ncid,id_stsource,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_stsource,'units','1/(meter^3 sec)'))


     base_label = '*  dudt : s dot, (# or  energy) /(meter**3*second)'
     CALL set_label(label,base_label,strlen,'ntot')
     CALL netcdf_err( nf90_def_var(ncid, "dudtsv", nf90_double,    &
          DIMIDS = (/idim_ntot,idim_rho/),VARID=id_dudtsv))
     CALL netcdf_err( nf90_put_att(ncid,id_dudtsv,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_dudtsv,'units','(1 or kev)/(meter^3 sec)'))

     !neutral  densities:
     base_label = '*  neutral density, #/meter**3, species: '
     CALL set_label(label,base_label,strlen,'nneu')
     CALL netcdf_err( nf90_def_var(ncid, "enn", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nneu/),VARID=id_enn))
     CALL netcdf_err( nf90_put_att(ncid,id_enn,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_enn,'units',       &
          'meter^3'))

     base_label = '*  neutral density,due to wall source , species: '
     CALL set_label(label,base_label,strlen,'nneu')
     CALL netcdf_err( nf90_def_var(ncid, "ennw", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nneu/),VARID=id_ennw))
     CALL netcdf_err( nf90_put_att(ncid,id_ennw,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_ennw,'units',       &
          'meter^3'))

     base_label = '*  neutral density,due to volume  source , species: '
     CALL set_label(label,base_label,strlen,'nneu')
     CALL netcdf_err( nf90_def_var(ncid, "ennv", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nneu/),VARID=id_ennv))
     CALL netcdf_err( nf90_put_att(ncid,id_ennv,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_ennv,'units',       &
          'meter^3'))

     base_label = '* volume  source of neutrals , species: '
     CALL set_label(label,base_label,strlen,'nneu')
     CALL netcdf_err( nf90_def_var(ncid, "volsn", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nneu/),VARID=id_volsn))
     CALL netcdf_err( nf90_put_att(ncid,id_volsn,'long_name',label  ))
     CALL netcdf_err(  nf90_put_att(ncid,id_volsn,'units',       &
          '#/(meter**3*second)'))


     !z and zsq:
     base_label = '* charge  z ,species : '
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "z", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_z))
     CALL netcdf_err( nf90_put_att(ncid,id_z,'long_name',label))
     CALL netcdf_err( nf90_def_var(ncid, "zsq", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_zsq))
     CALL netcdf_err( nf90_put_att(ncid,id_zsq,'long_name',label))


     CALL netcdf_err( nf90_def_var(ncid, "stfuse", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_stfuse))
     CALL netcdf_err( nf90_put_att(ncid,id_stfuse,'long_name',    &
          '*  stfuse : thermal fusion rate,#/(meter**3*second)'  ))
     CALL netcdf_err(  nf90_put_att(ncid,id_stfuse,'units',       &
          '#/(meter**3*second)'))


     CALL netcdf_err( nf90_def_var(ncid, "sbfuse", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_sbfuse))
     CALL netcdf_err( nf90_put_att(ncid,id_sbfuse,'long_name',    &
          '*  sbfuse : beam fusion rate,#/(meter**3*second)'  ))
     CALL netcdf_err(  nf90_put_att(ncid,id_sbfuse,'units',       &
          '#/(meter**3*second)'))


     CALL netcdf_err( nf90_def_var(ncid, "spellet", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_spellet))
     CALL netcdf_err( nf90_put_att(ncid,id_spellet,'long_name',    &
          '*  spellet : ion pellet source,#/(meter**3*second)'  ))
     CALL netcdf_err(  nf90_put_att(ncid,id_spellet,'units',       &
          '#/(meter**3*second)'))




     CALL netcdf_err( nf90_def_var(ncid, "sbion", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_sbeame))
     CALL netcdf_err( nf90_put_att(ncid,id_sbeame,'long_name',    &
          '*  sbeame : beam electron source,#/(meter**3*second)'  ))
     CALL netcdf_err(  nf90_put_att(ncid,id_sbeame,'units',       &
          '#/(meter**3*second)'))

     base_label = '*  sbeam : beam ion  source,#/(meter**3*second), species: '
     CALL set_label(label,base_label,strlen,'nf')
     CALL netcdf_err( nf90_def_var(ncid, "sbeam", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nbion/),VARID=id_sbeam))
     CALL netcdf_err( nf90_put_att(ncid,id_sbeam,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_sbeam,'units',       &
          '#/(meter**3*second)'))

     base_label = '*  fast ion density, #/meter**3, species: '
     CALL set_label(label,base_label,strlen,'nf')
     CALL netcdf_err( nf90_def_var(ncid, "enbeam", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nbion/),VARID=id_enbeam))
     CALL netcdf_err( nf90_put_att(ncid,id_enbeam,'long_name',label  ))
     CALL netcdf_err(  nf90_put_att(ncid,id_enbeam,'units',       &
          '#/(meter**3)'))



     CALL netcdf_err( nf90_def_var(ncid, "curden", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_curden))
     CALL netcdf_err(  nf90_put_att(ncid,id_curden,'long_name',    &
          '*  total torpodal current density,<Jtor R0/R>, amps/meter**2'))
     CALL netcdf_err(  nf90_put_att(ncid,id_curden,'units','amps/meter^2)'))

     CALL netcdf_err( nf90_def_var(ncid, "curohm", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_curohm))
     CALL netcdf_err(  nf90_put_att(ncid,id_curohm,'long_name',    &
          '* ohmic current density, amps/meter**2'))
     CALL netcdf_err(  nf90_put_att(ncid,id_curohm,'units',        &
          'amps/meter^2)'))

     CALL netcdf_err( nf90_def_var(ncid, "curboot", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_curboot))
     CALL netcdf_err(  nf90_put_att(ncid,id_curboot,'long_name',    &
          '* bootstrap current density, amps/meter**2'))
     CALL netcdf_err(  nf90_put_att(ncid,id_curboot,'units',        &
          'amps/meter^2)'))

     CALL netcdf_err( nf90_def_var(ncid, "curdbeam", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_curbeam))
     CALL netcdf_err(  nf90_put_att(ncid,id_curbeam,'long_name',    &
          '* beam driven  current density, amps/meter**2'))
     CALL netcdf_err(  nf90_put_att(ncid,id_curbeam,'units',        &
          'amps/meter^2)'))


     CALL netcdf_err( nf90_def_var(ncid, "currf", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_currf))
     CALL netcdf_err(  nf90_put_att(ncid,id_currf,'long_name',    &
          '* rf driven  current density, amps/meter**2'))
     CALL netcdf_err(  nf90_put_att(ncid,id_currf,'units',        &
          'amps/meter^2)'))


     CALL netcdf_err( nf90_def_var(ncid, "etor", nf90_double,     &
          DIMIDS = (/idim_rho/),VARID=id_etor))
     CALL netcdf_err(  nf90_put_att(ncid,id_etor,'long_name',    &
          '*  toroidal electric field profile, V/m'))
     CALL netcdf_err(  nf90_put_att(ncid,id_etor,'units',        &
          'Volts/meter)'))



     CALL netcdf_err( nf90_def_var(ncid, "rbp", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_rbp))
     CALL netcdf_err(  nf90_put_att(ncid,id_rbp,'long_name',    &
          '*  rho*bp0*fcap*gcap*hcap, tesla*meters'))
     CALL netcdf_err(  nf90_put_att(ncid,id_rbp,'units',        &
          'tesla*meters)'))

     CALL netcdf_err( nf90_def_var(ncid, "psivalnpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_psivalnpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_psivalnpsi,'long_name',    &
          '* psivalnpsi(npsi) grid,edge to amgnetic axis,volt sec/rad'))
     CALL netcdf_err(  nf90_put_att(ncid,id_psivalnpsi,'units','Volt sec/rad'))


     CALL netcdf_err( nf90_def_var(ncid, "ravgnpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_ravgnpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_ravgnpsi,'long_name',    &
          '* ravg: <R> avg radius on mhd grid,m'))
     CALL netcdf_err(  nf90_put_att(ncid,id_ravgnpsi,'units','meters'))

     CALL netcdf_err( nf90_def_var(ncid, "ravginpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_ravginpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_ravginpsi,'long_name',    &
          '* ravgi:<1/R>  on mhd grid,1/m'))
     CALL netcdf_err(  nf90_put_att(ncid,id_ravginpsi,'units','meters'))


     CALL netcdf_err( nf90_def_var(ncid, "fpsinpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_fpsinpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_fpsinpsi,'long_name',    &
          '* fpsi: f of psi (= R*Bt) on psivalnpsi(npsi) grid,tesla meters'))
     CALL netcdf_err(  nf90_put_att(ncid,id_fpsinpsi,'units','Tesla meters'))


     CALL netcdf_err( nf90_def_var(ncid, "pprim", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_pprim))
     CALL netcdf_err(  nf90_put_att(ncid,id_pprim,'long_name',    &
          '* pprim: dp/dpsi on transport grid,nt/(m**2-volt-sec) = amp/m**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_pprim,'units','nt/(m**2-volt-sec'))


     CALL netcdf_err( nf90_def_var(ncid, "ffprim", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_ffprim))
     CALL netcdf_err(  nf90_put_att(ncid,id_ffprim,'long_name',    &
          '* f*df/dpsi: on transport grid      kg/(A sec^2 '))
     CALL netcdf_err(  nf90_put_att(ncid,id_ffprim,'units', 'kg/(A sec^2)'))


     CALL netcdf_err( nf90_def_var(ncid, "bp", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_bp))
     CALL netcdf_err(  nf90_put_att(ncid,id_bp,'long_name',    &
          '*<Bp> : flux avg B poloidal on transport grid, Tesla'))
     CALL netcdf_err(  nf90_put_att(ncid,id_bp,'units',        &
          'tesla'))

     CALL netcdf_err( nf90_def_var(ncid, "bprmaj", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_bprmaj))
     CALL netcdf_err(  nf90_put_att(ncid,id_bprmaj,'long_name',    &
          '* B poloidal on rmaj (and transport rho)  grid, tesla'))
     CALL netcdf_err(  nf90_put_att(ncid,id_bprmaj,'units',        &
          'tesla)'))


     CALL netcdf_err( nf90_def_var(ncid, "btotrmaj", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_btotrmaj))
     CALL netcdf_err(  nf90_put_att(ncid,id_btotrmaj,'long_name',    &
          '* Btotal on transport rho grid, tesla'))
     CALL netcdf_err(  nf90_put_att(ncid,id_btotrmaj,'units',        &
          'tesla)'))



     CALL netcdf_err( nf90_def_var(ncid, "zeff", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_zeff))
     CALL netcdf_err(  nf90_put_att(ncid,id_zeff,'long_name',    &
          '*  zeff profile'))


     CALL netcdf_err( nf90_def_var(ncid, "angrot", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_angrot))
     CALL netcdf_err(  nf90_put_att(ncid,id_angrot,'long_name',    &
          '*  angular rotation speed profile, rad/sec'))
     CALL netcdf_err(  nf90_put_att(ncid,id_angrot,'units',        &
          'rad/sec'))

     !new 8/06/07:
     CALL netcdf_err( nf90_def_var(ncid, "d", nf90_double,    &
          DIMIDS = (/idim_ntot,idim_ntot,idim_rho/),VARID=id_d))
     CALL netcdf_err(  nf90_put_att(ncid,id_d,'long_name',    &
          '*  diffusion matrix(ntot,ntot,nj) , on half grid'))
     CALL netcdf_err(  nf90_put_att(ncid,id_d,'units',        &
          'mixed'))


     CALL netcdf_err( nf90_def_var(ncid, "chieinv", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_chieinv))
     CALL netcdf_err(  nf90_put_att(ncid,id_chieinv,'long_name',    &
          '*  electron thermal diffusivity, meters**2/sec, on half grid'))
     CALL netcdf_err(  nf90_put_att(ncid,id_chieinv,'units',        &
          'meters^2/sec'))

     CALL netcdf_err( nf90_def_var(ncid, "chiinv", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_chiinv))
     CALL netcdf_err(  nf90_put_att(ncid,id_chiinv,'long_name',    &
          '*  ion thermal diffusivity, meters**2/sec, on half grid'))
     CALL netcdf_err(  nf90_put_att(ncid,id_chiinv,'units',        &
          'meters^2/sec'))

     CALL netcdf_err( nf90_def_var(ncid, "xkineo", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_xkineo))
     CALL netcdf_err(  nf90_put_att(ncid,id_xkineo,'long_name',    &
          '*  ion neoclassical thermal conductivity'))
     CALL netcdf_err(  nf90_put_att(ncid,id_xkineo,'units',        &
          ' 1/(meter*second'))

     CALL netcdf_err( nf90_def_var(ncid, "xkeneo", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_xkeneo))
     CALL netcdf_err(  nf90_put_att(ncid,id_xkeneo,'long_name',    &
          '*  electron neoclassical thermal conductivity'))
     CALL netcdf_err(  nf90_put_att(ncid,id_xkeneo,'units',        &
          ' 1/(meter*second'))


     CALL netcdf_err( nf90_def_var(ncid, "dpedtc", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_dpedt))
     CALL netcdf_err(  nf90_put_att(ncid,id_dpedt,'long_name',    &
          '*  1.5*dpedt, power density due to electron pressure, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_dpedt,'units',        &
          ' watts/meter**3'))


     CALL netcdf_err( nf90_def_var(ncid, "dpidtc", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_dpidt))
     CALL netcdf_err(  nf90_put_att(ncid,id_dpidt,'long_name',    &
          '* 1.5*dpidt,power density due to ion pressure, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_dpidt,'units',        &
          ' watts/meter**3'))

     CALL netcdf_err( nf90_def_var(ncid, "qconde", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qconde))
     CALL netcdf_err(  nf90_put_att(ncid,id_qconde,'long_name',    &
          '*  electron conduction, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qconde,'units',        &
          ' watts/meter**3'))


     CALL netcdf_err( nf90_def_var(ncid, "qcondi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qcondi))
     CALL netcdf_err(  nf90_put_att(ncid,id_qcondi,'long_name',    &
          '*  iom conduction, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qcondi,'units',        &
          ' watts/meter**3'))

     CALL netcdf_err( nf90_def_var(ncid, "qconve", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qconve))
     CALL netcdf_err(  nf90_put_att(ncid,id_qconve,'long_name',    &
          '*  electron convection, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qconve,'units',        &
          ' watts/meter**3'))

     CALL netcdf_err( nf90_def_var(ncid, "qconvi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qconvi))
     CALL netcdf_err(  nf90_put_att(ncid,id_qconvi,'long_name',    &
          '*  ion convection, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qconvi,'units',        &
          ' watts/meter**3'))


     CALL netcdf_err( nf90_def_var(ncid, "qbeame", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qbeame))
     CALL netcdf_err(  nf90_put_att(ncid,id_qbeame,'long_name',    &
          '*  power to elec. from beam, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qbeame,'units',        &
          ' watts/meter**3'))

     CALL netcdf_err( nf90_def_var(ncid, "qbeami", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qbeami))
     CALL netcdf_err(  nf90_put_att(ncid,id_qbeami,'long_name',    &
          '*  power to ions  from beam, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qbeami,'units',        &
          ' watts/meter**3'))

     CALL netcdf_err( nf90_def_var(ncid, "qdelt", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qdelt))
     CALL netcdf_err(  nf90_put_att(ncid,id_qdelt,'long_name',    &
          '*  power to ions  from beam, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qdelt,'units',        &
          ' watts/meter**3'))

     CALL netcdf_err( nf90_def_var(ncid, "qexch", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qexch))
     CALL netcdf_err(  nf90_put_att(ncid,id_qexch,'long_name',    &
          '* qexch,anomalous electron-ion energy exchange term, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qexch,'units',        &
          ' watts/meter**3'))


     CALL netcdf_err( nf90_def_var(ncid, "qrfe", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qrfe))
     CALL netcdf_err(  nf90_put_att(ncid,id_qrfe,'long_name',    &
          '*  qrfe, RF electron heating, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qrfe,'units',        &
          ' watts/meter**3'))


     CALL netcdf_err( nf90_def_var(ncid, "qrfi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qrfi))
     CALL netcdf_err(  nf90_put_att(ncid,id_qrfi,'long_name',    &
          '*  qrfi, RF electron heating, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qrfi,'units',        &
          ' watts/meter**3'))


     CALL netcdf_err( nf90_def_var(ncid, "qione", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qione))
     CALL netcdf_err(  nf90_put_att(ncid,id_qione,'long_name',    &
          '*  qione, electron power density due to  recombination and impact ionization, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qione,'units',        &
          ' watts/meter**3'))


     CALL netcdf_err( nf90_def_var(ncid, "qioni", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qioni))
     CALL netcdf_err(  nf90_put_att(ncid,id_qioni,'long_name',    &
          '*  qioni, ion power density due to  recombination and impact ionization,watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qioni,'units',        &
          ' watts/meter**3'))


     CALL netcdf_err( nf90_def_var(ncid, "qcx", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qcx))
     CALL netcdf_err(  nf90_put_att(ncid,id_qcx,'long_name',    &
          '*  qcx, ion power density due to neutral-ion charge exchange, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qcx,'units',        &
          ' watts/meter**3'))



     CALL netcdf_err( nf90_def_var(ncid, "qe2d", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qe2d))
     CALL netcdf_err(  nf90_put_att(ncid,id_qe2d,'long_name',    &
          '*  2d electron heating, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qe2d,'units',        &
          ' watts/meter**3'))



     CALL netcdf_err( nf90_def_var(ncid, "qi2d", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qi2d))
     CALL netcdf_err(  nf90_put_att(ncid,id_qi2d,'long_name',    &
          '*  2d ion heating, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qi2d,'units',        &
          ' watts/meter**3'))



     CALL netcdf_err( nf90_def_var(ncid, "qfuse", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qfuse))
     CALL netcdf_err(  nf90_put_att(ncid,id_qfuse,'long_name',    &
          '* total  fusion electron heating, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qfuse,'units',        &
          ' watts/meter**3'))


     CALL netcdf_err( nf90_def_var(ncid, "qfusi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qfusi))
     CALL netcdf_err(  nf90_put_att(ncid,id_qfusi,'long_name',    &
          '* total  fusion ion heating, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qfusi,'units',        &
          ' watts/meter**3'))



     CALL netcdf_err( nf90_def_var(ncid, "qbfuse", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qbfuse))
     CALL netcdf_err(  nf90_put_att(ncid,id_qbfuse,'long_name',    &
          '*  beam fusion electron heating, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qbfuse,'units',        &
          ' watts/meter**3'))



     CALL netcdf_err( nf90_def_var(ncid, "qbfusi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qbfusi))
     CALL netcdf_err(  nf90_put_att(ncid,id_qfusi,'long_name',    &
          '* beam  fusion ion heating, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qbfusi,'units',        &
          ' watts/meter**3'))



     CALL netcdf_err( nf90_def_var(ncid,"qmag", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qmag))
     CALL netcdf_err(  nf90_put_att(ncid,id_qmag,'long_name',    &
          '*  qmag electron heating, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qmag,'units',        &
          ' watts/meter**3'))


     CALL netcdf_err( nf90_def_var(ncid,"qsawe", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qsawe))
     CALL netcdf_err(  nf90_put_att(ncid,id_qsawe,'long_name',    &
          '*  sawtooth electron heating, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qsawe,'units',        &
          ' watts/meter**3'))

     CALL netcdf_err( nf90_def_var(ncid,"qsawi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qsawi))
     CALL netcdf_err(  nf90_put_att(ncid,id_qsawi,'long_name',    &
          '*  sawtooth ion heating, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qsawi,'units',        &
          ' watts/meter**3'))


     CALL netcdf_err( nf90_def_var(ncid,"qrad", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qrad))
     CALL netcdf_err(  nf90_put_att(ncid,id_qrad,'long_name',    &
          '*  radiated power density, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qrad,'units',        &
          ' watts/meter**3'))



     CALL netcdf_err( nf90_def_var(ncid,"qohm", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qohm))
     CALL netcdf_err(  nf90_put_att(ncid,id_qohm,'long_name',    &
          '*  (electron) ohmic power density, watts/meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qohm,'units',        &
          ' watts/meter**3'))
     CALL netcdf_err( nf90_def_var(ncid,"rmajavnpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_rmajavnpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_rmajavnpsi,'long_name',    &
          '*  average major radius of each flux surface,meters, evaluated at elevation of magnetic axis'))
     CALL netcdf_err(  nf90_put_att(ncid,id_rmajavnpsi,'units',        &
          ' meters'))


     CALL netcdf_err( nf90_def_var(ncid,"rminavnpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_rminavnpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_rminavnpsi,'long_name',    &
          '*  average minor radius of each flux surface,meters, evaluated at elevation of magnetic axis'))
     CALL netcdf_err(  nf90_put_att(ncid,id_rminavnpsi,'units',        &
          ' meters'))


     CALL netcdf_err( nf90_def_var(ncid,"psivolp", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_psivolp))
     CALL netcdf_err(  nf90_put_att(ncid,id_psivolp,'long_name',    &
          '*  volume of each flux surface on rho  grid, meters**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_psivolp,'units',        &
          ' meters^3'))

     CALL netcdf_err( nf90_def_var(ncid,"psivolpnpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_psivolpnpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_psivolpnpsi,'long_name',    &
          '*  volume of each flux surface on npsi psi grid, meters**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_psivolpnpsi,'units',        &
          ' meters^3'))

     CALL netcdf_err( nf90_def_var(ncid,"elongxnpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_elongxnpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_elongxnpsi,'long_name',    &
          '*  elongation of each flux surface'))

     CALL netcdf_err( nf90_def_var(ncid,"elongx", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_elongx))
     CALL netcdf_err(  nf90_put_att(ncid,id_elongx,'long_name',    &
          '*  elongation of each flux surface on rho grid'))

     CALL netcdf_err( nf90_def_var(ncid,"triangnpsi_u", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_triangnpsi_u))
     CALL netcdf_err(  nf90_put_att(ncid,id_triangnpsi_u,'long_name',    &
          '*  upper triangularity of each flux surface'))

     CALL netcdf_err( nf90_def_var(ncid,"triangnpsi_l", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_triangnpsi_l))
     CALL netcdf_err(  nf90_put_att(ncid,id_triangnpsi_l,'long_name',    &
          '*  lower triangularity of each flux surface'))


     CALL netcdf_err( nf90_def_var(ncid,"triangnpsi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_triangnpsi_lcl))
     CALL netcdf_err(  nf90_put_att(ncid,id_triangnpsi_lcl,'long_name',    &
          '*  average  triangularity of each flux surface'))

     CALL netcdf_err( nf90_def_var(ncid,"pindentnpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_pindentnpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_pindentnpsi,'long_name',    &
          '*  indentation of each flux surface'))

     CALL netcdf_err( nf90_def_var(ncid,"sfareanpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_sfareanpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_sfareanpsi,'long_name',    &
          '*  surface area of each flux surface,this is 4*pi*pi*R0*hcap*rho*<ABS(grad rho)>'))
     CALL netcdf_err(  nf90_put_att(ncid,id_sfareanpsi,'units',        &
          ' meters^2'))

     CALL netcdf_err( nf90_def_var(ncid,"cxareanpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_cxareanpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_cxareanpsi,'long_name',    &
          '*  cross-sectional area of each flux'))
     CALL netcdf_err(  nf90_put_att(ncid,id_cxareanpsi,'units',        &
          ' meters^2'))


     CALL netcdf_err( nf90_def_var(ncid,"grho1npsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_grho1npsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_grho1npsi,'long_name',    &
          '*  flux surface average absolute grad rho'))

     CALL netcdf_err( nf90_def_var(ncid,"grho2npsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_grho2npsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_grho2npsi,'long_name',    &
          '*  flux surface average (grad rho)**2'))

     CALL netcdf_err( nf90_def_var(ncid, "nplasbdry", nf90_int, id_nplasbdry))
     CALL netcdf_err(  nf90_put_att(ncid,id_nplasbdry,'long_name',      &
                   '*  nplasdry : number of points on plasma boundary') )




     CALL netcdf_err( nf90_def_var(ncid,"rplasbdry", nf90_double,         &
          DIMIDS = (/idim_nplasbdry/),VARID=id_rplasbdry))
     CALL netcdf_err(  nf90_put_att(ncid,id_rplasbdry,'long_name',    &
          '*  r points for plasma boundary, meters'))
     CALL netcdf_err(  nf90_put_att(ncid,id_rplasbdry,'units',        &
          ' meters'))
     CALL netcdf_err( nf90_def_var(ncid,"zplasbdry", nf90_double,         &
          DIMIDS = (/idim_nplasbdry/),VARID=id_zplasbdry))
     CALL netcdf_err(  nf90_put_att(ncid,id_zplasbdry,'long_name',    &
          '*  z points for plasma boundary, meters'))
     CALL netcdf_err(  nf90_put_att(ncid,id_zplasbdry,'units',        &
          ' meters'))

     CALL netcdf_err( nf90_def_var(ncid,"rlimiter", nf90_double,         &
          DIMIDS = (/idim_nlimiter/),VARID=id_rlimiter))
     CALL netcdf_err(  nf90_put_att(ncid,id_rlimiter,'long_name',    &
          '*  R points for limiter, meters'))
     CALL netcdf_err(  nf90_put_att(ncid,id_rlimiter,'units',        &
          ' meters'))

     CALL netcdf_err( nf90_def_var(ncid,"zlimiter", nf90_double,         &
          DIMIDS = (/idim_nlimiter/),VARID=id_zlimiter))
     CALL netcdf_err(  nf90_put_att(ncid,id_zlimiter,'long_name',    &
          '*  Z points for limiter, meters'))
     CALL netcdf_err(  nf90_put_att(ncid,id_zlimiter,'units',        &
          ' meters'))

     CALL netcdf_err( nf90_def_var(ncid,"storqueb", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_storqueb))
     CALL netcdf_err(  nf90_put_att(ncid,id_storqueb,'long_name',    &
          '* beam   torque density, nt-m/m**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_storqueb,'units',        &
          ' NT-meters/meter^3'))


     CALL netcdf_err( nf90_def_var(ncid,"Kpol_c", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_kpol_c))
     CALL netcdf_err(  nf90_put_att(ncid,id_kpol_c,'long_name',    &
          '* CER impurity poloidal velocity/Bpol'))
     CALL netcdf_err(  nf90_put_att(ncid,id_kpol_c,'units',        &
          ' meters/second/Tesla'))


     CALL netcdf_err( nf90_def_var(ncid,"Kpol_d", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_kpol_d))
     CALL netcdf_err(  nf90_put_att(ncid,id_kpol_d,'long_name',    &
          '* main ion poloidal velocity/Bpol'))
     CALL netcdf_err(  nf90_put_att(ncid,id_kpol_d,'units',        &
          ' meters/second/Tesla'))


     CALL netcdf_err( nf90_def_var(ncid,"Kpol_exp", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_kpol_exp))
     CALL netcdf_err(  nf90_put_att(ncid,id_kpol_exp,'long_name',    &
          '* Measured CER impurity poloidal velocity/Bpol'))
     CALL netcdf_err(  nf90_put_att(ncid,id_kpol_exp,'units',        &
          ' meters/second/Tesla'))


     CALL netcdf_err( nf90_def_var(ncid,"angrot_c", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_angrot_c))
     CALL netcdf_err(  nf90_put_att(ncid,id_angrot_c,'long_name',    &
          '* CER impurity toroidal angular velocity'))
     CALL netcdf_err(  nf90_put_att(ncid,id_angrot_c,'units',        &
          ' radians/second'))


     CALL netcdf_err( nf90_def_var(ncid,"angrot_d", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_angrot_d))
     CALL netcdf_err(  nf90_put_att(ncid,id_angrot_d,'long_name',    &
          '* main ion toroidal angular velocity'))
     CALL netcdf_err(  nf90_put_att(ncid,id_angrot_d,'units',        &
          ' radians/second'))


     CALL netcdf_err( nf90_def_var(ncid,"ave_vpar_d", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_ave_vpar_d))
     CALL netcdf_err(  nf90_put_att(ncid,id_ave_vpar_d,'long_name',    &
          '* main ion average parallel velocity'))
     CALL netcdf_err(  nf90_put_att(ncid,id_ave_vpar_d,'units',        &
          ' meters/second'))


     CALL netcdf_err( nf90_def_var(ncid,"udia_d", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_udia_d))
     CALL netcdf_err(  nf90_put_att(ncid,id_udia_d,'long_name',    &
          '* main ion diamagnetic velocity/RBpol'))
     CALL netcdf_err(  nf90_put_att(ncid,id_udia_d,'units',        &
          '1/second'))

     CALL netcdf_err( nf90_def_var(ncid,"udia_c", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_udia_c))
     CALL netcdf_err(  nf90_put_att(ncid,id_udia_c,'long_name',    &
          '* CER ion diamagnetic velocity/RBpol'))
     CALL netcdf_err(  nf90_put_att(ncid,id_udia_c,'units',        &
          '1/second'))

     CALL netcdf_err( nf90_def_var(ncid,"ugrt", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_ugrt))
     CALL netcdf_err(  nf90_put_att(ncid,id_ugrt,'long_name',    &
          '* ion temperature gradient/(eRBpol)'))
     CALL netcdf_err(  nf90_put_att(ncid,id_ugrt,'units',        &
          '1/second'))

     CALL netcdf_err( nf90_def_var(ncid,"sqz_d", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_sqz_d))
     CALL netcdf_err(  nf90_put_att(ncid,id_sqz_d,'long_name',    &
          '* main ion orbit squeezing factor'))
     CALL netcdf_err(  nf90_put_att(ncid,id_sqz_d,'units',        &
          'dimensionless'))

     CALL netcdf_err( nf90_def_var(ncid,"Epsi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_epsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_epsi,'long_name',    &
          '* Er/RBp neoclassical'))
     CALL netcdf_err(  nf90_put_att(ncid,id_epsi,'units',        &
          '1/second'))

     CALL netcdf_err( nf90_def_var(ncid,"Epsi_exp", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_epsi_exp))
     CALL netcdf_err(  nf90_put_att(ncid,id_epsi_exp,'long_name',    &
          '* Er/RBp experimental'))
     CALL netcdf_err(  nf90_put_att(ncid,id_epsi_exp,'units',        &
          '1/second'))


     CALL netcdf_err( nf90_def_var(ncid,"cer_bp", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_cer_bp))
     CALL netcdf_err(  nf90_put_att(ncid,id_cer_bp,'long_name',    &
          '* Bp along Z=0 coord'))
     CALL netcdf_err(  nf90_put_att(ncid,id_cer_bp,'units',        &
          'Tesla'))


     CALL netcdf_err( nf90_def_var(ncid,"cer_btdr", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_cer_btdr))
     CALL netcdf_err(  nf90_put_att(ncid,id_cer_btdr,'long_name',    &
          '* Bt/R along Z=0 cord'))
     CALL netcdf_err(  nf90_put_att(ncid,id_cer_btdr,'units',        &
          'Tesla/meter'))

     CALL netcdf_err( nf90_def_var(ncid,"cer_btdr", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_cer_btdr))
     CALL netcdf_err(  nf90_put_att(ncid,id_cer_btdr,'long_name',    &
          '* Bt/R along Z=0 coord'))
     CALL netcdf_err(  nf90_put_att(ncid,id_cer_btdr,'units',        &
          'Tesla/meter'))

     CALL netcdf_err( nf90_def_var(ncid,"totcur_bc ", nf90_double,          &
          VARID=id_totcur_bc ))
     label = '* boundary condition total current,amps,at time '     &
          //bc_asc_time(1:LEN_TRIM(bc_asc_time))
     CALL netcdf_err(  nf90_put_att(ncid,id_totcur_bc ,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_totcur_bc ,'units',        &
          ' Amps'))


     CALL netcdf_err( nf90_def_var(ncid,"vloop_bc ", nf90_double,          &
          VARID=id_vloop_bc ))
     label = '* boundary condition loop voltage,Volts,at time '//bc_asc_time(1:LEN_TRIM(bc_asc_time))

     CALL netcdf_err(  nf90_put_att(ncid,id_vloop_bc ,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_vloop_bc ,'units',        &
          ' VOLTS'))




     CALL netcdf_err( nf90_def_var(ncid,"fix_edge_te_bc", nf90_double,          &
          VARID=id_fix_edge_te_bc))
     label ='* boundary condition rho flags at time '                           &
          //bc_asc_time(1:LEN_TRIM(bc_asc_time))
     CALL netcdf_err(  nf90_put_att(ncid,id_fix_edge_te_bc ,'long_name',label))

     CALL netcdf_err( nf90_def_var(ncid,"fix_edge_ti_bc ", nf90_double,          &
          VARID=id_fix_edge_ti_bc ))
     CALL netcdf_err(  nf90_put_att(ncid,id_fix_edge_ti_bc,'long_name',label))
     CALL netcdf_err( nf90_def_var(ncid,"fix_edge_rot_bc ", nf90_double,          &
          VARID=id_fix_edge_rot_bc ))
     CALL netcdf_err(  nf90_put_att(ncid,id_fix_edge_rot_bc,'long_name',label))

     CALL netcdf_err( nf90_def_var(ncid,"fix_edge_ni_bc ", nf90_double, &
          DIMIDS = (/idim_nion /),VARID=id_fix_edge_ni_bc ))
     CALL netcdf_err(  nf90_put_att(ncid,id_fix_edge_ni_bc,'long_name',label))




     CALL netcdf_err( nf90_def_var(ncid,"te_bc", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_te_bc ))
     label = '* boundary condition TE,kev,at time '//bc_asc_time(1:LEN_TRIM(bc_asc_time))
     CALL netcdf_err(  nf90_put_att(ncid,id_te_bc,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_te_bc ,'units',        &
          'Kev'))

     CALL netcdf_err( nf90_def_var(ncid,"ti_bc", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_ti_bc ))
     label = '* boundary condition TI,kev,at time '//bc_asc_time(1:LEN_TRIM(bc_asc_time))
     CALL netcdf_err(  nf90_put_att(ncid,id_ti_bc,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_ti_bc ,'units',        &
          'Kev'))

     CALL netcdf_err( nf90_def_var(ncid,"ene_bc", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_ene_bc ))
     label ='* bc profile: ene, 1/m**3 at time  '//bc_asc_time(1:LEN_TRIM(bc_asc_time))
     CALL netcdf_err(  nf90_put_att(ncid,id_ene_bc,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_ene_bc ,'units',        &
          '1./meter^3'))

     CALL netcdf_err( nf90_def_var(ncid,"zeff_bc", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_zeff_bc ))
     label = '* bc profile: zeff, at time  '//bc_asc_time(1:LEN_TRIM(bc_asc_time))
     CALL netcdf_err(  nf90_put_att(ncid,id_zeff_bc,'long_name',label))

     CALL netcdf_err( nf90_def_var(ncid,"angrot_bc", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_angrot_bc ))
     label = '*  angular rotation speed profile, rad/sec at time '        &
          //bc_asc_time(1:LEN_TRIM(bc_asc_time))
     CALL netcdf_err(  nf90_put_att(ncid,id_angrot_bc,'long_name',label))

     CALL netcdf_err( nf90_def_var(ncid,"wbeam", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_wbeam ))
     label = '* fast ion stored energy density KEV/m**3 '
     CALL netcdf_err(  nf90_put_att(ncid,id_wbeam,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_wbeam ,'units',        &
          'Kev/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid,"walp", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_walp))
     label = '* fast alpha stored energy density KEV/m**3 ' 
     CALL netcdf_err(  nf90_put_att(ncid,id_walp,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_walp ,'units',        &
          'Kev/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid,"enalp", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_enalp))
     label = '* fast alpha density 1/m**3'
     CALL netcdf_err(  nf90_put_att(ncid,id_enalp,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_enalp ,'units',        &
          '1/meter^3'))



     CALL netcdf_err( nf90_def_var(ncid,"eps", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_eps))
     label = '* horizontal inverse aspect ratio = (rmax-rmin)/(rmax+rmin)'
     CALL netcdf_err(  nf90_put_att(ncid,id_eps,'long_name',label))





     CALL netcdf_err( nf90_def_var(ncid,"rcap", nf90_double,            &
          DIMIDS = (/idim_rho/),VARID=id_rcap))
     label = '* rcap = < R>,m'
     CALL netcdf_err(  nf90_put_att(ncid,id_rcap,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_rcap ,'units',        &
          'meters'))



     CALL netcdf_err( nf90_def_var(ncid,"rcapi", nf90_double,            &
          DIMIDS = (/idim_rho/),VARID=id_rcapi))
     label = '* rcap i= <1/ R>,1/m'
     CALL netcdf_err(  nf90_put_att(ncid,id_rcapi,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_rcapi ,'units',        &
          '1/meters'))




     CALL netcdf_err( nf90_def_var(ncid,"r2cap", nf90_double,           &
          DIMIDS = (/idim_rho/),VARID=id_r2cap))
     label = '*r2cap = <R0**2/R**2>'
     CALL netcdf_err(  nf90_put_att(ncid,id_r2cap,'long_name',label))

     CALL netcdf_err( nf90_def_var(ncid,"r2capi", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_r2capi))
     label = '*r2capi = <R**2>       M**2'
     CALL netcdf_err(  nf90_put_att(ncid,id_r2capi,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_r2capi ,'units',            &
          'meters^2'))



     CALL netcdf_err( nf90_def_var(ncid,"xhm2", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_xhm2))
     label = '*xhm2 = < (B total/ B axis)**2 > (=1 for circular plasmas)'
     CALL netcdf_err(  nf90_put_att(ncid,id_xhm2,'long_name',label))

     CALL netcdf_err( nf90_def_var(ncid,"xi11", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_xi11))
     label ='*xi11 ( = 1.95 sqrt(eps)for circular plasmas)'
     CALL netcdf_err(  nf90_put_att(ncid,id_xi11,'long_name',label))


     CALL netcdf_err( nf90_def_var(ncid,"xi33", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_xi33))
     label ='*xi33 ( = 1.95 sqrt(eps)for circular plasmas)'
     CALL netcdf_err(  nf90_put_att(ncid,id_xi33,'long_name',label))

     CALL netcdf_err( nf90_def_var(ncid,"xips", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_xips))
     label = '*xips = <(Baxis/B)**2)> - 1./(<(B/Baxis)**2> )( = 2 eps**2 for circular plasmas)'
     CALL netcdf_err(  nf90_put_att(ncid,id_xips,'long_name',label))




     CALL netcdf_err( nf90_def_var(ncid,"dfdt", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_dfdt))
     label = '(d/dt)Fcap '
     CALL netcdf_err(  nf90_put_att(ncid,id_dfdt,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_dfdt,'units',        &
          '1/sec'))



     CALL netcdf_err( nf90_def_var(ncid,"dgdt", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_dgdt))
     label = '(d/dt)Gcap '
     CALL netcdf_err(  nf90_put_att(ncid,id_dgdt,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_dgdt ,'units',        &
          '1/sec'))



     CALL netcdf_err( nf90_def_var(ncid,"dhdt", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_dhdt))
     label = '(d/dt)Hcap '
     CALL netcdf_err(  nf90_put_att(ncid,id_dhdt,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_dhdt,'units',        &
          '1/sec'))








     !         ---------------------------------------------------------------

     CALL netcdf_err(nf90_enddef(ncid))   ! leave netcdf define mode

     !         ---------------------------------------------------------------






     !         put values of variables

     profile%te0     = get_element(profile%te,1)
     profile%ti0     = get_element(profile%ti,1)



     CALL netcdf_err( nf90_put_var(ncid,id_shot,shot_id%shot_nmbr),id_shot)

     CALL netcdf_err( nf90_put_var(ncid,id_time,time),id_time)

     CALL netcdf_err( nf90_put_var(ncid,id_tgcnmf,tGCNMF),id_tgcnmf)

     CALL netcdf_err( nf90_put_var(ncid,id_time_bc,time_bc),id_time_bc)

     CALL netcdf_err( nf90_put_var(ncid,id_nj,nj),id_nj)

     CALL netcdf_err(nf90_put_var(ncid,id_psiaxis,mhd_dat%psiaxis),id_psiaxis)

     CALL netcdf_err(nf90_put_var(ncid,id_psibdry,mhd_dat%psibdry),id_psibdry)

     CALL netcdf_err(nf90_put_var(ncid,id_rgeom, dischg%rgeom),id_rgeom)

     CALL netcdf_err(nf90_put_var(ncid,id_btgeom, dischg%btgeom),id_btgeom)

     CALL netcdf_err(nf90_put_var(ncid,id_rmag, dischg%rma),id_rmag)

     CALL netcdf_err(nf90_put_var(ncid,id_zma, dischg%zma),id_zma)

     CALL netcdf_err(nf90_put_var(ncid,id_rmajor, dischg%rmajor),id_rmajor)

     CALL netcdf_err(nf90_put_var(ncid,id_rsep, dischg%rsep),id_rsep)
     CALL netcdf_err(nf90_put_var(ncid,id_zsep, dischg%zsep),id_zsep)

     CALL netcdf_err(nf90_put_var(ncid,id_rplasmin, dischg%rplasmin),id_rplasmin)
     CALL netcdf_err(nf90_put_var(ncid,id_rplasmax, dischg%rplasmax),id_rplasmax)
     CALL netcdf_err(nf90_put_var(ncid,id_zplasmin, dischg%zplasmin),id_zplasmin)
     CALL netcdf_err(nf90_put_var(ncid,id_zplasmax, dischg%zplasmax),id_zplasmax)



     CALL netcdf_err(nf90_put_var(ncid,id_kappa, dischg%kappa),id_kappa)

     CALL netcdf_err(nf90_put_var(ncid,id_deltao, dischg%deltao),id_deltao)

     CALL netcdf_err(nf90_put_var(ncid,id_pindento, dischg%pindento),id_pindento)

     CALL netcdf_err(nf90_put_var(ncid,id_volume, dischg%volo),id_volume)

     CALL netcdf_err(nf90_put_var(ncid,id_circum, dischg%circum),id_circum)

     CALL netcdf_err(nf90_put_var(ncid,id_areao, dischg%areao),id_areao)

     CALL netcdf_err( nf90_put_var(ncid,id_nion,nion),id_nion)

     CALL netcdf_err( nf90_put_var(ncid,id_nprim,nprim),id_nprim)

     CALL netcdf_err( nf90_put_var(ncid,id_fd_thermal,fd_thermal),id_fd_thermal)

     CALL netcdf_err( nf90_put_var(ncid,id_nimp,nimp),id_nimp)

     CALL netcdf_err( nf90_put_var(ncid,id_nneu,nneu),id_nneu)

     CALL netcdf_err( nf90_put_var(ncid,id_nbion,nbion),id_nbion)

     CALL netcdf_err( nf90_put_var(ncid,id_fd_beam,fd_beam),id_fd_beam)

     CALL netcdf_err( nf90_put_var(ncid,id_namep,namep),id_namep)

     CALL netcdf_err( nf90_put_var(ncid,id_namei,namei),id_namei)

     CALL netcdf_err( nf90_put_var(ncid,id_namen,namen ),id_namen)

     CALL netcdf_err( nf90_put_var(ncid,id_nameb,nameb ),id_nameb)

     CALL netcdf_err( nf90_put_var(ncid,id_namepel,pellet%name),id_namepel)

     CALL netcdf_err( nf90_put_var(ncid,id_btor,   mhd_dat%btor ),id_btor)
     CALL netcdf_err( nf90_put_var(ncid,id_tot_cur,mhd_dat%tot_cur ),id_tot_cur)
     CALL netcdf_err( nf90_put_var(ncid,id_totohm_cur,mhd_dat%totohm_cur ),id_totohm_cur)
     CALL netcdf_err( nf90_put_var(ncid,id_totboot_cur,mhd_dat%totboot_cur ),id_totboot_cur)
     CALL netcdf_err( nf90_put_var(ncid,id_totbeam_cur,mhd_dat%totbeam_cur ),id_totbeam_cur)
     CALL netcdf_err( nf90_put_var(ncid,id_totrf_cur,  mhd_dat%totrf_cur ),id_totrf_cur)
     CALL netcdf_err( nf90_put_var(ncid,id_betap,  mhd_dat%betap ),id_betap)
     CALL netcdf_err( nf90_put_var(ncid,id_beta,  mhd_dat%beta ),id_beta)
     CALL netcdf_err( nf90_put_var(ncid,id_ali,    mhd_dat%ali ),id_ali)

     CALL netcdf_err( nf90_put_var(ncid,id_te0,profile%te0 ),id_te0)
     CALL netcdf_err( nf90_put_var(ncid,id_ti0,profile%ti0 ),id_ti0)







     ! ---     check that psir_grid is monotonic; it may not be in certain cases where
     ! ---     the current profile was evolved
     !

     work_nj(1:nj) = get_values(psir_grid)
     CALL check_monotonic (work_nj, nj, monotonic, 1)
     IF (.NOT. monotonic) THEN    ! psir_grid is not monotonic
        WRITE  (nlog, 5) time
        WRITE  (ncrt, 5) time
5       FORMAT (' the psi grid, calculated from the poloidal' / &
             ' B field evolution, is not monotonic'        / &
             ' therefore the data at this time (', 1pe12.6, &
             ') was not calculated')
        go to 2000
     END IF

     CALL netcdf_err( nf90_put_var(ncid,id_psir_grid,work_nj ),id_psir_grid)

     CALL netcdf_err( nf90_put_var(ncid,id_psi,mhd_dat%psi ),id_psi)

     work_npsi = get_values(rho_mhd_gridnpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_rho_mhd_gridnpsi,work_npsi ),id_rho_mhd_gridnpsi)

     work_nr = get_values(dischg%rmhdgrid)
     CALL netcdf_err( nf90_put_var(ncid,id_rmhdgrid,work_nr ),id_rmhdgrid)

     work_nz = get_values(dischg%zmhdgrid)
     CALL netcdf_err( nf90_put_var(ncid,id_zmhdgrid,work_nz ),id_zmhdgrid)

     work_nj = get_values(rho_grid)
     CALL netcdf_err( nf90_put_var(ncid,id_rho_grid,work_nj ),id_rho_grid)

     work_nj(1:nj) = get_values(mhd_dat%fcap)
     CALL netcdf_err( nf90_put_var(ncid,id_fcap,work_nj ),id_fcap)

     work_nj(1:nj) = get_values(mhd_dat%gcap)
     CALL netcdf_err( nf90_put_var(ncid,id_gcap,work_nj ),id_gcap)

     work_nj(1:nj) = get_values(mhd_dat%hcap)
     CALL netcdf_err( nf90_put_var(ncid,id_hcap,work_nj ),id_hcap)

     work_nj(1:nj) = get_values(profile%te)
     CALL netcdf_err( nf90_put_var(ncid,id_te,work_nj ),id_te)

     work_nj(1:nj) = get_values(profile%ti)
     CALL netcdf_err( nf90_put_var(ncid,id_ti,work_nj ),id_ti)

     work_nj(1:nj) = get_values(profile%press)
     CALL netcdf_err( nf90_put_var(ncid,id_press,work_nj ),id_press)


     work_nj(1:nj) = get_values(profile%pressb)
     CALL netcdf_err( nf90_put_var(ncid,id_pressb,work_nj ),id_pressb)

     work_nj = get_values(mhd_dat%q_value)
     CALL netcdf_err( nf90_put_var(ncid,id_q_value,work_nj ),id_q_value)

     work_nj(1:nj) = get_values(profile%ene)
     CALL netcdf_err( nf90_put_var(ncid,id_ene,work_nj ),id_ene)



     DO j=1,nion
        work_nj_nion(1:nj,j) = get_values(profile%en(j))
     ENDDO
     !the species are ordered as 1..nprim,1..nimp
     !where nion = nprim+nimp
     CALL netcdf_err( nf90_put_var(ncid,id_enp,work_nj_nion ),id_enp)

     DO j=1,nimp
        work_nj_nion(1:nj,j) = get_values(profile%en(nprim+j))
     ENDDO
     !the species are ordered as 1..nprim,1..nimp
     !where nion = nprim+nimp
     CALL netcdf_err( nf90_put_var(ncid,id_eni,work_nj_nion ),id_eni)

     IF( ALLOCATED(work_nj_ntot))DEALLOCATE(work_nj_ntot)
     ALLOCATE(work_nj_ntot(nj,ntot))

     DO j=1,nion
        work_nj_nion(1:nj,j) = get_values(profile%flux(j))
        CALL netcdf_err( nf90_put_var(ncid,id_pflux,work_nj_nion),id_pflux)
     ENDDO
 
     work_nj(:) = get_values(profile%flux(nion+1)) 
     CALL netcdf_err( nf90_put_var(ncid,id_efluxe,work_nj),id_efluxe)

     work_nj(:) = get_values(profile%flux(nion+2)) 
     CALL netcdf_err( nf90_put_var(ncid,id_efluxi,work_nj),id_efluxi)

     work_nj(:) = get_values(profile%flux(nion+3)) 
     CALL netcdf_err( nf90_put_var(ncid,id_fdyflux,work_nj),id_fdyflux)


     work_nj(:) = get_values(profile%flux(nion+4)) 
     CALL netcdf_err( nf90_put_var(ncid,id_rotflux,work_nj),id_rotflux)

     work_nj(:) = tglf_p_output(:,1)
     CALL netcdf_err( nf90_put_var(ncid,id_tglf_p_fluxe,work_nj ),id_tglf_p_fluxe)
     work_nj(:) = tglf_p_output(:,2)
     CALL netcdf_err( nf90_put_var(ncid,id_tglf_p_fluxp,work_nj ),id_tglf_p_fluxp)
     work_nj(:) = tglf_p_output(:,3)
     CALL netcdf_err( nf90_put_var(ncid,id_tglf_p_fluxi,work_nj ),id_tglf_p_fluxi)


     work_nj(:) = tglf_e_output(:,1)
     CALL netcdf_err( nf90_put_var(ncid,id_tglf_e_fluxe,work_nj ),id_tglf_e_fluxe)
     work_nj(:) = tglf_e_output(:,2)
     CALL netcdf_err( nf90_put_var(ncid,id_tglf_e_fluxp,work_nj ),id_tglf_e_fluxp)
     work_nj(:) = tglf_e_output(:,3)
     CALL netcdf_err( nf90_put_var(ncid,id_tglf_e_fluxi,work_nj ),id_tglf_e_fluxi)


 

     CALL netcdf_err( nf90_put_var(ncid,id_enn,enn ),id_enn)

     CALL netcdf_err( nf90_put_var(ncid,id_ennw,ennw ),id_ennw)

     CALL netcdf_err( nf90_put_var(ncid,id_ennv,ennv ),id_ennv)

     CALL netcdf_err( nf90_put_var(ncid,id_volsn,volsn ),id_volsn)

     work_nj = get_values(prtcl_src%stfuse) 
     CALL netcdf_err( nf90_put_var(ncid,id_stfuse,work_nj),id_stfuse)

     work_nj = get_values(prtcl_src%sbfuse) 
     CALL netcdf_err( nf90_put_var(ncid,id_sbfuse,work_nj),id_sbfuse)

     work_nj = get_values(prtcl_src%spellet) 
     CALL netcdf_err( nf90_put_var(ncid,id_spellet,work_nj),id_spellet)

     CALL netcdf_err( nf90_put_var(ncid,id_sbeame,sbeame),id_sbeame)

     CALL netcdf_err( nf90_put_var(ncid,id_sbeam,sbeam),id_sbeam)

     CALL netcdf_err( nf90_put_var(ncid,id_sion,sion),id_sion)

!     CALL netcdf_err( nf90_put_var(ncid,id_srecom,srecom),id_srecom)
     DO jj=1,nion
        work_nj_nion(:,jj) = get_values(prtcl_src%srecom(jj))
     ENDDO
     CALL netcdf_err( nf90_put_var(ncid,id_srecom,work_nj_nion),id_srecom)


     CALL netcdf_err( nf90_put_var(ncid,id_scx,scx),id_scx)

     CALL netcdf_err( nf90_put_var(ncid,id_sbcx,sbcx),id_sbcx)

     CALL netcdf_err( nf90_put_var(ncid,id_stsource,stsource),id_stsource)

     CALL netcdf_err( nf90_put_var(ncid,id_dudtsv,dudtsv),id_dudtsv)

     CALL netcdf_err( nf90_put_var(ncid,id_enbeam,enbeam),id_enbeam)

     work_nj(1:nj) = get_values(mhd_dat%curden)
     CALL netcdf_err( nf90_put_var(ncid,id_curden,work_nj),id_curden)

     work_nj(1:nj) = get_values(mhd_dat%curohm)
     CALL netcdf_err( nf90_put_var(ncid,id_curohm,work_nj),id_curohm)

     work_nj(1:nj) = get_values(mhd_dat%curboot)
     CALL netcdf_err( nf90_put_var(ncid,id_curboot,work_nj),id_curboot)

     CALL netcdf_err( nf90_put_var(ncid,id_curbeam,curbeam),id_curbeam)

     CALL netcdf_err( nf90_put_var(ncid,id_currf,currf),id_currf)

     work_nj(1:nj) = get_values(profile%etor)
     CALL netcdf_err( nf90_put_var(ncid,id_etor,work_nj),id_etor)

     work_nj(1:nj) = get_values(mhd_dat%rbp)
     CALL netcdf_err( nf90_put_var(ncid,id_rbp,work_nj),id_rbp)

     work_npsi(1:mhd_dat%npsi) = get_values(mhd_dat%ravgnpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_ravgnpsi,work_npsi),id_ravgnpsi)

     work_npsi(1:mhd_dat%npsi) = get_values(mhd_dat%ravginpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_ravginpsi,work_npsi),id_ravginpsi)

     work_npsi(1:mhd_dat%npsi) = get_values(mhd_dat%psivalnpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_psivalnpsi,work_npsi),id_psivalnpsi)

     work_npsi(1:mhd_dat%npsi) = get_values(mhd_dat%fpsinpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_fpsinpsi,work_npsi),id_fpsinpsi)

     work_nj(1:nj) = get_values(mhd_dat%pprim)
     CALL netcdf_err( nf90_put_var(ncid,id_pprim,work_nj),id_pprim)

     work_nj(1:nj) = get_values(mhd_dat%ffprim)
     CALL netcdf_err( nf90_put_var(ncid,id_ffprim,work_nj),id_ffprim)

     work_nj(1:nj) = get_values(mhd_dat%bp)
     CALL netcdf_err( nf90_put_var(ncid,id_bp,work_nj),id_bp)

     work_nj(1:nj) = get_values(mhd_dat%bprmaj)
     CALL netcdf_err( nf90_put_var(ncid,id_bprmaj,work_nj),id_bprmaj)

     work_nj(1:nj) = get_values(mhd_dat%btotrmaj)
     CALL netcdf_err( nf90_put_var(ncid,id_btotrmaj,work_nj),id_btotrmaj)

     work_nj(1:nj) = get_values(profile%angrot)
     CALL netcdf_err( nf90_put_var(ncid,id_angrot,work_nj),id_angrot)

     CALL netcdf_err( nf90_put_var(ncid,id_zeff,zeff),id_zeff)

     CALL netcdf_err( nf90_put_var(ncid,id_d,diffuse%dcoef),id_d)

     work_nj(1:nj) = get_values(diffuse%chieinv)
     CALL netcdf_err( nf90_put_var(ncid,id_chieinv,work_nj),id_chieinv)


     work_nj(1:nj) = get_values(diffuse%chiinv)
     CALL netcdf_err( nf90_put_var(ncid,id_chiinv,work_nj),id_chiinv)

     work_nj(1:nj) = get_values(diffuse%xkineo)
     CALL netcdf_err( nf90_put_var(ncid,id_xkineo,work_nj),id_xkineo)

     work_nj(1:nj) = get_values(diffuse%xkeneo)
     CALL netcdf_err( nf90_put_var(ncid,id_xkeneo,work_nj),id_xkeneo)


     work_nj(1:nj) = get_values(wpdot%dpedt)
     CALL netcdf_err( nf90_put_var(ncid,id_dpedt,work_nj),id_dpedt)

     DO jj = 1,nion
        work_nj_nion(1:nj,jj) = get_values(wpdot%dpidt(jj))
     ENDDO
     CALL netcdf_err( nf90_put_var(ncid,id_dpidt,work_nj_nion),id_dpidt)

     work_nj(1:nj) = get_values(pwrden%qconde)
     CALL netcdf_err( nf90_put_var(ncid,id_qconde,work_nj),id_qconde)


     work_nj(1:nj) = get_values(pwrden%qcondi)
     CALL netcdf_err( nf90_put_var(ncid,id_qcondi,work_nj),id_qcondi)

     work_nj(1:nj) = get_values(pwrden%qconve)
     CALL netcdf_err( nf90_put_var(ncid,id_qconve,work_nj),id_qconve)


     work_nj(1:nj) = get_values(pwrden%qconvi)
     CALL netcdf_err( nf90_put_var(ncid,id_qconvi,work_nj),id_qconvi)


     work_nj(1:nj) = get_values(pwrden%qbeame)
     CALL netcdf_err( nf90_put_var(ncid,id_qbeame,work_nj),id_qbeame)

     work_nj(1:nj) = get_values(pwrden%qbeami)
     CALL netcdf_err( nf90_put_var(ncid,id_qbeami,work_nj),id_qbeami)


     work_nj(1:nj) = get_values(pwrden%qdelt)
     CALL netcdf_err( nf90_put_var(ncid,id_qdelt,work_nj),id_qdelt)

     work_nj(1:nj) = get_values(pwrden%qexch)
     CALL netcdf_err( nf90_put_var(ncid,id_qexch,work_nj),id_qexch)

     work_nj(1:nj) = get_values(pwrden%qrfe)
     CALL netcdf_err( nf90_put_var(ncid,id_qrfe,work_nj),id_qrfe)

     work_nj(1:nj) = get_values(pwrden%qrfi)
     CALL netcdf_err( nf90_put_var(ncid,id_qrfi,work_nj),id_qrfi)

     work_nj(1:nj) = get_values(pwrden%qione)
     CALL netcdf_err( nf90_put_var(ncid,id_qione,work_nj),id_qione)

     work_nj(1:nj) = get_values(pwrden%qioni)
     CALL netcdf_err( nf90_put_var(ncid,id_qioni,work_nj),id_qioni)

     work_nj(1:nj) = get_values(pwrden%qcx)
     CALL netcdf_err( nf90_put_var(ncid,id_qcx,work_nj),id_qcx)

     work_nj(1:nj) = get_values(pwrden%qe2d)
     CALL netcdf_err( nf90_put_var(ncid,id_qe2d,work_nj),id_qe2d)
     work_nj(1:nj) = get_values(pwrden%qi2d)
     CALL netcdf_err( nf90_put_var(ncid,id_qi2d,work_nj),id_qi2d)


     work_nj(1:nj) = get_values(pwrden%qfuse)
     CALL netcdf_err( nf90_put_var(ncid,id_qfuse,work_nj),id_qfuse)
     work_nj(1:nj) = get_values(pwrden%qfusi)
     CALL netcdf_err( nf90_put_var(ncid,id_qfusi,work_nj),id_qfusi)

     !
     ! --- beam fusion electron ,ion heating profile
     ! --- fraction of beam fusion energy deposited on thermal electron and distributions:
     work_nj(1:nj) = get_values(pwrden%qbfuse)
     CALL netcdf_err( nf90_put_var(ncid,id_qbfuse,work_nj),id_qbfuse)
     work_nj(1:nj) = get_values(pwrden%qbfusi)
     CALL netcdf_err( nf90_put_var(ncid,id_qbfusi,work_nj),id_qbfusi)


     work_nj(1:nj) = get_values(pwrden%qmag)
     CALL netcdf_err( nf90_put_var(ncid,id_qmag,work_nj),id_qmag)

     work_nj(1:nj) = get_values(pwrden%qsawe)
     CALL netcdf_err( nf90_put_var(ncid,id_qsawe,work_nj),id_qsawe)
     work_nj(1:nj) = get_values(pwrden%qsawi)
     CALL netcdf_err( nf90_put_var(ncid,id_qsawi,work_nj),id_qsawi)

     work_nj(1:nj) = get_values(pwrden%qrad)
     CALL netcdf_err( nf90_put_var(ncid,id_qrad,work_nj),id_qrad)

     work_nj(1:nj) = get_values(pwrden%qohm)
     CALL netcdf_err( nf90_put_var(ncid,id_qohm,work_nj),id_qohm)

     work_npsi(1:mhd_dat%npsi) = get_values(dischg%rmajavnpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_rmajavnpsi,work_npsi),id_rmajavnpsi)

     work_npsi(1:mhd_dat%npsi) = get_values(dischg%rminavnpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_rminavnpsi,work_npsi),id_rminavnpsi)

     work_npsi(1:mhd_dat%npsi) = get_values(dischg%psivolpnpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_psivolpnpsi,work_npsi),id_psivolpnpsi)

     work_nj(1:nj) = get_values(dischg%psivolpnj)
     CALL netcdf_err( nf90_put_var(ncid,id_psivolp,work_nj),id_psivolp)


     work_npsi(1:mhd_dat%npsi) = get_values(dischg%elongxnpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_elongxnpsi,work_npsi),id_elongxnpsi)

     work_nj(1:nj) = get_values(dischg%elongxnj)
     CALL netcdf_err( nf90_put_var(ncid,id_elongx,work_nj),id_elongx)

     work_npsi(1:mhd_dat%npsi) = get_values(dischg%triangnpsi_u)
     CALL netcdf_err( nf90_put_var(ncid,id_triangnpsi_u,work_npsi),id_triangnpsi_u)

     work_npsi(1:mhd_dat%npsi) = get_values(dischg%triangnpsi_l)
     CALL netcdf_err( nf90_put_var(ncid,id_triangnpsi_l,work_npsi),id_triangnpsi_l)
      triangnpsi_lcl(:) = (dischg%triangnpsi_u%data(:) + dischg%triangnpsi_l%data(:))*0.5_DP
      tension = 0.0               ! don't use tension option of tspline
      tmax    = 0.0               ! max allowed tension
      bpar(1) = 0.0               ! set boundary conditions on spline
      bpar(2) = 0.0
      bpar(3) = 0.0
      bpar(4) = 0.0
      kpsi    = mhd_dat%npsi
      work_npsi(1:mhd_dat%npsi) = get_values(mhd_dat%psivalnpsi)
      CALL tspline (work_npsi,triangnpsi_lcl,kpsi,bpar,cs2spline,kpsi, &
                     ier,tension,aspline,bspline,cspline,dspline,       &
                     espline,fspline,tmax,psir_grid,work_nj,nj)

     work_npsi(1:mhd_dat%npsi) = get_values(dischg%pindentnpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_pindentnpsi,work_npsi),id_pindentnpsi)

     work_npsi(1:mhd_dat%npsi) = get_values(dischg%sfareanpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_sfareanpsi,work_npsi),id_sfareanpsi)

     work_npsi(1:mhd_dat%npsi) = get_values(dischg%cxareanpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_cxareanpsi,work_npsi),id_cxareanpsi)

     work_npsi(1:mhd_dat%npsi) = get_values(dischg%grho1npsi)
     CALL netcdf_err( nf90_put_var(ncid,id_grho1npsi,work_npsi),id_grho1npsi)
     work_npsi(1:mhd_dat%npsi) = get_values(dischg%grho2npsi)
     CALL netcdf_err( nf90_put_var(ncid,id_grho2npsi,work_npsi),id_grho2npsi)

     CALL netcdf_err( nf90_put_var(ncid,id_storqueb,storqueb),id_storqueb)

     CALL netcdf_err( nf90_put_var(ncid,id_kpol_c,kpol_c),id_kpol_c)
     CALL netcdf_err( nf90_put_var(ncid,id_kpol_d,kpol_d),id_kpol_d)
     CALL netcdf_err( nf90_put_var(ncid,id_kpol_exp,kpol_exp),id_kpol_exp)


     CALL netcdf_err( nf90_put_var(ncid,id_angrot_c,angrot_c),id_angrot_c)
     CALL netcdf_err( nf90_put_var(ncid,id_angrot_d,angrot_d),id_angrot_d)

     CALL netcdf_err( nf90_put_var(ncid,id_ave_vpar_d,ave_vpar_d),id_ave_vpar_d)
     CALL netcdf_err( nf90_put_var(ncid,id_udia_d,udia_d),id_udia_d)
    CALL netcdf_err( nf90_put_var(ncid,id_udia_c,udia_c),id_udia_c)
    CALL netcdf_err( nf90_put_var(ncid,id_ugrt,ugrt),id_ugrt)
    CALL netcdf_err( nf90_put_var(ncid,id_sqz_d,sqz_d),id_sqz_d)
    CALL netcdf_err( nf90_put_var(ncid,id_epsi,epsi),id_epsi)
    CALL netcdf_err( nf90_put_var(ncid,id_epsi_exp,epsi_exp),id_epsi_exp)
    CALL netcdf_err( nf90_put_var(ncid,id_cer_btdr,cer_btdr),id_cer_btdr)  
    CALL netcdf_err( nf90_put_var(ncid,id_cer_bp,cer_bp),id_cer_bp)
    CALL netcdf_err( nf90_put_var(ncid,id_cer_btdr,cer_btdr),id_cer_btdr)


     CALL netcdf_err( nf90_put_var(ncid,id_totcur_bc,totcur_bc),id_totcur_bc)

     CALL netcdf_err( nf90_put_var(ncid,id_vloop_bc,vloop_bc),id_vloop_bc)

     CALL netcdf_err( nf90_put_var(ncid,id_fix_edge_te_bc,fix_edge_te_bc),id_fix_edge_te_bc)
     CALL netcdf_err( nf90_put_var(ncid,id_fix_edge_ti_bc,fix_edge_ti_bc),id_fix_edge_ti_bc)
     CALL netcdf_err( nf90_put_var(ncid,id_fix_edge_rot_bc,fix_edge_rot_bc),id_fix_edge_rot_bc)
     CALL netcdf_err( nf90_put_var(ncid,id_fix_edge_ni_bc,fix_edge_ni_bc),id_fix_edge_ni_bc)

     CALL netcdf_err( nf90_put_var(ncid,id_te_bc,te_bc),id_te_bc)

     CALL netcdf_err( nf90_put_var(ncid,id_ti_bc,ti_bc),id_ti_bc)




     CALL netcdf_err( nf90_put_var(ncid,id_ene_bc,ene_bc),id_ene_bc)

     CALL netcdf_err( nf90_put_var(ncid,id_zeff_bc,zeff_bc),id_zeff_bc)

     CALL netcdf_err( nf90_put_var(ncid,id_angrot_bc,angrot_bc),id_angrot_bc)

     CALL netcdf_err( nf90_put_var(ncid,id_en_bc,en_bc),id_en_bc)


     CALL netcdf_err( nf90_put_var(ncid,id_flux_bc,flux_bc),id_flux_bc)

     CALL netcdf_err( nf90_put_var(ncid,id_z,z),id_z)
     CALL netcdf_err( nf90_put_var(ncid,id_zsq,zsq),id_zsq)

     CALL netcdf_err( nf90_put_var(ncid,id_wbeam,wbeam),id_wbeam)
     CALL netcdf_err( nf90_put_var(ncid,id_walp,walp),id_walp)
     CALL netcdf_err( nf90_put_var(ncid,id_enalp,enalp),id_enalp)

     work_nj(1:nj) = get_values(wpdot%dnedt)
     CALL netcdf_err( nf90_put_var(ncid,id_dnedt,work_nj),id_dnedt) 

     CALL netcdf_err( nf90_put_var(ncid,id_eps,eps),id_eps)
     work_nj(1:nj) = get_values(mhd_dat%rcap)
     CALL netcdf_err( nf90_put_var(ncid,id_rcap,work_nj),id_rcap)
     work_nj(1:nj) = get_values(mhd_dat%rcapi)
     CALL netcdf_err( nf90_put_var(ncid,id_rcapi,work_nj),id_rcapi)
     work_nj(1:nj) = get_values(mhd_dat%r2capi)
     CALL netcdf_err( nf90_put_var(ncid,id_r2capi,work_nj),id_r2capi)
     work_nj(1:nj) = get_values(mhd_dat%r2cap)
     CALL netcdf_err( nf90_put_var(ncid,id_r2cap,work_nj),id_r2cap)

     CALL netcdf_err( nf90_put_var(ncid,id_xi11,xi11),id_xi11)
     CALL netcdf_err( nf90_put_var(ncid,id_xi33,xi33),id_xi33)
     CALL netcdf_err( nf90_put_var(ncid,id_xips,xips),id_xips)
     CALL netcdf_err( nf90_put_var(ncid,id_xhm2,xhm2),id_xhm2)

     CALL netcdf_err( nf90_put_var(ncid,id_dfdt,dfdt),id_dfdt)
     CALL netcdf_err( nf90_put_var(ncid,id_dgdt,dgdt),id_dgdt)       
     CALL netcdf_err( nf90_put_var(ncid,id_dhdt,dhdt),id_dhdt)

     CALL netcdf_err( nf90_put_var(ncid,id_nplasbdry,dischg%nplasbdry),id_nplasbdry)
     work_nplasbdry(1:dischg%nplasbdry) = get_values(dischg%rplasbdry)
     CALL netcdf_err( nf90_put_var(ncid,id_rplasbdry,work_nplasbdry),id_rplasbdry)
     work_nplasbdry(1:dischg%nplasbdry) = get_values(dischg%zplasbdry)
     CALL netcdf_err( nf90_put_var(ncid,id_zplasbdry,work_nplasbdry),id_zplasbdry)


     work_nlimiter(1:dischg%nlimiter) = get_values(dischg%rlimiter)
     CALL netcdf_err( nf90_put_var(ncid,id_rlimiter,work_nlimiter),id_rlimiter)
     work_nlimiter(1:dischg%nlimiter) = get_values(dischg%zlimiter)
     CALL netcdf_err( nf90_put_var(ncid,id_zlimiter,work_nlimiter),id_zlimiter)


     !-------------------------------------------------------------------------
     !
     !     netcdf read file 
     ! ------------------------------------------------------------------------

  ELSE  rwncd
     !          READ netcdf file: Cant assume that dim ids are as set in this
     !           subroutine because netcdf file most likely is  foreign.



     CALL netcdf_err( nf90_inquire(ncid, nDimensions, nVariables, nAttributes, &
          unlimitedDimId))


     DO j=1,nDimensions  
        CALL netcdf_err( NF90_INQUIRE_DIMENSION(ncid, j, gen_name, idlen)) 
        CALL search_dims(gen_name,dim_names,nf_max_name,dim_size ,max_dims,index )
        label = "error in determination of required dimensions"
        Assign_dim : SELECT CASE(index)
        CASE(0) 
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
           lerrno = 33
           CALL  terminate(lerrno,nlog)
        CASE(1) 
           nj     = idlen
        CASE(2) 
           nion   = idlen
        CASE(3) 
           nprim  = idlen
        CASE(4) 
           IF(name_size .NE. idlen)THEN
              label ='error in character array declaration'
              WRITE(ncrt,FMT='(a)')label
              WRITE(nlog,FMT='(a)')label
              lerrno = 34
              CALL  terminate(lerrno,nlog)
           ENDIF
        CASE(5) 
           nimp  = idlen
        CASE(6) 
           nneu  = idlen
        CASE(7) 
           nbion = idlen
        CASE(8) 
           dischg%nplasbdry = idlen
        CASE(9) 
           ntot = idlen
        CASE(10)
           dischg%nr_mhd = idlen
        CASE(11)
           dischg%nz_mhd = idlen
        CASE(12)
           mhd_dat%npsi =  idlen
        CASE(13)
           dischg%nlimiter = idlen
        CASE DEFAULT
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
           lerrno = 33
           CALL  terminate(lerrno,nlog)
        END SELECT Assign_dim
     ENDDO




     IF(ALLOCATED(work_nj))DEALLOCATE(work_nj)
     IF(ALLOCATED(work_npsi))DEALLOCATE(work_npsi)
     IF(ALLOCATED(work_nplasbdry))DEALLOCATE(work_nplasbdry)
     IF(ALLOCATED(work_nlimiter))DEALLOCATE(work_nlimiter)
     IF(ALLOCATED(work_nr))DEALLOCATE(work_nr)
     IF(ALLOCATED(work_nz))DEALLOCATE(work_nz)


     IF(ASSOCIATED(mhd_dat%psi))DEALLOCATE(mhd_dat%psi)

     ALLOCATE(work_nj(nj))
     ALLOCATE(work_npsi(mhd_dat%npsi))
     ALLOCATE(work_nplasbdry(dischg%nplasbdry))
     ALLOCATE(work_nlimiter(dischg%nlimiter))
     ALLOCATE(work_nr(dischg%nr_mhd))
     ALLOCATE(work_nz(dischg%nz_mhd))


     ALLOCATE(mhd_dat%psi(dischg%nr_mhd,dischg%nz_mhd))

     !     in the future we may introduce time dependent metrics.
     !     then dfdt (= d/dt FCAP ,etc) will be calculated internally in this code
     !     by reading in fcap_bc, etc.
     !     in demo code just assume time independent for each
     !     slice that the solver is called.
     IF(.NOT. ALLOCATED(dfdt))THEN 
        ALLOCATE(dfdt(nj))
        dfdt(:) = zeroc
     ENDIF
     IF(.NOT. ALLOCATED(dgdt))THEN 
        ALLOCATE(dgdt(nj))
        dgdt(:) = zeroc
     ENDIF
     IF(.NOT. ALLOCATED(dhdt))THEN
        ALLOCATE(dhdt(nj))
        dhdt(:) = zeroc
     ENDIF


     label(:) = " "
     CALL netcdf_err(  &
          nf90_inquire_attribute(ncid,NF90_GLOBAL,"title", len = titlen))
     IF( titlen .LE. strlen) &
          CALL netcdf_err( nf90_get_att(ncid,NF90_GLOBAL,"title",title))

     CALL netcdf_err( nf90_inq_varid(ncid,"shot", id_shot)) !given name get id
     CALL netcdf_err( nf90_get_var(ncid,id_shot,shot_id%shot_nmbr)) !given id, get value(s) 
     CALL netcdf_err( nf90_inq_varid(ncid,"nion", id_nion),-id_nion) 
     CALL netcdf_err( nf90_get_var(ncid,id_nion,nion),-id_nion) 
     CALL netcdf_err( nf90_inq_varid(ncid,"nprim",id_nprim),-id_nprim) 
     CALL netcdf_err( nf90_get_var(ncid,id_nprim,nprim),-id_nprim) 
     CALL netcdf_err( nf90_inq_varid(ncid,"nprim",id_fd_thermal),-id_fd_thermal) 
     CALL netcdf_err( nf90_get_var(ncid,id_fd_thermal,fd_thermal),-id_fd_thermal)
     CALL netcdf_err( nf90_inq_varid(ncid,"nimp",id_nimp),-id_nimp)
     CALL netcdf_err( nf90_get_var(ncid,id_nimp,nimp),-id_nimp)
     CALL netcdf_err( nf90_inq_varid(ncid,"nneu",id_nneu),-id_nneu)
     CALL netcdf_err( nf90_get_var(ncid,id_nneu,nneu),-id_nneu)
     
     CALL netcdf_err( nf90_inq_varid(ncid,"ibion",id_nbion),-id_nbion)
     CALL netcdf_err( nf90_get_var(ncid,id_nbion,nbion),-id_nbion)
     CALL netcdf_err( nf90_inq_varid(ncid,"fd_beam",id_fd_beam),-id_fd_beam)
     CALL netcdf_err( nf90_get_var(ncid,id_fd_beam,fd_beam),-id_fd_beam)



     IF( .NOT. ASSOCIATED(namep))THEN
        ALLOCATE ( namep(nprim),namei(nimp),namen(nneu),nameb(nbion))
        ALLOCATE (profile%en(nion))   !2d array 
        ALLOCATE (profile%flux(ntot)) !2d array
     ENDIF


     ! nj,nion are now known  so allocate some arrays: 
     IF( ALLOCATED(work_nj_ntot))DEALLOCATE(work_nj_ntot)
     ALLOCATE(work_nj_ntot(nj,ntot))
     IF( ALLOCATED(work_nj_nion))DEALLOCATE(work_nj_nion)
     ALLOCATE(work_nj_nion(nj,nion))
     IF(.NOT. ALLOCATED(te_bc))ALLOCATE(te_bc(nj))
     IF(.NOT. ALLOCATED(ti_bc))ALLOCATE(ti_bc(nj))
     IF(.NOT. ALLOCATED(ene_bc))ALLOCATE(ene_bc(nj))
     IF(.NOT. ALLOCATED(zeff_bc))ALLOCATE(zeff_bc(nj))
     IF(.NOT. ALLOCATED(angrot_bc))ALLOCATE(angrot_bc(nj))
     IF(.NOT. ALLOCATED(eps))ALLOCATE(eps(nj))
     IF(.NOT. ALLOCATED(storqueb))ALLOCATE(storqueb(nj))
     IF(.NOT. ALLOCATED(wbeam))ALLOCATE(wbeam(nj))
     IF(.NOT. ALLOCATED(walp))ALLOCATE(walp(nj))
     IF(.NOT. ALLOCATED(enalp))ALLOCATE(enalp(nj))
!!$     IF(.NOT. ALLOCATED(dnedt))ALLOCATE(dnedt(nj))
!!$     IF(.NOT. ALLOCATED(dangrotdt))ALLOCATE(dangrotdt(nj))
     IF(.NOT. ASSOCIATED(zeff))ALLOCATE(zeff(nj))
     IF(.NOT. ASSOCIATED(z))ALLOCATE(z(nj,nion))
     IF(.NOT. ALLOCATED(en_bc))ALLOCATE(en_bc(nj,nion))
     IF(.NOT. ALLOCATED(flux_bc))ALLOCATE(flux_bc(nj,ntot))
!!$     IF(.NOT. ALLOCATED(dnidt))ALLOCATE(dnidt(nj,nion))
     IF(.NOT. ASSOCIATED(zsq))ALLOCATE(zsq(nj,nion))


     CALL netcdf_err( nf90_inq_varid(ncid,"namep",id_namep),-id_namep)
     CALL netcdf_err( nf90_get_var(ncid,id_namep,namep),-id_namep)
     CALL netcdf_err( nf90_inq_varid(ncid,"namei",id_namei))
     CALL netcdf_err( nf90_get_var(ncid,id_namei,namei),-id_namei)
     CALL netcdf_err( nf90_inq_varid(ncid,"namen",id_namen))
     CALL netcdf_err( nf90_get_var(ncid,id_namen,namen),-id_namen)
     CALL netcdf_err( nf90_inq_varid(ncid,"nameb",id_nameb))
     CALL netcdf_err( nf90_get_var(ncid,id_nameb,nameb),-id_nameb)
     CALL netcdf_err( nf90_inq_varid(ncid,"namepel",id_namepel))
     CALL netcdf_err( nf90_get_var(ncid,id_namepel,pellet%name),-id_namepel)

     CALL netcdf_err( nf90_inq_varid(ncid,"time", id_time),-id_time) 
     CALL netcdf_err( nf90_get_var(ncid,id_time,time),-id_time) 

     eqtime = time       !equilibirum time. It is assumed that the quatities
     !fcap,gcap,hcap,eps,xhm2,xi11,xi33,xips
     !were calculated at this time.

     CALL netcdf_err( nf90_inq_varid(ncid,"psiaxis", id_psiaxis),-id_psiaxis) 
     CALL netcdf_err( nf90_get_var(ncid,id_psiaxis,mhd_dat%psiaxis),-id_psiaxis) 

     CALL netcdf_err( nf90_inq_varid(ncid,"psibdry", id_psibdry),-id_psibdry) 
     CALL netcdf_err( nf90_get_var(ncid,id_psibdry,mhd_dat%psibdry),-id_psibdry)

     CALL netcdf_err( nf90_inq_varid(ncid,"rgeom", id_rgeom),-id_rgeom) 
     CALL netcdf_err( nf90_get_var(ncid,id_rgeom,dischg%rgeom),-id_rgeom) 

     CALL netcdf_err( nf90_inq_varid(ncid,"btgeom", id_btgeom),-id_btgeom) 
     CALL netcdf_err( nf90_get_var(ncid,id_btgeom,dischg%btgeom),-id_btgeom) 


     CALL netcdf_err( nf90_inq_varid(ncid,"rmag", id_rmag),-id_rmag) 
     CALL netcdf_err( nf90_get_var(ncid,id_rmag,dischg%rma),-id_rmag)    
 
     dischg%rmag = dischg%rma
     CALL netcdf_err( nf90_inq_varid(ncid,"zma", id_zma),-id_zma) 
     CALL netcdf_err( nf90_get_var(ncid,id_zma,dischg%zma),-id_zma) 


     CALL netcdf_err( nf90_inq_varid(ncid,"rsep", id_rsep),-id_rsep) 
     CALL netcdf_err( nf90_get_var(ncid,id_rsep,dischg%rsep),-id_rsep) 

     CALL netcdf_err( nf90_inq_varid(ncid,"zsep", id_zsep),-id_zsep) 
     CALL netcdf_err( nf90_get_var(ncid,id_zsep,dischg%zsep),-id_zsep) 

     CALL netcdf_err( nf90_inq_varid(ncid,"rmajor", id_rmajor),-id_rmajor) 
     CALL netcdf_err( nf90_get_var(ncid,id_rmajor,dischg%rmajor),-id_rmajor) 

     CALL netcdf_err( nf90_inq_varid(ncid,"rplasmin", id_rplasmin),-id_rplasmin) 
     CALL netcdf_err( nf90_get_var(ncid,id_rplasmin,dischg%rplasmin),-id_rplasmin) 
     CALL netcdf_err( nf90_inq_varid(ncid,"rplasmax", id_rplasmax),-id_rplasmax) 
     CALL netcdf_err( nf90_get_var(ncid,id_rplasmax,dischg%rplasmax),-id_rplasmax)

     CALL netcdf_err( nf90_inq_varid(ncid,"zplasmin", id_zplasmin),-id_zplasmin) 
     CALL netcdf_err( nf90_get_var(ncid,id_zplasmin,dischg%zplasmin),-id_zplasmin) 
     CALL netcdf_err( nf90_inq_varid(ncid,"zplasmax", id_zplasmax),-id_zplasmax) 
     CALL netcdf_err( nf90_get_var(ncid,id_zplasmax,dischg%zplasmax),-id_zplasmax)

     CALL netcdf_err( nf90_inq_varid(ncid,"kappa", id_kappa),-id_kappa) 
     CALL netcdf_err( nf90_get_var(ncid,id_kappa,dischg%kappa),-id_kappa) 
     CALL netcdf_err( nf90_inq_varid(ncid,"deltao", id_deltao),-id_deltao) 
     CALL netcdf_err( nf90_get_var(ncid,id_deltao,dischg%deltao),-id_deltao) 
     CALL netcdf_err( nf90_inq_varid(ncid,"pindento", id_pindento),-id_pindento) 
     CALL netcdf_err( nf90_get_var(ncid,id_pindento,dischg%pindento),-id_pindento) 
     CALL netcdf_err( nf90_inq_varid(ncid,"volo", id_volume),-id_volume) 
     CALL netcdf_err( nf90_get_var(ncid,id_volume,dischg%volo),-id_volume) 
     CALL netcdf_err( nf90_inq_varid(ncid,"circum", id_circum),-id_circum) 
     CALL netcdf_err( nf90_get_var(ncid,id_circum,dischg%circum),-id_circum) 
     CALL netcdf_err( nf90_inq_varid(ncid,"areao", id_areao),-id_areao) 
     CALL netcdf_err( nf90_get_var(ncid,id_areao,dischg%areao),-id_areao)
     CALL netcdf_err( nf90_inq_varid(ncid,"btor", id_btor),-id_btor) 
     CALL netcdf_err( nf90_get_var(ncid,id_btor,mhd_dat%btor),-id_btor)
     CALL netcdf_err( nf90_inq_varid(ncid,"totcur", id_tot_cur),-id_tot_cur) 
     CALL netcdf_err( nf90_get_var(ncid,id_tot_cur,mhd_dat%tot_cur),-id_tot_cur)
     CALL netcdf_err( nf90_inq_varid(ncid,"totohm", id_totohm_cur),-id_totohm_cur) 
     CALL netcdf_err( nf90_get_var(ncid,id_totohm_cur,mhd_dat%totohm_cur),-id_totohm_cur)
     CALL netcdf_err( nf90_inq_varid(ncid,"totboot", id_totboot_cur),-id_totboot_cur) 
     CALL netcdf_err( nf90_get_var(ncid,id_totboot_cur,mhd_dat%totboot_cur),-id_totboot_cur)
     CALL netcdf_err( nf90_inq_varid(ncid,"totbeam", id_totbeam_cur),-id_totbeam_cur) 
     CALL netcdf_err( nf90_get_var(ncid,id_totbeam_cur,mhd_dat%totbeam_cur),-id_totbeam_cur)
     CALL netcdf_err( nf90_inq_varid(ncid,"totrf", id_totrf_cur),-id_totrf_cur) 
     CALL netcdf_err( nf90_get_var(ncid,id_totrf_cur,mhd_dat%totrf_cur),-id_totrf_cur)
     ibcur =1_I4B ;irfc = 1_I4B
     IF(ABS(mhd_dat%totbeam_cur) .LT. 1.e-5)ibcur = 0
     IF(ABS(mhd_dat%totrf_cur)   .LT. 1.e-5)irfc  = 0

     CALL netcdf_err( nf90_inq_varid(ncid,"betap", id_betap),-id_betap) 
     CALL netcdf_err( nf90_get_var(ncid,id_betap,mhd_dat%betap),-id_betap)
     CALL netcdf_err( nf90_inq_varid(ncid,"beta", id_beta),-id_beta) 
     CALL netcdf_err( nf90_get_var(ncid,id_beta,mhd_dat%beta),-id_beta)
     CALL netcdf_err( nf90_inq_varid(ncid,"ali", id_ali),-id_ali) 
     CALL netcdf_err( nf90_get_var(ncid,id_ali,mhd_dat%ali),-id_ali)
     CALL netcdf_err( nf90_inq_varid(ncid,"te0", id_te0),-id_te0) 
     CALL netcdf_err( nf90_get_var(ncid,id_te0,profile%te0),-id_te0)
     CALL netcdf_err( nf90_inq_varid(ncid,"ti0", id_ti0),-id_ti0) 
     CALL netcdf_err( nf90_get_var(ncid,id_ti0,profile%ti0),-id_ti0)

     CALL netcdf_err( nf90_inq_varid(ncid,"psir", id_psir_grid),-id_psir_grid) 
     CALL netcdf_err( nf90_get_var(ncid,id_psir_grid,work_nj),-id_psir_grid)
     psir_grid = new_Vector(nj,work_nj) ! [volt*sec/rad]

     CALL check_monotonic( work_nj, nj, monotonic,one)


     IF (.NOT. monotonic) THEN    ! psir_grid is not monotonic
        WRITE  (ncrt, 17) time
17       FORMAT (' the psi grid,psir_grid, read in  from iterdb statefile' / &
             ' is not monotonic'        /                                    &
             ' therefore  we do not continue with this run',/,               &
             ' the time value is ', 1pe14.8)
        WRITE(nlog,17)time
        lerrno = 13_I4B
        CALL  terminate(lerrno,nlog)
     END IF



     CALL netcdf_err( nf90_inq_varid(ncid,"psi", id_psi),-id_psi) 
     CALL netcdf_err( nf90_get_var(ncid,id_psi,mhd_dat%psi),-id_psi) ! [volt*sec/rad]

     CALL netcdf_err( nf90_inq_varid(ncid,"rho_mhd_gridnpsi", id_rho_mhd_gridnpsi),-id_rho_mhd_gridnpsi) 
     CALL netcdf_err( nf90_get_var(ncid,id_rho_mhd_gridnpsi,work_npsi),-id_rho_mhd_gridnpsi)
     rho_mhd_gridnpsi = new_Vector(mhd_dat%npsi,work_npsi)


     CALL netcdf_err( nf90_inq_varid(ncid,"rmhdgrid", id_rmhdgrid),-id_rmhdgrid) 
     CALL netcdf_err( nf90_get_var(ncid,id_rmhdgrid,work_nr),-id_rmhdgrid)
     dischg%rmhdgrid = new_Vector(dischg%nr_mhd,work_nr)

     CALL netcdf_err( nf90_inq_varid(ncid,"zmhdgrid", id_zmhdgrid),-id_zmhdgrid) 
     CALL netcdf_err( nf90_get_var(ncid,id_zmhdgrid,work_nz),-id_zmhdgrid)
     dischg%zmhdgrid = new_Vector(dischg%nz_mhd,work_nz)

     CALL netcdf_err( nf90_inq_varid(ncid,"r", id_rho_grid),-id_rho_grid) 
     CALL netcdf_err( nf90_get_var(ncid,id_rho_grid,work_nj),-id_rho_grid)
     rho_grid = new_Vector(nj,work_nj)
     elmt = get_element(rho_grid,nj)
     elmt = 1.0_DP/elmt
     rho_gridn = real_mult_Vector (elmt,rho_grid)


     CALL netcdf_err( nf90_inq_varid(ncid,"fcap", id_fcap),-id_fcap) 
     CALL netcdf_err( nf90_get_var(ncid,id_fcap,work_nj), -id_fcap)
     mhd_dat%fcap = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"gcap", id_gcap),-id_gcap) 
     CALL netcdf_err( nf90_get_var(ncid,id_gcap,work_nj),-id_gcap)
     mhd_dat%gcap = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"hcap", id_hcap),-id_hcap) 
     CALL netcdf_err( nf90_get_var(ncid,id_hcap,work_nj),-id_hcap)
     mhd_dat%hcap = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"q", id_q_value),-id_q_value) 
     CALL netcdf_err( nf90_get_var(ncid,id_q_value,work_nj),-id_q_value)
     mhd_dat%q_value = new_Vector(nj,work_nj)


     CALL netcdf_err( nf90_inq_varid(ncid,"Te", id_te),-id_te) 
     CALL netcdf_err( nf90_get_var(ncid,id_te,work_nj),-id_te)
     profile%te = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"Ti", id_ti),-id_ti) 
     CALL netcdf_err( nf90_get_var(ncid,id_ti,work_nj),-id_ti)
     profile%ti = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"press", id_press),-id_press) 
     CALL netcdf_err( nf90_get_var(ncid,id_press,work_nj),-id_press)
     profile%press = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"pressb", id_pressb),-id_pressb) 
     CALL netcdf_err( nf90_get_var(ncid,id_pressb,work_nj),-id_pressb)
     profile%pressb = new_Vector(nj,work_nj)




     CALL netcdf_err( nf90_inq_varid(ncid,"ene", id_ene),-id_ene) 
     CALL netcdf_err( nf90_get_var(ncid,id_ene,work_nj),-id_ene)
     profile%ene = new_Vector(nj,work_nj)


     CALL netcdf_err( nf90_inq_varid(ncid,"enp", id_enp),-id_enp) 
     CALL netcdf_err( nf90_get_var(ncid,id_enp,work_nj_nion),-id_enp)
     !the species are ordered as 1..npriim,1..nimp
     !where nion = nprim+nimp
     DO jj=1,nprim
        profile%en(jj)  = new_Vector(nj,work_nj_nion(1,jj))
     ENDDO

     CALL netcdf_err( nf90_inq_varid(ncid,"eni", id_eni),-id_eni) 
     CALL netcdf_err( nf90_get_var(ncid,id_eni,work_nj_nion),-id_eni)
     !the species are ordered as 1..npriim,1..nimp
     !where nion = nprim+nimp
     DO jj=1,nimp
        profile%en(nprim+jj)  = new_Vector(nj,work_nj_nion(1,nprim+jj))
     ENDDO

     CALL netcdf_err( nf90_inq_varid(ncid,"p_flux", id_pflux),-id_pflux) 
     CALL netcdf_err( nf90_get_var(ncid,id_pflux,work_nj_nion),-id_pflux)
     DO jj=1,nion
        profile%flux(jj)  = new_Vector(nj,work_nj_nion(1,jj))
     ENDDO

     CALL netcdf_err( nf90_inq_varid(ncid,"e_fluxe", id_efluxe),-id_efluxe) 
     CALL netcdf_err( nf90_get_var(ncid,id_efluxe,work_nj),-id_efluxe)
     profile%flux(nion+1)  = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"e_fluxi", id_efluxi),-id_efluxi) 
     CALL netcdf_err( nf90_get_var(ncid,id_efluxi,work_nj),-id_efluxi)
     profile%flux(nion+2)  = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"fday_flux", id_fdyflux),-id_fdyflux) 
     CALL netcdf_err( nf90_get_var(ncid,id_fdyflux,work_nj),-id_fdyflux)
     profile%flux(nion+3)  = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"rot_flux", id_rotflux),-id_rotflux) 
     CALL netcdf_err( nf90_get_var(ncid,id_rotflux,work_nj),-id_rotflux)
     profile%flux(nion+4)  = new_Vector(nj,work_nj)



     IF(ALLOCATED(tglf_p_output))DEALLOCATE(tglf_p_output)
     ALLOCATE(tglf_p_output(nj,3))
     CALL netcdf_err( nf90_inq_varid(ncid,"tglf_elct_p_flux", id_tglf_p_fluxe),-id_tglf_p_fluxe) 
     CALL netcdf_err( nf90_get_var(ncid,id_tglf_p_fluxe,work_nj),-id_tglf_p_fluxe)
     tglf_p_output(:,1) = work_nj(:)  ! NOTE : to load tglf_p,e_flux* need to convert to zc grid
     CALL netcdf_err( nf90_inq_varid(ncid,"tglf_ion_p_flux",id_tglf_p_fluxp), id_tglf_p_fluxp) 
     CALL netcdf_err( nf90_get_var(ncid,id_tglf_p_fluxp,work_nj),-id_tglf_p_fluxp)
     tglf_p_output(:,2) = work_nj(:)
     CALL netcdf_err( nf90_inq_varid(ncid,"tglf_imp_p_flux", id_tglf_p_fluxi),-id_tglf_p_fluxi) 
     CALL netcdf_err( nf90_get_var(ncid,id_tglf_p_fluxi,work_nj),-id_tglf_p_fluxi)
     tglf_p_output(:,3) = work_nj(:)

     IF(ALLOCATED(tglf_e_output))DEALLOCATE(tglf_e_output)
     ALLOCATE(tglf_e_output(nj,3))
     CALL netcdf_err( nf90_inq_varid(ncid,"tglf_elc_e_flux", id_tglf_e_fluxe),-id_tglf_e_fluxe) 
     CALL netcdf_err( nf90_get_var(ncid,id_tglf_e_fluxe,work_nj),-id_tglf_e_fluxe)
     tglf_e_output(:,1) = work_nj(:)
     CALL netcdf_err( nf90_inq_varid(ncid,"tglf_ion_e_flux", id_tglf_e_fluxp),-id_tglf_e_fluxp) 
     CALL netcdf_err( nf90_get_var(ncid,id_tglf_e_fluxp,work_nj),-id_tglf_e_fluxp)
     tglf_e_output(:,2) = work_nj(:)
     CALL netcdf_err( nf90_inq_varid(ncid,"tglf_imp_e_flux", id_tglf_e_fluxi),-id_tglf_e_fluxi) 
     CALL netcdf_err( nf90_get_var(ncid,id_tglf_e_fluxi,work_nj),-id_tglf_e_fluxi)
     tglf_e_output(:,3) = work_nj(:)





     IF( .NOT. ALLOCATED(stsource))THEN
        ntot = nion+dp4
        ALLOCATE(stsource(nion,nj),sion(nj,nion),       &
             sbcx(nj,nion),scx(nj,nion),dudtsv(ntot,nj))
             sion(:,:) = zeroc ;scx(:,:) = zeroc
             sbcx(:,:) = zeroc  ; dudtsv(:,:) = zeroc
     ENDIF
     IF(.NOT. ASSOCIATED(prtcl_src%srecom))ALLOCATE(prtcl_src%srecom(nion))

     CALL netcdf_err( nf90_inq_varid(ncid,"sion", id_sion),-id_sion) 
     CALL netcdf_err( nf90_get_var(ncid,id_sion,sion),-id_sion)
     CALL netcdf_err( nf90_inq_varid(ncid,"srecom", id_srecom),-id_srecom) 
!     CALL netcdf_err( nf90_get_var(ncid,id_srecom,srecom),-id_srecom)
     CALL netcdf_err( nf90_get_var(ncid,id_srecom,work_nj_nion),-id_srecom)
     DO jj=1,nion
       prtcl_src%srecom(jj) = new_Vector(nj,work_nj_nion(1,jj))
     ENDDO
     CALL netcdf_err( nf90_inq_varid(ncid,"scx", id_scx),-id_scx) 
     CALL netcdf_err( nf90_get_var(ncid,id_scx,scx),-id_scx)
     CALL netcdf_err( nf90_inq_varid(ncid,"sbcx", id_sbcx),-id_sbcx) 
     CALL netcdf_err( nf90_get_var(ncid,id_sbcx,sbcx),-id_sbcx)
     CALL netcdf_err( nf90_inq_varid(ncid,"s", id_stsource),-id_stsource) 
     CALL netcdf_err( nf90_get_var(ncid,id_stsource,stsource),-id_stsource)
     CALL netcdf_err( nf90_inq_varid(ncid,"dudtsv", id_dudtsv),-id_dudtsv) 
     CALL netcdf_err( nf90_get_var(ncid,id_dudtsv,dudtsv),-id_dudtsv)
     CALL netcdf_err( nf90_inq_varid(ncid,"tGCNMf", id_tGCNMf),-id_tGCNMf) 
     CALL netcdf_err( nf90_get_var(ncid,id_tGCNMf,tGCNMf),-id_tGCNMf) 

     CALL netcdf_err( nf90_inq_varid(ncid,"time_bc", id_time_bc),-id_time_bc) 
     CALL netcdf_err( nf90_get_var(ncid,id_time_bc,time_bc),-id_time_bc) 


     CALL netcdf_err( nf90_inq_varid(ncid,"rmajavnpsi",id_rmajavnpsi),-id_rmajavnpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_rmajavnpsi,work_npsi),-id_rmajavnpsi)
     dischg%rmajavnpsi = new_Vector(mhd_dat%npsi,work_npsi)


     !
     ! --- fast ion density
     ! 
     IF( .NOT. ALLOCATED(enbeam)) ALLOCATE(enbeam(nj,nbion))
     IF( .NOT. ALLOCATED(enbeam_tot)) ALLOCATE(enbeam_tot(nj))
     CALL netcdf_err( nf90_inq_varid(ncid,"enbeam",id_enbeam),-id_enbeam)
     CALL netcdf_err( nf90_get_var(ncid,id_enbeam,enbeam),-id_enbeam)
     enbeam_tot(:) = zeroc

     DO jn =1,nbion
        enbeam_tot(:) = enbeam_tot(:) + enbeam(:,jn)
     ENDDO
     IF(.NOT. ALLOCATED(enn))ALLOCATE(enn(nj,nneu))
     CALL netcdf_err( nf90_inq_varid(ncid,"enn",id_enn),-id_enn)
     CALL netcdf_err( nf90_get_var(ncid,id_enn,enn),0)
     IF(.NOT. ALLOCATED(ennw))ALLOCATE(ennw(nj,nneu))
     CALL netcdf_err( nf90_inq_varid(ncid,"ennw",id_ennw))
     CALL netcdf_err( nf90_get_var(ncid,id_ennw,ennw),0)
     IF(.NOT. ALLOCATED(ennv))ALLOCATE(ennv(nj,nneu))
     CALL netcdf_err( nf90_inq_varid(ncid,"ennv",id_ennv))
     CALL netcdf_err( nf90_get_var(ncid,id_ennv,ennv),0)
     IF(.NOT. ALLOCATED(volsn))ALLOCATE(volsn(nj,nneu))
     volsn(:,:) = zeroc
     CALL netcdf_err( nf90_inq_varid(ncid,"volsn",id_volsn))
     CALL netcdf_err( nf90_get_var(ncid,id_volsn,volsn),0)

     CALL netcdf_err( nf90_inq_varid(ncid,"stfuse",id_stfuse),-id_stfuse)
     CALL netcdf_err( nf90_get_var(ncid,id_stfuse,work_nj),-id_stfuse)
     prtcl_src%stfuse  = new_vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"sbfuse",id_sbfuse),-id_sbfuse)
     CALL netcdf_err( nf90_get_var(ncid,id_sbfuse,work_nj),-id_sbfuse)
     prtcl_src%sbfuse  = new_vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"spellet",id_spellet),-id_spellet)
     CALL netcdf_err( nf90_get_var(ncid,id_spellet,work_nj),-id_spellet)
     prtcl_src%spellet  = new_vector(nj,work_nj)


     IF(.NOT. ALLOCATED(sbeame))ALLOCATE(sbeame(nj))
     sbeame = zeroc
     CALL netcdf_err( nf90_inq_varid(ncid,"sbion",id_sbeame),-id_sbeame)
     CALL netcdf_err( nf90_get_var(ncid,id_sbeame,sbeame),-id_sbeame)
     !
     ! --- thermal ion source due to beams
     !

     IF(.NOT. ALLOCATED(sbeam))ALLOCATE(sbeam(nj,nbion))
     sbeam(:,:) = zeroc
     CALL netcdf_err( nf90_inq_varid(ncid,"sbeam",id_sbeam),-id_sbeam)
     CALL netcdf_err( nf90_get_var(ncid,id_sbeam,sbeam),-id_sbeam)


     ! 
     ! ---  current profiles  :
     ! 

     CALL netcdf_err( nf90_inq_varid(ncid,"curden",id_curden),-id_curden)
     CALL netcdf_err( nf90_get_var(ncid,id_curden,work_nj),-id_curden)
     mhd_dat%curden = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"curohm",id_curohm),-id_curohm)
     CALL netcdf_err( nf90_get_var(ncid,id_curohm,work_nj),-id_curohm)
     mhd_dat%curohm = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"curboot",id_curboot),-id_curboot)
     CALL netcdf_err( nf90_get_var(ncid,id_curboot,work_nj),-id_curboot)
     mhd_dat%curboot = new_Vector(nj,work_nj)

     !beam current density
     !no mhd_dat%curbeam because beam current is fixed in gcnm code
     IF(.NOT. ALLOCATED(curbeam))THEN
        ALLOCATE(curbeam(nj))
        curbeam(:) = zeroc
     ENDIF
     CALL netcdf_err( nf90_inq_varid(ncid,"curdbeam",id_curbeam),-id_curbeam)
     CALL netcdf_err( nf90_get_var(ncid,id_curbeam,curbeam),-id_curbeam)

     ! RF current density
     ! no mhd_dat%currf  because rf  current is fixed in gcnm code
     IF(.NOT. ALLOCATED(currf))THEN
        ALLOCATE(currf(nj))
        currf(:) = zeroc
     ENDIF
     ! < Jrf dot B/Bt0>        
     CALL netcdf_err( nf90_inq_varid(ncid,"currf",id_currf),-id_currf)
     CALL netcdf_err( nf90_get_var(ncid,id_currf,currf),-id_currf)


     ! toroidal electric field < E dot B /Bt0 >
     CALL netcdf_err( nf90_inq_varid(ncid,"etor",id_etor),-id_etor)
     CALL netcdf_err( nf90_get_var(ncid,id_etor,work_nj),-id_etor)
     profile%etor = new_Vector(nj,work_nj)

     ! rho*bp0*fcap*gcap*hcap
     CALL netcdf_err( nf90_inq_varid(ncid,"rbp",id_rbp),-id_rbp)
     CALL netcdf_err( nf90_get_var(ncid,id_rbp,work_nj),-id_rbp)
     mhd_dat%rbp = new_Vector(nj,work_nj)

     !ravgnpsi (1d from edge to mag axis)
     CALL netcdf_err( nf90_inq_varid(ncid,"ravgnpsi",id_ravgnpsi),-id_ravgnpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_ravgnpsi,work_npsi),-id_ravgnpsi)
     mhd_dat%ravgnpsi = new_Vector(mhd_dat%npsi,work_npsi)

     !ravginpsi (1d from edge to mag axis)
     CALL netcdf_err( nf90_inq_varid(ncid,"ravginpsi",id_ravginpsi),-id_ravginpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_ravginpsi,work_npsi),-id_ravginpsi)
     mhd_dat%ravginpsi = new_Vector(mhd_dat%npsi,work_npsi)



     !psi grid (1d from edge to mag axis)
     CALL netcdf_err( nf90_inq_varid(ncid,"psivalnpsi",id_psivalnpsi),-id_psivalnpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_psivalnpsi,work_npsi),-id_psivalnpsi)
     mhd_dat%psivalnpsi = new_Vector(mhd_dat%npsi,work_npsi)

     !fpsi   on  psivalnpsi
     CALL netcdf_err( nf90_inq_varid(ncid,"fpsinpsi",id_fpsinpsi),-id_fpsinpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_fpsinpsi,work_npsi),-id_fpsinpsi)
     mhd_dat%fpsinpsi = new_Vector(mhd_dat%npsi,work_npsi)

     !pprim   on  rho grid
     CALL netcdf_err( nf90_inq_varid(ncid,"pprim",id_pprim),-id_pprim)
     CALL netcdf_err( nf90_get_var(ncid,id_pprim,work_nj),-id_pprim)
     mhd_dat%pprim = new_Vector(nj,work_nj)

     ! ffprim   on  rho grid
     CALL netcdf_err( nf90_inq_varid(ncid,"ffprim",id_ffprim),-id_ffprim)
     CALL netcdf_err( nf90_get_var(ncid,id_ffprim,work_nj),-id_ffprim)
     mhd_dat%ffprim = new_Vector(nj,work_nj)

     !<Bp>   on  rho grid
     CALL netcdf_err( nf90_inq_varid(ncid,"bp",id_bp),-id_bp)
     CALL netcdf_err( nf90_get_var(ncid,id_bp,work_nj),-id_bp)
     mhd_dat%bp = new_Vector(nj,work_nj)


     !bprmaj on rmaj corresponding to rho grid
     CALL netcdf_err( nf90_inq_varid(ncid,"bprmaj",id_bprmaj),-id_bprmaj)
     CALL netcdf_err( nf90_get_var(ncid,id_bprmaj,work_nj),-id_bprmaj)
     mhd_dat%bprmaj = new_Vector(nj,work_nj)

     !btotrmaj on rmaj corresponding to rho grid
     CALL netcdf_err( nf90_inq_varid(ncid,"btotrmaj",id_btotrmaj),-id_btotrmaj)
     CALL netcdf_err( nf90_get_var(ncid,id_btotrmaj,work_nj),-id_btotrmaj)
     mhd_dat%btotrmaj = new_Vector(nj,work_nj)

     !zeff profiles
     CALL netcdf_err( nf90_inq_varid(ncid,"zeff",id_zeff),-id_zeff)
     CALL netcdf_err( nf90_get_var(ncid,id_zeff,zeff),-id_zeff)

     ! angular rotation speed profile
     CALL netcdf_err( nf90_inq_varid(ncid,"angrot",id_angrot),-id_angrot)
     CALL netcdf_err( nf90_get_var(ncid,id_angrot,work_nj),-id_angrot)
     profile%angrot = new_Vector(nj,work_nj)


     ! diffusivity matrix 
     CALL netcdf_err( nf90_inq_varid(ncid,"d",id_d),-id_d)
     IF( ASSOCIATED (diffuse%dcoef))DEALLOCATE(diffuse%dcoef)
     ALLOCATE(diffuse%dcoef(ntot,ntot,nj))
     CALL netcdf_err( nf90_get_var(ncid,id_d,diffuse%dcoef),-id_d)


     ! thermal diff. profiles, electron and ion
     CALL netcdf_err( nf90_inq_varid(ncid,"chieinv",id_chieinv),-id_chieinv)
     CALL netcdf_err( nf90_get_var(ncid,id_chieinv,work_nj),-id_chieinv)
     diffuse%chieinv  = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"chiinv",id_chiinv),-id_chiinv)
     CALL netcdf_err( nf90_get_var(ncid,id_chiinv,work_nj),-id_chiinv)
     diffuse%chiinv  = new_Vector(nj,work_nj)

     ! --- ion neoclassical thermal conductivity
     CALL netcdf_err( nf90_inq_varid(ncid,"xkineo",id_xkineo),-id_xkineo)
     CALL netcdf_err( nf90_get_var(ncid,id_xkineo,work_nj),-id_xkineo)
     diffuse%xkineo  = new_Vector(nj,work_nj)

     ! --- electron  neoclassical thermal conductivity
     CALL netcdf_err( nf90_inq_varid(ncid,"xkeneo",id_xkeneo),-id_xkeneo)
     CALL netcdf_err( nf90_get_var(ncid,id_xkeneo,work_nj),-id_xkeneo)
     diffuse%xkeneo  = new_Vector(nj,work_nj)

     !---not  read in(these are pointers NOT vectors):
     !---diffuse%xkitot ,diffuse%xketot plus others

     !d(electron energy)/dt profile
     CALL netcdf_err( nf90_inq_varid(ncid,"dpedtc",id_dpedt),-id_dpedt)
     CALL netcdf_err( nf90_get_var(ncid,id_dpedt,work_nj),-id_dpedt)
     wpdot%dpedt = new_Vector(nj,work_nj)

     ! --- d(ion energy)/dt profile
     IF(.NOT. ASSOCIATED(wpdot%dpidt))ALLOCATE(wpdot%dpidt(nion))
     CALL netcdf_err( nf90_inq_varid(ncid,"dpidtc",id_dpidt),-id_dpidt)
     CALL netcdf_err( nf90_get_var(ncid,id_dpidt,work_nj_nion),-id_dpidt)
     DO jj = 1,nion
        wpdot%dpidt(jj)  = new_Vector(nj,work_nj_nion(:,jj))
     ENDDO

     ! --- electron conduction profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qconde",id_qconde),-id_qconde)
     CALL netcdf_err( nf90_get_var(ncid,id_qconde,work_nj),-id_qconde)
     pwrden%qconde = new_Vector(nj,work_nj)

     ! --- ion conduction profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qcondi",id_qcondi),-id_qcondi)
     CALL netcdf_err( nf90_get_var(ncid,id_qcondi,work_nj),-id_qcondi)
     pwrden%qcondi = new_Vector(nj,work_nj)

     ! --- electron convection profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qconve",id_qconve),-id_qconve)
     CALL netcdf_err( nf90_get_var(ncid,id_qconve,work_nj),-id_qconve)
     pwrden%qconve = new_Vector(nj,work_nj)

     ! --- ion convection profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qconvi",id_qconvi),-id_qconvi)
     CALL netcdf_err( nf90_get_var(ncid,id_qconvi,work_nj),-id_qconvi)
     pwrden%qconvi = new_Vector(nj,work_nj)

     ! --- beam electron profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qbeame",id_qbeame),-id_qbeame)
     CALL netcdf_err( nf90_get_var(ncid,id_qbeame,work_nj),-id_qbeame)
     pwrden%qbeame = new_Vector(nj,work_nj)

     ! --- electron ion equilibration profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qdelt",id_qdelt),-id_qdelt)
     CALL netcdf_err( nf90_get_var(ncid,id_qdelt,work_nj),-id_qdelt)
     pwrden%qdelt = new_Vector(nj,work_nj)

     ! --- beam ion profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qbeami",id_qbeami),-id_qbeami)
     CALL netcdf_err( nf90_get_var(ncid,id_qbeami,work_nj),-id_qbeami)
     pwrden%qbeami = new_Vector(nj,work_nj)

     ! --- anomalous electron ion energy exchange term,
     ! --- due for example to drift ballooning mode
     CALL netcdf_err( nf90_inq_varid(ncid,"qexch",id_qexch),-id_qexch)
     CALL netcdf_err( nf90_get_var(ncid,id_qexch,work_nj),-id_qexch)
     pwrden%qexch = new_Vector(nj,work_nj)


     ! --- RF electron heating profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qrfe",id_qrfe),-id_qrfe)
     CALL netcdf_err( nf90_get_var(ncid,id_qrfe,work_nj),-id_qrfe)
     pwrden%qrfe = new_Vector(nj,work_nj)

     ! --- RF ion heating profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qrfi",id_qrfi),-id_qrfi)
     CALL netcdf_err( nf90_get_var(ncid,id_qrfi,work_nj),-id_qrfi)
     pwrden%qrfi = new_Vector(nj,work_nj)

     ! --- qione heating profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qione",id_qione),-id_qione)
     CALL netcdf_err( nf90_get_var(ncid,id_qione,work_nj),-id_qione)
     pwrden%qione = new_Vector(nj,work_nj)

     ! --- qioni heating profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qioni",id_qioni),-id_qioni)
     CALL netcdf_err( nf90_get_var(ncid,id_qioni,work_nj),-id_qioni)
     pwrden%qioni = new_Vector(nj,work_nj)



     ! --- qxc, ion heating profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qcx",id_qcx),-id_qcx)
     CALL netcdf_err( nf90_get_var(ncid,id_qcx,work_nj),-id_qcx)
     pwrden%qcx   = new_Vector(nj,work_nj)

     ! --- 2d electron heating profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qe2d",id_qe2d),-id_qe2d)
     CALL netcdf_err( nf90_get_var(ncid,id_qe2d,work_nj),-id_qe2d)
     pwrden%qe2d   = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"qi2d",id_qi2d),-id_qi2d)
     CALL netcdf_err( nf90_get_var(ncid,id_qi2d,work_nj),-id_qi2d)
     pwrden%qi2d   = new_Vector(nj,work_nj)

     ! --- fusion electron heating  profile
     !      qfuse = qtfuse + qbfuse 
     !      (sum of thermal and fast ion progenated fusion heating to electrons
     !       fast part includes beam and alphas)
     !      (since qbfuse and qfuse are input, qtfuse is implicitely determined)
     CALL netcdf_err( nf90_inq_varid(ncid,"qfuse",id_qfuse),-id_qfuse)
     CALL netcdf_err( nf90_get_var(ncid,id_qfuse,work_nj),-id_qfuse)
     pwrden%qfuse  = new_Vector(nj,work_nj)

     ! --- fusion ion heating profile
     !      qfusi = qtfusi + qbfusi 
     !      (sum of thermal and fast ion progenated fusion heating to ions
     !       fast part includes beam and alphas)
     !      (since qbfusi and qfusi are input, qtfusi is implicitely determined)
     CALL netcdf_err( nf90_inq_varid(ncid,"qfusi",id_qfusi),-id_qfusi)
     CALL netcdf_err( nf90_get_var(ncid,id_qfusi,work_nj),-id_qfusi)
     pwrden%qfusi  = new_Vector(nj,work_nj)

     ! --- beam fusion electron heating profile
     ! --- (fraction of beam fusion energy deposited on 
     ! --- thermal electron distribution
     CALL netcdf_err( nf90_inq_varid(ncid,"qbfuse",id_qbfuse),-id_qbfuse)
     CALL netcdf_err( nf90_get_var(ncid,id_qbfuse,work_nj),-id_qbfuse)
     pwrden%qbfuse  = new_Vector(nj,work_nj)

     ! --- beam fusion ion heating profile
     ! --- (fraction of beam fusion energy deposited on thermal ion distribution)
     CALL netcdf_err( nf90_inq_varid(ncid,"qbfusi",id_qbfusi),-id_qbfusi)
     CALL netcdf_err( nf90_get_var(ncid,id_qbfusi,work_nj),-id_qbfusi)
     pwrden%qbfusi  = new_Vector(nj,work_nj)

     ! --- mag electron heating profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qmag",id_qmag),-id_qmag)
     CALL netcdf_err( nf90_get_var(ncid,id_qmag,work_nj),-id_qmag)
     pwrden%qmag  = new_Vector(nj,work_nj)

     ! --- sawtooth electron heating profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qsawe",id_qsawe),-id_qsawe)
     CALL netcdf_err( nf90_get_var(ncid,id_qsawe,work_nj),-id_qsawe)
     pwrden%qsawe  = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"qsawi",id_qsawi),-id_qsawi)
     CALL netcdf_err( nf90_get_var(ncid,id_qsawi,work_nj),-id_qsawi)
     pwrden%qsawi  = new_Vector(nj,work_nj)

     ! --- radiated power density
     CALL netcdf_err( nf90_inq_varid(ncid,"qrad",id_qrad),-id_qrad)
     CALL netcdf_err( nf90_get_var(ncid,id_qrad,work_nj),-id_qrad)
     pwrden%qrad  = new_Vector(nj,work_nj)

     ! --- qohm,ohmic heating profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qohm",id_qohm),-id_qohm)
     CALL netcdf_err( nf90_get_var(ncid,id_qohm,work_nj),-id_qohm)
     pwrden%qohm  = new_Vector(nj,work_nj)

     ! --- average minor radius
     CALL netcdf_err( nf90_inq_varid(ncid,"rminavnpsi",id_rminavnpsi),-id_rminavnpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_rminavnpsi,work_npsi),-id_rminavnpsi)
     dischg%rminavnpsi = new_Vector(mhd_dat%npsi,work_npsi)

     CALL netcdf_err( nf90_inq_varid(ncid,"psivolpnpsi",id_psivolpnpsi),-id_psivolpnpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_psivolpnpsi,work_npsi),-id_psivolpnpsi)
     dischg%psivolpnpsi = new_Vector(mhd_dat%npsi,work_npsi)

     CALL netcdf_err( nf90_inq_varid(ncid,"psivolp",id_psivolp),-id_psivolp)
     CALL netcdf_err( nf90_get_var(ncid,id_psivolp,work_nj),-id_psivolp)
     dischg%psivolpnj = new_Vector(nj,work_nj)



     ! ---  elongation
     CALL netcdf_err( nf90_inq_varid(ncid,"elongxnpsi",id_elongxnpsi),-id_elongxnpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_elongxnpsi,work_npsi),-id_elongxnpsi)
     dischg%elongxnpsi = new_Vector(mhd_dat%npsi,work_npsi)

     CALL netcdf_err( nf90_inq_varid(ncid,"elongx",id_elongx),-id_elongx)
     CALL netcdf_err( nf90_get_var(ncid,id_elongx,work_nj),-id_elongx)
     dischg%elongxnj  = new_Vector(nj,work_nj)



     ! --- triangularity
     CALL netcdf_err( nf90_inq_varid(ncid,"triangnpsi_u",id_triangnpsi_u),-id_triangnpsi_u)
     CALL netcdf_err( nf90_get_var(ncid,id_triangnpsi_u,work_npsi),-id_triangnpsi_u)
     dischg%triangnpsi_u = new_Vector(mhd_dat%npsi,work_npsi)
     CALL netcdf_err( nf90_inq_varid(ncid,"triangnpsi_l",id_triangnpsi_l),-id_triangnpsi_l)
     CALL netcdf_err( nf90_get_var(ncid,id_triangnpsi_l,work_npsi),-id_triangnpsi_l)
     dischg%triangnpsi_l = new_Vector(mhd_dat%npsi,work_npsi)
     CALL netcdf_err( nf90_inq_varid(ncid,"pindentnpsi",id_pindentnpsi),-id_pindentnpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_pindentnpsi,work_npsi),-id_pindentnpsi)
     dischg%pindentnpsi  = new_Vector(mhd_dat%npsi,work_npsi)

     ! --- surface area
     CALL netcdf_err( nf90_inq_varid(ncid,"sfareanpsi",id_sfareanpsi),-id_sfareanpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_sfareanpsi,work_npsi),-id_sfareanpsi)
     dischg%sfareanpsi = new_Vector(mhd_dat%npsi,work_npsi)

     ! --- cross-sectional area
     CALL netcdf_err( nf90_inq_varid(ncid,"cxareanpsi",id_cxareanpsi),-id_cxareanpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_cxareanpsi,work_npsi),-id_cxareanpsi)
     dischg%cxareanpsi = new_Vector(mhd_dat%npsi,work_npsi)

     ! --- flux surface average grad rho
     CALL netcdf_err( nf90_inq_varid(ncid,"grho1npsi",id_grho1npsi),-id_grho1npsi)
     CALL netcdf_err( nf90_get_var(ncid,id_grho1npsi,work_npsi),-id_grho1npsi)
     dischg%grho1npsi = new_Vector(mhd_dat%npsi,work_npsi)


     ! --- flux surface average (grad rho)**2
     CALL netcdf_err( nf90_inq_varid(ncid,"grho2npsi",id_grho2npsi),-id_grho2npsi)
     CALL netcdf_err( nf90_get_var(ncid,id_grho2npsi,work_npsi),-id_grho2npsi)
     dischg%grho2npsi = new_Vector(mhd_dat%npsi,work_npsi)

     ! --- beam torque density - new input added here so that backward
     !     compatibility is maintained
     CALL netcdf_err( nf90_inq_varid(ncid,"storqueb",id_storqueb),-id_storqueb)
     CALL netcdf_err( nf90_get_var(ncid,id_storqueb ,storqueb),-id_storqueb)

     CALL netcdf_err( nf90_inq_varid(ncid,"Kpol_c",id_kpol_c),-id_kpol_c)
     CALL netcdf_err( nf90_get_var(ncid,id_kpol_c,kpol_c),-id_kpol_c)
     CALL netcdf_err( nf90_inq_varid(ncid,"Kpol_d",id_kpol_d),-id_kpol_d)
     CALL netcdf_err( nf90_get_var(ncid,id_kpol_d,kpol_d),-id_kpol_d)
     CALL netcdf_err( nf90_inq_varid(ncid,"Kpol_exp",id_kpol_exp),-id_kpol_exp)
     CALL netcdf_err( nf90_get_var(ncid,id_kpol_exp,kpol_exp),-id_kpol_exp)

     CALL netcdf_err( nf90_inq_varid(ncid,"angrot_c",id_angrot_c),-id_angrot_c)
     CALL netcdf_err( nf90_get_var(ncid,id_angrot_c,angrot_c),-id_angrot_c)
     CALL netcdf_err( nf90_inq_varid(ncid,"angrot_d",id_angrot_d),-id_angrot_d)
     CALL netcdf_err( nf90_get_var(ncid,id_angrot_d,angrot_d),-id_angrot_d)
     CALL netcdf_err( nf90_inq_varid(ncid,"ave_vpar_d",id_ave_vpar_d),-id_ave_vpar_d)
     CALL netcdf_err( nf90_get_var(ncid,id_ave_vpar_d,ave_vpar_d),-id_ave_vpar_d)
     CALL netcdf_err( nf90_inq_varid(ncid,"udia_d",id_udia_d),-id_udia_d)
     CALL netcdf_err( nf90_get_var(ncid,id_udia_d,udia_d),-id_udia_d)

     CALL netcdf_err( nf90_inq_varid(ncid,"udia_c",id_udia_c),-id_udia_c)
     CALL netcdf_err( nf90_get_var(ncid,id_udia_c,udia_c),-id_udia_c)

     CALL netcdf_err( nf90_inq_varid(ncid,"ugrt",id_ugrt),-id_ugrt)
     CALL netcdf_err( nf90_get_var(ncid,id_ugrt,ugrt),-id_ugrt)

     CALL netcdf_err( nf90_inq_varid(ncid,"sqz_d",id_sqz_d),-id_sqz_d)
     CALL netcdf_err( nf90_get_var(ncid,id_sqz_d,ugrt),-id_sqz_d)

     CALL netcdf_err( nf90_inq_varid(ncid,"Epsi",id_epsi),-id_epsi)
     CALL netcdf_err( nf90_get_var(ncid,id_epsi,epsi),-id_epsi)

     CALL netcdf_err( nf90_inq_varid(ncid,"Epsi_exp",id_epsi_exp),-id_epsi_exp)
     CALL netcdf_err( nf90_get_var(ncid,id_epsi_exp,epsi_exp),-id_epsi_exp)

     CALL netcdf_err( nf90_inq_varid(ncid,"cer_bp",id_cer_bp),-id_cer_bp)
     CALL netcdf_err( nf90_get_var(ncid,id_cer_bp,cer_bp),-id_cer_bp)

     CALL netcdf_err( nf90_inq_varid(ncid,"cer_btdr",id_cer_btdr),-id_cer_btdr)
     CALL netcdf_err( nf90_get_var(ncid,id_cer_btdr,cer_btdr),-id_cer_btdr)

     !     in the future we may introduce time dependent metrics.
     !     then dfdt (= d/dt FCAP ,etc) will be calculated internally in this code
     !     by reading in fcap_bc, etc.
     !     for now just assume time independent for each
     !     slice that the solver is called.
     IF(.NOT. ALLOCATED(dfdt))ALLOCATE(dfdt(nj))
     IF(.NOT. ALLOCATED(dgdt))ALLOCATE(dgdt(nj))
     IF(.NOT. ALLOCATED(dhdt))ALLOCATE(dhdt(nj))

     !
     CALL netcdf_err( nf90_inq_varid(ncid,"dfdt",id_dfdt),-id_dfdt)
     rcode = nf90_get_var(ncid,id_dfdt ,dfdt )
     IF(rcode /= nf90_noerr)dfdt(:) = zeroc

     CALL netcdf_err( nf90_inq_varid(ncid,"dgdt",id_dgdt),-id_dgdt)
     rcode =  nf90_get_var(ncid,id_dgdt ,dgdt )
     IF(rcode /= nf90_noerr)dgdt(:) = zeroc

     CALL netcdf_err( nf90_inq_varid(ncid,"dhdt",id_dhdt),-id_dhdt)
     rcode =   nf90_get_var(ncid,id_dhdt ,dhdt )
     IF(rcode /= nf90_noerr)dhdt(:) = zeroc


     !**_bc means ** is for  boundary condtions 
     !we give the profile over the entire rho grid but most often
     !only the edge values are actually used to get the 
     !boundary conditions. SOme model such as glf23 may have
     !boundaries that are set as far in as rho 0.6 so the **_bc
     !profiles must have appropriate values at least in that range.

     !boundary condition total current :
     CALL netcdf_err( nf90_inq_varid(ncid,"totcur_bc",id_totcur_bc),-id_totcur_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_totcur_bc ,totcur_bc ),-id_totcur_bc)
     !boundary condition loop voltage :
     CALL netcdf_err( nf90_inq_varid(ncid,"vloop_bc",id_vloop_bc),-id_vloop_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_vloop_bc ,vloop_bc ),-id_vloop_bc)


     !boundary indicators (determine rho point at which bc is applied)
     CALL netcdf_err( nf90_inq_varid(ncid,"fix_edge_te_bc",id_fix_edge_te_bc),-id_fix_edge_te_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_fix_edge_te_bc ,fix_edge_te_bc ),-id_fix_edge_te_bc)

     CALL netcdf_err( nf90_inq_varid(ncid,"fix_edge_ti_bc",id_fix_edge_ti_bc),-id_fix_edge_ti_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_fix_edge_ti_bc ,fix_edge_ti_bc ),-id_fix_edge_ti_bc)

     CALL netcdf_err( nf90_inq_varid(ncid,"fix_edge_rot_bc",id_fix_edge_rot_bc),-id_fix_edge_rot_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_fix_edge_rot_bc ,fix_edge_rot_bc ),-id_fix_edge_rot_bc)

     IF(.NOT. ALLOCATED(fix_edge_ni_bc))ALLOCATE(fix_edge_ni_bc(nion))
     CALL netcdf_err( nf90_inq_varid(ncid,"fix_edge_ni_bc",id_fix_edge_ni_bc),-id_fix_edge_ni_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_fix_edge_ni_bc ,fix_edge_ni_bc ),-id_fix_edge_ni_bc)

     ! convert indecies input as grid points equivalen normalized rho values
     CALL convert_edge_indecies


     ! boundary condition TE,TI kev,ene 1/m**3
     CALL netcdf_err( nf90_inq_varid(ncid,"te_bc",id_te_bc),-id_te_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_te_bc,te_bc ),-id_te_bc)

     CALL netcdf_err( nf90_inq_varid(ncid,"ti_bc",id_ti_bc),-id_ti_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_ti_bc,ti_bc ),-id_ti_bc)


     CALL netcdf_err( nf90_inq_varid(ncid,"ene_bc",id_ene_bc),-id_ene_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_ene_bc,ene_bc ),-id_ene_bc)


     CALL netcdf_err( nf90_inq_varid(ncid,"zeff_bc",id_zeff_bc),-id_zeff_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_zeff_bc,zeff_bc ),-id_zeff_bc)

     ! boundary condition,toroidal rotation rad/sec:
     CALL netcdf_err( nf90_inq_varid(ncid,"angrot_bc",id_angrot_bc),-id_angrot_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_angrot_bc,angrot_bc ),-id_angrot_bc)

     ! boundary condition,ion densities,1/m**3:
     CALL netcdf_err( nf90_inq_varid(ncid,"en_bc",id_en_bc),-id_en_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_en_bc,en_bc ),-id_en_bc)

     !ion flux  boundary condition,1/(m^2 sec):
     CALL netcdf_err( nf90_inq_varid(ncid,"flux_bc",id_flux_bc),-id_flux_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_flux_bc,flux_bc ),-id_flux_bc)

     ! --- primary, impurity charge,charge square
     CALL netcdf_err( nf90_inq_varid(ncid,"z",id_z),-id_z)
     CALL netcdf_err( nf90_get_var(ncid,id_z,z ),-id_z)

     CALL netcdf_err( nf90_inq_varid(ncid,"zsq",id_zsq),-id_zsq)
     CALL netcdf_err( nf90_get_var(ncid,id_zsq,zsq ),-id_zsq)

     !fast ion stored energy density KEV/m**3 ' 
     CALL netcdf_err( nf90_inq_varid(ncid,"wbeam",id_wbeam),-id_wbeam)
     CALL netcdf_err( nf90_get_var(ncid,id_wbeam,wbeam ),-id_wbeam)

     !fast alpha particle stored energy density KEV/m**3 ' 
     CALL netcdf_err( nf90_inq_varid(ncid,"walp",id_walp),-id_walp)
     CALL netcdf_err( nf90_get_var(ncid,id_walp,walp ),-id_walp)

     !fast alpha density 1/m**3 ' 
     CALL netcdf_err( nf90_inq_varid(ncid,"enalp",id_enalp),-id_enalp)
     CALL netcdf_err( nf90_get_var(ncid,id_enalp,enalp ),-id_enalp)


     ! rate of change of ene  1/(m**3 sec)
     CALL netcdf_err( nf90_inq_varid(ncid,"dnedt",id_dnedt),-id_dnedt)
     CALL netcdf_err( nf90_get_var(ncid,id_dnedt,work_nj),-id_dnedt)
     wpdot%dnedt  =  new_Vector(nj,work_nj)


     !horizontal inverse aspect ratio = (rmax-rmin)/(rmax+rmin)
     CALL netcdf_err( nf90_inq_varid(ncid,"eps",id_eps),-id_eps)
     CALL netcdf_err( nf90_get_var(ncid,id_eps,eps),-id_eps)

     !rcap = < R>,m
     CALL netcdf_err( nf90_inq_varid(ncid,"rcap",id_rcap),-id_rcap)
     CALL netcdf_err( nf90_get_var(ncid,id_rcap,work_nj),-id_rcap)
     mhd_dat%rcap = new_Vector(nj,work_nj)


     !rcapi = < 1/R>,m
     CALL netcdf_err( nf90_inq_varid(ncid,"rcapi",id_rcapi),-id_rcapi)
     CALL netcdf_err( nf90_get_var(ncid,id_rcapi,work_nj),-id_rcapi)
     mhd_dat%rcapi = new_Vector(nj,work_nj)


     !r2cap  = <R0**2/R**2>
     CALL netcdf_err( nf90_inq_varid(ncid,"r2cap",id_r2cap),-id_r2cap)
     CALL netcdf_err( nf90_get_var(ncid,id_r2cap,work_nj),-id_r2cap)
     mhd_dat%r2cap = new_Vector(nj,work_nj)


     !r2capi = <R**2>       M**2'
     CALL netcdf_err( nf90_inq_varid(ncid,"r2capi",id_r2capi),-id_r2capi)
     CALL netcdf_err( nf90_get_var(ncid,id_r2capi,work_nj),-id_r2capi)
     mhd_dat%r2capi = new_Vector(nj,work_nj)

     !     some quantities required for trapped particle fraction, and
     !     neocl;assical transport

     IF(.NOT. ALLOCATED(xips))ALLOCATE(xips(nj))
     IF(.NOT. ALLOCATED(xips0))ALLOCATE(xips0(nj))
     IF(.NOT. ALLOCATED(xhm2))ALLOCATE(xhm2(nj))
     IF(.NOT. ALLOCATED(xhm20))ALLOCATE(xhm20(nj))
     IF(.NOT. ALLOCATED(xi11))ALLOCATE(xi11(nj))
     IF(.NOT. ALLOCATED(xi110))ALLOCATE(xi110(nj))
     IF(.NOT. ALLOCATED(xi33))ALLOCATE(xi33(nj))
     IF(.NOT. ALLOCATED(xi330))ALLOCATE(xi330(nj))


     !xhm2 = < (B total/ B axis)**2 > (=1 for circular plasmas)
     CALL netcdf_err( nf90_inq_varid(ncid,"xhm2",id_xhm2),-id_xhm2)
     CALL netcdf_err( nf90_get_var(ncid,id_xhm2,xhm2),-id_xhm2)

     !xi11  = 1.95 sqrt(eps)for circular plasmas (see doc for complete desc.)
     CALL netcdf_err( nf90_inq_varid(ncid,"xi11",id_xi11),-id_xi11)
     CALL netcdf_err( nf90_get_var(ncid,id_xi11,xi11),-id_xi11)

     !xi33 ( = 1.95 sqrt(eps)
     CALL netcdf_err( nf90_inq_varid(ncid,"xi33",id_xi33),-id_xi33)
     CALL netcdf_err( nf90_get_var(ncid,id_xi33,xi33),-id_xi33)

     !xips = <(Baxis/B)**2)> - 1./(<(B/Baxis)**2> )
     CALL netcdf_err( nf90_inq_varid(ncid,"xips",id_xips),-id_xips)
     CALL netcdf_err( nf90_get_var(ncid,id_xips,xips),-id_xips)

     xhm20(:) = xhm2(:)
     xi110(:) = xi11(:)
     xi330(:) = xi33(:)
     xips0(:) = xips(:)



     ! --- plasma boundary

      CALL netcdf_err( nf90_inq_varid(ncid,"nplasbdry",id_nplasbdry),-id_nplasbdry)
      CALL netcdf_err( nf90_get_var(ncid,id_nplasbdry,dischg%nplasbdry),-id_nplasbdry)
     CALL netcdf_err( nf90_inq_varid(ncid,"rplasbdry",id_rplasbdry),-id_rplasbdry)
     CALL netcdf_err( nf90_get_var(ncid,id_rplasbdry,work_nplasbdry),-id_rplasbdry)
     dischg%rplasbdry  = new_Vector(dischg%nplasbdry,work_nplasbdry)

     CALL netcdf_err( nf90_inq_varid(ncid,"zplasbdry",id_zplasbdry),-id_zplasbdry)
     CALL netcdf_err( nf90_get_var(ncid,id_zplasbdry,work_nplasbdry),-id_zplasbdry)
     dischg%zplasbdry  = new_Vector(dischg%nplasbdry,work_nplasbdry)


     ! ---  limiter for plasma:
     !          CALL netcdf_err( nf90_inq_varid(ncid,"nplasbdry",id_nlimiter))
     !          CALL netcdf_err( nf90_get_var(ncid,id_nplasbdry,dischg%nplasbdry))
     CALL netcdf_err( nf90_inq_varid(ncid,"rlimiter",id_rlimiter),-id_rlimiter)
     CALL netcdf_err( nf90_get_var(ncid,id_rlimiter,work_nlimiter),-id_rlimiter)
     dischg%rlimiter = new_Vector(dischg%nlimiter,work_nlimiter)

     CALL netcdf_err( nf90_inq_varid(ncid,"zlimiter",id_zlimiter),-id_zlimiter)
     CALL netcdf_err( nf90_get_var(ncid,id_zlimiter,work_nlimiter),-id_zlimiter)
     dischg%zlimiter  = new_Vector(dischg%nlimiter,work_nlimiter)

     ! -------------------------------------------------------------------------
     ! ---  END NETCDF file read/ write processing 
     ! -------------------------------------------------------------------------
  END IF  rwncd


  eqtime = time            !equilibirum time. It is assumed that the quatities
  !fcap,gcap,hcap,eps,xhm2,xi11,xi33,xips
  !were calcualted at this time. Hence subrotuine 
  !neointrp does nothing


2000 CONTINUE
  CALL netcdf_err(nf90_close(ncid) )         !close netcdf file

  IF(ALLOCATED(work_nj))DEALLOCATE(work_nj)
  IF(ALLOCATED(work_npsi))DEALLOCATE(work_npsi)
  IF(ALLOCATED(work_nplasbdry))DEALLOCATE(work_nplasbdry)
  IF(ALLOCATED(work_nlimiter))DEALLOCATE(work_nlimiter)
  IF(ALLOCATED(work_nr))DEALLOCATE(work_nr)
  IF(ALLOCATED(work_nz))DEALLOCATE(work_nz)
  IF(ALLOCATED(work_nj_nion))DEALLOCATE(work_nj_nion)
  IF(ALLOCATED(work_nj_ntot))DEALLOCATE(work_nj_ntot)
  RETURN
END SUBROUTINE iter_dbase_nc_compat


