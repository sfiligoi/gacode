!------------------------------------------------------
! gyro_write_timedata.F90
!
! PURPOSE:
!  This is the master file controlling output of
!  data on a per-timestep basis.  This file also 
!  contains the MPI IO routines 
!
!  - write_distributed_real
!  - write_distributed_complex
!  - write_local_real
!-----------------------------------------------------

subroutine gyro_write_timedata

  use gyro_globals
  use mpi
#ifdef HAVE_HDF5
  use hdf5_api
#endif

  !---------------------------------------------------
  implicit none
  !
  real, dimension(:,:,:), allocatable :: a3
  !
  complex, dimension(n_theta_plot,n_x,n_kinetic) :: n_plot
  complex, dimension(n_theta_plot,n_x,n_kinetic) :: e_plot
  complex, dimension(n_theta_plot,n_x,n_kinetic) :: v_plot
  !
  real :: pi=3.141592653589793
#ifdef HAVE_HDF5
  real, dimension(:), allocatable, save :: zeta_phi
  real, dimension(:,:), allocatable :: a2

  integer, parameter :: hr4=SELECTED_REAL_KIND(6,37)
  character(60) :: description
  character(64) :: step_name
  character(128) :: dumpfile
  character(128) :: meshfile
  character(20)   :: openmethod
  integer(HID_T) :: dumpGid,dumpFid,gid3D,fid3D
  integer(HID_T) :: dumpTGid,dumpTFid
  type(hdf5InOpts) :: h5in
  type(hdf5ErrorType) :: h5err
  integer :: number_label
  logical :: write_threed
  logical :: h5_rewind=.false.
#endif

  !---------------------------------------------------

  !---------------------------------------------------
  ! Timestep data:
  !
  if (io_method < 3) then
    if (i_proc == 0) then
       call gyro_write_step(trim(path)//'out.gyro.t',1)
    endif
  endif
  !---------------------------------------------------
#ifdef HAVE_HDF5
  if (io_method > 1 ) then
!---------------------------------------------------
      ! Determine if the 3D files need to be written 
      if (n_torangle_3d > 1 ) then
         write_threed = .true.
      else
         write_threed = .false.
      endif
     !---------------------------------------------------

      !---------------------------------------------------
      ! io_control: 0=no I/O, 1=open, 2=append, 3=rewind

      select case (io_control)
      case(0)
         return
      case(1)
         openmethod='overwr'
      case(2)
         openmethod='append'
      case(3)
         openmethod='append'
         h5_rewind = .true.
      end select
         
      if (i_proc == 0) then

         !---------------------------------------------------
         ! Initialization (not executed during timestepping)
         !---------------------------------------------------

         call vshdf5_inith5vars(h5in, h5err)
         h5in%comm=MPI_COMM_SELF
         h5in%info=MPI_INFO_NULL
         h5in%wrd_type=H5T_NATIVE_REAL
         !h5in%wrd_type=h5in%h5_kind_type_r4
         h5in%typeConvert=.true.
         h5in%doTranspose=.true.
         h5in%verbose=.true.
         h5in%debug=.false.
         h5in%wrVsTime=.true.
         h5in%vsTime=t_current
         h5in%vsStep=step


         if (debug_flag==1) h5in%debug=.true.

         !------------------------------------------------
         ! Open the monolithic timedata file (incremental)
         dumpfile=TRIM(path)//"out.gyro.timedata.h5" 
         description="Time-dependent GYRO data"
         call open_h5file(trim(openmethod),dumpfile,dumpTFid,description,dumpTGid,h5in,h5err)
         if (h5err%errBool) call catch_error(h5err%errorMsg)

          meshfile=TRIM(path)//"gyroMesh.h5"


            !---------------------------------------------------
            ! Call for each timestep
            !---------------------------------------------------

            number_label = nint(t_current/dt)
            ! Write number_label to step_name using write statement
            if (number_label > 999999) then
               write(step_name,fmt='(i7.7)') number_label
            else if (number_label > 99999) then
               write(step_name,fmt='(i6.6)') number_label
            else
               write(step_name,fmt='(i5.5)') number_label
            endif

            ! Open a new file at each data step     

            dumpfile    = trim(path)//"gyro"//trim(step_name)//".h5"
            description = "GYRO field file"
            call open_newh5file(dumpfile,dumpFid,description,dumpGid,h5in,h5err)

            if (write_threed) then
               dumpfile    = trim(path)//"gyro3D"//trim(step_name)//".h5"
               description = "GYRO 3D field file"
               call open_newh5file(dumpfile,fid3d,description,gid3D,h5in,h5err)
            endif

          ! make external links
          call make_external_link(TRIM(meshfile),"poloidalMesh", &
            dumpGid,"poloidalMesh", h5in,h5err)
          call make_external_link(TRIM(meshfile),"threeDMesh", &
            gid3D,"threeDMesh", h5in,h5err)

      endif !i_proc ==0
  endif ! io_method > 1
#endif
  !--------------------------------------------------
  ! Output of field-like quantities:
  !
  if (plot_n_flag+plot_e_flag+plot_v_flag > 0) then
     n_plot(:,:,:) = moments_plot(:,:,:,1)
     e_plot(:,:,:) = moments_plot(:,:,:,2)
     v_plot(:,:,:) = moments_plot(:,:,:,3)
  endif
  !
  if (plot_u_flag == 1) then

     ! POTENTIALS

    if(io_method < 3) then
       call write_distributed_complex(&
            trim(path)//'out.gyro.moment_u',&
            10,&
            size(phi_plot(:,:,1:n_field)),&
            phi_plot(:,:,1:n_field))
    endif
#ifdef HAVE_HDF5
    if(io_method >1 ) then
       h5in%units="dimensionless"
       h5in%mesh="/threeDMesh/cartMesh"
       call write_distributed_complex_sorf_h5("phi",&
            dumpGid,gid3D,&
            n_theta_plot*n_x*n_field,&
            n_theta_plot,n_x,n_field,&
            phi_plot(:,:,1:n_field),&
            write_threed,&
            .false., &
            h5in,h5err)
    endif
#endif
  endif !u_flag==1

  if (plot_epar_flag == 1) then

     ! PARALLEL ELECTRIC FIELD

    if(io_method < 3) then
     call write_distributed_complex(&
          trim(path)//'out.gyro.moment_epar',&
          10,&
          size(phi_plot(:,:,n_field+1)),&
          phi_plot(:,:,n_field+1))
    endif
#ifdef HAVE_HDF5
    if(io_method >1) then
       h5in%units="dimensionless"
       h5in%mesh="/cartMesh"
       call write_distributed_complex_sorf_h5("epar",&
            dumpGid,gid3D,&
            n_theta_plot*n_x*n_field,&
            n_theta_plot,n_x,n_field,&
            phi_plot(:,:,n_field+1),&
            write_threed,&
            .false., &
            h5in,h5err)
    endif
#endif
  endif !epar_flag==1

  if (plot_n_flag == 1) then

     ! DENSITY

    if(io_method < 3) then
     call write_distributed_complex(&
          trim(path)//'out.gyro.moment_n',&
          10,&
          size(n_plot),&
          n_plot)
    endif

#ifdef HAVE_HDF5
    if(io_method >1) then
     h5in%units="dimensionless"
     h5in%mesh="/cartMesh"
     call write_distributed_complex_sorf_h5("density",&
          dumpGid,gid3D,&
          n_theta_plot*n_x*n_kinetic,&
          n_theta_plot,n_x,n_kinetic,&
          moments_plot(:,:,:,1),&
          write_threed,&
          .false., &
          h5in,h5err)
    endif
#endif
  endif !n_flag ==1 

  if (plot_e_flag == 1) then

     ! ENERGY

    if(io_method < 3) then
    call write_distributed_complex(&
          trim(path)//'out.gyro.moment_e',&
          10,&
          size(e_plot),&
          e_plot)
    endif
#ifdef HAVE_HDF5
    if(io_method >1) then
    h5in%units="dimensionless"
    h5in%mesh="/cartMesh"
    call write_distributed_complex_sorf_h5("energy",&
          dumpGid,gid3D,&
          n_theta_plot*n_x*n_kinetic,&
          n_theta_plot,n_x,n_kinetic,&
          moments_plot(:,:,:,2),&
          write_threed,&
          .false., &
          h5in,h5err)
    endif
#endif
  endif !e_flag==1

  if (plot_v_flag == 1) then

     ! PARALLEL VELOCITY

    if(io_method < 3) then
     call write_distributed_complex(&
          trim(path)//'out.gyro.moment_v',&
          10,&
          size(v_plot),&
          v_plot)
    endif

#ifdef HAVE_HDF5
    if(io_method >1) then
     h5in%units="dimensionless"
     h5in%mesh="/cartMesh"
     call write_distributed_complex_sorf_h5("v_par",&
          dumpGid,gid3D,&
          n_theta_plot*n_x*n_kinetic,&
          n_theta_plot,n_x,n_kinetic,&
          moments_plot(:,:,:,3),&
          write_threed,&
          .false., &
          h5in,h5err)
    endif
#endif
  endif !v_flag==1

  !--------------------------------------------------

  !--------------------------------------------------
  ! Output of field at r=r0:
  !
  if (field_r0_flag == 1) then
    if(io_method < 3) then
     call write_distributed_complex(&
          trim(path)//'out.gyro.field_r0',&
          10,&
          size(field_r0_plot),&
          field_r0_plot)
    endif
!field_r0_grid,n_field  size of field_r0 (2d)
!    if(io_method >1) then
!     h5in%units="dimensionless"
!     h5in%mesh="/cartMesh"
!     call write_distributed_complex_2d_h5("field_r0",&
!          dumpGid,gid3D,&
!          n_theta_plot*n_x*n_kinetic,&
!          n_theta_plot,n_x,n_kinetic,&
!          field_r0_plot,&
!          write_threed,&
!          .false., &
!          h5in,h5err)
!    endif

  endif
  !--------------------------------------------------

  call gyro_kxky_spectrum


  if(io_method < 3) then
  call write_distributed_real(&
       trim(path)//'out.gyro.kxkyspec',&
       10,&
       size(kxkyspec),&
       kxkyspec)
  endif
#ifdef HAVE_HDF5
  if(io_method > 1) then
    h5in%units="dimensionless"
    h5in%mesh="/t_current"
    h5in%vsCentering="nodal"
    call write_distributed_real_h5("kxkyspec",dumpTGid,&
         n_x,1,1,&
         size(kxkyspec),&
         kxkyspec,&
         h5in,h5err)
    if(h5err%errBool) write(*,*) h5err%errorMsg

    if (i_proc == 0) then
       h5in%units="dimensionless"
       h5in%mesh="/t_current"
       call add_h5(dumpTGid,'k_perp_squared',k_perp_squared,h5in,data_step,h5err)
       if(h5err%errBool) write(*,*) h5err%errorMsg
    endif !i_proc == 0
  endif !io_method > 1
#endif

  if (i_proc == 0) then
     if(io_method < 3) then
     call write_local_real(&
          trim(path)//'out.gyro.k_perp_squared',&
          10,&
          size(k_perp_squared),&
          k_perp_squared)
     endif
  endif

  call gyro_field_fluxave

  !-------------------------------------------------------------------
  ! Calculation of fundamental nonlinear fluxes and related 
  ! diffusivities
  !
  call gyro_nonlinear_flux
  call gyro_diffusivity
  call gyro_gbflux
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  ! Output specific to linear/nonlinear operation:
  !
  if (nonlinear_flag == 0) then

     !=============
     ! BEGIN LINEAR 
     !=============
#ifndef HAVE_HDF5
     call gyro_write_freq(trim(path)//'out.gyro.freq',10)
#else
     h5in%mesh="/t_current"
     call gyro_write_freq(trim(path)//'out.gyro.freq',10,dumpTGid,h5in,h5err)
#endif
  
!srinathV...we'll have to revisit getting the linear run info 

     if (plot_u_flag == 1) then        

        ! PHI
#ifndef HAVE_HDF5
        call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_phi',10,1,0)
#else
        h5in%mesh="/t_current"
        call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_phi',10,1,0, &
            "ballon_phi",dumpTGid,h5in,h5err)
#endif

        if (n_field > 1) then
           ! A_PARALLEL 
#ifndef HAVE_HDF5
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_a',10,2,0)
#else
           h5in%mesh="/t_current"
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_a',10,2,0, &
              "balloon_a",dumpTGid,h5in,h5err)
#endif
        endif

        if (n_field > 2) then
           ! B_PARALLEL 
#ifndef HAVE_HDF5
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_aperp',10,3,0)
#else
           h5in%mesh="/t_current"
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_aperp', &
            10,3,0, "balloon_aperp",dumpTGid,h5in,h5err)
#endif

        endif

        ! E_PARALLEL
        if (eparallel_plot_flag == 1) then
#ifndef HAVE_HDF5
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_epar',10,n_field+1,0)
#else
           h5in%mesh="/t_current"
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_epar',10, &
            n_field+1,0,"balloon_epar",dumpTGid,h5in,h5err)
#endif
        endif

     endif

     if (plot_n_flag == 1) then

        ! DENSITY
        if (electron_method /= 3) then
#ifndef HAVE_HDF5
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_n_ion',10,5,1)
#else
           h5in%mesh="/t_current"
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_n_ion', &
            10,5,1,"n_ion",dumpTGid,h5in,h5err)
#endif

        endif
        if (electron_method > 1) then 
#ifndef HAVE_HDF5
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_n_elec',10,5,indx_e)
#else
           h5in%mesh="/t_current"
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_n_elec', &
            10,5,indx_e,"n_elec",dumpTGid,h5in,h5err)
#endif
        endif
     endif

     if (plot_e_flag == 1) then

        ! ENERGY
        if (electron_method /= 3) then
#ifndef HAVE_HDF5
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_e_ion',10,6,1)
#else
           h5in%mesh="/t_current"
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_e_ion', &
            10,6,1,"e_ion",dumpTGid,h5in,h5err)
#endif

        endif
        if (electron_method > 1) then 
#ifndef HAVE_HDF5
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_e_elec',10,6,indx_e)
#else
           h5in%mesh="/t_current"
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_e_elec', &
            10,6,indx_e,"e_elec",dumpTGid,h5in,h5err)
#endif
        endif
     endif

     if (plot_v_flag == 1) then

        ! ENERGY
        if (electron_method /= 3) then
#ifndef HAVE_HDF5
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_v_ion',10,7,1)
#else
           h5in%mesh="/t_current"
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_v_ion', &
              10,7,1,"v_ion",dumpTGid,h5in,h5err)
#endif
        endif
        if (electron_method > 1) then 
#ifndef HAVE_HDF5
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_v_elec',10,7,indx_e)
#else
           h5in%mesh="/t_current"
           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_v_elec', &
            10,7,indx_e,"v_elec",dumpTGid,h5in,h5err)
#endif
        endif
     endif

     !-----------------------------------------------------------------
     ! Distribution function data:
     !
     if (n_proc == 1 .and. n_n == 1 .and. dist_print == 1) then
        call gyro_write_h(trim(path)//'out.gyro.hp',trim(path)//'out.gyro.ht',10,11)
     endif
     !-----------------------------------------------------------------

     if (i_proc == 0 .and. lindiff_method > 1) then
        if(io_method <3) then
          call write_local_real( &
               trim(path)//'out.gyro.diff',10,size(diff),diff) 
          call write_local_real( &
               trim(path)//'out.gyro.diff_i',10,size(diff_i),diff_i)
          call write_local_real( &
               trim(path)//'out.gyro.gbflux',10,size(gbflux),gbflux)
          call write_local_real( &
               trim(path)//'out.gyro.gbflux_mom',10,size(gbflux_mom),gbflux_mom)
          call write_local_real( &
               trim(path)//'out.gyro.gbflux_i',10,size(gbflux_i),gbflux_i)

          if (trapdiff_flag == 1) then
             call write_local_real( &
                  trim(path)//'out.gyro.diff_trapped',&
                  10,size(diff_trapped),diff_trapped)
             call write_local_real( &
                  trim(path)//'out.gyro.diff_i_trapped',&
                  10,size(diff_i_trapped),diff_i_trapped)
             call write_local_real( &
                  trim(path)//'out.gyro.gbflux_trapped',&
                  10,size(gbflux_trapped),gbflux_trapped)
             call write_local_real( &
                  trim(path)//'out.gyro.gbflux_i_trapped',&
                  10,size(gbflux_i_trapped),gbflux_i_trapped)
          endif
        endif !io_method < 3
#ifdef HAVE_HDF5
        if(io_method > 1 ) then 
        h5in%mesh="/t_current"
        call add_h5(dumpTGid,'diff',diff,h5in,data_step,h5err)
        call add_h5(dumpTGid,'diff_i',diff_i,h5in,data_step,h5err)
        call add_h5(dumpTGid,'gbflux',gbflux,h5in,data_step,h5err)
        call add_h5(dumpTGid,'gbflux_mom',gbflux_mom,h5in,data_step,h5err)
        call add_h5(dumpTGid,'gbflux_i',gbflux_i,h5in,data_step,h5err)

        if (trapdiff_flag == 1) then
           h5in%mesh="/t_current"
           call add_h5(dumpTGid,'diff_trapped',diff_trapped,h5in,data_step,h5err)
           call add_h5(dumpTGid,'diff_i_trapped',diff_i_trapped,h5in,data_step,h5err)
           call add_h5(dumpTGid,'gbflux_trapped',gbflux_trapped,h5in,data_step,h5err)
           call add_h5(dumpTGid,'gbflux_i_trapped',gbflux_i_trapped,h5in,data_step,h5err)
        endif !trapdiff_flag == 1
      endif !io_method > 1
#endif
    endif !i_proc ==0 and lindiff >1 

     if (lindiff_method >= 4) then
       if(io_method < 3 ) then

        call write_distributed_real(&
             trim(path)//'out.gyro.diff_n',&
             10,&
             size(diff_n),&
             diff_n)

        call write_distributed_real(&
             trim(path)//'out.gyro.phi_squared_QL_n',&
             10,&
             size(phi_squared_QL_n),&
             phi_squared_QL_n)

        call write_distributed_real(&
             trim(path)//'out.gyro.g_squared_QL_n',&
             10,&
             size(g_squared_QL_n),&
             g_squared_QL_n)

        call write_distributed_real(&
             trim(path)//'out.gyro.gbflux_n',&
             10,&
             size(gbflux_n),&
             gbflux_n)
      endif !io_method < 3
#ifdef HAVE_HDF5
        if(io_method > 1 ) then
           h5in%mesh="/t_current"
           call write_distributed_real_h5('diff_n',dumpTGid,&
                 n_kinetic,n_field,2, &
                 size(diff_n),&
                 diff_n,&
                 h5in,h5err)

            call write_distributed_real_h5('phi_squared_QL_n',dumpTGid,&
                 1,1,1,&
                 size(phi_squared_QL_n),&
                 phi_squared_QL_n,&
                 h5in,h5err)
            call write_distributed_real_h5('g_squared_QL_n',dumpTGid,&
                 size(g_squared_QL_n),&
                 3,1,1,&
                 g_squared_QL_n,&
                 h5in,h5err)

            call write_distributed_real_h5("gbflux_n",dumpTGid,&
                 n_kinetic,n_field,4,&
                 size(gbflux_n),&
                 gbflux_n,&
                 h5in,h5err)
        endif !io_method > 1
#endif
     endif !lindiff_method >= 4

     !=============
     ! END LINEAR 
     !=============

  else !nonlinear !=0

     !================
     ! BEGIN NONLINEAR 
     !================
    if (io_method < 3) then
     call write_distributed_real(&
          trim(path)//'out.gyro.diff_n',&
          10,&
          size(diff_n),&
          diff_n)

     call write_distributed_real(&
          trim(path)//'out.gyro.gbflux_n',&
          10,&
          size(gbflux_n),&
          gbflux_n)

     if (lindiff_method >= 4) then
        call write_distributed_real(&
             trim(path)//'out.gyro.phi_squared_QL_n',&
             10,&
             size(phi_squared_QL_n),&
             phi_squared_QL_n)
        call write_distributed_real(&
             trim(path)//'out.gyro.g_squared_QL_n',&
             10,&
             size(g_squared_QL_n),&
             g_squared_QL_n)
     endif !lindiff_method >= 4

     if (nonlinear_transfer_flag == 1) then
        call write_distributed_real(&
             trim(path)//'out.gyro.nl_transfer',&
             10,&
             size(nl_transfer),&
             nl_transfer)
     endif !nonlinear_transfer_flag ==1 
    endif !io_method <3
#ifdef HAVE_HDF5
  if (io_method > 1 ) then
        h5in%units="diff units"
        h5in%mesh="/t_current"
     call write_distributed_real_h5("diff_n",dumpTGid,&
          n_kinetic,n_field,2,&
          size(diff_n),&
          diff_n,&
          h5in,h5err)
     if(h5err%errBool) write(*,*) h5err%errorMsg

     call write_distributed_real_h5("gbflux_n",dumpTGid,&
          n_kinetic,n_field,4,&
          size(gbflux_n),&
          gbflux_n,&
          h5in,h5err)
     if(h5err%errBool) write(*,*) h5err%errorMsg

     if (lindiff_method >= 4) then
        call write_distributed_real_h5('phi_squared_QL_n',dumpTGid,&
             1,1,1,&
             size(phi_squared_QL_n),&
             phi_squared_QL_n,&
             h5in,h5err)
        if(h5err%errBool) write(*,*) h5err%errorMsg
        call write_distributed_real_h5('g_squared_QL_n',dumpTGid,&
             size(g_squared_QL_n),&
             3,1,1,&
             g_squared_QL_n,&
             h5in,h5err)
        if(h5err%errBool) write(*,*) h5err%errorMsg
     endif !lindiff_method >=4

     if (nonlinear_transfer_flag == 1) then
        call write_distributed_real_h5('nl_transfer',dumpTGid,&
             n_x,2,1,&
             size(nl_transfer),&
             nl_transfer,&
             h5in,h5err)
        if(h5err%errBool) write(*,*) h5err%errorMsg
     endif !nonlinear_transfer_flag==1
    endif !io_method > 1
#endif
    

     if (i_proc == 0) then
        if(io_method < 3) then
        call write_local_real(trim(path)//'out.gyro.field_rms',10,size(ave_phi),ave_phi)

        call write_local_real( &
             trim(path)//'out.gyro.diff',10,size(diff),diff)
        call write_local_real( &
             trim(path)//'out.gyro.diff_i',10,size(diff_i),diff_i)

        call write_local_real( &
             trim(path)//'out.gyro.gbflux',10,size(gbflux),gbflux)
        call write_local_real( &
             trim(path)//'out.gyro.gbflux_mom',10,size(gbflux_mom),gbflux_mom)
        call write_local_real( &
             trim(path)//'out.gyro.gbflux_i',10,size(gbflux_i),gbflux_i)

        if (trapdiff_flag == 1) then
           call write_local_real( &
                trim(path)//'out.gyro.diff_trapped',&
                10,size(diff_trapped),diff_trapped)
           call write_local_real( &
                trim(path)//'out.gyro.diff_i_trapped',&
                10,size(diff_i_trapped),diff_i_trapped)
           call write_local_real( &
                trim(path)//'out.gyro.gbflux_trapped',10,&
                size(gbflux_trapped),gbflux_trapped)
           call write_local_real( &
                trim(path)//'out.gyro.gbflux_i_trapped',10,&
                size(gbflux_i_trapped),gbflux_i_trapped)
        endif !trapdiff_flag == 1

        call write_local_real( &
             trim(path)//'out.gyro.zerobar',10,&
             size(field_fluxave),transpose(field_fluxave))

       endif !io_method < 3
#ifdef HAVE_HDF5
      if(io_method > 1 ) then
        h5in%units="diff units"
        h5in%mesh="/t_current"
        call add_h5(dumpTGid,'field_rms',ave_phi,h5in,data_step,h5err)
        if(h5err%errBool) write(*,*) h5err%errorMsg
        call add_h5(dumpTGid,'diff',diff,h5in,data_step,h5err)
        if(h5err%errBool) write(*,*) h5err%errorMsg
        call add_h5(dumpTGid,'diff_i',diff_i,h5in,data_step,h5err)
        if(h5err%errBool) write(*,*) h5err%errorMsg
        call add_h5(dumpTGid,'gbflux',gbflux,h5in,data_step,h5err)
        if(h5err%errBool) write(*,*) h5err%errorMsg
        call add_h5(dumpTGid,'gbflux_mom',gbflux_mom,h5in,data_step,h5err)
        if(h5err%errBool) write(*,*) h5err%errorMsg
        call add_h5(dumpTGid,'gbflux_i',gbflux_i,h5in,data_step,h5err)
        if(h5err%errBool) write(*,*) h5err%errorMsg

        if (trapdiff_flag == 1) then
           call add_h5(dumpTGid,'diff_trapped',diff_trapped,h5in,data_step,h5err)
           if(h5err%errBool) write(*,*) h5err%errorMsg
           call add_h5(dumpTGid,'diff_i_trapped',diff_i_trapped,h5in,data_step,h5err)
           if(h5err%errBool) write(*,*) h5err%errorMsg
           call add_h5(dumpTGid,'gbflux_trapped',gbflux_trapped,h5in,data_step,h5err)
           if(h5err%errBool) write(*,*) h5err%errorMsg
           call add_h5(dumpTGid,'gbflux_i_trapped',gbflux_i_trapped,h5in,data_step,h5err)
           if(h5err%errBool) write(*,*) h5err%errorMsg
        endif

        call add_h5(dumpTGid,'zerobar',transpose(field_fluxave),h5in,data_step,h5err)
        if (h5err%errBool) write(*,*) h5err%errorMsg
      endif !io_method > 1
#endif

        allocate(a3(n_kinetic,4,n_x))
        do i=1,n_x
           a3(:,1,i) = h0_n(:,i)
           a3(:,2,i) = h0_e(:,i)
           a3(:,3,i) = source_n(:,i)
           a3(:,4,i) = source_e(:,i)
        enddo

        if (io_method < 3) then
          call write_local_real( &
               trim(path)//'out.gyro.source',10,size(a3),a3)

          call write_local_real( &
               trim(path)//'out.gyro.moment_zero',10,&
               size(moments_zero_plot),moments_zero_plot)
        endif !io_method <3)
#ifdef HAVE_HDF5
      if (io_method > 1 ) then
        h5in%mesh="/t_current"
        call add_h5(dumpTGid,'source',a3,h5in,data_step,h5err)
        if(h5err%errBool) write(*,*) h5err%errorMsg

        call add_h5(dumpTGid,'moments_zero',moments_zero_plot,h5in,data_step,h5err)
        if(h5err%errBool) write(*,*) h5err%errorMsg
      endif !io_method > 1
#endif
       deallocate(a3)
     endif! i_proc ==0

     !================
     ! END NONLINEAR 
     !================

  endif
  !-------------------------------------------------------------------

  call gyro_write_error(trim(path)//'out.gyro.error',10)

  !------------------------------------------------------------
  ! Entropy diagnostics
  !
  if (entropy_flag == 1) then
     call gyro_entropy 
     if (i_proc == 0) then 
        if(io_method < 3) then
          call write_local_real(&
            trim(path)//'out.gyro.entropy.out',10,size(entropy),entropy)
        endif
#ifdef HAVE_HDF5
        if(io_method > 1 ) then
          h5in%mesh="/t_current"
          call add_h5(dumpTGid,'entropy',entropy,h5in,data_step,h5err)
          if(h5err%errBool) write(*,*) h5err%errorMsg
        endif
#endif
     endif
  endif
  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Velocity-space diagnostics
  !
  if (velocity_output_flag == 1) then
     call gyro_nonlinear_flux_velocity
     call write_distributed_real(&
          trim(path)//'out.gyro.flux_velocity',&
          10,&
          size(nonlinear_flux_velocity),&
          nonlinear_flux_velocity)
  endif

!srinathv: need and hdf5 implementiotion for above
!flux_velocity size:n_energy,n_lambda,n_kinetic,n_field,n_moment

  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Write precision-monitoring data
  !
  call gyro_write_precision(10,sum(abs(gbflux)))
  !------------------------------------------------------------

  !------------------------------------------------------------
  call gyro_write_timers(trim(path)//'out.gyro.timing',10)
  !------------------------------------------------------------

  if (i_proc == 0 .and. debug_flag == 1) print *,'[gyro_write_timedata done]'


#ifdef HAVE_HDF5
  if(io_method >1 ) then
    if (i_proc == 0) then
       h5in%mesh="mesh-structured"
       h5in%units="time"
       skip_at_tzero : if (.not. hdf5_skip) then
         call add_h5(dumpTGid,'data_step',data_step,h5in,data_step,h5err)
          if(h5err%errBool) write(*,*) h5err%errorMsg
         call add_h5(dumpTGid,'t_current',t_current,h5in,data_step,h5err)
          if(h5err%errBool) write(*,*) h5err%errorMsg
        endif skip_at_tzero

       h5in%mesh='mesh-structured'
       ! dump in the field
        call dump_h5(dumpGid,'data_step',data_step,h5in,h5err)
        if(h5err%errBool) write(*,*) h5err%errorMsg
        call dump_h5(dumpGid,'t_current',t_current,h5in,h5err)
        if(h5err%errBool) write(*,*) h5err%errorMsg
        call dump_h5(dumpGid,'n_proc',n_proc,h5in,h5err)
        if(h5err%errBool) write(*,*) h5err%errorMsg

      ! dump 3D
      if (write_threed) then
        h5in%mesh='mesh-structured'
        call dump_h5(gid3D,'data_step',data_step,h5in,h5err)
        if(h5err%errBool) write(*,*) h5err%errorMsg
        call dump_h5(gid3D,'t_current',t_current,h5in,h5err)
        if(h5err%errBool) write(*,*) h5err%errorMsg
        call dump_h5(gid3D,'n_proc',n_proc,h5in,h5err)
        if(h5err%errBool) write(*,*) h5err%errorMsg
      endif
 
      call close_h5file(dumpTFid,dumpTGid,h5err)
         if(h5in%debug) write(*,*) "closing dumpTFid"
         if(h5err%errBool) then
               write(*,*) h5err%errorMsg, " for dumpTFid"
            return
         endif
         call close_h5file(dumpFid,dumpGid,h5err)
         if(h5in%debug) write(*,*) "closing dumpFid"
         if(h5err%errBool) then
               write(*,*) h5err%errorMsg, " for dumpFid"
            return
         endif
         if (write_threed) then
           call close_h5file(fid3d,gid3d,h5err)
             if(h5in%debug) write(*,*) "closing fid3d"
             if(h5err%errBool) then
               write(*,*) h5err%errorMsg, " for fid3d"
               return
            endif
         endif
      endif
  endif !io_method > 1
#endif

10 format(t2,a,t24,es9.3)

end subroutine gyro_write_timedata

#ifdef HAVE_HDF5
!===========================================================================
subroutine write_hdf5_restart

  use mpi
  use gyro_globals
  !use hdf5
  use hdf5_api

  !---------------------------------------------------
  implicit none
  !
  character(60) :: description
  character(64) :: step_name, tempVarName
  character(128) :: dumpfile
  integer(HID_T) :: dumpGid,dumpFid
  type(hdf5ErrorType) :: errval
  character(4) :: iname
  type(hdf5InOpts) :: h5in
  type(hdf5ErrorType) :: h5err
  integer :: number_label

  !---------------------------------------------------
  if (i_proc == 0) then
     call vshdf5_inith5vars(h5in, h5err)
     h5in%comm=MPI_COMM_SELF
     h5in%info=MPI_INFO_NULL
     h5in%wrd_type=H5T_NATIVE_DOUBLE
     !h5in%wrd_type=h5in%h5_kind_type_r8
     h5in%doTranspose=.true.
     h5in%vsTime=t_current
     h5in%wrVsTime=.true.
     h5in%verbose=.true.

     !---------------------------------------------------
     ! Timestep data:
     !
     number_label=NINT(t_current/dt)
     if (number_label>999999) THEN
        write(step_name,fmt='(i7.7)') number_label
     else if (data_step>99999) THEN
        write(step_name,fmt='(i6.6)') number_label
     else
        write(step_name,fmt='(i5.5)') number_label
     endif

     
     dumpfile=TRIM(path)//"gyroRestart"//TRIM(step_name)//".h5"
     description="GYRO restart file"
     call open_newh5file(dumpfile,dumpFid,description,dumpGid,h5in,h5err)

     h5in%mesh=' '
     call write_attribute(dumpGid,"data_step",data_step,errval)
     call write_attribute(dumpGid,"n_proc",n_proc,errval)
     call write_attribute(dumpGid,"i_restart",i_restart,errval)
  endif

  if (n_proc-1 == 0) then
     call dump_h5(dumpGid,'h_0_real',REAL(h_0),h5in,h5err)
     call dump_h5(dumpGid,'h_0_imag',AIMAG(h_0),h5in,h5err)
     call close_h5file(dumpFid,dumpGid,h5err)
     return
  endif

  do i_proc_w=1,n_proc-1

     if (i_proc == 0) then
        call MPI_RECV(h_0,&
             size(h_0),&
             MPI_DOUBLE_COMPLEX,&
             i_proc_w,&
             i_proc_w,&
             GYRO_COMM_WORLD,&
             recv_status,&
             i_err)

        WRITE(iname,fmt='(i4.4)') i_proc_w
        tempVarName="h_0_real"//iname
        call dump_h5(dumpGid,tempVarName,REAL(h_0),h5in,h5err)
        tempVarName="h_0_imag"//iname
        call dump_h5(dumpGid,tempVarName,AIMAG(h_0),h5in,h5err)

     else if (i_proc == i_proc_w) then
        call MPI_SEND(h,&
             size(h),&
             MPI_DOUBLE_COMPLEX,&
             0,&
             i_proc_w,&
             GYRO_COMM_WORLD,&
             i_err)
     endif
  enddo

  if (i_proc == 0) call close_h5file(dumpFid,dumpGid,h5err)

end subroutine write_hdf5_restart


subroutine myhdf5_close

  use gyro_globals
  !use hdf5
  use hdf5_api

  implicit none
  integer :: ierr

  call h5close_f(ierr)

end subroutine myhdf5_close
#endif

!===========================================================================
!------------------------------------------------------
! write_distributed_real.f90
!
! PURPOSE:
!  Control merged output of distributed real array.
!------------------------------------------------------

subroutine write_distributed_real(datafile,io,n_fn,fn)

  use mpi
  use gyro_globals, only : &
       n_n,&
       n_n_1,&
       n_proc_1,&
       recv_status,&
       data_step,&
       GYRO_COMM_WORLD,&
       i_proc,&
       i_err,&
       io_control, &
       fmtstr

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io
  integer, intent(in) :: n_fn
  real, intent(in) :: fn(n_fn)
  !
  integer :: data_loop
  integer :: i_group_send
  integer :: i_send
  integer :: in
  !
  real :: fn_recv(n_fn)
  !------------------------------------------------------

  select case (io_control)

  case(0)

     return 

  case(1)
 
     ! Initial open

     if (i_proc == 0) then
        open(unit=io,file=datafile,status='replace')
        close(io)
     endif

  case(2)

     ! Output

     if (i_proc == 0) &
          open(unit=io,file=datafile,status='old',position='append')

     do in=1,n_n

        !-----------------------------------------
        ! Subgroup collector:
        !
        i_group_send = (in-1)/n_n_1

        if (i_group_send /= 0) then

           i_send = i_group_send*n_proc_1

           if (i_proc == 0) then

              call MPI_RECV(fn_recv,&
                   n_fn,&
                   MPI_DOUBLE_PRECISION,&
                   i_send,&
                   in,&
                   GYRO_COMM_WORLD,&
                   recv_status,&
                   i_err)

           else if (i_proc == i_send) then

              call MPI_SEND(fn,&
                   n_fn,&
                   MPI_DOUBLE_PRECISION,&
                   0,&
                   in,&
                   GYRO_COMM_WORLD,&
                   i_err)

           endif

        else

           fn_recv(:) = fn(:)

        endif
        !
        !-----------------------------------------

        if (i_proc == 0) then

           write(io,fmtstr) fn_recv(:)

        endif

     enddo ! in

     if (i_proc == 0) close(io)

  case(3)

     ! Reposition after restart

     if (i_proc == 0) then

        open(unit=io,file=datafile,status='old')
        do data_loop=0,data_step

           do in=1,n_n
              read(io,fmtstr) fn_recv(:)
           enddo

        enddo ! data_loop

        endfile(io)
        close(io)

     endif

  end select

end subroutine write_distributed_real

!===========================================================================
#ifdef HAVE_HDF5
subroutine write_distributed_real_h5(varName,rGid,n1,n2,n3,n_fn,fn,h5in,h5err)

  use mpi
  use hdf5_api
  use gyro_globals, only : &
       n_n,&
       n_n_1,&
       n_proc_1,&
       recv_status,&
       data_step,&
       GYRO_COMM_WORLD,&
       i_proc,&
       i_err,&
       electron_method!, &
      ! io_control

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: varName
  integer(HID_T), intent(in) :: rGid
  integer, intent(in) :: n_fn,n1,n2,n3
  real, intent(in) :: fn(n_fn)
  type(hdf5InOpts), intent(inout) :: h5in
  type(hdf5ErrorType), intent(inout) :: h5err
  !
  integer :: i_group_send
  integer :: i_send
  integer :: ifld,ikin,imom,i
  integer :: in

  real, dimension(:,:,:,:), allocatable :: buffn
  character(128) :: tempVarName 
  character(128), dimension(:,:,:),allocatable :: vnameArray
  character(1) :: ikin_name
  !
  real :: fn_recv(n_fn)
  !------------------------------------------------------


  !------------------------------------------------------
  !  set up the names for setting species and "field" name
  ALLOCATE(vnameArray(n3,3,4))
  vnameArray=" "
  do ikin=1,n3
     if(electron_method==2 .and. ikin==n3 ) THEN
        tempVarName=trim(varName)//"_electron"
     elseif(electron_method==3 .or. (electron_method==4.and.ikin==n3)) THEN
        tempVarName=trim(varName)//"_electron"
     else
        write(ikin_name,fmt='(i1.1)') ikin-1
        tempVarName=trim(varName)//"_ion"//ikin_name
     endif
     do i=1,4
        vnameArray(ikin,1,i)=trim(tempVarName)//"_phi"
        vnameArray(ikin,2,i)=trim(tempVarName)//"_Apar"
        vnameArray(ikin,3,i)=trim(tempVarName)//"_Bpar"
     enddo
     do i=1,3
        vnameArray(ikin,i,1)=trim(vnameArray(ikin,i,1))//"_density"
        vnameArray(ikin,i,2)=trim(vnameArray(ikin,i,2))//"_energy"
        vnameArray(ikin,i,3)=trim(vnameArray(ikin,i,3))//"_momentum"
        vnameArray(ikin,i,4)=trim(vnameArray(ikin,i,4))//"_energyExchange"
     enddo
  enddo

  !Rectangular array (n_kinetic,n_field,2=i,n_x,n_time)
  !  allocate(kxkyspec(n_x))
  !  allocate(gbflux_n(n_kinetic,n_field,p_moment))
  !  allocate(diff_n(n_kinetic,n_field,n_moment))
  !  real, dimension(1) :: phi_squared_QL_n
  !  real, dimension(3) :: g_squared_QL_n

  ! n1 = n_kinetic; n2 = n_field, n3=n_moment, n4=n_n
  if (i_proc==0) then
     allocate(buffn(n1,n2,n3,n_n)); buffn=0.
  endif
  !-----------------------------------------


  !------------------------------------------------------
  do in=1,n_n

     !-----------------------------------------
     ! Subgroup collector:
     !
     i_group_send = (in-1)/n_n_1

     if (i_group_send /= 0) then

        i_send = i_group_send*n_proc_1

        if (i_proc == 0) then

           call MPI_RECV(fn_recv,&
                n_fn,&
                MPI_DOUBLE_PRECISION,&
                i_send,&
                in,&
                GYRO_COMM_WORLD,&
                recv_status,&
                i_err)

        else if (i_proc == i_send) then

           call MPI_SEND(fn,&
                n_fn,&
                MPI_DOUBLE_PRECISION,&
                0,&
                in,&
                GYRO_COMM_WORLD,&
                i_err)

        endif

     else

        fn_recv(:) = fn(:)

     endif
     !
     !-----------------------------------------

     if (i_proc == 0) then
        buffn(:,:,:,in)=reshape(fn_recv,(/n1,n2,n3/))
     endif
  enddo ! in

  !-----------------------------------------
  if (i_proc /= 0) return


  if (n3==1) then
     if(n2==1) then
        call add_h5(rGid,trim(varName),buffn(:,1,1,:),h5in,data_step,h5err)
     else  
        call add_h5(rGid,trim(varName),buffn(:,:,1,:),h5in,data_step,h5err)
     endif
  else
     ! n1 = n_kinetic; n2 = n_field, n3=n_moment, n4=n_n
     do ikin=1,n1
        do ifld=1,n2
           do imom=1,n3
              tempVarName=trim(vnameArray(ikin,ifld,imom))
              call add_h5(rGid,trim(tempVarName),buffn(ikin,ifld,imom,:),h5in,data_step,h5err)
           enddo
        enddo
     enddo
  endif

  deallocate(buffn)
  deallocate(vnameArray)

end subroutine write_distributed_real_h5
#endif
!===========================================================================

!------------------------------------------------------
! write_distributed_complex.f90
!
! PURPOSE:
!  Control merged output of complex distributed array.
!------------------------------------------------------

subroutine write_distributed_complex(datafile,io,n_fn,fn)

  use mpi
  use gyro_globals, only : &
       n_n,&
       n_n_1,&
       n_proc_1,&
       recv_status,&
       data_step,&
       GYRO_COMM_WORLD,&
       i_proc,&
       i_err, &
       io_control, &
       fmtstr

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io
  integer, intent(in) :: n_fn
  complex, intent(in) :: fn(n_fn)
  !
  integer :: data_loop
  integer :: i_group_send
  integer :: i_send
  integer :: in
  !
  complex :: fn_recv(n_fn)
  !------------------------------------------------------


  select case (io_control)

  case(0)

     return

  case(1)

     ! Open

     if (i_proc == 0) then
        open(unit=io,file=datafile,status='replace')
        close(io)
     endif

  case(2)

     ! Append

     if (i_proc == 0) &
          open(unit=io,file=datafile,status='old',position='append')

     do in=1,n_n

        !-----------------------------------------
        ! Subgroup collector:
        !
        i_group_send = (in-1)/n_n_1

        if (i_group_send /= 0) then

           i_send = i_group_send*n_proc_1

           if (i_proc == 0) then

              call MPI_RECV(fn_recv,&
                   n_fn,&
                   MPI_DOUBLE_COMPLEX,&
                   i_send,&
                   in,&
                   GYRO_COMM_WORLD,&
                   recv_status,&
                   i_err)

           else if (i_proc == i_send) then

              call MPI_SEND(fn,&
                   n_fn,&
                   MPI_DOUBLE_COMPLEX,&
                   0,&
                   in,&
                   GYRO_COMM_WORLD,&
                   i_err)

           endif

        else

           fn_recv(:) = fn(:)

        endif
        !
        !-----------------------------------------

        if (i_proc == 0) then

           write(io,fmtstr) fn_recv(:)

        endif

     enddo ! in

     if (i_proc == 0) close(io)

  case(3)

     ! Rewind

     if (i_proc == 0) then

        open(unit=io,file=datafile,status='old')
        do data_loop=0,data_step

           do in=1,n_n
              read(io,fmtstr) fn_recv(:)
           enddo

        enddo ! data_loop

        endfile(io)
        close(io)

     endif

  end select

end subroutine write_distributed_complex
!===========================================================================
#ifdef HAVE_HDF5
subroutine write_distributed_complex_sorf_h5(vname,rGid,r3Did,&
     n_fn,n1,n2,n3,fn,plot3d,plotwedge,h5in,h5err)
! _sorf_ stands for species and fields. 
! this allows for the unrolling of the arrays 
! and indicating the species OR  field quantities.
  use mpi
  use hdf5_api
  use gyro_globals, only : &
       q, &
       n0,&
       d_n,&
       n_n,&
       n_n_1,&
       n_proc_1,&
       n_theta_plot,&
       t_current,&
       debug_flag,&
       recv_status,&
       data_step,&
       w0_s,&
       ir_norm,&
       GYRO_COMM_WORLD,&
       i_proc,&
       i_err, &
       electron_method,&
       n_torangle_wedge,&
       n_torangle_3d, &
       io_control, &
       alpha_phi, &
       alpha_phi_wedge

  !------------------------------------------------------
  !   mom: n1,n2,n3=n_theta_plot,n_x,n_kinetic
  !   phi: n1,n2,n3=n_theta_plot,n_x,n_field
  !------------------------------------------------------

  implicit none
  !
  real :: pi=3.141592653589793
  character(*), intent(in) :: vname
  integer(HID_T), intent(in) :: rGid,r3Did
  integer, intent(in) :: n_fn,n1,n2,n3
  complex, intent(in) :: fn(n_fn)
  logical, intent(in) :: plot3d, plotwedge
  character(128) :: tempVarName , tempVarNameGr
  character(128), dimension(:),allocatable :: vnameArray
  character(1) :: ikin_name
  integer(HID_T) :: grGid
  type(hdf5InOpts), intent(inout) :: h5in
  type(hdf5ErrorType), intent(inout) :: h5err
  !
  integer :: i_group_send, ispcs
  integer :: i_send, iphi, istart,nn,ikin,in, ix,nphi
  !
  complex :: fn_recv(n_fn), c_i
  complex, dimension(:,:,:,:), allocatable :: buffn
  real, dimension(:,:,:,:), allocatable:: real_buff
  real, dimension(:,:), allocatable:: alpha_loc
  real :: omega_exp
  logical :: iswedge

  !------------------------------------------------------

  if (io_control < 2) return

  c_i = (0.0,1.0)

  omega_exp = w0_s(ir_norm)

  if (.not. plotwedge) then
     iswedge=.false.
     allocate(buffn(0:n1,n2,n3,n_n)); buffn=0.
  else
     iswedge=.true.
     allocate(buffn(0:n1-1,n2,n3,n_n)); buffn=0.
  endif

  !when n3=n_kinetic
  !electron_method =1 => n3=n_ion (gk ions and addiabtic electrons )
  !electron_method =2 => n3=n_spec (gk ions and drift electrons)
  !electron_method =3 => n3=1 (gk electrons and addiabtic ions)
  !electron_method =4 => n3=n_ion (gk electrons and gk ions)

  if (trim(vname) /= "phi") then 
     ALLOCATE(vnameArray(n3))
     vnameArray=" "
     do ikin=1,n3
        if(electron_method==2 .and. ikin==n3 ) THEN
           tempVarName=trim(vname)//"_electron"
        elseif(electron_method==3 .or. (electron_method==4.and.ikin==n3)) THEN
           tempVarName=trim(vname)//"_electron"
        else
           write(ikin_name,fmt='(i1.1)') ikin-1
           tempVarName=trim(vname)//"_ion"//ikin_name
        endif
        vnameArray(ikin)=tempVarName
     enddo
  else
     ALLOCATE(vnameArray(3))
     vnameArray=" "
     vnameArray(1)="phi"
     vnameArray(2)="A_par"
     vnameArray(3)="B_par"
  endif

  do in=1,n_n
     !-----------------------------------------
     ! Subgroup collector:
     !
     i_group_send = (in-1)/n_n_1
     if (i_group_send /= 0) then
        i_send = i_group_send*n_proc_1
        if (i_proc == 0) then
           call MPI_RECV(fn_recv,&
                n_fn, MPI_DOUBLE_COMPLEX, i_send, in,&
                GYRO_COMM_WORLD, recv_status, i_err)
        else if (i_proc == i_send) then
           call MPI_SEND(fn,&
                n_fn, MPI_DOUBLE_COMPLEX, 0, in,&
                GYRO_COMM_WORLD, i_err)
        endif
     else
        fn_recv(:) = fn(:)
     endif

     if (i_proc == 0) then
        buffn(0:n1-1,:,:,in)=reshape(fn_recv,(/n1,n2,n3/))
     endif


  enddo ! in

  !-----------------------------------------
  if (i_proc /= 0) return
  !-----------------------------------------
  !-----------------------------------------
  ! Apply boundary conditions
  !-----------------------------------------
  if (.not.iswedge) then
     do in=1,n_n
        nn=n0+(in-1)*d_n
        !apply theta BC: z_n(r,,2*pi) = z_n(r,0)exp(I*n*(nu(r,2*pi)-nu(r,0)))
        !with nu(r,2*pi) - nu(r,0) = -2*pi*q by definition
        ! phase[*] = EXP(-2*!PI*C_I*n[i_n]*profile_data.q[*])
        do ix=1,n2
           buffn(n1,ix,:,in)=buffn(0,ix,:,in)*exp(-2.*pi*c_i*nn*q(ix))
        enddo
     enddo ! in
  endif

  !-----------------------------------------
  ! Dump each species independently
  !-----------------------------------------
  if (.not.iswedge) then
    h5in%mesh="/poloidalMesh/cartMesh"
  else
    h5in%mesh="/wedgeMesh/cartMesh"
  endif

  do ispcs=1,n3
     tempVarNameGr=trim(vnameArray(ispcs))//"_modes"
     call make_group(rGid,trim(tempVarNameGr),grGid,h5in,h5err)
     !call make_group(rGid,trim(tempVarNameGr),grGid,"",h5err)
     tempVarName=trim(vnameArray(ispcs))//"_real"
     call dump_h5(grGid,trim(tempVarName),real(buffn(:,:,ispcs,:)),h5in,h5err)
     tempVarName=trim(vnameArray(ispcs))//"_imag"
     call dump_h5(grGid,trim(tempVarName),aimag(buffn(:,:,ispcs,:)),h5in,h5err)
     CALL close_group(trim(tempVarNameGr),grGid,h5err)
  enddo ! in
  if(.not.plot3d) then
     deallocate(buffn)
     return
  endif

  !-----------------------------------------
  ! Tranform into real space
  !   mom: n1,n2,n3=n_theta_plot,n_x,n_kinetic
  !   phi: n1,n2,n3=n_theta_plot,n_x,n_field
  !-----------------------------------------
  nphi = 1 
  if (.not. iswedge) then
     if (n_torangle_3d > 0) nphi=n_torangle_3d
     allocate(real_buff(0:n1,n2,n3,nphi))
     allocate(alpha_loc(0:n1,n2))
  else
     if (n_torangle_wedge > 0) nphi=n_torangle_wedge
     allocate(real_buff(0:n1-1,n2,n3,nphi))
     allocate(alpha_loc(0:n1-1,n2))
  endif
  if (n0==0) then
     istart=2
     do iphi=1,nphi
        real_buff(:,:,:,iphi)=real(buffn(:,:,:,1))
     enddo
  else
     istart=1
     real_buff(:,:,:,:)=0.
  endif
  do iphi=1,nphi
     !Get alpha coordinate on either the coarse or wedge mesh.
     ! Include doppler shift here

     if (.not. iswedge) then
        alpha_loc=alpha_phi(:,:,iphi)+omega_exp*t_current
     else
        alpha_loc=alpha_phi_wedge(:,:,iphi)+omega_exp*t_current
     endif
     do in=istart,n_n
        nn=n0+(in-1)*d_n
        do ikin=1,n3
           real_buff(:,:,ikin,iphi)=real_buff(:,:,ikin,iphi)&
                +2.*real(buffn(:,:,ikin,in)*exp(-c_i*nn*alpha_loc(:,:)))
        enddo
     enddo
  enddo

  deallocate(buffn,alpha_loc)

  ! Mapping of the variable names to array indices depends on input types


  if (.not. iswedge) then
     h5in%mesh="/threeDMesh/cartMesh"
     do ikin=1,n3
        tempVarName=trim(vnameArray(ikin))
        call dump_h5(r3Did,trim(tempVarName),real_buff(:,:,ikin,:),h5in,h5err)
     enddo
  else
     h5in%mesh="/wedgeMesh/cartMesh"
     ! Dump each phi slice as a separate variable
     do ikin=1,n3
        tempVarNameGr=trim(vnameArray(ikin))//"_toroidal"
        call make_group(r3Did,trim(tempVarNameGr),grGid,h5in,h5err)
        call dump_h5(grGid,trim(vnameArray(ikin)),real_buff(:,:,ikin,:),h5in,h5err)
        call close_group(trim(tempVarNameGr),grGid,h5err)
     enddo
  endif

  deallocate(real_buff)
  deallocate(vnameArray)

end subroutine write_distributed_complex_sorf_h5
#endif
!===========================================================================

!------------------------------------------------------
! write_local_real.f90
!
! PURPOSE:
!  This routine write a vector of nondistributed reals.
!------------------------------------------------------

subroutine write_local_real(datafile,io,n_fn,fn)

  use gyro_globals, only : &
       data_step, &
       io_control, &
       fmtstr

  !---------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io
  integer, intent(in) :: n_fn
  real, intent(in) :: fn(n_fn)
  !
  integer :: data_loop
  real :: dummy(n_fn)
  !---------------------------------------------------

  select case (io_control)

  case(0)

     return

  case(1)

     ! Open

     open(unit=io,file=datafile,status='replace')
     close(io)

  case(2)

     ! Append

     open(unit=io,file=datafile,status='old',position='append')
     write(io,fmtstr)  fn(:)
     close(io)

  case(3)

     ! Rewind

     open(unit=io,file=datafile,status='old')

     do data_loop=0,data_step
        read(io,fmtstr) dummy(:)
     enddo

     endfile(io)
     close(io)

  end select

end subroutine write_local_real

