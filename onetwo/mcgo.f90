!
      MODULE MCGO
      USE param,   only : kj,kb,ke,le_mcgo,kz,kjeb

      IMPLICIT NONE
!
      real*8    rho_mcgo(kj,kb),prespar_mcgo(kj,ke,kb), &
                presprp_mcgo(kj,ke,kb), &
                density_mcgo(kj,le_mcgo,kb),curdens_mcgo(kj,le_mcgo,kb), &
                dtneut_mcgo(kj,le_mcgo,kb),ddneut_mcgo(kj,le_mcgo,kb), &
                enere_mcgo(kj,le_mcgo,kb),eneri_mcgo(kj,le_mcgo,kb), &
                psiz_mcgo(kz,kb),fbth_mcgo(kj,ke,kb), &
                eng_cx_mcgo(kj,le_mcgo,kb),eng_cx_12_mcgo(kj,ke,kb), &
                eng_forb_mcgo(kj,le_mcgo,kb),eng_forb_12_mcgo(kj,ke,kb), &
                enere_12_mcgo(kj,ke,kb),eneri_12_mcgo(kj,ke,kb), &
                density_12_mcgo(kj,ke,kb),prespar_12_mcgo(kj,ke,kb), &
                presprp_12_mcgo(kj,ke,kb),dtneut_12_mcgo(kj,ke,kb), &
                ddneut_12_mcgo(kj,ke,kb),curdens_12_mcgo(kj,ke,kb), &
                fbth_12_mcgo(kj,ke,kb),beam_net_cur_mcgo(kz,kb), &
                beam_net_cur_12_mcgo(kj,kb),forbit_mcgo(ke,kb), &
                fpcx_mcgo(ke,kb),fene_mcgo(kz,ke,kb), &
                feni_mcgo(kz,ke,kb),fbe_12_mcgo(kj,ke,kb), &
                fbi_12_mcgo(kj,ke,kb)
!
      integer   npart_mcgo,npart_mcgo_max, &
                npart_full(kb),nparts_12(ke), &
                npart_half(kb),npart_third(kb),onetwo_beam, &
                spawn_mcgo,read_mcgo_file(kb),mcgo_info_available, &
                mcgo_fast_ion_target,use_Callen
      integer,private :: j,ic,ib
!
      character* 24 mcgo_input_file(kb), mcgo_output_file(kb), &
                    mcgo_input_file2(kb), mcgo_output
      character* 16 mcgo_env_path_name ! stores name of environment..
!                                      ..variable
      character*256 mcgo_path          ! stores path that..
!                                      ..mcgo_env_path_name points to



     data   mcgo_input_file(1)      ,      &
            mcgo_input_file(2)      ,      &
            mcgo_output_file(1)     ,      &
            mcgo_output_file(2)     ,      &
            mcgo_input_file2(1)     ,      &
            mcgo_input_file2(2)     ,      &
            mcgo_env_path_name      ,      &
            npart_mcgo_max          ,      &
            read_mcgo_file(1)       ,      &
            read_mcgo_file(2)       ,      &
            spawn_mcgo              ,      &
            mcgo_info_available            &
           /'mcgo_input_12_beam1'   ,      &
             'mcgo_input_12_beam2'  ,      &
            ' '                     ,      &
            ' '                     ,      &
            'mcgo_input_freya_beam1',      &
            'mcgo_input_freya_beam2',      &
            'MCGO_PATH'             ,      &
            10000                   ,      &
            0                       ,      &
            0                       ,      &
            0                       ,      &
            0                        /
!

      data presprp_12_mcgo(:,:,:) /kjeb*0.0 /
      data prespar_12_mcgo(:,:,:) /kjeb*0.0 /


!
!      common /mcgo/ &
!                    rho_mcgo,psiz_mcgo,prespar_mcgo,presprp_mcgo, &
!                    curdens_mcgo,dtneut_mcgo,ddneut_mcgo, &
!                    enere_mcgo,eneri_mcgo,density_mcgo, &
!                    enere_12_mcgo,eneri_12_mcgo, &
!                    density_12_mcgo,prespar_12_mcgo, &
!                    presprp_12_mcgo,dtneut_12_mcgo, &
!                    ddneut_12_mcgo,curdens_12_mcgo, &
!                    beam_net_cur_12_mcgo,eng_cx_mcgo, &
!                    forbit_mcgo,fpcx_mcgo,fene_mcgo,feni_mcgo, &
!                    fbth_12_mcgo,fbth_mcgo,beam_net_cur_mcgo, &
!                    fbe_12_mcgo,fbi_12_mcgo,eng_cx_12_mcgo, &
!                    eng_forb_12_mcgo,eng_forb_mcgo, &
!                    npart_mcgo_max,npart_mcgo,npart_full, &
!                    npart_half,npart_third,onetwo_beam,spawn_mcgo, &
!                    read_mcgo_file,mcgo_info_available,nparts_12, &
!                    mcgo_fast_ion_target,use_Callen
!
!      common /mcgo_files/ &
!                    mcgo_input_file , mcgo_output_file, &
!                    mcgo_input_file2, mcgo_output
!
!      common /mcgo_path_info/ &
!                    mcgo_path, &
!                    mcgo_env_path_name
!
! FILE USAGE:
!
!   mcgo_input_file      = 'mcgo_input_12_beam1'
!                          'mcgo_input_12_beam2'
!                           an input file for MCGO created by ONETWO
!
!   mcgo_input_file2(ib) = 'mcgo_input_freya_beam1'
!                          'mcgo_input_freya_beam2'
!                           an input file for MCGO created by ONETWO
!                           The fast ion deposition determined by FREYA.
!                           This file is unformatted for speed.
!                           There is one such file for each beam.
!
!   mcgo_output_file(ib) = 'mcgo_12_output_beam1'
!                          'mcgo_12_output_beam2'
!                           The output from MCGO using the above two
!                           files as input. This output file is read
!                           by ONETWO and is also unformatted, one
!                           file for each beam.
!
    END MODULE MCGO
