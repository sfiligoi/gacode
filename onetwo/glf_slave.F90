SUBROUTINE glf_slave(nj)

   USE mpi12
   IMPLICIT NONE
   
   INTEGER nj,i
   INTEGER jin,jout
   
   INTEGER IPARAM(32)
   REAL*8  RPARAM(32)
   
   INTEGER jeigen_iglf
   INTEGER jroot_iglf
   INTEGER glf23_iglf
   INTEGER jshoot
   INTEGER jmm
   INTEGER jmaxm
   INTEGER itport_pt(5)
   INTEGER irotstab
   INTEGER igrad
   INTEGER idengrad_glf
   INTEGER ibtflag_glf
   INTEGER i_delay
   
   REAL*8 zpte_glf(1)
   REAL*8 zpti_glf(1)
   REAL*8 zpne_glf(1)
   REAL*8 zpni_glf(1)
   REAL*8 btor_exp
   REAL*8 arho_iglf
   REAL*8 rmajor_glf
   REAL*8 zimp_glf
   REAL*8 amassimp_glf
   REAL*8 amassgas_glf
   REAL*8 exbmult_glf
   REAL*8 x_alpha_glf
   
   REAL*8 te_glf(nj)
   REAL*8 ti_glf(nj)
   REAL*8 ne_glf(nj)
   REAL*8 ni_glf(nj)
   REAL*8 ns_glf(nj)
   REAL*8 angrotp_exp(nj)
   REAL*8 egamma_exp(nj)
   REAL*8 gamma_p_exp(nj)
   REAL*8 vphim_glf(nj)
   REAL*8 vparm_glf(nj)
   REAL*8 vperm_glf(nj)
   REAL*8 zeff_exp(nj)
   REAL*8 rho_glf(nj)
   REAL*8 grho1_glf(nj)
   REAL*8 grho2_glf(nj)
   REAL*8 rmin_glf(nj)
   REAL*8 rmaj_glf(nj)
   REAL*8 q_glf(nj)
   REAL*8 shat_exp(nj)
   REAL*8 alpha_exp(nj)
   REAL*8 elong_glf(nj)
   
   REAL*8 diffnem
   REAL*8 chietem
   REAL*8 chiitim
   REAL*8 etaphim
   REAL*8 etaparm
   REAL*8 etaperm
   REAL*8 exchm
   
   REAL*8 diff_m(nj)
   REAL*8 chie_m(nj)
   REAL*8 chii_m(nj)
   REAL*8 etaphi_m(nj)
   REAL*8 etapar_m(nj)
   REAL*8 etaper_m(nj)
   REAL*8 exch_m(nj)
   REAL*8 egamma_m(nj)
   REAL*8 egamma_d(nj)
   REAL*8 gamma_p_m(nj)
   REAL*8 anrate_m(nj)
   REAL*8 anrate2_m(nj)
   REAL*8 anfreq_m(nj)
   REAL*8 anfreq2_m(nj)
   
   INTEGER ierr
   REAL*8 start_time, end_time
   REAL*8 glf_load_buf
   
   !print *,'SLAVE STARTING'
   
   start_time = mpi_wtime()
   
   CALL mpi12_recv_int(iparam,32,1)
   
   jeigen_iglf   = iparam( 1) 
   jroot_iglf    = iparam( 2) 
   glf23_iglf    = iparam( 3) 
   jshoot        = iparam( 4) 
   jmm           = iparam( 5) 
   jmaxm         = iparam( 6) 
   itport_pt(1)  = iparam( 7) 
   itport_pt(2)  = iparam( 8) 
   itport_pt(3)  = iparam( 9) 
   itport_pt(4)  = iparam(10) 
   itport_pt(5)  = iparam(11) 
   irotstab      = iparam(12) 
   igrad         = iparam(13) 
   idengrad_glf  = iparam(14) 
   ibtflag_glf   = iparam(15) 
   i_delay       = iparam(16) 
   
   !print *,'INT',iparam(1),iparam(2)
   
   CALL mpi12_recv_real(rparam,32,2)
   
   zpte_glf(1)   = rparam( 1)
   zpti_glf(1)   = rparam( 2)
   zpne_glf(1)   = rparam( 3)
   zpni_glf(1)   = rparam( 4)
   btor_exp      = rparam( 5)
!  rho_glf       = rparam( 6)
   arho_iglf     = rparam( 7)
   rmajor_glf    = rparam( 8)
   zimp_glf      = rparam( 9)
   amassimp_glf  = rparam(10)
   amassgas_glf  = rparam(11)
   exbmult_glf   = rparam(12)
   x_alpha_glf   = rparam(13)
   
   end_time = mpi_wtime()
   mpi_etimes(1) = mpi_etimes(1) + end_time-start_time
   
   !print *,'REAL',rparam(1),rparam(2)

   start_time = mpi_wtime()
   
   CALL mpi12_d_recv(te_glf     ,3)
   CALL mpi12_d_recv(ti_glf     ,3)
   CALL mpi12_d_recv(ne_glf     ,3)
   CALL mpi12_d_recv(ni_glf     ,3)
   CALL mpi12_d_recv(ns_glf     ,3)
   CALL mpi12_d_recv(angrotp_exp,3)
   CALL mpi12_d_recv(egamma_exp ,3)
   CALL mpi12_d_recv(gamma_p_exp,3)
   CALL mpi12_d_recv(vphim_glf  ,3)
   CALL mpi12_d_recv(vparm_glf  ,3)
   CALL mpi12_d_recv(vperm_glf  ,3)
   CALL mpi12_d_recv(zeff_exp   ,3)
   CALL mpi12_d_recv(rho_glf    ,3)  
   CALL mpi12_d_recv(grho1_glf  ,3)
   CALL mpi12_d_recv(grho2_glf  ,3)
   CALL mpi12_d_recv(rmin_glf   ,3)
   CALL mpi12_d_recv(rmaj_glf   ,3)
   CALL mpi12_d_recv(q_glf      ,3)
   CALL mpi12_d_recv(shat_exp   ,3)
   CALL mpi12_d_recv(alpha_exp  ,3)
   CALL mpi12_d_recv(elong_glf  ,3)
   
   !CALL MPI_BCAST(te_glf     ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   !CALL MPI_BCAST(ti_glf     ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   !CALL MPI_BCAST(ne_glf     ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   !CALL MPI_BCAST(ni_glf     ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   !CALL MPI_BCAST(ns_glf     ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   !CALL MPI_BCAST(angrotp_exp,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   !CALL MPI_BCAST(egamma_exp ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   !CALL MPI_BCAST(gamma_p_exp,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   !CALL MPI_BCAST(vphim_glf  ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   !CALL MPI_BCAST(vparm_glf  ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   !CALL MPI_BCAST(vperm_glf  ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   !CALL MPI_BCAST(zeff_exp   ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   !CALL MPI_BCAST(rho_glf    ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   !CALL MPI_BCAST(grho1_glf  ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   !CALL MPI_BCAST(grho2_glf  ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   !CALL MPI_BCAST(rmin_glf   ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   !CALL MPI_BCAST(rmaj_glf   ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   !CALL MPI_BCAST(q_glf      ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   !CALL MPI_BCAST(shat_exp   ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   !CALL MPI_BCAST(alpha_exp  ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   !CALL MPI_BCAST(elong_glf  ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   
   end_time = mpi_wtime()
   mpi_etimes(2) = mpi_etimes(2) + end_time-start_time
   
   start_time = mpi_wtime()
   jin = js_mpi12(myid)-1
   jout = js_mpi12(myid+1)-2
   mpi_counts(1) = jout-jin+1
   jmm = -1
   CALL callglf2d(    & ! INPUT------------------------------------
        jin,jout,     &
        jeigen_iglf,  & ! 0
        jroot_iglf,   & ! 8
        glf23_iglf,   & ! 1
        jshoot,       & ! 0
        jmm,          & ! 0
        jmaxm,        & ! nj-1
        itport_pt,    & ! dimension 5
        irotstab,     & ! 1
        te_glf,       &
        ti_glf,       &
        ne_glf,       &
        ni_glf,       & 
        ns_glf,       & 
        igrad,        & ! 0 
        idengrad_glf, & ! 2
        zpte_glf(1),  & ! not used since igrad = 0
        zpti_glf(1),  & ! not used since igrad = 0
        zpne_glf(1),  & ! not used since igrad = 0
        zpni_glf(1),  & ! not used since igrad = 0
        angrotp_exp,  & 
        egamma_exp,   &
        gamma_p_exp,  &
        vphim_glf,    &
        vparm_glf,    &
        vperm_glf,    &
        zeff_exp,     &
        btor_exp,     & ! constant
        ibtflag_glf,  & ! 1
        rho_glf,      & ! normalized
        arho_iglf,    & ! r(nj)
        grho1_glf,    & 
        grho2_glf,    &
        rmin_glf,     &
        rmaj_glf,     &
        rmajor_glf,   & ! constant
        zimp_glf,     & ! constant
        amassimp_glf, & ! constant
        q_glf,        &
        shat_exp,     &
        alpha_exp,    &
        elong_glf,    &
        amassgas_glf, & ! constant
        exbmult_glf,  & ! constant
        x_alpha_glf,  & ! constant 
        i_delay,      & ! 0
                      !
                      ! OUTPUT-----------------------------------
                      !                   
        diffnem,      & ! constat not used
        chietem,      & ! constat not used
        chiitim,      & ! constat not used
        etaphim,      & ! constat not used
        etaparm,      & ! constat not used
        etaperm,      & ! constat not used
        exchm,        & ! constat not used  
        diff_m,       & ! not used
        chie_m,       & ! ***************************
        chii_m,       & ! ***************************
        etaphi_m,     & ! ***************************
        etapar_m,     & ! not used
        etaper_m,     & ! not used
        exch_m,       & ! not used 
        egamma_m,     & ! not used
        egamma_d,     & ! not used
        gamma_p_m,    & ! not used
        anrate_m,     & ! not used
        anrate2_m,    & ! not used
        anfreq_m,     & ! not used
        anfreq2_m     & ! not used
        )
   end_time = mpi_wtime()
   mpi_etimes(4) = mpi_etimes(4) + end_time-start_time
   glf_load_buf = end_time-start_time
   
   start_time = mpi_wtime()
   
   IF(iloadbal.eq.1) THEN
      CALL dynamic_balance(glf_load_buf)
   ENDIF   
   
   !do i=js_mpi12(myid),js_mpi12(myid+1)-1
   !	  print *,'**h**',i,chii_m(i)
   !end do	  
   
   CALL mpi12_g_send(diff_m     ,4)
   CALL mpi12_g_send(chie_m     ,4)
   CALL mpi12_g_send(chii_m     ,4)
   CALL mpi12_g_send(etaphi_m   ,4)
   
   end_time = mpi_wtime()
   
   mpi_etimes(3) = mpi_etimes(3) + end_time-start_time

END SUBROUTINE