!---------------------------------------------------------------------------------------------
! include file nbnamelist_post_201112.inc
! Use this module with nubeam version 201112  and P_nfreya
! allows inclusion of namelist in varios places
!------------------------------------------------------------------------------HSJ-2-14-2012--


  NAMELIST/nbdrive_naml/                                                                              & 
      abeama                  ,      adiff_0                 ,      adiff_a                 ,  &
      adiff_ntime             ,      adiff_time              ,      adiff_xpin              ,  &
      adiff_xpout             ,      aimp                    ,      asrd                    ,  &
      blk_mpi_bcast_int       ,      blk_mpi_bcast_r8        ,      bmwidra                 ,  &
      bmwidza                 ,      bsrd                    ,      cxpcon                  ,  &
      cxsplt                  ,      divra                   ,      divza                   ,  &
      dn0out                  ,      dtn                     ,      dt_acc                  ,  &
      dxbsmoo                 ,      dxbsm_nc                ,      ebdmax                  ,  &
      edbfac                  ,      einja                   ,      erngfi                  ,  &
      fbemax                  ,      fbemin                  ,      fbltim                  ,  &
      fbtrap_depth            ,      fdtnxy                  ,      ffulla                  ,  &
      fhalfa                  ,      foclra                  ,      foclza                  ,  &
      fporcelli               ,      fppcon                  ,      frac_depmax             ,  &
      frac_depmin             ,      frac_dep_lim            ,      frac_orbrr              ,  &
      fshper                  ,      fshwid                  ,      fvpvmn                  ,  &
      fvpvmx                  ,      gflr_ll                 ,      gflr_min                ,  &
      gflr_op                 ,      gflr_rl                 ,      gflr_xv                 ,  &
      goocon                  ,      levmod_halo             ,      lev_nbidep              ,  &
      lmidtube                ,      lunnbx                  ,      lunres                  ,  &
      mrstrt                  ,      nbapsh2                 ,      nbapsha                 ,  &
      nbbcal                  ,      nbbcx_avg               ,      nbbcx_bb                ,  &
      nbbox                   ,      nbeam                   ,      nbebox                  ,  &
      nbetube                 ,      nbfallgr                ,      nbsbox                  ,  &
      nbshapa                 ,      nbstube                 ,      nbtube                  ,  &
      nchdvp                  ,      nclass                  ,      ncx0                    ,  &
      ndep0                   ,      ndep0_max               ,      ndepbox                 ,  &
      ndeptube                ,      ndep_set_beam           ,      ndep_set_ien            ,  &
      ndtorb                  ,      nerngfi                 ,      nfbon_species           ,  &
      nfbon_vpvopt            ,      ngoocon_ebin            ,      ngoocon_rbin            ,  &
      ngoocon_vpvbin          ,      ngradd_opt              ,      ngyro                   ,  &
      nkdifb                  ,      nlbbcx                  ,      nlbcde                  ,  &
      nlbcoh                  ,      nlbcpa                  ,      nlbdat                  ,  &
      nlbeamcx                ,      nlbflr                  ,      nlbfpp                  ,  &
      nlbgflr                 ,      nlbout                  ,      nlbox                   ,  &
      nlco                    ,      nlcprb                  ,      nldep0_gather           ,  &
      nlebei                  ,      nlfatom                 ,      nlfbmflr                ,  &
      nlfbon                  ,      nlfdep                  ,      nlhvion                 ,  &
      nlminsv                 ,      nlorbo                  ,      nlpsirz                 ,  &
      nlsawb                  ,      nlsawf                  ,      nlsym2b                 ,  &
      nltest_output           ,      nltrk_dep0              ,      nlusf3                  ,  &
      nlusfa                  ,      nlusfp                  ,      nlusft                  ,  &
      nmcurb                  ,      nmimp                   ,      nmix_kdsaw              ,  &
      nonlin                  ,      nper_cx                 ,      nptclf                  ,  &
      nptclh                  ,      nptcls                  ,      nptcl_max               ,  &
      nrip                    ,      nsdbgb                  ,      nseed                   ,  &
      nsegtube                ,      nsigexc                 ,      nth0                    ,  &
      nubeam_dt               ,      nxbox                   ,      nybox                   ,  &
      nznbma                  ,      nznbme                  ,      nzones                  ,  &
      nzone_fb                ,      only_io                 ,      orbrzv_option           ,  &
      orbrzv_rf_option        ,      orbrzv_toric_frnm       ,      orbrzv_toric_prftot0    ,  &
      orbrzv_zzerr_con        ,      phitube                 ,      pinja                   ,  &
      plfhe3                  ,      plfhe4                  ,      plfsp                   ,  &
      plfst                   ,      quasi_check             ,      rapedg2                 ,  &
      rapedga                 ,      ref_namelist            ,      rhotube                 ,  &
      rtcena                  ,      rtube                   ,      sawflag                 ,  &
      taurip                  ,      tboffa                  ,      tbona                   ,  &
      tfshof                  ,      tfshon                  ,      thetatube               ,  &
      wghta                   ,      xboxhw                  ,      xbzeta                  ,  &
      xcfafus                 ,      xcfanbi                 ,      xdatsf3                 ,  &
      xdatsfa                 ,      xdatsfp                 ,      xdatsft                 ,  &
      xdepmod                 ,      xdfafus                 ,      xdfanbi                 ,  &
      xefafus                 ,      xefanbi                 ,      xfishmax                ,  &
      xfishmin                ,      xl1tube                 ,      xl2tube                 ,  &
      xlbapa                  ,      xlbapa2                 ,      xlbox1                  ,  &
      xlbox2                  ,      xlbtna                  ,      xpand_nptcl             ,  &
      xplasma_in_memory       ,      xrapoff2                ,      xrapoffa                ,  &
      xswfrac_allfast         ,      xswfrac_beam            ,      xswfrac_fusn            ,  &
      xybapa                  ,      xybsca                  ,      xzapoff2                ,  &
      xzapoffa                ,      xzbeama                 ,      xzetatube               ,  &
      xzimp                   ,      xzpedg2                 ,      xzpedga                 ,  &
      yboxhw                  ,      ytube                      
