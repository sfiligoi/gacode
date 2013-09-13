  MODULE source_terms_gcnmp
    USE nrtype, ONLY : DP,I4B
    IMPLICIT NONE
    SAVE
    PUBLIC
    REAL(DP),DIMENSION(:,:),ALLOCATABLE :: stsource,sion,srecom,scx,  &
         dudtsv,sbcx,siadd,s
    REAL(DP),DIMENSION(:),ALLOCATABLE ::      sbeame,qohm ,           &
                     qione,qfuse,qioni,qcx,qe2d,qi2d,qfusi,qbfuse,    &
                     qbfusi,qbeame,qbeami,qrfe,qrfi,                  &
                     sscxl,curbe,curbi,sdfuse,sdtfuse,qtfuse,         &
                     curbet,sprbeame,sprbeami,ssprcxl,qtfusi,         &
                     sprcxre,spreimpt,sprcx,sprcxe,sprcxree,          &
                     spreimpe,smagtorque,flxadd,sngas,xnue,           &
                     xnue_frozen,qdimpl,sione,stfuse,scurdri,         &
                     etap,scurfast,qbeami_rot,qbeame_rot,             &
                     sbfuse,sfus,qrad,prad,prad_imp,qexch,            &
                     storque,spr2d,brems_tot,                         &
                     qexsav,qieneo,dneo,xkeneo,xkineo,                &
                     chiineo,xden,xdet,qdelt,spellet,                 & 
                     spaleo,qepaleo, qipaleo

   REAL(DP),DIMENSION(:,:),ALLOCATABLE :: s2d,xdin,xdit,sbeam,        &
                                          brems_nions

   REAL(DP) beam_th_mult, stfuse_mult, sbfuse_mult,qrad_mult,         &
            qfuse_mult_e,qfuse_mult_i,qrf_mult_e,qrf_mult_i,          &
            qbeam_mult_e,qbeam_mult_i,qrad_tot

   INTEGER(I4B) k
   LOGICAL imp_called


   CONTAINS

     SUBROUTINE allocate_source_arrays

      USE grid_class,                    ONLY : nj

      USE plasma_properties,             ONLY : pwrden,prtcl_src,pellet

      USE vector_class,                  ONLY : get_values

      USE ions_gcnmp,                    ONLY : nimp,nion,dzdtim
      
      USE solcon_gcnmp,                  ONLY : ntot

      USE common_constants,              ONLY : zeroc,kevpjou

      USE vector_class,                  ONLY : zero_vector
      IMPLICIT NONE
  
      INTEGER(I4B) jj

        ! resistivity, ohm m on full grid (1:nj)
      IF(.NOT. ALLOCATED(etap))ALLOCATE(etap(nj)) 
      etap(:) = zeroc


      !electron ion radiation
      IF(.NOT. ALLOCATED(qrad))    ALLOCATE(qrad(nj))
      qrad(:) =  get_values(pwrden%qrad) * kevpjou  ! kev/(m**3 sec),qrad_mult is applied in source
      IF(.NOT. ALLOCATED(prad))ALLOCATE(prad(nj))
        prad(:)   = zeroc  
      IF(.NOT. ALLOCATED(prad_imp))ALLOCATE(prad_imp(nj))
        prad_imp(:)   = zeroc

     IF( .NOT. ASSOCIATED(pwrden%brems_nions))THEN
          ALLOCATE(pwrden%brems_nions(nion))
          DO jj =1,nion
             pwrden%brems_nions(jj)  = zero_vector(nj)
          ENDDO

     ENDIF

      IF(.NOT. ALLOCATED(brems_nions))ALLOCATE(brems_nions(nj,nion))
        DO k =1,nion
           brems_nions(:,k) = pwrden%brems_nions(k)%data(:)*kevpjou  ! kev/(m^3 sec)
        ENDDO





        !electron power density due to recombination and 
        !impact ionization, watts/meter**3 ==> kev/(m**3 sec)
      IF(.NOT. ALLOCATED(qione))ALLOCATE(qione(nj))
        qione(1:nj) = get_values(pwrden%qione)*kevpjou

        !ion power density due to recombination and 
        !impact ionization, watts/meter**3  ==> kev/(m**3 sec)
      IF(.NOT. ALLOCATED(qioni))ALLOCATE(qioni(nj))
        qioni(1:nj) = get_values(pwrden%qioni)*kevpjou


        !ion power density due to neutral-ion charge exchange, watts/meter**3
      IF(.NOT. ALLOCATED(qcx))ALLOCATE(qcx(nj))
        qcx(1:nj)  = get_values(pwrden%qcx)*kevpjou 

        !2d ion heating, watts/meter**3
      IF(.NOT. ALLOCATED(qe2d))ALLOCATE(qe2d(nj))
        qe2d(1:nj) = get_values(pwrden%qe2d)*kevpjou

        !2d ion heating, watts/meter**3
      IF(.NOT. ALLOCATED(qi2d))ALLOCATE(qi2d(nj))
        qi2d(1:nj) = get_values(pwrden%qi2d)*kevpjou

       IF(.NOT. ALLOCATED(qohm))ALLOCATE(qohm(nj))  
                 !ohmic power density, watts/meter**3
                 qohm(1:nj) = get_values(pwrden%qohm)*kevpjou


      IF(.NOT. ALLOCATED(qexch))ALLOCATE(qexch(nj))
         qexch(:) = get_values(pwrden%qexch)       
         !anomalous e-i energy exchange not implemented 
         !drift balloning model for example could be used
         !note that pwrden%qexch = 0 in rw_iterdb

        !fusion electron heating, watts/meter**3
      IF(.NOT. ALLOCATED(qfuse))ALLOCATE(qfuse(nj))
        qfuse(1:nj) = get_values(pwrden%qfuse)*kevpjou
      IF(.NOT. ALLOCATED(qtfuse))ALLOCATE(qtfuse(nj))
        qtfuse(1:nj) = zeroc

        !fusion ion heating, watts/meter**3
      IF(.NOT. ALLOCATED(qfusi))ALLOCATE(qfusi(nj))
        qfusi(1:nj) = get_values(pwrden%qfusi)*kevpjou
      IF(.NOT. ALLOCATED(qtfusi))ALLOCATE(qtfusi(nj))
        qtfusi(1:nj) = zeroc

        !beam fusion electron heating, watts/meter**3
      IF(.NOT. ALLOCATED(qbfuse))ALLOCATE(qbfuse(nj))
      IF(.NOT. ALLOCATED(qbfusi))ALLOCATE(qbfusi(nj))
        qbfuse(1:nJ) = get_values(pwrden%qbfuse)*kevpjou        
        qbfusi(1:nj) = get_values(pwrden%qbfusi)*kevpjou


        !power to elec.,resp ions  from beam, watts/meter**3
        IF(.NOT. ALLOCATED(qbeame))ALLOCATE(qbeame(nj))
        qbeame(1:nj) = get_values(pwrden%qbeame) *kevpjou
        IF(.NOT. ALLOCATED(qbeami))ALLOCATE(qbeami(nj))
        qbeami(1:nj) = get_values(pwrden%qbeami) *kevpjou

        !qrfe, RF electron,resp ion  heating, watts/meter**3
      IF(.NOT. ALLOCATED(qrfi))ALLOCATE(qrfi(nj))
      IF(.NOT. ALLOCATED(qrfe))ALLOCATE(qrfe(nj))
        qrfe(1:nj) = get_values(pwrden%qrfe) *kevpjou
        qrfi(1:nj) = get_values(pwrden%qrfi) *kevpjou

      IF(.NOT. ALLOCATED(sione))ALLOCATE(sione(nj))
        sione(:)  = zeroc  ! calculated in sub source

      IF(.NOT. ALLOCATED(qdimpl))ALLOCATE(qdimpl(nj))
        qdimpl(:)  = zeroc 


      IF( .NOT. ALLOCATED(srecom))ALLOCATE(srecom(nj,nion))
        DO k =1,nion
           srecom(:,k) = prtcl_src%srecom(k)%data(:)
        ENDDO

      IF(.NOT. ALLOCATED(spellet))ALLOCATE(spellet(nj))
        spellet(1:nj) = get_values(prtcl_src%spellet) ! load from  iterdb file input
                                                 ! recalc in source if called for
      IF(.NOT. ALLOCATED(stfuse))ALLOCATE(stfuse(nj))
        stfuse(:)  = get_values(prtcl_src%stfuse)   ! t burnup rate

      IF(.NOT. ALLOCATED(sdfuse))ALLOCATE(sdfuse(nj))
        sdfuse(:)  = get_values(prtcl_src%stfuse)   ! d burnup rate (set to input t burnup rate on entry)

      IF(.NOT. ALLOCATED(sdtfuse))ALLOCATE(sdtfuse(nj))
        sdtfuse(:)  =  zeroc   ! d burnup rate (set to input t burnup rate on entry)

      IF(.NOT. ALLOCATED(sbfuse))ALLOCATE(sbfuse(nj))
        sbfuse(:)  = get_values(prtcl_src%sbfuse)

      IF(.NOT. ALLOCATED(sfus))ALLOCATE(sfus(nj))
        sfus(:)   = stfuse(:) + sbfuse(:)
        ! sfus = stfuse  + sbfuse (contribution from thermal and fast ion fusion)
      
      ! Paleoclassical term 1/m^3, should be subtracted from sources
      IF(.NOT. ALLOCATED(spaleo))ALLOCATE(spaleo(nj))
      IF(.NOT. ALLOCATED(qepaleo))ALLOCATE(qepaleo(nj))
      IF(.NOT. ALLOCATED(qipaleo))ALLOCATE(qipaleo(nj))
      spaleo = zeroc  ; qepaleo = zeroc  ; qipaleo = zeroc

      IF(.NOT. ALLOCATED(qexch))ALLOCATE(qexch(nj))
        qexch(:)   = zeroc  !anomalous electron ion eng exchange
      IF(.NOT. ALLOCATED(qbeami_rot))ALLOCATE(qbeami_rot(nj))
      IF(.NOT. ALLOCATED(qbeame_rot))ALLOCATE(qbeame_rot(nj))
        qbeami_rot(:) = zeroc      ! energy sources due to rotation
        qbeame_rot(:) = zeroc
  

      IF(.NOT. ALLOCATED(siadd))ALLOCATE(siadd(nj,2))
        siadd(:,1) =  zeroc
        siadd(:,2) =  zeroc
      IF(.NOT. ALLOCATED(sscxl))ALLOCATE(sscxl(nj))
        sscxl(:)  = zeroc


      IF(.NOT. ALLOCATED(curbe))ALLOCATE(curbe(nj))
        curbe(:)  = zeroc !electron drag correction to curbi
      IF(.NOT. ALLOCATED(curbi))ALLOCATE(curbi(nj))
        curbi(:)  = zeroc ! unshielded bema ion current
        IF(.NOT. ALLOCATED(curbet))ALLOCATE(curbet(nj))
        curbet(:) = zeroc  ! trapped electron correction to curbe
        ! note that curbi,curbe,curbet are not used, only the total
        ! beam driven current is needed here.


      IF(.NOT. ALLOCATED(storque))ALLOCATE(storque(nj))
        storque(:) = zeroc
      IF(.NOT. ALLOCATED(spr2d))ALLOCATE(spr2d(nj))
        spr2d(:) = zeroc 
      IF(.NOT. ALLOCATED(s2d))ALLOCATE(s2d(nj,nion))
        s2d(:,:) = zeroc 

      IF(.NOT. ASSOCIATED(dzdtim))ALLOCATE(dzdtim(nj,nimp))
        dzdtim(:,:) = zeroc 


      IF(.NOT. ALLOCATED(scurdri))ALLOCATE(scurdri(nj))
        scurdri(:) = zeroc
      IF(.NOT. ALLOCATED(sprbeame))                                      &
            ALLOCATE(sprbeame(nj),sprbeami(nj),ssprcxl(nj),                &
                     sprcxre(nj),spreimpt(nj),sprcx(nj),                   &
                     sprcxe(nj),sprcxree(nj),spreimpe(nj),smagtorque(nj))
 

        IF(.NOT. ALLOCATED(flxadd))ALLOCATE(flxadd(2))
        IF(.NOT. ALLOCATED(sngas))ALLOCATE(sngas(2))
        IF(.NOT. ALLOCATED(xnue))ALLOCATE(xnue(nj),xnue_frozen(nj))
        IF(.NOT. ALLOCATED(s))ALLOCATE(s(ntot,nj))
        IF(.NOT. ALLOCATED(scurfast))ALLOCATE(scurfast(nj))
        IF(.NOT. ALLOCATED(qdelt))ALLOCATE(qdelt(nj))
     
        IF(.NOT. ALLOCATED(dneo))ALLOCATE(dneo(nj-1))
        IF(.NOT. ALLOCATED(xkeneo))ALLOCATE(xkeneo(nj-1))
        IF(.NOT. ALLOCATED(xkineo))ALLOCATE(xkineo(nj-1))
        IF(.NOT. ALLOCATED(chiineo))ALLOCATE(chiineo(nj-1))
        IF(.NOT. ALLOCATED(qieneo))ALLOCATE(qieneo(nj-1))


        s(:,:)      = zeroc     ; scurfast(:)    = zeroc
        sprbeame(:) = zeroc     ; sprcx(:)       = zeroc
        sprbeami(:) = zeroc     ; ssprcxl(:)     = zeroc 
        sprcxre(:)  = zeroc     ; spreimpt(:)    = zeroc 
        sprcxe(:)   = zeroc     ; sprcxree(:)    = zeroc
        spreimpe(:) = zeroc     ; smagtorque(:)  = zeroc
        xnue(:)     = zeroc     ; xnue_frozen(:) = zeroc
        qdelt(:)    = zeroc     ; dneo(:)        = zeroc
        xkeneo(:)   = zeroc     ; xkineo(:)      = zeroc
        chiineo(:)  = zeroc     ; qieneo(:)      = zeroc   

 
     END SUBROUTINE allocate_source_arrays

  END MODULE source_terms_gcnmp
