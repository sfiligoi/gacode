     MODULE rf
        USE param, only : krf,kj,kprim,krt,kzrf,kcrm
        implicit none  
!
        save
        character*8  :: no_rf ='no rf'  
        character(len=80),dimension(:) ::  extcurrf_id(krf),            &
                                 extqerf_id(krf), extqirf_id(krf)
        character(len =8), dimension(:)   ::  irfmodel(krf)
        character(len =8), dimension(:)   ::  rfmode(krf)
        CHARACTER *256 :: ech_input = 'none'

!
        integer,parameter :: krfrad =  kj
        integer,dimension(:)   ::                                        &
                   irf(krf), idamp(krf), nkfcd(krf),                     &
                   lifw(kprim), nihfw(kprim),nray(krf),   &
                   extcurrf_nj(krf), extqerf_nj(krf), extqirf_nj(krf)
        integer                                                          &
                  iswchfw, nfwsmth, jresonrf,beam_spec(4),               &
                  impath, necsmth, nzrffw, nrayptrt,                     &
                  irfplt,  irfech, ifb, lmode, ifbprof,nrfrad, iside,    &
                   jrf1min, jrf1max, jrf2min, jrf2max, jrfmin, jrfmax,   &
                   nhigh, navg, ichmod,active_rf_models,save_curray_input, &
	           nmodel
        integer  ::  case_no = 0
        real *8                                                        &
                wrfx, betalm, ylaunch, rfcur,gafsep,                   &
                totrtpe, totrtpi, totrtc, curray_fi,                   &
                totecpe, totecc, prfes, prfis ,relrf, xkpar, ykperp,   &
                relrf_pow
!
        real*8, dimension(:)  ::                                       &
             rfon(krf), rfoff(krf), rftime(krf), irfcur(krf),          &
             turnonp(krf), turnonc(krf),rfmodel_power_e(krf),          &
             rframp_timeup(krf),rframp_timedown(krf),                  &
             rfpow(krf), freq(krf), wrfo(krf), rfmodel_cd(krf),        &
             rfrad(krfrad), rfrow1(krfrad), rfmodel_power_i(krf),      &
             rfrow2(krfrad), rfrow(krfrad), omoder(krfrad),            &
             omodei(krfrad), xmoder(krfrad), xmodei(krfrad),           &
             xec(krf), zec(krf),                                       &
             thetec(krf), phaiec(krf), hlwec(krf), ratwec(krf),        &
             timrfp(5), rf1(krfrad),                                   &
             rf2(krfrad), rf3(krfrad), rf4(krfrad),                    &
             rnp(krf), alphaf(3),                                      &
             fbscrch(kj),  zrffw(kzrf),                                &
             pzrffw(kzrf), htsfw(kzrf), freqfw(kzrf),                  &
             rnpfw(kzrf), totfwpe(kzrf), totfwpi(kzrf),                &
             totfwc(kzrf),                                             &
             extcurrf(krf), extqerf(krf), extqirf(krf),                &
             extcurrf_amps(krf), extqerf_watts(krf),                   &
             extqirf_watts(krf),                                       &
             rf_ext_curtot(krf), rf_ext_qetot(krf),                    &
             rf_ext_qitot(krf),                                        &
             gamloss(krf) 
!
!
        real *8, dimension(:,:)   ::                                   &
             enrf(kj,krf), terf(kj,krf), tirf(kj,krf),                 &
             qrfes (kj,krf), qrfis(kj, krf),                           &
             extcurrf_rho(kj,krf), extqerf_rho(kj,krf),                &
             extqirf_rho(kj,krf), currfs(kj,krf),                      &
             extcurrf_curr(kj,krf), extqerf_qe(kj,krf),                &
             extqirf_qi(kj,krf),xntor(kzrf,krf), rpant(kzrf,krf),      &
             currf_save(kj)


! --- GENRAY related items
      real *8,dimension(:),allocatable :: rgenray,pgre,pgrit,pgrc  !size nj
      real *8,dimension(:),allocatable :: totgrps
      real *8,dimension(:,:),allocatable :: pgri
      real *8,dimension(:),allocatable :: charge_nc,dmass_nc  !size nspecgr
      real *8,dimension(:,:),allocatable :: en_nc,temp_nc  !size nj,nspecgr
      integer save_genray_io,nspecgr
      real *8  totgrpe, totgrpi, totgrc, totgrp, genray_fi
      character(len=128), dimension(:) :: genraydat(krf)
      character(len=128) :: gfilename
      character(len=1) :: char1,char2
      character(len=2) :: char_name1,char_name2
      character(len=5) :: char_name

! --- end GENRAY related items
     
! 
! --- curray related items:
!
         integer, parameter   :: kspec  = 11
         character *256 runid_curray,curray_file1,curray_file2,        &
                        curray_in_spawn,trxplout_spawn
         real *8 psistep,pkexpnt,epserr,epser1
         integer  igraph,iprint,idcur,indvar,ichois,igrill,modcd,      &
                  idmpsw,nminor,kalfa,nspect,irayiort,icurdrrt,        &
                  incrt,bmaxrt,nicsmth
         integer  indx_curray(krf),nthinrt(kcrm)
         integer,  dimension(:,:) :: nnkpolrt(krt,kcrm), nnkparrt(krt,kcrm)
         integer,dimension(:) :: nnkpol(krt),nnkpar(krt),ichoisrt(kcrm), &
                                 idmpswrt(kcrm)
         real *8,dimension(:),allocatable ::  rcurray,pwe,pwc,pwit !size kj
!
         real *8,dimension(:,:),allocatable :: pwi
         ! note allocatable arrays cant be in namelists
        real*8, dimension(:,:)  ::                                     &
             powersrt(krt,kcrm), anzinfrt(krt,kcrm),                   &
             anzsuprt(krt,kcrm),                    &
             anpinfrt(krt,kcrm), anpsuprt(krt,kcrm)                  
              
         real *8, dimension(:) :: powers(krt),                         &
                                  anzinf(krt),anzsup(krt),             &
                                  anpinf(krt),                         &
                                  anpsup(krt),heightrt(kcrm),          &
                                  thgrilrt(kcrm),islofart(kcrm),       &
                                  psi_startrt(kcrm),maxrefrt(kcrm)
         contains
!
         subroutine wrt_curray_in(nprim,nion,nimp,ptot0,ebkev,         &
                                  nalfa1,atmf1,azf1,freqcy,nthin,         &
                                 islofa,psi0,maxref,irayio,inc,bmax,   &
                                 idrive,nprofs,thgril,nminor,heigth,   &
                                 time_eqdsk)   
!---------------------------------------------------------------------------
!  creates curray input file, curray_in,  from information 
!  available in Onetwo.
!  Note that curray_in contains soem information which is also passed
!  in file trxpl.out. In particular curray only accepts one fast ion
!  distribution. In the case of beams and fusion alphas this becomes
!  probematic if we try to run strictly with curray_in (and not using 
!  trxpl.out) . In this file we always set nalfa,atmf,azf =0 (indicating no
!  fast ions) and then pass the information about fast ions to curray
!  through trxpl.out.
!
!
! ---------------------------------------------------------------HSJ--------
       USE param,only:  kion
       USE ename ,only : eqdskfilename,eqfile,eqdsk_tdem
       use io, only :eqdskin,ncrt,nout
       USE ions,only :  atw,atomno
          character (len = 64) eqdsk_name,old_eqdsk_name,eqdsk_name_prev
          character (len =9) :: curray_in ='curray_in'
          character (len = 256) :: command
          integer iprof,ioread,nprim,nspec,nimp,kboot,nion,             &
                  irayio,icurdr,nraypts,io12,nalfa,nthin,bmax,          &
                  islofa,maxref,inc,idrive,nminor,nprofs,lenge,         &
                  isimm,ISHELL,nalfa1
          integer,dimension(6) :: iounit = (/6, 12, 13, 14, 15, 16 /) 
          real *8,dimension(:) :: atm(kion),azi(kion)
          real *8 ptot0,heigth,ebkev,atmf,azf,freqcy,psi0,              &
                  thgril,time_eqdsk,atmf1,azf1
!
          integer nlim_eqd_t,np_eqd_t,nw_t,nh_t,npset,nspectl,j
          namelist /input/eqdsk_name,nprofs,nraypts,nthin,islofa,inc,   &
                         maxref,bmax,                                   &
                         psi0,psistep,iprof,pkexpnt,freqcy,epserr,      &
                         epser1,igraph,iprint,idcur,indvar,ioread,      &
                         ichois,igrill,thgril,modcd,idmpsw,             &
                         nprim,nspec,nimp,nminor,atm,azi,               &
                         kboot,kalfa,nalfa,isimm,nspect,ptot0,powers,   &
                         nnkpar,anzinf,anzsup,nnkpol,anpinf,anpsup,     &
                         heigth,irayio,icurdr,idrive,iounit,            &
                         ebkev,atmf,azf 
          
!
          nalfa =0   !see comments above
          atmf  =0
          azf   =0
!
!
        eqfile = eqdskfilename
        nspec = nion 
        nspectl = 0
        do j=1,krt
           if(powers(j) .gt. 0.0)nspectl = nspectl+1
        enddo
!        if(nspectl .ne. nspect)then
!           print *,'error , at least one toroidal wave number bin has zero power'
!           call STOP("wrt_curray_in, zeo power in bin",1)
!        endif
        nspect =nspectl
        isimm = nspect
        icurdr = icurdrrt
        if(eqdsk_tdem .ne. 'tdem' ) then
             lenge = LEN_TRIM(eqdskfilename)
             eqdsk_name  = eqdskfilename(1:lenge)
        else
           call wrt_tdem_eqdsk(time_eqdsk,eqdsk_name)
           eqfile = eqdsk_name
           print *,'eqdsk file created in sub rf_mhddat'
        endif
        print *,'eqdsk_name =', eqdsk_name
        print *,'eqfile =',eqfile
        print *,'eqdsk_tdem =',eqdsk_tdem
        if (eqfile .eq. 'none')  eqdsk_name =  'eqdskin'
          atm(:) = atw(:)
          azi(:) = atomno(:)
          kalfa = 0           ! no fast ion slowing down ditribution effects
          kboot  = 0          ! no bootstrap current calcs in curray
          iprof  = 1          ! input profiles in file trxpl.out
          ioread = 1          ! means curray called from Onetwo
                              ! NOTE: curray no longer reads  raytrin
                              ! ioread =1 still required to write curray output
                              ! file raytrout however
          nraypts = nrayptrt  ! input as nrayptrt in inone 
!
          if(save_curray_input .eq. 1 .and. case_no .eq. 0) then  
              case_no = case_no + 1            !initialization case
              old_eqdsk_name = eqdsk_name
              call save_cfile(curray_in,case_no,eqdsk_name)
              curray_in_spawn = curray_file1
              !save unique eqdsk name in curray_file2:
              command = 'cp '//eqdsk_name//' '//curray_file2
              if (ISHELL (command) .lt. 0)   &
                call STOP ('sub wrt_curray: failure of spawned cp  command', 67)
              eqdsk_name = '"'//curray_file2(1:LEN_TRIM(curray_file2))//'"'
              eqdsk_name_prev =  eqdsk_name
          elseif(save_curray_input .eq. 1 )then
              case_no = case_no+1
              if(old_eqdsk_name(1:LEN_TRIM(old_eqdsk_name)) ==  &
                            eqdsk_name(1:LEN_TRIM(eqdsk_name)))  then
                   call save_cfile(curray_in,case_no) !eqdsk same as previous
                   curray_in_spawn = curray_file1
                   eqdsk_name = eqdsk_name_prev
              else
                   old_eqdsk_name = eqdsk_name !new eqdsk not yet saved
                   call save_cfile(curray_in,case_no,eqdsk_name)
                   curray_in_spawn = curray_file1
!
                   !save unique eqdsk name in curray_file2:
                   command = 'cp '//eqdsk_name//' '//curray_file2
                   if (ISHELL (command) .lt. 0)   &
                   call STOP ('sub wrt_curray: failure of spawned cp  command', 67)
                   eqdsk_name = '"'//curray_file2(1:LEN_TRIM(curray_file2))//'"'
                   eqdsk_name_prev = eqdsk_name 
              endif
          else
             curray_in_spawn ='curray_in'
             eqdsk_name = '"'//eqdsk_name(1:LEN_TRIM(eqdsk_name))//'"'
             eqdsk_name_prev = eqdsk_name 
          endif
          eqdsk_name = ADJUSTL(eqdsk_name)
!
          !write the curray_in namelist to file curray_in_spawn
          call getioun(io12,42)
          open  (unit = io12, file = curray_in_spawn, status = 'UNKNOWN')
!
          write (unit = io12, fmt = '(3x, a)') runid_curray
          write (unit = io12, nml = input)
!
          close (unit = io12)
          call giveupus(io12)
          print *,'curray_in =',curray_in
!
          return
         end subroutine wrt_curray_in
!
!
      subroutine wrt_trxpl_out(nrho,nion,r)
! ------------------------------------------------------------------------------
! creates curray input file trxpl.out for curray from information in Onetwo:
! ----------------------------------------------------------------HSJ-----------
!  Curray description of quantitites put out by this rotuine:
!  'trxpl.out' : input plasma profile parameters from TRANSP runs.
!                called when iprof = 1
!
!  nrho    : number radial grid points
!
!  nspece   : nspec+1    number species including electrons,excluding fast ions
!                        (see below)
!
!  xrho(i) : radial grid points in normalized rho [=sqrt(phi)]
!
!  aze     : electron charge number (not used in CURRAY)
!  atme    : electron mass number (not used in CURRAY)
!  ne_tr(i) : electron density profile (m^-3)
!  te_tr(i) : electron temoerature profile (keV)
!
!  azi(k)  : ion charge number
!  atm(k)  : ion mass number
!  ni_tr(i,k) : ion density profile (m^-3)
!  ti_tr(i,k) : ion temperature profile (keV)
!
!  nfast      : number of fast ion species (not passed to trxpl.out)
! -----------------------------------------------------------------------------
          USE param,only:  kion
          USE fusion,only : enalp,walp,nalp_thresh
          USE ions,only :  atw,atomno
          USE soln, only : ene,en,te,ti
          USE nub, only : neg_ion_source,ibion,nbeams
          USE nub2,only : enb,nb_thresh,wb,enbmin_curray,tmin_curray, &
                        enbmin
          USE numbrs , only : nj
          integer nrho,nspec,nin,nfast,ik,jb,ic,jup,     &
                  nspece,nion,j,ipres,inc_alpha,nbeamlc
          real *8, dimension(:) :: r(nrho),xrho(nrho),            &
                                   ni_tr(nrho),                   &
                                   ne_tr(nrho),                   &
                                   ti_tr(nrho),el_tr(nrho),elp_tr(nrho)  
          real *8,dimension(:) :: atm(kion),azi(kion)
          real *8 aze,atme,atmb,azb ,ppsum,densum
          character (len = 9) :: trxplout = 'trxpl.out'
!
!
!
!
           !determine nfast, the number of effective fast ion species.
           ! an effective species is a fast ion species of a particular
           ! type, particular birth energy, particular injector source
           ! nfast = nbeamlc + nalpha
           !if beam is off nbeamlc = 0 
!
           nbeamlc = 0
           if(ibion .gt. 0 )then       ! check the beams
              atmb = atw(ibion) 
              azb = 1.0
!              if(beam_spec(1) .eq. -1)then       eliminated this part 7/12/04 HSJ
!                 nbeamlc = 1  need to checkk enbeam here 
!              else
                 do jb = 1,nbeams
                    jup = 3
                    if(neg_ion_source(jb) .gt. 0 )jup = 1
                    do ic = 1,jup
                        ipres =0
                        do j =1 ,nrho 
                           ! above threshold density means beam present :
                           if(enb(j,ic,jb) .gt. nb_thresh ) ipres  = 1
                        enddo
                        if(ABS(beam_spec(ic)) .eq. 1 .and. ipres .eq. 1)then
                         nbeamlc = nbeamlc + 1    !add beam components as  species
                        endif
                    enddo
                 enddo
!              endif
           endif
!
!
!
           nfast = nbeamlc
           ! same for fusion alphas:
           ipres = 0
           inc_alpha = 0
           do j=1,nrho
              if(enalp(j) .gt. nalp_thresh)ipres = 1
           enddo
           if(beam_spec(4) .eq. 1 .and. ipres .eq. 1)then
                 nfast = nfast +1
                 inc_alpha = 1
           endif
           !due to the above nfast (and hence nspec) may vary during
           !the simulation period
!
!
          call getioun(nin,42)
          open (unit = nin, file = trxplout, status = 'UNKNOWN')
           nspec = nion
           nspece = nspec+1+nfast
           atm(:) = atw(:)
           azi(:) = atomno(:)
           aze =0.0
           atme =0.0
           elp_tr(:) = 0.0                         !apparently not used in curray 
           xrho(1:nrho) = r(1:nrho)/r(nrho)                   
           ne_tr(1:nrho) = ene(1:nrho)*1.e+6
           write(nin,*) nrho
           write(nin,*) nspece
           write(nin,*)
           write(nin,8201) xrho(:nrho)
           write(nin,*)
           write(nin,'(f19.12)') aze
           write(nin,'(f19.12)') atme
           write(nin,*)
           write(nin,8201) ne_tr(:nrho)  ! m-3
           write(nin,*)
           write(nin,8201) te(:nrho)  ! keV
           do ik=1,nspec                  ! loop over thermal ions, nspec = nion
              write(nin,*)
              write(nin,'(f19.12)') azi(ik)
              write(nin,'(f19.12)') atm(ik)
              write(nin,*) 
              ni_tr(1:nrho) = en(1:nrho,ik)*1.e6
              write(nin,8201) ni_tr(:nrho)  ! m-3
              write(nin,*)
              write(nin,8201) ti(:nrho)  ! keV, all thermal species have same Ti
           enddo
           nspec = nion+1   !count up the species, nion thermal +1 electron species
!
!
           if_fast_ions: if ( nfast > 0)then  !beam and fusion contributions:
!
              
              if_beam_ions:        if(nbeamlc .gt. 0)then
                 atmb = atw(ibion) 
                 azb = 1.0
                 if(beam_spec(1) .gt. -1 )then
                    do jb = 1,nbeams
                       jup = 3
                       if(neg_ion_source(jb) .gt. 0 )jup = 1
                       do ic = 1,jup
                          ipres  = 0
                          do j =1 ,nrho 
                             ! above threshold density means beam present :
                             if(enb(j,ic,jb) .gt. nb_thresh ) ipres  = 1
                          enddo
                          if(beam_spec(ic) .eq. 1 .and. ipres .eq. 1)then
                             ni_tr(1:nrho) = enb(1:nrho,ic,jb) * 1.e6
                             ti_tr(1:nrho) = 0.667* 6.2415097e+15 *   &
                                wb(1:nrho,ic,jb)/enb(1:nrho,ic,jb)   !kev 
                             write(nin,*)
                             write(nin,'(f19.12)') azb
                             write(nin,'(f19.12)') atmb
                             write(nin,*) 
                             write(nin,8201) ni_tr(:nrho)  ! m-3
                             write(nin,*)
                             write(nin,8201) ti_tr(:nrho)  ! keV
                          endif
                       enddo
                    enddo
                 else       !beam_spec(1) = -1 ===> combine all beamlets
                   do jb = 1,nbeams
                      jup = 3
                      if(neg_ion_source(jb) .gt. 0 )jup = 1
                      do j =1,nrho
                       ppsum  =0.0
                       densum =0.0
                          do ic = 1,jup
                             ppsum = ppsum +wb(j,ic,jb)    !sum partial pressures
                             densum = densum + enb(j,ic,jb)
                          enddo !beam energies
                       ni_tr(j) = densum
                       ti_tr(j) = ppsum
                       enddo       !nrho loop
                       ti_tr(:) = 0.667* 6.2415097e+15 *ti_tr(:)/ni_tr(:)  !kev 
                       ni_tr(:) = ni_tr(:) * 1.e6
                       do j = 1,nrho   !beams can be aimed so that there are regions
                         !where there are no  fast ions. this will
                         !cause curray to choke:
                         if(ti_tr(j) .le. 0.0d0)ti_tr(j) =  tmin_curray
                         if(ni_tr(j) .le. 1.1*jb*enbmin*1.e6)ni_tr(j) = &
                                                       enbmin_curray*1.e6
                       enddo
                       write(nin,*)
                       write(nin,'(f19.12)') azb
                       write(nin,'(f19.12)') atmb
                       write(nin,*) 
                       write(nin,8201) ni_tr(:nrho)  ! m-3
                       write(nin,*)
                       write(nin,8201) ti_tr(:nrho)  ! keV
                     enddo    !nbeams
                 endif        !beam_spec
              endif if_beam_ions
!
!
!
              ! fusion alphas:
              if(beam_spec(4) .eq. 1 .and. inc_alpha .eq. 1 )then
                 atmb  = 4.0
                 azb   = 2.0
                 ni_tr(1:nrho) = enalp(1:nrho) * 1.e6
                 ti_tr(1:nrho) = 0.667* walp(1:nrho)/enalp(1:nrho)
                 write(nin,*)
                 write(nin,'(f19.12)') azb
                 write(nin,'(f19.12)') atmb
                 write(nin,*) 
                 write(nin,8201) ni_tr(:nrho)  ! m-3
                 write(nin,*)
                 write(nin,8201) ti_tr(:nrho)  ! keV
              endif
!
           endif if_fast_ions           
!
!
!
           write(nin,*)
           write(nin,8201) elp_tr(:nrho)
!
 8201      format(5(1x,e16.9))  
           close (unit = nin)
           call giveupus(nin)    
 
           if(save_curray_input .eq. 1) then             
               call save_cfile(trxplout,case_no)
                trxplout_spawn = curray_file1
           else
                trxplout_spawn = 'trxpl.out'
           endif
        
!
          return
      end subroutine wrt_trxpl_out
!
!
!
      subroutine save_cfile(file1,case_no,file2)
!  ------------------------------------------------------------------
!    save files 1 and optionally file 2 by appending time stamp
!    and case number to name
!  -------------------------------------------------------------------
       USE solcon,only : time
       integer ,intent(in) :: case_no
       character(len =*) ,intent(in) :: file1
       character(len =*) ,intent(in),optional :: file2
       character(len = 9) :: cht 
       character(len = 5) :: chi
       character(len=18)  :: storesp
!
!
!
       write(storesp,'(i5)')case_no
       read(storesp,'(a)')chi
       write(storesp,'(f14.6)')time
       read(storesp,'(a)')cht
       chi = ADJUSTL(chi)     !move leading blanks to end
       cht = ADJUSTL(cht)
!
       curray_file1 = file1(1:LEN_TRIM(file1))//'_'//cht(1:LEN_TRIM(cht))    &
                                       //'_'//chi(1:LEN_TRIM(chi))
       curray_file1 = ADJUSTL(curray_file1)
       print *,'curray_file1 =',curray_file1
       if(present(file2))then
           curray_file2 = file2(1:LEN_TRIM(file2))//'_'//cht(1:LEN_TRIM(cht))    &
                                           //'_'//chi(1:LEN_TRIM(chi))
           curray_file2 = ADJUSTL(curray_file2)
           print *,'curray_file2 =',curray_file2
       endif
!
      return
      end subroutine save_cfile
!
!
!
! --- end Curray related items
!
!
!--- external RF heating and current drive profiles are defined by
!                 extcurrf, extqerf, extqirf
!                 extcurrf_amps, extqerf_watts, extqirf_watts
!                 extcurrf_id, extqerf_id, extqirf_id
!                 extcurrf_nj, extqerf_nj, extqirf_nj
!                 extcurrf_rho(kj), extqerf_rho, extqirf_rho
!                 extcurrf_curr(kj), extqerf_qe, extqerf_qi
!                 rf_ext_curtot, rf_ext_qetot, rf_ext_qitot
!
!--- see subroutine INIT for definitions, and the subroutines
!                 GET_EXTERNAL_RFCUR,
!                 GET_EXTERNAL_RFQE,
!                 GET_EXTERNAL_RFQI to see what is done ----------- HSJ
!


   END MODULE rf
 
