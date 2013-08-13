       subroutine writeglf2d(         !all these are inputs to callglf2d
     &   leigen, nroot, iglf          !execept litport,nj
     & , jshoot, jmm, jmaxm, itport_pt
     & , irotstab, te_m, ti_m, ne_m, ni_m, ns_m
     & , igrad, idengrad, zpte_in, zpti_in, zpne_in, zpni_in
     & , angrotp_exp, egamma_exp, gamma_p_exp, vphi_m, vpar_m, vper_m
     & , zeff_exp, bt_exp, bt_flag, rho
     & , arho_exp, gradrho_exp, gradrhosq_exp
     & , rmin_exp, rmaj_exp, rmajor_exp, zimp_exp, amassimp_exp
     & , q_exp, shat_exp, alpha_exp, elong_exp, amassgas_exp
     & , alpha_e, x_alpha, i_delay,lprint,litport,jpdm1)

c-------------------------------------------------------------------------
c       subroutine writes glf namelist file  nlglf
c       which will be used to get input into the stand alone glf23 code
c       The stand alone code uses mpi whcih is otherwise not possible
c       in Onetwo at this time.
c-------------------------------------------------------------HSJ-2/25/03
      implicit none
         character(len =12) ::  glf_namelst = 'glf_namelist'
         integer,parameter :: jpd = 50
         integer,parameter :: litport1 = 5
         integer
     .           leigen, lprint, nroot, iglf, jshoot, jmm, jmaxm,   
     .            irotstab,i_delay,igrad, idengrad, bt_flag,iounit,
     .           litport,iostat,jpdm1,i
         integer, dimension(:) :: itport_pt(1:litport1)

         real *8                                                       
     .        zpte_in, zpti_in, zpne_in, zpni_in,drho,bt_exp, arho_exp, 
     .        rmajor_exp, zimp_exp, amassimp_exp, amassgas_exp,        
     .        alpha_e, x_alpha              

         real *8, dimension(:) ::     !jpd =nj-1,allocated length in cray102.f
     .        te_m(0:jpd), ti_m(0:jpd), ne_m(0:jpd),                   
     .        ni_m(0:jpd), ns_m(0:jpd),                                
     .        angrotp_exp(0:jpd), egamma_exp(0:jpd),                   
     .        gamma_p_exp(0:jpd), vphi_m(0:jpd),                       
     .        vpar_m(0:jpd), vper_m(0:jpd), zeff_exp(0:jpd),           
     .        rho(0:jpd), gradrho_exp(0:jpd), gradrhosq_exp(0:jpd),    
     .        rmin_exp(0:jpd), rmaj_exp(0:jpd),  q_exp(0:jpd),         
     .        shat_exp(0:jpd),alpha_exp(0:jpd), elong_exp(0:jpd)

      namelist /nlglf/ leigen, lprint, nroot, iglf, jshoot, jmm, jmaxm
     & ,  irotstab, te_m, ti_m, ne_m, ni_m, ns_m,itport_pt
     & , igrad, idengrad, zpte_in, zpti_in, zpne_in, zpni_in
     & , angrotp_exp, egamma_exp, gamma_p_exp, vphi_m, vpar_m, vper_m
     & , zeff_exp, bt_exp, bt_flag, rho, arho_exp
     & , gradrho_exp, gradrhosq_exp
     & , rmin_exp, rmaj_exp, rmajor_exp, q_exp, shat_exp
     & , alpha_exp, elong_exp, zimp_exp, amassimp_exp, amassgas_exp
     & , alpha_e, x_alpha
       

!       Namelist dump requires that size of arrays(jpd) be
!       explicitely defined (jpd passed through argument list doesnt work??)
         if(jpd .ne. jpdm1)then
            write(6,FMT ='("jpd must equal nj-1 for namelist write")')
            write(6,FMT ='("jpd =",i4," nj =",i4)')jpd,jpdm1+1
            return
         endif
         if(litport .ne. litport1)
     .  Call STOP('sub writeglf2d,litport .ne. litport1',1)
!      create a new input file each time this routine is called:
       iounit = 110
       call getioun(iounit,iounit)
       OPEN(unit = iounit, file = glf_namelst, status = 'UNKNOWN', 
     .                                           iostat = iostat)


!       write(iounit,nlglf) !this is just too easy, so it doesnt work right
!       instead we have to do it manually:
       write(iounit,FMT ='(4x,
     .          "Namelist input for Glf23 written by Onetwo")')
       write(iounit,FMT = '("$nlglf")')
       write(iounit,2)leigen, lprint, nroot
 2     format(4x,'leigen = ',i3,2x,'lprint = ',i3,2x,'nroot = ',i3)
       write(iounit,3)iglf, jshoot, jmm, jmaxm
 3     format(4x,'iglf =',i3,2x,'jshoot =',i3,2x,'jmm =',i3,2x,
     .         'jmaxm =',i3)
       write(iounit,4)(itport_pt(i),i=1,litport)
 4     format(4x,'itport_pt =',5(i3,2x))
       write(iounit,5)igrad,idengrad,bt_flag,irotstab
 5     format(4x,'igrad = ',i3,2x,'idengrad =',i3,2x,'bt_flag =',i3,/,
     .        4x, "irotstab = ",i3)
       write(iounit,6)zpte_in, zpti_in, zpne_in, zpni_in
 6     format(4x,'zpte_in = ',1pe14.6,4x,'zpti_in = ',1pe14.6,4x,
     .        'zpne_in = ',1pe14.6, 2x,'zpni_in = ',1pe14.6)
       write(iounit,FMT = '(4x,"bt_exp =",1pe14.6,4x,
     . "arho_exp =",1pe14.6,4x,"rmajor_exp =",1pe14.6)')
     .                            bt_exp,arho_exp,rmajor_exp
       write(iounit,FMT = '(4x,"zimp_exp =",1pe14.6,4x,
     . "amassimp_exp =",1pe14.6,4x,"amassgas_exp =",1pe14.6)')
     .  zimp_exp,amassimp_exp ,amassgas_exp
       write(iounit,FMT='(4x,"alpha_e = ",1pe14.6,4x,
     .          "x_alpha = ",1pe14.6)')alpha_e, x_alpha

       write(iounit,FMT = '(4x,"rho(0) =")')
       write(iounit,400)(rho(i),i=0,jpd)
       write(iounit,FMT = '(4x,"te_m(0) =")')
       write(iounit,400)(te_m(i),i=0,jpd)
       write(iounit,FMT = '(4x,"ti_m(0) =")')
       write(iounit,400)(ti_m(i),i=0,jpd)
       write(iounit,FMT = '(4x,"ne_m(0) =")')
       write(iounit,400)(ne_m(i),i=0,jpd)
       write(iounit,FMT = '(4x,"ni_m(0) =")')
       write(iounit,400)(ni_m(i),i=0,jpd)
       write(iounit,FMT = '(4x,"ns_m(0) =")')
       write(iounit,400)(ns_m(i),i=0,jpd)
       write(iounit,FMT = '(4x,"angrotp_exp(0) =")')
       write(iounit,400)(angrotp_exp(i),i=0,jpd)
       write(iounit,FMT = '(4x,"egamma_exp(0) =")')
       write(iounit,400)(egamma_exp(i),i=0,jpd)
       write(iounit,FMT = '(4x,"gamma_p_exp(0) =")')
       write(iounit,400)(gamma_p_exp(i),i=0,jpd)
       write(iounit,FMT = '(4x,"vphi_m(0) =")')
       write(iounit,400)(vphi_m(i),i=0,jpd)
       write(iounit,FMT = '(4x,"vpar_m(0) =")')
       write(iounit,400)(vpar_m(i),i=0,jpd)
       write(iounit,FMT = '(4x,"vper_m(0) =")')
       write(iounit,400)(vper_m(i),i=0,jpd)
       write(iounit,FMT = '(4x,"zeff_exp(0) =")')
       write(iounit,400)(zeff_exp(i),i=0,jpd)
       write(iounit,FMT = '(4x,"gradrho_exp(0) =")')
       write(iounit,400)(gradrho_exp(i),i=0,jpd)
       write(iounit,FMT = '(4x,"gradrhosq_exp(0) =")')
       write(iounit,400)(gradrhosq_exp(i),i=0,jpd)
       write(iounit,FMT = '(4x,"rmin_exp(0) =")')
       write(iounit,400)(rmin_exp(i),i=0,jpd)
       write(iounit,FMT = '(4x,"rmaj_exp(0) =")')
       write(iounit,400)(rmaj_exp(i),i=0,jpd)
       write(iounit,FMT = '(4x,"q_exp(0) =")')
       write(iounit,400)(q_exp(i),i=0,jpd)
       write(iounit,FMT = '(4x,"shat_exp(0) =")')
       write(iounit,400)(shat_exp(i),i=0,jpd)
       write(iounit,FMT = '(4x,"alpha_exp(0) =")')
       write(iounit,400)(alpha_exp(i),i=0,jpd)
       write(iounit,FMT = '(4x,"elong_exp(0) =")')
       write(iounit,400)(elong_exp(i),i=0,jpd)
 400   format((4x,5(2x,1pe14.6)))
       write(iounit,FMT = '(1x,"$END")')
       close(iounit)

       return
       end

