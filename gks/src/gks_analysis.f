c**************************************************************************
       subroutine g_vs_logk(imax,dx,ilin)
c**************************************************************************
c    This function makes a plot of the  growthrate/aky1
c    verses log(aky1).Setting ilin=1 forces a linear spectrum.
c    The quasilinear fluxes are saved and plotted also.
c    GKS code is used for gamma.
c    set by user: shat,pk,eps,epsa,uprim1,fprim1,temp3,teti
c                 zeff,beta,eps1,fprim3, tprim1,tprim3,other gks inputs
c**************************************************************************
       use gks_var
       implicit none
       include 'input.m'
       include 'glf.m'
       include 'gks_out.m'
       include 'data_exp.m'
       real*8 xk
       real*8 dx,save_aky1
       integer i,j,imax,ilin
c
c       if(imax .gt. 50)imax = 50
       if (uprim1 .ne. 0.0) ipar=0
       save_aky1 = aky1
       aky1 = abs(aky1)
       icontinue=0
c  calculate spectum
c
       write(*,*)" g_vs_logk "
       write(*,100)
       j=0
       if(dx.lt.0.0)j=imax
       do i=0,imax
          xk = ABS(REAL(j-i))*dx
          aky1=10.D0**xk
          if(ilin.eq.1) aky1=0.7D0+xk
          call gstotal
          write(*,200) kys(1),aky1,agammas(1),dgammas(1)
     >    ,afreqs(1),dfreqs(1)
          if(dgammas(1) .gt. tol_f*agammas(1))then
             agammas(1)= 0.0
             icontinue = 0
           else 
cgms             icontinue = 1
            icontinue=0
          endif
          aky_m(i)=aky1
          anrate_m(i) = agammas(1)
          dnrate_m(i) = dgammas(1)
          anfreq_m(i) = afreqs(1)
          ky_m(i) = kys(1)
          peflx_m(i) = peflxa/kys(1)
          qeflx_m(i) = eeflxa/kys(1)
          qiflx_m(i) = eiflxa/kys(1)
          if(igks_model.eq.1)then  ! tglf only output
            gamma_ion(i)= agammas(2)
            freq_ion(i) = afreqs(2)
            gamma_electron(i) = agammas(3)
            freq_electron(i) = afreqs(3)
            phi_bar_m(i) = phi_bar_k
            ne_bar_m(i) = ne_bar_k
            te_bar_m(i) = te_bar_k
            ti_bar_m(i) = ti_bar_k
            ne_te_phase_ion(i) = ne_te_phase_k(2)
            ne_te_phase_electron(i) = ne_te_phase_k(3)
          endif        
       enddo
       aky1 = save_aky1
c
 100   format(/,6x,'ky_m',10x,'aky1',10x,'gamma',
     &        9x,'dgamma',8x,'freq',10x,'dfreq')
 200   format(1p6e14.4)
       return
       end
c@gks_max.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c    This subroutine adjusts aky1 to find the maximum growth rate
c    from the GKS code using a golden mean section search.
c       input: all gks inputs set externally, the first call to
c              gstotal must have been done and aky1 set
c              tol_f*gammas(1) is the maximum dgammas(1) accepted
c              as converged
c       
c       if chimax = 1 the maximum of gamma/kys**2 is found
c
c       output: agammas(1),kys(1) set to
c               maximum growth rate case
c               afreqs(1),dgammas(1),dfreqs(1) will be near the
c               maximum case, dgammas(1) is negative if no maximum is
c               found for xkymin_gf < aky1 < xkymax_gf
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine gks_max
c
      use gks_var
      implicit none
      include 'glf.m'
      include 'data_exp.m'
c
c      integer chimax
      real*8 kmax, kmin, xl, xr, xm, xn, xp, dx
      real*8 gl, gr, gm, gn, gp, ww, xptest
      real*8 g_first, dg_first, f_first, kys_first, gmgn
c
c set some constants
c
      xptest = 0.2D0
c      chimax = 0
c      kmax = xkymax_gf
c      kmin = xkymin_gf
      kmax = aky1*4.D0
      kmin = aky1/10.D0
      ww = (3.D0 - DSQRT(5.D0))/2.D0
c
c first search for bracketed maximum
c
            xm = aky1
            gm = agammas(1)
            g_first  = agammas(1)
            dg_first = dgammas(1)
            f_first  = afreqs(1)
            kys_first = kys(1)
            gl = -1.D0
            gn = -1.D0
            gr = -1.D0
            if (dgammas(1) .gt. tol_f*gm) gm = 0.D0
c            if(chimax .eq. 1) then
c              gm = gm/(kys(1)*kys(1))
c              g_first = g_first/(kys(1)*kys(1))
c            endif
c            dx = ww*(1.D0-ww)*xm
            dx = (1.D0-ww)*xm
            xl = xm - ww*dx
            aky1 = xl
            call gstotal
            gl = agammas(1)
            if (dgammas(1) .gt. tol_f*gl) gl = 0.D0
c            if(chimax .eq. 1) gl = gl/(kys(1)*kys(1))
            write(*,*) ' dx = ',dx
            write(*,*) ' xm = ',xm
            write(*,*) ' xl = ',xl
            write(*,*) ' gm = ',gm
            write(*,*) ' gl = ',gl
            if (gm .ge. gl ) then   !  search right
              write(*,*) ' right'
              do
               if (gr .gt. gn) then !  shift brackets right
                 xl = xn
                 gl = gn
                 xm = xr
                 gm = gr
                 gn = -1.D0
                 gr = -1.D0
               endif
               xn = xl + (1.D0-ww)*dx
               aky1 = xn
               call gstotal
               gn = agammas(1)
               if (dgammas(1).gt.tol_f*gn) gn=0.D0
c               if(chimax .eq. 1) gn = gn/(kys(1)*kys(1))
               if( gn .gt. gm) then  !  keep going right
                xr = xm + (1.D0-ww)*dx
                aky1 = xr
                call gstotal
                gr = agammas(1)
                if (dgammas(1).gt.tol_f*gr) gr=0.D0
c                if(chimax .eq. 1) gr = gr/(kys(1)*kys(1))
               endif
               write(*,*) ' xl = ',xl
               write(*,*) ' xm = ',xm
               write(*,*) ' xn = ',xn
               write(*,*) ' xr = ',xr
               write(*,*) ' gl = ',gl
               write(*,*) ' gm = ',gm
               write(*,*) ' gn = ',gn
               write(*,*) ' gr = ',gr
               if (xl.lt.kmin.or.xr.gt.kmax.or.
     &              gl+gm+gn.eq.0.D0.or.gn.gt.gr) exit
c             until(xl.lt.kmin.or.xr.gt.kmax.or.
c    &              gl+gm+gn.eq.0.D0.or.gn.gt.gr)
              enddo
            else
c search left
              write(*,*) ' left'
              do
c shift brackets left
               xr = xm
               gr = gm
               xn = xl
               gn = gl
               gm = -1.D0
               gl = -1.D0
               xm = xr - (1.D0-ww)*dx
               aky1=xm
               call gstotal
               gm = agammas(1)
               if (dgammas(1).gt.tol_f*gm) gm=0.D0
c               if(chimax .eq. 1) gm = gm/(kys(1)*kys(1))
               if (gm .gt. gn) then    !   keep going left
                 xl = xn - (1.D0-ww)*dx
                 aky1 = xl
                 call gstotal
                 gl = agammas(1)
                 if (dgammas(1).gt.tol_f*gl) gl=0.D0
c                 if(chimax .eq. 1) gl = gl/(kys(1)*kys(1))
               endif
               write(*,*) ' xl = ',xl
               write(*,*) ' xm = ',xm
               write(*,*) ' xn = ',xn
               write(*,*) ' xr = ',xr
               write(*,*) ' gl = ',gl
               write(*,*) ' gm = ',gm
               write(*,*) ' gn = ',gn
               write(*,*) ' gr = ',gr
               if (xl.lt.kmin.or.xr.gt.kmax.or.
     &               gm+gn+gr.eq.0.D0.or.gm.gt.gl) exit
              enddo
c             until (xl.lt.kmin.or.xr.gt.kmax.or.
c    &               gm+gn+gr.eq.0.D0.or.gm.gt.gl)
            endif
            if (xl.ge.kmin .and. xr.le.kmax .and.
     &          gm+gn .ne.0.D0) then
c  found bracketed maximum, shrink brackets
             if (gm .ge. gn) then
                xp = (xn-xl)/xm
                if(xp .gt. xptest)then
                 xr = xn
                 gr = gn
                 xn = xm
                 gn = gm
                 xm = xl + ww*(xr-xl)
                 gm = -1.D0
                endif
             else
                xp = (xr-xm)/xn
                if(xp .gt. xptest)then
                 xl = xm
                 gl = gm
                 xm = xn
                 gm = gn
                 xn = xl + (1.D0-ww)*(xr-xl)
                 gn = -1.D0
                endif
             endif
             write(*,*) ' xp = ',xp
             write(*,*) ' xptest = ',xptest
             if(xp .gt. xptest) then
               write(*,*) ' shrink'
              do
               if (gn .eq. -1.D0)then
                aky1=xn
               else
                aky1=xm
               endif
               call gstotal
               gp = agammas(1)
               if (dgammas(1) .gt. tol_f*gp) gp=0.D0
c               if (chimax .eq. 1) gp = gp/(kys(1)*kys(1))
               if (gn .eq. -1.D0) then
                gn = gp
               else
                gm = gp
               endif
               xp = (xr-xl)/xm
               if (xp .gt. xptest) then
                if (gn .ge. gm) then    !  shrink bracket right
                 write(*,*) ' s right'
                 xl = xm
                 gl = gm
                 xm = xn
                 gm = gn
                 xn = xl + (1.D0-ww)*(xr-xl)
                 gn = -1.0
                else                    !  shrink bracket left
                 write(*,*) ' s Left'
                 xr = xn
                 gr = gn
                 xn = xm
                 gn = gm
                 xm = xl + ww*(xr-xl)
                endif
               endif
               write(*,*) ' xl = ',xl
               write(*,*) ' xm = ',xm
               write(*,*) ' xn = ',xn
               write(*,*) ' xr = ',xr
               write(*,*) ' xp = ',xp
               write(*,*) ' gl = ',gl
               write(*,*) ' gm = ',gm
               write(*,*) ' gn = ',gn
               write(*,*) ' gr = ',gr
               if (xp .le. xptest) exit
             enddo
            endif
            if (gn .ge. gm) then
              agammas(1) = gn
              kys(1) = kys(1)*xn/aky1
            else
              agammas(1) = gm
              kys(1) = kys(1)*xm/aky1
            endif
          else
             if(xl.lt.kmin .or. xr.gt.kmax) then
c exceeded bounds tag kys negative
               kys(1)=-kys(1)             
c failed to find a maximum tag dgammas(1) negative and leave ky unchanged
             else
               agammas(1)=g_first
               afreqs(1) = f_first
               kys(1) = kys_first
               dgammas(1) = -dg_first
             endif
          endif
c
      return
      end
c@gks_profile.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c    This subroutine calls gstotal in GKS to find the growth rates
c    at the gridpoints from jstart_prof to jend_prof.
c    It is assumed that profile data already exists.
c
c    The ky is determined as follows:
c           if kys0 = 0 aky1 is used
c           if kys0 = 0 and aky1 < 0 aky1 is assumed to be wavenumber (1/cm)
c           if kys0 > 0 kys0 is used
c           if kys0 < 0 the maximum growth rate is found
c    Output:
c           anrate_m(j) the growth rate profile
c           dnrate_m(j) the error in anrate_m
c           anfreq_m(j) the frequency profile
c           dnfreq_m(j) the error in the anfreq_m
c           ky_m(j)     the wavenumber profile in rho_s units
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine gks_profile(jstart_prof,jend_prof)
c
      use gks_var
      implicit none
      include 'input.m'
      include 'data_exp.m'
      include 'glf.m'
      include 'gks_out.m'
c
      real*8 zpimp, shat_min, dshat, save_ls
      real*8 gamma1, freq1, dgamma1
      real*8 gamma2, freq2, dgamma2
      real*8 save_aky1, kys_last
      real*8 tem,tim,nem,nim,nzm
      real*8 qm,rminm,rmajm,zeffm
      integer j, dj, k, iprint_gks
      integer jstart_prof,jend_prof
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      iprint_gks = 0 ! no diagnostic printout
c
c        write(*,*) 'inside gks_profile ...'
c        write (*,'(a10,i3)') 'nstep  = ',nstep
c        write (*,'(a10,i2)') 'igeo_m = ',igeo_m
c        write (*,'(a10,i2)') 'igeo   = ',igeo
c        write (*,'(a10,i2)') 'iptot  = ',iptot
c        write (*,'(a13,1p2e14.6)') 'delt   =    ',delt
c        write (*,'(a13,1p2e14.6)') 'kys0   =    ',kys0
c        if(kys0 .lt. 0.0) then
c          write (*,'(a13,1p2e14.6)') 'xkymin_gf = ',xkymin_gf
c          write (*,'(a13,1p2e14.6)') 'xkymax_gf = ',xkymax_gf
c        endif
c        write (*,'(a13,1p2e14.6)') 'tol_f  =    ',tol_f
c        write (*,'(a13,1p2e14.6)') 'xalpha =    ',abs(xalpha)
c        write (*,'(a13,1p2e14.6)') 'cbetae =    ',cbetae
c        write (*,'(a13,1p2e14.6)') 'alpha_p=    ',alpha_p
c        write (*,'(a13,1p2e14.6)') 'alpha_p_cur=    ',alpha_p_cur
c
c... set some constants
c
      shat_min = 0.10
      shat_min = 0.0
      if(ABS(alpha_p)+ABS(alpha_p_cur) .ne. 0.0) ipar = 0
      bt_s = bt_exp
      amass3_s = (5.446D-4)/amassgas_exp
      rmin_s = arho_exp
      kys_last =  abs(kys0)
      save_aky1 = aky1
c
c... begin loop over grid
c
      if (jstart_prof .le. jend_prof) then
         dj = 1
      else
         dj = -1
      endif
c
      do j=jstart_prof,jend_prof,dj
         tem=0.5*(te_exp(j+1)+te_exp(j))
         tim=0.5*(ti_exp(j+1)+ti_exp(j))
         nem=0.5*(ne_exp(j+1)+ne_exp(j))
         nim=0.5*(ni_exp(j+1)+ni_exp(j))
         nzm=0.5*(nz_exp(j+1)+nz_exp(j))
         write(*,*)"nem=",nem
         write(*,*)"nim=",nim
         write(*,*)"nzm=",nzm
         qm=0.5*(q_exp(j+1)+q_exp(j))
         rminm=0.5*(rmin_exp(j+1)+rmin_exp(j))
         rmajm=0.5*(rmaj_exp(j+1)+rmaj_exp(j))
         zeffm=0.5*(zeff_exp(j+1)+zeff_exp(j))
c
         zeff_s = zeffm
         rmaj_s = rmajm
         rsurf_s = rminm
         ne_s = nem
         ni_s = nim
         nz_s = nzm
         nfast_s = nem - nim - z2*nzm
         if (nfast_s .le. 0.0) nfast_s = 0.0
         te_s = tem
         ti_s = tim
         qgks_s = qm
         lne_s = arho_exp/zpne_exp(j)
         lni_s = arho_exp/zpni_exp(j)
         lnz_s = arho_exp/zpnz_exp(j)
         lti_s = arho_exp/zpti_exp(j)
         lte_s = arho_exp/zpte_exp(j)
         ls_s = rmajm*qm/shat_exp(j)
         if(j.eq.35)write(*,*)"drhodr=",drhodr(j)
         if (igeo_m .ge. 2) then
            lne_s = lne_s/drhodr(j)
            lni_s = lni_s/drhodr(j)
            lnz_s = lnz_s/drhodr(j)
            lte_s = lte_s/drhodr(j)
            lti_s = lti_s/drhodr(j)
         else
            lne_s = lne_s/sqrt(elong_exp(j))
            lni_s = lni_s/sqrt(elong_exp(j))
            lnz_s = lnz_s/sqrt(elong_exp(j))
            lte_s = lte_s/sqrt(elong_exp(j))
            lti_s = lti_s/sqrt(elong_exp(j))
         endif
         if (igeo_m .eq. 2) ls_s = ls_s/drhodrrrho(j)
         shift_s = abs(xalpha)*alpha_exp(j)
         uprim1_s = sqrt(2.0*te_s/ti_s)*gamma_p_i_exp(j)
         uprim3_s = sqrt(2.0*te_s/ti_s)*gamma_p_e_exp(j)
         egamma = alpha_e*egamma_exp(j)
         mach1 = mach_i_exp(j)
         mach3 = mach_e_exp(j)
         delt_s = delt
c... choose aky1
         aky1_s = save_aky1
         if (save_aky1 .lt. 0) then
           aky1_s = -aky1_s*rhosda_exp(j)*arho_exp*100.D0
c            aky1_s = (-1.0e5)*sqrt(2.0*ti_s/te_s)*aky1_s*te_s/
c     &               (bt_s*csda_exp(j)*arho_exp)
         endif
         if (kys0 .ne. 0) aky1_s = kys_last*sqrt(2.0*ti_s/te_s)
c
c... for weak mag. shear, interpolate
c
         if(abs(shat_exp(j)).lt.shat_min .and. kys0.ge.0.0) then
            save_ls = ls_s
            ls_s = -shat_exp(j)*save_ls/shat_min
            call sharegks
            if (igeo .eq. 0) call ifwrite
            if (igeo .eq. 1) call ifwritegeo
            fprim1 = rmin_s/lni_s
            fprim3 = (rmin_s/lne_s)/fprim1
c            fprim2 = fprim3
            fprim2=zpnz_exp(j)/zpni_exp(j)
            fprim4 = 1.0e-6
            beta = cbetae*beta
            call gstotal
            gamma1 = agammas(1)
            freq1  = afreqs(1)
            dgamma1 = dgammas(1)
            ls_s = shat_exp(j)*save_ls/shat_min
            call sharegks
            if (igeo .eq. 0) call ifwrite
            if (igeo .eq. 1) call ifwritegeo
            fprim1 = rmin_s/lni_s
            fprim3 = (rmin_s/lne_s)/fprim1
c            fprim2 = fprim3
            fprim2 = zpnz_exp(j)/zpni_exp(j)
            fprim4 = 1.0e-6
            beta = cbetae*beta
c            shat=shat*drhodrrrho(j)
            call gstotal
            gamma2 = agammas(1)
            freq2  = afreqs(1)
            dgamma2 = dgammas(1)
            dshat = (shat_exp(j) + shat_min)/(2.0*shat_min)
            agammas(1) = gamma1 + dshat*(gamma2-gamma1)
            afreqs(1) = freq1 + dshat*(freq2-freq1)
            dgammas(1) = dgamma1 + dshat*(dgamma2-dgamma1)
            shat = rmin_s*qgks_s/(save_ls*epsa)
         else
           call sharegks
            if (igeo .eq. 0) call ifwrite
            if (igeo .eq. 1) call ifwritegeo
            fprim1 = rmin_s/lni_s
            fprim3 = (rmin_s/lne_s)/fprim1
c            fprim2 = fprim3
            fprim2 = zpnz_exp(j)/zpni_exp(j)
            fprim4 = 1.0e-6
            beta = cbetae*beta
c... run gks code
            call gstotal
            write(*,*)"an1=",an1,nim/nem
            write(*,*)"an2=",an2,nzm/nem
         endif
c
c... find maximum growth rate by golden section search
         if (kys0 .lt. 0.0) then
c           write(*,*) 'calling gks_max ...'
c           write(*,*)
           call gks_max
         endif
c
c... store the results for this gridpoint
c
         anrate_m(j) = agammas(1)
         dnrate_m(j) = dgammas(1)
         anfreq_m(j) = afreqs(1)
         dnfreq_m(j) = dfreqs(1)
         ky_m(j)     = kys(1)
         kys_last    = kys(1)
         if(igks_model.eq.1)then  ! tglf only output
           gamma_ion(j)= agammas(2)
           freq_ion(j) = afreqs(2)
           gamma_electron(j) = agammas(3)
           freq_electron(j) = afreqs(3)
           phi_bar_m(j) = phi_bar_k
           ne_bar_m(j) = ne_bar_k
           te_bar_m(j) = te_bar_k
           ti_bar_m(j) = ti_bar_k
           ne_te_phase_ion(j) = ne_te_phase_k(2)
           ne_te_phase_electron(j) = ne_te_phase_k(3)
         endif        
         write(*,*)j
         write(*,*) 'agammas(1) = ',agammas(1)
         write(*,*) 'dgammas(1) = ',dgammas(1)
         write(*,*) 'afreqs(1) = ',afreqs(1)
         write(*,*) 'kys(1) = ',kys(1)
c         write(*,*) 'peflxa = ',peflxa
c         write(*,*) 'eeflxa = ',eeflxa
c         write(*,*) 'eiflxa = ',eiflxa
c
c... diagnostic printout
c
         if (iprint_gks.gt.0) then
           write(*,*) 'ntheta = ',ntheta
           write(*,*) 'nperiod = ',nperiod
           write(*,*) 'ipar = ',ipar
           write(*,*) 'ngauss = ',ngauss
           write(*,*) 'icv = ',icv
           write(*,*) 'negrid = ',negrid
           write(*,*) 'nspec = ',nspec
           write(*,*) 'eps = ',eps
           write(*,*) 'shift = ',shift
           write(*,*) 'dbeam = ',dbeam
           write(*,*) 'shat = ',shat
           write(*,*) 'pk = ',pk
           write(*,*) 'epsa = ',epsa
           write(*,*) 'epsl = ',epsl
           write(*,*) 'width0 = ',width0
           write(*,*) 'beta = ',beta
           write(*,*) 'zeff = ',zeff
           write(*,*) 'teti = ',teti
           write(*,*) 'zeff5 = ',zeff5
           write(*,*) 'fprim1 = ',fprim1
           write(*,*) 'fprim2 = ',fprim2
           write(*,*) 'fprim3 = ',fprim3
           write(*,*) 'fprim4 = ',fprim4
           write(*,*) 'fprim5 = ',fprim5
           write(*,*) 'tprim1 = ',tprim1
           write(*,*) 'tprim2 = ',tprim2
           write(*,*) 'tprim3 = ',tprim3
           write(*,*) 'tprim4 = ',tprim4
           write(*,*) 'tprim5 = ',tprim5
           write(*,*) 'vnewk1 = ',vnewk1
           write(*,*) 'vnewk2 = ',vnewk2
           write(*,*) 'vnewk3 = ',vnewk3
           write(*,*) 'vnewk4 = ',vnewk4
           write(*,*) 'vnewk5 = ',vnewk5
           write(*,*) 'bakdif1 = ',bakdif1
           write(*,*) 'bakdif2 = ',bakdif2
           write(*,*) 'bakdif3 = ',bakdif3
           write(*,*) 'bakdif4 = ',bakdif4
           write(*,*) 'bakdif5 = ',bakdif5
           write(*,*) 'amass2 = ',amass2
           write(*,*) 'amass3 = ',amass3
           write(*,*) 'amass4 = ',amass4
           write(*,*) 'amass5 = ',amass5
           write(*,*) 'temp2 = ',temp2
           write(*,*) 'temp3 = ',temp3
           write(*,*) 'temp4 = ',temp4
           write(*,*) 'temp5 = ',temp5
           write(*,*) 'z2 = ',z2
           write(*,*) 'z3 = ',z3
           write(*,*) 'z4 = ',z4
           write(*,*) 'z5 = ',z5
         endif
c
c        j = j + dj
c        if (j .eq. jend_prof + dj) goto 900
      enddo
c
c        write(*,*)
c        do k=jstart_prof,jend_prof
c          write(*,100) k,rho(k),ky_j(k),anrate_m(k),dnrate_m(k),
c     &                 anfreq_m(k),dnfreq_m(k)
c        enddo
c
  50  format(1p2e14.6)
 100  format(2x,i2,2x,0p2f9.5,1p6e14.6)
 900  continue
      return
      end
c@geogks_profile.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c  routine goes over an input experimental profile to
c  calculate the dimensionless variables which drive ifwritegeo
c  for Miller's triangulated shifted ellipse local equillibrium
c  input to gyrokinetic stability code gks
c
c  ...the usual ifwrite variables of specific density, temperature,
c  and b_field are avoided in vafor of the dimensionless variables
c                                                            
c  jgeoin (inside); jgeoout (outside); jgeodel (spaceing) grid
c
c
c   a_mult: aspectratio_loc_s(j)*a_mult
c                 shift_loc_s(j)*a_mult, 
c                rmaj_mag_center*a_mult
c   ab_div:     beta_loc_0_s(j)/ab_div
c                shift_loc_s(j)/ab_div
c         ....ab_div=a_mult keeps alpha_mhd_loc_s(j) fixed
c                           and shift_loc_s(j) fixed
c   b_div:         bt_mag_center/b_div
c         .... but beta_loc_0 same ie n(j)*t(j)/b_div**2
c         .... and rhosda_loc_s(j)*b_div
c              ie only confinement but not stability is changed
c         .....b_div=a_mult  for constant center post current
c   e_mult:        (kappa_loc_s(j)-1.)*e_mult+1.
c         arho_exp*(((kappa_loc_s(j)-1.)*e_mult+1.)/kappa_loc_s(j))**.5
c   d_mult:        delta_loc_s(j)*d_mult
c
c   this function calls gstotal to find growth rates
c   at the gridpoints from jstart_prof to jend_prof
c   it is assumed that runstack_exp has been run
c   The ky is determined as follows:
c           if kys0 = 0 aky1 is used
c           if kys0 = 0 and aky1 < 0 aky1 is assumed to be wavenumber in 1/cm
c           if kys0 > 0 kys0 is used
c           if kys0 < 0 the maximum growth rate is found
c   Output:
c           anrate_m(j) the growth rate profile
c           dnrate_m(j) the error in anrate_m
c           anfreq_m(j) the frequency profile
c           ky_m(j)     the wavenumber profile in rho_s units
c   for TGLF only
c           phi_bar_m(j) the incremental potential fluctuation intensity at ky_m
c           ne_bar_m(j) the incremental electron density fluctuation intensity at ky_m
c           te_bar_m(j) the incremental electron temperature fluctuation intensity at ky_m
c           ti_bar_m(j) the incremental ion temperature fluctuation intensity at ky_m
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine geogks_profile(jgeoin,jgeoout,jgeodel)
c
      use gks_var
      implicit none
      include 'input.m'
      include 'glf.m'
      include 'data_exp.m'
      include 'gks_out.m'
c
      real*8 rmaj_exp_mult(0:jmaxm)
      real*8 a_mult,ab_div,b_div,e_mult,d_mult
      real*8 kys_last, save_aky1
      integer jcount,jgeoin,jgeoout,jgeodel,jgeostep
      integer j, jp, k
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c      write(*,*) 'inside geogks_profile ...'
c      if (igeo.ne.1) then
c        write(*,*) 'Warning: igeo not set to 1, exiting'
c        return
c      endif
c      write (*,'(a10,i3)') 'nstep  = ',nstep
c      write (*,'(a10,i2)') 'igeo_m = ',igeo_m
c      write (*,'(a10,i2)') 'igeo   = ',igeo
c      write (*,'(a10,i2)') 'iptot  = ',iptot
c      write (*,'(a13,1p2e14.6)') 'delt   =    ',delt
c      write (*,'(a13,1p2e14.6)') 'kys0   =    ',kys0
c      if(kys0 .lt. 0.0) then
c        write (*,'(a13,1p2e14.6)') 'xkymin_gf = ',xkymin_gf
c        write (*,'(a13,1p2e14.6)') 'xkymax_gf = ',xkymax_gf
c      endif
c      write (*,'(a13,1p2e14.6)') 'tol_f  =    ',tol_f
c      write (*,'(a13,1p2e14.6)') 'xalpha =    ',abs(xalpha)
c      write (*,'(a13,1p2e14.6)') 'cbetae =    ',cbetae
c      write (*,'(a13,1p2e14.6)') 'alpha_p_i=    ',alpha_p_i
c      write (*,'(a13,1p2e14.6)') 'alpha_p_e=    ',alpha_p_e
c
c set scalings
c
      a_mult=1.D0
      ab_div=1.D0
      b_div=1.D0
      e_mult=1.D0
cgms 3/13/06      d_mult=-1.D0
      d_mult=1.D0
c
c      igeo_print=0    ! turn off diagnostic printout
      p_prime_zero=-1.D0  
      jp=1
      kys_last=abs(kys0)
      save_aky1=aky1
c
      if(i_bpar.eq.0) then
        p_prime_zero=1.D-9
        x_bpar=1.D-9
        y_bpar=1.D-9
        xy_bpar=1.D-9
        y_mhd=0.D0
      endif
c
      if(i_bpar.eq.1) then
        p_prime_zero=-1.D0
        x_bpar=1.D0
        y_bpar=1.D0
        xy_bpar=1.D0
        y_mhd=0.D0
      endif
c
      rmaj_exp_mult(jmaxm)=rmaj_exp(jmaxm)*a_mult
      do j=jmaxm-1,0,-1
        rmaj_exp_mult(j)=rmaj_exp_mult(j+1)+
     &    (rmaj_exp(j)-rmaj_exp(j+1))*a_mult/ab_div
      enddo
c 
      write(*,*)"setunits"
      call geogks_setunits(a_mult,ab_div,b_div,e_mult)
c
      jgeostep = jgeodel
      if(jgeoin.gt.jgeoout)jgeostep=-abs(jgeodel)
      jcount=0
c
      do j=jgeoin,jgeoout,jgeostep
       if(kys_last.ge.0.D0)then
        jcount=jcount+1
c        aspectratio_loc=rmaj_exp_mult(j)/rmin_exp(j)
c        shift_loc=(rmaj_exp_mult(j+1)-rmaj_exp_mult(j-1))/
c     &            (rmin_exp(j+1)-rmin_exp(j-1))
c
      write(*,*)"setup"
        call geogks_setup(j,ab_div,b_div,e_mult,d_mult)
c
c convert_theta_theta_c
c
        if(save_aky1.lt.0.0)then
          aky1 = abs(save_aky1)*100.0*rhosda_loc_s(j)*a_unit_exp
          aky1 = aky1*sqrt(2.0*tiote_loc)
        endif
        if(kys0.ne.0.0) aky1 = kys_last*sqrt(2.0*tiote_loc)
c        aky1 = kys_last*sqrt(2.)*sqrt(tiote_loc)
      write(*,*)"gstotal",kys0
        call gstotal
        if(kys0.lt.0.0) call gks_max
c
      write(*,*)"return from gstotal"
        anrate_m(j) = agammas(1)
        dnrate_m(j) = dgammas(1)
        anfreq_m(j) = afreqs(1)
        dnfreq_m(j) = dfreqs(1)
        ky_m(j) = kys(1)
        if(igks_model.eq.1)then  ! tglf only output
          gamma_ion(j)= agammas(2)
          freq_ion(j) = afreqs(2)
          gamma_electron(j) = agammas(3)
          freq_electron(j) = afreqs(3)
          phi_bar_m(j) = phi_bar_k
          ne_bar_m(j) = ne_bar_k
          te_bar_m(j) = te_bar_k
          ti_bar_m(j) = ti_bar_k
          ne_te_phase_ion(j) = ne_te_phase_k(2)
          ne_te_phase_electron(j) = ne_te_phase_k(3)
          mhd_DR_m(j) = mhd_DR_k
        endif
c        write(*,*)
        write(*,*) 'agammas(1) = ',agammas(1)
        write(*,*) 'dgammas(1) = ',dgammas(1)
        write(*,*) 'afreqs(1) = ',afreqs(1)
        write(*,*) 'kys(1) = ',kys(1)
c        write(*,*) 'peflxa = ',peflxa
c        write(*,*) 'eeflxa = ',eeflxa
c        write(*,*) 'eiflxa = ',eiflxa
c
        if(dgammas(1).gt.0.D0) kys_last=kys(1)
       endif
      enddo
c
c        write(*,*)
c        do k=jgeoin,jgeoout,jgeostep
c          write(*,100) k,rho(k),ky_j(k),anrate_m(k),dnrate_m(k),
c     &                 anfreq_m(k),dnfreq_m(k)
c        enddo
c
 100  format(2x,i2,2x,0p2f9.5,1p6e14.6)
c
      return
      end
c@geogks_setunits.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c   a_mult: aspectratio_loc_s(j)*a_mult
c                 shift_loc_s(j)*a_mult, 
c                rmaj_mag_center*a_mult
c   ab_div:     beta_loc_0_s(j)/ab_div
c                shift_loc_s(j)/ab_div
c         ....ab_div=a_mult keeps alpha_mhd_loc_s(j) fixed
c                           and shift_loc_s(j) fixed
c   b_div:         bt_mag_center/b_div
c         .... but beta_loc_0 same ie n(j)*t(j)/b_div**2
c         .... and rhosda_loc_s(j)*b_div
c              ie only confinement but not stability is changed
c         .....b_div=a_mult  for constant center post current
c   e_mult:        (kappa_loc_s(j)-1.)*e_mult+1.
c         arho_exp*(((kappa_loc_s(j)-1.)*e_mult+1.)/kappa_loc_s(j))**.5
c
c  global variables set in this function:
c        real a_unit_exp
c        real bt_mag_center
c        real rmaj_mag_center 
c        real b_unit_loc_s(0:jmaxm)
c       real rhosda_loc_s(0:jmaxm),csda_loc_s(0:jmaxm)
c
c  notes :
c      bt_mag_center=bt_exp is field on magnetic axis without plasma
c      rmaj_mag_center=magnetic axis without plasma=rmaj_exp(jmaxm)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine geogks_setunits(a_mult,ab_div,b_div,e_mult)
c
      use gks_var
      implicit none
      include 'input.m'
      include 'glf.m'
      include 'data_exp.m'
c
      real*8 a_mult,ab_div,b_div,e_mult
      real*8 tem,rhom,rminm
      integer j
c
      a_unit_exp=rmin_exp(jmaxm)
      bt_mag_center=bt_exp/b_div
      rmaj_mag_center=rmaj_exp(jmaxm)*a_mult/a_unit_exp
c
      do j=1,jmaxm-1
cgms        b_unit_loc_s(j)=bt_mag_center*
cgms     &    ((elong_exp(j)-1.)*e_mult+1.)/elong_exp(j)*
cgms     &    rho(j)*arho_exp/rmin_exp(j)*(rho(j+1)-rho(j-1))*arho_exp/
cgms     &    (rmin_exp(j+1)-rmin_exp(j-1))
         rhom=0.5*(rho(j+1)+rho(j))*arho_exp
         rminm=0.5*(rmin_exp(j+1)+rmin_exp(j))
         b_unit_loc_s(j) = bt_mag_center*(rhom/rminm)*
     >   arho_exp*(rho(j+1)-rho(j))/(rmin_exp(j+1)-rmin_exp(j))
      enddo
c
      b_unit_loc_s(0)=b_unit_loc_s(1) 
      b_unit_loc_s(jmaxm)=b_unit_loc_s(jmaxm-1)
cgms      b_unit_loc_s(jmaxm)=bt_mag_center*
cgms     &  ((elong_exp(jmaxm)-1.0)*e_mult+1.0)/elong_exp(jmaxm)*
cgms     &  rho(jmaxm)*arho_exp/rmin_exp(jmaxm)*
cgms     &  (rho(jmaxm)-rho(jmaxm-1))*arho_exp/
cgms     &  (rmin_exp(jmaxm)-rmin_exp(jmaxm-1))
c
      do j=0,jmaxm
        tem=0.5*(te_exp(j+1)+te_exp(j))
        rhosda_loc_s(j)=((1.02D2*(tem*1.D3)**0.5)/
     &                  b_unit_loc_s(j)/1.D4)*
     &                  (amassgas_exp)**0.5/(a_unit_exp*100.D0)
        csda_loc_s(j)=9.79D5*(tem*1.D3)**0.5/
     &                (a_unit_exp*100.D0)/(amassgas_exp)**0.5
      enddo
c
      return
      end
c@geogks_setup.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c   e_mult:        (kappa_loc_s(j)-1.)*e_mult+1.
c         arho_exp*(((kappa_loc_s(j)-1.)*e_mult+1.)/kappa_loc_s(j))**.5
c   d_mult:        delta_loc_s(j)*d_mult
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine geogks_setup(jloc,ab_div,b_div,e_mult,d_mult)
c
      use gks_var
      implicit none
      include 'input.m'
      include 'glf.m'
      include 'data_exp.m'
c
      real*8 kappa_loc_p,kappa_loc_m
      real*8 alpha_c_s(0:jmaxm)
      real*8 alpha_c_exp(0:jmaxm)
      real*8 b2_ave_geo_s(0:jmaxm)
      real*8 dbeta_loc_0drho_s(0:jmaxm)
      real*8 local_vnewk3 
      real*8 beta_loc_0_sav
      real*8 alpha_mhd_loc_p
      real*8 ab_div,b_div,e_mult,d_mult
      real*8 lnlamda,taue,mui,vthi,xnuei
      real*8 tem,tim,nem,nim,nzm,ptotm,zeffm
      real*8 rminm,rmajm,qm,deltam,elongm
      integer j, jp, jloc
c
       pi=atan2(0.0,-1.0)
c
      j = jloc
      jp = 1
      tem=0.5*(te_exp(j+1)+te_exp(j))
      tim=0.5*(ti_exp(j+1)+ti_exp(j))
      nem=0.5*(ne_exp(j+1)+ne_exp(j))
      nim=0.5*(ni_exp(j+1)+ni_exp(j))
      nzm=0.5*(nz_exp(j+1)+nz_exp(j))
      ptotm=0.5*(ptot_exp(j+1)+ptot_exp(j))
      zeffm=0.5*(zeff_exp(j+1)+zeff_exp(j))
      rminm=0.5*(rmin_exp(j+1)+rmin_exp(j))
      rmajm=0.5*(rmaj_exp(j+1)+rmaj_exp(j))
      qm=0.5*(q_exp(j+1)+q_exp(j))
      deltam=0.5*(delta_exp(j+1)+delta_exp(j))
      elongm=0.5*(elong_exp(j+1)+elong_exp(j))
c
      rmin_loc=rminm/a_unit_exp
      aspectratio_loc=(rmajm/a_unit_exp)/rmin_loc
      shift_loc = (rmaj_exp(j+1)-rmaj_exp(j))
     >            /(rmin_exp(j+1)-rmin_exp(j))
      delta_loc=deltam*d_mult
c      write(*,*)"geogks_setup",j,d_mult,delta_loc
      kappa_loc=(elongm-1.D0)*e_mult+1.D0
      s_delta_loc=rminm/(1.D0-deltam**2*d_mult**2)**0.5*
     &            d_mult*(delta_exp(j+1)-delta_exp(j))/
     &            (rmin_exp(j+1)-rmin_exp(j))
      kappa_loc_p=(elong_exp(j+1)-1.D0)*e_mult+1.D0
      kappa_loc_m=(elong_exp(j)-1.D0)*e_mult+1.D0
      s_kappa_loc=rminm/kappa_loc*
     &           (kappa_loc_p-kappa_loc_m)/
     &           (rmin_exp(j+1)-rmin_exp(j)) 
      q_loc=qm
      shat_loc=(rminm/qm)*(q_exp(j+1)-q_exp(j))/
     &         (rmin_exp(j+1)-rmin_exp(j))
      tiote_loc=tim/tem
c
      nione_loc=nim/nem
      fastionfrac_loc=1.0-nione_loc-z2*nzm/nem
      if(fastionfrac_loc.lt.0.0)fastionfrac_loc=0.0
      zeff_loc=zeffm
c
cgms      dlntidr_loc=-(LOG(ti_exp(j-1))-LOG(ti_exp(j+jp)))/
cgms     &            (rmin_exp(j-1)-rmin_exp(j+jp))*a_unit_exp+1.D-10
cgms      dlntedr_loc=-(LOG(te_exp(j-1))-LOG(te_exp(j+jp)))/
cgms     &            (rmin_exp(j-1)-rmin_exp(j+jp))*a_unit_exp+1.D-10  
cgms      dlnnidr_loc=-(LOG(ni_exp(j-1))-LOG(ni_exp(j+jp)))/
cgms     &            (rmin_exp(j-1)-rmin_exp(j+jp))*a_unit_exp+1.D-10  
cgms      dlnnedr_loc=-(LOG(ne_exp(j-1))-LOG(ne_exp(j+jp)))/
cgms     &            (rmin_exp(j-1)-rmin_exp(j+jp))*a_unit_exp+1.D-10 
cgms      dlnnimpdr_loc=-(LOG(nz_exp(j-1))-LOG(nz_exp(j+jp)))/
cgms     &              (rmin_exp(j-1)-rmin_exp(j+jp))*a_unit_exp+1.D-10 
c 
      dlntidr_loc = -(a_unit_exp/tim)*(ti_exp(j+1)-ti_exp(j))/
     >                (rmin_exp(j+1)-rmin_exp(j)) 
      dlntedr_loc = -(a_unit_exp/tem)*(te_exp(j+1)-te_exp(j))/
     >                (rmin_exp(j+1)-rmin_exp(j)) 
      dlnnidr_loc = -(a_unit_exp/nim)*(ni_exp(j+1)-ni_exp(j))/
     >                (rmin_exp(j+1)-rmin_exp(j)) 
      dlnnedr_loc = -(a_unit_exp/nem)*(ne_exp(j+1)-ne_exp(j))/
     >                (rmin_exp(j+1)-rmin_exp(j)) 
      dlnnimpdr_loc = -(a_unit_exp/nzm)*(nz_exp(j+1)-nz_exp(j))/
     >                (rmin_exp(j+1)-rmin_exp(j)) 
      if(ABS(dlntidr_loc).lt.1.D-10)dlntidr_loc=1.D-10
      if(ABS(dlntedr_loc).lt.1.D-10)dlntedr_loc=1.D-10
      if(ABS(dlnnidr_loc).lt.1.D-10)dlnnidr_loc=1.D-10
      if(ABS(dlnnedr_loc).lt.1.D-10)dlnnedr_loc=1.D-10
      if(ABS(dlnnimpdr_loc).lt.1.D-10)dlnnimpdr_loc=1.D-10
c
c      write(*,*) 'geogks_setup'
      write(*,*) 'j = ',jloc
c      write(*,*) 'aspectratio_loc = ',aspectratio_loc
c      write(*,*) 'rmin_loc = ',rmin_loc
c      write(*,*) 'kappa_loc = ',kappa_loc
c      write(*,*) 's_kappa_loc = ',s_kappa_loc
c      write(*,*) 'delta_loc = ',delta_loc
c      write(*,*) 's_delta_loc = ',s_delta_loc
c      write(*,*) 'shift_loc = ',shift_loc
c      write(*,*) 'q_loc = ',q_loc
c      write(*,*) 'shat_loc = ',shat_loc
c      write(*,*) 'tiote_loc = ',tiote_loc
c      write(*,*) 'nione_loc = ',nione_loc
c      write(*,*) 'dlntidr_loc = ',dlntidr_loc
c      write(*,*) 'dlntedr_loc = ',dlntedr_loc
c      write(*,*) 'dlnnidr_loc = ',dlnnidr_loc
c      write(*,*) 'dlnnedr_loc = ',dlnnedr_loc
c      write(*,*) 'dlnnimpdr_loc =',dlnnimpdr_loc
c      write(*,*) 'nione_loc =',nione_loc
c      write(*,*) 'fastionfrac_loc =',fastionfrac_loc
c      write(*,*) 'nimp/ne =',nzm/nem
c
cgms changed data reader so that ptot_exp is filled even if iptot=0
cgms      dlnpdr_loc=0.0
cgms      if(iptot.gt.0) then
cgms        dlnpdr_loc=-(LOG(ptot_exp(j-1))-LOG(ptot_exp(j+jp)))/
cgms     &             (rmin_exp(j-1)-rmin_exp(j+jp))*a_unit_exp+1.D-10
        dlnpdr_loc=-(a_unit_exp/ptotm)*(ptot_exp(j+1)-ptot_exp(j))/
     >              (rmin_exp(j+1)-rmin_exp(j))
        if(ABS(dlnpdr_loc).lt.1.D-10)dlnpdr_loc=1.D-10      
cgms      endif
      write(*,*) 'dlnpdr_loc =',dlnpdr_loc
c
c JK      alpha_mhd_loc=0.0 (defaulted in gks_defaults)
c
cgms      beta_loc=1.0
cgms      beta_loc_0=403.D0*(nem*tem+nim*tim)/
cgms     &           (1.D5*bt_mag_center**2*b_div**2)/ab_div
cgms      if(iptot.gt.0) then
        beta_loc_0=0.0
        beta_loc=(1.6022D-4)*(8.0*pi)*ptotm/
     &     (b_unit_loc_s(j)**2*b_div**2)/ab_div 
cgms      endif
c      write(*,*) 'beta_loc = ',beta_loc
c
c toroidal field with no plasma
c  
c gks collisionality (xnu/w_star_i)*(ky*rho_i)
      local_vnewk3=0.117D0*ne_exp(j)*te_exp(j)**(-1.5)/
     &             (ti_exp(j)**0.5)*(a_unit_exp)*
     &             (amassgas_exp/2.)**0.5
cgms correction should be (xnu/w_star_i)*(ky*rho_i)/2=xnu*a/sqrt(2*ti/mi)
cgms put in missing 1/2 
      local_vnewk3=local_vnewk3/2.
c
cgms new calculation of vnewk3 with lnlamda factor 
cgms lnlamda and taue from NRL formulary
cgms note: for Te=1Kev, ne=10**13 lnlamda = 15.94 and taue=1.088D-3/lnlamda
      lnlamda = 15.94D0-0.5*LOG(nem)+LOG(tem)
      taue = 1.088D-3*(tem**1.5)/(nem*lnlamda)
cgms xnuei = 3/4 sqrt(pi)/taue
      xnuei = 1.329D0/taue
      mui = amassgas_exp
      vthi = 9.79D5*(2.0*tim*1.D3/mui)**0.5
      local_vnewk3 = xnuei*a_unit_exp*100.D0/vthi
c
      local_vnewk3=cnewk3*local_vnewk3
c
cgms      xnu_loc =local_vnewk3/(2.*te_exp(j)/ti_exp(j))**0.5
cgms correction: should have 
      xnu_loc  =local_vnewk3*(2.0*tiote_loc)**0.5
      write(*,*) 'xnu_loc = ',xnu_loc
c      write(*,*) 'vnewk3 = ',local_vnewk3
      local_vnewk3=0.
c note: xnu_loc differs from xnu_exp because of zeff's and 
c       unit_length a_unit_exp=rmin_exp(jmaxm) vs arho_exp
c       xnu_exp(j)=xnu_loc*(zeff_exp(j)+zeff_e)*(a_unit_exp/arho_exp)
c
c
c overrides
      ne_loc=0.
      te_loc=0.
      b00_loc=0.  
c
c main call to ifwritegeo
c
      call ifwritegeo
c
c driver variables
c       
      aspectratio_loc_s(j)=aspectratio_loc
      kappa_loc_s(j)=kappa_loc
      s_kappa_loc_s(j)=s_kappa_loc
      delta_loc_s(j)=delta_loc
      s_delta_loc_s(j)=s_delta_loc
c
      shift_loc_s(j)=shift_loc
      beta_loc_s(j)=beta_loc_out
      q_loc_s(j)=q_loc
      shat_loc_s(j)=shat_loc_out
c
      dlntidr_loc_s(j)=dlntidr_loc
      dlntedr_loc_s(j)=dlntedr_loc
      dlnnidr_loc_s(j)=dlnnidr_loc
      dlnnedr_loc_s(j)=dlnnedr_loc
c
      tiote_loc_s(j)=tiote_loc
      nione_loc_s(j)=nione_loc
      xnu_loc_s(j) = xnu_loc
c
c diagnostics
c
      shat_mhd_loc_s(j)=shat_mhd_loc_out 
      alpha_mhd_loc_s(j)=alpha_mhd_loc_out
      alpha_c_s(j)=alpha_c
      alpha_c_exp(j)=drhodr(j)*
     &               q_exp(j)**2*rmaj_exp(j)/arho_exp*
     &               betae_exp(j)*((ti_exp(j)*ni_exp(j)/
     &               te_exp(j)/ne_exp(j))*
     &               (zpni_exp(j)+zpti_exp(j))+
     &               zpne_exp(j)+zpte_exp(j))
c
c should be same as alpha_c_s
c same as "old" alpha_exp(j)/(sqrt(elong_exp(j))/((1.+elong_exp(j)**2)/2.))
c *drhodr(j)
c ...and before 2/27/98 mlt alpha_exp used ne for ni
c
      beta_loc_0_s(j)=beta_loc_0
      volume_loc_s(j)=volume_loc
      b2_ave_geo_s(j)=b2_ave_geo
      b_norm_loc_s(j)=b_norm*b_unit_loc_s(j)
c
      delta_dor_s(j)=delta_dor
      d_prime_dor_s(j)=d_prime_dor
      k_prime_dor_s(j)=k_prime_dor
c
c overwrite ifwritegeo version of beta with  
      beta=403.D-5*(nem*tim/b_unit_loc_s(j)**2)
      beta=cbetae*beta
c   debye length/rhos
      debyelorhos=cdebye*(7.43D2/(nem*1.D13)**0.5)/
     & ((1.02D2*amassgas_exp**0.5)/(b_unit_loc_s(j)*1.D4))            
c
      uprim1=sqrt(2.0*tem/tim)*gamma_p_i_exp(j)/
     &       (arho_exp/a_unit_exp)
      uprim2=sqrt(2.0*tem/tim)*gamma_p_i_exp(j)/
     &       (arho_exp/a_unit_exp)
      uprim3=sqrt(2.0*tem/tim)*gamma_p_e_exp(j)/
     &       (arho_exp/a_unit_exp)
      egamma = alpha_e*egamma_exp(j)
      mach1 = mach_i_exp(j)
      mach2 = mach_i_exp(j)
      mach3 = mach_e_exp(j)
c
      return
      end
c@gks_max.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c     This subroutine varies tprim3 to find the marginal
c     stability point. One call to gstotal must be made before.
c     Must set xkymax_gf,xkymin_gf to maximim and minimum
c     range over which a maximim in k will be sought.
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine find_gradte_crit
c
      use gks_var
      implicit none
      include 'glf.m'
      include 'data_exp.m'
c
      integer i, ierr
      real*8 tprim_bot, tprim_top, tprim_save, aky1_save
c
      tprim_save=tprim3
      aky1_save=aky1
      tprim_top=-999.D0
      tprim_bot=-999.D0
      ierr=0
c
      write(*,*) ' start'
      write(*,*) ' agammas(1) = ',agammas(1)
      write(*,*) ' dgammas(1) = ',dgammas(1)
      write(*,*) ' tprim3 = ',tprim3
c
c   find an unstable upper bound
c
      if (dgammas(1).gt. tol_f*agammas(1)) then
        do i=1,10
           tprim3=(1.0+0.2)*tprim3
           call gstotal
           write(*,*) '   up, i = ',i
           write(*,*) '   agammas(1) = ',agammas(1)
           write(*,*) '   dgammas(1) = ',dgammas(1)
           write(*,*) '   tprim3 = ',tprim3
           write(*,*) '   tprim_bot = ',tprim_bot
           if (dgammas(1).lt.tol_f*agammas(1)) exit
c       until((i.eq.6).or.(dgammas(1).lt.tol_f*agammas(1)))
        enddo
        if (dgammas(1).lt.tol_f*agammas(1)) tprim_top=tprim3
      else
          tprim_top=tprim3
      endif
      if (tprim_top.eq.-999) ierr=1
c
c   maximize over k
c
      if(ierr.eq.0) call gks_max
      write(*,*) '   maximum'
      write(*,*) '   agammas(1) = ',agammas(1)
      write(*,*) '   dgammas(1) = ',dgammas(1)
      write(*,*) '   tprim3 = ',tprim3
      write(*,*) '   tprim_top = ',tprim_top
      write(*,*) '   tprim_bot = ',tprim_bot
      if (dgammas(1).lt.0.0) ierr=1
c
c  find a stable lower bound
c
      if (ierr.eq.0) then
        do i=1,10
          tprim_top=tprim3
          tprim3=tprim3/(1.0+0.3)
          call gstotal
          write(*,*) '   down, i = ',i
          write(*,*) '   agammas(1) = ',agammas(1)
          write(*,*) '   dgammas(1) = ',dgammas(1)
          write(*,*) '   tprim3 = ',tprim3
          write(*,*) '   tprim_bot = ',tprim_bot
          if (dgammas(1).gt.tol_f*agammas(1)) exit
c       until((i.eq.6).or.(dgammas(1).gt.tol_f*agammas(1)))
        enddo
        if (dgammas(1).gt.tol_f*agammas(1)) tprim_bot=tprim3
      endif
c
      if (tprim_bot.eq.-999) ierr=1
c
c   find critical point
c
      if (ierr.eq.0) then
        do i=1,8
          tprim3=(tprim_top+tprim_bot)/2
          call gstotal
          write(*,*) '   find crit, i = ',i
          write(*,*) '   agammas(1) = ',agammas(1)
          write(*,*) '   dgammas(1) = ',dgammas(1)
          write(*,*) '   tprim3 = ',tprim3
          write(*,*) '   tprim_top = ',tprim_top
          write(*,*) '   tprim_bot = ',tprim_bot
          if (dgammas(1).gt.tol_f*agammas(1)) then
            tprim_bot=tprim3
          else
            tprim_top=tprim3
          endif
          if (abs(tprim_top-tprim_bot).le.0.02*abs(tprim3)) exit
c       until((i.eq.8).or.(abs(tprim_top-tprim_bot).le.0.02*abs(tprim3)))
        enddo
      endif
c
      if(ierr.ne.0)then
        tprim3=tprim_save
        aky1=aky1_save
        kys(1)=-kys(1)
      endif
c
      return
      end
c@gradte_crit.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c    This subroutine finds the critical grad-Te profile
c    The experimental profile must be loaded first.
c    Here, kys0 is used as described in gks_profile
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine gradte_crit(jstart_prof,jend_prof,jdel)
c
      use gks_var
      implicit none
      include 'input.m'
      include 'glf.m'
      include 'gks_out.m'
      include 'data_exp.m'
c
      integer i, j, dj, jstart_prof, jend_prof, jdel
      real*8 save_tprim3, last_tprim3
c
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c  begin loop over grid
c
      write(*,*) 'Finding gradte-crit ...'
c
      if (jstart_prof .le. jend_prof) then
         dj = jdel
      else
         dj = -jdel
      endif
      last_tprim3 = -1.D0
      kys0 = 0.42/SQRT(amass3)  ! start at ky*rho_e = 0.42
c
      do j=jstart_prof,jend_prof,dj
       write(*,*)"igeo=",igeo
       if (igeo.eq.0) call gks_profile(j,j)
       if (igeo.eq.1) call geogks_profile(j,j,1)
       if(kys0*debyelorhos.ge.0.2)then
c adjust starting wavenumber down due to debye cut-off
         kys0 = 0.2/debyelorhos
         write(*,*)"kys0=",kys0
         if (igeo.eq.0) call gks_profile(j,j)
         if (igeo.eq.1) call geogks_profile(j,j,1)
       endif
       save_tprim3 = tprim3
       if (last_tprim3.gt.0) tprim3=last_tprim3
       if (save_tprim3.eq.0.0) save_tprim3 = 1.0e-10
       call find_gradte_crit
       write(*,*) ' tprim3 = ', tprim3
       write(*,*) ' kys(1) = ', kys(1)
       write(*,*) ' save_tprim3 = ', save_tprim3
       if (kys(1).gt.0.0) kys0=kys(1)
       zpte_crit_m(j) = (tprim3/save_tprim3)*zpte_exp(j)
       ky_m(j) = kys(1)
       last_tprim3 = tprim3
       if (kys(1).lt.0) last_tprim3=-1
      enddo
c
      write(*,*) ' shot = ', shot
c     write(*,*) ' time = ', time_d
c      write(*,*)
c      do j=jstart_prof,jend_prof,dj
c        write(*,100) j, rho(j), zpte_exp(j), 
c     &               zpte_crit_m(j), ky_m(j)
c      enddo
c
 100  format(2x,i2,2x,0p2f9.5,1p6e14.6)
c
      return
      end

c***********************************************************
      subroutine find_gradti_crit
c***********************************************************
c     This function varies tprim1 to find the marginal
c     stability point.One call to gstotal must be made before.
c     Must set xkymax_gf,xkymin_gf to maximim and minimum
c     range over which a maximim in k will be sought.
c     if z3 > 0 finds critical gradient for adiabatic electrons and zero beta
c     if z3 < 0 finds gradient at which frequency becomes positive or zero
c     this sould be the place where the ITG becomes the TEM 
c     Note that tprim2=1 sets impurity temperature gradient to that of main ion.
c***********************************************************
      use gks_var
      implicit none
      include 'glf.m'
      include 'gks_out.m'
      include 'data_exp.m'
       real*8 tprim_bot,tprim_top
       real*8 tprim_save,aky1_save,sign_freq
       integer i,ierr
       tprim_save=tprim1
       aky1_save=aky1
       tprim_top=-999.D0
       tprim_bot=-999.D0
       ierr=0
       sign_freq=-1.D0
       write(*,*)"start",agammas(1),dgammas(1),afreqs(1),tprim1
       if(z3.lt.0.D0.and.afreqs(1).ge.0.D0)then
        dgammas(1)= 2.0*tol_f*abs(agammas(1))
        sign_freq=1.0
       endif
       if(dgammas(1).gt. tol_f*agammas(1)) then
c   find an unstable upper bound
        do i=1,10
         tprim1=(1.2D0)*tprim1
         call gstotal
         if(z3.lt.0.D0.and.afreqs(1).ge.0.D0)then
          dgammas(1)= 2.0*tol_f*abs(agammas(1))
          sign_freq=1.0
         endif
        write(*,*) "up",i,agammas(1),dgammas(1),afreqs(1)
        write(*,*)tprim1,tprim_bot
c        until((i.eq.10).or.(dgammas(1).lt.tol_f*agammas(1)))
          if(dgammas(1).lt.tol_f*agammas(1))exit
        enddo
        if(dgammas(1).lt.tol_f*agammas(1))tprim_top=tprim1
       else
        tprim_top=tprim1
       endif
       if(tprim_top.eq.-999)ierr=1
c   maximize over k
       if(ierr.eq.0)call gks_max
       write(*,*)"maximum",agammas(1),dgammas(1),afreqs(1)
       write(*,*)tprim1,tprim_top,tprim_bot
       if(dgammas(1).lt.0.0)ierr=1
       if(ierr.eq.0)then
        if(z3.lt.0.0.and.afreqs(1).ge.0.0)then
         write(*,*)"maximum k found electron mode"
         dgammas(1)= 2.0*tol_f*abs(agammas(1))
         sign_freq=1.0
c refind the upper bound with a negative frequency
         tprim_top=-999.0
         do i=1,10
          tprim1=(1.2D0)*tprim1
          call gstotal
          if(z3.lt.0.0.and.afreqs(1).ge.0.0)then
           dgammas(1)= 2.0*tol_f*abs(agammas(1))
           sign_freq=1.0
          endif
         write(*,*)"up",i,agammas(1),dgammas(1),afreqs(1)
         write(*,*)tprim1,tprim_bot
c         until((i.eq.10).or.(dgammas(1).lt.tol_f*agammas(1)))
          if(dgammas(1).lt.tol_f*agammas(1))exit
         enddo
         if(dgammas(1).lt.tol_f*agammas(1))tprim_top=tprim1
         if(tprim_top.eq.-999)ierr=1
         endif
        endif
        if(ierr.eq.0)then
c  find a stable lower bound
        do i=1,10
         tprim_top=tprim1
         tprim1=tprim1/(1.0+0.3)
         call gstotal
         if(z3.lt.0.0.and.afreqs(1).ge.0.0)then
          dgammas(1)= 2.0*tol_f*abs(agammas(1))
          sign_freq=1.0
         endif
        write(*,*)"down",i,agammas(1),dgammas(1),afreqs(1)
        write(*,*) tprim1,tprim_bot
c        until((i.eq.10).or.(dgammas(1).gt.tol_f*agammas(1)))
        if(dgammas(1).gt.tol_f*agammas(1))exit
        enddo
        if(dgammas(1).gt.tol_f*agammas(1))tprim_bot=tprim1
       endif
       if(tprim_bot.eq.-999)ierr=1
       if(ierr.eq.0)then
c   find critical point
        do i=1,10
         tprim1=(tprim_top+tprim_bot)/2
         call gstotal
         if(z3.lt.0.0.and.afreqs(1).ge.0.0)then
          dgammas(1)= 2.0*tol_f*abs(agammas(1))
          sign_freq=1.0
         endif
         write(*,*)"find crit",i,agammas(1),dgammas(1),afreqs(1)
         write(*,*) tprim1,tprim_top,tprim_bot
         if(dgammas(1).gt. tol_f*agammas(1))then
          tprim_bot=tprim1
         else
          tprim_top=tprim1
         endif
         if(abs(tprim_top-tprim_bot).le.0.02*abs(tprim1))exit
c        until((i.eq.10).or.(abs(tprim_top-tprim_bot).le.0.02*abs(tprim1)))
        enddo
       endif
       afreqs(1)=sign_freq
       if(ierr.ne.0)then
        tprim1=tprim_save
        aky1=aky1_save
        kys(1)=-kys(1)
       endif
       return
       end
c***********************************************************
       subroutine gradti_crit(jstart_prof,jend_prof,jdel)
c***********************************************************
c         this function finds the critical gradti profile for
c         electrostatic ITG modes with adiabatic electrons.
c         The experimental profiles must be loaded
c         kys0 is used as described in gks_profile.
c***********************************************************
      use gks_var
      implicit none
      include 'input.m'
      include 'glf.m'
      include 'gks_out.m'
      include 'data_exp.m'
      integer i,j,dj
      integer jstart_prof,jend_prof,jdel
      real*8 save_tprim1,last_tprim1
      real*8 save_cbetae,save_z3,save_nstep
      if(z3 .gt. 0.0)then
c  GKS uses adiabatic electrons for z3 > 0
c        save_z3 = z3
c        z3 = 1
c  Electrostatic
        save_cbetae = cbetae
        cbetae = 1.0D-33
      endif
      save_nstep=nstep
c  begin loop over grid
      if (jstart_prof .le. jend_prof) then
         dj = jdel
      else
         dj = -jdel
      endif
      last_tprim1 = -1.D0
      kys0=0.D0    ! use aky1 for starting wavenumber
      j = jstart_prof
      do j = jstart_prof,jend_prof,dj
        if(igeo.eq.0)call gks_profile(j,j)
        if(igeo.eq.1)call geogks_profile(j,j,1)
        save_tprim1 = tprim1
        if(last_tprim1 .gt. 0.0)tprim1=last_tprim1
        if(save_tprim1 .eq. 0.0)save_tprim1 = 1.0e-10
        nstep=save_nstep
        call gstotal
        call find_gradti_crit
        write(*,*) 'tprim1 = ',tprim1
        write(*,*) 'kys(1) = ',kys(1)
        write(*,*) 'save_tprim1 = ',save_tprim1
        write(*,*) 'afreqs(1) = ',afreqs(1)
        if(kys(1).gt.0.0)kys0=kys(1)
        zpti_crit_m(j) = (tprim1/save_tprim1)*zpti_exp(j)
        ky_m(j) = kys(1)
        anfreq_m(j) = afreqs(1)
        last_tprim1 = tprim1
        if(kys(1).lt.0.D0)last_tprim1=-1
       enddo
c      until (dj*(j-jend_prof) .gt. 0.0)
c      output "gradti_crit.dat"
c      shot,time_d,tok,arho_exp,z3,cbetae,zpti_crit,ky_m
c      anfreq_m
c      output tty
c      z3 = save_z3
      if(z3.gt.0.D0)cbetae = save_cbetae
      return
      end

