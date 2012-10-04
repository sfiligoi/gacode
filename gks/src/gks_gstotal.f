      subroutine gks
************************************************************************
!      use gks_var
!      implicit none
cgms      Use(Ifwcom)
cgms      Use(Basiscom)
cgms      Use(Common)
crew_______________geo__________________________________________________
cgms      Use(Geo)
       

c 1/21/98 bpar added for i_bpar=1 flag.
c note phi= e*Phi/T_1   apar=e*A_par/(cT_1)*vth_1 
c      bpar = e**2*B(theta)*sigma/(c**2*M_1*T_1) 
c      where B_par=k_perp**2*sigma
c    see Antonsen and Lane, PF 23(6),(1980) 1205 gyprkinetic equation
c    with Phi, A_par, and B_par...here assume T_par=T_perp

c 9/01/98 new coding bpar=B_par/B(theta)
c   in new coding zJ1(z) is replaced by J1(z)/z
c also
c   in new coding B_par is lumped with phi not apar
c   so there is a polarization part associated with B_par in the QN
c and there is added polarization term in the perp amperes law 
c like J1**2 in addition to the J1*J0 part

c  i_bpar=0 go through old coding only
c  i_bpar=-1 go through old coding except go through timeadv_bpar
c     with bpar in mhd limit directly related to phi if y_mhd=1.
c  i_bpar=1 go through new coding timeadv_bpar, matrix_bpar, getfield_bpar
c     which uses ordering 
c     phi(ntgridl:ntgrid),apar(ntgridl:ntgrid),bpar(ntgridl:ntgrid) 
c  i_bpar=2 go through new coding timeadv_bpar, matrix_bpar2, getfield_bpar2
c     which uses ordering
c     [phi,apar,bpar](ntgridl:ntgrid)
c  x_bpar=1. full  x_bpar-> 0 forces bpar->0
c  y_bpar=1. full  y_bpar-> 0 drops bpar from eq. of motion
c  xy_bpar=1. full xy_bpar->0 kills bpar polarization terms
c  y_bpar->-1 and xy_bpar->-1 flips bpar in Antonsen/Lane EoM
c  x_bpar->-1 flips perp ampers law in Antonsen/Lane
c  y_mhd=1. b_par->del_bpar=bpar_mhd-b_bpar(true)


       character*20 gy0inp,gyinp,gyout
crew       real t(3,3),ratio(3,3),akyarr(16)
crew       complex alpha,alpha1,alp1,al1a1
crew       complex phim(0:ntm),aparm(0:ntm)

cnopar
       complex amphi
       real test,phi_max
       
       
crews       open(1,file='name',status='old')
crews       read(1,5) gy0inp
crews       read(1,5) gyinp
crews       read(1,5) gyout
crews 5     format(20A)
     
crews       open(10,file=gy0inp,status="old")
crews       open(9, file=gyinp, status="old")
crews       open(8, file=gyout, status="unknown")
       open(99, file='tempdat', status="unknown")

    
crews       if (noread.eq.0) then
       
cnopar iphi
  
crews       read(10,*) nscreen, nout, nstep, ngamstep, isvar, nwrite, iphi
crews       read(10,*) delt, ominst, phydif, power, tol
crews       read(10,*) gamma, absom, fcv, fv, gridfac
crews       read(10,*) fexp1, fexp3, test1, test2
crews       read(10,*) alr, ali, al1r, al1i
crews       read(10,*) alar, alai, al1ar, al1ai
crews       read(10,*) ncspec1,ncspec2,ncspec3,ncspec4,ncspec5
crews       read(10,*)    nce1,   nce2,   nce3,   nce4,   nce5
crews       read(10,*) naky
crews       read(10,*) aky1, aky2, aky3, aky4
crews       read(10,*) aky5, aky6, aky7, aky8
crews       read(10,*) aky9, aky10,aky11,aky12
crews       read(10,*) aky13,aky14,aky15,aky16
crews       read(10,*) omegacom     
cnopar 
crews       read(10,*) nt0, thetamin, thetamax, nt0min, nt0max
crews       read(10,*) ecut

cnopar negrid, uprim1,2,3,4,5 

crews       read(9,*) ntheta, nperiod, ngauss, icv, negrid
crews       read(9,*) eps, shift, dbeam, shat, pk
crews       read(9,*) epsl, width0, beta, zeff, teti, zeff5
crews       read(9,*) uprim1, uprim2, uprim3, uprim4, uprim5
crews       read(9,*) fprim1, fprim2, fprim3, fprim4, fprim5
crews       read(9,*) tprim1, tprim2, tprim3, tprim4, tprim5
crews       read(9,*) vnewk1, vnewk2, vnewk3, vnewk4, vnewk5
crews       read(9,*) bakdif1,bakdif2,bakdif3,bakdif4,bakdif5
crews       read(9,*)         amass2, amass3, amass4, amass5
crews       read(9,*)          temp2,  temp3,  temp4,  temp5
crews       read(9,*)             z2,     z3,     z4,     z5
  
crews       endif

       time=0.
       omegacom = 0.0
       zi=cmplx(0.,1.)
cgms       pi= 3.1415926
       pi=atan2(0.0,-1.0)
cjc       c= cexp( -zi*omegacom*time )
       c = exp(-zi*omegacom*time)

       nce(1)=nce1
       nce(2)=nce2
       nce(3)=nce3
       nce(4)=nce4
       nce(5)=nce5

       ncspec(1)=ncspec1
       ncspec(2)=ncspec2
       ncspec(3)=ncspec3
       ncspec(4)=ncspec4
       ncspec(5)=ncspec5

       z(1)=1.
       z(2)=z2
       z(3)=z3
       z(4)=z4
       z(5)=z5

       amass(1)=1.
       amass(2)=amass2
       amass(3)=amass3
       amass(4)=amass4
       amass(5)=amass5

       temp(1)=1.
       temp(2)=temp2
       temp(3)=temp3
       temp(4)=temp4
       temp(5)=temp5

       bakdif(1)=bakdif1
       bakdif(2)=bakdif2
       bakdif(3)=bakdif3
       bakdif(4)=bakdif4
       bakdif(5)=bakdif5
       
cnopar
       uprim(1)=uprim1
       uprim(2)=uprim2/SQRT(temp2/amass2)
       uprim(3)=uprim3/SQRT(temp3/amass3)
       uprim(4)=uprim4/SQRT(temp4/amass4)
       uprim(5)=uprim5/SQRT(temp5/amass5)

       fprim(1)=fprim1
       tprim(1)=tprim1
       fprim(2)=fprim2*fprim1
       tprim(2)=tprim2*tprim1
       fprim(3)=fprim3*fprim1
       tprim(3)=tprim3
       fprim(4)=fprim4
       tprim(4)=tprim4
       fprim(5)=fprim5*fprim1
       tprim(5)=tprim5*tprim1

       alpha=cmplx(alr,ali)
       alpha1=cmplx(al1r,al1i)
       ala=cmplx(alar,alai)
       al1a=cmplx(al1ar,al1ai)

       do 10 is=1,nspec
          fprimm(is)=fprim(is)
 10    continue

       akyarr(1)=aky1
       akyarr(2)=aky2
       akyarr(3)=aky3
       akyarr(4)=aky4
       akyarr(5)=aky5
       akyarr(6)=aky6
       akyarr(7)=aky7
       akyarr(8)=aky8
       akyarr(9)=aky9
       akyarr(10)=aky10
cshort       akyarr(11)=aky11
cshort       akyarr(12)=aky12 
cshort       akyarr(13)=aky13
cshort       akyarr(14)=aky14
cshort       akyarr(15)=aky15
cshort       akyarr(16)=aky16

crew       do 10 ii=1,3
crew       do 10 jj=1,3
crew              t(ii,jj)= 0.
crew          ratio(ii,jj)= 0.
crew 10    continue


crews       if(noread.eq.0) call openfile
       dtheta = 0.D0      
       if (nt0.ne.0) dtheta= (thetamax -thetamin)*pi / nt0


c------------------------------------
c   starting aky loop and gamma loop
c------------------------------------
       do 900 iaky=1,16
         if(iaky.gt.naky) goto 1000
	 
crew________________________________
crew    kys in units of rhos(1) not rhoi(1)

         ikys=iaky
       kysi=akyarr(iaky)/sqrt(2.)*sqrt(temp(3))
       kys(ikys)=kysi
	 
  
       do 11 ii=1,3
       do 11 jj=1,3
          t(ii,jj)= 0.
          ratio(ii,jj)= 0.
 11    continue
crew________________________________

       do 800 it0= nt0min, nt0max
          theta0= thetamin*pi + it0 * dtheta

       do 700 igam=1,ngamstep

         time=0.
	 c= cexp( -zi*omegacom*time )

         if (igam.eq.1) gammause= gamma
         if (igam.eq.2) gammause= gamma/2.
         
         aky=akyarr(iaky) 


crews         write(8,*) ' ======OUTPUT FOR aky(',iaky,')=',aky,'========'
crews         write(8,*) ' -------OUTPUT FOR theta0=',theta0,'---------'

         vnew0(1)=vnewk1/aky
         vnew0(2)=vnewk2/aky
         vnew0(3)=vnewk3/aky
         vnew0(4)=vnewk4/aky
         vnew0(5)=vnewk5/aky

cgms note that this is really vnew0=(1/2)vnewk3*(2/aky)
cgms the factor (2/aky) makes vnewk3*2/aky=nu_ei/wstar_i
cgms the extra 1/2 out front is needed for the pitch angle 
cgms scattering operator. 

         pp=pk/aky
	 
crew________geo________________________________________________________

         do j=-ntgridl,ntgrid
	      ppgeo(j)=pk_geo(j)/aky
	     enddo 
 
crew________geo________________________________________________________
	 

c    set up coordinates and integration factors
         call gridset

         icount=0

         tcomp=0.0
cgms         omm1=cmplx(1.0D50,1.0D50)         
         fexp(1)=fexp1
         fexp(2)=fexp1
         fexp(3)=fexp3
         fexp(4)=fexp3
         fexp(5)=fexp1
         
         if(nscreen.eq.0) goto 1000
	 
cnpar   reverses calls to ginitial and plasmacoef?
	 
c   set up initial conditions
         if (icount.eq.0.and.icontinue.ne.2) call ginitial

c   set up plasma coefficients
         if(icontinue.ne.2) call plasmacoef



         
         icount= 1
         alp1= 0.
         if(z(3).lt.0.) alp1=alpha

         istep=0
          
         do mi=1,ntmp2p1t3
          do mj=1,ntmp2p1t3
           am(mi,mj)=0.
          enddo
         enddo

         if(icontinue.ne.2.and.i_bpar.le.0) call matrix
         if(icontinue.ne.2.and.i_bpar.eq.1) call matrix_bpar
         if(icontinue.ne.2.and.i_bpar.eq.2) call matrix_bpar2
         
         do 31 is=1,nspec
            if(ncspec(is).eq.0) go to  31

            do 30 isign=1,2
            do 30 ie=1,negrid
            do 30 il=1,nlambda
            do 30 j=-ntgridl,ntgrid
               gnew(j,il,ie,isign,is)= g(j,il,ie,isign,is)
 30         continue

 31      continue


         do 40 j=-ntgridl,ntgrid
             phig(j)= 0.
            aparg(j)= 0.
            bparg(j)= 0.
  40     continue

         alphal= alp1
         alphapl= 0.

         if(i_bpar.le.0) call getfield
         if(i_bpar.eq.1) call getfield_bpar
         if(i_bpar.eq.2) call getfield_bpar2

	 call getterm
	 
cnopar loop 51
         do 51 is= 1,nspec
            if(ncspec(is).eq.0) goto 51

            do 50 ie= 1,negrid
               termi(ie,is)= term(ie,is)
               
               do 45 ig= -ntgridl,ntgrid
                  ginti(ig,ie,is)= gint(ig,ie,is)
 45            continue

 50         continue

 51      continue
 
         do 60 j=-ntgridl,ntgrid
             phig(j)=  phi(j)
            aparg(j)= apar(j)
            bparg(j)= bpar(j)

             phim(j)=  phig(j)
            aparm(j)= aparg(j)
            bparm(j)= bparg(j)
 60      continue

crew
      if( icontinue.ge.1) then
cgms added rescaling factor 1.0D-10 for restarts
        phi_max = 1.0
        do j=-ntgridl,ntgrid
          test = abs(real(phi(j)*conjg(phi(j))))
          if(test.gt.phi_max)phi_max=test
        enddo
        phi_max = sqrt(phi_max)
        write(*,*)"phi_max=",phi_max
c
       do is=1,nspec
         do ie=1,negrid
          do i=1,nlambda
           do j=-ntgridl,ntgrid
              g(j,i,ie,1,is)=gsav(j,i,ie,1,is)/phi_max
              g(j,i,ie,2,is)=gsav(j,i,ie,2,is)/phi_max
              gnew(j,i,ie,1,is)=g(j,i,ie,1,is)
              gnew(j,i,ie,2,is)=g(j,i,ie,2,is)
           enddo
          enddo
         enddo
        enddo
        do j=-ntgridl,ntgrid

             phi(j)=phisav(j)/phi_max
            apar(j)=aparsav(j)/phi_max
            bpar(j)=bparsav(j)/phi_max

             phig(j)=  phi(j)
            aparg(j)= apar(j)
            bparg(j)= bpar(j)

             phim(j)=  phig(j)
            aparm(j)= aparg(j)
            bparm(j)= bparg(j)
        enddo
      endif
crew
         ict= 0

         if (igam.eq.1) then
            tguess1= 2.5*tprim(1)
            tguess3= 2.5*tprim(3)
         endif

         if (igam.eq.2) then
           tguess1= t(1,1)
           tguess3= t(1,3)
         endif
 
        alphal= alpha
        alphapl= alpha1
             
****************main loop*************
crew
         timest(0)=0.
         do 500 istep=1,nstep
            
           icount=icount+1
           time=delt+time
cnpar
           c= cexp( -zi*omegacom*time )
           
	   
crew
           timest(istep)=time
 
           if(power.gt.0.) then 
              call wstarch

crew 3/2/98

            tp1_rec(istep,igam)=tp1
            tp3_rec(istep,igam)=tp3 
            if(istep.eq.nstep) then
             nstepd4=nstep/4
             tpave=0.
             tpcount=0
              do istep_cnt=nstepd4,nstep
               tpcount=tpcount+1
               tpave=tpave+tp1_rec(istep_cnt,igam)
              enddo
              tp1=tpave/tpcount
             tpave=0.
             tpcount=0
              do istep_cnt=nstepd4,nstep
               tpcount=tpcount+1
               tpave=tpave+tp3_rec(istep_cnt,igam)
              enddo
              tp3=tpave/tpcount 
            endif
crew
              t(igam,1)= tp1
              t(igam,3)= tp3
           endif


crew added omegaprev
      omega0=clog((phig(0)+aparg(0)+bparg(0))/
     &    (phim(0)+aparm(0)+bparm(0)))
     &            *cmplx(0.,1.)/delt
     >            +omegaprev

           if( ((( (istep-1)/nscreen)*nscreen).eq.(istep-1)).and.
     &         (istep.ne.1) )   absom=cabs(omega0)*fcv+(1.-fcv)*absom

           ivnewout=0
           if((( (istep-1)/nscreen)*nscreen).eq.(istep-1)) ivnewout=1
 
crew__________________________________________________
      do j=0,ntgrid-1
crew added omegaprev
      omega=clog((phig(j)+aparg(j)+bparg(j))/
     &   (phim(j)+aparm(j)+bparm(j)))
     &             *cmplx(0.,1.)/delt
     >             +omegaprev

              phim(j)= phig(j)
             aparm(j)=aparg(j)
             bparm(j)=bparg(j)

crew freqs  units of ws=c_s/a  (+) for electron (-) for ion
crew gammas units os ws=c_s/a
crew c_s=v_th_i/(sqrt(2*T_i/T_e)

       freqsj(j,ikys)=real(omega)*kysi/temp(3)
       gammasj(j,ikys)=aimag(omega)*kysi/temp(3)

      enddo
      jcntst=0
      jcntskip=ntgrid/10+1
      do j=0,ntgrid-1,jcntskip
       freqst(istep,jcntst,ikys)=freqsj(jcntst,ikys)
       gammast(istep,jcntst,ikys)=gammasj(jcntst,ikys)
       jcntst=jcntst+1
      enddo
  
      stkgammast(istep,iistk,ikys)=gammast(istep,0,ikys)
      stknstep(ikys,iistk)=istep
      omm=1.e6
      if(istep.gt.nteststep) omm=
     >abs(gammast(istep,0,ikys)-gammast(istep-nteststep,0,ikys))
     > /abs(gammast(istep,0,ikys))            

      if(omm.lt.tol) goto 700
crew_____________________________________________________

           do 100 j=0,ntgrid-1
crew             omega=clog((phig(j)+aparg(j))/(phim(j)+aparm(j)))
crew     &             *cmplx(0.,1.)/delt
crew             omm= abs(real(omega)-real(omm1))
crew     &           +abs(aimag(omega)-aimag(omm1))
     

crew              phim(j)= phig(j)
crew             aparm(j)=aparg(j)

             if((((istep-1)/nscreen)*nscreen).ne.(istep-1)) goto 100

crews             write(8,*) " j=",j,"  phig=",phig(j)," aparg=",aparg(j)
crews             write(8,*) " j=",j," omega=",omega

  100 continue

         if((((istep-1)/nscreen)*nscreen).ne.(istep-1)) goto 300



crews         call flux
crews         write(8,*) " qheat=",qheat(1),qheat(2)
crews         write(8,*) " pflux=",pflux(1),pflux(2)

         call kperp
         dmixing=aimag(omega0)*aky/akperp**2
         diffren=diff/aky
crews         write(8,*) " dmixing=",dmixing,"diffren=",diffren
	 
cnopar add loop 151
           do 151 is= 1, nspec
              if(ncspec(is).eq.0) goto 151

              do 150 ie= 1,negrid
                 if(nce(ie).eq.0) goto 150
                 
crew                 write(8,*) " term(",ie,",",is,")=",term(ie,is)
 150          continue

 151       continue
crew         write(8,*) " tp1=",tp1," time=",time

crewcc         if(omm.lt.tol) goto 9000
crew         omm1=omega

 300     continue
         
         do 400 j=-ntgridl,ntgrid
             phi(j)= 0.
             apar(j)= 0.
             bpar(j)= 0.
 400     continue
 
        if(i_bpar.le.0) then
         if(i_bpar.eq.0) call timeadv
         if(i_bpar.eq.-1) call timeadv_bpar
         call getfield
         if(i_bpar.eq.0) call timeadv
         if(i_bpar.eq.-1) call timeadv_bpar
        endif
        if(i_bpar.eq.1) then
         call timeadv_bpar
         call getfield_bpar
         call timeadv_bpar
        endif
        if(i_bpar.eq.2) then
         call timeadv_bpar
         call getfield_bpar2
         call timeadv_bpar
        endif            
        
         call replace
cnopar
         call getterm
crew         call output

         do 420 j=-ntgridl,ntgrid
             phig(j)= phi(j)
            aparg(j)=apar(j)
            bparg(j)=bpar(j)
 420     continue
        if(istep.gt.10) call current
crew____________________________________________________

 500     continue


        do 600 ig=0,ntgrid
           phig(ig)= phi(ig)
           aparg(ig)=apar(ig)
 600      bparg(ig)=bpar(ig)
 
 700   continue


 800   continue   


       if (power.gt.0.) then 

         if ((ngamstep.eq.1).and.(isvar.eq.1)) then
            ratio(1,1)= tprim(1)/t(1,1)
         endif

         if ((ngamstep.eq.1).and.(isvar.eq.3)) then
            ratio(1,3)= tprim(3)/t(1,3)
         endif

         if ((ngamstep.eq.1).and.(isvar.eq.13)) then
            ratio(1,1)= tprim(1)/t(1,1)
            ratio(1,3)= tprim(3)/t(1,3)
         endif


         if ((ngamstep.eq.2).and.(isvar.eq.1)) then
            t(3,1)= t(2,1) + (t(2,1)-t(1,1))

            ratio(1,1)= tprim(1)/t(1,1)
            ratio(2,1)= tprim(1)/t(2,1)
            ratio(3,1)= tprim(1)/t(3,1)
         endif

         if ((ngamstep.eq.2).and.(isvar.eq.3)) then
            t(3,3)= t(2,3) + (t(2,3)-t(1,3))

            ratio(1,3)= tprim(3)/t(1,3)
            ratio(2,3)= tprim(3)/t(2,3)
            ratio(3,3)= tprim(3)/t(3,3)
         endif

         if ((ngamstep.eq.2).and.(isvar.eq.13)) then
            t(3,1)= t(2,1) + (t(2,1)-t(1,1))
            t(3,3)= t(2,3) + (t(2,3)-t(1,3))

            ratio(1,1)= tprim(1)/t(1,1)
            ratio(2,1)= tprim(1)/t(2,1)
            ratio(3,1)= tprim(1)/t(3,1)

            ratio(1,3)= tprim(3)/t(1,3)
            ratio(2,3)= tprim(3)/t(2,3)
            ratio(3,3)= tprim(3)/t(3,3)
         endif

crew___________________________________________

      tp1crt(ikys)=2.*t(2,1)-t(1,1)
      tp3crt(ikys)=2.*t(2,3)-t(1,3)
      stktp1crt(ikys,iistk)=tp1crt(ikys)
      stktp3crt(ikys,iistk)=tp3crt(ikys)
crew___________________________________________

       endif

c       write( 6,*) "*****************"
c       write( 6,*) "t(1,1)=",t(1,1),"  t(1,3)=",t(1,3)
c       write( 6,*) "t(2,1)=",t(2,1),"  t(2,3)=",t(2,3)
c       write( 6,*) "t(3,1)=",t(3,1),"  t(3,3)=",t(3,3)
c       write( 6,*) "-----------------"
c       write( 6,*) "ratio(1,1)=",ratio(1,1),"  ratio(1,3)=",ratio(1,3)
c       write( 6,*) "ratio(2,1)=",ratio(2,1),"  ratio(2,3)=",ratio(2,3)
c       write( 6,*) "ratio(3,1)=",ratio(3,1),"  ratio(3,3)=",ratio(3,3)

crews       write(8,*) "*****************"
crews       write(8,*) "t(1,1)=",t(1,1),"  t(1,3)=",t(1,3)
crews       write(8,*) "t(2,1)=",t(2,1),"  t(2,3)=",t(2,3)
crews       write(8,*) "t(3,1)=",t(3,1),"  t(3,3)=",t(3,3)
crews       write(8,*) "-----------------"
crews       write(8,*) "ratio(1,1)=",ratio(1,1),"  ratio(1,3)=",ratio(1,3)
crews       write(8,*) "ratio(2,1)=",ratio(2,1),"  ratio(2,3)=",ratio(2,3)
crews       write(8,*) "ratio(3,1)=",ratio(3,1),"  ratio(3,3)=",ratio(3,3)
       
crew________________end of k-loop_________________
  
crew____________________________________________________
        do is=1,nspec
         do ie=1,negrid
          do i=1,nlambda
           do j=-ntgridl,ntgrid
              gsav(j,i,ie,1,is)=g(j,i,ie,1,is)
              gsav(j,i,ie,2,is)=g(j,i,ie,2,is)
           enddo
          enddo
         enddo
        enddo

         do j=-ntgridl,ntgrid
             phisav(j)=phi(j)
            aparsav(j)=apar(j)
            bparsav(j)=bpar(j)
         enddo

crew____________________________________________________
       call normalmode
       call quasilinear
       call normalmode
crew______________________________________________


 900   continue   


1000   continue

crews       close(10)
crews       close(11)
crews       close(12)
crew______________________________________________
crews       close(1)
       close(99)
crews       close(9)
crews       close(8)

       noread=1

       return
       end subroutine gks


************************************************************************
       subroutine egridset
************************************************************************

cgms      Use(Common)

       real gaus(10,10), w1(10,10)


       data gaus(8,1), gaus(8,2), gaus(8,3), gaus(8,4), gaus(8,5),
     >      gaus(8,6), gaus(8,7), gaus(8,8)
     >      / 0.01985,   0.10166,   0.23723,   0.40828,   0.59171,
     >        0.76276,   0.89833,   0.98014/

       data gaus(7,1), gaus(7,2), gaus(7,3), gaus(7,4), gaus(7,5),
     >      gaus(7,6), gaus(7,7)
     >      / 0.02544,   0.12923,   0.29707,   0.50000,   0.70292,
     >        0.87076,   0.97455/

       data gaus(6,1), gaus(6,2), gaus(6,3), gaus(6,4), gaus(6,5),
     >      gaus(6,6)
     >      / 0.03376,   0.16939,   0.38069,   0.61930,   0.83060,
     >        0.96623/

       data gaus(5,1), gaus(5,2), gaus(5,3), gaus(5,4), gaus(5,5)
     >      / 0.04691,   0.23076,   0.50000,   0.76923,   0.95308/

       data gaus(4,1), gaus(4,2), gaus(4,3), gaus(4,4)
     >      / 0.06943,   0.33001,   0.66999,   0.93057/

       data gaus(3,1), gaus(3,2), gaus(3,3)
     >      / 0.11270,   0.50000,   0.88729/

       data gaus(2,1), gaus(2,2)
     >      / 0.21132,   0.78867/

       data gaus(1,1) 
     >      / 0.50000/

 
       data w1(8,1), w1(8,2), w1(8,3), w1(8,4), w1(8,5),
     >      w1(8,6), w1(8,7), w1(8,8)
     >      /.05061,  .11119,  .15685,  .18134,  .18134,
     >       .15685,  .11119,  .05061/

       data w1(7,1), w1(7,2), w1(7,3), w1(7,4), w1(7,5),
     >      w1(7,6), w1(7,7)
     >      /.06474,  .13985,  .19091,  .20897,  .19091,
     >       .13985,  .06474/

       data w1(6,1), w1(6,2), w1(6,3), w1(6,4), w1(6,5),
     >      w1(6,6)
     >      /.08566,  .18038,  .23395,  .23395,  .18038,
     >       .08566/

       data w1(5,1), w1(5,2), w1(5,3), w1(5,4), w1(5,5)
     >      /.11846,  .23931,  .28444,  .23931,  .11846/

       data w1(4,1), w1(4,2), w1(4,3), w1(4,4)
     >      /.17399,  .32607,  .32607,  .17393/

       data w1(3,1),  w1(3,2),  w1(3,3)
     >      /.277777, .444444,  .277777/

       data w1(2,1), w1(2,2)
     >      /.50000,  .50000/

       data w1(1,1)
     >      /1.000/

       

       if (negrid.gt.negmax) then
          write(6,*) "negrid IS BIGGER THAN ",negmax,"   !!!"
          stop
       endif


       tt=2./3.
       ng1= max0( negrid-2, 0 )
       if (negrid.le.2)  ecut= 0.

       e(ng1+1)= ecut +0.58578
       w(ng1+1)= exp(-ecut) * sqrt( e(ng1+1) ) * 0.853553
       
       e(ng1+2)= ecut +3.41421
       w(ng1+2)= 0.146446 * exp(-ecut) * sqrt( e(ng1+2) )

       do 10 i= 1,ng1
          e(i)= ecut*gaus(ng1,i)**tt
          w(i)= tt*exp(-e(i))*w1(ng1,i)*ecut**1.5
 10    continue

       do 20 i=1,negrid
          w(i)=w(i)*0.5/sqrt(pi)
 20    continue


       return
       end subroutine egridset


************************************************************************
       subroutine lgridset
************************************************************************

cgms      Use(Common)

       real wx(nlm),xx(nlm)
       real xgauss(5,10),wgauss(5,10)

       data xgauss(3,1),xgauss(3,2),xgauss(3,3)/.23861,.66120,.93246/
       data wgauss(3,1),wgauss(3,2),wgauss(3,3)/.46791,.36076,.17132/
       data xgauss(4,1),xgauss(4,2),xgauss(4,3),xgauss(4,4)/.18343,
     &                                          .52553,.79666,.96028/
       data wgauss(4,1),wgauss(4,2),wgauss(4,3),wgauss(4,4)/.36268,
     &                                          .31370,.22238,.10122/
       data xgauss(5,1),xgauss(5,2),xgauss(5,3),xgauss(5,4),xgauss(5,5)
     &                             /.14887,.43339,.67940,.86506,.97390/
       data wgauss(5,1),wgauss(5,2),wgauss(5,3),wgauss(5,4),wgauss(5,5)
     &                             /.29552,.26926,.21908,.14945,.06667/



       if (abs(eps).le.1.e-33) then
         nlambda= ng2
       else
         nlambda= ng2 +ntheth+1
       endif

c    zero out weight matrix
       do 5 i=1,nlambda
       do 5 j=-ntgridl,ntgrid
          wl(j,i)=0.
 5     continue

c    set up u grid and weights
       do 10 i=1,ngauss
         xx(i)=.5*(1.+xgauss(ngauss,ngauss+1-i))
         xx(ngauss+i)=.5*(1.-xgauss(ngauss,i))
         wx(i)=.5*wgauss(ngauss,ngauss+1-i)
         wx(ngauss+i)=.5*wgauss(ngauss,i)
10     continue

c    change u grid and weights to al grid and weights
       do 20 i=1,ng2
         al(i)=(1.-xx(i)**2)/bmax
       do 20 j= -nthethl,ntheth
         wl(j,i)=wx(i)*2.*sqrt( (bmag(j)/bmax)
     >                         *( (1./bmax-al(i))/(1./bmag(j)-al(i))) )
20     continue

c     trapped particle region computed here
       if (abs(eps).le.1.e-33) goto 100

       do 30 i=1,ntheth+1
         al(i+ng2)=1./bmag(ntheth+1-i)
30     continue

c   set up weights
       do 40 i=1,ntheth
       do 40 j=-ntheth+i,ntheth-i
         ww=sqrt(abs(1.-al(ng2+i  )*bmag(j)))-
     &      sqrt(abs(1.-al(ng2+i+1)*bmag(j)))
         wl(j,ng2+i  )=wl(j,ng2+i  )+ww
         wl(j,ng2+i+1)=wl(j,ng2+i+1)+ww
40     continue

c     copy these coeficients to the various periods

100    continue

       do 200 iperiod=nperiod11,nperiod2
       do 200 i=1,nlambda
       do 200 j=-nthethl,ntheth
         wl(iperiod*ntheta+j,i)=wl(j,i)
crew
         wl(iperiod*ntheta-j,i)=wl(j,i)
200    continue

       return
       end subroutine lgridset


************************************************************************
       subroutine eintegrate
************************************************************************

cgms      Use(Common)



       do  91 is=1,nspec
         if(ncspec(is).eq.0) goto 91
         
         do 90 j=-ntgridl,ntgrid
90         geint(j,is)=0.

91     continue

       do 101 is=1,nspec
         if(ncspec(is).eq.0) goto 101
         
         do 100 i=1,negrid
         do 100 j=-ntgridl,ntgrid
           geint(j,is)=geint(j,is)+w(i)*gint(j,i,is)
100      continue  

101    continue

       return
       end subroutine eintegrate



************************************************************************
       subroutine lintegrate
************************************************************************

cgms      Use(Common)


       do 11 is=1,nspec
         if(ncspec(is).eq.0) goto 11

         do 10 i=1,negrid
         do 10 j=-ntgridl,ntgrid
            gint(j,i,is)=0.
 10      continue

11     continue


       do 21 is=1,nspec
         if(ncspec(is).eq.0) goto 21
         
         do 20 ie=1,negrid
         do 20 i=1,nlambda
         do 20 j=-ntgridl,ntgrid
           gint(j,ie,is)=gint(j,ie,is)+wl(j,i)*
     &                   (g1(j,i,ie,1,is)+g1(j,i,ie,2,is))
20       continue

21     continue

       return
       end subroutine lintegrate


************************************************************************
       subroutine ginitial
************************************************************************

cgms      Use(Common)
crew
cgms      Use(Geo)


c   set up initial conditions for g

       if (iphi.eq.0) then
            xk= 0.5065 
          do 10 j=-ntgridl,ntgrid
             phi(j)=exp(-( (theta(j)-theta0) /width0)**2)*cmplx(1.,1.)
cc             xlamd= 2.*atan(xk*delthet(j)/2.)
cc             phi(j)= cos(xlamd*j)

crew6.95 correction factor (phi_cor) to make actual phi initial gausian
crew            phi_cor=0.
crew            do is=1,nspec
crew             phi_cor=phi_cor+an(is)*z(is)**2*poldrift(j,is)/temp(is)
crew            enddo
crew            phi(j)=phi(j)*phi_cor
             phi_cor=1.+shat**2*theta(j)**2
             phi(j)=phi(j)*phi_cor
          if(ioddpar.eq.1) phi(j)=phi(j)*theta(j)

 10       continue

       endif

       if (iphi.eq.1) then
          do 20 j=-ntgridl,ntgrid
             phi(j)=exp(-(theta(j)/width0)**2)
             
cc             if ( abs(theta(j)-theta0).lt.pi ) then
cc                phi(j)= (theta(j)-theta0)**2 -pi*pi
cc             else
cc                phi(j)= 0.0
cc             endif

cc             phi(j)= cos(theta(j))
 20       continue
       endif

       if (ipar.eq.0) phi(-ntgrid)=0.
       phi( ntgrid)=0.


       if((iaky.eq.1).and.(igam.eq.1)) then
       
crew
       if(icontinue.eq.0) then


         do 101 is=1,nspec
           if(ncspec(is).eq.0) goto 101

           do 100 ie=1,negrid
           do 100 i=1,nlambda
           do 100 j=-ntgridl,ntgrid
              g(j,i,ie,1,is)=-phi(j)*z(is)
              g(j,i,ie,2,is)=g(j,i,ie,1,is)
cc              g(j,i,ie,1,is)= sqrt( abs(1. -al(i)*bmag(j)) )
cc              g(j,i,ie,2,is)=-g(j,i,ie,1,is)
cc              g(j,i,ie,1,is)= 1. -3.*(1.-al(i)*bmag(j))
cc              g(j,i,ie,2,is)=g(j,i,ie,1,is)
cc              g(j,i,ie,1,is)= 1.0
cc              g(j,i,ie,2,is)= 1.0
cc              g(j,i,ie,1,is)= 0.
cc              g(j,i,ie,2,is)= 0.
 100       continue

 101     continue
crew   end icontinue
       endif

       endif


       an(5)= zeff5/(z(5)**2)
cgms       an(2)= ( zeff -1. -z(5)*(z(5)-1.)*an(5) )/( (z(2)-1.)*z(2) )
       an(2)= an2
       
       anlt= 1. -an(2)*z(2) -an(5)*z(5)
       an(4)= dbeam
cgms       an(1)= anlt -an(4)
       an(1)=an1
       an(3)= 0.
       xteti=teti

       if(abs(z(3)+1.).le.1.e-33) an(3)=1.
       if(abs(z(3)+1.).le.1.e-33) xteti=0.
         

       do 201 is=1,nspec
         if(ncspec(is).eq.0) goto 201

         do 200 ie=1,negrid
         do 200 i=1,nlambda
         do 200 j=-ntgridl,ntgrid
            sqab=max(0.,(1.-al(i)*bmag(j)))

            g1(j,i,ie,1,is)=sqab
            g1(j,i,ie,2,is)=sqab

            sq(j,i,1)=sqrt(sqab)
            sq(j,i,2)=-sqrt(sqab)
200      continue

201    continue


       call lintegrate

       do 251 is=1,nspec
          if(ncspec(is).eq.0) goto 251

          do 250 ie=1,negrid
          do 250 j=-ntgridl,ntgrid
             g3int(j,ie,is)= gint(j,ie,is)
 250      continue

 251   continue

             
       do 301 is=1,nspec
         if(ncspec(is).eq.0) goto 301

         do 300 ie=1,negrid
         do 300 i=1,nlambda
         do 300 j=-ntgridl,ntgrid
           g1(j,i,ie,1,is)=1.
           g1(j,i,ie,2,is)=1.
300      continue

301    continue

       call integrate


       do 311 is=1,nspec
         if(ncspec(is).eq.0) goto 311

         do 310 j=-ntgridl,ntgrid
            aintnorm(j,is)=1./real(geint(j,is))
 310     continue

311    continue


       do 500 isign=1,2
       do 500 ie=1,negrid
       do 500 il=1,nlambda
       do 500 j=-ntgridl,ntgrid
         vpa(j,il,ie,isign)=sq(j,il,isign)*sqrt(e(ie))
         vper2(j,il,ie)=e(ie)*al(il)*bmag(j)
500    continue


       do 510 ie=1,negrid
       do 510 il=1,nlambda
       do 510 j=-ntgridl,ntgrid-1
          bmid=0.5*(bmag(j)+bmag(j+1))
          xtemp= 1.-al(il)*bmid
          xtemp= 0.5*( xtemp + abs(xtemp) )
          vpac(j,il,ie,1)= sqrt( xtemp*e(ie) )
          vpac(j,il,ie,2)=-vpac(j,il,ie,1)
 510   continue
           

c     set up vparallel coef matrix 
       do 521 is=1,nspec
         if(ncspec(is).eq.0) goto 521

         zstm=z(is)/sqrt(temp(is)*amass(is))
         
         do 520 ie=1,negrid
         do 520 il=1,nlambda
         do 520 j=-ntgridl,ntgrid-1
            vpar(j,il,ie,1,is)=zstm*pp*vpac(j,il,ie,1)*delt
     >                         /delthet(j)/(1.)
crrd      shifted circle factor for vpar eliminated in above
crew_______________geo__________________________________________________
     
          if (igeo.eq.1)
     >    vpar(j,il,ie,1,is)=zstm*ppgeo(j)*vpac(j,il,ie,1)*
     &                        delt/delthet(j)
     
crew_______________geo__________________________________________________
            vpar(j,il,ie,2,is)=-vpar(j,il,ie,1,is)
 520     continue
            
 521   continue


       return
       end subroutine ginitial





************************************************************************
       subroutine integrate
************************************************************************

cgms      Use(Common)


       call lintegrate
       call eintegrate

       return
       end subroutine integrate


************************************************************************
       subroutine gridset
************************************************************************

cgms      Use(Common)


       ntheth=ntheta/2
crew add ntheth
       nthethl=ntheth
       if (ipar.eq.1) nthethl=0
       
       ntgrid=ntheta/2+(nperiod-1)*ntheta
crew add ntgridl
       ntgridl=ntgrid
       if (ipar.eq.1) ntgridl=0
       
       ng2=2*ngauss
       
       nperiod2= nperiod-1
crew replaced nperiod1 with nperiod11 and nperiod10
       nperiod11=-nperiod+1
       nperiod10=-nperiod+1
       if (ipar.eq.1) then
        nperiod11=1
        nperiod10=0
       endif
       

       do 100 i=-ntgridl,ntgrid
100      theta(i)=real(i)*2.*pi/real(ntheta)
         
       do 110 i=-ntgridl,ntgrid-1
110      delthet(i)=theta(i+1)-theta(i)
         
       delthet(ntgrid)=delthet(ntgrid-1)
       if (ipar.eq.0) delthet(-ntgrid)=delthet(-ntgrid+1)
       

       call magcoord
       call egridset
       call lgridset

       return
       end subroutine gridset


************************************************************************
       subroutine magcoord
************************************************************************

cgms      Use(Common)
crew________geo________________________________________________________
cgms      Use(Geo)

       integer ibmin

c     set up magnetic field
       if(igeo.eq.0)then
         do 7 i=-ntgridl,ntgrid
           bmag(i)=1./(1.+eps*cos(theta(i)))
7        continue
         bmin=1./(1.+eps)
         bmax=1./(1.-eps)
       endif
       ntheth=ntheta/2
       
crew________geo________________________________________________________
       if (igeo.eq.1) then
       
           do i=-ntgridl,ntgrid
	     bmag(i)=b_geo(i)
	   enddo
c check that bmin is at theta=0
           bmin = b_geo(0)
           ibmin=0
           do i=0,ntheth
              if(b_geo(i).lt.bmin)then
                bmin=b_geo(i)
                ibmin=i
              endif
           enddo
	   if(ibmin.ne.0)then
              write(*,*)" fixing double bmag well",ibmin
c fix double well problem by making minimum at theta=0
c this only impacts the pitch angle weights for igeo=1
              do j=0,ibmin-1
                bmag(j) = bmin*(1.0-0.001*REAL(ibmin-j)/REAL(ibmin))
                if(nperiod.gt.1)then
                  do k=2,nperiod
                    jk=(k-1)*ntheta
                    bmag(jk+j)=bmag(j)
                    bmag(jk-j)=bmag(j)
                  enddo
                endif
              enddo
              if(ntgridl.ne.0)then
                do i=-ntgridl,-1
                  bmag(i)=bmag(-i)
                enddo
              endif
            endif
      	  bmin=bmag(0)
	  bmax=bmag(ntheth)

       endif
crew________geo________________________________________________________
       

       return
       end subroutine magcoord
       
************************************************************************
      subroutine matrix
************************************************************************

cgms      Use(Common)


       do 11 is=1,nspec
         if(ncspec(is).eq.0) goto 11

         do 10 isign=1,2
         do 10 ie=1,negrid
         do 10 il=1,nlambda
         do 10 j=-ntgridl,ntgrid
           g0(j,il,ie,isign,is)= 0. 
10       continue

11     continue


       
       do 100 mj= 1,2*(ntgrid+ntgridl+1)
         
          do 20 j= -ntgridl,ntgrid
             phig (j)= 0.
             aparg(j)= 0.
             

             phi (j)= 0.
             apar(j)= 0.
             
20        continue

            mjjudge= (mj-mj/2*2) + 2*(mj+1-(mj+1)/2*2)

          if (mjjudge.eq.1) then
            j= mj/2 -ntgridl
            phi (j)= cmplx(1.,0.)
          endif

          if (mjjudge.eq.2) then
            j= (mj-1)/2 -ntgridl
            apar(j)= cmplx(1.,0.)
          endif


          if(i_bpar.eq.0) call timeadv
          if(i_bpar.eq.-1) call timeadv_bpar
          call getan


             miodd = -1
             mieven= 0
          do 50 j= -ntgridl,ntgrid
             miodd = miodd  +2
             mieven= mieven +2
             
             am(miodd ,mj)= antot (j)  
             am(mieven,mj)= antota(j)  
50        continue

100    continue

c       write(6,*) "called matrix"


       return
       end subroutine matrix
************************************************************************
      subroutine matrix_bpar
************************************************************************
c      1/21/98 modification of matrix to include bpar
cgms      Use(Common)



       do 11 is=1,nspec
         if(ncspec(is).eq.0) goto 11

         do 10 isign=1,2
         do 10 ie=1,negrid
         do 10 il=1,nlambda
         do 10 j=-ntgridl,ntgrid
           g0(j,il,ie,isign,is)= 0. 
10       continue

11     continue

       mspan=(ntgrid+ntgridl+1)
       mspanp1=mspan+1
       mspan2=2*mspan
       mspan2p1=mspan2+1
       mspan3=3*mspan
  
c       write(6,*) 'mspan = ',mspan,ntgrid,ntgridl
       do mj=1,3*mspan
        do mi=1,3*mspan
         am(mi,mj)=0.
        enddo
       enddo
       
       nspan=2
       if(i_bpar.eq.1) nspan=3
       do 100 mj= 1,nspan*mspan
         
          do 20 j= -ntgridl,ntgrid
             phig (j)= 0.
             aparg(j)= 0.
             bparg(j)= 0.

             phi (j)= 0.
             apar(j)= 0.
             bpar(j)= 0.
20        continue


          if(mj.ge.1.and.mj.le.mspan)         mjjudge=1
          if(mj.ge.mspanp1.and.mj.le.mspan2)  mjjudge=2
          if(mj.ge.mspan2p1.and.mj.le.mspan3) mjjudge=3
                    
          if (mjjudge.eq.1) then
            j= mj-(mspan-ntgrid-1)-1
            phi (j)= cmplx(1.,0.)
          endif

          if (mjjudge.eq.2) then
            j= mj-mspan-(mspan-ntgrid-1)-1
            apar(j)= cmplx(1.,0.)
          endif
          
          if (mjjudge.eq.3) then
            j= mj-2*mspan-(mspan-ntgrid-1)-1
            bpar(j)= cmplx(1.,0.)
          endif


          call timeadv_bpar
          call getan

 
             mi1=0
             mi2=mspan
          do 50 jj= -ntgridl,ntgrid
             mi1 = mi1+1
             mi2 = mi2+1
             am(mi1,mj)= antot (jj)  
             am(mi2,mj)= antota(jj)
50        continue

          if(i_bpar.eq.1) then
             mi3=mspan2
          do 51 jj= -ntgridl,ntgrid
             mi3 = mi3+1
             am(mi3,mj)= antotb(jj)
51        continue
          endif

100    continue


c       write(6,*) "called matrix_bpar"

       return
       end subroutine matrix_bpar

************************************************************************
      subroutine matrix_bpar2
************************************************************************
c      8/29/98 modification of matrix to include bpar
c      this used 1*j-> phi, 2*j-> apar, 3*j-> bpar
cgms      Use(Common)



       do 11 is=1,nspec
         if(ncspec(is).eq.0) goto 11

         do 10 isign=1,2
         do 10 ie=1,negrid
         do 10 il=1,nlambda
         do 10 j=-ntgridl,ntgrid
           g0(j,il,ie,isign,is)= 0.
10       continue

11     continue

       mspan=(ntgrid+1)


       do 100  mj=1,3*(ntgrid+1)

          do 20 jj= 0,ntgrid
             phig (jj)= 0.
             aparg(jj)= 0.
             bparg(jj)= 0.

             phi (jj)= 0.
             apar(jj)= 0.
             bpar(jj)= 0.
20        continue


        mprejudge=(mj-mj/3*3)+3*(mj+1-(mj+1)/3*3)+3*(mj+2-(mj+2)/3*3)
c          7,5,9
          if(mprejudge.eq.7) then
            mjjudge=1
            j=(mj+3-mjjudge)/3-1 
            phi(j)= cmplx(1.,0.)
          endif
          if(mprejudge.eq.5) then
            mjjudge=2
            j=(mj+3-mjjudge)/3-1
            apar(j)= cmplx(1.,0.)
          endif
          if(mprejudge.eq.9) then
            mjjudge=3
            j=(mj+3-mjjudge)/3-1
            bpar(j)= cmplx(1.,0.)
          endif
c          [mj,mjjudge,j]  



          call timeadv_bpar
          call getan

          mione=-2
          mitwo=-1
          mithree=0
          do 50 jj= 0,ntgrid
           mione=mione+3
           mitwo=mitwo+3
           mithree=mithree+3
             am(mione,mj)= antot (jj)
             am(mitwo,mj)= antota(jj)
             am(mithree,mj)= antotb(jj)
50        continue


100    continue


c       write(6,*) "called matrix_bpar2"

       return
       end subroutine matrix_bpar2  
************************************************************************
       subroutine plasmacoef
************************************************************************

cgms      Use(Common)
crew_______________geo__________________________________________________
cgms      Use(Geo)

       real gam(-ntml:ntm,nspec)


       do 501 is=1,nspec
         if(ncspec(is).eq.0) goto 501

         tz=temp(is)/z(is)
         zstm=z(is)/sqrt(temp(is)*amass(is))
         smz=sqrt(amass(is)/abs(z(is)))*temp(is)
         if(z(is).lt.0.) smz=0.


         do 100 ie=1,negrid
         do 100 il=1,nlambda
         do 100 j=-ntgridl,ntgrid
crew           zzz=smz*aky*sqrt( e(ie)*al(il)*bmag(j)
crew     >                      *(1.+(shat*(theta(j)-theta0))**2) )
crew CORRECTION on smz factor
         gyro_fix=1.
         if (igyro_fix.eq.1) gyro_fix=bmag(j)
           zzz=sqrt(amass(is)*temp(is))/abs(z(is))
     >          *(aky/gyro_fix)*sqrt( e(ie)*al(il)*bmag(j)
     >         *(1.+(shat*(theta(j)-theta0) 
     >         - shift*sin(theta(j)))**2))
crrd       non circular shift correction added to above
     
crew_______________geo__________________________________________________

       b_geo_here=b_geo(j)
       if (abs(igeot).eq.1.or.abs(igeot).eq.4) b_geo_here=1.
       if (igeo.eq.1) 
     >     zzz=sqrt(amass(is)*temp(is))/abs(z(is))
     >      *aky*qrat_geo(j)/b_geo_here*sqrt( e(ie)*al(il)*bmag(j)
     >                         *(1.+kxoky_geo(j)**2) )

crew_______________geo__________________________________________________
           aj0(j,il,ie,is)=aj0f(zzz)
           azj1(j,il,ie,is)=azj1f(zzz)

         if(igyro_e.eq.1.and.is.eq.3) then
crew  06.08.01 rho_e -> 0  no gyro aveage in electrons
          aj0(j,il,ie,is) = 1.
          azj1(j,il,ie,is) = 0.
         endif

100      continue


c   set up wstar
         do 300 ie=1,negrid
crew    signomega=1 ion modes + flow to -x;   =-1 ion modes + flow to +x

crew________geo________________________________________________________
crew  note wstar does not change!!!!!!!
crew________geo________________________________________________________

        wstar(ie,is)=signomega*delt*(fprim(is)+tprim(is)*(e(ie)-1.5))

 300     continue

c   set up omega drift matrix
         if(icv.ne.0) then

           do 350 ie=1,negrid
           do 350 il=1,nlambda
           do 350 j=-ntgridl,ntgrid
             shfac=1.-shift*cos(theta(j))
             shfac=1.
crrd      non circular factor set to unity

crew    signomega
             wdrift(j,il,ie,is)=
     >           signomega*delt*epsl*e(ie)*(1.-al(il)*bmag(j)*0.5)*
     &                          (-xwell+cos(theta(j))/shfac+
     &                          (shat*(theta(j)-theta0) - shift*
     &                          sin(theta(j)) )*sin(theta(j))*shfac)
crew error fixed above shatxshift 4/4/96
crew xwell added
     
crew________geo________________________________________________________
         
	    if(igeo.eq.1)
     >       wdrift(j,il,ie,is)=
     >       signomega*delt*epsl_geo(j)*(e(ie)*(1.-al(il)*bmag(j)*0.5)*
     &              (costheta_geo(j)+kxoky_geo(j)*sintheta_geo(j))
     &           +e(ie)*(1.-al(il)*bmag(j))*costheta_p_geo(j))
     
crew________geo________________________________________________________

350        continue

         else

           do 360 ie=1,negrid
           do 360 il=1,nlambda
           do 360 j=-ntgridl,ntgrid
crew   signomega
             wdrift(j,il,ie,is)=
     >          signomega*delt*epsl*e(ie)*(1.-al(il)*bmag(j)*0.5)
360        continue

         endif

         do 370 ie=1,negrid
         do 370 il=1,nlambda
         do 370 j=-ntgridl,ntgrid-1
            wdrift(j,il,ie,is)=0.5*( wdrift(j  ,il,ie,is)
     >                              +wdrift(j+1,il,ie,is) )
370      continue


         fac=1.-fexp(is)

c    set up time advancement matrices
         facz=fac*tz
         fexpz=fexp(is)*tz

         do 400 ie=1,negrid
         do 400 il=1,nlambda
         do 400 j=-ntgridl,ntgrid-1
           vptemp= vpar(j,il,ie,1,is)
           wdtemp= wdrift(j,il,ie,is)
crew added omegaprev terms
          ainv(j,il,ie,is)=1./(1.+bakdif(is)+zi*facz *wdtemp
     >                           -zi*facz/tz*delt*omegaprev
     &                           +2.*facz*vptemp)
             r(j,il,ie,is)=   (1.-bakdif(is)+zi*facz *wdtemp
     >                           -zi*facz/tz*delt*omegaprev
     &                           -2.*facz*vptemp)*ainv(j,il,ie,is)
             a(j,il,ie,is)=   (1.+bakdif(is)-zi*fexpz*wdtemp
     >                           +zi*fexpz/tz*delt*omegaprev
     &                           -2.*fexpz*vptemp)
             b(j,il,ie,is)=   (1.-bakdif(is)-zi*fexpz*wdtemp
     >                           +zi*fexpz/tz*delt*omegaprev
     &                           +2.*fexpz*vptemp)
400      continue



c  zero out advancement matrices for unphysical regions of phase space
         if (abs(eps).le.1.e-33) goto 480

         do 430 iperiod=nperiod10,nperiod2
            j0= iperiod*ntheta
            jplus2=  j0 +ntheth
            jminus1= j0 -ntheth
         do 430 ie=1,negrid
         do 430 il=2,ntheth+1
           iln=il +ng2
           jplus= jplus2 -il +2
           do 410 j=jplus,jplus2
             ainv(j,iln,ie,is)=0.
                r(j,iln,ie,is)=0.
 410        continue

           ainv(jplus-1,iln,ie,is)=1.
              r(jplus-1,iln,ie,is)=0.

           jminus=jminus1 +il -2

           do 420 j=jminus1,jminus
             ainv(j,iln,ie,is)=0.
                r(j,iln,ie,is)=0.
 420       continue

           ainv(jminus,iln,ie,is)=1.

 430      continue


         do 440 iperiod=nperiod10,nperiod2
           j0=iperiod*ntheta
         do 440 ie=1,negrid
            wdtemp= wdrift(j0,nlambda,ie,is)
           ainv(j0,nlambda,ie,is)=1./(1.+zi*facz*wdtemp)
              a(j0,nlambda,ie,is)=1.-zi*fexpz*wdtemp
 440       continue


480      continue


501   continue



       do 601 is=1,nspec
         if(ncspec(is).eq.0) goto 601

         alpp=4./(3.*sqrt(pi))
         
         if(is.eq.3) then

           do 520 ie=1,negrid
              vnew(ie,is)=vnew0(is)/e(ie)**1.5*( zeff+alpp*e(ie)/
     &                                     sqrt(1.+(alpp*e(ie))**2) )
520        continue

         else
            
           do 540 ie=1,negrid
             vnew(ie,is)=vnew0(is)*alpp*e(ie)/
     &                   (sqrt(1.+(alpp*e(ie))**2)*e(ie)**1.5)
 540       continue

         endif

 601   continue


c  compute polarization factor needed for quasinuetrality
c 1/21/98 gamtotb i_bpar=1 added. 
c  note: gamtota does not enter equations?

       do 620 j= -ntgridl,ntgrid
          gamtot (j)= 0.0
          gamtotpb(j)=0.0
          gamtota(j)= 0.0
          gamtotb(j)= 0.0
          gamtotbb(j)=0.0
 620    continue
 
c get gamtot


       do 651 is=1,nspec
         if(ncspec(is).eq.0) goto 651

         do 650 ie=1,negrid
         do 650 il=1,nlambda
         do 650 j=-ntgridl,ntgrid
           g1(j,il,ie,1,is)=(1.-aj0(j,il,ie,is)**2)
           g1(j,il,ie,2,is)=g1(j,il,ie,1,is  )
650      continue

651    continue

       call integrate


       do 661 is=1,nspec
         if(ncspec(is).eq.0) goto 661

         do 660 j=-ntgridl,ntgrid
           gam(j,is)=real(geint(j,is))
crew added
           poldrift(j,is)=gam(j,is)
660      continue

661    continue


       do 671 is= 1,nspec
          if(ncspec(is).eq.0) goto 671

          do 670 j=-ntgridl,ntgrid
             gamtot(j)= gamtot(j) +z(is)*z(is)*an(is)*gam(j,is)/temp(is)
 670      continue

 671   continue
 
c rew 2/3/97 add charge separation lamda_Debye term assuming k_par<<k_perp
c   debyelorhos=debye_length_e / rhos = ((Me/Mi)*omega_pe/omega_ce)**.5 
c common parameter is de=1.+(omega_pe/omega_ce)**2
c                       =1.+(debyelorhos*Mi/Me**.5)**2
c aky=ky*rho_i  with rho_i=(2.*Ti/Mi)**.5*(cMi/eB)
c note the debye term adds to pol. drift for both species so for ETG modes
c at high k_perp when ions are adiabatic
c    Te/Ti+ k_perp**2*rho_e**2 -> Te/Ti+ de*k_perp**2*rho_e**2
c   see Horton,Hong, and Tang  PF 31 (1988) 2971
c   typical DIII-D deuterium 2.1T 6.e13  debyelorhos=1.39642E-02
c    so de=1.+(60.x1.4e-2)**2 which is order 2 instead of 1.
c    thus debyelorhos is small for ITG but significant for ETG
c 
       write(*,*)"debyelorhos=",debyelorhos,"aky=",aky
        if(igeo.eq.0) then
          do  j=-ntgridl,ntgrid
           gamtot(j)=gamtot(j)+(aky**2)*(temp3/2.)*(debyelorhos**2)*
     >    (1.+(shat*(theta(j)-theta0) - shift*sin(theta(j)))**2)
          enddo
        endif
        if( igeo.eq.1) then
          do  j=-ntgridl,ntgrid
           gamtot(j)=gamtot(j)+(aky**2)*(temp3/2.)*(debyelorhos**2)*
     >      ((qrat_geo(j)/b_geo(j))**2)*(1.+kxoky_geo(j)**2)
          enddo
        endif
  
c get gamtotpb


       do 656 is=1,nspec
         if(ncspec(is).eq.0) goto 656

         do 655 ie=1,negrid
         do 655 il=1,nlambda
         do 655 j=-ntgridl,ntgrid
           g1(j,il,ie,1,is)=
     >        (-aj0(j,il,ie,is)*azj1(j,il,ie,is)*vper2(j,il,ie))
           g1(j,il,ie,2,is)=g1(j,il,ie,1,is  )
655      continue

656    continue

       call integrate


       do 666 is=1,nspec
         if(ncspec(is).eq.0) goto 666

         do 665 j=-ntgridl,ntgrid
           gam(j,is)=real(geint(j,is))
665      continue

666    continue


       do 676 is= 1,nspec
          if(ncspec(is).eq.0) goto 676

          do 675 j=-ntgridl,ntgrid
             gamtotpb(j)= gamtotpb(j) +xy_bpar*
     > z(is)*z(is)*an(is)*gam(j,is)/temp(is)*2.*temp(is)/(z(is))
 675      continue
c 9/3/98 factor 2. was missed
c 9/18/98 abs(z)->(z)

 676   continue 

c get gamtota (not used ??????)

       do 681 is=1,nspec
         if(ncspec(is).eq.0) goto 681

         do 680 ie=1,negrid
         do 680 il=1,nlambda
         do 680 j=-ntgridl,ntgrid
           g1(j,il,ie,1,is)=(aj0(j,il,ie,is)**2)*(vpa(j,il,ie,1)**2)
           g1(j,il,ie,2,is)=g1(j,il,ie,1,is  )
680      continue

681    continue

       call integrate


       do 691 is=1,nspec
         if(ncspec(is).eq.0) goto 691

         do 690 j=-ntgridl,ntgrid
           gam(j,is)=real(geint(j,is))
690      continue

691    continue


       do 701 is= 1,nspec
          if(ncspec(is).eq.0) goto 701

          do 700 j=-ntgridl,ntgrid
             gamtota(j)= gamtota(j) +2.*beta*an(is)*z(is)*z(is)*
     >                                  gam(j,is)/amass(is)
 700      continue

 701   continue
 
 
c get gamtotb

       do 851 is=1,nspec
         if(ncspec(is).eq.0) goto 851

         do 850 ie=1,negrid
         do 850 il=1,nlambda
         do 850 j=-ntgridl,ntgrid
c note (-) sign like (1-aj0**2)
           g1(j,il,ie,1,is)=-aj0(j,il,ie,is)*azj1(j,il,ie,is)
     >             *vper2(j,il,ie)
           g1(j,il,ie,2,is)=g1(j,il,ie,1,is  )
850      continue

851    continue

       call integrate


       do 861 is=1,nspec
         if(ncspec(is).eq.0) goto 861

         do 860 j=-ntgridl,ntgrid
           gam(j,is)=real(geint(j,is))

860      continue

861    continue


       do 871 is= 1,nspec
          if(ncspec(is).eq.0) goto 871
 
          sign_lhs_rhs=1.
          do 870 j=-ntgridl,ntgrid
             gamtotb(j)=gamtotb(j) + sign_lhs_rhs*x_bpar*
     >(-1.)*beta*z(is)*an(is)*gam(j,is)
 870      continue
c 9/18/98 abs(z) removed

 871   continue

c get gamtotbb

       do 856 is=1,nspec
         if(ncspec(is).eq.0) goto 856

         do 855 ie=1,negrid
         do 855 il=1,nlambda
         do 855 j=-ntgridl,ntgrid
c note (-) sign like (1-aj0**2)
           g1(j,il,ie,1,is)=-azj1(j,il,ie,is)**2
     >             *vper2(j,il,ie)**2
           g1(j,il,ie,2,is)=g1(j,il,ie,1,is  )
855      continue

856    continue

       call integrate


       do 866 is=1,nspec
         if(ncspec(is).eq.0) goto 866

         do 865 j=-ntgridl,ntgrid
           gam(j,is)=real(geint(j,is))

865      continue

866    continue


       do 876 is= 1,nspec
          if(ncspec(is).eq.0) goto 876

          sign_lhs_rhs=1.
c sign for moving from lhs to rhs same as for -j0**2 and already accounted
          do 875 j=-ntgridl,ntgrid
             gamtotbb(j)=gamtotbb(j) + sign_lhs_rhs*x_bpar*xy_bpar*
     >(-1.)*beta*z(is)*an(is)*2.*temp(is)/(z(is))*gam(j,is)
 875      continue
c 9/3/98 factor 2. was missed
c 9/18/98  z**2/abs(z)-> z/z
 876   continue   
  


       return
       end subroutine plasmacoef


************************************************************************
       real function aj0f(xx)
************************************************************************

	   real x,y,f0,t0
           x=abs(xx)

	   if ( x.gt.3. ) goto 100

	   y=(x/3.)**2
	   aj0f=1.-2.2499997*y+1.2656208*y*y-.3163866*y**3+.0444479*y**4
     &  -.0039444*y**5 + .0002100*y**6

	   return

100	  y=3./x
	  f0=.79788456-.0000077*y-.00552740*y*y-.00009512*y**3
     & +y**4*(.00137237-.00072805*y+.00014476*y*y)
	  t0=x-.78539816-.04166397*y-.00003954*y*y+.00262573*y**3
     & +y**4*(-.00054125-.00029333*y+.00013558*y*y)
	  aj0f=f0*cos(t0)/sqrt(x)

	  return
	  end function aj0f
************************************************************************
       real function azj1f(xx)
************************************************************************
crew 1/21/98
c    this function is needed for i_bpar=1
c    azj1(z)=z*J1(z) 
c 9/1/98 azj1(z)=J1(z)/z
c  we use identity J1(z)=-d/dz J0(z)

       real x,f0,xp,f0p
       x=abs(xx)
       
       f0=aj0f(x)
       xp=x+.01*x
       f0p=aj0f(xp)
       
c 9/1/98       azj1f=-x*(f0p-f0)/(.01*x)
       azj1f=-(f0p-f0)/(.01*x)/x
          
	  
	  return
	  end function azj1f


************************************************************************
        subroutine replace
************************************************************************

cgms      Use(Common)


       do 101 is=1,nspec
         if(ncspec(is).eq.0) goto 101

         do 100 isign=1,2
         do 100 ie=1,negrid
         do 100 il=1,nlambda
         do 100 j=-ntgridl,ntgrid
           g(j,il,ie,isign,is)=gnew(j,il,ie,isign,is)
100      continue

101    continue

       return
       end subroutine replace


************************************************************************
       subroutine kperp
************************************************************************

cgms      Use(Common)
crew________geo________________________________________________________
cgms      Use(Geo)



       real akperp,phinorm
       akperp=0.
       phinorm=0.

       do 100 j=-ntgridl,ntgrid
         akperp=akperp+cabs((phig(j)+aparg(j)+bparg(j))**2)*
     >          (1.+shat**2*(theta(j)-theta0)**2)
     
crew________geo________________________________________________________
      if (igeo.eq.1)
     >         akperp=akperp+cabs((phig(j)+aparg(j)+bparg(j))**2)*
     > qrat_geo(j)**2/b_geo(j)**2*(1.+kxoky_geo(j)**2)
     
crew________geo________________________________________________________
         phinorm=phinorm+cabs(phig(j)+aparg(j)+bparg(j))**2
100    continue

       akperp=sqrt(akperp/phinorm)*aky

       return
       end subroutine kperp

************************************************************************
       subroutine flux
************************************************************************

cgms      Use(Common)


c   compute particle flux
       do 101 is=1,nspec
         if(ncspec(is).eq.0) goto 101

         do 100 isign=1,2
         do 100 ie=1,negrid
         do 100 il=1,nlambda
         do 100 j=-ntgridl,ntgrid
           g1(j,il,ie,isign,is)=g(j,il,ie,isign,is)*aj0(j,il,ie,is)
100      continue

101    continue

       call integrate

c   compute normalizing phi**2 integral
       anorm=0.
       do 130 j=-ntgridl,ntgrid
          anorm=anorm+cabs(phig(j)+aparg(j))**2
 130   continue

c   compute particle flux
       do 201 is=1,nspec
         if(ncspec(is).eq.0) goto 201

         pflux(is)=0.
         do 200 j=-ntgridl,ntgrid
          pflux(is)=pflux(is)-aimag(geint(j,is)*conjg(phig(j)+aparg(j)))
200      continue

201    continue

       do 250 is=1,nspec
         if(ncspec(is).eq.0) goto 250

         pflux(is)=pflux(is)/anorm
 250   continue

c   compute q flux
       do 301 is=1,nspec
         if(ncspec(is).eq.0) goto 301

         do 300 isign=1,2
         do 300 ie=1,negrid
         do 300 il=1,nlambda
         do 300 j=-ntgridl,ntgrid
           g1(j,il,ie,isign,is)=g1(j,il,ie,isign,is)*(e(ie)-1.5)
300      continue

301    continue

       call integrate

c   compute q flux
       do 401 is=1,nspec
         if(ncspec(is).eq.0) goto 401

         qheat(is)=0.
         do 400 j=-ntgridl,ntgrid
           qheat(is)=qheat(is)-aimag(geint(j,is)*conjg(phig(j)))
400      continue

401    continue

       do 450 is=1,nspec
         if(ncspec(is).eq.0) goto 450

          qheat(is)=qheat(is)/anorm
450    continue

       return
       end subroutine flux
************************************************************************
       subroutine solfp1
************************************************************************

cgms      Use(Common)

       integer ivb(-ntml:ntm)


       if(vnew(2,3).ge.-1.e-33) then

         do 31 is=1,nspec
           if(ncspec(is).eq.0) go to 31

           do 10 isign=1,2
           do 10 ie=1,negrid
           do 10 j=1,nlambda
           do 10 i=-ntgridl,ntgrid
             g1(i,j,ie,isign,is)=0.
 10        continue

           if(abs(ncspec(is)*vnew(2,is)).le.1.e-33) goto 31
           
           do 20 isign=1,2
           do 20 ie=1,negrid
           do 20 j=1,nlambda
           do 20 i=-ntgridl,ntgrid
             g1(i,j,ie,isign,is)=gnew(i,j,ie,isign,is)*sq(i,j,isign)
20         continue

 31      continue

         call lintegrate

         do 51 is=1,nspec
            if(ncspec(is).eq.0) goto 51

            do 50 ie=1,negrid
            do 50 j=-ntgridl,ntgrid
               g1int(j,ie,is)= gint(j,ie,is)
 50         continue

 51      continue

       else

         do 61 is=1,nspec
           if(ncspec(is).eq.0) go to 61

           do 60 isign=1,2
           do 60 ie=1,negrid
           do 60 j=1,nlambda
           do 60 i=-ntgridl,ntgrid
             g1(i,j,ie,isign,is)= gnew(i,j,ie,isign,is)
 60        continue

 61     continue

         call integrate

         do 91 is=1,nspec
            if(ncspec(is).eq.0) goto 91

            do 90 j=-ntgridl,ntgrid
               g0eint(j,is)= geint(j,is)
 90         continue

 91      continue

       endif


       if(vnew(2,3).ge.-1.e-33) then
         if(abs(delt-tcomp).le.1.e-33) goto 200

c      calculate tri-diagonal coefficient matrix
           ivb(0)=ntheta/2

         do 100 i=1,ntheta/2+1
            ivb(i)=ntheta/2+2-i
 100     continue

         do 110 i=ntheta/2+2,ntheta
            ivb(i)=i-ntheta/2
 110     continue


         do 120 i=0,ntgrid
           ipointer=mod(i+1,ntheta)
           jend(i)=ng2+ivb(ipointer)
 120     continue
         
         do 130 i=-ntgridl,0
           jend(i)=jend(-i)
 130     continue

         
         do 161 is=1,nspec
           if(abs(ncspec(is)*vnew(2,is)).le.1.e-33) goto 161
           
           do 160 ie=1,negrid
              
cgms           do 140 j=2,nlambda
           do 140 i=-ntgridl,ntgrid
            jendm1 = jend(i)-1
            do 140 j=2,jendm1
              slb0= sqrt( abs(1.-bmag(i)*al(j-1)) )
              slb1= sqrt( abs(1.-bmag(i)*al(j  )) )
              slb2= sqrt( abs(1.-bmag(i)*al(j+1)) )
              
              slbl= (slb1+slb0)/2.
              slbr= (slb1+slb2)/2.

              cc(i,j,ie,is)=-vnew(ie,is)*delt*(1.-slbr*slbr)
     &                      /(slbr-slbl)/(slb2-slb1)
              aa(i,j,ie,is)=-vnew(ie,is)*delt*(1.-slbl*slbl)
     &                      /(slbr-slbl)/(slb1-slb0)
              bb(i,j,ie,is)=1.0-(aa(i,j,ie,is)+cc(i,j,ie,is))
140        continue


c  apply boundary conditions

           do 150 i=-ntgridl,ntgrid
             slb1= sqrt( abs(1.-bmag(i)*al(1)) )
             slb2= sqrt( abs(1.-bmag(i)*al(2)) )

             slbr= (slb1+slb2)/2.

             aa(i,1,ie,is)=0.0
             bb(i,1,ie,is)=vnew(ie,is)*delt*(-1.-slbr)
     &                     /(slb2-slb1)
             cc(i,1,ie,is)=-bb(i,1,ie,is)
             bb(i,1,ie,is)= bb(i,1,ie,is)+1.0


             slb0= sqrt( abs(1.-bmag(i)*al(jend(i)-1)) )

             slbl= slb0/2.
             
             aa(i,jend(i),ie,is)=-0.5*vnew(ie,is)*delt*(1.-slbl*slbl)
     &                             /slbl/slb0
             bb(i,jend(i),ie,is)=-2.0*aa(i,jend(i),ie,is)+1.
             cc(i,jend(i),ie,is)=     aa(i,jend(i),ie,is)
 150       continue

 160     continue

 161     continue


         tcomp=delt


 200     continue


c   solve for gnew(n+1) row by row

         do 301 is=1,nspec
           if(abs(ncspec(is)*vnew(2,is)).le.1.e-33) goto 301
           
           do 210 isign=1,2
           do 210 ie=1,negrid
           do 210 il=1,nlambda
           do 210 j=-ntgridl,ntgrid
              g2(j,il,ie,isign,is)= 0.
 210       continue


           do 300 ie=1,negrid
           do 300 ii=-ntgridl,ntgrid
              
             do 220 k=1,jend(ii)
               a1(k)=  aa(ii,k,ie,is)
               b1(k)=  bb(ii,k,ie,is)
               c1(k)=  cc(ii,k,ie,is)
                d(k)=gnew(ii,k,ie,1,is)
 220         continue

             do 230 k=1,jend(ii)-1
               a1(jend(ii)+k)=  cc(ii,jend(ii)-k,ie,is)
               b1(jend(ii)+k)=  bb(ii,jend(ii)-k,ie,is)
               c1(jend(ii)+k)=  aa(ii,jend(ii)-k,ie,is)
                d(jend(ii)+k)=gnew(ii,jend(ii)-k,ie,2,is)
 230         continue

             mn=2*jend(ii)-1

c   solve tri-diagonal system by tomas algorithm
 
             betaa(1)=b1(1)
             delta(1)= d(1)

             do 240 i=1,mn-1
                  ql(i+1)=a1(i+1)/betaa(i)
               betaa(i+1)=b1(i+1)-ql(i+1)*c1(i)
               delta(i+1)= d(i+1)-ql(i+1)*delta(i)
 240         continue

             x(mn)=delta(mn)/betaa(mn)
             do 250 i=mn-1,1,-1
                x(i)=(delta(i)-c1(i)*x(i+1))/betaa(i)
 250         continue

             do 260 k=1,jend(ii)
                g2(ii,k,ie,1,is)=x(k)
 260         continue

             do 270 k=1,jend(ii)-1
                g2(ii,k,ie,2,is)=x(2*jend(ii)-k)
 270         continue
             g2(ii,jend(ii),ie,2,is)=g2(ii,jend(ii),ie,1,is)
                

 300      continue

 301      continue


       else


         do 501 is=1,nspec
           if(ncspec(is).eq.0) goto 501

           do 500 ie=1,negrid
           do 500 j=1,nlambda

             if((vnew0(is).gt.1.e-33).or.(vnew0(is).lt.-1.e-33)) then
               vnewtemp=delt*vnewfe(e(ie),al(j),vnew0(is))
             else
               vnewtemp=0.
             endif

             factor=(1.-vnewtemp*fv)/(1.+vnewtemp*(1.-fv))
             
             do 500 i=-ntgridl,ntgrid
               g2(i,j,ie,1,is)=gnew(i,j,ie,1,is)*factor
               g2(i,j,ie,2,is)=gnew(i,j,ie,2,is)*factor
 500         continue

501      continue

       endif


       if(vnew(2,3).ge.-1.e-33) then

         do 601 is=1,nspec
           if(ncspec(is).eq.0) go to 601

           do 520 isign=1,2
           do 520 ie=1,negrid
           do 520 j=1,nlambda
           do 520 i=-ntgridl,ntgrid
             g1(i,j,ie,isign,is)=0.
 520       continue

           if(abs(ncspec(is)*vnew(2,is)).le.1.e-33) goto 601
           
           do 560 isign=1,2
           do 560 ie=1,negrid
           do 560 j=1,nlambda
           do 560 i=-ntgridl,ntgrid
             g1(i,j,ie,isign,is)=g2(i,j,ie,isign,is)*sq(i,j,isign)
 560       continue

 601     continue
         
         call lintegrate

         do 651 is=1,nspec
            if(ncspec(is).eq.0) goto 651

            do 650 ie=1,negrid
            do 650 j=-ntgridl,ntgrid
               g2int(j,ie,is)= gint(j,ie,is)
 650        continue

 651     continue

         
         do 701 is=1,nspec
           if(ncspec(is).eq.0) goto 701

           if(abs(vnew(2,is)).le.1.e-33) then
           
              do 670 isign=1,2
              do 670 ie=1,negrid
              do 670 j=1,nlambda
              do 670 i=-ntgridl,ntgrid
                 g1(i,j,ie,isign,is)=gnew(i,j,ie,isign,is)
 670          continue

           else
              
              do 690 isign=1,2
              do 690 ie=1,negrid
              do 690 j=1,nlambda
              do 690 i=-ntgridl,ntgrid
ccc               g1(i,j,ie,isign,is)=g2(i,j,ie,isign,is)+sq(i,j,isign)*
ccc     &                           (g1int(i,ie,is)-g2int(i,ie,is))
ccc     &                           /g3int(i,ie,is)

                 g1(i,j,ie,isign,is)=g2(i,j,ie,isign,is)
 690          continue

           endif

701      continue


       else


         do 751 is=1,nspec
           if(ncspec(is).eq.0) goto 751

           do 750 isign=1,2
           do 750 ie=1,negrid
           do 750 j=1,nlambda
           do 750 i=-ntgridl,ntgrid
             g1(i,j,ie,isign,is)=g2(i,j,ie,isign,is)
 750       continue

 751     continue

         call integrate

         do 761 is=1,nspec
            if(ncspec(is).eq.0) goto 761

            do 760 j=-ntgridl,ntgrid
               g1eint(j,is)= geint(j,is)
 760        continue

 761     continue


         do 801 is=1,nspec
           if(abs(ncspec(is)*vnew(2,is)).le.1.e-33) goto 801

           do 800 isign=1,2
           do 800 ie=1,negrid
           do 800 j=1,nlambda
           do 800 i=-ntgridl,ntgrid
             add=(g0eint(i,is)-g1eint(i,is))*aintnorm(i,is)
             g1(i,j,ie,isign,is)=g2(i,j,ie,isign,is)+add
 800       continue

 801     continue

       endif


       return 
       end subroutine solfp1                      
       

************************************************************************
       subroutine wstarch
************************************************************************

cgms      Use(Common)


       ict=ict+1

       if(ict.eq.1) delt=time

       if(ict.lt.4) return

       anorm=0.
       do 10 j=-ntgridl,ntgrid
         anorm=anorm+cabs(phig(j)+aparg(j))**2
 10    continue
       anorm=sqrt(anorm/(2*ntgrid+1))

       if(ict.eq.4) an1=anorm

       if ((isvar.eq.1).or.(isvar.eq.13)) then
         tp1=tguess1*(an1*exp(gammause*time)/anorm)**power

         do 20 is=1,2
         do 20 ie=1,negrid
crew      signomega
           wstar(ie,is)=signomega*delt*(fprimm(is)+tp1*(e(ie)-1.5))
20       continue
       endif

       if ((isvar.eq.3).or.(isvar.eq.13)) then
         tp3=tguess3*(an1*exp(gammause*time)/anorm)**power

         do 30 ie=1,negrid
crew      signomega
           wstar(ie,3)=signomega*delt*(fprimm(3)+tp3*(e(ie)-1.5))
 30      continue
       endif

       return
       end subroutine wstarch 



************************************************************************
       function vnewfe(edata,aldata,vnew0data)
************************************************************************

cgms      Use(Common)

       real vhat,vestar,vnewf1,vnewfh
       vhat=sqrt(edata)
       vep=2.*abs(vnew0data)

       vestar=vep*zeff*sqrt(amass(3)*temp(1)/temp(3))/(eps**1.5*pp)

       vnewfl=vep*eps*(zeff+hee(vhat))/(vhat**3*((1.-eps-aldata)**2
     &           +1.e-8))*(.111*(absom/vep*zeff)+1.31)/(11.79*absom
     &           /vep*zeff+1)

       vnewfh=1.08*vep*(zeff+0.6289)/1.6289

       vnewfe=(vnewfl+vnewfh*vestar/(zeff*vhat**3)*(zeff+hee(vhat)))
     &         /(1.+vestar/(zeff*vhat**3)*(zeff+hee(vhat)))


       return
       end  function vnewfe   


************************************************************************
       function hee(vhat)
************************************************************************

       real pi,hee
cgms       pi= 3.1415926
       pi=atan2(0.0,-1.0)
       hee=1./pi/vhat*exp(-vhat**2)+(1-1/2./vhat**2)*erf(vhat)

       return
       end function hee


************************************************************************
       function erf(x)
************************************************************************
        real a1,a2,a3,a4,a5,a6
        a1=0.0705230784
        a2=0.0422820123
        a3=0.0092705272
        a4=0.0001520143
        a5=0.0002765672
        a6=0.0000430638

        erf=1.-1./(1.+a1*x+a2*x*x+a3*x**3+a4*x**4+
     &                            a5*x**5+a6*x**6)**16

        return
        end function erf
************************************************************************
       subroutine timeadv
************************************************************************

cgms      Use(Common)

       real dthet(-ntml:ntm)
       complex source(-ntml:ntm,nlm,negrid,2,nspec)


       do 1000 is=1,nspec
         if(ncspec(is).eq.0) goto 1000
         
         do 50 isign=1,2
         do 50 ie=1,negrid
         do 50 il=1,nlambda
         do 50 j= -ntgridl,ntgrid
            source(j,il,ie,isign,is)= 0.
                g2(j,il,ie,isign,is)= 0.
 50      continue

         
         stm=sqrt(temp(is)/amass(is))
         zstm= z(is)/sqrt(temp(is)*amass(is))

c     set up source term
         do 60 isign=1,2
         do 60 ie=1,negrid
         do 60 il=1,nlambda
         do 60 j=-ntgridl,ntgrid
            gnew(j,il,ie,isign,is)=0.
 60      continue


         if (abs(eps).le.1.e-33) then
            lmax=nlambda
         else
            lmax=nlambda-1
         endif

c     compute gyrophase averaged phi=phia, apar=apara
c        (using g1 as their tempary storage)
         do 80 isign=1,2
         do 80 ie=1,negrid
         do 80 il=1,nlambda
         do 80 j=-ntgridl,ntgrid
            g1(j,il,ie,isign,1)=(   fexp(is) *phig(j) +
     >                          (1.-fexp(is))* phi(j) )
     >                          *aj0(j,il,ie,is)

            g1(j,il,ie,isign,2)=(   fexp(is) *aparg(j) +
     >                          (1.-fexp(is))* apar(j) )
     >                          *aj0(j,il,ie,is)
 80      continue


c     compute source term in finite difference equations
         do 100 isign=1,2
         do 100 ie=1,negrid
         do 100 il=1,lmax
         do 100 j=-ntgridl,ntgrid-1
crew omegaprev added
crew fixed error uprim-> delt*uprim
            transfer=-2.*vpar(j,il,ie,isign,is)*
     &                (g1(j+1,il,ie,isign,1)-g1(j,il,ie,isign,1))
     &               -zstm*vpac(j,il,ie,isign)*aj0(j,il,ie,is)*
     &               ( (apar(j+1)+apar(j)-aparg(j+1)-aparg(j))
     >      -zi*delt*omegaprev/2.*
     >                 (apar(j+1)+apar(j)+aparg(j+1)+aparg(j)) )
     &        +zi*(wstar(ie,is)+vpac(j,il,ie,isign)*delt*uprim(is))*
     &                ( g1(j+1,il,ie,isign,1)+ g1(j,il,ie,isign,1)
     &                 -stm*vpac(j,il,ie,isign)*
     &                  (g1(j+1,il,ie,isign,2)+g1(j,il,ie,isign,2)) ) 
     &            -zi*wdrift(j,il,ie,is)*
     &                ( g1(j+1,il,ie,isign,1)+ g1(j,il,ie,isign,1) )

            source(j,il,ie,isign,is)= transfer
 100     continue


         
         if (istep.eq.0) then

            
           do 110 ie=1,negrid
           do 110 il=1,lmax
           do 110 j=-ntgridl,ntgrid-1
              source(j,il,ie,1,is)=source(j  ,il,ie,1,is)
     >                             +   g0(j  ,il,ie,1,is)*b(j,il,ie,is)     
     >                             +   g0(j+1,il,ie,1,is)*a(j,il,ie,is)
 110       continue

           do 115 ie=1,negrid
           do 115 il=1,lmax
           do 115 j=-ntgridl,ntgrid-1
              source(j,il,ie,2,is)=source(j  ,il,ie,2,is)
     >                             +   g0(j  ,il,ie,2,is)*a(j,il,ie,is)     
     >                             +   g0(j+1,il,ie,2,is)*b(j,il,ie,is)     
115        continue

c      special source term for totally trapped particles
           if(abs(eps).le.1.e-33) goto 210

           do 120 iperiod=nperiod10,nperiod2
              j0=iperiod*ntheta
           do 120 isign=1,2
           do 120 ie=1,negrid
crew fixed error uprim-> delt*uprim
              transfer= g0(j0,nlambda,ie,isign,is)*a(j0,nlambda,ie,is)
     &                 -zstm*vpac(j0,nlambda,ie,isign)
     &                       *aj0(j0,nlambda,ie,is)*(apar(j0)-aparg(j0))
     &          +zi*(wstar(ie,is)+vpac(j0,nlambda,ie,isign)*
     &                  delt*uprim(is))*( g1(j0,nlambda,ie,isign,1) 
     &                                   -stm*vpac(j0,nlambda,ie,isign)*
     &                                    g1(j0,nlambda,ie,isign,2) ) 
     &             -zi*wdrift(j0,nlambda,ie,is)
     &                    *g1(j0,nlambda,ie,isign,1)
      
              source(j0,nlambda,ie,isign,is)= transfer
 120       continue


        else

            
           do 130 ie=1,negrid
           do 130 il=1,lmax
           do 130 j=-ntgridl,ntgrid-1
              source(j,il,ie,1,is)= source(j  ,il,ie,1,is)
     >                             +     g(j  ,il,ie,1,is)*b(j,il,ie,is)     
     >                             +     g(j+1,il,ie,1,is)*a(j,il,ie,is)
 130       continue


           do 135 ie=1,negrid
           do 135 il=1,lmax
           do 135 j=-ntgridl,ntgrid-1
              source(j,il,ie,2,is)= source(j  ,il,ie,2,is)
     >                             +     g(j  ,il,ie,2,is)*a(j,il,ie,is)     
     >                             +     g(j+1,il,ie,2,is)*b(j,il,ie,is)     
135        continue


c      special source term for totally trapped particles
           if(abs(eps).le.1.e-33) goto 210

           do 140 iperiod=nperiod10,nperiod2
              j0=iperiod*ntheta
           do 140 isign=1,2
           do 140 ie=1,negrid
crew fixed error uprim->delt*uprim
              transfer= g(j0,nlambda,ie,isign,is)*a(j0,nlambda,ie,is)
     &                 -zstm*vpac(j0,nlambda,ie,isign)
     &                       *aj0(j0,nlambda,ie,is)*(apar(j0)-aparg(j0))
     &              +zi*(wstar(ie,is)+vpac(j0,nlambda,ie,isign)*
     &                  delt*uprim(is))*( g1(j0,nlambda,ie,isign,1) 
     &                                   -stm*vpac(j0,nlambda,ie,isign)*
     &                                    g1(j0,nlambda,ie,isign,2) ) 
     &              -zi*wdrift(j0,nlambda,ie,is)
     &                    *g1(j0,nlambda,ie,isign,1)
      
              source(j0,nlambda,ie,isign,is)= transfer
 140       continue

        endif


c   start time advance using ainv and r maticesdelt*
c   start by initializing trapped regions
         
         do 200 iperiod=nperiod10,nperiod2
            jplus2=iperiod*ntheta +ntheth
         do 200 ie=1,negrid
         do 200 il=2,ntheth
            iln=il+ng2
            jplus=jplus2 -il +1
            source(jplus,iln,ie,2,is)=0.
                g2(jplus,il ,ie,2,is)=1.
 200     continue


 210     continue 


         do 250 isign=1,2
         do 250 ie=1,negrid
         do 250 il=1,nlambda
         do 250 j=-ntgrid, ntgrid
            g1(j,il,ie,isign,is)= 0.
250      continue


c   time advance particles with vpar<0 using r,ainv and source
         do 300 ie=1,negrid
         do 300 il=1,nlambda
         do 300 j=ntgrid-1,-ntgridl,-1
            gnew(j,il,ie,2,is)=  -r(j,il,ie,is)*gnew(j+1,il,ie,2,is)+
     &                         ainv(j,il,ie,is)*source(j,il,ie,2,is)
 300     continue



         if (abs(eps).le.1.e-33) goto 315

         do 310 ie=1,negrid
         do 310 il=2,ntheth+1
         do 310 j=ntgrid-1,-ntgridl,-1
            iln=il+ng2
            g1(j,il,ie,2,is)=-r(j,iln,ie,is)*g1(j+1,il,ie,2,is)+
     &                       g2(j,il,ie,2,is)
 310     continue


315      continue


c  time advance gnew with vpar>0

         if(ipar.eq.0.and.abs(eps).le.1.e-33) goto 425
         
crew  add theta=0 reflection         
         
         if(ipar.eq.1) then
         
c  set up gnew values at theta=0
         do 400 ie=1,negrid
         do 400 il=1,nlambda
400         gnew(0,il,ie,1,is)=gnew(0,il,ie,2,is)


         if(abs(eps).le.1.e-33) goto 425

         do 405 ie=1,negrid
         do 405 il=2,ntheth
            g1(0,il,ie,1,is)=g1(0,il,ie,2,is)
405      continue  

         endif       


c  set up gnew in trapped regions
         do 410 iperiod=nperiod11,nperiod2
            jminus1=iperiod*ntheta -ntheth
         do 410 ie=1,negrid
         do 410 il=2,ntheth
            jminus=jminus1 +il -2
            iln=il+ng2
            source(jminus,iln,ie,1,is)=gnew(jminus+1,iln,ie,2,is)
                g2(jminus,il ,ie,1,is)=  g1(jminus+1,il ,ie,2,is)
410      continue


         do 420 iperiod=nperiod11,nperiod2
            j0=iperiod*ntheta
         do 420 ie=1,negrid
            gnew(j0,nlambda,ie,1,is)=gnew(j0,nlambda,ie,2,is)
420      continue


425      continue


c   time advance using r,ainv and source
         do 500 ie=1,negrid
         do 500 il=1,lmax
         do 500 j=-ntgridl,ntgrid-1
            gnew(j+1,il,ie,1,is)=  -r(j,il,ie,is)*gnew(j,il,ie,1,is)+
     &                           ainv(j,il,ie,is)*source(j,il,ie,1,is)
500      continue


         do 510 ie=1,negrid
         do 510 il=2,ntheth
         do 510 j=-ntgridl,ntgrid-1
            iln=il+ng2
            g1(j+1,il,ie,1,is)=-r(j,iln,ie,is)*g1(j,il,ie,1,is)+
     &                         g2(j,il,ie,1,is)
510      continue



c   add correct amount of homogeneous solution g2
         if(abs(eps).le.1.e-33) goto 580

         do 520 iperiod=nperiod10,nperiod2
            jminus1= iperiod*ntheta -ntheth
            jplus2 = iperiod*ntheta +ntheth
         do 520 ie=1,negrid
         do 520 il=2,ntheth
            iln=il+ng2
            jplus= jplus2 -il +1
            beta1=(gnew(jplus,iln,ie,1,is)-gnew(jplus,iln,ie,2,is))/
     &            (1.-g1(jplus,il,ie,1,is)+1.e-10)
crew3/4/98     &            (1.-g1(jplus,il,ie,1,is))
         do 520 j=jminus1,jplus2
          gnew(j,iln,ie,1,is)=gnew(j,iln,ie,1,is)+beta1*g1(j,il,ie,1,is)
          gnew(j,iln,ie,2,is)=gnew(j,iln,ie,2,is)+beta1*g1(j,il,ie,2,is)
520      continue


c  zero out spurious gnew outside trapped boundary
         do 570 iperiod=nperiod10,nperiod2
            jplus2= iperiod*ntheta +ntheth
         do 570 ie=1,negrid
         do 570 il=2,ntheth+1
            iln=il+ng2
            jplus= jplus2 -il +2
            gnew(jplus,iln,ie,1,is)=0.
 570     continue


         do 575 iperiod=nperiod11,nperiod2
            jminus1= iperiod*ntheta -ntheth
         do 575 ie=1,negrid
         do 575 il=2,ntheth+1
            iln=il+ng2
            jminus= jminus1 +il -2
            gnew(jminus,iln,ie,2,is)=0.
575      continue


580      continue


         go to 1000
crew
crew
crewcc         if (istep.eq.0) goto 1000
crew
crew         diff=0.
crew
crew         do 900 j=-ntgridl,ntgrid
crew            diff=diff+cabs(phig(j))+cabs(aparg(j))
crew 900     continue
crew
crew         diff=(diff/((2*ntgrid+1)*phydif))**power
crew
crew         do 910 j=-ntgridl,ntgrid
crew            dthet(j)=exp(-diff*delt*(1.+shat**2*(theta(j)
crew     >           -theta0)**2))
crew 910     continue
crew
crew         do 930 isign=1,2
crew         do 930 ie=1,negrid
crew         do 930 il=1,nlambda
crew         do 930 j=-ntgridl,ntgrid
crew            gnew(j,il,ie,isign,is)=gnew(j,il,ie,isign,is)*dthet(j)
crew 930     continue
crew
crew
1000   continue


       call solfp1


       do 1101 is=1,nspec
         if(ncspec(is).eq.0) goto 1101

         do 1100 isign=1,2
         do 1100 ie=1,negrid
         do 1100 il=1,nlambda
         do 1100 j=-ntgridl,ntgrid
           gnew(j,il,ie,isign,is)=g1(j,il,ie,isign,is)
1100     continue

1101   continue


       return
       end subroutine timeadv

************************************************************************
       subroutine timeadv_bpar
************************************************************************

cgms      Use(Common)
       real dthet(-ntml:ntm)
       complex source(-ntml:ntm,nlm,negrid,2,nspec)

       

       do 1000 is=1,nspec
         if(ncspec(is).eq.0) goto 1000
         
         do 50 isign=1,2
         do 50 ie=1,negrid
         do 50 il=1,nlambda
         do 50 j= -ntgridl,ntgrid
            source(j,il,ie,isign,is)= 0.
                g2(j,il,ie,isign,is)= 0.
 50      continue

         
         stm=sqrt(temp(is)/amass(is))
         zstm= z(is)/sqrt(temp(is)*amass(is))
         tz=temp(is)/abs(z(is))

c     set up source term
         do 60 isign=1,2
         do 60 ie=1,negrid
         do 60 il=1,nlambda
         do 60 j=-ntgridl,ntgrid
            gnew(j,il,ie,isign,is)=0.
 60      continue


         if (abs(eps).le.1.e-33) then
            lmax=nlambda
         else
            lmax=nlambda-1
         endif

c     compute gyrophase averaged phi=phia, apar=apara
c      1/21/98 added bpar=bpara
c        (using g1 as their tempary storage)
c       1/21/98 must add g1(....,3)
         do 80 isign=1,2
         do 80 ie=1,negrid
         do 80 il=1,nlambda
         do 80 j=-ntgridl,ntgrid
            g1(j,il,ie,isign,1)=(   fexp(is) *phig(j) +
     >                          (1.-fexp(is))* phi(j) )
     >                          *aj0(j,il,ie,is)
     >     +y_bpar*(fexp(is)*bparg(j)+(1.-fexp(is))*bpar(j))
     >             *2.*temp(is)/(z(is))*azj1(j,il,ie,is)*vper2(j,il,ie)
c   the y_bpar line only addition
c 9/3/98 found factor 2. was missed
c 9/18/98 abs(z)->z
            g1(j,il,ie,isign,2)=(   fexp(is) *aparg(j) +
     >                          (1.-fexp(is))* apar(j) )
     >                          *aj0(j,il,ie,is)
     
 80      continue


c     compute source term in finite difference equations
         do 100 isign=1,2
         do 100 ie=1,negrid
         do 100 il=1,lmax
         do 100 j=-ntgridl,ntgrid-1
crew omegaprev added
crew fixed error uprim-> delt*uprim
            transfer=-2.*vpar(j,il,ie,isign,is)*
     &                (g1(j+1,il,ie,isign,1)-g1(j,il,ie,isign,1))
     &               -zstm*vpac(j,il,ie,isign)*aj0(j,il,ie,is)*
     &               ( (apar(j+1)+apar(j)-aparg(j+1)-aparg(j))
     >      -zi*delt*omegaprev/2.*
     >                 (apar(j+1)+apar(j)+aparg(j+1)+aparg(j)) ) 
     &        +zi*(wstar(ie,is)+vpac(j,il,ie,isign)*delt*uprim(is))*
     &                ( g1(j+1,il,ie,isign,1)+ g1(j,il,ie,isign,1)
     &                 -stm*vpac(j,il,ie,isign)*
     &                  (g1(j+1,il,ie,isign,2)+g1(j,il,ie,isign,2)))
     &            -zi*wdrift(j,il,ie,is)*
     &                ( g1(j+1,il,ie,isign,1)+ g1(j,il,ie,isign,1)) 

            source(j,il,ie,isign,is)= transfer
 100     continue


         
         if (istep.eq.0) then

            
           do 110 ie=1,negrid
           do 110 il=1,lmax
           do 110 j=-ntgridl,ntgrid-1
              source(j,il,ie,1,is)=source(j  ,il,ie,1,is)
     >                             +   g0(j  ,il,ie,1,is)*b(j,il,ie,is)     
     >                             +   g0(j+1,il,ie,1,is)*a(j,il,ie,is)
 110       continue

           do 115 ie=1,negrid
           do 115 il=1,lmax
           do 115 j=-ntgridl,ntgrid-1
              source(j,il,ie,2,is)=source(j  ,il,ie,2,is)
     >                             +   g0(j  ,il,ie,2,is)*a(j,il,ie,is)     
     >                             +   g0(j+1,il,ie,2,is)*b(j,il,ie,is)     
115        continue

c      special source term for totally trapped particles
           if(abs(eps).le.1.e-33) goto 210

           do 120 iperiod=nperiod10,nperiod2
              j0=iperiod*ntheta
           do 120 isign=1,2
           do 120 ie=1,negrid
crew fixed error uprim-> delt*uprim
              transfer= g0(j0,nlambda,ie,isign,is)*a(j0,nlambda,ie,is)
     &               -zstm*vpac(j0,nlambda,ie,isign)
     &                   *aj0(j0,nlambda,ie,is)*(apar(j0)-aparg(j0))
     &               +zi*(wstar(ie,is)+vpac(j0,nlambda,ie,isign)*
     &                  delt*uprim(is))*( g1(j0,nlambda,ie,isign,1) 
     &                                -stm*vpac(j0,nlambda,ie,isign)*
     &                                 g1(j0,nlambda,ie,isign,2))
     &             -zi*wdrift(j0,nlambda,ie,is)
     &                    *g1(j0,nlambda,ie,isign,1)
      
              source(j0,nlambda,ie,isign,is)= transfer
 120       continue


        else

            
           do 130 ie=1,negrid
           do 130 il=1,lmax
           do 130 j=-ntgridl,ntgrid-1
              source(j,il,ie,1,is)= source(j  ,il,ie,1,is)
     >                             +     g(j  ,il,ie,1,is)*b(j,il,ie,is)     
     >                             +     g(j+1,il,ie,1,is)*a(j,il,ie,is)
 130       continue


           do 135 ie=1,negrid
           do 135 il=1,lmax
           do 135 j=-ntgridl,ntgrid-1
              source(j,il,ie,2,is)= source(j  ,il,ie,2,is)
     >                             +     g(j  ,il,ie,2,is)*a(j,il,ie,is)     
     >                             +     g(j+1,il,ie,2,is)*b(j,il,ie,is)     
135        continue


c      special source term for totally trapped particles
           if(abs(eps).le.1.e-33) goto 210

           do 140 iperiod=nperiod10,nperiod2
              j0=iperiod*ntheta
           do 140 isign=1,2
           do 140 ie=1,negrid
crew fixed error uprim->delt*uprim
              transfer= g(j0,nlambda,ie,isign,is)*a(j0,nlambda,ie,is)
     &                 -zstm*vpac(j0,nlambda,ie,isign)
     &                     *aj0(j0,nlambda,ie,is)*(apar(j0)-aparg(j0))
     &              +zi*(wstar(ie,is)+vpac(j0,nlambda,ie,isign)*
     &                  delt*uprim(is))*( g1(j0,nlambda,ie,isign,1) 
     &                                 -stm*vpac(j0,nlambda,ie,isign)*
     &                                  g1(j0,nlambda,ie,isign,2))
     &              -zi*wdrift(j0,nlambda,ie,is)
     &                    *g1(j0,nlambda,ie,isign,1)
      
              source(j0,nlambda,ie,isign,is)= transfer
 140       continue

        endif


c   start time advance using ainv and r maticesdelt*
c   start by initializing trapped regions
         
         do 200 iperiod=nperiod10,nperiod2
            jplus2=iperiod*ntheta +ntheth
         do 200 ie=1,negrid
         do 200 il=2,ntheth
            iln=il+ng2
            jplus=jplus2 -il +1
            source(jplus,iln,ie,2,is)=0.
                g2(jplus,il ,ie,2,is)=1.
 200     continue


 210     continue 


         do 250 isign=1,2
         do 250 ie=1,negrid
         do 250 il=1,nlambda
         do 250 j=-ntgrid, ntgrid
            g1(j,il,ie,isign,is)= 0.
250      continue


c   time advance particles with vpar<0 using r,ainv and source
         do 300 ie=1,negrid
         do 300 il=1,nlambda
         do 300 j=ntgrid-1,-ntgridl,-1
            gnew(j,il,ie,2,is)=  -r(j,il,ie,is)*gnew(j+1,il,ie,2,is)+
     &                         ainv(j,il,ie,is)*source(j,il,ie,2,is)
 300     continue



         if (abs(eps).le.1.e-33) goto 315

         do 310 ie=1,negrid
         do 310 il=2,ntheth+1
         do 310 j=ntgrid-1,-ntgridl,-1
            iln=il+ng2
            g1(j,il,ie,2,is)=-r(j,iln,ie,is)*g1(j+1,il,ie,2,is)+
     &                       g2(j,il,ie,2,is)
 310     continue


315      continue


c  time advance gnew with vpar>0

         if(ipar.eq.0.and.abs(eps).le.1.e-33) goto 425
         
crew  add theta=0 reflection         
         
         if(ipar.eq.1) then
         
c  set up gnew values at theta=0
         do 400 ie=1,negrid
         do 400 il=1,nlambda
400         gnew(0,il,ie,1,is)=gnew(0,il,ie,2,is)


         if(abs(eps).le.1.e-33) goto 425

         do 405 ie=1,negrid
         do 405 il=2,ntheth
            g1(0,il,ie,1,is)=g1(0,il,ie,2,is)
405      continue  

         endif       


c  set up gnew in trapped regions
         do 410 iperiod=nperiod11,nperiod2
            jminus1=iperiod*ntheta -ntheth
         do 410 ie=1,negrid
         do 410 il=2,ntheth
            jminus=jminus1 +il -2
            iln=il+ng2
            source(jminus,iln,ie,1,is)=gnew(jminus+1,iln,ie,2,is)
                g2(jminus,il ,ie,1,is)=  g1(jminus+1,il ,ie,2,is)
410      continue


         do 420 iperiod=nperiod11,nperiod2
            j0=iperiod*ntheta
         do 420 ie=1,negrid
            gnew(j0,nlambda,ie,1,is)=gnew(j0,nlambda,ie,2,is)
420      continue


425      continue


c   time advance using r,ainv and source
         do 500 ie=1,negrid
         do 500 il=1,lmax
         do 500 j=-ntgridl,ntgrid-1
            gnew(j+1,il,ie,1,is)=  -r(j,il,ie,is)*gnew(j,il,ie,1,is)+
     &                           ainv(j,il,ie,is)*source(j,il,ie,1,is)
500      continue


         do 510 ie=1,negrid
         do 510 il=2,ntheth
         do 510 j=-ntgridl,ntgrid-1
            iln=il+ng2
            g1(j+1,il,ie,1,is)=-r(j,iln,ie,is)*g1(j,il,ie,1,is)+
     &                         g2(j,il,ie,1,is)
510      continue



c   add correct amount of homogeneous solution g2
         if(abs(eps).le.1.e-33) goto 580

         do 520 iperiod=nperiod10,nperiod2
            jminus1= iperiod*ntheta -ntheth
            jplus2 = iperiod*ntheta +ntheth
         do 520 ie=1,negrid
         do 520 il=2,ntheth
            iln=il+ng2
            jplus= jplus2 -il +1
            beta1=(gnew(jplus,iln,ie,1,is)-gnew(jplus,iln,ie,2,is))/
     &            (1.-g1(jplus,il,ie,1,is)+1.e-10)
crew3/3/98     &            (1.-g1(jplus,il,ie,1,is))
         do 520 j=jminus1,jplus2
          gnew(j,iln,ie,1,is)=gnew(j,iln,ie,1,is)+beta1*g1(j,il,ie,1,is)
          gnew(j,iln,ie,2,is)=gnew(j,iln,ie,2,is)+beta1*g1(j,il,ie,2,is)
520      continue


c  zero out spurious gnew outside trapped boundary
         do 570 iperiod=nperiod10,nperiod2
            jplus2= iperiod*ntheta +ntheth
         do 570 ie=1,negrid
         do 570 il=2,ntheth+1
            iln=il+ng2
            jplus= jplus2 -il +2
            gnew(jplus,iln,ie,1,is)=0.
 570     continue


         do 575 iperiod=nperiod11,nperiod2
            jminus1= iperiod*ntheta -ntheth
         do 575 ie=1,negrid
         do 575 il=2,ntheth+1
            iln=il+ng2
            jminus= jminus1 +il -2
            gnew(jminus,iln,ie,2,is)=0.
575      continue


580      continue


         go to 1000

cc         if (istep.eq.0) goto 1000

crew         diff=0.
crew
crew         do 900 j=-ntgridl,ntgrid
crew            diff=diff+cabs(phig(j))+cabs(aparg(j))
crew 900     continue
crew
crew         diff=(diff/((2*ntgrid+1)*phydif))**power
crew
crew         do 910 j=-ntgridl,ntgrid
crew            dthet(j)=exp(-diff*delt*(1.+shat**2*(theta(j)
crew     >           -theta0)**2))
crew 910     continue
crew
crew         do 930 isign=1,2
crew         do 930 ie=1,negrid
crew         do 930 il=1,nlambda
crew         do 930 j=-ntgridl,ntgrid
crew            gnew(j,il,ie,isign,is)=gnew(j,il,ie,isign,is)*dthet(j)
crew 930     continue
crew
crew
1000   continue
  


       call solfp1


       do 1101 is=1,nspec
         if(ncspec(is).eq.0) goto 1101

         do 1100 isign=1,2
         do 1100 ie=1,negrid
         do 1100 il=1,nlambda
         do 1100 j=-ntgridl,ntgrid
           gnew(j,il,ie,isign,is)=g1(j,il,ie,isign,is)
1100     continue

1101   continue


       return
       end subroutine timeadv_bpar



************************************************************************
       subroutine getfield
************************************************************************

cgms      Use(Common)
crew________________geo_______________________________________________
cgms      Use(Geo)
       complex gk(2*(ntml+ntm+1),2*(ntml+ntm+1)), gf(2*(ntml+ntm+1)),
     >                                       u(2*(ntml+ntm+1))


       if (iphi.eq.0) then

          call getan

          do 100 mi= 1,2*(ntgridl+ntgrid+1),2
             j= mi/2 -ntgridl

             gridfac1= 1.
crew          if ((j.eq.-ntgrid).or.(j.eq.ntgrid))  gridfac1= gridfac
            if (ipar.eq.0) then
             if ((j.eq.-ntgridl).or.(j.eq.ntgrid))  gridfac1= gridfac
            endif
            if (ipar.eq.1) then
             if (j.eq.ntgrid) gridfac1= gridfac
            endif
                
             do 50 mj= 1,2*(ntgridl+ntgrid+1)    
                gk(mi,mj)= -am(mi,mj)*test1
 50          continue

  
             gk(mi,mi)= (xteti +gamtot(j)+ alphal)*gridfac1
     >                 -am(mi,mi)*test1 

             gf(mi)= (antot(j) +alphal*phig(j))*test1
 100      continue

          do 200 mi= 2,2*(ntgridl+ntgrid+1),2
             j= (mi-1)/2 -ntgridl

             gridfac1= 1.
crew          if ((j.eq.-ntgrid).or.(j.eq.ntgrid))  gridfac1= gridfac
            if (ipar.eq.0) then
             if ((j.eq.-ntgridl).or.(j.eq.ntgrid))  gridfac1= gridfac
            endif
            if (ipar.eq.1) then
             if (j.eq.ntgrid) gridfac1= gridfac
            endif


             do 150 mj= 1,2*(ntgridl+ntgrid+1)    
                gk(mi,mj)= -am(mi,mj)*test2
 150         continue
             gyro_fix=1.0
cgms             if(igyro_fix.eq.1)gyro_fix=bmag(j)
             gk(mi,mi)= (((aky/gyro_fix)**2) * 
     >                  (1. + (shat*(theta(j)-theta0)
     >                   - shift*sin(theta(j)) )**2) 
     >                   +ala)*gridfac1 -am(mi,mi)*test2
crew error fixed above shatxshift 4/4/96
crew________________geo_______________________________________________
     
          if(igeo.eq.1)
     >    gk(mi,mi)= ((aky**2)*qrat_geo(j)**2*b_unit**2/b2_ave_geo*
     >               (1. +  kxoky_geo(j)**2) 
     >               +ala)*gridfac1 -am(mi,mi)*test2
c note: "beta" in gstotal is 8.*pi*ne*ti/b2_ave_geo
c so we do not divide by b_geo(j)**2
crew________________geo_______________________________________________
          
             gf(mi)= (antota(j) +ala*aparg(j))*test2
 200      continue
  
c         if(ntgridl.eq.0.and.ntgrid.le.48) then
c          do mi=1,147
c           gfss(mi)=gf(mi)
c           do mj=1,147
c             gkss(mi,mj)=gk(mi,mj)
c           enddo
c          enddo
c         endif


          ndim= 2*(ntml+ntm+1)
          neqn= 2*(ntgridl+ntgrid+1)
          call gauss(gk,ndim,gf,neqn,u)
  
c        if(ntgridl.eq.0.and.ntgrid.le.48) then
c         do mi=1,147
c           uss(mi)=u(mi)
c         enddo  
c        endif

             miphi = -1
             miapar= 0 
          do 300 j= -ntgridl,ntgrid
             miphi = miphi  +2
             miapar= miapar +2
          
             phi (j)= u(miphi )
             apar(j)= u(miapar)
 300      continue


       endif



       if (iphi.eq.1) then

          do 400 j= -ntgridl,ntgrid
             phi (j)= exp(-((theta(j)-theta0)/width0)**2)*c
cc             phi(j)= cos(theta(j))*c
             
cc             if ( abs(theta(j)-theta0).lt.pi ) then
cc                phi(j)= (theta(j)-theta0)**2 -pi*pi
cc             else
cc                phi(j)= 0.0
cc             endif

             apar(j)= 0.
 400      continue

       endif



       return
       end subroutine getfield
************************************************************************
       subroutine getfield_bpar
************************************************************************

cgms      Use(Common)
crew________________geo_______________________________________________
cgms      Use(Geo)
       complex gk(3*(ntml+ntm+1),3*(ntml+ntm+1)), gf(3*(ntml+ntm+1)),
     >                                       u(3*(ntml+ntm+1))
       real ak4(-ntml:ntm)
       real wkspce(3*(ntml+ntm+1))


       if (iphi.eq.0) then

          call getan
          
          
          mspan=(ntgrid+ntgridl+1)
          mspanp1=mspan+1
          mspan2=2*mspan
          mspan2p1=mspan2+1
          mspan3=3*mspan
  
          nspan=2
          if(i_bpar.eq.1) nspan=3
          
c          mjjudge.eq.1 then
c            j= mj-(mspan-ntgrid-1)-1           

c          mjjudge.eq.2 then
c            j= mj-mspan-(mspan-ntgrid-1)-1
           
c           mjjudge.eq.3 then
c            j= mj-2*mspan-(mspan-ntgrid-1)-1
           
          
c   mjjudge 1,(1,2,3)

          do 100 mi= 1,mspan
             j= mi-(mspan-ntgrid-1)-1

             gridfac1= 1.
crew          if ((j.eq.-ntgrid).or.(j.eq.ntgrid))  gridfac1= gridfac
            if (ipar.eq.0) then
             if ((j.eq.-ntgridl).or.(j.eq.ntgrid))  gridfac1= gridfac
            endif
            if (ipar.eq.1) then
             if (j.eq.ntgrid) gridfac1= gridfac
            endif
                
             do 50 mj= 1, nspan*mspan   
                gk(mi,mj)= -am(mi,mj)*test1
 50          continue

             gk(mi,mi)= (xteti +gamtot(j)+ alphal)*gridfac1 
     >                 -am(mi,mi)*test1
 
             mj=mi+mspan2
             gk(mi,mj)= (gamtotpb(j)-am(mi,mj))*test2   

          
             gf(mi)= (antot(j) +alphal*phig(j))*test1
 100      continue
  
c            write(6,*) "mjjudge=1"
 
c mjjudge 2,(1,2,3)

          do 200 mi= mspanp1,mspan2
             j= mi-mspan-(mspan-ntgrid-1)-1
             
             gridfac1= 1.
crew          if ((j.eq.-ntgrid).or.(j.eq.ntgrid))  gridfac1= gridfac
            if (ipar.eq.0) then
             if ((j.eq.-ntgridl).or.(j.eq.ntgrid))  gridfac1= gridfac
            endif
            if (ipar.eq.1) then
             if (j.eq.ntgrid) gridfac1= gridfac
            endif


             do 150 mj= 1,nspan*mspan   
                gk(mi,mj)= -am(mi,mj)*test2
 150         continue
             gyro_fix=1.0
cgms             if(igyro_fix.eq.1)gyro_fix=bmag(j)
             gk(mi,mi)= (((aky/gyro_fix)**2) * 
     >                  (1. + ( (shat*(theta(j)-theta0)
     >                   - shift*sin(theta(j)) ) )**2) 
     >                   +ala)*gridfac1 -am(mi,mi)*test2
crew error fixed above shatxshift 4/4/96
crew________________geo_______________________________________________
     
          if(igeo.eq.1)
     >    gk(mi,mi)= ((aky**2)*qrat_geo(j)**2*b_unit**2/b2_ave_geo*
     >               (1. +  kxoky_geo(j)**2) 
     >               +ala)*gridfac1 -am(mi,mi)*test2
c note: "beta" in gstotal is 8.*pi*ne*ti/b2_ave_geo
c so we do not divide by b_geo(j)**2
crew________________geo_______________________________________________
          
             gf(mi)= (antota(j) +ala*aparg(j))*test2
 200      continue
  
c          write(6,*) "mjjudge=2"
 
          if(nspan.eq.3) then
c mjjudge 3,(1,2,3)

          do 300 mi= mspan2p1,mspan3
              j= mi-2*mspan-(mspan-ntgrid-1)-1
             
             gridfac1= 1.
crew          if ((j.eq.-ntgrid).or.(j.eq.ntgrid))  gridfac1= gridfac
            if (ipar.eq.0) then
             if ((j.eq.-ntgridl).or.(j.eq.ntgrid))  gridfac1= gridfac
            endif
            if (ipar.eq.1) then
             if (j.eq.ntgrid) gridfac1= gridfac
            endif

             do 240 mj= 1,3*mspan   
                gk(mi,mj)= -am(mi,mj)*test2
 240         continue
  
             mj=mi-mspan2
             gk(mi,mj)= gamtotb(j)-am(mi,mj)*test2
             gyro_fix=1.0
cgms             if(igyro_fix.eq.1)gyro_fix=bmag(j)
             gk(mi,mi)=(bmag(j)/gyro_fix)**2
     >                +gamtotbb(j) -am(mi,mi)*test2
crew________________geo_______________________________________________
     
          if(igeo.eq.1) 
     > gk(mi,mi)= b_geo(j)**2*b_unit**2/b2_ave_geo+gamtotbb(j)
     >        -am(mi,mi)*test2

c note: "beta" in gstotal is 8.*pi*ne*ti/b2_ave_geo
c so we do not divide by b_geo(j)**2
crew________________geo_______________________________________________

             ak4(j)=gk(mi,mi)+am(mi,mi)*test2
          
             gf(mi)= (antotb(j) +ala*aparg(j))*test2
 300      continue
  
c         write(6,*) "mjjudge=3"

         endif
  
c         if(ntgridl.eq.0.and.ntgrid.le.48) then
c          do mi=1,147
c           gfss(mi)=gf(mi)
c           do mj=1,147
c             gkss(mi,mj)=gk(mi,mj)
c           enddo
c          enddo
c         endif


          ndim= 3*(ntml+ntm+1)
          neqn= nspan*mspan
c          write(6,*) "ndim=",ndim,"neqn=",neqn
cgms overwrite i_solve option
       i_solve = 0
          if(i_solve.eq.0) call gauss(gk,ndim,gf,neqn,u)
  
cgms         if(i_solve.eq.1)
cgms    >     call f04ade(gk,ndim,gf,ndim,neqn,1,u,ndim,wkspce,ifail)

c        if(ntgridl.eq.0.and.ntgrid.le.48) then
c         do mi=1,147
c           uss(mi)=u(mi)
c         enddo  
c        endif 

             
             mi1=0
             mi2=mspan
          do 400 j= -ntgridl,ntgrid
             mi1 = mi1+1
             mi2 = mi2+1
             phi (j)= u(mi1 )
             apar(j)= u(mi2)
             bpar(j)= 0.
 400      continue
  
          if(nspan.eq.3) then
             mi3=mspan2
          do 401 j= -ntgridl,ntgrid
             mi3 = mi3+1
             bpar(j)= u(mi3)
 401      continue
          endif
      



       endif



       if (iphi.eq.1) then

          do 500 j= -ntgridl,ntgrid
             phi (j)= exp(-((theta(j)-theta0)/width0)**2)*c
cc             phi(j)= cos(theta(j))*c
             
cc             if ( abs(theta(j)-theta0).lt.pi ) then
cc                phi(j)= (theta(j)-theta0)**2 -pi*pi
cc             else
cc                phi(j)= 0.0
cc             endif

             apar(j)= 0.
             bpar(j)= 0.
 500      continue

       endif



       return
       end subroutine getfield_bpar
************************************************************************
       subroutine getfield_bpar2
************************************************************************

cgms      Use(Common)
crew________________geo_______________________________________________
cgms      Use(Geo)
       complex gk(3*(ntm+1),3*(ntm+1)), gf(3*(ntm+1)),
     >                                       u(3*(ntm+1))
       real ak4(0:ntm)

       if(ntgridl.ne.0) write(6,*) 'ntgridl.ne.0 do not use i_bpar=2'

       if (iphi.eq.0) then

          call getan


          mspan=(ntgrid+1)

          mijudge=1
          do 100  mi=1,3*(ntgrid+1),3
            j=(mi+3-mijudge)/3-1

             gridfac1= 1.
crew          if ((j.eq.-ntgrid).or.(j.eq.ntgrid))  gridfac1= gridfac
            if (ipar.eq.0) then
             if ((j.eq.-ntgridl).or.(j.eq.ntgrid))  gridfac1= gridfac
            endif
            if (ipar.eq.1) then
             if (j.eq.ntgrid) gridfac1= gridfac
            endif

           do 50 mj=1,3*(ntgrid+1)
             gk(mi,mj)=-am(mi,mj)*test1
50         continue


             gk(mi,mi)= (xteti +gamtot(j)+ alphal)*gridfac1
     >                 -am(mi,mi)*test1

            do 60 mj=3,3*(ntgrid+1),3
             gk(mi,mj)=gk(mi,mj)+gamtotpb(j)*gridfac1*test2
60         continue 

             gf(mi)= (antot(j) +alphal*phig(j))*test1
100      continue

          mijudge=2
          do 200  mi=2,3*(ntgrid+1),3
           j=(mi+3-mijudge)/3-1
 
            gridfac=1.
crew          if ((j.eq.-ntgrid).or.(j.eq.ntgrid))  gridfac1= gridfac
            if (ipar.eq.0) then
             if ((j.eq.-ntgridl).or.(j.eq.ntgrid))  gridfac1= gridfac
            endif
            if (ipar.eq.1) then
             if (j.eq.ntgrid) gridfac1= gridfac
            endif
  
           do 150 mj=1,3*(ntgrid+1)
             gk(mi,mj)=-am(mi,mj)*test2
150         continue
             gyro_fix=1.0
cgms             if(igyro_fix.eq.1)gyro_fix=bmag(j)
             gk(mi,mi)= (((aky/gyro_fix)**2) *
     >                  (1. + ( (shat*(theta(j)-theta0)
     >                   - shift*sin(theta(j)) ) )**2)
     >                   +ala)*gridfac1 -am(mi,mi)*test2
crew error fixed above shatxshift 4/4/96
crew________________geo_______________________________________________

          if(igeo.eq.1)
     >    gk(mi,mi)= ((aky**2)*qrat_geo(j)**2*b_unit**2/b2_ave_geo*
     >               (1. +  kxoky_geo(j)**2)
     >               +ala)*gridfac1 -am(mi,mi)*test2
c note: "beta" in gstotal is 8.*pi*ne*ti/b2_ave_geo
c so we do not divide by b_geo(j)**2
crew________________geo_______________________________________________ 

             gf(mi)= (antota(j) +ala*aparg(j))*test2           
200        continue

          mijudge=3
          do 300 mi=3,3*(ntgrid+1),3
           j=(mi+3-mijudge)/3-1 
 
            gridfac=1.
crew          if ((j.eq.-ntgrid).or.(j.eq.ntgrid))  gridfac1= gridfac
            if (ipar.eq.0) then
             if ((j.eq.-ntgridl).or.(j.eq.ntgrid))  gridfac1= gridfac
            endif
            if (ipar.eq.1) then
             if (j.eq.ntgrid) gridfac1= gridfac
            endif 
         
c            do 350 mj=1,3*(ntgrid+1)
c             gk(mi,mj)= gamtotb(j)-am(mi,mj)*test2
c350         continue

            do 350 mj=1,3*(ntgrid+1)
             gk(mi,mj)= -am(mi,mj)*test2
350         continue

            do 360 mj=1,3*(ntgrid+1),3
             gk(mi,mj)=gk(mi,mj)+gamtotb(j)*test2
360         continue
             gyro_fix=1.0
cgms             if(igyro_fix.eq.1)gyro_fix=bmag(j)
             gk(mi,mi)= (bmag(j)/gyro_fix)**2
     >             +gamtotbb(j) -am(mi,mi)*test2
crew________________geo_______________________________________________

          if(igeo.eq.1)
     > gk(mi,mi)=b_geo(j)**2*b_unit**2/b2_ave_geo
     >   +gamtotbb(j) -am(mi,mi)*test2

c note: "beta" in gstotal is 8.*pi*ne*ti/b2_ave_geo
c so we do not divide by b_geo(j)**2
crew________________geo_______________________________________________

          ak4(j)=gk(mi,mi)+am(mi,mi)*test2

             gf(mi)= (antotb(j) +ala*aparg(j))*test2  

300      continue

c         if(ntgridl.eq.0.and.ntgrid.le.48) then
c          do mi=1,147
c           gfss(mi)=gf(mi)
c           do mj=1,147
c             gkss(mi,mj)=gk(mi,mj)
c           enddo
c          enddo
c         endif



          ndim= 3*(ntm+1)
          neqn= 3*(ntgrid+1)
c          write(6,*) "ndim=",ndim,"neqn=",neqn
          call gauss(gk,ndim,gf,neqn,u)

c        if(ntgridl.eq.0.and.ntgrid.le.48) then
c          do mi=1,147
c           uss(mi)=u(mi)
c          enddo 
c        endif

          mione=-2
          mitwo=-1
          mithree=0
          do 400 j= 0,ntgrid
             mione=mione+3
             mitwo=mitwo+3
             mithree=mithree+3
             phi (j)= u(mione)
             apar(j)= u(mitwo)
             bpar(j)= u(mithree)
 400      continue

          do 402 j=0,ntgrid
           bpar_mhd(j)=-gamtotb(j)*phi(j)/ak4(j)
 402      continue
          if(i_mhd.eq.1) then
          do 403 j=-ntgridl,ntgrid
           bpar(j)=bpar_mhd(j)
 403      continue
          endif 


       endif



       if (iphi.eq.1) then

          do 500 j= -ntgridl,ntgrid
             phi (j)= exp(-((theta(j)-theta0)/width0)**2)*c
cc             phi(j)= cos(theta(j))*c

cc             if ( abs(theta(j)-theta0).lt.pi ) then
cc                phi(j)= (theta(j)-theta0)**2 -pi*pi
cc             else
cc                phi(j)= 0.0
cc             endif

             apar(j)= 0.
             bpar(j)= 0.
 500      continue

       endif



       return
       end subroutine getfield_bpar2     

************************************************************************
      subroutine getan
************************************************************************
c 1/21/98 antotb i_bpar=1 section added

cgms      Use(Common)



       do 10 j=-ntgridl,ntgrid
          antot (j)= 0.0
          antota(j)= 0.0
          antotb(j)= 0.0
 10    continue

c get antot

       do 51 is=1,nspec
         if(ncspec(is).eq.0) goto 51

         do 50 isign=1,2
         do 50 ie=1,negrid
         do 50 il=1,nlambda
         do 50 j=-ntgridl,ntgrid
           g1(j,il,ie,isign,is)=gnew(j,il,ie,isign,is)*aj0(j,il,ie,is)
50       continue

51     continue

       call integrate

       
       do 101 is= 1,nspec
         if(ncspec(is).eq.0) goto 101

         do 100 j=-ntgridl,ntgrid
            antot(j)= antot(j) +z(is)*an(is)*geint(j,is)
100      continue

101    continue

c get antota

       do 161 is=1,nspec
         if(ncspec(is).eq.0) goto 161

         stm=sqrt(temp(is)/amass(is))
         do 160 isign=1,2
         do 160 ie=1,negrid
         do 160 il=1,nlambda
         do 160 j=-ntgridl,ntgrid
           g1(j,il,ie,isign,is)=gnew(j,il,ie,isign,is)*aj0(j,il,ie,is)
     >                          *stm*vpa(j,il,ie,isign)
160      continue

161    continue

       call integrate

       
       do 201 is= 1,nspec
         if(ncspec(is).eq.0) goto 201

         do 200 j=-ntgridl,ntgrid
            antota(j)= antota(j) +2.*beta*z(is)*an(is)*geint(j,is)
200      continue

201    continue

c get antotb

      if(i_bpar.ge.1) then

       do 301 is=1,nspec
         if(ncspec(is).eq.0) goto 301

         stm=abs(z(is))/amass(is)
         do 300 isign=1,2
         do 300 ie=1,negrid
         do 300 il=1,nlambda
         do 300 j=-ntgridl,ntgrid
         g1(j,il,ie,isign,is)=gnew(j,il,ie,isign,is)*temp(is)*
     >     vper2(j,il,ie)*azj1(j,il,ie,is)

300      continue

301    continue

       call integrate

       
       do 401 is= 1,nspec
         if(ncspec(is).eq.0) goto 401

         do 400 j=-ntgridl,ntgrid
          antotb(j)= antotb(j) +(-1.)*beta*x_bpar*an(is)*geint(j,is)
400      continue
c 9/18/98  z factor removed

401    continue

      endif
      
       return
       end subroutine getan


c-----------------------------------------------------------------------
c
c   routine name       - gausse
c
c-----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   purpose            - solution of a full matrix linear system 
c                        of equations using gauss elimination
c                        without pivoting
c
c   usage              - call gauss(gk,igk,gf,n,u)
c
c   arguments in: gk   - n by n complex matrix
c                 igk  - row dimension of matrix gk exactly as
c                        specified in the dimension statement in the
c                        calling program
c                 gf   - complex vector of length n
c                 n    - number of equations
c             out:u    - complex vector of length n containing
c                        the solution
c
c   required routines  - rhsub,tri
c
c-----------------------------------------------------------------------

      subroutine gauss(gk,igk,gf,n,u)

      integer igk,n
      complex gk(igk,*),gf(*),u(*)


      call tri(gk,igk,n)
      call rhsub(gk,gf,igk,n,  u)

      return
      end subroutine gauss

c-----------------------------------------------------------------------
c
c   routine name       - tri
c
c-----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   purpose            - forward gauss elimination
c
c   usage              - call tri(a,n,m)
c
c   arguments in: a    - m by m matrix
c                 n    - row dimension of matrix a exactly as
c                        specified in the dimension statement in the
c                        calling program
c                 m    - number of equations
c             out:a    - a after gauss factorization
c
c   required routines  - none
c
c-----------------------------------------------------------------------

      subroutine tri(a,n,m)

      integer n,m,m1,j1,i,k
      complex a(n,n)
      complex fac

      m1=m-1
      tiny=1.e-30
      zero= 0.d0

c.....eliminate degree of freedom i
      do 30  i=1,m1
c.....check for excessively small pivot
      if(abs(a(i,i)).lt.tiny) go to 99
      j1=i+1
c.....modify row j
      do 20  j=j1,m
      if(a(j,i).eq.zero) go to 20
      fac=a(j,i)/a(i,i)
      do 10 k=j1,m
 10    a(j,k)=a(j,k)-a(i,k)*fac
 20    continue
 30    continue
      return
 99    write (6,100)i,abs(a(i,i))
100   format(/' reduction failed due to small pivot'/
     >   ' equation no.',i5,' pivot ',e10.3)
      stop
      end subroutine tri


c-----------------------------------------------------------------------
c
c   routine name       - rhsub
c
c-----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - april 2, 1988
c
c   purpose            - backward gauss substitution
c
c   usage              - call rhsub(a,x,b,n,m)
c
c   arguments in: a    - m by m matrix
c                 b    - vector of length m containing the right
c                        hand side
c                 n    - row dimension of matrix a exactly as
c                        specified in the dimension statement in the
c                        calling program
c                 m    - number of equations
c             out:x    - vector of length m containing the solution
c
c   required routines  - none
c
c-----------------------------------------------------------------------

      subroutine rhsub(a,b,n,m,x)

c
c   purpose            - backward gauss substitution
c
c   usage              - call rhsub(a,x,b,n,m)
c
c   arguments in: a    - m by m matrix


      integer n,m,i,j1,m1,j,ib
      complex a(n,*),x(*),b(*)

      m1=m-1
c.....begin forward reduction of right hand side
      do 20  i=1,m1
      j1=i+1
      do 10 j=j1,m
10    b(j)=b(j)-b(i)*a(j,i)/a(i,i)
20    continue
c.....begin back substitution
      x(m)=b(m)/a(m,m)
      do 40 i=1,m1
      ib=m-i
      j1=ib+1
      do 30  j=j1,m
30    b(ib)=b(ib)-a(ib,j)*x(j)
      x(ib)=b(ib)/a(ib,ib)
40    continue

      return
      end subroutine rhsub

************************************************************************
      subroutine output
************************************************************************

cgms      Use(Common)
      complex gtime

      
      if (istep.ge.nwrite) then

         do 500 is= 1,nspec
            if (ncspec(is).eq.0) goto 500
         
            do 400 ie= 1,negrid
               if (nce(ie).eq.0) goto 400
        
               ichannel= is*10 +ie

               write(ichannel,*) '***** is= ',is,' ie= ',ie,' *****'
               
               write(ichannel,*) '------------------'
               write(ichannel,*) 'istep= ',istep
               write(ichannel,*) ' '

               do 100 ig= 0,ntgrid
               do 100 il= 1,nlambda
               do 100 isign= 1,2
cc                  gtime= (3.-2*isign)*sqrt(abs(1.-al(il)*bmag(ig)))*
cc     >                 exp( -2*vnew(ie,is)*time )
cc                  gtime= (1. -3.*(1.-al(il)*bmag(ig)))*
cc     >                  exp( -6.*vnew(ie,is)*time )

                  write(ichannel,*) 'ig=',ig,' il=',il,' isign=',isign,
     >                              '  g=',g(ig,il,ie,isign,is)
cc                  write(ichannel,*) '        ',' gt=',gtime/c

 100           continue

 400     continue

 500  continue


      do 551 is= 1,nspec
      do 551 ie= 1,negrid
         ncse= ncspec(is) * nce(ie)
         if (ncse.eq.0) goto 551

         write(7,*) '***** is= ',is,' ie= ',ie,' *****'
               
         write(7,*) '------------------'
         write(7,*) 'istep= ',istep
         write(7,*) ' '

cc         gtime= termi(ie,is)*exp(-2.*vnew(ie,is)*time)
cc         gtime= termi(ie,is)*exp(-6.*vnew(ie,is)*time)

         write(7,*) 'termi=',termi(ie,is)
         write(7,*) 'term =',term (ie,is)
cc         write(7,*) '  gtime =',gtime
         write(7,*) ' '

         write(7,*) '------------------'
         write(7,*) 'istep= ',0
         write(7,*) ' '

         write(7,*) ' '
         write(7,*) '***** is= ',is,' ie= ',ie,' *****'
         write(7,*) ' '

         do 550 ig= -ntgridl,ntgrid
            write(7,*) 'ig=',ig,'  ginti=',ginti(ig,ie,is)
 550     continue

 551  continue


      do 601 is= 1,nspec
      do 601 ie= 1,negrid
         ncse= ncspec(is) * nce(ie)
         if (ncse.eq.0) goto 601
               
         write(7,*) '------------------'
         write(7,*) 'istep= ',istep
         write(7,*) ' '

         write(7,*) ' '
         write(7,*) '***** is= ',is,' ie= ',ie,' *****'
         write(7,*) ' '

         do 600 ig= -ntgridl,ntgrid
            write(7,*) 'ig=',ig,'  gint=',gint(ig,ie,is)
 600     continue

 601  continue


      endif


      return
      end subroutine output
      
*********************************************************************
       subroutine getterm
*********************************************************************

cgms      Use(Common)
       complex term1, term2


       do 21 is= 1,nspec
          if(ncspec(is).eq.0) goto 21

          do 20 isign= 1,2
          do 20 ie= 1,negrid
          do 20 il= 1,nlambda
          do 20 ig= -ntgridl,ntgrid
             g1(ig,il,ie,isign,is)= g(ig,il,ie,isign,is)
 20       continue

 21    continue

       call lintegrate

       do 71 is=1,nspec
          if(ncspec(is).eq.0) goto 71

          do 70 ie= 1,negrid
             term1= 0.0
             term2= 0.0

             do 60 ig= -ntgridl,ntgrid
                term1= term1 +conjg( phi(ig) )*gint(ig,ie,is)
                term2= term2 +conjg( phi(ig) )*phi(ig)
 60          continue

crwew             term(ie,is)= term1/term2
              term(ie,is)= term1/(1.e-10+term2)

 70       continue

 71    continue


       do 81 is=1,nspec
          if(ncspec(is).eq.0) goto 81

          do 80 ie= 1,negrid
          do 80 ig= -ntgridl,ntgrid
             gint(ig,ie,is)= gint(ig,ie,is)/c
 80       continue

 81    continue


       return
       end subroutine getterm

************************************************************************
      subroutine openfile
************************************************************************

cgms      Use(Common)
      character*20 gstep

      
      do 200 is= 1,nspec
         if(ncspec(is).eq.0) goto 200

         do 100 ie= 1,negrid
            if(nce(ie).eq.0) goto 100

            ichannel= is*10 +ie

            if (ichannel.eq.11) gstep= 'gstep11.dat'
            if (ichannel.eq.12) gstep= 'gstep12.dat'
            if (ichannel.eq.13) gstep= 'gstep13.dat'
            if (ichannel.eq.14) gstep= 'gstep14.dat'
            if (ichannel.eq.15) gstep= 'gstep15.dat'

            if (ichannel.eq.21) gstep= 'gstep21.dat'
            if (ichannel.eq.22) gstep= 'gstep22.dat'
            if (ichannel.eq.23) gstep= 'gstep23.dat'
            if (ichannel.eq.24) gstep= 'gstep24.dat'
            if (ichannel.eq.25) gstep= 'gstep25.dat'

            if (ichannel.eq.31) gstep= 'gstep31.dat'
            if (ichannel.eq.32) gstep= 'gstep32.dat'
            if (ichannel.eq.33) gstep= 'gstep33.dat'
            if (ichannel.eq.34) gstep= 'gstep34.dat'
            if (ichannel.eq.35) gstep= 'gstep35.dat'

            if (ichannel.eq.41) gstep= 'gstep41.dat'
            if (ichannel.eq.42) gstep= 'gstep42.dat'
            if (ichannel.eq.43) gstep= 'gstep43.dat'
            if (ichannel.eq.44) gstep= 'gstep44.dat'
            if (ichannel.eq.45) gstep= 'gstep45.dat'

            if (ichannel.eq.51) gstep= 'gstep51.dat'
            if (ichannel.eq.52) gstep= 'gstep52.dat'
            if (ichannel.eq.53) gstep= 'gstep53.dat'
            if (ichannel.eq.54) gstep= 'gstep54.dat'
            if (ichannel.eq.55) gstep= 'gstep55.dat'

            open (ichannel,file=gstep,status='unknown')
 100     continue

 200  continue


      return
      end subroutine openfile



      
************************************************************************
      subroutine normalmode
************************************************************************
cgms      Use(Basiscom)
cgms      Use(Common)
cgms      Use(Ifwcom)
      
      complex phi00
      
      agammasprev(ikys)=agammas(ikys)
      dgammasprev(ikys)=dgammas(ikys)
      istepk(ikys)=istep
      

      phi00=phi(0)
      
      do j=-ntgridl,ntgrid
       phinorm(j,ikys)=phi(j)/phi00
       aparnorm(j,ikys)=apar(j)/phi00
       bparnorm(j,ikys)=bpar(j)/phi00
      enddo
      
      rnorm=0
      afreqs(ikys)=0.
      agammas(ikys)=0.
      do j=0,ntgrid-1
       rnorm=rnorm+cabs(phinorm(j,ikys))
       afreqs(ikys)=afreqs(ikys)
     >     +freqsj(j,ikys)*cabs(phinorm(j,ikys))
       agammas(ikys)=agammas(ikys)
     >     +gammasj(j,ikys)*cabs(phinorm(j,ikys))
      enddo
      afreqs(ikys)=afreqs(ikys)/rnorm
      agammas(ikys)=agammas(ikys)/rnorm

      nsteps(ikys)=nstep
  
      afreq(ikys)=afreqs(ikys)*temp(3)/kys(ikys)/tprim1
      agamma(ikys)=agammas(ikys)*temp(3)/kys(ikys)/tprim1
      
      dfreqs(ikys)=0.
      dgammas(ikys)=0.
      do j=0,ntgrid-1
       dfreqs(ikys)=dfreqs(ikys)
     >   +abs(freqsj(j,ikys)-afreqs(ikys))*cabs(phinorm(j,ikys))
       dgammas(ikys)=dgammas(ikys)
     >   +abs(gammasj(j,ikys)-agammas(ikys))*cabs(phinorm(j,ikys))
      enddo
      dfreqs(ikys)=dfreqs(ikys)/rnorm
      dgammas(ikys)=dgammas(ikys)/rnorm
  
      dtgammas(ikys)=0.
      rnorm=0.
      do nstepa=nstep/2,nstep
      rnorm=rnorm+1.
       dtgammas(ikys)=dtgammas(ikys)+
     >      (gammast(nstepa,0,ikys)-agammas(ikys))**2
      enddo
      dtgammas(ikys)=(dtgammas(ikys)/rnorm)**0.5
 
      rnorm=0.
      tcenter=0.
      do j=1-ntgridl,ntgrid-1
       rnorm=rnorm+(cabs(phinorm(j,ikys)))**2
       tcenter=tcenter+theta(j)*(cabs(phinorm(j,ikys)))**2
      enddo
       tcenter=tcenter/rnorm
      do j=1-ntgridl,ntgrid-1
       twidths(ikys)=twidths(ikys)
     >   +(theta(j)-tcenter)**2*(cabs(phinorm(j,ikys)))**2 
      enddo
       twidths(ikys)=sqrt(twidths(ikys)/rnorm)
       dmxls(ikys)=agammas(ikys)/(shat*kys(ikys)*twidths(ikys))**2
       dmxlps(ikys)=agammas(ikys)/kys(ikys)**2
     > /(1.+shat**2*twidths(ikys)**2)
       
       
cmnt compute d_eff and chi_eff
cmnt The {theta ave} flows pflx and eflx are normed to 
cmnt    vth(1)*ne0*cabs(phi(0))**2 and temp(is)*vth(1)*ne0*cabs(phi(0))**2
cmnt x phisqnorm = {theta ave} cabs(phi)**2/cabs(phi(0))**2
cmnt phi=(q(1)*phi_phys/temp(1))
cmnt We use the mixing rule from balancing simple EXB nonlinearity against the 
cmnt complex time derivative 
cmnt   phi^=(e*phi/T_e)/(rhos/a)=cabs(gammas+zi*freqs)/(kxs*kys)
cmnt alternatively we use
cmnt   phi^=(e*phi/T_e)/(rhos/a)=(cabs(gammas))**0.5*gamma_damp**0.5/(kxs*kys)

cmnt d_eff=plasma_flux/(-dn/dx); 
cmnt chi_eff=(energy_flux)/(-ndT/dx)
cmnt d_eff and chi_eff are left in units cs/a*rhos**2
cmnt d_effc,diff_effs,d_effw correspond to
cmnt (c)  (agammas)**0.5*gamma_damp**0.5
cmnt (w)eak (agammas), and  (s)trong (afreqs)  saturation rules

      akys=kys(ikys)
      akxs=shat*kys(ikys)*twidths(ikys)
      akxs=sqrt(akys**2+akxs**2)

    
      qq_loc=2.*epsa/pk
      agamma_damps=sqrt(temp(1)/temp(3))/(qq_loc/epsa)*3./2.*
     > sqrt(pi/2.)/(2.+1/qq_loc**2)

c better formula is (temp(1)/temp(3))*(akys)*epsa using kxs=kys

      agamma_damps=(temp(1)/temp(3))*akys*epsa
  
      agamma_drive=abs(agammas(ikys))-dgammas(ikys)
      if(agamma_drive.lt.0.) agamma_drive=0.


      is=1 
      
       pflx=pflxea(is)+pflxma(is)
       qflx=eflxea(is)+eflxma(is)
       
       di_effc(ikys)=pflx*(temp(3)/temp(1))**2*
     > ((agamma_drive)**(0.5-adamp)*agamma_damps**(0.5+adamp)/
     > (akys*akxs))**2/(an(is)*fprim(is)*sqrt(temp(3)/temp(1)/2.))
       di_effw(ikys)=pflx*(temp(3)/temp(1))**2*
     > (abs(agammas(ikys))/
     > (akys*akxs))**2/(an(is)*fprim(is)*sqrt(temp(3)/temp(1)/2.))
       di_effs(ikys)=pflx*(temp(3)/temp(1))**2*
     > (abs(afreqs(ikys))/
     > (akys*akxs))**2/(an(is)*fprim(is)*sqrt(temp(3)/temp(1)/2.))
     
       chii_effc(ikys)=qflx*(temp(3)/temp(1))**2*
     > ((agamma_drive)**(0.5-adamp)*agamma_damps**(0.5+adamp)/
     > (akys*akxs))**2/(an(is)*tprim(is)*sqrt(temp(3)/temp(1)/2.))
       chii_effw(ikys)=qflx*(temp(3)/temp(1))**2*
     > (abs(agammas(ikys))/
     > (akys*akxs))**2/(an(is)*tprim(is)*sqrt(temp(3)/temp(1)/2.))
       chii_effs(ikys)=qflx*(temp(3)/temp(1))**2*
     > (abs(afreqs(ikys))/
     > (akys*akxs))**2/(an(is)*tprim(is)*sqrt(temp(3)/temp(1)/2.))
     
      if (an(3).ne.0.) then
      is=3 
      
       pflx=pflxea(is)+pflxma(is)
       qflx=eflxea(is)+eflxma(is)
       
       de_effc(ikys)=pflx*(temp(3)/temp(1))**2*
     > ((agamma_drive)**(0.5-adamp)*agamma_damps**(0.5+adamp)/
     > (akys*akxs))**2/(an(is)*fprim(is)*sqrt(temp(3)/temp(1)/2.))
       de_effw(ikys)=pflx*(temp(3)/temp(1))**2*
     > (abs(agammas(ikys))/
     > (akys*akxs))**2/(an(is)*fprim(is)*sqrt(temp(3)/temp(1)/2.))
       de_effs(ikys)=pflx*(temp(3)/temp(1))**2*
     > (abs(afreqs(ikys))/
     > (akys*akxs))**2/(an(is)*fprim(is)*sqrt(temp(3)/temp(1)/2.))
     
       chie_effc(ikys)=qflx*(temp(3)/temp(1))**2*
     > ((agamma_drive)**(0.5-adamp)*agamma_damps**(0.5+adamp)/
     > (akys*akxs))**2/(an(is)*tprim(is)*sqrt(temp(3)/temp(1)/2.))
       chie_effw(ikys)=qflx*(temp(3)/temp(1))**2*
     > (abs(agammas(ikys))/
     > (akys*akxs))**2/(an(is)*tprim(is)*sqrt(temp(3)/temp(1)/2.))
       chie_effs(ikys)=qflx*(temp(3)/temp(1))**2*
     > (abs(afreqs(ikys))/
     > (akys*akxs))**2/(an(is)*tprim(is)*sqrt(temp(3)/temp(1)/2.))
 
      endif


      di_effgbc(ikys)=dgyrobohmnorm*di_effc(ikys)
      di_effgbs(ikys)=dgyrobohmnorm*di_effs(ikys)
      di_effgbw(ikys)=dgyrobohmnorm*di_effw(ikys)

      de_effgbc(ikys)=dgyrobohmnorm*de_effc(ikys)
      de_effgbs(ikys)=dgyrobohmnorm*de_effs(ikys)
      de_effgbw(ikys)=dgyrobohmnorm*de_effw(ikys)

      chii_effgbc(ikys)=dgyrobohmnorm*chii_effc(ikys)
      chii_effgbs(ikys)=dgyrobohmnorm*chii_effs(ikys)
      chii_effgbw(ikys)=dgyrobohmnorm*chii_effw(ikys)

      chie_effgbc(ikys)=dgyrobohmnorm*chie_effc(ikys)
      chie_effgbs(ikys)=dgyrobohmnorm*chie_effs(ikys)
      chie_effgbw(ikys)=dgyrobohmnorm*chie_effw(ikys)


  
      stkafreqs(ikys,iistk)=afreqs(ikys)
      stkagammas(ikys,iistk)=agammas(ikys)
      stknstep(ikys,iistk)=nsteps(ikys) 
      stkdfreqs(ikys,iistk)=dfreqs(ikys)
      stkdgammas(ikys,iistk)=dgammas(ikys)   
      stktwidths(ikys,iistk)=twidths(ikys)
      stkdmxls(ikys,iistk)=dmxls(ikys)
      stkdmxlps(ikys,iistk)=dmxlps(ikys)
      
      stkdi_effc(ikys,iistk)=di_effc(ikys)
      stkdi_effs(ikys,iistk)=di_effs(ikys)
      stkdi_effw(ikys,iistk)=di_effw(ikys)
      
      stkde_effc(ikys,iistk)=de_effc(ikys)
      stkde_effs(ikys,iistk)=de_effs(ikys)
      stkde_effw(ikys,iistk)=de_effw(ikys)
      
      stkchii_effc(ikys,iistk)=chii_effc(ikys)
      stkchii_effs(ikys,iistk)=chii_effs(ikys)
      stkchii_effw(ikys,iistk)=chii_effw(ikys)
      
      stkchie_effc(ikys,iistk)=chie_effc(ikys)
      stkchie_effs(ikys,iistk)=chie_effs(ikys)
      stkchie_effw(ikys,iistk)=chie_effw(ikys)

      
      stkdi_effgbc(ikys,iistk)=dgyrobohmnorm*di_effc(ikys)
      stkdi_effgbs(ikys,iistk)=dgyrobohmnorm*di_effs(ikys)
      stkdi_effgbw(ikys,iistk)=dgyrobohmnorm*di_effw(ikys)

      stkde_effgbc(ikys,iistk)=dgyrobohmnorm*de_effc(ikys)
      stkde_effgbs(ikys,iistk)=dgyrobohmnorm*de_effs(ikys)
      stkde_effgbw(ikys,iistk)=dgyrobohmnorm*de_effw(ikys)

      stkchii_effgbc(ikys,iistk)=dgyrobohmnorm*chii_effc(ikys)
      stkchii_effgbs(ikys,iistk)=dgyrobohmnorm*chii_effs(ikys)
      stkchii_effgbw(ikys,iistk)=dgyrobohmnorm*chii_effw(ikys)

      stkchie_effgbc(ikys,iistk)=dgyrobohmnorm*chie_effc(ikys)
      stkchie_effgbs(ikys,iistk)=dgyrobohmnorm*chie_effs(ikys)
      stkchie_effgbw(ikys,iistk)=dgyrobohmnorm*chie_effw(ikys)
      
       freqdel(ikys)=delt*afreqs(ikys)*temp(3)/kys(ikys)
       gammadel(ikys)=delt*agammas(ikys)*temp(3)/kys(ikys)


  
      do j=-ntgridl,ntgrid
       stkphinorm(j,ikys,iistk)=real(phinorm(j,ikys))
      enddo

      return
  
      end subroutine normalmode
      
************************************************************************
      subroutine chkweight
************************************************************************
cgms      Use(Basiscom)
cgms      Use(Common)
      
      ng2=2*ngauss
      ng2p1=ng2+1

      do i=-ntgridl,ntgrid
       do ie=1,negrid
        do il=1,ng2
	     g1(i,il,ie,1,1)=1.
	     g1(i,il,ie,2,1)=1.
        enddo
        do il=ng2p1,nlambda
         g1(i,il,ie,1,1)=0.
         g1(i,il,ie,2,1)=0.  
        enddo
       enddo
    
       call integrate
      
       weightchkp(i)=geint(i,1)
      enddo
       
      do i=-ntgridl,ntgrid
       do ie=1,negrid
        do il=1,ng2
         g1(i,il,ie,1,1)=0.
         g1(i,il,ie,2,1)=0.
        enddo
        do il=ng2p1,nlambda
         g1(i,il,ie,1,1)=1.
         g1(i,il,ie,2,1)=1.
        enddo
       enddo

       call integrate

       weightchkt(i)=geint(i,1)
       weightchk(i)=weightchkt(i)+weightchkp(i)
      enddo  
      return
      end subroutine chkweight
      
************************************************************************
      subroutine quasilinear
************************************************************************
cgms      Use(Basiscom)
cgms      Use(Common)
      
      complex phi00

      

      phi00=phi(0)
      
      ng2p1=2*ngauss+1
      ng2=2*ngauss
     
      
c***********************************************************************      
c   compute total EXB particle flux pflxe
c  in units of ne0*vth(1)*(cabs(phi00))**2
c*****
      do  is=1,nspec
       if(ncspec(is).ne.0) then

       do  isign=1,2
        do  ie=1,negrid
         do  il=1,nlambda
          do j=-ntgridl,ntgrid
           g1(j,il,ie,isign,is)=g(j,il,ie,isign,is)*aj0(j,il,ie,is)
          enddo
         enddo
        enddo
       enddo
      
       endif
      enddo

      call integrate


      do is=1,nspec
       
       do j=-ntgridl,ntgrid
        pflxe(j,is)=0.
       enddo
       if(ncspec(is).ne.0) then
       do j=-ntgridl,ntgrid
          pflxe(j,is)=
     >     +an(is)*real(conjg(geint(j,is))*(-zi*aky*phi(j)))
     >     /(cabs(phi00))**2
       enddo
 
c   computing density and response functions

       do j=-ntgridl,ntgrid
          dengyronorm(j,is)=an(is)*geint(j,is)/phi00
          dennorm(j,is)=
     > an(is)*(geint(j,is)-z(is)*poldrift(j,is)/temp(is)*phi(j))/phi00
          respfunc(j,is)=dennorm(j,is)*phi00/(-z(is)*phi(j)/temp(is))
       enddo
c   computing fourier transforms and response functions

      thetamaxft=real(ntgrid)/real(ntheta)*2.*pi
      delpark=pi/thetamaxft
      do j=-ntgrid,ntgrid
       park(j)=real(j)*delpark
      enddo
 
      do m=-ntgrid+1,ntgrid
       phf(m)=0.
       denf(m,is)=0.
       do j=-ntgrid+1,-1
        phf(m)=phf(m)
     >   +exp(-zi*real(j)*real(m)*pi/real(ntgrid))*phi(-j)/phi00
        denf(m,is)=denf(m,is)
     >   +exp(-zi*real(j)*real(m)*pi/real(ntgrid))*dennorm(-j,is)
        enddo
       do j=0,ntgrid
        phf(m)=phf(m)
     >   +exp(-zi*real(j)*real(m)*pi/real(ntgrid))*phi(j)/phi00
        denf(m,is)=denf(m,is)
     >   +exp(-zi*real(j)*real(m)*pi/real(ntgrid))*dennorm(j,is)
        enddo

       enddo
  
       phf(-ntgrid)=phf(ntgrid)
       denf(-ntgrid,is)=denf(ntgrid,is)
  
       phf00=phf(0)
       do m=-ntgrid,ntgrid
         phf(m)=phf(m)/phf00
         denf(m,is)=denf(m,is)/phf00
         respff(m,is)=denf(m,is)/(-z(is)*phf(m)/temp(is))
       enddo
 
       
       endif
      enddo
     
     
c***********************************************************************      
c   compute total magnetic flutter particle flux pflxm 
c   in units of ne0*vth(1)*(cabs(phi00))**2
c*****
      do  is=1,nspec
       if(ncspec(is).ne.0) then
       stm=sqrt(temp(is)/amass(is))

       do  isign=1,2
        do  ie=1,negrid
         do  il=1,nlambda
          do j=-ntgridl,ntgrid
           g1(j,il,ie,isign,is)=g(j,il,ie,isign,is)*aj0(j,il,ie,is)
     >       *stm*vpa(j,il,ie,isign)
           enddo
          enddo
        enddo
       enddo
      
       endif
      enddo

      call integrate


      do is=1,nspec
       
       do j=-ntgridl,ntgrid
        pflxm(j,is)=0.
       enddo
       if(ncspec(is).ne.0) then
       do j=-ntgridl,ntgrid
          pflxm(j,is)=
     >     +an(is)*real(conjg(geint(j,is))*(zi*aky*apar(j)))
     >     /(cabs(phi00))**2
       enddo
       endif
      enddo
     
c***********************************************************************      
c   compute  EXB trapped particle flux pflxte 
c  in units of ne0*vth(1)*(cabs(phi00))**2
c*****
      do  is=1,nspec
       if(ncspec(is).ne.0) then

       do  isign=1,2
        do  ie=1,negrid
         do  il=ng2p1,nlambda
          do j=-ntgridl,ntgrid
           g1(j,il,ie,isign,is)=g(j,il,ie,isign,is)*aj0(j,il,ie,is)
	      enddo
	     enddo
	     do  il=1,ng2
          do j=-ntgridl,ntgrid
           g1(j,il,ie,isign,is)=0.0
          enddo
         enddo
        enddo
       enddo
      
       endif
      enddo

      call integrate


      do is=1,nspec
       
       do j=-ntgridl,ntgrid
        pflxte(j,is)=0.
       enddo
       if(ncspec(is).ne.0) then
       do j=-ntgridl,ntgrid
         pflxte(j,is)=
     >     +an(is)*real(conjg(geint(j,is))*(-zi*aky*phi(j)))
     >     /(cabs(phi00))**2
       enddo
       endif
      enddo
     
c***********************************************************************      
c   compute total magnetic flutter trapped particle flux pflxtm
c   in units of ne0*vth(1)*(cabs(phi00))**2
c*****
      do  is=1,nspec
       if(ncspec(is).ne.0) then
       stm=sqrt(temp(is)/amass(is))

       do  isign=1,2
        do  ie=1,negrid
         do  il=ng2p1,nlambda
          do j=-ntgridl,ntgrid
           g1(j,il,ie,isign,is)=g(j,il,ie,isign,is)*aj0(j,il,ie,is)
     >       *stm*vpa(j,il,ie,isign)
	      enddo
	     enddo
	     do  il=1,ng2
          do j=-ntgridl,ntgrid
           g1(j,il,ie,isign,is)=0.0
          enddo
         enddo
        enddo
       enddo
      
       endif
      enddo

      call integrate


      do is=1,nspec

       do j=-ntgridl,ntgrid
        pflxtm(j,is)=0.
       enddo
       if(ncspec(is).ne.0) then
       do j=-ntgridl,ntgrid
          pflxtm(j,is)=
     >     +an(is)*real(conjg(geint(j,is))*(zi*aky*apar(j)))
     >     /(cabs(phi00))**2
       enddo
       endif
      enddo

c***********************************************************************      
c   compute total EXB energy flux eflxe 
c  in units of ne0*vth(1)*temp(is)*(cabs(phi00))**2
c*****
      do  is=1,nspec
       if(ncspec(is).ne.0) then

       do  isign=1,2
        do  ie=1,negrid
         do  il=1,nlambda
          do j=-ntgridl,ntgrid
           g1(j,il,ie,isign,is)=
     >       e(ie)*g(j,il,ie,isign,is)*aj0(j,il,ie,is)
          enddo
         enddo
        enddo
       enddo
      
       endif
      enddo

      call integrate


      do is=1,nspec
     
       do j=-ntgridl,ntgrid
        eflxe(j,is)=0.
       enddo
       if(ncspec(is).ne.0) then
       do j=-ntgridl,ntgrid
          eflxe(j,is)=
     >     +an(is)*real(conjg(geint(j,is))*(-zi*aky*phi(j)))
     >     /(cabs(phi00))**2
       enddo
       endif
      enddo
      
c***********************************************************************     
c   compute total magnetic flutter energy flux 
c   in units of ne0*vth(1)*temp(is)*(cabs(phi00))**2
c*****
      do  is=1,nspec
       if(ncspec(is).ne.0) then
       stm=sqrt(temp(is)/amass(is))

       do  isign=1,2
        do  ie=1,negrid
         do  il=1,nlambda
          do j=-ntgridl,ntgrid
           g1(j,il,ie,isign,is)=g(j,il,ie,isign,is)*aj0(j,il,ie,is)
     >       *e(ie)*stm*vpa(j,il,ie,isign)
          enddo
         enddo
        enddo
       enddo
      
       endif
      enddo

      call integrate


      do is=1,nspec
       
       do j=-ntgridl,ntgrid
        eflxm(j,is)=0.
       enddo
       if(ncspec(is).ne.0) then
       do j=-ntgridl,ntgrid
        eflxm(j,is)=
     >     +an(is)*real(conjg(geint(j,is))*(zi*aky*apar(j)))
     >     /(cabs(phi00))**2
       enddo
       endif
      enddo

c***********************************************************************     
c   compute  EXB trapped energy flux 
c  in units of ne0*vth(1)*temp(is)*(cabs(phi00))**2
c*****
      do  is=1,nspec
       if(ncspec(is).ne.0) then

       do  isign=1,2
        do  ie=1,negrid
         do  il=ng2p1,nlambda
          do j=-ntgridl,ntgrid
           g1(j,il,ie,isign,is)=
     >       e(ie)*g(j,il,ie,isign,is)*aj0(j,il,ie,is)
	      enddo
	     enddo
	     do  il=1,ng2
          do j=-ntgridl,ntgrid
           g1(j,il,ie,isign,is)=0.0
          enddo
         enddo
        enddo
       enddo
      
       endif
      enddo

      call integrate


      do is=1,nspec
    
       do j=-ntgridl,ntgrid
        eflxte(j,is)=0.
       enddo
       if(ncspec(is).ne.0) then
       do j=-ntgridl,ntgrid
          eflxte(j,is)=
     >     +an(is)*real(conjg(geint(j,is))*(-zi*aky*phi(j)))
     >     /(cabs(phi00))**2
       enddo
       endif
      enddo

c***********************************************************************     
c   compute trapped magnetic flutter energy flux 
c   in units of ne0*vth(1)*temp(is)*(cabs(phi00))**2
c*****
      do  is=1,nspec
       if(ncspec(is).ne.0) then
       stm=sqrt(temp(is)/amass(is))
 
       do  isign=1,2
        do  ie=1,negrid
         do  il=ng2p1,nlambda
          do j=-ntgridl,ntgrid
           g1(j,il,ie,isign,is)=g(j,il,ie,isign,is)*aj0(j,il,ie,is)
     >       *e(ie)*stm*vpa(j,il,ie,isign)
	      enddo
	     enddo
	     do  il=1,ng2
          do j=-ntgridl,ntgrid
           g1(j,il,ie,isign,is)=0.0
          enddo
         enddo
        enddo
       enddo
      
       endif
      enddo

      call integrate


      do is=1,nspec
       
       do j=-ntgridl,ntgrid
        eflxtm(j,is)=0.
       enddo
       if(ncspec(is).ne.0) then
       do j=-ntgridl,ntgrid
          eflxtm(j,is)=
     >     +an(is)*real(conjg(geint(j,is))*(zi*aky*apar(j)))
     >     /(cabs(phi00))**2
       enddo
       endif
      enddo
     
c***********************************************************************      
c   compute total electrostatic momentum flux mflxe (parallel stress) 
c   in units of mi0*ne0*vth(1)*vth(1)*(cabs(phi00))**2
c*****
      do  is=1,nspec
       if(ncspec(is).ne.0) then
       stm=sqrt(temp(is)/amass(is))

       do  isign=1,2
        do  ie=1,negrid
         do  il=1,nlambda
          do j=-ntgridl,ntgrid
           g1(j,il,ie,isign,is)=g(j,il,ie,isign,is)*aj0(j,il,ie,is)
     >       *stm*vpa(j,il,ie,isign)
           enddo
          enddo
        enddo
       enddo
      
       endif
      enddo

      call integrate


      do is=1,nspec
       
       do j=-ntgridl,ntgrid
        mflxe(j,is)=0.
       enddo
       if(ncspec(is).ne.0) then
! note negative sign in order to make this positive for vprim > 0
       do j=-ntgridl,ntgrid
          mflxe(j,is)=
     >     -amass(is)*an(is)*real(conjg(geint(j,is))*(-zi*aky*phi(j)))
     >     /(cabs(phi00))**2
       enddo
       endif
      enddo
c______________________________________________________________________
c  compute the flux surface average flow ie integrate over all theta
c*****
        normpi=(ntgridl+ntgrid)/ntgrid 
        anormpi=normpi
        phisqrnorm=0.
      do j=-ntgridl,ntgrid
        phisqrnorm=phisqrnorm+
     >   delthet(j)/(anormpi*pi)*conjg(phi(j))*phi(j)/(cabs(phi00))**2
      enddo

      do is=1,nspec
       
        pflxea(is)=0.
        pflxma(is)=0.
        pflxtea(is)=0.
        pflxtma(is)=0.
	
	eflxea(is)=0.
        eflxma(is)=0.
        eflxtea(is)=0.
        eflxtma(is)=0.
 
        mflxea(is) = 0.0
	
      enddo
      
      do is=1,nspec
       do j=-ntgridl,ntgrid 
        pflxea(is)=pflxea(is)+delthet(j)*pflxe(j,is)/(anormpi*pi)
     >   /phisqrnorm
        pflxma(is)=pflxma(is)+delthet(j)*pflxm(j,is)/(anormpi*pi)
     >   /phisqrnorm
        pflxtea(is)=pflxtea(is)+delthet(j)*pflxte(j,is)/(anormpi*pi)
     >   /phisqrnorm
        pflxtma(is)=pflxtma(is)+delthet(j)*pflxtm(j,is)/(anormpi*pi)
     >   /phisqrnorm
	
        eflxea(is)=eflxea(is)+delthet(j)*eflxe(j,is)/(anormpi*pi)
     >   /phisqrnorm
        eflxma(is)=eflxma(is)+delthet(j)*eflxm(j,is)/(anormpi*pi)
     >   /phisqrnorm
        eflxtea(is)=eflxtea(is)+delthet(j)*eflxte(j,is)/(anormpi*pi)
     >   /phisqrnorm
        eflxtma(is)=eflxtma(is)+delthet(j)*eflxtm(j,is)/(anormpi*pi)
     >   /phisqrnorm

        mflxea(is)=mflxea(is)+delthet(j)*mflxe(j,is)/(anormpi*pi)
     >   /phisqrnorm
       enddo
      enddo
      
      do j=-ntgridl,ntgrid
       peflx(j)=0.
       eeflx(j)=0.
      enddo
      is=3
      do j=-ntgridl,ntgrid
       peflx(j)=peflx(j)+pflxe(j,is)+pflxm(j,is)
       eeflx(j)=eeflx(j)+eflxe(j,is)+eflxm(j,is)
      enddo

      do j=-ntgridl,ntgrid
       piflx(j)=0.
       eiflx(j)=0.
       miflx(j)=0.
      enddo
      do is=1,nspec
       if(is.ne.3) then
       do j=-ntgridl,ntgrid
        piflx(j)=piflx(j)+pflxe(j,is)+pflxm(j,is)
        eiflx(j)=eiflx(j)+eflxe(j,is)+eflxm(j,is)
        miflx(j)=miflx(j)+mflxe(j,is)
       enddo
       endif
      enddo
 
       peflxa=pflxea(3)+pflxma(3)
       eeflxa=eflxea(3)+eflxma(3)
      
      piflxa=0.
      eiflxa=0.
      miflxa=0.
      do is=1,nspec
       if(is.ne.3) then
        piflxa=piflxa+pflxea(is)+pflxma(is)
        eiflxa=eiflxa+eflxea(is)+eflxma(is)
        miflxa=miflxa+mflxea(is)
       endif
      enddo
       
       

      return
      end subroutine quasilinear
************************************************************************
      subroutine  current
************************************************************************
cgms      Use(Basiscom)
cgms      Use(Common)

      complex phi00
      complex phf00


      phi00=phi(0)

      ng2p1=2*ngauss+1
      ng2=2*ngauss


c***********************************************************************
c compute normalized trapped electron current
c*****
       is=3
       stm=sqrt(temp(is)/amass(is))

       do  isign=1,2
        do  ie=1,negrid
         do  il=ng2p1,nlambda
          do j=-ntgridl,ntgrid
           g1(j,il,ie,isign,is)=g(j,il,ie,isign,is)*aj0(j,il,ie,is)
     >       *stm*vpa(j,il,ie,isign)
          enddo
         enddo
         do  il=1,ng2
          do j=-ntgridl,ntgrid
           g1(j,il,ie,isign,is)=0.0
          enddo
         enddo
        enddo
      enddo

      call integrate



       is=3
       do j=-ntgridl,ntgrid
        curtrap(istep,j)=
     >     +an(is)*geint(j,is)/cabs(phi00)
       enddo
c***********************************************************************
c compute normalized passing electron current
c*****
       is=3
       stm=sqrt(temp(is)/amass(is))

       do  isign=1,2
        do  ie=1,negrid
         do  il=ng2p1,nlambda
          do j=-ntgridl,ntgrid
           g1(j,il,ie,isign,is)=0.0
          enddo
         enddo
         do  il=1,ng2
          do j=-ntgridl,ntgrid
           g1(j,il,ie,isign,is)=g(j,il,ie,isign,is)*aj0(j,il,ie,is)
     >       *stm*vpa(j,il,ie,isign)   
          enddo
         enddo
        enddo
      enddo

      call integrate



       is=3
       do j=-ntgridl,ntgrid
        curpass(istep,j)=
     >     +an(is)*geint(j,is)/cabs(phi00)
       enddo

c***********************************************************************
c  compute normalized wave
c*****  
              
       wave(istep)=cwave*phi00/cabs(phi00)
               
       return
       end subroutine  current
************************************************************************
      subroutine  colrate
************************************************************************
cgms      Use(Basiscom)
cgms      Use(Common)
 
c returns g1 as the collision rate relative to vnew_eff



      integer ivb(0:ntm)
      complex entrap(0:ntm),enpass(0:ntm)

      call chkweight 

      ng2=2*ngauss
      ng2p1=ng2+1

c*****
      do isig=1,2
       do ie=1,negrid
        do j=1,nlambda
         do i=-ntgridl,ntgrid
          colratedet(i,j,ie,isig)=0.
          colratedetm(i,j,ie,isig)=0.
          colratedetcm(i,j,ie,isig)=0.
          do is=1,nspec
           g1(i,j,ie,isig,is)=0.
           g2(i,j,ie,isig,is)=0.
          enddo
         enddo
        enddo
       enddo
      enddo

      is=3
c compute trapped electron density

      do isig=1,2
       do ie=1,negrid
        do j=ng2p1,nlambda
         do i=-ntgridl,ntgrid
          g1(i,j,ie,isig,is)=
     >       g(i,j,ie,isig,is)*aj0(i,j,ie,is)
         enddo
        enddo
       enddo
      enddo

      do isig=1,2
       do ie=1,negrid
        do j=1,ng2
         do i=-ntgridl,ntgrid
          g1(i,j,ie,isig,is)=0.
         enddo
        enddo
       enddo
      enddo


      call integrate

      do i=-ntgridl,ntgrid
       entrap(i)=geint(i,is)
      enddo

c compute passing electron  density

      do isig=1,2
       do ie=1,negrid
        do j=ng2p1,nlambda
         do i=-ntgridl,ntgrid
          g1(i,j,ie,isig,is)=0
         enddo
        enddo
       enddo
      enddo

      do isig=1,2
       do ie=1,negrid
        do j=1,ng2
         do i=-ntgridl,ntgrid
          g1(i,j,ie,isig,is)=
     >       g(i,j,ie,isig,is)*aj0(i,j,ie,is)
         enddo
        enddo
       enddo
      enddo


      call integrate

      do i=-ntgridl,ntgrid
       enpass(i)=geint(i,is)
      enddo
c set up pitch angle scattering operator 
           ivb(0)=ntheta/2

         do 2 i=1,ntheta/2+1
 2         ivb(i)=ntheta/2+2-i

         do 4 i=ntheta/2+2,ntheta
 4         ivb(i)=i-ntheta/2

         do 5 i=0,ntgrid
           ipointer=mod(i+1,ntheta)
           jend(i)=ng2+ivb(ipointer)
 5       continue
 
         do 6 i=-ntgridl,0
	   jend(i)=jend(-i)
 6       continue



cgms         alpp=4./(3.*sqrt(3.14159265))
         pi=atan2(0.0,-1.0)
         alpp=4./(3.*sqrt(pi))

      do ie=1,negrid
       vnew(ie,is)=vnew0(is)/e(ie)**1.5*( zeff+alpp*e(ie)/
     &          sqrt(1.+(alpp*e(ie))**2) )
      enddo

           do 40 ie=1,negrid
           do 40 i=-ntgridl,ntgrid
             jendm1=jend(i)-1

             slbnorm(1,i)=1.
             slbnorm(jend(i),i)=0.
             do 20 j=2,jendm1
                slb0= sqrt( abs(1.-bmag(i)*al(j-1)) )
                slb1= sqrt( abs(1.-bmag(i)*al(j  )) )
                slb2= sqrt( abs(1.-bmag(i)*al(j+1)) )
 
                slbnorm(j,i)=slb0

                slbl= (slb1+slb0)/2.
                slbr= (slb1+slb2)/2.


                cc(i,j,ie,is)=vnew(ie,is)*(1.-slbr*slbr)
     &                        /(slbr-slbl)/(slb2-slb1)
                aa(i,j,ie,is)=vnew(ie,is) *(1.-slbl*slbl)
     &                        /(slbr-slbl)/(slb1-slb0)

c norm to vnew_eff=vnew0(is)*zeff/eps
             
                cc(i,j,ie,is)=cc(i,j,ie,is)*eps/vnew0(is)/zeff
                aa(i,j,ie,is)=aa(i,j,ie,is)*eps/vnew0(is)/zeff

                bb(i,j,ie,is)=-(aa(i,j,ie,is)+cc(i,j,ie,is))
 
      colratedet(i,j,ie,1)=
     >               cc(i,j,ie,is)*g(i,j+1,ie,1,is)/g(i,j,ie,1,is)
     >              +bb(i,j,ie,is)
     >              +aa(i,j,ie,is)*g(i,j-1,ie,1,is)/g(i,j,ie,1,is)
 
      colratedet(i,j,ie,2)=
     >               cc(i,j,ie,is)*g(i,j+1,ie,2,is)/g(i,j,ie,2,is)
     >              +bb(i,j,ie,is)
     >              +aa(i,j,ie,is)*g(i,j-1,ie,2,is)/g(i,j,ie,2,is)  

      g2(i,j,ie,1,is)=colratedet(i,j,ie,1)*g(i,j,ie,1,is)
      g2(i,j,ie,2,is)=colratedet(i,j,ie,2)*g(i,j,ie,2,is)

      colratedetm(i,j,ie,1)=-e(ie)**(-1.5)*(1.-phi(i)/g(i,j,ie,1,is))
      colratedetm(i,j,ie,2)=-e(ie)**(-1.5)*(1.-phi(i)/g(i,j,ie,2,is))
 
      colratedetcm(i,j,ie,1)=
     >           -(1.-(entrap(i)+enpass(i))/g(i,j,ie,1,is))
      colratedetcm(i,j,ie,2)=
     >           -(1.-(entrap(i)+enpass(i))/g(i,j,ie,2,is))

 
 20          continue


 40        continue

c compute trapped electron collisions 

      do isig=1,2
       do ie=1,negrid
        do j=ng2p1,nlambda
         do i=-ntgridl,ntgrid
          g1(i,j,ie,isig,is)=
     >colratedet(i,j,ie,isig)*g(i,j,ie,isig,is)*aj0(i,j,ie,is)
         enddo
        enddo
       enddo
      enddo

      do isig=1,2
       do ie=1,negrid
        do j=1,ng2
         do i=-ntgridl,ntgrid
          g1(i,j,ie,isig,is)=0.
         enddo
        enddo
       enddo
      enddo    


      call integrate

      do i=-ntgridl,ntgrid
       colratetrap(i)=geint(i,is)
      enddo

c compute passing electron collisions 

      do isig=1,2
       do ie=1,negrid
        do j=ng2p1,nlambda
         do i=-ntgridl,ntgrid
          g1(i,j,ie,isig,is)=0
         enddo
        enddo
       enddo
      enddo

      do isig=1,2
       do ie=1,negrid
        do j=1,ng2
         do i=-ntgridl,ntgrid
          g1(i,j,ie,isig,is)=
     >colratedet(i,j,ie,isig)*g(i,j,ie,isig,is)*aj0(i,j,ie,is)     
         enddo
        enddo
       enddo
      enddo


      call integrate

      do i=-ntgridl,ntgrid
       colratepass(i)=geint(i,is)
      enddo    


c compute total electron collisions

      do isig=1,2
       do ie=1,negrid
        do j=1,nlambda
         do i=-ntgridl,ntgrid
          g1(i,j,ie,isig,is)=
     >colratedet(i,j,ie,isig)*g(i,j,ie,isig,is)*aj0(i,j,ie,is)
         enddo
        enddo
       enddo
      enddo


      call integrate

      do i=-ntgridl,ntgrid
       colratetot(i)=geint(i,is)
      enddo


c norm colratetrap and colratepass to 
c   vnew_pass*(enpass_0*entrap-entrap_0*enpass)

c*****
      do i=-ntgridl,ntgrid
       colratetrapn(i)=(colratetrap(i)+1.e-33)/
     > (weightchkp(i)*entrap(i)-weightchkt(i)*enpass(i)+1.e-33)
       colratepassn(i)=(colratepass(i)+1.e-33)/
     > (weightchkp(i)*entrap(i)-weightchkt(i)*enpass(i)+1.e-33)
       colratetotn(i)=(colratetot(i)+1.e-33)/
     > (weightchkp(i)*entrap(i)-weightchkt(i)*enpass(i)+1.e-33)
      enddo

       ntheth=ntheta/2
      
       colratetrapn(ntheth)=0.
       colratetrapn(ntheth+ntheta)=0.
       colratepassn(ntheth)=0.
       colratepassn(ntheth+ntheta)=0.
       colratetotn(ntheth)=0.
       colratetotn(ntheth+ntheta)=0.
  
c compute model trapped electron collisions

      do isig=1,2
       do ie=1,negrid
        do j=ng2p1,nlambda
         do i=-ntgridl,ntgrid
          g1(i,j,ie,isig,is)=
     >colratedetm(i,j,ie,isig)*g(i,j,ie,isig,is)*aj0(i,j,ie,is)
         enddo
        enddo
       enddo
      enddo

      do isig=1,2
       do ie=1,negrid
        do j=1,ng2
         do i=-ntgridl,ntgrid
          g1(i,j,ie,isig,is)=0.
         enddo
        enddo
       enddo
      enddo


      call integrate

      do i=-ntgridl,ntgrid
       colratetrapm(i)=geint(i,is)
      enddo

c compute model passing electron collisions

      do isig=1,2
       do ie=1,negrid
        do j=ng2p1,nlambda
         do i=-ntgridl,ntgrid
          g1(i,j,ie,isig,is)=0
         enddo
        enddo
       enddo
      enddo

      do isig=1,2
       do ie=1,negrid
        do j=1,ng2
         do i=-ntgridl,ntgrid
          g1(i,j,ie,isig,is)=
     >colratedetm(i,j,ie,isig)*g(i,j,ie,isig,is)*aj0(i,j,ie,is)
         enddo
        enddo
       enddo
      enddo


      call integrate

      do i=-ntgridl,ntgrid
       colratepassm(i)=geint(i,is)
      enddo

c norm to model

      do i=-ntgridl,ntgrid
       colratetrapnm(i)=
     >   (colratetrap(i)+1.e-33)/(colratetrapm(i)+1.e-33)
       colratepassnm(i)=
     >   (colratepass(i)+1.e-33)/(colratepassm(i)+1.e-33)
      enddo
  
       ntheth=ntheta/2

       colratetrapnm(ntheth)=0.
       colratetrapnm(ntheth+ntheta)=0.
       colratepassnm(ntheth)=0.
       colratepassnm(ntheth+ntheta)=0.
  
c compute conservative model trapped electron collisions

      do isig=1,2
       do ie=1,negrid
        do j=ng2p1,nlambda
         do i=-ntgridl,ntgrid
          g1(i,j,ie,isig,is)=
     >colratedetcm(i,j,ie,isig)*g(i,j,ie,isig,is)*aj0(i,j,ie,is)
         enddo
        enddo
       enddo
      enddo

      do isig=1,2
       do ie=1,negrid
        do j=1,ng2
         do i=-ntgridl,ntgrid
          g1(i,j,ie,isig,is)=0.
         enddo
        enddo
       enddo
      enddo


      call integrate

      do i=-ntgridl,ntgrid
       colratetrapcm(i)=geint(i,is)
      enddo

c compute conservative model passing electron collisions

      do isig=1,2
       do ie=1,negrid
        do j=ng2p1,nlambda
         do i=-ntgridl,ntgrid
          g1(i,j,ie,isig,is)=0
         enddo
        enddo
       enddo
      enddo

      do isig=1,2
       do ie=1,negrid
        do j=1,ng2
         do i=-ntgridl,ntgrid
          g1(i,j,ie,isig,is)=
     >colratedetcm(i,j,ie,isig)*g(i,j,ie,isig,is)*aj0(i,j,ie,is)
         enddo
        enddo
       enddo
      enddo


      call integrate

      do i=-ntgridl,ntgrid
       colratepasscm(i)=geint(i,is)
      enddo

c norm to model

      do i=-ntgridl,ntgrid
       colratetrapncm(i)=
     >   (colratetrap(i)+1.e-33)/(colratetrapcm(i)+1.e-33)
       colratepassncm(i)=
     >   (colratepass(i)+1.e-33)/(colratepasscm(i)+1.e-33)
      enddo

       ntheth=ntheta/2

       colratetrapnm(ntheth)=0.
       colratetrapnm(ntheth+ntheta)=0.
       colratepassnm(ntheth)=0.
       colratepassnm(ntheth+ntheta)=0.
 
      return
      end subroutine  colrate



************************************************************************
       subroutine ifwrite
************************************************************************
cgms      Use(Ifwcom)
cgms      Use(Basiscom)
cgms      Use(Common)
  
       real lfast
       real lnlamda,taue,mui,vthi,xnuei

c      converts physical plasma machine variables to "gypinp" inputs
c      bt,zeff,zeff5
c      ne,nfast,nfastp,te,ti,q,lne,lte,lti,ls
c      rmin(a),rsurf(r),rmaj
c      read from ifwin input file
c
c      rmin="a"
c      rsurf="r"
c      unit of length: meter
c      if set rmin=1 unit of length is "a"
c         must then have input lengths in units of "a"
c      unit density: 10**13 cm**-3
c      unit temperature: keV
c      unit bt: Tesla
c
       bakdif1=0.0
       bakdif2=0.0
       bakdif3=0.0
       bakdif4=0.0
       bakdif5=0.0

crew6/30/97       z2=6.
crew       z3=-1.
       z4=1.
       z5=27.

crew6/30/97       amass2=6.
crew       amass3=2.7e-4
       amass4=1.0
       amass5=28.6

       epsa=rmin/rmaj
       eps=rsurf/rmaj
c       ntheta=32
c       nperiod=2
c       ngauss=5

cgms       an5=zeff5/(z5**2)
       an5 = 0.0
cgms       an2=(zeff-1.-z5*(z5-1.)*an5)/((z2-1.)*z2)
       an2=nz/ne
       anlt=1.-an2*z2-an5*z5
       dbeam=nfast/ne
       an4=dbeam
cgms       an1=anlt-an4
       an1=ni/ne
       anb=dbeam
       if (anlt.lt.0.) write(6,*) "anlt.lt.0"

       shat=rmin*q/(ls*epsa)
       lfast=-nfast/(nfastp+1.e-10)
      
       if (lfast.eq.0.)  lfast=1.e10

cgms       fprim1=((rmin/lne)*(1.-z2*an2-z5*an5)-anb*rmin/lfast)/an1
       fprim1=rmin/lni

crew       shift=q**2*rmaj*4.e5*ne*(te*(1./lne+1./lte)+an1*ti*
crew     &  (fprim1/rmin+1./lti)+anb*50.0/lfast+an2*ti*(1./lne+1./lti))
crew     &       /(shat*(bt*10**4)**2)
  
crew  BEWARE shift comes from gksin as an input

       pk=2.*epsa/q
       tprim1=rmin/lti
       epsl=2.*rmin/rmaj
       icv=1
       width0=3.
       fprim3=(rmin/lne)/fprim1
       fprim2=(rmin/lnz)/fprim1
       tprim2=1.
       tprim3=rmin/lte
       tprim4=0.0
       fprim4=rmin/lfast
       tprim5=1.
       fprim5=(rmin/lne)/fprim1
       temp2=1.
       temp3=te/ti
       temp4=50./ti
       temp5=temp2
       teti=ti/te
c
crew       beta=403.*ne*te/(1.e5*bt**2)
crew CORRECTION
       beta=403.*ne*ti/(1.e5*bt**2)
       debyelorhos=cdebye*
     > 7.43e2/1.02e2*(bt*1.e4)/(ne*1.e13)**.5*(amass3/.0005446)**.5
c cdebyelorhos defaults to zero on input 
c      write(*,*)"debyelorhos=",debyelorhos,"bt=",bt,"ne=",ne
       vnewk4=0.0
crew   vnewk3=0.051*ne*te**(-1.5)/(ti**0.5)
crew   this close to correct for hardwire rmin=.8 deuterium with
crew   no 1/2 factor in solfp collision operator
crew  CORRECTION (used coulomb log=15.)
crew   vnewk3=0.063*ne*te**(-1.5)/(ti**0.5)*rmin*(5.4e-4/amass3)**0.5
crew forgot to multiply by 3/4*pi**.5 see HH 4.36
       vnewk3=0.117*ne*te**(-1.5)/(ti**0.5)*rmin*(2.723014D-4/amass3)**0.5
cgms   need to divide this by 2 to get the right units
cgms   vnewk3 = nu_ei*a/vth_i where vth_i=sqrt(2*ti/mi)
       vnewk3=vnewk3/2.
cgms new calculation of vnewk3 with lnlamda factor 
cgms lnlamda and taue from NRL formulary
cgms note: for Te=1Kev, ne=10**13 lnlamda = 15.94 and taue=1.088D-3/lnlamda
       lnlamda = 15.94D0-0.5*LOG(ne)+LOG(te)
       taue = 1.088D-3*(te**1.5)/(ne*lnlamda)
cgms xnuei = 3/4 sqrt(pi)/taue
       xnuei = 1.329D0/taue
       mui = (5.446D-4)/amass3
       vthi = 9.79D5*(2.0*ti*1.D3/mui)**0.5
       vnewk3 = xnuei*rmin*100.D0/vthi
       write(*,*)"debug",temp3,mui,vnewk3*(2.0/temp3)**0.5
c
       vnewk3=cnewk3*vnewk3
c
crew this is normed to deuterium. agreement with MK valid
crew  have not traced vnewk3 CORRECTION to vnewk2 and vnewk5 ?
       vnewk2=vnewk3*0.291*(zeff-1.)*(te/ti)**1.5
       vnewk2=cnewk2*vnewk2
       vnewk1=0.0
crew       vnewk5=vnewk2*zeff5*((z5/z2)**2)/(zeff-1)
       vnewk5=vnewk2*zeff5*((z5/z2)**2)/(zeff-.999999999)

crew add vnewstar calc.
       vnewstare=zeff*2.91e-6*ne*1.e13*15./(te*1.e3)**2*rmaj*100.*q
     >     /(rsurf/rmaj)**1.5/4.19e7
       vnewstari=4.78e-8*ne*1.e13*15./(ti*1.e3)**2*rmaj*100.*q
     >     /(rsurf/rmaj)**1.5/9.79e5

c  ifwrite produces the gyinp input physics variables

c        ntheta,nperiod,ngauss,icv

c        eps,shift,dbeam,shat,pk

c        epsl,width0,beta,zeff,teti,zeff5

c        fprim1,fprim2,fprim3,fprim4,fprim5

c        tprim1,tprim2,tprim3,tprim4,tprim5

c        vnewk1,vnewk2,vnewk3,vnewk4,vnewk5

c        bakdif1,bakdif2,bakdif3,bakdif4,bakdif5

c               amass2,amass3,amass4,amass5

c                temp2, temp3, temp4, temp5

c                   z2,    z3,    z4,    z5
  
c     compute gyrobohm norm

      dgyrobohmnorm=(9.79e5*(te*1.e3)**.5)/rmin/100.
     >  *((1.02e2*(te*1.e3)**.5)/bt/1.e4)**2
     >  *(5.446e-4/amass3)**.5
c
       return
       end subroutine ifwrite
************************************************************************
       subroutine ifwritegeo
************************************************************************
cgms      Use(Ifwcom)
cgms      Use(Basiscom)
cgms      Use(Common)
cgms      Use(Geo)
  
       real lfast
       real lnlamda,taue,mui,vthi,xnuei
       integer i_s_delta_loc, i_s_kappa_loc
       
crew_________geo____geo_____geo________________________________________
c     ifwritegeo is analog of ifwrite for real geometry
c     12/4/97 replace brizard's input
c     It is a single flux surface routine 
c     It feeds output directly to GSTOTAL and must be used with igeo=1
c     it replaces the theta shat_alpha modeltheta dependence
c        bmag(j)=1./(1.+eps*cos(theta(j))
c        and curvature drift cos(theta(j))+kx/ky*sin(theta(j))
c        kx/ky = (theta-theta0)shat+alpha*sin(theta)
c     with a real geometry theta dependence in "_geo(j)" quanities 
c    
c     INPUTS:
c
c::::::shape veriables::::::::::
c        aspectratio_loc
c        shift_loc
c        delta_loc
c        kappa_loc
c        s_kappa_loc: degenerate as (kappa_loc-kappa0)/kappa_loc
c        s_delta_loc: degenerate as  delta_loc/(1-delta_loc**2)**0.5
c::::::current profile::::::::::
c        q_loc
c        shat_loc = rmin_loc/q_loc*d q_loc/d rmin_loc     
c::::::pure plasma profile::::::
c        tiote_loc
c        nione_loc   = 1 for pure plasma
c
c        dlntidr_loc   = -(d ti/ d rmin_loc)/ti
c        dlntedr_loc
c        dlnnidr_loc
c        dlnnedr_loc
c
c        beta_loc     
c        xnu_loc
c
c::::::impurity profile::::::
c         
c        zeff_loc        : used in ei-collizionality,
c                           ion dilution and impurity dynamics            
c        dlnnimpdr_loc      : only used in impurity dynamics
c        fastionfrac_loc : used in ion dilution and fast ion dynamics
c             
c:::::overrides:::::::::if supplied with non zero values for
c         ne_loc, te_loc, b00_loc (flux surface center field)
c         xnu_loc and beta_loc will be overridden  and
c         b_unit will be calculated
c         
c         ne_loc in 10**13 cm**-3   te_loc in keV
c     MUST BE CAREFUL that ne_loc and te_loc are 0 otherwise
c
c:::::overrides:::::::::
c         if beta_loc.eq.0 then  calulated from alpha_mhd_loc
c            and shat_loc calculated from shat_mhd_loc
c         if alpha_mhd_loc.eq.0 then calulated from beta_loc 
c            and shat_mhd_loc calulated from shat_loc
c         if beta_loc_0.neq.0 the beta_loc is recalulated
c           from beta_loc_0=p/(bt_mag_center**2/(8.*pi)
c         if dlnpdr_loc   = -(d p/ d rmin_loc)/p.neq.0 where p is the 
c           total MHD pressure, "geometric" MHD beta and gradient of beta
c           are calucalted from p.....only for iptot.gt.0
c         
c         "name"_loc_out is output value of overridden variable
c
c         if b00_loc.neq.0 b_unit is not 1.0 and b_norm = f/rmaj0_loc
c           is b00_loc
c
c         rmin_loc is set to 1 unless overridden
c
c    MHD OUTPUTS:
c        Chu's CAMINO shat_mhd and alpha_mhd are calculated
c
c
c::::::shape:::::::::::
c     This setup assumes updown symmetric triangulated elongated ellipse
c     shaped equibibria from Miller GA-A22365
c     
c     The equiibrium is centered at major radius = "rmaj0_loc" and 
c     minor radius = "rmin_loc"  "rmaj0_loc = aspectratio_loc*rmaj0_loc"
c   
c     The Shafranov shift is (d rmaj0_loc /d rmin_loc) = "shift_loc"
c   
c     The elongation is "kappa_loc" and the triangulrity is "delta_loc"
c
c     The variation of elongation is "s_kappa_loc" = 
c           (rmin_loc/kappa_loc) d kappa_loc / d rmin_loc
c
c     The variation of triabgularity is "s_delta_loc"  = 
c           (rmin_loc/(1-delta_loc**2)**0.5 d delta_loc / d rmin_loc
c
c     A good first model for s_kappa_loc is (kappa_loc-kappa_loc0)/kappa_loc
c        we will use kappa_loc0=1.
c  
c     A good first model for s_delta_loc  is delta_loc/(1-delta_loc**2)**0.5
c
c     The flux surface major radius "rmaj_theta" is given by Miller Eq 34
c     The flux surface height          "z_theta" is given by Miller Eq 35
c
c     The poloidal field is 
c               "bp_theta" = (rmin_loc/rmaj_theta)*b_unit/q_loc*grad_r_theta
c     grad_r_theta is given by combining Miller Eq 36 and Eq37
c     
c     The toroidal field is
c               "bt_theta" = f/rmaj_theta = (rmin_loc/rmaj_theta)*b_unit/
c     (integral_dtheta/2*pi*d length_1/d theta/
c                       ((rmaj_theta/rmin_loc)*grad_r_theta))
c      
c      length_1 is poloidal length at rmin_loc=1.
c      dlength_1/dtheta=
c           (1/rmin_loc)*(d rmaj_theta/d theta)**2+d z_theta/d theta)**2)**0.5
c      
c      theta is not physical poloidal angle but 0=outside pi/2=top  pi=inside
c    
c::::::::units::::::::::::::::
c     
c      a_unit is the arbitary unit of length. 
c           a_units is usually thought of as the plasma minor
c           radius "a". In gstotal this is "rmin"
c           We think of rmin_loc as the variable "rsurf" 
c             i.e.eps=rsurf/rmaj0
c             but one can think of a_unit as rmin_loc also
c           rmin_loc and rmaj_loc are in units of a_unit
c           The gradients (gradti_loc, gradne_loc, etc) are in units of 1/a_unit
c           All rates are in unnit of cs_loc/a_init where 
c             cs_loc=sqrt(te_loc/mass_ion)
c 
c      b_unit is the "unit field"  = b0*rho_loc*d rho_loc/(rmin_loc*d rmin_loc)
c          where  toroidal flux is chi_loc=b0*rho_loc**2/2 and b0 is field on 
c          on magnetic axis without plasma
c    
c          b_unit obtained in terms toroidal field on flux surface center 
c          from b00_loc= f/rmaj0_loc
c          
c          The growth rates (gamma) do not depend on b_unit since b_unit comes 
c          in only in  mode label kyrhos_loc= n*q_loc/rmin_loc*rhos_unit_loc
c          where rhos_unit_loc= cs_loc/(e*b_unit/c*mass_ion)           
c          however the diffusion gamma/kperp**2 is in gyrobohm units of
c             cs_loc/a*rhos_unit_loc**2  proportional to b_unit**2
c
c          Nonlinear simulations are need to assess confinement versus shape 
c          at FIXED MINOR RADIUS radius, we can guess n*T*tau proportional to 
c          beta*<bt_theta**2>/(rhos_unit_loc**2*cs_loc/a_unit)*rmin_loc**2
c            /{max_gamma/kperp**2_in_gyrobohm_units*<|grad_r_theta|**2>}
c          At the same beta and collisionality wrt rmin_loc n*T*tau is 
c           proportional to b_unit**3=(b0*kappa_glob)**3   
c           kappa_glob=d (rho_loc)**2/ d (rmin_loc)**2
c
c           Thus confinement at fixed b0 is proportional to
c
c      kappa_glob**3/{max_gamma/kperp**2_in_gyrobohm_units*<|grad_r_theta|**2>}
c      
c     Note: At delta_loc=0 1/<|grad_r_theta|**2>=kappa_loc**2/(1+kappa_loc**2)
c         which gives the max factor 2  volume/surface benefit over a circle
c         and it is kappa_glob which more important.
c
c     Note: We should seek to minimize gamma/(cs_loc/a)/(kperp*rhos_unit_loc)**2
c             and not just gamma
c
c      OUTPUTS:

c      b2_ave_geo                 real
c      epsl_geo(-ntml:ntm)        real
c      costheta_geo(-ntml:ntm)    real
c      sintheta_geo(-ntml:ntm)    real
c      kxoky_geo(-ntml:ntm)       real
c      pk_geo(-ntml:ntm)          real
c      b_geo(-ntml:ntm)           real
c      qrat_geo(-ntml:ntm)        real

c   diagnostic:
c      rmaj0_loc                  real
c      shat_mhd_loc               real
c      alpha_mhd_loc              real

ccccccccccccc

c      converts physical plasma machine variables to "gypinp" inputs
c      bt,zeff,zeff5
c      ne,nfast,nfastp,te,ti,q,lne,lte,lti,ls
c      rmin(a),rsurf(r),rmaj
c      read from ifwin input file
c
c      rmin="a"
c      rsurf="r"
c      unit of length: meter
c      if set rmin=1 unit of length is "a"
c         must then have input lengths in units of "a"
c      unit density: 10**13 cm**-3
c      unit temperature: keV
c      unit bt: Tesla
c


c      ntheta=32
c      nperiod=2

c      ntheta=32
c     check that 96 is comparable with 32
c      nperiod=2
c      nperiod = "nround" interface
c      ngauss=5

c**********************************************************************
       
       igeo=1
       
cgms       pi= 3.14159265
       pi=atan2(0.0,-1.0)
       
c set up theta grid   -(1+nperiod)*pi to  (1+nperiod)*pi
c j=0 is theta=0; j=ntheta/2 is theta=pi/2; j=ntheta is theta=pi

      

      ntgrid=ntheta/2+(nperiod-1)*ntheta
      ntheth=ntheta/2
      
      if(ipar.eq.0) ntgridl=ntgrid
      if(ntgridl.ne.0) ntgridl=ntgrid      

      do j=-ntgridl,ntgrid
cgms       theta(j)=real(j)*2.*3.14159265/real(ntheta)
       theta(j)=real(j)*2.*pi/real(ntheta)
      enddo
 
cgms      j_0=theta0/(2.*3.14159265/real(ntheta))
      j_0=theta0/(2.*pi/real(ntheta))
      
c quantities are even or odd in theta so we only need 
c to compute 0 to ntgrid
c to speed up calulation we could make use of 2pi perodic properties

      x_delta=asin(delta_loc) 
      
c compute rmaj_theta and z_theta in a_units

      if(rmin_loc.eq.0.) rmin_loc=rsurf

      rmaj0_loc=aspectratio_loc*rmin_loc
      kappa0=1.
      i_s_kappa_loc=0
      if(s_kappa_loc.eq.0) i_s_kappa_loc=1
      if(s_kappa_loc.eq.0)
     >       s_kappa_loc= (kappa_loc-kappa0)/kappa_loc
      i_s_delta_loc=0
      if(s_delta_loc.eq.0) i_s_delta_loc=1 
      if(s_delta_loc.eq.0.)
     >       s_delta_loc=  delta_loc/(1.-delta_loc**2)**0.5  
c      igeo_print=0
      if(igeo_print.ne.0) then
      write(6,*) " rmin_loc =", rmin_loc , " rmaj0_loc  =", rmaj0_loc
      write(6,*) " kappa_loc=", kappa_loc, " s_kappa_loc=", s_kappa_loc
      write(6,*) " delta_loc=", delta_loc, " s_delta_loc=", s_delta_loc
      write(6,*) " x_delta  =", x_delta  
      write(6,*) " shift_loc=", shift_loc 
      write(6,*) " q_loc    =", q_loc,     " shat_loc   =", shat_loc
      endif

      do j=0,ntgrid
      
       rmaj_theta(j)=
     >    rmaj0_loc+rmin_loc*cos(theta(j)+x_delta*sin(theta(j)))
     
       z_theta(j)=kappa_loc*rmin_loc*sin(theta(j))
       
       grad_r_theta(j)=
     > 1./kappa_loc*((sin(theta(j)+x_delta*sin(theta(j))))**2*
     > (1.+x_delta*cos(theta(j)))**2+
     >      kappa_loc**2*cos(theta(j))**2)**0.5/
     > (cos(x_delta*sin(theta(j)))+shift_loc*cos(theta(j))+
     > (s_kappa_loc-s_delta_loc*cos(theta(j))+(1.+s_kappa_loc)*
     > x_delta*cos(theta(j)))*sin(theta(j))*
     >   sin(theta(j)+x_delta*sin(theta(j))))
     
      enddo
      
c if grad_r_theta goes to infinity because denominator goes through
c  zero, magnetic filed lines are intersecting and 
c  the magnetic surfaces are not "nicely nested"
 
 
c compute derivative of poloidal arc lenght normed to rmin_loc
c   dlength_1/dtheta, z_l, and  r_curv=1./(rmaj_l*z_ll-z_l*rmaj_ll)
   
   
      do j=0,ntgrid
      
        dldtheta=rmin_loc*(kappa_loc**2*(cos(theta(j)))**2+
     >  (sin(theta(j)+x_delta*sin(theta(j)))*
     >        (1.+x_delta*cos(theta(j))))**2)**0.5
      
       
       ddlddtheta=rmin_loc**2*0.5/dldtheta*(
     >2.*kappa_loc**2*cos(theta(j))*(-sin(theta(j)))+
     >2.*sin(theta(j)+x_delta*sin(theta(j)))*(1.+x_delta*cos(theta(j)))
     >*(cos(theta(j)+x_delta*sin(theta(j)))*
     >       (1.+x_delta*cos(theta(j)))**2+
     > 2.*(sin(theta(j)+x_delta*sin(theta(j))))**2*
     >(1.+x_delta*cos(theta(j)))*x_delta*(-sin(theta(j)))))

     
       rmaj_l(j)=rmin_loc*(-sin(theta(j)+x_delta*sin(theta(j))))*
     >        (1.+x_delta*cos(theta(j)))/dldtheta
     
       z_l(j)=rmin_loc*kappa_loc*cos(theta(j))/dldtheta
     
       z_ll(j)=1./dldtheta*(
     >  kappa_loc*rmin_loc*(-sin(theta(j)))/dldtheta+
     >  kappa_loc*rmin_loc*cos(theta(j))*(-1./dldtheta**2)*ddlddtheta)
     
       rmaj_ll(j)=1./dldtheta*(
     >  -rmin_loc*cos(theta(j)+x_delta*sin(theta(j)))*
     >  (1.+x_delta*cos(theta(j)))**2/dldtheta
     >  -rmin_loc*sin(theta(j)+x_delta*sin(theta(j)))*
     >  x_delta*(-sin(theta(j)))/dldtheta
     >  -rmin_loc*sin(theta(j)+x_delta*sin(theta(j)))*
     >  (1.+x_delta*cos(theta(j)))*(-1./dldtheta**2)*ddlddtheta)
     
      dl1dtheta(j)=dldtheta/rmin_loc
      
      r_curv(j)=1./(rmaj_l(j)*z_ll(j)-z_l(j)*rmaj_ll(j))
   
c      write(6,*) j,rmin_loc,dldtheta,ddlddtheta
c      write(6,*) rmaj_l(j),z_l(j),rmaj_ll(j),z_ll(j),r_curv(j)
     
       
      enddo
c 
c  find the volume of shape
      delta_z=kappa_loc*rmin_loc/real(ntgrid)
      volume_loc=0.
      do j=1,ntgrid
       zint=real(j-1)*kappa_loc*rmin_loc/real(ntgrid)
       asin_arg=zint/kappa_loc/rmin_loc
        if(asin_arg.ge.1.) asin_arg=1.
        if(asin_arg.le.-1.) asin_arg=-1. 
       theta_z_out=asin(asin_arg)
       theta_z_in=pi-theta_z_out
       rmaj_out=rmaj0_loc+
     >     rmin_loc*cos(theta_z_out+x_delta*sin(theta_z_out))
       rmaj_in =rmaj0_loc+
     >     rmin_loc*cos(theta_z_in+x_delta*sin(theta_z_in))
       volume_loc=volume_loc+delta_z*0.5*2.*2.*pi*
     >  1./2.*(rmaj_out**2-rmaj_in**2)
       zint=real(j)*kappa_loc*rmin_loc/real(ntgrid)
       asin_arg=zint/kappa_loc/rmin_loc
        if(asin_arg.ge.1.) asin_arg=1.
        if(asin_arg.le.-1.) asin_arg=-1.
       theta_z_out=asin(asin_arg) 
       theta_z_in=pi-theta_z_out
       rmaj_out=rmaj0_loc+
     >     rmin_loc*cos(theta_z_out+x_delta*sin(theta_z_out))
       rmaj_in =rmaj0_loc+
     >     rmin_loc*cos(theta_z_in+x_delta*sin(theta_z_in))
       volume_loc=volume_loc+delta_z*0.5*2.*2.*pi*
     >  1./2.*(rmaj_out**2-rmaj_in**2) 
       enddo
       
 

      
c compute "f" in units of b_units vis integration 0 to 2pi
c  f -> rmaj in circle

c     bt_theta=f/rmaj_theta
c     f=b_unit*f1

      f1=0.
      dtheta_loc=theta(1)-theta(0)
      do j=0,ntheta-1
c       f1=f1+dtheta_loc/2./pi*dl1dtheta(j)/
c     &    (rmaj_theta(j)*grad_r_theta(j))
       f1=f1+dtheta_loc/2./pi*dl1dtheta(j)*(
     &    0.5/(rmaj_theta(j)*grad_r_theta(j))
     &   + 0.5/(rmaj_theta(j+1)*grad_r_theta(j+1)))
      enddo
      f1=1./f1
      
      b_unit=1.
      
c     b00_loc=f*b_unit/rmaj0_loc
      if(b00_loc.ne.0.) b_unit=b00_loc*rmaj0_loc/f1
      
      if(igeo_print.ne.0) then
      write(6,*) "b00_loc= ", b00_loc
      write(6,*) "b_unit = ", b_unit
      write(6,*) "f1 = ", f1
      endif
      f=b_unit*f1 
 
      b_norm=f/rmaj0_loc
      if(igeo_print.ne.0) then
      write(6,*) "b_norm = ", b_norm
      endif

      
c compute bt_theta, bp_theta, b_theta in unints of b_units
c bt_theta -> bt
c bp -> bt*rmin_loc/(q*rmaj)
c b_theta -> bt/(1.+eps*cos(theta))
      do j=0,ntgrid
      
       bt_theta(j)=f/rmaj_theta(j)
       
       bp_theta(j)=b_unit*rmin_loc/q_loc*grad_r_theta(j)/rmaj_theta(j)
       
       b_theta(j)=(bt_theta(j)**2+bp_theta(j)**2)**0.5

      enddo

      b_theta_2_ave=0.
      do j=1,ntheta
       b_theta_2_ave=b_theta_2_ave+0.5*dtheta_loc/(2.*pi)*
     >   (b_theta(j-1)**2+b_theta(j)**2)
      enddo
  
c compute beta_loc from beta_loc_0 if beta_loc_0 is not zero
c beta_loc_0 is p/((bt_mag_center)**2/8.pi)
c rmaj_mag_center must be in same units as rmaj0, ie a_unit
c which we take to be rmin at outermost surface

      if(beta_loc_0.ne.0.and.rmaj_mag_center.ne.0.) then
       bt_mag_center=f/rmaj_mag_center
       beta_loc=beta_loc_0/b_theta_2_ave*(bt_mag_center)**2
      endif

c compute volume_prime
      
      volume_prime=0.
      do j=1,ntheta
       volume_prime=volume_prime+0.5*dtheta_loc*(2.*pi)*rmin_loc*(
     >  dl1dtheta(j-1)/bp_theta(j-1)+dl1dtheta(j)/bp_theta(j)) 
      enddo
   
      angml=0
      do j=1,ntheta
       angml=angml+0.5*dtheta_loc*(-f)*rmin_loc*(
     > dl1dtheta(j-1)/rmaj_theta(j-1)**2/bp_theta(j-1)+
     > dl1dtheta(j)/rmaj_theta(j)**2/bp_theta(j)) 
      enddo
      q_check=-angml/(2.*pi)
  
      if(igeo_print.ne.0) then
      write(6,*)  "q_check = ",q_check
      endif


c compute Miller's "D_0", "D_p", and "D_ff'"  quantities needed for kx

      
      d_0(0)=0.
      d_p(0)=0.
      d_ffp(0)=0.
    
      sign_flip=-1.
      sign_curv=1.
      do j=1,ntgrid
  
       d_0(j)=d_0(j-1)+0.5*dtheta_loc*(f*sign_flip)*rmin_loc*(
     >   dl1dtheta(j-1)/rmaj_theta(j-1)**3/bp_theta(j-1)**2*
     >     (-2.*z_l(j-1)/rmaj_theta(j-1)+sign_curv*2./r_curv(j-1))+
     >   dl1dtheta(j)/rmaj_theta(j)**3/bp_theta(j)**2*
     >     (-2.*z_l(j)/rmaj_theta(j)+sign_curv*2./r_curv(j)))
     
       d_p(j)=d_p(j-1)+0.5*dtheta_loc*(-4.*pi*f)*rmin_loc*(
     > dl1dtheta(j-1)/rmaj_theta(j-1)**2/bp_theta(j-1)**3+
     > dl1dtheta(j)/rmaj_theta(j)**2/bp_theta(j)**3)
     
c       d_ffp(j)=d_ffp(j-1)+0.5*dtheta_loc*(-1./f)*rmin_loc*(
c     >dl1dtheta(j-1)*b_theta(j-1)**2/rmaj_theta(j-1)**2/bp_theta(j-1)**3+
c     >dl1dtheta(j)*b_theta(j)**2/rmaj_theta(j)**2/bp_theta(j)**3)
 
       d_ffp(j)=d_ffp(j-1)+0.5*dtheta_loc*(-1./f)*rmin_loc*(
     >dl1dtheta(j-1)*(1./bp_theta(j-1)**3/rmaj_theta(j-1)**2)*
     >   ((f/rmaj_theta(j-1))**2+bp_theta(j-1)**2)+
     >dl1dtheta(j)*(1./bp_theta(j)**3/rmaj_theta(j)**2)*
     >   ((f/rmaj_theta(j))**2+bp_theta(j)**2))
  
     
      enddo
  
c arclength check
      arclength=0.
      do j=1,ntheta
       arclength=arclength+0.5*dtheta_loc*rmin_loc*(
     >  dl1dtheta(j-1)+dl1dtheta(j))
      enddo
      twopiol=2.*pi/arclength
   

c begin computing gstotal ballooning mode quantities

c b_geo replaces bmaj(j)=b_theta(j)/b_unit
      do j=0,ntgrid
       b_geo(j)=b_theta(j)/b_unit
      enddo
      
c  pk_geo close to pk=2.*rmin/(rmaj*q), the coefficient of d/dtheta

       do j=0,ntgrid
        pk_geo(j)=2./rmin_loc/dl1dtheta(j)*bp_theta(j)/b_theta(j)
       enddo
     
c  the mode label k_theta=n*q_loc/rmin_loc and
c  we use kyrhos_loc  for n*q_loc/rmin_loc*rhos_unit_loc
c  which is exactly the same as for the circle 
      
c  "wstar" remains unchanged from circle with logarithmic density
c   gradients along minor axis

c   "ky*rhos" in Bessel function is kyrhos_loc*qrat_geo(j)/b_geo(j)
c   "kx*rhos"  = kxoky_geo(j)* ky*rhos
c    qrat_geo -> 1 in circle

        do j=0,ntgrid
         qrat_geo(j)=rmin_loc*b_theta(j)/rmaj_theta(j)/bp_theta(j)/q_loc
        enddo
        
c compute kxoky_geo(j)
c
c -xi*delt_perp_y=-n*gradS*b^Xrho^="ky"=n/rmaj_theta*b_theta/bp_theta
c
c in Miller's notation
c -xi*delt_perp_x=-n*gradS*rho^="kx"=
c       -n*rmaj_theta(j)*bp_theta(j)*((S_rho/rmaj_theta*bp_theta) + 
c         kx0/(n*rmaj_theta(j)*bp_theta(j))
c
c  (S_rho/rmaj_theta*bp_theta)=d_0(j)+d_p(j)*p_prime+
c         d_ffp(j)*f*f_prime
c
c  Miller uses  2*pi*q_prime=
c    -d_0(ntheta)-d_p(ntheta)*p_prime-d_ffp(ntheta)*(f*f_prime) 
c
c  (S_rho/rmaj_theta*bp_theta)=
c  1./d_ffp(ntheta)*(d_0(j)*d_ffp(ntheta)-d_0(ntheta)*d_ffp(j))+
c  1./d_ffp(ntheta)*(d_p(j)*d_ffp(ntheta)-d_p(ntheta)*d_ffp(j))*p_prime+
c  1./d_ffp(ntheta)*(-2*pi*q_prime)
 
      b2_ave_geo=b_theta_2_ave
      
c   "(ky*rhos)**2" in Ampere's law is 
c    (kyrhos_loc*qrat_geo(j))**2/b2_ave_geo(j)    

c    d psi_loc / drmin_loc =  rmin_loc/q_loc

      q_prime=1./(rmin_loc/q_loc)*shat_loc*(q_loc/rmin_loc)
     
cgms      p_prime=1./(rmin_loc/q_loc)*(-1.)*
cgms     >  ((dlntedr_loc+dlnnedr_loc)+(dlntidr_loc+dlnnidr_loc)*
cgms     >    nione_loc*tiote_loc)/(1.+nione_loc*tiote_loc) 
cgms     >  *b_theta_2_ave/(8.*pi)*beta_loc/b_unit
  
cgms      if(iptot.gt.0)then
      if(beta_loc_0.ne.0)then 
        p_prime=1./(rmin_loc/q_loc)*(-1.)*
     >  (dlnpdr_loc)*b_theta_2_ave/(8.*pi)*beta_loc/b_unit
      else
        p_prime=(q_loc/rmin_loc)*(-1.)*
     >  (dlnpdr_loc)*(beta_loc/(8.0*pi))/b_unit
      endif  
      i_alpha_mhd_loc=0
      if(alpha_mhd_loc.eq.0.) i_alpha_mhd_loc=1
      if(alpha_mhd_loc.eq.0.) then  
       shat_mhd_loc=2.*volume_loc*q_prime/volume_prime/q_loc

       alpha_mhd_loc=-2.*volume_prime/(2.*pi)**2*
     >  (volume_loc/(2.*pi**2*rmaj0_loc))**0.5*4.*pi*p_prime     
      endif
c      write(*,*)"debug",alpha_mhd_loc,p_prime
c      write(*,*)"debug",volume_loc,volume_prime
      i_beta_loc=0
      if(beta_loc.eq.0.)i_beta_loc=1
      if(beta_loc.eq.0.) then
       p_prime=1./(-2.*volume_prime/(2.*pi)**2*
     >   (volume_loc/(2.*pi**2*rmaj0_loc))**0.5*4.*pi)*alpha_mhd_loc
      beta_loc=1./(1./(rmin_loc/q_loc)*(-1.)*
     >  ((dlntedr_loc+dlnnedr_loc)+(dlntidr_loc+dlnnidr_loc)*
     >    nione_loc*tiote_loc)/(1.+nione_loc*tiote_loc) 
     >  *b_theta_2_ave/(8.*pi))*p_prime*b_unit
       if(iptot.gt.0) 
     >  beta_loc=1./(1./(rmin_loc/q_loc)*(-1.)*
     >  (dlnpdr_loc)
     >  *b_theta_2_ave/(8.*pi))*p_prime*b_unit 
      q_prime=1./(2.*volume_loc/volume_prime/q_loc)*shat_mhd_loc
      shat_loc=1./(1./(rmin_loc/q_loc)*(q_loc/rmin_loc))*q_prime*b_unit
      endif
 
      if(igeo_print.eq.2) then
       write(6,*) "alpha_mhd_loc =",alpha_mhd_loc," beta_loc =",beta_loc
       write(6,*) "shat_mhd_loc =",shat_mhd_loc, " shat_loc =",shat_loc
      endif
      if(beta_loc_0.ne.0..and.rmaj_mag_center.ne.0.)
     &bp_c=f/rmaj_mag_center*rmin_loc/q_loc/rmaj0_loc
      if(beta_loc_0.eq.0..and.rmaj_mag_center.eq.0.)
     &bp_c=f/rmaj0_loc*rmin_loc/q_loc/rmaj0_loc
cgms there is a problem here when neither case is true so bp_c is undefined
cgms added the following to fix this   
      if(rmaj_mag_center.ne.0)then
        bp_c=f/rmaj_mag_center*rmin_loc/q_loc/rmaj0_loc
      else
        bp_c=f/rmaj0_loc*rmin_loc/q_loc/rmaj0_loc
      endif
c
      alpha_c=2.*4.*pi*rmin_loc**2/rmaj0_loc/
     >  bp_c**2*
     >  ((dlntedr_loc+dlnnedr_loc)+(dlntidr_loc+dlnnidr_loc)*
     >    nione_loc*tiote_loc)/(1.+nione_loc*tiote_loc)
     >  *b_theta_2_ave/(8.*pi)*beta_loc
 
      if(iptot.gt.0) 
     > alpha_c=2.*4.*pi*rmin_loc**2/rmaj0_loc/
     >  bp_c**2*
     >  (dlnpdr_loc)
     >  *b_theta_2_ave/(8.*pi)*beta_loc 
  
      volume_c=2.*pi*rmaj0_loc*pi*rmin_loc**2
  
      ff_prime=-(2.*pi*q_prime+d_0(ntheta)+d_p(ntheta)*p_prime)
     >     /d_ffp(ntheta)
c       write(*,*)"debug d_0,d_p,d_ffp",
c     >  d_0(ntheta),d_p(ntheta),d_ffp(ntheta)
       delta_dor=x_delta
       d_prime_dor=s_delta_loc/rmin_loc
       k_prime_dor=kappa_loc/rmin_loc*s_kappa_loc 
c
       q_prime_loc = q_prime
       p_prime_loc = p_prime
c
c      write(*,*) "        f= ", f,       " ff_prime= ",ff_prime 
      write(*,*) " q_prime = ", q_prime, "  p_prime= ",p_prime       
      if(igeo_print.ne.0) then
      write(6,*)
      write(6,*) " rmin_loc= ", rmin_loc, " rmaj0_loc   = ", rmaj0_loc
      write(6,*) "arclength= ",arclength, " twopiol     = ", twopiol
      write(6,*) " volume_c= ", volume_c, " volume_loc  = ", volume_loc
      write(6,*) " shat_loc= ", shat_loc,"shat_mhd_loc= ", shat_mhd_loc
      write(6,*) "alpha_c = ", alpha_c,"alpha_mhd_loc= ", alpha_mhd_loc
      write(6,*) " beta_loc= ", beta_loc
      write(6,*) " bt_mag_center= ", bt_mag_center
      write(6,*) " rmaj_mag_center=", rmaj_mag_center
      write(6,*) " beta_loc_0=", beta_loc_0
      write(6,*) "        f= ", f,       " ff_prime= ",ff_prime 
      write(6,*) " q_prime = ", q_prime, "  p_prime= ",p_prime      
      
      write(6,*) " delta_dor =", delta_dor, " d_prime_dor=", d_prime_dor
      write(6,*) " k_prime_dor=", k_prime_dor
      endif
  
       ai3=-1./(2.*pi)*d_p(ntheta)/(2.*f)
       ai2=-1./(2.*pi)*d_0(ntheta)/(f/pi)
       ai1=(-1./(2.*pi)*d_ffp(ntheta)-q_loc/f**2)/(f/(2.*pi))

       q_prime_test =ff_prime*(q_loc/f**2+f/2./pi*ai1)+
     >                             f/pi*ai2+2.*f*ai3*p_prime


        ff_prime_plus=(q_prime-(f/pi*ai2+2.*f*ai3*p_prime))
     >   /(q_loc/f**2+f/2./pi*ai1)
        ff_prime_minus=(q_prime-(-f/pi*ai2+2.*f*ai3*p_prime))
     >   /(q_loc/f**2+f/2./pi*ai1)
  
       alimit=4.*pi*p_prime*rmaj0_loc/(ff_prime/rmaj0_loc) 

      if(igeo_print.ne.0) then
        write(6,*) " i1=",ai1," i2=",ai2," i3=",ai3
        write(6,*) "q_prime_test= ",q_prime_test 
 
        write(6,*) " ff_prime_plus= ",ff_prime_plus
        write(6,*) " ff_prime_minus= ",ff_prime_minus
        write(6,*) " alimit= ",alimit
      endif
     
      do j=0,ntgrid

        s1_check(j)=d_0(j)+d_p(j)*p_prime+d_ffp(j)*ff_prime
        kxoky_geo(j)=
     >   -(rmaj_theta(j)*bp_theta(j))**2/b_theta(j)*(
     >   1./d_ffp(ntheta)*(d_0(j)*d_ffp(ntheta)-d_0(ntheta)*d_ffp(j))+
     >   1./d_ffp(ntheta)*(d_p(j)*d_ffp(ntheta)-d_p(ntheta)*d_ffp(j))*
     >      p_prime+
     >   1./d_ffp(ntheta)*(-2*pi*d_ffp(j)*q_prime))
        kxoky_0_geo(j)=
     >   -(rmaj_theta(j)*bp_theta(j))**2/b_theta(j)*(
     >   1./d_ffp(ntheta)*(d_0(j)*d_ffp(ntheta)-d_0(ntheta)*d_ffp(j)))
        kxoky_p_geo(j)=
     >   -(rmaj_theta(j)*bp_theta(j))**2/b_theta(j)*(
     >   1./d_ffp(ntheta)*(d_p(j)*d_ffp(ntheta)-d_p(ntheta)*d_ffp(j))*
     >      p_prime)
        kxoky_s_geo(j)=
     >   -(rmaj_theta(j)*bp_theta(j))**2/b_theta(j)*(
     >   1./d_ffp(ntheta)*(-2*pi*d_ffp(j)*q_prime))
      enddo
 
      if(abs(igeot).eq.3) then
      do j=0,ntgrid
      if(igeot.gt.0)
     & kxoky_geo(j)=shat_loc*theta(j)-alpha_c*sin(theta(j))
      if(igeot.lt.0)
     & kxoky_geo(j)=shat_mhd_loc*theta(j)-alpha_mhd_loc*sin(theta(j))
      enddo
      endif
      
c kxody_geo(j) -> 0 - shift*sin(theta) + shat*(theta-theta0)
c with 1st, 2nd, and 3rd term d_0, d_p, and d_ffp terms respectively
c The d_ffp*q_prime term is only secular part in any "_geo(j)"
      
c theta(j_0)=theta0 so that the j_0 term 
c   represents kx0 which must be theta independent

c compute curvature terms
c  In the circular limit we have 
c    ((1.-lamda*bmag(j)/2)*cos(theta(j))+
c     (1.-lamda*bmag(j)/2)*kxoky*sin(theta(j)))
c        *epsl
c   where epsl=2.*rmin/rmaj and kxody=-shift*sin(theta)+shat*theta
c
c  This gets replaced by 
c     ((1.-lamda*bmag(j)/2)*costheta_geo(j)+
c      (1.-lamda*bmag(j))*costheta_p_geo(j)+
c       (1.-lamda*bmag(j)/2)*kxody_geo(j)*sintheta_geo(j))
c        *epsl_geo(j)
c
c   where epsl_geo(j)=2./rmaj0_loc*qrat_geo(j)/b_geo(j)
c
c  Note (1.-lamda*bmag(j)) is curv_drift and 
c          lamda*bmag(j)/2. is grad_b_drift
c       (1.-lamda*bmag(j)/2.) adds curv_drift and grad_b_drift 
c  costheta_p_geo is a very small O(1%) difference between
c    curve_drift and grad_b_drift


      sign2_curv=1.
      if(igeot.eq.2.or.abs(igeot).eq.3) p_prime_zero=0.
      p_prime_zero=0.
      if(i_bpar.eq.1)p_prime_zero=1.0
c      write(*,*)"debug: p_prime_zero=",p_prime_zero
      do j=0,ntgrid
      
       epsl_geo(j)=2./rmaj0_loc*qrat_geo(j)/b_geo(j)
       
       costheta_p_geo(j)=-p_prime_zero*rmaj0_loc*(
     >    4.*pi*p_prime*rmaj_theta(j)*bp_theta(j)/b_theta(j)**2)

       costheta_geo(j)=rmaj0_loc*bp_theta(j)**2/
     >    b_theta(j)**2/r_curv(j)*sign2_curv
     >   -rmaj0_loc*f**2/b_theta(j)**2/
     >    rmaj_theta(j)**3*(-z_l(j))
     >   -costheta_p_geo(j)
          
c note: (1.-lamda*bmag(j)/2.)*rmaj0_loc*(
c   -4.*pi*p_prime*rmaj_theta(j)*bp_theta(j)/b_theta(j)**2)
c          +lamda*bmag(j)/2.*rmaj0_loc*(
c    4.*pi*p_prime*rmaj_theta(j)*bp_theta(j)/b_theta(j)**2)=
c     lamda*bmaj(j)*rmaj0_loc*(
c    4.*pi*p_prime*rmaj_theta(j)*bp_theta(j)/b_theta(j)**2)     
c    

c diagnostics for dorland
      if(igeo_print.eq.1) then 
       dbpdrho1(j)=rmin_loc/r_curv(j)
       dbpdrho2(j)=rmin_loc*(-4.*pi*rmaj_theta(j)*p_prime/bp_theta(j))
       dbpdrho3(j)=rmin_loc*(-ff_prime/rmaj_theta(j)/bp_theta(j))
       write(6,*) theta(j),dbpdrho1(j)
       write(6,*) theta(j),dbpdrho2(j)
       write(6,*) theta(j),dbpdrho3(j)
      endif 

      enddo
      
      sintheta_geo(0)=0.
      sintheta_geo(ntgrid)=0.
c     j=0 corresponds to 0 where b_theta(j) is min
c      ntgrid corrsponds to pi+(nperiod-1)*2*pi where b_theta in max
      do j=1,ntgrid-1
       sintheta_geo(j)=rmaj0_loc/
     >  (rmin_loc*dl1dtheta(j))*f/(rmaj_theta(j)*b_theta(j))*
     >  (b_theta(j+1)**2-b_theta(j-1)**2)/(2.*dtheta_loc)/
     >  (2.*b_theta(j)**2)
       enddo
       
c   compute -j parts
      if(ntgridl.ne.0) then
       do j=1,ntgridl
        b_geo(-j)=b_geo(j)
        qrat_geo(-j)=qrat_geo(j)
        pk_geo(-j)=pk_geo(j)
        epsl_geo(-j)=epsl_geo(j)
        costheta_geo(-j)=costheta_geo(j)
        costheta_p_geo(-j)=costheta_p_geo(j)
        sintheta_geo(-j)=-sintheta_geo(j)
        kxoky_geo(-j)=-kxoky_geo(j)
        kxoky_0_geo(-j)=-kxoky_0_geo(j)
        kxoky_p_geo(-j)=-kxoky_p_geo(j) 
        kxoky_s_geo(-j)=-kxoky_s_geo(j)
       enddo
      endif
      
c compute MHD code check qunatities
c 
c    compute theta_bar the equal arclength angle


      arc_length=0.
      do j=1,ntheta
       arc_length=arc_length+0.5*dtheta_loc*rmin_loc*
     > (dl1dtheta(j-1)+dl1dtheta(j))
      enddo
      
      theta_bar(0)=0.
      do j=1,ntgrid
       theta_bar(j)=theta_bar(j-1)+0.5*dtheta_loc*rmin_loc*
     > (dl1dtheta(j-1)+dl1dtheta(j))*2.*pi/arc_length
      enddo
      
c compute amhd, bmhd, cmhd in 
c     amhd d/dtheta_bar (bmhd) d/dtheta_bar[phi]  + cmhd [phi]
      
      do j=0,ntgrid
       dtheta_bardtheta=rmin_loc*dl1dtheta(j)*2.*pi/arc_length
       amhd(j)=pk_geo(j)/2.*dtheta_bardtheta
       bmhd(j)=amhd(j)*qrat_geo(j)**2*(b_unit**2/b2_ave_geo)*
     > (1.+kxoky_geo(j)**2)
c introduce divergence factor
       amhd(j)=amhd(j)*b_geo(j)
       bmhd(j)=bmhd(j)/b_geo(j) 
       fac_mhd=1.
c    apparently ideal MHD acts to set grad_B=curvature so fac_mhd=2.
c    this results in p_prime term in costheta_p_geo completely
c    cancelling p_prime term in costheta_geo
       fac_mhd=2.
       cmhd(j)=beta_loc/2.*2./rmaj0_loc*
     >  ((dlntedr_loc+dlnnedr_loc)+(dlntidr_loc+dlnnidr_loc)*
     >    nione_loc*tiote_loc)/(1.+nione_loc*tiote_loc)*
     >  qrat_geo(j)*(1./b_geo(j))*
     >  (costheta_geo(j) + fac_mhd*0.5*costheta_p_geo(j)+
     >   kxoky_geo(j)*sintheta_geo(j))
       if(iptot.gt.0) then
       cmhd(j)=beta_loc/2.*2./rmaj0_loc*
     >  (dlnpdr_loc)*
     >  qrat_geo(j)*(1./b_geo(j))*
     >  (costheta_geo(j) + fac_mhd*0.5*costheta_p_geo(j)+
     >   kxoky_geo(j)*sintheta_geo(j))             
       endif
      enddo
  
        c=0.
       do j=1,ntheta
        c=c+(2./pk_geo(j-1)+2./pk_geo(j))/2.*(theta(j)-theta(j-1))
       enddo
        c=2.*pi/c
       do j=0,ntheta
        dtheta_cdtheta(j)=c*2./pk_geo(j)
       enddo
       theta_c(0)=0.
       do j=1,ntheta
        theta_c(j)=theta_c(j-1)+dtheta_cdtheta(j)*(theta(j)-theta(j-1))
       enddo 

       do j=ntheta,ntgrid
        jj=j-ntheta
        theta_c(j)=2.*pi+theta_c(jj)
        dtheta_cdtheta(j)=dtheta_cdtheta(jj)
       enddo

       do j=0,-ntgrid
        theta_c(j)=-theta_c(-j)
        dtheta_cdtheta(j)=dtheta_cdtheta(-j)
       enddo
  
c  x_geo_c are to be used for equally spaced theta_c grid
c in which pk_geo_c(j) is a constant and the kpar variable
c which is -xi*pk_geo_c*d/dtheta_c can be FFT'd
c if x_geo is replace by x_geo_c then theta effectively becomes 
c theta_c inside gstotal

       do j=0,ntgrid
        theta_c_test=theta(j)
        do jj=0,ntgrid-1
         if(theta_c_test.gt.theta_c(jj).and.
     >        theta_c_test.le.theta_c(jj+1)) then

        b_geo_c(j)=b_geo(jj)+(b_geo(jj+1)-b_geo(jj))/
     >     (theta_c(jj+1)-theta_c(jj))*(theta_c_test-theta_c(jj))
  
        qrat_geo_c(j)=qrat_geo(jj)+(qrat_geo(jj+1)-qrat_geo(jj))/
     >     (theta_c(jj+1)-theta_c(jj))*(theta_c_test-theta_c(jj)) 

        pk_geo_c(j)=pk_geo(jj)*dtheta_cdtheta(jj)+
     >    (pk_geo(jj+1)*dtheta_cdtheta(jj+1)
     >                   -pk_geo(jj)*dtheta_cdtheta(jj))/
     >     (theta_c(jj+1)-theta_c(jj))*(theta_c_test-theta_c(jj))

        epsl_geo_c(j)=epsl_geo(jj)+(epsl_geo(jj+1)-epsl_geo(jj))/
     >     (theta_c(jj+1)-theta_c(jj))*(theta_c_test-theta_c(jj))  

        costheta_geo_c(j)=costheta_geo(jj)+
     >       (costheta_geo(jj+1)-costheta_geo(jj))/
     >     (theta_c(jj+1)-theta_c(jj))*(theta_c_test-theta_c(jj))

        costheta_p_geo_c(j)=costheta_p_geo(jj)+
     >       (costheta_p_geo(jj+1)-costheta_p_geo(jj))/
     >     (theta_c(jj+1)-theta_c(jj))*(theta_c_test-theta_c(jj))

        sintheta_geo_c(j)=sintheta_geo(jj)+
     >       (sintheta_geo(jj+1)-sintheta_geo(jj))/
     >     (theta_c(jj+1)-theta_c(jj))*(theta_c_test-theta_c(jj))

        kxoky_geo_c(j)=kxoky_geo(jj)+
     >       (kxoky_geo(jj+1)-kxoky_geo(jj))/
     >     (theta_c(jj+1)-theta_c(jj))*(theta_c_test-theta_c(jj)) 

         go to 2000
         endif
        enddo
 2000 continue
       enddo
       b_geo_c(0)=b_geo(0)
       qrat_geo_c(0)=qrat_geo(0)
       pk_geo_c(0)=pk_geo(0)*dtheta_cdtheta(0)
       epsl_geo_c(0)=epsl_geo(0)
       costheta_geo_c(0)=costheta_geo(0) 
       costheta_p_geo_c(0)=costheta_p_geo(0)
       sintheta_geo_c(0)=sintheta_geo(0)
       kxoky_geo_c(0)=kxoky_geo(0)

c   compute -j parts
      if(ntgridl.ne.0) then
       do j=1,ntgridl
        b_geo_c(-j)=b_geo_c(j)
        qrat_geo_c(-j)=qrat_geo_c(j)
        pk_geo_c(-j)=pk_geo_c(j)
        epsl_geo_c(-j)=epsl_geo_c(j)
        costheta_geo_c(-j)=costheta_geo_c(j)
        costheta_p_geo_c(-j)=costheta_p_geo_c(j)
        sintheta_geo_c(-j)=-sintheta_geo_c(j)
        kxoky_geo_c(-j)=-kxoky_geo_c(j)
       enddo
      endif


c convert to ifwrite quantities

c      rmin=rmin
       rmin=1.
c should be no dependence on rmin which is the unit length
       rmaj=rmaj0_loc
       rsurf=rmin_loc
       bt=b00_loc
       ne=ne_loc
       te=te_loc
       ti=tiote_loc*te_loc
       shat=shat_loc
       q=q_loc
       lne=1./dlnnedr_loc
       lti=1./dlntidr_loc
       lte=1./dlntedr_loc
       zeff=zeff_loc            


       
crew________________geo____________________________________________________________

c      compute all other ifwrite quantities by repeating ifwrite calculations

       bakdif1=0.0
       bakdif2=0.0
       bakdif3=0.0
       bakdif4=0.0
       bakdif5=0.0

c       z2=6.
c       z3=-1.
c       z4=1.
c       z5=27.

c       amass2=6.
c  amass3=2.7e-4 is deuterium  5.4e-4 is hydrogen
c       amass3=2.7e-4
       amass4=1.0
       amass5=28.6

       epsa=rmin/rmaj
       eps=rsurf/rmaj
c       ntheta=32
c       nperiod=2
c       ngauss=5

cgms       an5=zeff5/(z5**2)
       an5 = 0.0
       an2=(zeff-1.-z5*(z5-1.)*an5)/((z2-1)*z2)
       anlt=1.-an2*z2-an5*z5
c       dbeam=nfast/ne
       dbeam=fastionfrac_loc
       an4=dbeam
cgms       an1=anlt-an4
       an1 = nione_loc
       anb=dbeam
cgms        if(an1.ne.nione_loc) then
cgms          dbeam=an1-nione_loc
cgms        endif
c       if (an1.ne.nione_loc) then
c        write(6,*) "an1=", an1, "nione_loc=", nione_loc
c       endif
       if (anlt.lt.0.) write(6,*) "anlt.lt.0"
       
       lfast=-dbeam*ne/(nfastp+1.e-10)
  
       if(lfast.eq.0.) lfast=1.e10

       
crew       shat=rmin*q/(ls*epsa)
 

crew       shift=q**2*rmaj*4.e5*ne*(te*(1./lne+1./lte)+an1*ti*
crew     &  (fprim1/rmin+1./lti)+anb*5
crew     &       /(shat*(bt*10**4)**2)

c       shift=shift
       
crew......reset to agree with brizards code
c      fprim1=((rmin/lne)*(1.-z2*an2-z5*an5)-anb*rmin/lfast)/an1

       fprim1=dlnnidr_loc

       pk=2.*epsa/q
       tprim1=dlntidr_loc
       epsl=2.*rmin/rmaj
       icv=1
       width0=3.
       
c         beta=beta_loc*ne*ti/(ne*te+ni*ti) hence

        beta=beta_loc*tiote_loc/(1.+nione_loc*tiote_loc)
c note there is a compensating factor b2_ave_geo  

crew       beta=400.*ne*te/(1.e5*bt**2)
crew CORRECTION
c       beta=400.*ne*ti/(1.e5*bt**2)   circular version
c       "beta" in gstotal is 8*pi*ne*ti/<b**2>
        if(bt.ne.0..and.ne.ne.0.) 
     >  beta=400.*ne*tiote_loc*te/(1.e5*bt**2)
        debyelorhos=0.0
        if(bt.ne.0.0.and.ne.ne.0.0) debyelorhos=cdebye*
     > 7.43e2/1.02e2*(bt*1.e4)/(ne*1.e13)**.5*(amass3/5.446D-4)**.5
c  cdebye defaults to zero on input  
       write(*,*)"debyelorhos=",debyelorhos,"bt=",bt,"ne=",ne

       fprim3=(dlnnedr_loc)/fprim1
cab    fprim2=fprim3
       fprim2=(dlnnimpdr_loc)/fprim1
       tprim2=1.
       tprim3=dlntedr_loc
       tprim4=0.0
       fprim4=rmin/lfast
       tprim5=1.
       fprim5=(dlnnedr_loc)/fprim1

       if(ne.ne.0.and.te.ne.0.) then
       vnewk4=0.0
crew   vnewk3=0.051*ne*te**(-1.5)/(ti**0.5)
crew   this close to correct for hardwire rmin=.8 deuterium with
crew   no 1/2 factor in solfp collision operator
crew  CORRECTION (used coulomb log=15.)
crew   vnewk3=0.063*ne*te**(-1.5)/(ti**0.5)*rmin*(5.4e-4/amass3)**0.5
crew forgot to multiply by 3/4*pi**.5 see HH 4.36

c      3/4*pi**.5=1.329  makes nu_ei pitch angle rate 
c      0.063*(1.329)*(1.414)=0.117
c      (2.91*10**-6)*10**13/(10**3)**1.5=0.92*10**3
c      0.92*10**3*1.3*29*15*root2*root2*100./(9.79*10**5*(10**3)**0.5)
c      0.118 (rmin/m)*(A/A_D)**0.5 where  (.00027/amass3)=A/A_D

       vnewk3=0.117*ne*te**(-1.5)/(ti**0.5)*rmin*(2.723014D-4/amass3)**0.5
cgms   need to divide by 2 to get
cgms   vnewk3 = nu_ei*a/vth_i where vth_i=sqrt(2*ti/mi)
       vnewk3=vnewk3/2.
c
cgms new calculation of vnewk3 with lnlamda factor 
cgms lnlamda and taue from NRL formulary
cgms note: for Te=1Kev, ne=10**13 lnlamda = 15.94 and taue=1.088D-3/lnlamda
       lnlamda = 15.94D0-0.5*LOG(ne)+LOG(te)
       taue = 1.088D-3*(te**1.5)/(ne*lnlamda)
cgms xnuei = 3/4 sqrt(pi)/taue
       xnuei = 1.329D0/taue
       mui = (5.446D-4)/amass3
       vthi = 9.79D5*(2.0*ti*1.D3/mui)**0.5
       vnewk3 = xnuei*rmin*100.D0/vthi
       write(*,*)"debug",temp3,mui,vnewk3*(2.0/temp3)**0.5

       vnewk3=cnewk3*vnewk3
crew this is normed to deuterium. agreement with MK valid
crew  have not traced vnewk3 CORRECTION to vnewk2 and vnewk5 ?
       vnewk2=vnewk3*0.291*(zeff-1.)*(te/ti)**1.5
       vnewk2=cnewk2*vnewk2
       vnewk1=0.0
crew       vnewk5=vnewk2*zeff5*((z5/z2)**2)/(zeff-1)
       vnewk5=vnewk2*zeff5*((z5/z2)**2)/(zeff-.999999999)

crew add vnewstar calc.
       vnewstare=zeff*2.91e-6*ne*1.e13*15./(te*1.e3)**2*rmaj*100.*q
     >     /(rsurf/rmaj)**1.5/4.19e7
       vnewstari=4.78e-8*ne*1.e13*15./(ti*1.e3)**2*rmaj*100.*q
     >     /(rsurf/rmaj)**1.5/9.79e5
      endif
      if(ne.eq.0.and.te.eq.0.) then
      
c  need relation of xnu_loc to vnewk3
c  xnu_loc=nu_ei/(sqrt(te/mi)/a_units)
c  vnewk3=[nu_ei*(ky*rho_i)]/w*_i=[xnu_loc)](2Te/Ti)**0.5
cgms correction should be 
cgms  vnewk3=[nu_ei*(ky*rho_i/2)]/w*_i=[xnu_loc)](Te/(2Ti))**0.5
cgms this agrees with GS2 documentation.
c  note: theis is pure ion case. The gstotal multiplies this but
c   (zeff + ee_collision) 
       
cgms       vnewk3=xnu_loc*(2./tiote_loc)**0.5
cgms correction should be
       vnewk3=xnu_loc/(2.*tiote_loc)**0.5
       vnewk2=vnewk3*0.291*(zeff-1.)*(1./tiote_loc)**1.5
       vnewk2=cnewk2*vnewk2
       vnewk1=0.0
       vnewk5=vnewk2*zeff5*((z5/z2)**2)/(zeff-.999999999)
       vnewk4=0.
       
      endif
c       write(*,*)"vnewk1 =",vnewk1
c       write(*,*)"vnewk2 =",vnewk2
c       write(*,*)"vnewk3 =",vnewk3
c       write(*,*)"vnewk4 =",vnewk4
c       write(*,*)"vnewk5 =",vnewk5      
c
       temp2=1.
       temp3=1./tiote_loc
       if(ti.ne.0.) temp4=50./ti
       if(ti.eq.0.) temp4=10.
c arbitray fast temp at 10.*ti
       temp5=temp2

       teti=tiote_loc

c  ifwritegeo produces the gyinp input physics variables

c        ntheta,nperiod,ngauss,icv

c        eps,shift,dbeam,shat,pk

c        epsl,width0,beta,zeff,teti,zeff5

c        fprim1,fprim2,fprim3,fprim4,fprim5

c        tprim1,tprim2,tprim3,tprim4,tprim5

c        vnewk1,vnewk2,vnewk3,vnewk4,vnewk5

c        bakdif1,bakdif2,bakdif3,bakdif4,bakdif5

c               amass2,amass3,amass4,amass5

c                temp2, temp3, temp4, temp5

c                   z2,    z3,    z4,    z5

c        PLUS the "geo" variables
  
c     compute gyrobohm norm

        
        dgyrobohmnorm=1.
      if(bt.ne.0.and.ne.ne.0.) 
     >  dgyrobohmnorm=(9.79e5*(te*1.e3)**.5)/rmin/100.
     >  *((1.02e2*(te*1.e3)**.5)/bt/1.e4)**2
     >  *(5.4e-4/amass3)**.5
  
c experiment with error correction in gstotal
c       beta=beta*2.

crrd    test of non-circular geometric factors;abs(igeot)=1 activates option
c      reverts back to ifwrite output
c     abs(igeot)=4 also suts off trapped particles

      if (abs(igeot).eq.1.or.abs(igeot).eq.4) then

      a_unit=1.
      do j=0,ntgrid

       b_geo(j)=1./(1.+rmin_loc/rmaj0_loc*cos(theta(j)))
       if(abs(igeot).eq.4) b_geo(j)=1./(1.+1.e-3*cos(theta(j))) 
       if(abs(igeot).eq.4) eps=1.e-34
       qrat_geo(j)=1.
       epsl_geo(j)=2.*a_unit/rmaj0_loc
       costheta_geo(j)=cos(theta(j))
       costheta_p_geo(j)=0.
       sintheta_geo(j)=sin(theta(j))
       if(igeot.gt.0)
     >  kxoky_geo(j)=shat_loc*theta(j)-
     >  alpha_c*sin(theta(j))
       if(igeot.lt.0)
     >  kxoky_geo(j)=shat_mhd_loc*theta(j)-
     >  alpha_mhd_loc*sin(theta(j))
       pk_geo(j)=2.*a_unit/(rmaj0_loc*q_loc)
  
       b2_ave_geo=b_unit**2
       
      enddo


c   compute -j parts
      if(ntgridl.ne.0) then
       do j=1,ntgridl
        b_geo(-j)=b_geo(j)
        qrat_geo(-j)=qrat_geo(j)
        pk_geo(-j)=pk_geo(j)
        epsl_geo(-j)=epsl_geo(j)
        costheta_geo(-j)=costheta_geo(j)
        costheta_p_geo(-j)=costheta_p_geo(j)
        sintheta_geo(-j)=-sintheta_geo(j)
        kxoky_geo(-j)=-kxoky_geo(j)
        kxoky_0_geo(-j)=-kxoky_0_geo(j)
        kxoky_p_geo(-j)=-kxoky_p_geo(j)
        kxoky_s_geo(-j)=-kxoky_s_geo(j)
       enddo
      endif  
      
      endif
   
       beta_loc_out=beta_loc
       alpha_mhd_loc_out=alpha_mhd_loc
       shat_loc_out=shat_loc
       shat_mhd_loc_out=shat_mhd_loc
       s_kappa_loc_out=s_kappa_loc
       s_delta_loc_out=s_delta_loc
       xnu_loc_out=xnu_loc
       
      
      
       if(i_s_kappa_loc.eq.1) s_kappa_loc=0.
       if(i_s_delta_loc.eq.1) s_delta_loc=0.
       if(i_alpha_mhd_loc.eq.1) alpha_mhd_loc=0.
       if(i_beta_loc.eq.1) beta_loc=0. 
   
       return
       end subroutine ifwritegeo

************************************************************************
       subroutine  mhd_crit_shoot
************************************************************************
cgms      Use(Basiscom)
cgms      Use(Common)
cgms      Use(Geo)
      
       if(abs(igeot).eq.1) then
        do j=0,ntgrid
        amhd(j)=1.

      if(igeot.gt.0)
     >rkxoky=shat_loc*theta_bar(j)-alpha_c*sin(theta_bar(j))
      if(igeot.lt.0)
     >rkxoky=shat_mhd_loc*theta_bar(j)-alpha_mhd_loc*sin(theta_bar(j))
        bmhd(j)=1.+rkxoky**2
      if(igeot.gt.0)
     >  cmhd(j)=
     >    alpha_c*(cos(theta_bar(j))+rkxoky*sin(theta_bar(j)))
      if(igeot.lt.0)
     >  cmhd(j)=
     >    alpha_mhd_loc*(cos(theta_bar(j))+rkxoky*sin(theta_bar(j)))
        enddo
       endif
     
       pmhd(0)=1.
       pmhd(1)=pmhd(0)-theta_bar(1)**2/amhd(0)*cmhd(0)*pmhd(0)
     >    /(bmhd(1)+bmhd(0))
       do j=1,ntgrid-1
        pmhd(j+1)=
     >2./(bmhd(j+1)+bmhd(j))*(theta_bar(j+1)-theta_bar(j))*(
     >     -cmhd(j)*pmhd(j)/amhd(j)*(theta_bar(j+1)-theta_bar(j-1))/2.
     > -(bmhd(j)+bmhd(j-1))/2./(theta_bar(j)-theta_bar(j-1))*pmhd(j-1)
     > +(bmhd(j+1)+bmhd(j))/2./(theta_bar(j+1)-theta_bar(j))*pmhd(j)
     > +(bmhd(j)+bmhd(j-1))/2./(theta_bar(j)-theta_bar(j-1))*pmhd(j))        
       enddo
       gamma_mhd_ave_fin=-real(pmhd(ntgrid))
               
       return
       end subroutine  mhd_crit_shoot
************************************************************************
       subroutine  mhd_crit
************************************************************************
cgms      Use(Basiscom)
cgms      Use(Common)
cgms      Use(Geo)
   

       complex spmhd(0:ntgrid)
       complex fmhd(0:ntgrid),gmhd(0:ntgrid)
       real hmhd(0:ntgrid)
       complex d_pmhd(0:ntgrid),d_fmhd(0:ntgrid)
       complex sd_pmhd(0:ntgrid),sd_fmhd(0:ntgrid)  
       complex sfmhd(0:ntgrid),sgmhd(0:ntgrid) 
       real omega_mhd(0:ntgrid),gamma_mhd(0:ntgrid)
 
       if(ntgridr.eq.0) ntgridr=ntgrid
       theta_bar_max=theta_bar(ntgridr)
  
       do j=0,ntgridr
        pmhd(j)=cos(theta_bar(j)/theta_bar(ntgridr)*pi/2.)
        gmhd(j)=0.
        fmhd(j)=0.
        hmhd(j)=amhd(j)*bmhd(j)/amhd(0)/bmhd(0)*theta_bar_max**2
       enddo
       pmhd(0)=pmhd(1)
       pmhd(ntgridr)=0.
       
       do n=1,nstep_mhd
  
        do j=0,ntgridr
         spmhd(j)=pmhd(j)
         sfmhd(j)=fmhd(j)
         sgmhd(j)=gmhd(j)
        enddo

        do j=0,ntgridr
         gmhd(j)=sgmhd(j)+delt_mhd*(pmhd(j)+spmhd(j))/2.
     >    -delt_mhd*damp_mhd*sgmhd(j)*.5/(1.+.5*delt_mhd*damp_mhd)
        enddo

        pmhd(0)=pmhd(1)
        pmhd(ntgridr)=0.
        do j=1,ntgridr-1
         d_pmhd(j)=
     >      (pmhd(j+1)-pmhd(j-1))/(theta_bar(j+1)-theta_bar(j-1))
        enddo
        d_pmhd(0)=0.
        d_pmhd(ntgridr)=
     >    (pmhd(ntgridr-1)-pmhd(ntgridr))/
     >         (theta_bar(ntgridr-1)-theta_bar(ntgridr))
  
        do j=1,ntgridr-1
         sd_pmhd(j)=
     >      (spmhd(j+1)-spmhd(j-1))/(theta_bar(j+1)-theta_bar(j-1))
        enddo
        sd_pmhd(0)=0.
        sd_pmhd(ntgridr)=
     >    (spmhd(ntgridr-1)-spmhd(ntgridr))/
     >         (theta_bar(ntgridr-1)-theta_bar(ntgridr))  

         do j=0,ntgridr
          fmhd(j)=sfmhd(j)+delt_mhd*
     >     theta_bar_max*bmhd(j)/bmhd(0)*(d_pmhd(j)+sd_pmhd(j))/2.
     >    -delt_mhd*damp_mhd*sfmhd(j)*.5/(1.+.5*delt_mhd*damp_mhd)
         enddo
     
        do j=1,ntgridr-1
         d_fmhd(j)=
     >      (fmhd(j+1)-fmhd(j-1))/(theta_bar(j+1)-theta_bar(j-1))
        enddo
        d_fmhd(0)=2.*fmhd(1)/2./(theta_bar(1))
        d_fmhd(ntgridr)=
     >    (fmhd(ntgridr-1)-fmhd(ntgridr))/
     >         (theta_bar(ntgridr-1)-theta_bar(ntgridr))
  
        do j=1,ntgridr-1
         sd_fmhd(j)=
     >      (sfmhd(j+1)-sfmhd(j-1))/(theta_bar(j+1)-theta_bar(j-1))
        enddo
        sd_fmhd(0)=2.*sfmhd(1)/2./(theta_bar(1)) 
        sd_fmhd(ntgridr)=
     >    (sfmhd(ntgridr-1)-sfmhd(ntgridr))/
     >         (theta_bar(ntgridr-1)-theta_bar(ntgridr)) 
  
        do j=1,ntgridr-1
         pmhd(j)=spmhd(j)+delt_mhd/hmhd(j)*
     >    (s_mhd*theta_bar_max*amhd(j)/amhd(0)*
     >      (d_fmhd(j)+sd_fmhd(j))/2.+
     >    theta_bar_max**2/amhd(0)/bmhd(0)*cmhd(j)*
     >      (gmhd(j)+sgmhd(j))/2.)
     >    -delt_mhd*damp_mhd*spmhd(j)*.5/(1.+.5*delt_mhd*damp_mhd)
     >    +delt_mhd*xmu_mhd*(spmhd(j+1)-2*spmhd(j)+spmhd(j-1))
     >    /(theta_bar(j+1)-theta_bar(j))
     >    /(theta_bar(j)-theta_bar(j-1))
        enddo
        pmhd(0)=pmhd(1)
        pmhd(ntgridr)=0.
  
        do j=0,ntgridr-1
         omega_mhd(j)=aimag(
     >   (pmhd(j)-spmhd(j))*2./(pmhd(j)+spmhd(j))/delt_mhd)
         gamma_mhd(j)=real(
     >   (pmhd(j)-spmhd(j))*2./(pmhd(j)+spmhd(j))/delt_mhd) 
     >   +damp_mhd
        enddo
        do j=0,ntgridr
         pmhd(j)=pmhd(j)/pmhd(0)
        enddo
  
    
c        nn=mod(n,100)

        nn=n
        omega_mhd_ave(nn)=0.
        gamma_mhd_ave(nn)=0.
        pmhd2_ave=0.
        do j=0,ntgridr-1
         omega_mhd_ave(nn)=omega_mhd_ave(nn)+
     >        omega_mhd(j)*(cabs(pmhd(j)))**2
         gamma_mhd_ave(nn)=gamma_mhd_ave(nn)+
     >        gamma_mhd(j)*(cabs(pmhd(j)))**2 
         pmhd2_ave=pmhd2_ave+(cabs(pmhd(j)))**2
        enddo
        omega_mhd_ave(nn)=omega_mhd_ave(nn)/pmhd2_ave
        gamma_mhd_ave(nn)=gamma_mhd_ave(nn)/pmhd2_ave

        omega_mhd_dev=0.
        do j=0,ntgridr-1
         omega_mhd_dev=omega_mhd_dev+
     >    (omega_mhd(j)-omega_mhd_ave(nn))**2*(cabs(pmhd(j)))**2
         gamma_mhd_dev=gamma_mhd_dev+
     >    (gamma_mhd(j)-gamma_mhd_ave(nn))**2*(cabs(pmhd(j)))**2
        enddo
        omega_mhd_dev=sqrt(omega_mhd_dev/pmhd2_ave)
        gamma_mhd_dev=sqrt(gamma_mhd_dev/pmhd2_ave)
 
        gamma_mhd_ave_fin=gamma_mhd_ave(nn)

c    end time step loop
       if(n.gt.10) then
        test=abs((gamma_mhd_ave(nn)-gamma_mhd_ave(nn-1))/
     >      gamma_mhd_ave(nn))
        if(test.lt.eps_mhd) go to 100
       endif
       enddo
 100  continue
       return
       end subroutine  mhd_crit
************************************************************************
       
************************************************************************
       subroutine profile
************************************************************************
cgms      Use(Profcom)
cgms      Use(Ifwcom)
cgms      Use(Basiscom)
cgms      Use(Common)

       
c      "profile" produces input to "ifwrite"
  
c      must read ifwin
c      must read prfin
c      must call makeprofile

c      ifwrite converts physical plasma machine variables to "gypinp"
c inputs
c      bt,zeff,zeff5
c      ne,nfast,nfastp,te,ti,q,lne,lte,lti,ls
c      rmin(a),rsurf(r),rmaj
c      read from ifwin input file
c
c      rmin="a"
c      rsurf="r"
c      unit of length: meter
c      if set rmin=1 unit of length is "a"
c         must then have input lengths in units of "a"
c      unit density: 10**13 cm**-3
c      unit temperature: keV
c      unit bt: Tesla
c
      
      do ip=ipmin,ipmax
      
       te=teprf(ip)
       ti=tiprf(ip)
       ne=neprf(ip)
       
       lte=rmin/rlteprf(ip)
       lti=rmin/rltiprf(ip)
       lne=rmin/rlneprf(ip)
       q=qprf(ip)
       ls=rmaj*qprf(ip)/shatprf(ip)
       rsurf=rmin*rprf(ip)
       
       call ifwrite

       betaprf(ip)=beta
       vnewstareprf(ip)=vnewstare
       vnewstariprf(ip)=vnewstari
       tetiprf(ip)=teti
       dgyrobohms(ip)=9.79e5*(te*1.e3)**.5/rmin/100.
     >  *(1.02e2*(te*1.e3)**.5/bt/1.e4)**2
     >  *(5.446e-4/amass3)**.5
 
      if(deltprf0.gt.0.) delt=
     >  deltprf0/(tprim1/3.+tprim3/3.+2.*fprim1/1.)

      if(deltprfgam.gt.0.) delt=
     > deltprfgam/(abs(maxagammasprf(ip))
     > /kys(iakymax(ip))/teti)

      deltprf(ip)=delt

       
!       call gstotal
       call gks
  
       
       maxagammasprfprev(ip)=maxagammasprf(ip)
       maxdgammasprfprev(ip)=maxdgammasprf(ip)

     
       do iaky=1,naky
        istepkprf(iaky,ip)=istepk(iaky)
        agammasprf(iaky,ip)=agammas(iaky)
        dgammasprf(iaky,ip)=dgammas(iaky)
        afreqsprf(iaky,ip)=afreqs(iaky)
        dfreqsprf(iaky,ip)=dfreqs(iaky)
	    dmxlsprf(iaky,ip)=dmxls(iaky)
        dmxlpsprf(iaky,ip)=dmxlps(iaky)
	
	    di_effcsprf(iaky,ip)=di_effc(iaky)
	    di_effssprf(iaky,ip)=di_effs(iaky)
	    di_effwsprf(iaky,ip)=di_effw(iaky)
	    de_effcsprf(iaky,ip)=de_effc(iaky)
	    de_effssprf(iaky,ip)=de_effs(iaky)
	    de_effwsprf(iaky,ip)=de_effw(iaky)
	    chii_effcsprf(iaky,ip)=chii_effc(iaky)
	    chii_effssprf(iaky,ip)=chii_effs(iaky)
	    chii_effwsprf(iaky,ip)=chii_effw(iaky)
	    chie_effcsprf(iaky,ip)=chie_effc(iaky)
	    chie_effssprf(iaky,ip)=chie_effs(iaky)
	    chie_effwsprf(iaky,ip)=chie_effw(iaky)
	
        gammadelprf(iaky,ip)=gammadel(iaky)
        freqdelprf(iaky,ip)=freqdel(iaky)
        tp1crtprf(iaky,ip)=tp1crt(iaky)
        tp3crtprf(iaky,ip)=tp3crt(iaky)
       enddo
       
       maxagammasprf(ip)=-1.e6
       iakymax(ip)=0
       do iaky=1,naky
        if(agammas(iaky).gt.maxagammasprf(ip))
     >    iakymax(ip)=iaky
        maxagammasprf(ip)=agammas(iakymax(ip))
       enddo
       
       maxdgammasprf(ip)=dgammas(iakymax(ip))
       maxafreqsprf(ip)=afreqs(iakymax(ip))
       maxdmxlsprf(ip)=dmxls(iakymax(ip))
       maxdmxlpsprf(ip)=dmxlps(iakymax(ip))
       
       maxdi_effcsprf(ip)=di_effc(iakymax(ip))
       maxdi_effssprf(ip)=di_effs(iakymax(ip))
       maxdi_effwsprf(ip)=di_effw(iakymax(ip))
       maxde_effcsprf(ip)=de_effc(iakymax(ip))
       maxde_effssprf(ip)=de_effs(iakymax(ip))
       maxde_effwsprf(ip)=de_effw(iakymax(ip))
       maxchii_effcsprf(ip)=chii_effc(iakymax(ip))
       maxchii_effssprf(ip)=chii_effs(iakymax(ip))
       maxchii_effwsprf(ip)=chii_effw(iakymax(ip))
       maxchie_effcsprf(ip)=chie_effc(iakymax(ip))
       maxchie_effwsprf(ip)=chie_effs(iakymax(ip))
       maxchie_effssprf(ip)=chie_effw(iakymax(ip))
       
       dmxlsgyrobohm(ip)=dmxls(iakymax(ip))*dgyrobohms(ip)
       dmxlpsgyrobohm(ip)=dmxlps(iakymax(ip))*dgyrobohms(ip)
       
       di_effcgyrobohm(ip)=di_effc(iakymax(ip))*dgyrobohms(ip)
       di_effsgyrobohm(ip)=di_effs(iakymax(ip))*dgyrobohms(ip)
       di_effwgyrobohm(ip)=di_effw(iakymax(ip))*dgyrobohms(ip)
       de_effcgyrobohm(ip)=de_effc(iakymax(ip))*dgyrobohms(ip)
       de_effsgyrobohm(ip)=de_effs(iakymax(ip))*dgyrobohms(ip)
       de_effwgyrobohm(ip)=de_effw(iakymax(ip))*dgyrobohms(ip)
       chii_effcgyrobohm(ip)=chii_effc(iakymax(ip))*dgyrobohms(ip)
       chii_effsgyrobohm(ip)=chii_effs(iakymax(ip))*dgyrobohms(ip)
       chii_effwgyrobohm(ip)=chii_effw(iakymax(ip))*dgyrobohms(ip)
       chie_effcgyrobohm(ip)=chie_effc(iakymax(ip))*dgyrobohms(ip)
       chie_effsgyrobohm(ip)=chie_effs(iakymax(ip))*dgyrobohms(ip)
       chie_effwgyrobohm(ip)=chie_effw(iakymax(ip))*dgyrobohms(ip)



      enddo


      return 
      end subroutine profile
************************************************************************
       subroutine makeprofile
************************************************************************
cgms      Use(Profcom)
cgms      Use(Ifwcom)
cgms      Use(Basiscom)
cgms      Use(Common)


c      "profile" produces input to "ifwrite"

c      must read ifwin
c      must read prfin

      aipm=ipmaxmax
      do ip=0,ipmaxmax
       zerop(ip)=0.

       aip=ip
       rprf(ip)=aip/aipm
       xxq=(rprf(ip))**2
       if(ip.eq.0) xxq=1.e-4
       if(ip.eq.ipmaxmax) xxq=1.-1.e-4

       teprf(ip)=teprf0*(1.-fteprfa)*(1.-xxq)**alfprfte
     >    +teprf0*fteprfa
       tiprf(ip)=tiprf0*(1.-ftiprfa)*(1.-xxq)**alfprfti
     >    +tiprf0*ftiprfa
       neprf(ip)=neprf0*(1.-fneprfa)*(1.-xxq)**alfprfne
     >    +neprf0*fneprfa

       rlteprf(ip)=(1.-fteprfa)*2.*xxq**.5*alfprfte
     >   *(1.-xxq)**(alfprfte-1.)/
     >   ((1.-fteprfa)*(1.-xxq)**alfprfte+fteprfa)
       rltiprf(ip)=(1.-ftiprfa)*2.*xxq**.5*alfprfti
     >   *(1.-xxq)**(alfprfti-1.)/
     >   ((1.-ftiprfa)*(1.-xxq)**alfprfti+ftiprfa)
       rlneprf(ip)=(1.-fneprfa)*2.*xxq**.5*alfprfne
     >   *(1.-xxq)**(alfprfne-1.)/
     >   ((1.-fneprfa)*(1.-xxq)**alfprfne+fneprfa)

       qprf(ip)=1./((1./qprf0-1./qprfa)/xxq/alfprfj*
     >     ((1.-xxq)-(1.-xxq)**(alfprfj+1))+1./qprfa)

       shatprf(ip)=-2.*xxq*qprf(ip)*
     > (1./qprf0-1./qprfa)/alfprfj*(-1./xxq**2*((1.-xxq)
     >  -(1-xxq)**(alfprfj+1.))
     >  +1./xxq*(-1.+(alfprfj+1.)*(1.-xxq)**alfprfj))

c note: alfprf=(qprfa/qprf0-1.) corresponds to ja=0 and shatprfa=2.


      enddo  



      return
      end subroutine makeprofile
c************************************************
c@sharegks.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine sharegks
c
c This routine communicates with subroutine ifwrite
c
!      use gks_var
!      implicit none
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c input
c 
      bt=bt_s
      zeff=zeff_s
      zeff5=0.0
      rsurf=rsurf_s
      ne=ne_s
      ni=ni_s
      nz=nz_s
      nfast=nfast_s
      nfastp=-.0001
      te=te_s
      ti=ti_s
      q=qgks_s
      lne=lne_s
      lni=lni_s
      lnz=lnz_s
      lte=lte_s
      lti=lti_s
      ls=ls_s
      amass3=amass3_s
      rmin=rmin_s
      rmaj=rmaj_s
      aky1=aky1_s
      shift=shift_s
      uprim1=uprim1_s
      uprim3=uprim3_s
      teti=ti/te
      delt=delt_s
c      icontinue=icontinue_s
c     
c output
c
       agammas_s=agammas(1)
       dgammas_s=dgammas(1)
       dtgammas_s=dtgammas(1)
       afreqs_s=afreqs(1)
       dfreqs_s=dfreqs(1)
       chie_eff_s= chie_effc(1)
       chii_eff_s= chii_effc(1)
       diff_eff_s= de_effc(1)
c
      return
      end subroutine sharegks

