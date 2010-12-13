      SUBROUTINE get_xgrid_functions
!***********************************************************
!
!***********************************************************
      USE tglf_internal_interface
!
        new_width=.FALSE.
!      write(*,*)"get_xgrid_functions"
!
      if(igeo.eq.0)call xgrid_functions_sa
      if(igeo.ne.0)call xgrid_functions_geo
!
      END SUBROUTINE get_xgrid_functions
!
!
      SUBROUTINE xgrid_functions_sa
!***********************************************************
!
!***********************************************************
      USE tglf_dimensions
      USE tglf_internal_interface
      USE tglf_species
      USE tglf_hermite
      USE tglf_xgrid
!
      IMPLICIT NONE
      INTEGER :: i,is
      REAL :: thx,dthx,sn,cn,eps,Rx,Rx1,Rx2,kx0,f_gam
!
! debug
!      write(*,*)"shat_sa=",shat_sa,"alpha_sa=",alpha_sa
!      write(*,*)"theta0_sa=",theta0_sa,"xwell_sa=",xwell_sa
!      write(*,*)"rmin_sa=",rmin_sa,"rmaj_sa=",rmaj_sa
!      write(*,*)"q_sa=",q_sa,"b_model_sa=",b_model_sa
! debug
!
! set the units used by the tglf equations
      R_unit=rmaj_sa
      q_unit=q_sa
! fill the x-grid eikonal function arrays wdx and b0x 
      eps = rmin_sa/rmaj_sa
      kx0=0.0
      if(vexb_shear_in.ne.0.0)then
!        f_gam=get_GAM_freq()
!        kx0 = -alpha_kx0_in*(vexb_shear_in/ABS(vexb_shear_in))*(ABS(vexb_shear_in)/f_gam)**0.2
!        gamma_reference_GQ=get_GAM_freq()
        kx0 = alpha_kx0_in*TANH(vexb_shear_in/gamma_reference_GQ)
      endif
!
      do i=1,nx
        thx = width_in*x(i)
        sn = sin(thx)
        cn = cos(thx)
        Rx = 1.0 + eps*cn
        Bx(i) = 1.0/Rx
        kxx(i) = -shat_sa*(thx-theta0_sa) + alpha_sa*sn + kx0
!        wdx(i) =  -xwell_sa*MIN(1.0,alpha_sa)+ &
!         cn+sn*(shat_sa*(thx-theta0_sa) - alpha_sa*sn)
        wdx(i) =  -xwell_sa*MIN(1.0,alpha_sa)+ cn - sn*kxx(i)
!        b0x(i) = 1.0+(shat_sa*(thx-theta0_sa) - alpha_sa*sn)**2
        b0x(i) = 1.0+(kxx(i))**2
        if(b_model_sa.eq.1)then
! put B dependence into k_per**2
          b0x(i) = b0x(i)/Bx(i)**2
        endif
        if(b_model_sa.eq.2)then
! put 1/R(theta) factor into wd
          wdx(i) = wdx(i)/Rx
! put B dependence into k_per**2
          b0x(i) = b0x(i)/Bx(i)**2
        endif
! debug
!      write(*,*)i,thx,Bx(i),wdx(i),b0x(i)
!
! momentum equation stress projection coefficients
!
       cx_tor_par(i) = rmaj_sa
       cx_tor_per(i) = -rmin_sa/q_sa
       cx_par_par(i) = 1.0
      enddo
!
! compute the effective trapped fraction
!
      call get_ft_sa
!
! compute the flux surface averages
!
      B2_ave_out = 0.0
      R2_ave_out = 0.0
      RBt_ave_out = rmaj_sa
      dthx = pi_2/REAL(2*nx)
      thx = 0.0
      do i=1,2*nx
        Rx1 = 1.0 + eps*COS(thx)
        thx = thx + dthx
        Rx2 = 1.0 + eps*COS(thx)
        R2_ave_out  = R2_ave_out + dthx*(Rx1**2 + Rx2**2)/2.0
        B2_ave_out = B2_ave_out + dthx*(0.5/Rx1**2 + 0.5/Rx2**2)
      enddo
      R2_ave_out = (R2_ave_out*rmaj_sa**2)/pi_2
      B2_ave_out = B2_ave_out/pi_2
!
      END SUBROUTINE xgrid_functions_sa
!
      SUBROUTINE xgrid_functions_geo
!******************************************************************************!************************
!
!   PURPOSE: compute the geometric coefficients on the x-grid
!   
!
!******************************************************************************!************************
!
      USE tglf_dimensions
      USE tglf_internal_interface
      USE tglf_species
      USE tglf_hermite
      USE tglf_xgrid
      USE tglf_sgrid
      IMPLICIT NONE
!
      INTEGER :: i,j,k,m,is
      INTEGER :: m1,m2
      REAL :: y_x,s_x
      REAL :: Ly,dy
      REAL :: norm,ww,thx
      REAL :: sign_theta,loops
      REAL :: dkxky1,dkxky2
      REAL :: wd1,wd2
      REAL :: b1,b2
      REAL :: y1,y2
      REAL :: kx0,f_gam
      REAL :: kxx1,kxx2
      REAL :: cxtorper1,cxtorper2
      REAL :: B2x1,B2x2,R2x1,R2x2
!
!
!  find length along magnetic field y
!
      y(0)=0.D0
! pk_geo = 2 Bp/B = 2 ds/dy
      do m=1,ms
!        y(m) = y(m-1)+0.5D0*(2.D0/pk_geo(m)+2.D0/pk_geo(m-1))*ds
        y(m) = y(m-1)+s_p(m)*ds*4.0/(pk_geo(m)+pk_geo(m-1))
!      y(1)=y(0)+0.5D0*(2.D0*ds/(0.25D0*(pk_geo(1)+pk_geo(ms-1))))
!      do m=1,ms-1
!        y(m+1)=y(m-1)+2.D0*ds/(0.25D0*(pk_geo(m+1)+pk_geo(m-1)))
      enddo
! set the global units
      Ly=y(ms)
      R_unit = Rmaj_s/((qrat_geo(0)/b_geo(0))*costheta_geo(0))
      q_unit = Ly/(pi_2*R_unit)
! save f for output
       RBt_ave_out = f
!
!      write(*,*)"qrat_geo(0)=",qrat_geo(0),"b_geo(0)=",b_geo(0)
!      write(*,*)"costheta_geo(0)=",costheta_geo(0)
!      write(*,*)"sintheta_geo(0)=",sintheta_geo(0)
!      write(*,*)"Ly = ",Ly,"R_unit = ",R_unit,"Rmaj_s =",Rmaj_s
!      write(*,*)"q_unit=",q_unit
!      open(2,file='y_s.dat')
!      do m=0,ms
!        write(2,*)m,pi_2*y(m)/Ly
!      enddo
!      close(2)
!
! model for kx0 = kr/k_theta induced by vexb_shear
!
        kx0 = 0.0
        if(vexb_shear_in.ne.0.0)then
!          f_gam = get_GAM_freq()
!          kx0 = -alpha_kx0_in*(vexb_shear_in/ABS(vexb_shear_in))*(ABS(vexb_shear_in)/f_gam)**0.01
          kx0 = alpha_kx0_in*TANH(vexb_shear_in/gamma_reference_GQ)
!        write(*,*)"kx0=",kx0
        endif
!
!*************************************************************
!  begin calculation of wdx and b0x
!*************************************************************
!  compute drift wdx and pependicular wavenumber squared b0x
!  at the Hermite nodes x(i)
!  thx is the ballooning angle = 2 pi y/Ly
!  x is the argument of the Hermite basis functions = thx/width_in
!
      do i=1,nx
        thx = width_in*x(i)
	sign_theta=1.D0
	if(thx.lt.0.D0) sign_theta=-1.D0 
!
!  use quasi-periodic property of S_prime to evaluate it when thx > 2pi
!  S_prime(0)=0
!  S_prime(t+2pi*m) = m*S_prime(2pi) + S_prime(t)
!
	loops = REAL(INT(ABS(thx/pi_2)))
	y_x = Ly*(ABS(thx) - loops*pi_2)/pi_2
	if(thx.lt.0.D0)then
!
!  not always up/down symmetric so can't just use symmetry
!  for negative theta = -(t+m*2pi) use
!  S_prime(-t-m*2pi) = -(m+1)*S_prime(2pi)+S_prime(2pi-t)
!
          y_x=Ly-y_x
	  loops=loops+1.0
	endif
	do m=1,ms
	  if(y(m).ge.y_x) exit 
	enddo
!      write(*,*)"exit at",m,ms,y_x,y(ms)
        if(m.gt.ms)m=ms
        m1=m-1
        m2=m
!
! dkxky is the offset for S_prime due to the number of loops
!
        dkxky1 = sign_theta*loops*S_prime(ms)
        y1=y(m1)
        dkxky2 = sign_theta*loops*S_prime(ms)
        y2=y(m2)
!        write(*,*)"check interpolation",m1,m2
!        write(*,*)"y=",y1,y2
!        write(*,*)S_prime(m1)+dkxky1,S_prime(m2)+dkxky2
!
! note that costheta_geo and sintheta_geo are periodic so
! we can use f(-t-m*2pi) = f(2pi-t) if f(0)=f(2pi) and 0<t<2pi
!
! intepolate kxx
        kxx1 = (kx_factor(m1)*(S_prime(m1)+dkxky1)+kx0*b_geo(m1)/qrat_geo(m1)**2)*qrat_geo(m1)/b_geo(m1)
        kxx2 = (kx_factor(m2)*(S_prime(m2)+dkxky1)+kx0*b_geo(m2)/qrat_geo(m2)**2)*qrat_geo(m2)/b_geo(m2)
        kxx(i) = kxx1 +(kxx2-kxx1)*(y_x-y1)/(y2-y1)
! interpolate wdx        
	wd1 = (qrat_geo(m1)/b_geo(m1))*(costheta_geo(m1) &
        +(kx_factor(m1)*(S_prime(m1)+dkxky1)+kx0*b_geo(m1)/qrat_geo(m1)**2)*sintheta_geo(m1))
	wd2 = (qrat_geo(m2)/b_geo(m2))*(costheta_geo(m2) &
        +(kx_factor(m2)*(S_prime(m2)+dkxky2)+kx0*b_geo(m2)/qrat_geo(m2)**2)*sintheta_geo(m2))
!        write(*,*)"wd1,,wd2=",wd1,wd2
        wdx(i) = wd1 +(wd2-wd1)*(y_x-y1)/(y2-y1)
        wdx(i) = (R_unit/Rmaj_s)*wdx(i)
!        write(*,*)i,"wdx = ",wdx(i),y_x,x(i),y1,y2
! interpolate b0x
        b1 = (1.D0+(kx_factor(m1)*(S_prime(m1)+dkxky1)+kx0*b_geo(m1)/qrat_geo(m1)**2)**2) &
             *(qrat_geo(m1)/b_geo(m1))**2	
        b2 = (1.D0+(kx_factor(m2)*(S_prime(m2)+dkxky2)+kx0*b_geo(m2)/qrat_geo(m2)**2)**2) &
             *(qrat_geo(m2)/b_geo(m2))**2		
!        write(*,*)"b1,b2,b3,b4=",b1,b2,b3,b4
        b0x(i) = b1 +(b2-b1)*(y_x-y1)/(y2-y1)
        if(b0x(i).lt.0.0)then
         write(*,*)"interpolation error b0x < 0",i,b0x(i),b1,b2
         b0x(i)=(b1+b2)/2.0
        endif
!
! interpolate viscous stress projection coefficients
!
       cxtorper1 = -R(m1)*Bp(m1)/b_geo(m1)
       cxtorper2 = -R(m2)*Bp(m2)/b_geo(m2)
       cx_tor_par(i) = f/b_geo(m1) + (f/b_geo(m2)-f/b_geo(m1))*(y_x-y1)/(y2-y1)
       cx_tor_per(i) = cxtorper1 + (cxtorper2-cxtorper1)*(y_x-y1)/(y2-y1)
       cx_par_par(i) = b_geo(m1) + (b_geo(m2)-b_geo(m1))*(y_x-y1)/(y2-y1)
!
! debug
!        write(*,*)"b0x = ",b0x(i)
!        write(*,*)" y_x =",y_x," m = ",m
!        write(*,*)i," loops = ",loops
      enddo
!
!  compute flux surface averages
!
       B2_ave_out = 0.0
       R2_ave_out = 0.0
       dy = 1.0/REAL(ms)
       do i=1,ms
         B2x1 = b_geo(i-1)**2
         B2x2 = b_geo(i)**2
         B2_ave_out = B2_ave_out + dy*(B2x1+B2x2)/2.0
         R2x1 = R(i-1)**2
         R2x2 = R(i)**2
         R2_ave_out = R2_ave_out + dy*(R2x1+R2x2)/2.0
       enddo
!             
!      do m=0,ms
!        write(*,*)m,s_prime(ms-m),s_prime(ms)-s_prime(m)
!      enddo
!*************************************************************
!  end of calculation of wdx and b0x
!*************************************************************
!
! compute the effective trapped fraction
      call get_ft_geo
!
      END SUBROUTINE xgrid_functions_geo
!
!
      SUBROUTINE get_ft_sa
!
!  shifted circle version of get_ft
!
      USE tglf_dimensions
      USE tglf_internal_interface
      USE tglf_hermite
!
      IMPLICIT NONE
      INTEGER :: i,is
      REAL :: norm,ww
      REAL :: eps,theta_max,theta_eff,vshear_eff
      REAL :: sn,cn,thx,ftx,Bmax,Bmin
      REAL :: Rx,Bx,pol
!
!   compute pitch angle at bounce average boundary
!
       eps = rmin_sa/rmaj_sa
!       vshear_eff = 2.0*(ky/R_unit)*ABS(vpar_in)*sqrt_two*width_in
!       theta_eff = width_in/MAX(1.0-vshear_eff,0.1)
!       theta_max = theta_trapped_in*theta_eff*pi/sqrt_two
       theta_max = theta_trapped_in*width_in*pi/sqrt_two
       if(theta_max.gt.pi)theta_max=pi
       Bmax = 1.D0/(1.D0 + eps*COS(theta_max))
!       write(*,*)"theta_max = ",theta_max," Bmax = ",Bmax
       Bmin = 1.D0/(1.D0 + eps)
!
       if(ft_model_sa.eq.0)then
!       Gauss-Chebyshev integration of trappped fraction
!       over the bounce angle weighted by gaussian wavefunction
        ft = 0.D0
        norm = 0.D0
        do i=1,nx
         ww = COS(pi*REAL(2*i-1)/REAL(2*nx))
         thx = ww*theta_max 
         cn = COS(thx)
         Rx = (1.D0 + eps*cn)
         Bx = 1.D0/Rx
         ftx = 1.D0 - Bx/Bmax
         ftx = SQRT(ftx)
         ww = SQRT(1.D0 - ww**2)*EXP(-(thx/width_in)**2)          
         ft = ft + ww*ftx
         norm = norm + ww
        enddo
        ft = ft/norm
!        write(*,*)"ft = ",ft
       endif
!
       if(ft_model_sa.eq.1)ft = SQRT(1.D0-Bmin/Bmax)
!
       if(ft_model_sa.eq.2)then
!       Gauss-Chebyshev integration of trappped fraction
!       over the bounce angle
        ft = 0.D0
        norm = 0.D0
        do i=1,nx
         ww = COS(pi*REAL(2*i-1)/REAL(2*nx))
         thx = ww*theta_max 
         cn = COS(thx)
         Rx = (1.D0 + eps*cn)
         Bx = 1.D0/Rx
         ftx = 1.D0 - Bx/Bmax
         ftx = SQRT(ftx)
         ww = SQRT(1.D0 - ww**2)
         ft = ft + ww*ftx
         norm = norm + ww
        enddo
        ft = ft/norm
       endif
!
       if(ft_model_sa.eq.3)then
! hermite integration averaged ft with Gaussian envelope of wavefunction
         ft = 0.0
         do i=1,nx
           ww = wx(i)*h(1,i)*h(1,i)
           thx = x(i)*width_in + theta0_sa
           Rx = 1.0 + eps*COS(thx)
           Bx = 1.0/Rx
           ftx = SQRT(MAX(1.0 - Bx/Bmax,0.0))
!           write(*,*)i,"thx = ",thx,"ftx=",ftx
           ft = ft + ww*ftx
         enddo
!         write(*,*)"ft = ",ft,"ft0=",SQRT(1.D0-Bmin/Bmax)
       endif
!      write(*,*)"ft = ",ft,"ft_model_sa =",ft_model_sa
!
      END SUBROUTINE get_ft_sa
!
      SUBROUTINE get_ft_geo
!
! general geometry version of get_ft
!
      USE tglf_dimensions
      USE tglf_internal_interface
      USE tglf_sgrid
!
      IMPLICIT NONE
!
      INTEGER,PARAMETER :: nb_grid=25
      INTEGER :: i,j,k,m,is,m_max,m_min,j_max
      INTEGER :: pm(2,0:nb_grid),qm
      REAL :: Bmax,Bmin,By(0:nb_grid),delta_y(0:nb_grid)
      REAL :: Ly
      REAL :: B_bounce,kpar
      REAL :: db,test1,test2,test,bounce_y
      REAL :: dmin
!
!*************************************************************
!  begin trapped fraction model
!*************************************************************
!  find global Bmax and Bmin
!
! test case for debug
!      do m=0,ms/2
!       b_geo(m)= 1.0 +2.0*(REAL(ABS(m))/REAL(ms/2))**2 
!       b_geo(ms-m)=b_geo(m) +0.01*REAL(ABS(ms/2-m)*m)/REAL(ms/2)
!       b_geo(m)=b_geo(m) -0.01*REAL(ABS(ms/2-m)*m)/REAL(ms/2)
!      enddo
!
      Ly=y(ms)
      Bmax = b_geo(0)
      Bmin = b_geo(0)
      m_max = 0
      m_min = 0
      do m=1,ms
        if(b_geo(m).gt.Bmax)then
          Bmax = b_geo(m)
          m_max = m
	endif
        if(b_geo(m).lt.Bmin)then
          Bmin = b_geo(m)
          m_min = m
        endif
      enddo
!      write(*,*)"global Bmax = ",Bmax,"at m =",m_max
!      write(*,*)"global Bmin = ",Bmin,"at m =",m_min
!
!  make a table of Bmin=< By <= Bmax
!
      By(0)=Bmin
      dB = (Bmax - Bmin)/REAL(nb_grid)
      do i=1,nb_grid
	By(i) = By(i-1) + db
!      write(*,*)i,"By=",By(i)
      enddo
      
!
!  find pairs of m's at the same B starting at Bmin and taking the farthest pair
!
      pm(1,0)=m_min
      pm(2,0)=m_min
      do i=1,nb_grid-1
! go clockwise
        j_max = m_max - m_min
        if(j_max.lt.0)j_max=m_max+ms-m_min
	do m=1,j_max
! find the farthest gridpoint where b_geo=By
          j=m_min+m
          if(j.gt.ms)j=j-ms
	  test1=b_geo(j-1)-By(i)
	  test2=b_geo(j)-By(i)
	  if(test1*test2.le.0.D0)then
	    if(ABS(test1).lt.ABS(test2))then
	      qm=j-1
	    else
              qm=j
	    endif
!            write(*,*)i,j,"qm =",qm,"test1 =",test1,"test2=",test2
	  endif
        enddo
        pm(1,i)=qm
! go counterclockwise
        j_max = m_max - m_min
        if(j_max.lt.0)then
          j_max=ABS(j_max)
        else
          j_max=ms-j_max
        endif
	do m=1,j_max
! find the farthest gridpoint where b_geo=By
          j=m_min-m
          if(j.lt.0)j=j+ms
	  test1=b_geo(j+1)-By(i)
	  test2=b_geo(j)-By(i)
	  if(test1*test2.le.0.D0)then
	    if(ABS(test1).lt.ABS(test2))then
	      qm=j+1
	    else
              qm=j
	    endif
!            write(*,*)i,j,"qm =",qm,"test1 =",test1,"test2=",test2
	  endif
        enddo
        pm(2,i)=qm
       enddo
       pm(1,nb_grid)=m_max
       pm(2,nb_grid)=m_max
       do i=0,nb_grid
         if(pm(1,i).gt.ms.or.pm(2,i).gt.ms)then
           write(*,*)"error in get_ft_geo: pm out of bounds",pm(1,i),pm(2,i),ms
         endif
!         write(*,*)i,"pm(1,i) =",pm(1,i)," pm(2,i)=",pm(2,i)
       enddo
!		
!  now pm contains the pairs of m's with the same value of B
!
!  make a table of distances along field lines between the turning points for lookup
!
       delta_y(0)=0.D0
       do i=1,nb_grid
         if(y(pm(1,i)).gt.y(pm(2,i)))then
	  delta_y(i) = y(pm(1,i))-y(pm(2,i))
         else
           delta_y(i)=Ly+y(pm(1,i))-y(pm(2,i))
         endif
!         write(*,*)i,"delta_y =",delta_y(i)
       enddo
!
! compute trapped fraction
!
        kpar= pi_2/(Ly*sqrt_two*width_in)
	bounce_y = MIN(Ly,pi*theta_trapped_in/kpar)
!        kpar = pi_2/(Ly*sqrt_two*width_in*theta_trapped_in)  &
!               +xnu_factor_in*xnuei_in*(Ly/pi)*SQRT(mass_in(1))
!        bounce_y = MIN(Ly,pi/kpar)
!        write(*,*)"bounce_y =",bounce_y,"kpar =",kpar,"Ly=",Ly
!        write(*,*)"pi*theta_trapped/kpar =",pi*theta_trapped_in/kpar
	B_bounce = Bmax
	if(bounce_y.lt.Ly)then
	  do i=1,nb_grid
	    if(delta_y(i).gt.bounce_y)exit
          enddo
	  B_bounce = By(i-1)+(By(i)-By(i-1))* &
          (bounce_y-delta_y(i-1))/(delta_y(i)-delta_y(i-1))
!          write(*,*)i,"B_bounce =",B_bounce,Bmax
	endif
	ft = SQRT(1.D0 - Bmin/B_bounce)
!        write(*,*)"ft = ",ft
!*************************************************************
! end of trapped fraction model
!*************************************************************
 999    CONTINUE
        END SUBROUTINE get_ft_geo
!
!*************************************************************
!
!---------------------------------------------------------------
! mercier_luc.f [various callers]
!
! PURPOSE: Compute ballooning mode eikonal form of poloidally varying terms 
! in the TGLF equations using the mercier-luc local equilibrium expansion 
!
! DEFINITIONS:
!  Length is in units of rmin at the last closed 
!  flux surface (a_unit).
! NOTES:
!  To follow what's going on in this routine, it is
!  necessary to have the following papers at hand:
!
!  (1) R.L. Miller, M.S. Chu, J.M. Greene, Y.R. Lin-liu 
!      and R.E. Waltz, 'Noncircular, finite aspect ratio, 
!      local equilibrium model', Phys. Plasmas 5 (1998) 973.
!
!  (2) R.E. Waltz and R.L. Miller, 'Ion temperature gradient 
!      turbulence simulations and plasma flux surface shape',
!      Phys. Plasmas 6 (1999) 4265.
!   
!       **See MILLER_inputs for control variables** 
!
! REVISIONS:
! 13 Aug 02: jc
!  Documentation after MILLER de-spaghetti.
! 18 Nov 03: jc
!  Removed reset of s_delta and s_kappa when they 
!  are input as zero.
! 24 Nov 03: jc
!  Entire code has been rewritten for greater 
!  efficiency and clarity.  Documentation greatly 
!  improved as well.
!---------------------------------------------------------------

	SUBROUTINE mercier_luc

!-------------------------------------------
! the following must be defined from a previous call to one of the
! geometry routines miller_geo, ELITE_geo and stored in tglf_sgrid:
!      ms              ! the number of points in the s-grid (flux surface contour)
!      ds              ! the arc length differential on a flux surface
!      R(ms)           ! the major radius on the  s-grid
!      Z(ms)           ! the vertical coordinate on the s-grid
!      Bp(ms)          ! the poloidal magnetic field on the s-grid normalized to B_unit
!      q_s = local flux surface safety factor 
!      q_prime_s = dq/dpsi
!      p_prime_s = dp/dpsi 
!
        USE tglf_dimensions
        USE tglf_internal_interface
        USE tglf_sgrid
!
	IMPLICIT NONE
	INTEGER :: m,m1,m2,m3,m4
!
	REAL :: sin_u(0:ms),r_curv(0:ms)
        REAL :: psi_x(0:ms)
        REAL :: Bt(0:ms),B(0:ms)
	REAL :: d_0(0:ms),d_p(0:ms),d_ffp(0:ms)
	REAL :: delta_s, ds2
!
	REAL :: R_s, Z_s
	REAL :: R_ss, Z_ss
	REAL :: error_check
        REAL :: dq1,dq2
        REAL :: d0_s1,d0_s2
        REAL :: dp_s1,dp_s2
        REAL :: dffp_s1,dffp_s2
!
!-----------------------
!
! compute the first and second derivatives of R,Z on the s-grid
! and the local radius of curvature.
! Note that for the Mercier-Luc coordinate dR/ds = cos(u), dZ/ds = -sin(u)
! so (dR/ds)**2+(dZ/ds)**2 = 1, error_check compute the error in this relation
! to make sure that the input flux surface coordinates R(s), Z(s) are ok.
!
	delta_s = 12.0*ds 
	ds2 = 12.0*ds**2
!        write(*,*)"ms = ",ms," ds = ",ds
!        do m=0,ms
!          write(*,*)m,"R = ",R(m)," Z = ",Z(m)
!        enddo
! note that the point 0 and ms are the same so m+1->1 and m-1->ms-1 at m=0
        error_check=0.0
	do m=0,ms
          m1=m-2
          m2=m-1
          m3=m+1
          m4=m+2
          if(m.lt.2)then
           m1=ms+m-2
           m2=ms+m-1
          endif
          if(m.gt.ms-2)then
           m3=m+1-ms
           m4=m+2-ms
          endif
	  R_s = (R(m1)-8.0*R(m2)+8.0*R(m3)-R(m4))/delta_s
	  Z_s = (Z(m1)-8.0*Z(m2)+8.0*Z(m3)-Z(m4))/delta_s
          s_p(m) = SQRT(R_s**2 + Z_s**2)
	  R_ss = (-R(m1)+16.0*R(m2)-30.0*R(m)+16.0*R(m3)-R(m4))/ds2
	  Z_ss = (-Z(m1)+16.0*Z(m2)-30.0*Z(m)+16.0*Z(m3)-Z(m4))/ds2
	  r_curv(m) = (s_p(m)**3)/(R_s*Z_ss - Z_s*R_ss)
	  sin_u(m) = -Z_s/s_p(m)
!          write(*,*)m," r_curv =",r_curv(m),"sin_u =",sin_u(m)
!          write(*,*)"error=",m,ABS(s_p(m) -1.D0)
!	  error_check = error_check + (s_p(m)**2 -1.D0)**2
	enddo
!        open(3,file='curv.dat',status='replace')
!        do m=0,ms
!          write(3,*)r_curv(m),sin_u(m),s_p(m)
!        enddo
!        close(3)
!
	error_check = SQRT(error_check/DBLE(ms))
	if(error_check .gt. 1.00) &
          write(*,*)"error in s-grid derivative = ",error_check
!---------------------------------------------------------------
! Compute f=R*Bt such that the eikonal S which solves
! B*Grad(S)=0 has the correct quasi-periodicity S(s+Ls)=S(s)-2*pi*q_s
!
!      Ls
!       /  ds    
! 1/f = | ---------------
!       /   R**2 Bp 2pi q    
!       0
!
!   f -> R for a circular flux-surface
!
! Define psi_x 
       do m=0,ms
         psi_x(m)=R(m)*Bp(m)
       enddo
!
! First compute 2 pi q/f: 
!
        f = 0.0
	do m=1,ms
	  f = f &
          +s_p(m)*ds*2.0/(R(m-1)*psi_x(m-1)+R(m)*psi_x(m))
	enddo
!
	f = pi_2*q_s/f
!        write(*,*)"f = ",f,q_s
!        write(*,*)"ds=",ds
!
!-----------------------------------------------------------
!-----------------------------------------------------------
! Compute toroidal and total fields:
!
!  2    2    2
! B  = B  + B
!       p    t
!
! B  = f/R
!  t
!
	do m=0,ms
	  Bt(m) = f/R(m)
	  B(m)  = SQRT(Bt(m)**2 + Bp(m)**2)
!          write(*,*)m,"Bt =",Bt(m)," B =",B(m)
	enddo
!----------------------------------------------------------
!---------------------------------------------------------------
! Compute Miller's D  , D  and D    needed for kx.
!                   0    p      ff'
!
	d_0(0)   = 0.0
	d_p(0)   = 0.0
	d_ffp(0) = 0.0
! 
	dq1 = ds*s_p(0)*f/(R(0)*psi_x(0)**2) 
	d0_s1 = -dq1*(2.D0/r_curv(0)+2.D0*sin_u(0)/R(0))
	dp_s1 = dq1*4.0*pi*R(0)/Bp(0)
	dffp_s1 = dq1*(R(0)/Bp(0))*(B(0)/f)**2 
!
	do m=1,ms
	  dq2 = ds*s_p(m)*f/(R(m)*psi_x(m)**2) 
	  d0_s2 = -dq2*(2.D0/r_curv(m)+2.D0*sin_u(m)/R(m))
	  dp_s2 = dq2*4.D0*pi*R(m)/Bp(m)
	  dffp_s2 = dq2*(R(m)/Bp(m))*(B(m)/f)**2 
!
	  d_0(m)   = d_0(m-1)+0.5*(d0_s1+d0_s2)
	  d_p(m)   = d_p(m-1)+0.5*(dp_s1+dp_s2) 
	  d_ffp(m) = d_ffp(m-1)+0.5*(dffp_s1+dffp_s2)
!
	  d0_s1 = d0_s2
	  dp_s1 = dp_s2
	  dffp_s1 = dffp_s2
!
	enddo
!
!        d_p(0)=0.0
!        do m=1,ms
!           d_p(m)=d_p(m-1)+d_0(ms-(m-1))-d_0(ms-m)
!          write(*,*)t_s(m),d_0(m),d_p(m),d_ffp(m)
!          write(*,*)m,d_p(m),d_0(ms-m)-d_0(ms),d_0(m)
!        enddo
!
!---------------------------------------------------------------


!---------------------------------------------------------------
! Begin computing geometric quantities required for solution 
! of gyrokinetic equation:
!
! - b_geo replaces bmaj(j)=b_theta(j)/b_unit
! 
! - pk_geo is close to pk=2*rmin/(rmaj*q), the coefficient 
!   of d/dtheta
!  
! - qrat_geo -> 1 in a circle
!
! Note that for the physical quantity 
!
!         k_theta = nq/r 
! 
! we use 
! 
!         kyrhos_s = n*q_s/rmin_s*rhos_unit_s
!
! which is exactly the same as for the circle.
!
! Also, "omega_star" remains unchanged from circle with 
! logarithmic density gradients along minor axis.
!
! - "ky*rhos" in Bessel function is kyrhos_s*qrat_geo(j)/b_geo(j)
!
! - "kx*rhos" is kxoky_geo(j)*ky*rhos
!
	do m=0,ms
	  b_geo(m)    = B(m)
	  pk_geo(m)   = 2.D0*Bp(m)/B(m)
	  qrat_geo(m) = (rmin_s/R(m))*(B(m)/Bp(m))/q_s
	enddo
!---------------------------------------------------------------
!
! Determine ff_prime from:
!  
! 2 pi q_prime =  d_0(ms)
!                +d_p(ms)*p_prime
!                +d_ffp(ms)*ff_prime
!
	ff_prime = (pi_2*q_prime_s-d_0(ms)-d_p(ms)*p_prime_s) &
                   /d_ffp(ms)
!        write(*,*)"ff_prime=",ff_prime,"f=",f
!        write(*,*)"d_0,d_p,d_ffp",d_0(ms),d_p(ms),d_ffp(ms)
!---------------------------------------------------------------

!--------------------------------------------------------------
! Compute [[kx/ky]] (herein, kxoky_geo) from Waltz-Miller [2] 
! paper.  
!                       2
!                (R B_p)     S1
! kxoky_geo =   ---------- ------
!                   B       R B_p
! 
!   S1
! ------ = -(d_0(theta)+d_p(theta)*p_prime+d_ffp(theta)*ff_prime) 
!  R B_p
!
	do m=0,ms

          S_prime(m) = -(d_0(m)+d_p(m)*p_prime_s+d_ffp(m)*ff_prime)
          kx_factor(m) = (psi_x(m)**2)/B(m)
	  kxoky_geo(m) = S_prime(m)*kx_factor(m)
!        write(*,*)"check s_prime",S_prime(m)

	enddo
!---------------------------------------------------------------

!---------------------------------------------------------------
! Compute drift coefficients:
!
!
! p_prime_zero forces grad-B-curvature to zero to compensates 
! for b_par =0
!
        p_prime_zero_s = 0.0
        if(use_bpar_in)p_prime_zero_s = 1.0
!
!        write(*,*)"debug p_prime_zero",p_prime_zero_s

	do m=0,ms

	  epsl_geo(m) = 2.0/rmaj_s*qrat_geo(m)/b_geo(m)

! Waltz/Miller [[cos_p]]

	  costheta_p_geo(m) =  -p_prime_zero_s* &
          rmaj_s*(Bp(m)/B(m)**2)*(4.D0*pi*R(m)*p_prime_s)

! Waltz/Miller [[cos]

	  costheta_geo(m) =  &
      	  -rmaj_s*(Bp(m)/B(m)**2)*(Bp(m)/r_curv(m) &
          -(f**2/(Bp(m)*R(m)**3))*sin_u(m)) &
          -costheta_p_geo(m)

	enddo

!----------------------------------------------------------
! Functions which require theta-derivatives:
!
	do m=0,ms

! Waltz/Miller [[sin]]
          m1=m-2
          m2=m-1
          m3=m+1
          m4=m+2
          if(m.lt.2)then
           m1=ms+m-2
           m2=ms+m-1
          endif
          if(m.gt.ms-2)then
           m3=m+1-ms
           m4=m+2-ms
          endif
	  sintheta_geo(m) = -rmaj_s*(f/(R(m)*B(m)**2))* &
          (B(m1)-8.0*B(m2)+8.0*B(m3)-B(m4))/(delta_s*s_p(m))

	enddo
!
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!
!        call mercier_write
!        open(2,file='functions_geo.dat',status='replace')
!        do m=ms/2,0,-1
!        write(2,*)t_s(m),costheta_geo(m),sintheta_geo(m),kxoky_geo(m)
!        enddo        
!        do m=ms-1,ms/2,-1
!        write(2,*)t_s(m)+pi_2,costheta_geo(m),sintheta_geo(m) &
!          ,kx_factor(m)*(S_prime(m)-S_prime(ms))
!        enddo
!        close(2)
!
	END SUBROUTINE mercier_luc
!
!
!---------------------------------------------------------------
! mercier_write.f 
!
! PURPOSE:
!  Output of some internal MILLER variables, as well 
!  as selected values of required external variables 
!  like costheta, etc.
!
! REVISIONS:
! 30 june 2005: gms
!  Based on MILLER_write.f90 from GYRO
!---------------------------------------------------------------

      SUBROUTINE mercier_write
!
      USE tglf_internal_interface
      USE tglf_sgrid
!
      IMPLICIT NONE
      INTEGER :: j,jj
      character*7 ang(0:4)
      character*11 names(0:12)
      REAL :: theta_write(0:4)
!---------------------------------------------------

      open(unit=3,file='mercier.out',status='replace')
      write(3,10) 'rmin',rmin_loc
      write(3,10) 'rmaj0',rmaj_loc
      write(3,10) 'kappa',kappa_loc
      write(3,10) 's_kappa',s_kappa_loc
      write(3,10) 'delta',delta_loc
      write(3,10) 's_delta',s_delta_loc
      write(3,10) 'shift',shift_loc
      write(3,10) 'q',q_loc
      write(3,10) 'f',f
      write(3,10) 'ff_prime',ff_prime
      write(3,10) 'q_prime',q_prime_s
      write(3,10) 'p_prime',p_prime_s

      ang(0) = '0.00*pi'
      ang(1) = '0.50*pi'
      ang(2) = '1.00*pi'
      ang(3) = '1.50*pi'
      ang(4) = '2.00*pi'
      names(0)  = 'theta  '
      names(1)  = 't_s        '
      names(2)  = 'b_geo      '
      names(3)  = 'grad_r     '
      names(4)  = 'R          '
      names(5)  = 'Z          '
      names(6)  = 'pk_geo     '
      names(7)  = 'qrat_geo   '
      names(8)  = 'cos        '
      names(9)  = 'cos_p      '
      names(10) = 'sin        '
      names(11) = 'epsl       '
      names(12) = 'kxoky_geo  '
!
 
      theta_write(0) = 0.0
      theta_write(1) = 0.5*pi
      theta_write(2) = 1.0*pi
      theta_write(3) = 1.5*pi
      theta_write(4) = 2.0*pi
      write(3,*) 
      write(3,30)names(0),names(1),(names(j),j=2,7)
      write(3,*) 
 !
      do jj=0,4
          do j=ms,2,-1
           if(t_s(j).gt.theta_write(jj)-pi_2)exit
         enddo
         if(ABS(t_s(j-1)+pi_2-theta_write(jj)) &
           .lt.ABS(t_s(j)+pi_2-theta_write(jj)))j=j-1
!
        write(3,20)ang(jj),t_s(j),b_geo(j),q_s*R(j)*Bp(j)/rmin_loc, &
         R(j),Z(j),pk_geo(j),qrat_geo(j)
!
      enddo
!
      write(3,*) 
      write(3,30)names(0),names(1),(names(j),j=8,12)
      write(3,*)
!
      do jj=0,4
         do j=ms,2,-1
           if(t_s(j).gt.theta_write(jj)-pi_2)exit
         enddo
         if(ABS(t_s(j-1)+pi_2-theta_write(jj)) &
           .lt.ABS(t_s(j)+pi_2-theta_write(jj)))j=j-1
!
        write(3,20) ang(jj),t_s(j),costheta_geo(j),costheta_p_geo(j), &
        sintheta_geo(j),epsl_geo(j),kxoky_geo(j), &
       (S_prime(j)-S_prime(ms))*kx_factor(j)
!
      enddo
      close(3)

 10   format(t2,a,':',t20,f10.6)
 20   format(t2,a,1x,10(es10.3,1x)) 
 30   format(t2,10a)
!
      END SUBROUTINE mercier_write
!
!
!---------------------------------------------------------------
! miller_geo.f [various callers]
!
! PURPOSE:
!  Core routine for calculation of MILLER shaped flux surface
!  quantities R(theta), Z(theta), grad_r(theta)
!
! DEFINITIONS:
!  Length is in units of rmin at the last closed 
!  flux surface (a_unit).
!
!  R(r,theta) = R0(r) + r cos[theta + x_delta*sin(theta)]
!  Z(r,theta) = kappa r sin(theta)
!
!
! NOTES:
!  To follow what's going on in this routine, it is
!  necessary to have the following papers at hand:
!
!  (1) R.L. Miller, M.S. Chu, J.M. Greene, Y.R. Lin-liu 
!      and R.E. Waltz, 'Noncircular, finite aspect ratio, 
!      local equilibrium model', Phys. Plasmas 5 (1998) 973.
!
!  (2) R.E. Waltz and R.L. Miller, 'Ion temperature gradient 
!      turbulence simulations and plasma flux surface shape',
!      Phys. Plasmas 6 (1999) 4265.
!   
!       **See MILLER_inputs for control variables** 
!
! REVISIONS:
! 13 Aug 02: jc
!  Documentation after MILLER de-spaghetti.
! 18 Nov 03: jc
!  Removed reset of s_delta and s_kappa when they 
!  are input as zero.
! 24 Nov 03: jc
!  Entire code has been rewritten for greater 
!  efficiency and clarity.  Documentation greatly 
!  improved as well.
! 24 June 05: gms
!  produced this version which only computes R,Z,Bp for input 
!  into mercier_luc.f which completes the calculation of the Waltz-Miller
!  functions [[sin]],[[cos]], etc.
!---------------------------------------------------------------

      SUBROUTINE miller_geo
!
      USE tglf_internal_interface
      USE tglf_sgrid
!
      IMPLICIT NONE
!
!-------------------------------------------
!
      INTEGER,PARAMETER :: nzmax = 148, mts=5
      INTEGER :: i, j, k, m, l_theta
!-----------------------------------------------
!
      REAL :: theta, x_delta
      REAL :: dtheta
      REAL :: arg,darg
      REAL :: R_t,Z_t
      REAL :: R_r, Z_r
      REAL :: l_t, grad_r, det
      REAL :: scale_max, l_t1, arclength
      REAL :: test0, test1, test2
      REAL :: r1,r2,z1,z2
      REAL :: dx,theta_out,theta_in, rmaj_out, rmaj_in
      REAL :: theta1,theta2,dtheta1,dtheta2,error
      REAL :: arg1,arg2,save_theta2
!
!-------------------------------------------
!
!      write(*,*)"miller_geo"
      x_delta = ASIN(delta_loc)
!      write(*,*)"pi = ",pi," x_delta = ",x_delta
!
! set the flux surface constants needed for mercier-luc
!
      Rmaj_s  = rmaj_loc
      q_s     = q_loc
      if(rmin_loc.lt.0.00001)rmin_loc=0.00001
      rmin_s = rmin_loc
      p_prime_s = p_prime_loc
      q_prime_s = q_prime_loc
!
!--------------------------------------------------------------
!
! compute the arclength around the flux surface 
!
      theta = 0.D0
      arg   = theta+x_delta*sin(theta)
      darg = 1.D0+x_delta*cos(theta)
      r_t = -rmin_loc*sin(arg)*darg
      z_t = kappa_loc*rmin_loc*cos(theta) 
      l_t = SQRT(r_t**2+z_t**2)
! scale dtheta by l_t to keep mts points in each ds interval of size pi_2/ms
      dtheta = pi_2/(REAL(mts*ms)*l_t)
!
      l_t1 = l_t
      scale_max=l_t
      arclength = 0.D0
      do while(theta.lt.pi_2)
	theta = theta + dtheta
	if(theta.gt.pi_2)then
          theta=theta-dtheta
          dtheta=pi_2-theta
          theta = pi_2
        endif
!        write(*,*)"theta = ",theta,"dtheta=",dtheta
	arg   = theta+x_delta*sin(theta)
! d(arg)/dtheta
	darg = 1.D0+x_delta*cos(theta)
! dR/dtheta
	r_t = -rmin_loc*sin(arg)*darg
! dZ/dtheta
	z_t = kappa_loc*rmin_loc*cos(theta) 
! dl/dtheta
	l_t = SQRT(r_t**2+z_t**2)
! arclength along flux surface in poloidal direction
	arclength = arclength + 0.5D0*(l_t + l_t1)*dtheta
! save maximum expansion scale for later
	if(l_t.gt.scale_max) scale_max = l_t
	l_t1 = l_t			
      enddo
      Ls = arclength
!
! debug
!      write(*,*)"arclength = ", arclength
!      write(*,*)"scale_max =",scale_max
!
! Find the theta points which map to an equally spaced s-grid of ms points along the arclength
! by searching for the theta where dR**2 + dZ**2 >= ds**2 for a centered difference df=f(m+1)-f(m-1).
! This keeps the finite difference error of dR/ds, dZ/ds on the s-grid small
!
      ds = arclength/REAL(ms)
!      write(*,*)"ds=",ds
      t_s(0)=0.0
      t_s(ms)=-pi_2
!  make a first guess based on theta=0.0
      theta=0.0
      arg   = theta+x_delta*sin(theta)
      darg = 1.D0+x_delta*cos(theta)
      r_t = -rmin_loc*sin(arg)*darg
      z_t = kappa_loc*rmin_loc*cos(theta) 
      l_t = SQRT(r_t**2+z_t**2)
      dtheta = -ds/l_t
      theta=dtheta
      l_t1=l_t
!
      do m=1,ms/2
	arg   = theta+x_delta*sin(theta)
	darg = 1.D0+x_delta*cos(theta)
	r_t = -rmin_loc*sin(arg)*darg
	z_t = kappa_loc*rmin_loc*cos(theta) 
	l_t = SQRT(r_t**2+z_t**2)
        dtheta = -ds/(0.5*(l_t+l_t1))
        t_s(m)=t_s(m-1)+dtheta
        theta = t_s(m) +dtheta
        l_t1=l_t
       enddo
! distribute endpoint error over interior points
      dtheta = (t_s(ms/2)-(-pi))/REAL(ms/2)
!      write(*,*)"enpoint error =",dtheta
!      dtheta=0.0
!      t_s(ms/2)=-pi
      do m=1,ms/2
       t_s(m) = t_s(m)-REAL(m)*dtheta
       t_s(ms-m)=-pi_2 - t_s(m)
      enddo
!      write(*,*)"t_s(ms/2)+pi=",t_s(ms/2)+pi
!
!      open(2,file='t_s.dat',status='replace')
!      do m=0,ms
!        write(2,*)m,t_s(m) 
!      enddo
!      close(2)
!
!--------------------------------------------------------------
!
!
!---------------------------------------------------------------
! all equilibrium functions satisfy 
!
!                 f(0) = f(l_theta)
!
! Loop to compute most geometrical quantities needed for Mercie-Luc expansion
! R, Z, R*Bp on flux surface s-grid
!
! NOTES:
!  If grad_r_theta diverges because denominator goes 
!  through zero, magnetic field lines are intersecting 
!  and the magnetic surfaces are not nested.
!
      do m=0,ms

        theta = t_s(m)
        arg   = theta+x_delta*sin(theta)
	darg = 1.D0+x_delta*cos(theta)
	
! R(theta)
! Z(theta)

        R(m) = rmaj_loc+rmin_loc*cos(arg)
        Z(m) = kappa_loc*rmin_loc*sin(theta)
		
! dR/dtheta
! dZ/dtheta

        R_t = -rmin_loc*sin(arg)*darg
        Z_t = kappa_loc*rmin_loc*cos(theta) 

! dl/dtheta

        l_t = SQRT(R_t**2+Z_t**2)
! dR/dr
! dZ/dr
        R_r = shift_loc + cos(arg) &
              -sin(arg)*s_delta_loc*sin(theta)
        Z_r = kappa_loc*sin(theta)*(1.D0 +s_kappa_loc)
! Jacobian
        det = R_r*z_t - R_t*Z_r
! grad_r
        grad_r = ABS(l_t/det)


! Bp = (r/q) B_unit*Abs(grad_r)/R

        Bp(m) = (rmin_loc/(q_loc*R(m)))*grad_r
!
      enddo
!
!      write(*,*)"rmin=",rmin_loc,"rmaj=",rmaj_loc,"q_loc=",q_loc
!      write(*,*)shift_loc,kappa_loc,delta_loc,s_kappa_loc,s_delta_loc
!      open(3,file='Bp.dat',status='replace')
!      open(2,file='RZ.dat',status='replace')
!      do m=0,ms
!        write(3,*)m,Bp(m)
!        write(2,*)m,R(m),Z(m)
!      enddo
!      close(3)
!      close(2)
!
! Note the definitions:
!
! q_prime -> dq/dpsi = dq/dr dr/dpsi
! p_prime -> dp/dpsi = dp/dr dr/dpsi
!
!                R B_p
! and dpsi/dr = -------- = b_unit (r/q)
!               |grad r|
!
! So, we can write:
!                2
! q_prime = (q/r) s/b_unit
! p_prime = (q/r) (1/b_unit) dp/dr
!         = (q/r) (1/b_unit) p * dlnpdr
!
!
!      write(*,*)"p_prime_s =",p_prime_s,"q_prime_s =",q_prime_s
!
      END SUBROUTINE miller_geo
!
!---------------------------------------------------------------
!
      SUBROUTINE ELITE_geo
!
! interpolates the input R_ELITE,Z_ELITE,Bp_ELITE onto the s-grid
! used by mercier_luc and sets the arclength Ls and differential ds
! 
     USE tglf_internal_interface
     USE tglf_sgrid
!
     IMPLICIT NONE
!   
     INTEGER ::  i,j,k,imax,im,ip
     REAL :: arclength,drde,dzde,de
     REAL :: Rmax,Zmax,Bpmax,emax
     REAL :: B_unit,sign
     REAL :: e_length,s_length
     REAL :: da
     REAL :: a,b,c
!
! compute the arclength and find the maximum and average major radius
!
     Rmax = R_ELITE(0) 
     imax=0
     drde = (R_ELITE(1)-R_ELITE(n_ELITE))/2.0
     dzde = (Z_ELITE(1)-Z_ELITE(n_ELITE))/2.0
     da = SQRT(drde**2 + dzde**2)
     arclength = da
     do i=1,n_ELITE-1
       if(R_ELITE(i).gt.Rmax)then
         Rmax = R_ELITE(i)
         imax = i
       endif
       drde = (R_ELITE(i+1)-R_ELITE(i-1))/2.0
       dzde = (Z_ELITE(i+1)-Z_ELITE(i-1))/2.0
       da =  SQRT(drde**2 + dzde**2)
       arclength = arclength + da
     enddo 
     write(*,*)"imax = ",imax
     de = arclength/REAL(n_ELITE)
     ds = arclength/REAL(ms)
!
!  find maximum of quardratic fit through R(e):  R(e) = a + b e + c e^2
!
     im = imax-1
     if(im.lt.0)im=n_ELITE-1
     ip = imax + 1
     a = R_ELITE(im)
     b = (2.0*R_ELITE(imax)-0.5*R_ELITE(ip)-1.5*R_ELITE(im))/de
     c = (0.5*R_ELITE(ip)+0.5*R_ELITE(im)-R_ELITE(imax))/de**2
! dR/de=0 at emax=-b/(2 c)
     emax= -b/(2.0*c)
     write(*,*)"emax/de=",emax/de
     Rmax = a + b*emax + c*emax**2
! interpolate Z,Bp onto the point emax
     a = Z_ELITE(im)
     b = (2.0*Z_ELITE(imax)-0.5*Z_ELITE(ip)-1.5*Z_ELITE(im))/de
     c = (0.5*Z_ELITE(ip)+0.5*Z_ELITE(im)-Z_ELITE(imax))/de**2
     Zmax = a + b*emax + c*emax**2
     a = Bp_ELITE(im)
     b = (2.0*Bp_ELITE(imax)-0.5*Bp_ELITE(ip)-1.5*Bp_ELITE(im))/de
     c = (0.5*Bp_ELITE(ip)+0.5*Bp_ELITE(im)-Bp_ELITE(imax))/de**2
     Bpmax = a + b*emax + c*emax**2
!
! compute B_unit = d psi/dr *(q/r)
!
     rmin_s = arclength/pi_2
     B_unit = Rmax*Bpmax*q_ELITE/rmin_s
!
! interpolate R,Z,Bp onto the s-grid
!
      if(emax/de.lt.1.0)then
        imax=im
      else
        emax = emax - de
      endif
      e_length = -emax
      s_length = ds
      R(0) = Rmax
      Z(0) = 0.0
      Bp(0) = Bpmax/B_unit
      k = 0
      do i=1,n_ELITE-1
        j = imax + i
        if(j .gt. n_ELITE)then
! wrap around the flux surface 
          j = imax + i - n_ELITE
        endif
        e_length = e_length + de
        if(e_length.gt.s_length)then
          k = k+1
          R(k) = R_ELITE(j-1)+(R_ELITE(j)-R_ELITE(j-1))*(e_length-s_length)/de
          Z(k) = Z_ELITE(j-1)+(Z_ELITE(j)-Z_ELITE(j-1))*(e_length-s_length)/de - Zmax
          Bp(k) = (Bp_ELITE(j-1)+(Bp_ELITE(j)-Bp_ELITE(j-1))*(e_length-s_length)/de)/B_unit
          s_length = s_length + ds
        endif         
      enddo
      R(ms) = R(0)
      Z(ms) = Z(0)
      Bp(ms) = Bp(0)
      write(*,*)"k = ",k
!
! set remaining variables for tglf_sgrid input
!
       Ls = arclength
       Rmaj_s = Rmax
       p_prime_s = p_prime_ELITE
       q_s = q_ELITE
       q_prime_s = q_prime_ELITE
! debug
       write(*,*)"Rmax = ",Rmax
       write(*,*)"Zmax = ",Zmax
       write(*,*)"Bpmax = ",Bpmax
       write(*,*)"B_unit = ",B_unit
       write(*,*)"Ls = ",Ls
       write(*,*)"ds = ",ds
       write(*,*)"Rmaj_s =",Rmaj_s
       write(*,*)"rmin_s = ",rmin_s
       write(*,*)"p_prime_s = ",p_prime_s
       write(*,*)"q_s = ",q_s
       write(*,*)"q_prime_s = ",q_prime_s
!      open(3,file='Bp.dat',status='replace')
!      open(2,file='RZ.dat',status='replace')
!      do i=0,ms
!        write(3,*)i,Bp(i)
!        write(2,*)i,R(i),Z(i)
!      enddo
!      close(3)
!      close(2)
!
       END SUBROUTINE ELITE_geo

