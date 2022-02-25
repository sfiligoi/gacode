!         
      SUBROUTINE get_matrix
!
      USE gftm_dimensions
      USE gftm_global
      USE gftm_coeff
!
      IMPLICIT NONE
      INTEGER :: i,j,k 
!
      new_matrix=.FALSE.
!      write(*,*)"get_matrix"
!
      call get_gyro_average_xgrid
!
      call get_gyro_average_matrix
!
      call ave_theta
!
!      write(*,*)"ave_bp(1,1)=",ave_bp(1,1)
      call ave_inv0(ave_p0,ave_p0inv)
      call ave_inv0(ave_b0,ave_b0inv)
!          write(*,*)"ave_b0inv(i,j)=",ave_b0inv(1,1)
! debug
!      do i=1,nbasis
!      do j=1,nbasis
!        p0 = 0.0
!        do k=1,nbasis
!          p0 = p0 + ave_b0(i,k)*ave_b0inv(k,j)
!        enddo
!        write(*,*)"check b0/b0",i,j,p0
!      enddo
!      enddo
! debug
!
!
      END SUBROUTINE get_matrix
!
      SUBROUTINE get_gyro_average_matrix
!***************************************************************
!
!   compute the average h-bessel functions
!
!***************************************************************
      USE gftm_dimensions
      USE gftm_hermite
      USE gftm_xgrid
      USE gftm_coeff
      USE gftm_gyro_average
!
      IMPLICIT NONE
      INTEGER :: is,ie,ix,ib,jb


      REAL :: ww, zero_cut
!
      zero_cut = 1.E-12
!
      do is = ns0,ns
      do ie = 1,ne
      do ib=1,nbasis
      do jb=ib,nbasis
        mat_pe1j0s(is,ie,ib,jb) = 0.0
        mat_pe2j0s(is,ie,ib,jb) = 0.0
        mat_pe1j1s(is,ie,ib,jb) = 0.0
        mat_pe2j1s(is,ie,ib,jb) = 0.0
!
        do ix=1,nx
         ww = wx(ix)*h(ib,ix)*h(jb,ix)
!        write(*,*)i,j,k,ww,wx(k),h(i,k),h(j,k)
         mat_pe1j0s(is,ie,ib,jb) = mat_pe1j0s(is,ie,ib,jb) + ww*pe1j0sx(is,ie,ix)
         mat_pe2j0s(is,ie,ib,jb) = mat_pe2j0s(is,ie,ib,jb) + ww*pe2j0sx(is,ie,ix)
         mat_pe1j1s(is,ie,ib,jb) = mat_pe1j1s(is,ie,ib,jb) + ww*pe1j1sx(is,ie,ix)/SQRT(b2x(ix))
         mat_pe2j1s(is,ie,ib,jb) = mat_pe2j1s(is,ie,ib,jb) + ww*pe2j1sx(is,ie,ix)/SQRT(b2x(ix))
        enddo
! symmetrize
         mat_pe1j0s(is,ie,jb,ib) = mat_pe1j0s(is,ie,ib,jb)
         mat_pe2j0s(is,ie,jb,ib) = mat_pe2j0s(is,ie,ib,jb)
         mat_pe1j1s(is,ie,jb,ib) = mat_pe1j1s(is,ie,ib,jb)
         mat_pe2j1s(is,ie,jb,ib) = mat_pe2j1s(is,ie,ib,jb)
      enddo
      enddo
      enddo
      enddo
!
      END SUBROUTINE get_gyro_average_matrix
!
      SUBROUTINE ave_theta
!***********************************************************
!  compute  k-independent hermite basis averages
!***********************************************************
      USE gftm_dimensions
      USE gftm_global
      USE gftm_species
      USE gftm_hermite
      USE gftm_xgrid
      USE gftm_coeff
      USE gftm_sgrid
!
      IMPLICIT NONE
      INTEGER :: i,j,k,is
      REAL :: ww
      REAL :: gradz(nb,nb)
      REAL :: zero_cut
      REAL :: debye
!
      zero_cut = 1.E-12

! fill the Poisson equation phi multiplier px0 x-grid array
! note that pol is set in get_species
!
      do i=1,nx
        debye = debye_factor_in*b0x(i)*(ky*debye_s)**2  ! (k_per*debye length unit)^2
        p0x(i) = debye + pol
      enddo
!debug
!      write(*,*)"check p0",nx
!      do i=1,nx
!        write(*,*)i,b0x(i)
!      enddo
!      write(*,*)"debye_s = ",debye_s
!      write(*,*)"pol = ",pol
!      write(*,*)"ky = ", ky
!      write(*,*)"width_in =",width_in
!debug
!
!  compute the guass-hermite intregrals
!
       do i=1,nbasis
       do j=1,nbasis
         gradz(i,j) = 0.0
       enddo
       enddo
       do i=1,nbasis        
        if(i.lt.nbasis)then
         gradz(i,i+1)=SQRT(REAL(i)/2.0)
! gradz is an odd function
          gradz(i+1,i)=-gradz(i,i+1)
        endif
       enddo
!
       do i=1,nbasis
       do j=1,nbasis
!   initialize the averages
        ave_kpar(i,j) = 0.0
        ave_wdpar(i,j) = 0.0
        ave_wdper(i,j)=0.0
        ave_b0(i,j) = 0.0
        ave_lnB(i,j) = 0.0
        ave_p0inv(i,j) = 0.0
        ave_p0(i,j) = 0.0
        ave_kx(i,j) = 0.0
        ave_c_tor_par(i,j) = 0.0
        ave_c_tor_per(i,j) = 0.0
        ave_c_par_par(i,j) = 0.0
!
        do k=1,nx
         ww=wx(k)*h(i,k)*h(j,k)
           ave_kpar(i,j) = ave_kpar(i,j) + 0.5*wx(k)*   &
           (h(i+1,k)*SQRT(2.0*REAL(i))*h(j,k)          &
           -h(i,k)*SQRT(2.0*REAL(j))*h(j+1,k))/SQRT(b2x(k))
         ave_wdpar(i,j)    = ave_wdpar(i,j)  + ww*wdx(k)
         ave_wdper(i,j)   = ave_wdper(i,j) + ww*(wdx(k)+wdpx(k))
         ave_b0(i,j)    = ave_b0(i,j)  + ww*b0x(k)*ky*ky
         ave_lnB(i,j)   = ave_lnB(i,j) + ww*LOG(SQRT(b2x(k)))
         ave_p0inv(i,j) = ave_p0inv(i,j) + ww/p0x(k)
         ave_p0(i,j) = ave_p0(i,j) + ww*p0x(k)
         ave_kx(i,j) = ave_kx(i,j) + ww*kxx(k)
         ave_c_tor_par(i,j) = ave_c_tor_par(i,j) + ww*cx_tor_par(k)
         ave_c_tor_per(i,j) = ave_c_tor_per(i,j) + ww*cx_tor_per(k)
         ave_c_par_par(i,j) = ave_c_par_par(i,j) + ww*cx_par_par(k)
        enddo
        if(ABS(ave_kpar(i,j)).lt.zero_cut)ave_kpar(i,j) = 0.0
        if(ABS(ave_wdpar(i,j)).lt.zero_cut)ave_wdpar(i,j) = 0.0
        if(ABS(ave_wdper(i,j)).lt.zero_cut)ave_wdper(i,j) = 0.0
        if(ABS(ave_b0(i,j)).lt.zero_cut)ave_b0(i,j) = 0.0
        if(ABS(ave_lnB(i,j)).lt.zero_cut)ave_lnB(i,j) = 0.0
        if(ABS(ave_p0inv(i,j)).lt.zero_cut)ave_p0inv(i,j) = 0.0
        if(ABS(ave_p0(i,j)).lt.zero_cut)ave_p0(i,j) = 0.0
        if(ABS(ave_kx(i,j)).lt.zero_cut)ave_kx(i,j) = 0.0
        if(ABS(ave_c_tor_par(i,j)).lt.zero_cut)ave_c_tor_par(i,j) = 0.0
        if(ABS(ave_c_tor_per(i,j)).lt.zero_cut)ave_c_tor_per(i,j) = 0.0
        if(ABS(ave_c_par_par(i,j)).lt.zero_cut)ave_c_par_par(i,j) = 0.0
! symmetrize
!        ave_wdpar(j,i)    = ave_wdpar(i,j)
!        ave_wdper(j,i)    = ave_wdper(i,j)
!        ave_b0(j,i)    = ave_b0(i,j)
!        ave_lnB(j,i)   = ave_lnB(i,j)
!        ave_p0inv(j,i) = ave_p0inv(i,j)
!        ave_p0(j,i) = ave_p0(i,j)
!        ave_kx(j,i) = ave_kx(i,j)
!        ave_c_tor_par(j,i) = ave_c_tor_par(i,j)
!        ave_c_tor_per(j,i) = ave_c_tor_per(i,j)
!        ave_c_par_par(j,i) = ave_c_par_par(i,j)
        enddo
       enddo
       ave_p0_out = ave_p0(1,1)
!       write(*,*)"ave_wdpar(1,1)=",ave_wdpar(1,1)
!       write(*,*)"ave_wdper(1,1)=",ave_wdper(1,1)
!       write(*,*)"ave_p0inv(1,1) = ",ave_p0inv(1,1)
!       write(*,*)"ave_p0(1,1) = ",ave_p0(1,1)
!       write(*,*)"ave_b0(1,1) = ",ave_b0(1,1)
!       write(*,*)"ave_kx(1,1) = ",ave_kx(1,1)
!       write(*,*)"ave_tor_par(1,1) = ",ave_c_tor_par(1,1)
!       write(*,*)"ave_tor_per(1,1) = ",ave_c_tor_per(1,1)
!       write(*,*)"ave_par_par(1,1) = ",ave_c_par_par(1,1)
!
       do i=1,nbasis
       do j=1,nbasis
         gradB = 0.0
         do k=1,nbasis
           gradB = gradB + gradz(i,k)*ave_lnB(k,j) &
           - ave_lnB(i,k)*gradz(k,j)
         enddo
         ave_gradB(i,j) = gradB
       enddo
       enddo
!
!
!
!  debug
!       write(*,*)"check ave_p0inv"
!       do i=1,nbasis
!         write(*,*)"ave_p0inv(i,j) = ",(ave_p0inv(i,j),j=1,nbasis)
!       enddo
!       write(*,*)"check ave_kpar "
!       do i=1,nbasis
!       do j=1,nbasis
!        write(*,*)i,j,ave_kpar(i,j)
!       enddo
!       enddo
!
!       write(*,*)"check ave_wdpar "
!       do i=1,nbasis
!       do j=1,nbasis
!        write(*,*)i,j,ave_wdpar(i,j)
!       enddo
!       enddo
!       write(*,*)"check ave_wdper "
!       do i=1,nbasis
!       do j=1,nbasis
!        write(*,*)i,j,ave_wdper(i,j)
!       enddo
!       enddo
!       write(*,*)"check ave_gradB "
!       do i=1,nbasis
!       do j=1,nbasis
!        write(*,*)i,j,ave_gradB(i,j)
!       enddo
!       enddo
!
      END SUBROUTINE ave_theta
!
!
      SUBROUTINE ave_inv(ave_m,ave_minv)
!***************************************************************
!
!   compute the inverse matrix ave_minv
!   of the symmetric real matrix ave_m
!
!***************************************************************
!
      USE gftm_dimensions
! 
      IMPLICIT NONE
      INTEGER :: nm,is,i,j,k
      INTEGER :: lwork,info
      REAL,INTENT(IN),DIMENSION(ns,nbasis_max,nbasis_max) :: ave_m
      REAL,INTENT(OUT),DIMENSION(ns,nbasis_max,nbasis_max) :: ave_minv
      REAL :: detm,zero
!      REAL :: check
      REAL :: a(nbasis,nbasis)
      REAL :: w(nbasis)
      REAL,ALLOCATABLE,DIMENSION(:) :: work
!
      nm = nbasis
!
      if(nm.eq.1)then
        do is=ns0,ns
          ave_minv(is,1,1) = 1.0/ave_m(is,1,1)
        enddo
        go to 100
      endif
!
      if(nm.eq.2)then
        do is=ns0,ns
         detm = ave_m(is,1,1)*ave_m(is,2,2)-ave_m(is,1,2)*ave_m(is,2,1)
         if(detm.eq.0.0)detm = 1.E-12
         ave_minv(is,1,1) = ave_m(is,2,2)/detm
         ave_minv(is,2,2) = ave_m(is,1,1)/detm
         ave_minv(is,1,2) = -ave_m(is,1,2)/detm
         ave_minv(is,2,1) = -ave_m(is,2,1)/detm
        enddo
        go to 100
      endif
!
!  for nm > 2 need to use eigensytem solver
!
!  find the eigenvalues of ave_m
!
      zero = 1.0E-10
      do is = ns0,ns
       do i=1,nm
       do j=i,nm
         a(i,j) = ave_m(is,i,j)
       enddo
       enddo
       lwork=34*nm
       ALLOCATE(work(lwork))
! call LAPACK routine for symmetric real eigenvalue problem DSYEV
       call DSYEV('V','U',nm,a,nm,w,work,lwork,info)
       if(info.ne.0)CALL gftm_error(1,"DSYEV failed in ave_inv")
       DEALLOCATE(work)

!
! debug
!       write(*,*)"ave_m eigenvalues"
!       do i=1,nm
!       write(*,*)i,rr(i),ri(i)
!       enddo
!       write(*,*)"ave_m eigenvectors"
!       do i=1,nm
!       write(*,*)"******",i
!       do j=1,nm
!       write(*,*)j,vr(j,i),vi(j,i)
!       enddo
!       enddo
!
! regularize the eigenvalues
       do k=1,nm
           if(ABS(w(k)).lt.zero)then
             if(w(k).ge.0.0)then
               w(k)=zero
             else
               w(k)=-zero
             endif
           endif
       enddo
! compute ave_inv
       do i=1,nm
       do j=1,nm
         ave_minv(is,i,j) = 0.0
         do k=1,nm
           ave_minv(is,i,j) = ave_minv(is,i,j) + a(i,k)*a(j,k)/w(k)
         enddo
       enddo
       enddo
      enddo ! is loop
 100  continue
! debug check minv
!      do is=ns0,ns
!       write(*,*)"check minv",is
!        do i=1,nm
!        do j=1,nm
!         check = 0.0
!         do k=1,nm
!           check = check + ave_m(is,i,k)*ave_minv(is,k,j)
!         enddo
!         write(*,*)i,j,check
!        enddo
!        enddo
!      enddo
!
!
      END SUBROUTINE ave_inv
!
      SUBROUTINE ave_inv0(ave_m,ave_minv)
!***************************************************************
!
!   compute the inverse matrix ave_minv
!   of the symmetric real matrix ave_m
!
!***************************************************************
      USE gftm_dimensions
!
      IMPLICIT NONE
      INTEGER :: nm,i,j,k
      INTEGER :: lwork,info
      REAL :: detm,zero
!      REAL :: check
      REAL :: a(nbasis,nbasis)
      REAL :: w(nbasis)
      REAL :: work(34*nbasis)
      REAL,INTENT(IN),DIMENSION(nbasis_max,nbasis_max) :: ave_m
      REAL,INTENT(OUT),DIMENSION(nbasis_max,nbasis_max) :: ave_minv
!
      nm = nbasis
      zero = 1.0E-12
!
      if(nm.eq.1)then
        ave_minv(1,1) = 1.0/ave_m(1,1)
        go to 100
      endif
!
      if(nm.eq.2)then
        detm = ave_m(1,1)*ave_m(2,2)-ave_m(1,2)*ave_m(2,1)
        if(detm.eq.0.0)detm = zero
        ave_minv(1,1) = ave_m(2,2)/detm
        ave_minv(2,2) = ave_m(1,1)/detm
        ave_minv(1,2) = -ave_m(1,2)/detm
        ave_minv(2,1) = -ave_m(2,1)/detm
        go to 100
      endif
!
!  for nm > 2 need to use eigensytem solver
!
!
!  find the eigenvalues of ave_m
!
       do i=1,nm
       do j=i,nm
         a(i,j) = ave_m(i,j)
       enddo
       enddo
       lwork=34*nm
! call LAPACK routine for symmetric real eigenvalue problem DSYEV
       call DSYEV('V','U',nm,a,nm,w,work,lwork,info)
       if(info.ne.0)CALL gftm_error(1,"DSYEV failed in ave_inv0")
!
! debug
!       write(*,*)"ave_m eigenvalues"
!       do i=1,nm
!       write(*,*)i,rr(i),ri(i)
!       enddo
!       write(*,*)"ave_m eigenvectors"
!       do i=1,nm
!       write(*,*)"******",i
!       do j=1,nm
!       write(*,*)j,vr(j,i),vi(j,i)
!       enddo
!       enddo
!
! regularize the eigenvalues
       do k=1,nm
           if(ABS(w(k)).lt.zero)then
             if(w(k).ge.0.0)then
               w(k)=zero
             else
               w(k)=-zero
             endif
           endif
       enddo
! compute ave_inv
       do i=1,nm
       do j=1,nm
         ave_minv(i,j) = 0.0
         do k=1,nm
           ave_minv(i,j) = ave_minv(i,j) + a(i,k)*a(j,k)/w(k)
         enddo
       enddo
       enddo
 100   continue
! debug check minv
!       write(*,*)"check minv"
!        do i=1,nm
!        do j=1,nm
!         check = 0.0
!         do k=1,nm
!           check = check + ave_m(i,k)*ave_minv(k,j)
!         enddo
!         write(*,*)i,j,check
!        enddo
!        enddo
!
      END SUBROUTINE ave_inv0
! **************** begin gftm-2.0 velocity space matricies
!
      SUBROUTINE get_velocity_matrix
!
!  This subroutine initializes the parallel velocity and perpendicular energy
!  matricies and the corresponding derivative matricies needed for the
!  mirror force operator
!  G. M. Staebler Jan 6, 2022
!
      USE gftm_velocity_matrix
      IMPLICIT NONE
      INTEGER :: i
!
!  upar matrix
      mat_upar(:,:) = 0.0
      mat_upar(1,2) = SQRT(0.5)
      mat_upar(2,1) = SQRT(0.5)
      mat_upar(2,3) = 1.0
      mat_upar(3,2) = 1.0
      mat_upar(3,4) = SQRT(1.5)
      mat_upar(4,3) = SQRT(1.5)
      mat_upar(4,5) = SQRT(2.0)
      mat_upar(5,4) = SQRT(2.0)
      mat_upar(5,6) = SQRT(2.5)
      mat_upar(6,5) = SQRT(2.5)
      mat_upar(6,7) = SQRT(3.0)
      mat_upar(7,6) = SQRT(3.0)
      mat_upar(7,8) = SQRT(3.5)
      mat_upar(8,7) = SQRT(3.5)
      mat_upar(8,9) = 2.0
      mat_upar(9,8) = 2.0
      mat_upar(9,10) = SQRT(4.5)
      mat_upar(10,9) = SQRT(4.5)
      mat_upar(10,11) = SQRT(5.0)
      mat_upar(11,10) = SQRT(5.0)
      mat_upar(11,12) = SQRT(5.5)
      mat_upar(12,11) = SQRT(5.5)
      mat_upar(12,13) = SQRT(6.0)
      mat_upar(13,11) = SQRT(6.0)
! upar derivative matrix  d/dupar
      mat_dupar(:,:) = 0.0
      mat_dupar(1,2) = SQRT(0.5)
      mat_dupar(2,1) = -SQRT(0.5)
      mat_dupar(2,3) = 1.0
      mat_dupar(3,2) = -1.0
      mat_dupar(3,4) = SQRT(1.5)
      mat_dupar(4,3) = -SQRT(1.5)
      mat_dupar(4,5) = SQRT(2.0)
      mat_dupar(5,4) = -SQRT(2.0)
      mat_dupar(5,6) = SQRT(2.5)
      mat_dupar(6,5) = -SQRT(2.5)
      mat_dupar(6,7) = SQRT(3.0)
      mat_dupar(7,6) = -SQRT(3.0)
      mat_dupar(7,8) = SQRT(3.5)
      mat_dupar(8,7) = -SQRT(3.5)
      mat_dupar(8,9) = 2.0
      mat_dupar(9,8) = -2.0
      mat_dupar(9,10) = SQRT(4.5)
      mat_dupar(10,9) = -SQRT(4.5)
      mat_dupar(10,11) = SQRT(5.0)
      mat_dupar(11,10) = -SQRT(5.0)
      mat_dupar(11,12) = SQRT(5.5)
      mat_dupar(12,11) = -SQRT(5.5)
      mat_dupar(12,13) = SQRT(6.0)
      mat_dupar(13,12) = -SQRT(6.0)
!  eper matrix
      mat_eper(:,:) = 0.0
      mat_eper(1,1) = 1.0
      mat_eper(1,2) = 1.0
      mat_eper(2,1) = 1.0
      mat_eper(2,2) = 3.0
      mat_eper(2,3) = 2.0
      mat_eper(3,2) = 2.0
      mat_eper(3,3) = 5.0
      mat_eper(3,4) = 3.0
      mat_eper(4,3) = 3.0
      mat_eper(4,4) = 7.0
      mat_eper(4,5) = 4.0
      mat_eper(5,4) = 4.0
      mat_eper(5,5) = 9.0
      mat_eper(5,6) = 5.0
      mat_eper(6,5) = 5.0
      mat_eper(6,6) = 11.0
      mat_eper(6,7) = 6.0
      mat_eper(7,6) = 6.0
      mat_eper(7,7) = 13.0
!  eper derivative matrix  eper*d/deper
      mat_deper(:,:) = 0.0
      mat_deper(1,2) = 0.5
      mat_deper(2,1) = -0.5
      mat_deper(2,3) = 1.0
      mat_deper(3,2) = -1.0
      mat_deper(3,4) = 1.5
      mat_deper(4,3) = -1.5
      mat_deper(4,5) = 2.0
      mat_deper(5,4) = -2.0
      mat_deper(5,6) = 2.5
      mat_deper(6,5) = -2.5
      mat_deper(6,7) = 3.0
      mat_deper(7,6) = -3.0
!
      END SUBROUTINE get_velocity_matrix
!
!     gyro_average integrals
!
      SUBROUTINE get_gyro_average(b,n)
!
! Soubroutine to evaluate the perpendicular energy moments of the
! Bessel functions J0, J1. 
! Note: This routine is called many times to compute the polidal angle matricies
! so it needs to be reduce the output to just the n specified moments
!
      USE gftm_gyro_average
      IMPLICIT NONE
!
      INTEGER,INTENT(IN) :: n
      REAL,INTENT(IN)  :: b
      REAL :: w, b2, b4, b6, b8, b10, b12
!
      b2=b*b
      b4 = b2*b2
      b6 = b4*b2
      b8 = b6*b2
      b10 = b8*b2
      b12 = b10*b2
      w = EXP(- 0.25*b2)
      if(n.ge.1)then
        pe1j0(1)= w
        pe2j0(1) = -b2*w/4.0
        pe1j1(1) = w/2.0
        pe2j1(1) = -(b2 - 4.0)*w/8.0
      endif
      if(n.ge.2)then
        pe1j0(2) = -b2*w/4.0
        pe2j0(2) = ((b2 -4.0)**2)*w/16.0
        pe1j1(2) = -(b2 -4.0)*w/8.0
        pe2j1(2) = (b4 -16.0*b2 +48.0)*w/32.0
      endif
      if(n.ge.3)then
        pe1j0(3) = b4*w/32.0
        pe2j0(3) = -((b2 -8.0)**2)*b2*w/128.0
        pe1j1(3) = (b2 -8.0)*b2*w/64.0
        pe2j1(3) = -(b2 -8.0)*(b4 -20.0*b2 +32.0)*w/256.0
      endif
      if(n.ge.4)then
        pe1j0(4) = -b6*w/384.0
        pe2j0(4) = ((b2 -12.0)**2)*b4*w/1536.0
        pe1j1(4) = -(b2 -12.0)*b4*w/768.0
        pe2j1(4) = (b2 -4.0)*(b2 -12.0)*(b2 -24.0)*b2*w/3072.0
      endif
      if(n.ge.5)then
        pe1j0(5) = b8*w/6144.0
        pe2j0(5) = -((b2 -16.0)**2)*b6*w/24576.0
        pe1j1(5) = (b2 -16.0)*b6*w/12288.0
        pe2j1(5) = -(b2 -16.0)*(b4 -36.0*b2 +192.0)*b4*w/49152.0
      endif
      if(n.ge.6)then
        pe1j0(6) = -b10*w/122880.0
        pe2j0(6) = ((b2 -20.0)**2)*b8*w/491520.0
        pe1j1(6) = -(b2 -20.0)*b8*w/245760.0
        pe2j1(6) = (b2 - 20.0)*(b4 -44.0*b2 +320.0)*b6*w/983040.0
      endif
      if(n.ge.7)then
        pe1j0(7) = b12*w/2949120.0
        pe2j0(7) = -((b2 -24.0)**2)*b10*w/11796480.0
        pe1j1(7) = (b2 -24.0)*b10*w/5898240.0
        pe2j1(7) = -(b2 -40.0)*(b2 -24.0)*(b2 -12.0)*b8*w/23592960.0
      endif
!
      END SUBROUTINE get_gyro_average
!
      SUBROUTINE get_gyro_average_xgrid
!
!*************************************************************
!  Compute the gyro-average integrals at the hermite quadrature nodes
!*************************************************************
!
      USE gftm_dimensions
      USE gftm_global
      USE gftm_species
      USE gftm_xgrid
      USE gftm_gyro_average
!
      IMPLICIT NONE
!
      REAL :: b
      INTEGER :: ix, is, ie
!
      do ix=1,nx
      do is=ns0,ns
        b= ky*sqrt_two*vs(is)*mass(is)*SQRT(b0x(ix)/b2x(ix))/ABS(zs(is))
        CALL get_gyro_average(b,ne)
        do ie=1,ne
          pe1j0sx(is,ie,ix) = pe1j0(ie)
          pe2j0sx(is,ie,ix) = pe2j0(ie)
          pe1j1sx(is,ie,ix) = pe1j1(ie)
          pe2j1sx(is,ie,ix) = pe2j1(ie)
        enddo
      enddo
      enddo
!
      END SUBROUTINE get_gyro_average_xgrid
!
! **************** end 
