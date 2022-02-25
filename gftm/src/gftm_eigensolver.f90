      SUBROUTINE gftm_eigensolver
!*********************************************
!
!*********************************************
      USE gftm_dimensions
      USE gftm_global
      USE gftm_species
      USE gftm_eigen
      USE gftm_coeff
      USE gftm_GFS
      USE gftm_sgrid
      USE gftm_velocity_matrix
      USE gftm_gyro_average
!
      IMPLICIT NONE
!
      CHARACTER(1) :: rightvectors
      INTEGER :: i,j,i1,j1,i2,j2
      INTEGER :: is,js
      INTEGER :: iu,ju,ku
      INTEGER :: ie,je
      INTEGER :: iue,jue
      INTEGER :: itot,jtot
      INTEGER :: ib,jb,kb
      REAL :: one
      REAL :: k_par0, w_d0, betapsi, betasig
      REAL :: gcut = 0.001
      REAL :: beta2,max_freq
      REAL,DIMENSION(nsm) :: vpar, vpar_shear
      REAL :: collision1, collision2
      REAL :: nuei1,nuei2
! zggev
      INTEGER :: lwork, info, ifail
      REAL,ALLOCATABLE,DIMENSION(:) :: rwork
      COMPLEX,ALLOCATABLE,DIMENSION(:,:) :: at,bt
      COMPLEX,ALLOCATABLE,DIMENSION(:,:) :: vleft,vright
      COMPLEX,ALLOCATABLE,DIMENSION(:) :: work
!
       ifail = 0
       lwork = 8*ntot
!       write(*,*)" nu = ",nu," ne = ",ne," nbasis = ",nbasis
!       write(*,*)" ntot = ",ntot," nune = ",nune," nphase = ",nphase
       ALLOCATE(rwork(lwork))
       ALLOCATE(at(ntot,ntot))
       ALLOCATE(bt(ntot,ntot))
       ALLOCATE(vleft(ntot,ntot))
       ALLOCATE(vright(ntot,ntot))
       lwork = 33*ntot
       ALLOCATE(work(lwork))
!
       k_par0 = sqrt_two/(R_unit*q_unit*width_in)
!       write(*,*)"R_unit = ",R_unit," q_unit = ",q_unit," width_in = ",width_in
       w_d0 = ky/R_unit
       betapsi = 0.0
       betasig = 0.0
       if(use_bper_in)betapsi = betae_in
       if(use_bpar_in)betasig = betae_in
!
! load velocity subspace matricies
!
     call get_velocity_matrix
     do iu = 1,nu
     do ju = 1,nu
       matu(iu,ju) = mat_upar(iu,ju)
       matdu(iu,ju) = mat_dupar(iu,ju)
     enddo
     enddo
     do iu = 1,nu
     do ju = 1,nu
       matuu(iu,ju) = 0.0
       do ku = 1,nu
         matuu(iu,ju) =  matuu(iu,ju) + matu(iu,ku)*matu(ku,ju)
       enddo
     enddo
     enddo
     do ie = 1,ne
     do je = 1,ne
       mate(ie,je) = mat_eper(ie,je)
       matde(ie,je) = mat_deper(ie,je)
     enddo
     enddo
!
!  load the mirror force matrix
!
     do ie = 1,ne
     do je = 1,ne
     do iu = 1,nu
     do ju = 1,nu
       iue = iu + nu*(ie-1)
       jue = ju + nu*(je-1)
       matmirror(iue,jue) = matu(iu,ju)*matde(ie,je)-0.5*mate(ie,je)*matdu(iu,ju)
!       write(*,*)"matmirror",iue,jue,matmirror(iue,jue)
     enddo
     enddo
     enddo
     enddo
!
! copmute the field operator dot products projectors
!
     do is = ns0,ns
     do ie = 1,ne
     do iu = 1,nu
     do ib = 1,nbasis
     do jb = 1,nbasis
       p0invpe1j0s(is,ie,jb,ib) = 0.0
       b0invpe1j0s(is,ie,jb,ib) = 0.0
       do kb = 1,nbasis
       p0invpe1j0s(is,ie,jb,ib) = p0invpe1j0s(is,ie,jb,ib) + ave_p0inv(ib,kb)*mat_pe1j0s(is,ie,kb,jb)
       b0invpe1j0s(is,ie,jb,ib) = b0invpe1j0s(is,ie,jb,ib) + ave_b0inv(ib,kb)*mat_pe1j0s(is,ie,kb,jb)
       enddo
     enddo
     enddo
     enddo
     enddo
     enddo

!
! copmute the Maxwell field projectors
!
     do is = ns0,ns
     do ie = 1,ne
     do iu = 1,nu
     do ib = 1,nbasis
     do jb = 1,nbasis
       itot = ib+nbasis*(iu-1)+nbasis*nu*(ie-1)+nbasis*nu*ne*(is-ns0)
       phib(jb,itot) = as(is)*zs(is)*one(1,iu)*p0invpe1j0s(is,ie,jb,ib)
       psib(jb,itot) = betapsi*0.5*as(is)*zs(is)*vs(is)*one(2,iu)*b0invpe1j0s(is,ie,jb,ib)
       sigb(jb,itot)= -betasig*as(is)*taus(is)*one(1,iu)*mat_pe1j1s(is,ie,jb,ib)
     enddo
     enddo
     enddo
     enddo
     enddo
!
! debug
     go to 10
     do itot=1,ntot
       do ib=1,nbasis
         write(*,*)itot,ib,"phib = ",phib(ib,itot)
       enddo
     enddo
     do is = ns0,ns
     do ie = 1,ne
     do ib = 1,nbasis
     do jb = 1,nbasis
       write(*,*)is,ie,jb,ib,"pe2j0s = ",mat_pe2j0s(is,ie,jb,ib)
     enddo
     enddo
     enddo
     enddo
     do ib = 1,nbasis
     do jb = 1,nbasis
       write(*,*)ib,jb,"ave_wdpar = ",ave_wdpar(ib,jb)
       write(*,*)ib,jb,"ave_wdper = ",ave_wdper(ib,jb)
     enddo
     enddo
     do ib = 1,nbasis
     do jb = 1,nbasis
       write(*,*)ib,jb,"ave_kpar = ",ave_kpar(ib,jb)
       write(*,*)ib,jb," ave_gradb = ",ave_gradb(ib,jb)
     enddo
     enddo
10   continue
!
! copmute the gyro-average Maxwell field projectors
!
! first do the dot product on the subspace without nu for speed
!
     do is = 1,ns
     do ie = 1,ne
     do ib = 1,nbasis
     do jtot = 1,ntot
       pe1j0phib(is,ie,ib,jtot) = 0.0
       pe1j0psib(is,ie,ib,jtot) = 0.0
       pe1j1sigb(is,ie,ib,jtot) = 0.0
       pe2j0phib(is,ie,ib,jtot) = 0.0
       pe2j0psib(is,ie,ib,jtot) = 0.0
       pe2j1sigb(is,ie,ib,jtot) = 0.0
       do jb = 1,nbasis
         pe1j0phib(is,ie,ib,jtot) = pe1j0phib(is,ie,ib,jtot) + mat_pe1j0s(is,ie,ib,jb)*phib(jb,jtot)
         pe1j0psib(is,ie,ib,jtot) = pe1j0psib(is,ie,ib,jtot) + mat_pe1j0s(is,ie,ib,jb)*psib(jb,jtot)
         pe1j1sigb(is,ie,ib,jtot) = pe1j1sigb(is,ie,ib,jtot) + mat_pe1j1s(is,ie,ib,jb)*sigb(jb,jtot)
         pe2j0phib(is,ie,ib,jtot) = pe2j0phib(is,ie,ib,jtot) + mat_pe2j0s(is,ie,ib,jb)*phib(jb,jtot)
         pe2j0psib(is,ie,ib,jtot) = pe2j0psib(is,ie,ib,jtot) + mat_pe2j0s(is,ie,ib,jb)*psib(jb,jtot)
         pe2j1sigb(is,ie,ib,jtot) = pe2j1sigb(is,ie,ib,jtot) + mat_pe2j1s(is,ie,ib,jb)*sigb(jb,jtot)
       enddo
     enddo
     enddo
     enddo
     enddo
!
! now summ the gyro-average Maxwell fields
!
     do is = 1,ns
     do ie = 1,ne
     do iu = 1,nu
     do ib = 1,nbasis
       itot = ib+nbasis*(iu-1)+nbasis*nu*(ie-1)+nbasis*nu*ne*(is-1)
       do jtot = 1,ntot
         pu1pe1PSI(itot,jtot) = one(1,iu)*pe1j0phib(is,ie,ib,jtot)                    &
               - vs(is)*one(2,iu)*pe1j0psib(is,ie,ib,jtot)                            &
               + 2.0*(taus(is)/zs(is))*one(1,iu)*pe1j1sigb(is,ie,ib,jtot)
         pu1pe2PSI(itot,jtot) = one(1,iu)*pe2j0phib(is,ie,ib,jtot)                    &
               - vs(is)*one(2,iu)*pe2j0psib(is,ie,ib,jtot)                            &
               + 2.0*(taus(is)/zs(is))*one(1,iu)*pe2j1sigb(is,ie,ib,jtot)
         pu3pe1PSI(itot,jtot) = one(3,iu)*pe1j0phib(is,ie,ib,jtot)                    &
               - vs(is)*sqrt_two*matu(3,iu)*pe1j0psib(is,ie,ib,jtot)                  &
               + 2.0*(taus(is)/zs(is))*one(3,iu)*pe1j1sigb(is,ie,ib,jtot)
         pu2pe1PSI(itot,jtot) = one(2,iu)*pe1j0phib(is,ie,ib,jtot)                    &
               - vs(is)*sqrt_two*matu(2,iu)*pe1j0psib(is,ie,ib,jtot)                  &
               + 2.0*(taus(is)/zs(is))*one(2,iu)*pe1j1sigb(is,ie,ib,jtot)
       enddo
     enddo
     enddo
     enddo
     enddo
!
! compute mateq:  the linear GK operator tensor product
!
      do is = ns0,ns
      do js = ns0,ns
      do ie = 1,ne
      do je = 1,ne
      do iu = 1,nu
      do ju = 1,nu
      do ib = 1,nbasis
      do jb = 1,nbasis
        itot = ib+nbasis*(iu-1)+nbasis*nu*(ie-1)+nbasis*nu*ne*(is-ns0)
        jtot = jb+nbasis*(ju-1)+nbasis*nu*(je-1)+nbasis*nu*ne*(js-ns0)
        iue = iu + nu*(ie-1)
        jue = ju + nu*(je-1)
        mateq(itot,jtot) = -xi*gcut*vs(is)*k_par0*one(itot,jtot)               &
        + (-xi*ave_kpar(ib,jb)*matu(iu,ju)*one(ie,je)*vs(is)*k_par0              &
        - 2.0*w_d0*(taus(is)/zs(is))*ave_wdpar(ib,jb)*matuu(iu,ju)*one(ie,je)    &
        - w_d0*(taus(is)/zs(is))*ave_wdper(ib,jb)*mate(ie,je)*one(iu,ju)         &
        - xi*vs(is)*k_par0*ave_gradb(ib,jb)*matmirror(iue,jue)                     &
        - xi*one(1,is)*one(ib,jb)*xnue_in*(zeff_in*collision1(ie,iu,je,ju)     &
        + collision2(ie,iu,je,ju)))*one(is,js)
!        if((is.eq.1.and.js.eq.1).and.(ib.eq.1.and.jb.eq.1))then
!        nuei1 = collision1(ie,iu,je,ju)
!        write(*,*)"collision1",iue,jue,nuei1
!        nuei2 = collision2(ie,iu,je,ju)
!        write(*,*)"collision2",iue,jue,nuei2
!        endif
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!
! compute matb and mats:  the time derivative and source matricies
!
      do is = ns0,ns
      do js = ns0,ns
      do ie = 1,ne
      do je = 1,ne
      do iu = 1,nu
      do ju = 1,nu
      do ib = 1,nbasis
      do jb = 1,nbasis
        itot = ib+nbasis*(iu-1)+nbasis*nu*(ie-1)+nbasis*nu*ne*(is-ns0)
        jtot = jb+nbasis*(ju-1)+nbasis*nu*(je-1)+nbasis*nu*ne*(js-ns0)
        mats(itot,jtot) = (ky*rlns(is))*pu1pe1PSI(itot,jtot)                          &
        + (ky*rlts(is)/sqrt_two)*pu3pe1PSI(itot,jtot)                                 &
        + (ky*rlts(is))*pu1pe2PSI(itot,jtot)                                          &
        + (ky*vpar_shear_s(is)/vs(is))*pu2pe1PSI(itot,jtot)                           
 !
        matb(itot,jtot) = one(itot,jtot) - (zs(is)/taus(is))*pu1pe1PSI(itot,jtot)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!
!
!..find the eigenvalues and eigenvectors
!
!
      lwork=33*ntot
      rightvectors="N"
      if(iflux_in)rightvectors="V"
!      write(*,*)"iflux = ",iflux_in
      do i=1,ntot
      do j=1,ntot
        at(i,j) = mateq(i,j) + mats(i,j)
        bt(i,j) = matb(i,j)
        amat(i,j) = at(i,j)
        bmat(i,j) = bt(i,j)
!        write(*,*)i,j," mateq = ",mateq(i,j)
!        write(*,*)" matb = ",bt(i,j)
!        write(*,*)" mats = ",mats(i,j)
      enddo
      enddo
!      call system_clock(cpucount1)
      call zggev("N",rightvectors,ntot,at,ntot,bt,ntot,    &
                  alpha,beta,vleft,ntot,vright,ntot,work,lwork,rwork,info)
!      call system_clock(cpucount2,cpurate)
!      write(*,*)"cputime for zggev =",REAL(cpucount2-cpucount1)/REAL(cpurate)
!      write(*,*)"jmax = ",jmax,alpha(jmax)/beta(jmax)
!      write(*,*)"work(1)",work(1)
       if (info /= 0) then
         call gftm_error(1,"ERROR: ZGGEV failed in gftm_eigensolver.f90")
       endif
      if(iflux_in) then
          hetot(:,:) = vright(:,:)
      endif
!      cputime2=MPI_WTIME()
      do j1=1,ntot
        beta2=REAL(CONJG(beta(j1))*beta(j1))
        if(beta2.ne.0.0)then  ! zomega = -xi*(frequency+xi*growthrate)
!          zomega(j1)=alpha(j1)*CONJG(beta(j1))/beta2
          zomega(j1)=alpha(j1)/beta(j1)
        else
          zomega(j1)=0.0
        endif
!        write(*,*)j1,"zomega=",zomega(j1)
        rr(j1) = REAL(zomega(j1))
        ri(j1) = AIMAG(zomega(j1))
! filter out numerical instabilities that sometimes occur 
! with high mode frequency
!        if(filter_in.gt.0.0)then
 !         if(rr(j1).gt.0.0.and.ABS(ri(j1)).gt.max_freq)then
 !           rr(j1)=-rr(j1)
!            write(*,*)"filtered mode ",j1,rr(j1),ri(j1)
!            write(*,*)"beta = ",beta(j1),"alpha =",alpha(j1)
!            write(*,*)(vright(j1,j2),j2=1,ntot)
!          endif
!        endif
        do j2=1,ntot
!          if(iflux_in)then
!            vr(j1,j2) = REAL(vright(j1,j2))
!            vi(j1,j2) = AIMAG(vright(j1,j2))
!          else
            vr(j1,j2) = 0.0
            vi(j1,j2) = 0.0
!          endif
        enddo
      enddo
!
      if(ALLOCATED(rwork))DEALLOCATE(rwork)
      if(ALLOCATED(at))DEALLOCATE(at)
      if(ALLOCATED(bt))DEALLOCATE(bt)
      if(ALLOCATED(vleft))DEALLOCATE(vleft)
      if(ALLOCATED(vright))DEALLOCATE(vright)
      if(ALLOCATED(work))DEALLOCATE(work)
!
      END SUBROUTINE gftm_eigensolver
! ___________________________________________
!
      REAL FUNCTION one(i,j)
!
!     Kroneker delta function
!
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i,j
!
      if(i.eq.j)then
        one = 1.0
      else
        one = 0.0
      endif
      RETURN
!
      END FUNCTION ONE

