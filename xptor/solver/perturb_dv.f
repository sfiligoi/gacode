
      subroutine perturb_dv(n,ierr)
c
c     this subroutine find the eigenvalues and selected eigenvectors of amat(n,n)
c     starting from a good initial guess. 
c
c      inputs:
c        nmax dimension of amat
c        n the number of equations
c        zomega(n,n) first guess of eigenvectors
c        zevec(n) first guess of eigenvalues
c      output: 
c        zomega(n,n) new eigenvalues
c        zevec(n) new normalized eigenvectors
c        ierr error flag. 
c        if ierr is not 0 then zevec and zomega are unchanged
c       
c
      implicit none
      integer neq 
      parameter ( neq = 12 )
      COMPLEX*16 zevec(neq,neq), zevecl(neq,neq)
     &       ,zomega(neq), amat(1:neq,1:neq)
      COMMON /pert/ zevec,zevecl,zomega,amat
      INTEGER j1,j2,j3,k,n,nmax,m,ierr,i,iter,its,jm(4)
      REAL*8 tol,az(n),adz(n),relax,dzmax
      REAL*8 gmax,azmin,scale
      COMPLEX*16 xi,z(n,n),v(n,n,0:n),c(0:n),zg(n),vg(n,n),dz
      COMPLEX*16 adotv,zp, av(n)
      
c compute norms
      xi = DCMPLX(0.D0,1.D0)
      ierr = 0
      relax = 1.D-1
      tol=1.D-4
      iter = 1
      its=0
      azmin = DSQRT(DREAL(DCONJG(zomega(1))*zomega(1)))
      do j1=1,n
        az(j1) = DSQRT(DREAL(DCONJG(zomega(j1))*zomega(j1)))
        if(az(j1).lt.azmin)azmin=az(j1)
        av(j1) = DCMPLX(0.D0,0.D0)
        do j2=1,n
          av(j1) = av(j1) +
     >    zevecl(j2,j1)*zevec(j2,j1)
        enddo
      enddo
c        write(6,*)(av(j1),j1=1,n)
c initialize eigenvalues and eigenvectors
c
      do j1=1,n
        zg(j1) = zomega(j1)
        do j2=1,n
         vg(j2,j1) = zevec(j2,j1)
         v(j2,j1,0) = vg(j2,j1)
        enddo
      enddo
c compute moments
      do j1=1,n
       do m=1,n
        z(j1,m) = DCMPLX(0.D0,0.D0)
        do j2=1,n
         adotv = DCMPLX(0.D0,0.D0)
         do k=1,n
           adotv = adotv + amat(j2,k)*v(k,j1,m-1)
         enddo
         v(j2,j1,m) = xi*adotv
         z(j1,m) = z(j1,m) + zevecl(j2,j1)*v(j2,j1,m)
        enddo
         z(j1,m) = z(j1,m)/av(j1)
       enddo
      enddo
c       
c compute inital error in eigenvalues
      dzmax=0.D0
      do j1=1,n
        zp = z(j1,1)
        dz = zp - zg(j1)
        adz(j1) = DSQRT(DREAL(dz*DCONJG(dz)))
        if(dzmax.lt.adz(j1))dzmax = adz(j1)
      enddo
      dzmax = dzmax/azmin
      write(6,*)its,"  dzmax = ",dzmax
      write(6,*)(adz(j1)/az(j1),j1=1,n)
c      if(dzmax.lt.tol) go to 100
c find corrected eigenvalues using zg as a guess
      c(0) =  DCMPLX(1.D0,0.D0)
 1    its = its+1
      dzmax=0.D0
      do j1=1,n
        do m=1,n
         c(m)= z(j1,m)
        enddo
        k=0
        do j2=1,n
         if(j2.ne.j1)then
          k=k+1
          do i=n,k,-1
           c(i) = c(i)-zg(j2)*c(i-1)
          enddo
         endif  
        enddo
        zp = c(n)/c(n-1)
        dz = zp-zg(j1)
        zg(j1)=zp
        adz(j1) = DSQRT(DREAL(dz*DCONJG(dz)))
        if(dzmax.lt.adz(j1))dzmax = adz(j1)
c        scale=relax
c       if(adz(j1)/az(j1).lt.tol)scale = 0.D0
c        if(j1.gt.2)scale=0.D0
c        zg(j1) = zg(j1) + scale*(zp-zg(j1))
        az(j1) =DSQRT(REAL(DCONJG(zg(j1))*zg(j1)))
      enddo
      dzmax = dzmax/azmin
      if(dzmax.lt.tol) go to 100
      write(6,*)its,"  dzmax = ",dzmax
      write(6,*)(adz(j1)/az(j1),j1=1,n)
      if(its.lt.iter)go to 1
c maximum number of iterations exceeded without convergence
        ierr = 2
        return
 100  continue 
c      find four most unstable eigenvalues
      do j2=1,4
        jm(j2)=0
      enddo
      do j2=1,4
        gmax = 0.D0
        do j1=1,n
         if(j1.ne.jm(1).and.j1.ne.jm(2).and.j1.ne.jm(3))then
          if(DIMAG(zg(j1)).gt.gmax)then
           jm(j2) = j1
           gmax = DIMAG(zg(j1))
          endif
         endif
        enddo
      enddo
c compute four selected eigenvectors
      do j3=1,n          
       c(0) = DCMPLX(1.D0,0.D0)
       do j1=1,4
        do m=1,n
         c(m)= v(j3,jm(j1),m)
        enddo
        k=0
        do j2=1,n
         if(j2.ne.jm(j1))then
          k=k+1
          do i=n-1,k,-1
            c(i) = c(i)-zg(j2)*c(i-1)
          enddo
         endif
        enddo
        vg(j3,jm(j1)) = c(n-1)
       enddo
      enddo
c write results to output 
      do j1=1,n
        zomega(j1) = zg(j1)
        do j2=1,n
          zevec(j2,j1) = vg(j2,j1)
        enddo
      enddo
      return
      end
            
          
        
      
