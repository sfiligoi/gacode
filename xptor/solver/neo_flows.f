      SUBROUTINE neo_flows(kmax,upol_neo,udia_neo)
c
      IMPLICIT NONE
c
      include 'mpif.h'
      INCLUDE '../inc/tport.m'
      INCLUDE '../inc/glf.m'
      INCLUDE '../inc/ptor.m'
c
      INTEGER kgrid
      REAL*8 upol_neo(3,0:jmaxmt),udia_neo(3,0:jmaxmt)
c
      REAl*8 upol_sum(3,0:jmaxmt),udia_sum(3,0:jmaxmt)
      INTEGER  i,k,kmax
      INTEGER MPI_status(MPI_STATUS_SIZE)
c
      do i=1,nspecies
        do k=0,jmaxmt
          upol_neo(i,k)=0.0
          udia_neo(i,k)=0.0
          upol_sum(i,k)=0.0
          udia_sum(i,k)=0.0
        enddo
      enddo
      do k=1+i_proc,kmax-1,n_proc
        nem = (ne_m(k+1)+ne_m(k))/2.D0
        tim = (ti_m(k+1)+ti_m(k))/2.D0
        tem = (te_m(k+1)+te_m(k))/2.D0
        fim = (fi_m(k+1)+fi_m(k))/2.D0
        fzm = (fz_m(k+1)+fz_m(k))/2.D0
        nim = fim*nem
        nzm = fzm*nem
        vexbm= 0.0
        vpolm= 0.0
        gradnem = (ne_m(k+1)-ne_m(k))/dr(k,2)
        gradtim = (ti_m(k+1)-ti_m(k))/dr(k,2)
        gradtem = (te_m(k+1)-te_m(k))/dr(k,2)
        gradvexbm = 0.0
        gradvpolm = 0.0
        gradfim = (fi_m(k+1)-fi_m(k))/dr(k,2)
        gradfzm = (fz_m(k+1)-fz_m(k))/dr(k,2)
        gradnim = fim*gradnem + nem*gradfim
        gradnzm = fzm*gradnem + nem*gradfzm
        jm=k
        call neoclassical
        do i=1,nspecies
c these are on the 1/2 grid
          upol_neo(i,k) = vneo(i)
          udia_neo(i,k) = vdia(i)
        enddo
      enddo
c
      call MPI_BARRIER(MPI_COMM_WORLD,i_err)
      call MPI_REDUCE(upol_neo,upol_sum,3*(jmaxmt+1)
     >  ,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,i_err)
      call MPI_BCAST(upol_sum,3*(jmaxmt+1)
     >  ,MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(udia_neo,udia_sum,3*(jmaxmt+1)
     >  ,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,i_err)
      call MPI_BCAST(udia_sum,3*(jmaxmt+1)
     >  ,MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, i_err)
c
      do i=1,nspecies
c use zero gradient condition at boundaries
        upol_sum(i,kmax) = upol_sum(i,kmax-1)
        udia_sum(i,kmax) = udia_sum(i,kmax-1)
        upol_sum(i,0) = upol_sum(i,1)
        udia_sum(i,0) = udia_sum(i,1)
      enddo
c
      do i=1,nspecies
        do k=1,kmax-1
c interpolate onto the full grid
          upol_neo(i,k) = 0.5*(upol_sum(i,k)+upol_sum(i,k-1))
          udia_neo(i,k) = 0.5*(udia_sum(i,k)+udia_sum(i,k-1))
        enddo      
        upol_neo(i,kmax) = upol_neo(i,kmax-1)
        udia_neo(i,kmax) = udia_neo(i,kmax-1)
        upol_neo(i,0) = upol_neo(i,1)
        udia_neo(i,0) = udia_neo(i,1)
      enddo
c
      RETURN
      END ! subroutine neo_flows
