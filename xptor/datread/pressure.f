c@pressure.f
c jek 20-Jan-11
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c... Compute pressures using _d data
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c  
      subroutine pressure
c
      implicit none
      include 'mpif.h'
      include '../inc/data.m'
      include '../inc/input.m'
      include '../inc/glf.m'
c
      integer j
c
c... If total pressure ptot read, then construct pfast
c
      if(iptotr.eq.1) then
        do j=1,nj_d
           pfast_d(j)=ptot_d(j)-ene_d(j)*te_d(j)-en_d(j,1)*ti_d(j)
           if(ipfst.eq.0) pfast_d(j)=0.D0
        enddo
c
c... If not, then compute total thermal pressure instead
c
      elseif(iptotr.eq.0) then
         if(i_proc.eq.0)
     >   write(*,'(a32)') 'Computing total thermal pressure'
         do j=1,nj_d
            ptot_d(j) = ene_d(j)*(te_d(j) + ti_d(j))
            pfast_d(j) = 0.D0
         enddo
         if(ismooth_all.ne.0)call average7_1d(ptot_d,nj_d)
      endif
c
      end
