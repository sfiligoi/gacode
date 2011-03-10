      subroutine delay
************************************************************************
cmnt delays egamma_m by i_delay time steps (i_delay < or = 10) 
************************************************************************
       implicit none
c
       include '../inc/input.m'
       include '../inc/tport.m'
       include '../inc/model.m'
       include '../inc/data.m'
       include '../inc/share.m'
       include '../inc/sharegk.m'
       include '../inc/glf.m'
c
       integer ii 
c
       if(i_delay.eq.0) return
c
       egamma_d(jm,i_delay)=egamma_m(jm)
c  
       if(i_delay.eq.1) return
c
       do ii=1,i_delay-1
       egamma_d(jm,ii)=egamma_d(jm,ii+1)
       enddo
c  
       egamma_d(jm,i_delay)=egamma_m(jm)
c
       return
       end
