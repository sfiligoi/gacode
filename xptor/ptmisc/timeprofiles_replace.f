       subroutine timeprofiles_replace(ktime,ktime_mini,ktime_chk,Tvec)
************************************************************************
cmnt   fills  time profiles
************************************************************************
       include '../inc/input.m'
       include '../inc/tport.m'
       include '../inc/model.m'
       include '../inc/ptor.m'
c
       integer ktime,ktime_mini,ktime_chk
       real*8 aktime_mini,aktime_chk
       integer k
       real*8 Tvec(mxfields,mxgrid)
c
       if(ktime.eq.0) return
c
      aktime_mini=ktime_mini+1
      aktime_chk=ktime_chk
      if(ktime_chk.eq.0) aktime_chk=1.D0
      if(ktime_chk.eq.0) aktime_mini=1.D0
c
      if(xte_t(0,ktime).ne.0.) then
      do k=1,ngrid
       Tvec(2,k)=xte_t(k,ktime)+
     >    (xte_t(k,ktime+1)-xte_t(k,ktime))*aktime_mini/aktime_chk
      enddo 
      endif
c
      if(xti_t(0,ktime).ne.0.) then
      do k=1,ngrid
       Tvec(3,k)=xti_t(k,ktime)+
     >    (xti_t(k,ktime+1)-xti_t(k,ktime))*aktime_mini/aktime_chk
      enddo
      endif 
c
       return
       end 
