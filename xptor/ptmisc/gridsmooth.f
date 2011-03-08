************************************************************************
      subroutine gridsmooth
************************************************************************
cmnt  does gaussian smooth for egamma_m
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
       integer j, jj
       real*8 rho_c, weight
c
       do j=1,jout_m
        rho_c=rho(j)
        egamma_g(j)=0.D0
        weight=0.D0
        do jj=1,jout_m
         egamma_g(j)=egamma_g(j)+
     >   exp(-((rho(jj)-rho_c)/rho_g)**2.D0)*egamma_m(jj)
         weight=weight+exp(-((rho(jj)-rho_c)/rho_g)**2.D0)
        enddo
        egamma_g(j)=egamma_g(j)/weight
       enddo
c
       return
       end
