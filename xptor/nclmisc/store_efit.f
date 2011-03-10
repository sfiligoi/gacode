      subroutine store_efit(nr_er,k_edotb_e,k_potato_e,
     &                 a0_e,r0_e,bt0_e,
     &                 izim1_e,izim2_e, amuim1_e,amuim2_e,
     &                 izi0_e,amui0_e,rhop_er,rhot_er,
     &                 f_er,q_er,xj_er,xj_nb_ex_er,
     &                 rin_er,rout_er,elong_er,vol_er,vp_er,
     &                 phit_er,fm_er,grth_er,gph_er,gth_er,
     &                 b2_er,bm2_er,grho1_er,grho2_er,
     &                 gr2bm2_er,ftrap_er,fhat_er,
     &                 bpout_er,btout_er,psi_er,psip_er,
     &                 rm1_er,rm2_er)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c  nr_r-number of radial points
c  rhot_r(i)-normalized toroidal flux grid proportional to (Phi)**0.5
c  rhop_r-normalized poloidal flux grid proportional to psi
c  f_r(i)-2*pi*R*B_t/mu0 (A)
c  q_r(i)-safety factor
c  rhot_r(i)-normalized toroidal flux grid proportional to (Phi)**0.5
c  rin_r(i)-major radius grid on inside of torus in axis plane (m)
c  rout_r(i)-major radius grid on outside of torus in axis plane (m)
c  elong_r(i)-elongation
c  vol_r(i)-volume enclosed (m**3)
c  vp_r(i)-d vol_r/d rhot_r/a0 (m**2)
c  phit_r(i)-toroidal flux (Wb)
c  fm_r(3,i)-geometric factor

c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      implicit none
c
c ... already in common in efitdat.m:
c     rhop_r, rhot_r, q_r, rin_r, rout_r, elong_r, vol_r, phit_r
c     grho1_r, grho2_r, bpout_r, btout_r, psi_r
c ... need to add to common in efitdat.m:
c     vp_r, fm_r, grth_r, gph_r, gth_r, b2_r, bm2_r, gr2bm2_r
c     ftrap_r, fhat_r
c
      include '../inc/efitdat.m'
c
      integer j, nr_er, izim1_e, izim2_e,
     &     k_edotb_e, k_potato_e
      integer izi0_e(*)
      real a0_e, r0_e, bt0_e
      real amui0_e(*), amuim1_e, amuim2_e
      real rhop_er(*), rhot_er(*),
     &     f_er(*), q_er(*), xj_er(*), xj_nb_ex_er(*),
     &     rin_er(*), rout_er(*), elong_er(*), vol_er(*),
     &     vp_er(*), phit_er(*), fm_er(3,*),
     &     grth_er(*), gph_er(*), gth_er(*), b2_er(*), bm2_er(*),
     &     grho1_er(*), grho2_er(*), gr2bm2_er(*), ftrap_er(*),
     &     fhat_er(*), bpout_er(*), btout_er(*), psi_er(*),
     &     psip_er(*), rm1_er(*),rm2_er(*)
c
       nr_r=nr_er
       a0=a0_e
       r0=r0_e
       bt0=bt0_e
       izim1=izim1_e
       izim2=izim2_e
       amuim1=amuim1_e
       amuim2=amuim2_e
       k_edotb=k_edotb_e
       k_potato=k_potato_e
       do j=1,mx_ni
         izi0(j)=izi0_e(j)
         amui0(j)=amui0_e(j)
       enddo
       do j=1,nr_er
         rhop_r(j)=rhop_er(j)
         f_r(j)=f_er(j)
         q_r(j)=q_er(j)
         xj_r(j)=xj_er(j)
         xj_nb_ex_r(j)=xj_nb_ex_er(j)
         rhot_r(j)=rhot_er(j)
         rin_r(j)=rin_er(j)
         rout_r(j)=rout_er(j)
         elong_r(j)=elong_er(j)
         vol_r(j)=vol_er(j)
         vp_r(j)=vp_er(j)
         phit_r(j)=phit_er(j)
         grth_r(j)=grth_er(j)
         gph_r(j)=gph_er(j)
         gth_r(j)=gth_er(j)
         grho1_r(j)=grho1_er(j)
         grho2_r(j)=grho2_er(j)
         gr2bm2_r(j)=gr2bm2_er(j)
         b2_r(j)=b2_er(j)
         bm2_r(j)=bm2_er(j)
         fhat_r(j)=fhat_er(j)
         ftrap_r(j)=ftrap_er(j)
         fm_r(1,j)=fm_er(1,j)
         fm_r(2,j)=fm_er(2,j)
         fm_r(3,j)=fm_er(3,j)
         bpout_r(j)=bpout_er(j)
         btout_r(j)=btout_er(j)
         psi_r(j)=psi_er(j)
         psip_r(j)=psip_er(j)
c         write(*,100) j, rhot_r(j),xj_r(j),'store_efit'
       enddo
       q0_efit=q_r(1)
       el0=elong_r(1)
c
 100  format(i2,2x,0p1f8.6,1pe13.6,a12)
       end
