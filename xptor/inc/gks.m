c@gks.m 27-May-04 J. Kinsey
c common block for GKS growth rates, frequencies
cc
      integer jrmx, izrmx
      parameter (izrmx=20,jrmx=100)
      real*8 xp_gk(jrmx),yp_gk(izrmx,jrmx,3),
     &       zp_gk(jrmx)
c
      common /gks/ xp_gk, yp_gk, zp_gk
