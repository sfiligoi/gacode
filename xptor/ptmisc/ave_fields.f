ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ave_fields
c
c     average the  fields
c
      implicit none
c
      include '../inc/input.m'
      include '../inc/ptor.m'
      include '../inc/tport.m'
      integer i
      do i=1,ngrid-1
         if(itport_pt(1).ne.0)
     >   ne_m(i) = (1.D0-ave_dv)*ne_m(i)
     >    +ave_dv*(ne_m(i-1)+ne_m(i+1))/2.D0
         if(itport_pt(2).ne.0)
     >   te_m(i) =(1.D0-ave_dv) *te_m(i)
     >    +ave_dv*(te_m(i-1)+te_m(i+1))/2.D0
         if(itport_pt(3).ne.0)
     >   ti_m(i) = (1.D0-ave_dv)*ti_m(i)
     >     +ave_dv*(ti_m(i-1)+ti_m(i+1))/2.D0
         if(itport_pt(4).ne.0)
     >   vphi_m(i) = (1.D0-ave_dv)*vphi_m(i)
     >     +ave_dv*(vphi_m(i-1)+vphi_m(i+1))/2.D0
c         if(itport_pt(5).ne.0)
c     >   vper_m(i) = (1.D0-ave_dv)*vper_m(i)
c     >     +ave_dv*(vper_m(i-1)+vper_m(i+1))/2.D0
      enddo
      do i=1,ngrid-1
        ni_m(i) = fi_m(i)*ne_m(i)
      enddo
      ni_m(0)=ni_m(1)
      ti_m(0)=ti_m(1)
      te_m(0)=te_m(1)
      vphi_m(0)=vphi_m(1)
c      vper_m(0)=vper_m(1)
      return
      end
