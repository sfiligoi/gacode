
          program read_scratch_file
!reads scratch1 file output by onetwo
!and does some processing of contours:

          IMPLICIT NONE
          INTEGER,PARAMETER ::         N=1000
          INTEGER,PARAMETER ::         NW=129

          CHARACTER*4 CDEV,CDSP,C*1,filetype
          REAL*8 xp(n,nw),yp(n,nw),arclen(n,nw),bpinv(n,nw), &
                 rmin(nw),rmax(nw),zmin(nw),zmax(nw),rzmax(nw), &
                 rzmin(nw),zrmax(nw),zrmin(nw),psivlcpy(nw), &
                 rho(nw),psirr(nw),xmagax,ymagax
          INTEGER mpp(nw),ndisk,j,i,k,npsi,mp,nj
!          CALL METAFL('CONS')
          filetype ='EPS' ! only one page in eps mode
!          filetype='PS'
          CALL METAFL(filetype)
          ndisk = 10

          open (unit = ndisk, file = 'scratch1', status = 'UNKNOWN')
          j=0
 100      read (ndisk, 8100)  mp
          if (mp .lt. n  .and. mp .gt. 0)then
                    j = j+1
                     read (ndisk, 8120)  (xp(i,j),yp(i,j), i=1,mp)
                     read (ndisk, 8120)  (bpinv(i,j),i=1,mp)
                     read (ndisk, 8120)  (arclen(i,j),i=1,mp)
                     mpp(j) = mp
cray204.f8100                 format (i6,2x,i6)
8120                 format (6e12.5)
                     rmin(j) = xp(1,j)
                     rmax(j) = xp(1,j)
                     zmin(j) = yp(1,j)
                     zmax(j) = yp(1,j)
                     do k=1,mp
                        rmin(j) = MIN(xp(k,j),rmin(j))
                        rmax(j) = MAX(xp(k,j),rmax(j))
                        if(rmax(j) .eq. xp(k,j))zrmax(j) = yp(k,j)
                        if(rmin(j) .eq. xp(k,j))zrmin(j) = yp(k,j)
                        zmin(j) = MIN(yp(k,j),zmin(j))
                        zmax(j) = MAX(yp(k,j),zmax(j))
                        if(zmax(j) .eq. yp(k,j))rzmax(j) = xp(k,j)
                        if(zmin(j) .eq. yp(k,j))rzmin(j) = xp(k,j)
                     enddo
                     go to 100
           else
              !read remaining info
                     read (ndisk, 8100)  npsi,nj
                     if(npsi .gt. nw)then
                        print *,'error, recompile with nw .gt. ',npsi
                        call exit(1)
                     endif
                     read (ndisk, 8120)  (psivlcpy(i),i=1,npsi)
                     read(ndisk,  8120) xmagax,ymagax
                     read (ndisk, 8120)  (rho(i),i=1,npsi)
                     read (ndisk, 8120)  (psirr(j),j=1,nj)
!                    call close(ndisk)


                     call process(xp,yp,arclen,rmin,rmax,zmin,zmax, &
                               rzmax,rzmin,zrmax,zrmin,rho,psirr, &
                           psivlcpy,xmagax,ymagax,n,mpp,nw,npsi,nj,filetype)
                     call exit(1)

            endif


            END





         subroutine process(r,z,s,rmin,rmax,zmin,zmax,rzmax, &
                   rzmin,zrmax,zrmin,rho,psirr, &
                   psivlcpy,xmagax,ymagax,n,mpp,nw,npsi,nj,filetype)


         USE smthspline,               ONLY : smooth_values



         IMPLICIT NONE
         CHARACTER*4 filetype
         INTEGER n,mp,nw,mpp(nw),npsi,nj,j,ir,il,i,no_weight,ibc
         REAL*8 R(n,nw),z(n,nw),s(n,nw),rmin(nw),rmax(nw),            &
                    rzmax(nw),rzmin(nw),zrmax(nw),zrmin(nw),          &
                    zmax(nw),zmin(nw),psivlcpy(npsi),                 &
                     triang(npsi),elong(npsi),rho(npsi),psirr(nw),    &
                     pindentnpsi(nw),rminavnpsi(nw),deriv1(nw),       &
                     deriv2(nw),weight(nw),derivl,derivr,             &
                     triangsmth(nw),elongsmth(nw),pindentsmth(nw)
         REAL*8 psimin,psimax,triangmin,triangmax,elongmax,elongmin,  &
                pindentmax,pindentmin,rzmb,rzms,zms,zmb,rsc,zsc,      &
                zmb1,zmb2,zms1,zms2,xmagax,ymagax,armin,almin,advert, &
                xstep,ystep,dtriang,dpsi,delong,dpindent, xl,xr,xs,   &
                yb,yt,ys,sdevfctr

         psimin = 1.d35
         psimax = -1.d35
         triangmin = 1.d35
         triangmax = -1.d35
         elongmax = -1.d35
         elongmin = 1.d35
         pindentmin = 1.d35
         pindentmax = -1.d35
         do j=1,npsi - 1
            almin = 1.d12
            armin = almin

            psimin = MIN(psimin,psivlcpy(j))
            psimax = MAX(psimax,psivlcpy(j))
            ir = -1
            il = -1
            do i=1,mpp(j)-1
              advert = ABS (ymagax-Z(i,j))
              if (R(i,j) .le. xmagax) then
                almin = MIN (almin, advert)
                if (almin .eq. advert )il = i
              else
                armin = MIN (armin, advert)
                if ( armin  .eq.  advert)ir = i
              end if
            end do

            triang(j) = (0.5*(rmax(j)+rmin(j)) - rzmax(j))/ &
                        (0.5*(R(ir,j)-R(il,j)))
            triangmin = MIN(triangmin,triang(j))
            triangmax = MAX(triangmax,triang(j))
            elong(j) = (zmax(j)-zmin(j))/(rmax(j)-rmin(j))
            elongmax = MAX(elong(j),elongmax)
            elongmin = MIN(elong(j),elongmin)



!            rmajavnpsi (j) =  0.5 * (xp(ir) + xp(il))
            rminavnpsi (j) =  0.5 * (R(ir,j) - R(il,j))
            ! note multiply by 100 to get scale for plot:
            pindentnpsi(j) = 1000.*(R(il,j) - rmin(j)) / (2.0 * rminavnpsi(j))
            pindentmin = MIN(pindentnpsi(j),pindentmin)
            pindentmax = MAX(pindentnpsi(j),pindentmax)
         enddo
         no_weight = 1
         ibc = 2
         sdevfctr = 0.0005
         Call smooth_values(triang,psivlcpy,npsi-1,triangsmth,deriv1,deriv2, &
         weight,no_weight,sdevfctr,ibc,derivl,derivr)

         Call smooth_values(elong,psivlcpy,npsi-1,elongsmth,deriv1,deriv2, &
         weight,no_weight,sdevfctr,ibc,derivl,derivr)
         sdevfctr = 0.1
         Call smooth_values(pindentnpsi,psivlcpy,npsi-1,pindentsmth,deriv1,deriv2, &
         weight,no_weight,sdevfctr,ibc,derivl,derivr)

      CALL DISINI
      CALL PAGERA
      CALL COMPLX


      go to 50

      CALL AXSPOS(450,1800)
      CALL AXSLEN(800,1700)

      CALL NAME('R,cm ','X')
      CALL NAME('Z,cm ','Y')

      CALL LABDIG(-1,'X')
      CALL TICKS(10,'XY')

!      CALL TITLIN('Demonstration of CURVE',1)
      CALL TITLIN('  ',3)
      xl = 100.d0 ; xr = 250.d0 ;xs = xl
      xstep = (250.-100.)/5.
      ystep = (150. -(-150.))/5.
      yb = -125.d0 ; yt = 125.d0 ; ys = yb
      CALL GRAF(xl,xr,xs,xstep,yb,yt,ys,ystep)
      CALL TITLE
      CALL COLOR('GREEN')
      do j=1,nw,6
          CALL CURVE(R(1,j),z(1,j),mpp(j))
      enddo
      CALL COLOR('WHITE')
      CALL CURVE(rzmax,zmax,nw-1)
      CALL CURVE(rzmin,zmin,nw-1)
      CALL CURVE(rmax,zrmax,nw-1)
      CALL CURVE(rmin,zrmin,nw-1)
      CALL ENDGRF

   
!--------------------------------------------
      CALL AXSPOS(1650,1800)
      CALL AXSLEN(1000,700)
      CALL NAME('R at zmax,zmin','X')
      CALL NAME('Zmax,Zmin','Y')
      CALL LABDIG(-1,'X')
      CALL TICKS(2,'XY')
      CALL TITLIN('  ',3)
      rzmb = rzmax(1)
      rzms =rzmb
      zms = zmin(1)
      zmb = zms
      do j =1,nw-1
         rzmb = MAX(rzmax(j),rzmb)
         rzms = MIN(rzmax(j),rzms)
         zmb = MAX(zmax(j),zmb)
         zms = MIN(zmax(j),zms)
      enddo
      rsc = (rzmb-rzms)/4.
      zsc = (zmb-zms)/10.
      CALL GRAF(rzms,rzmb,rzms,rsc,zms,zmb,zms,zsc)
      CALL TITLE

      CALL CURVE(rzmax,zmax,nw-1)
      CALL CURVE(rzmin,zmin,nw-1)
      CALL ENDGRF




      CALL AXSPOS(1650,900)
      CALL AXSLEN(1000,700)
      CALL NAME('Rmin,Rmax','X')
      CALL NAME('z at Rmin,Rmax','Y')
      CALL LABDIG(-1,'X')
      CALL TICKS(2,'XY')
      CALL TITLIN('  ',3)
      rzmb = -1.
      rzms =1.d35
      zms1 = 1.d35
      zmb1 = -1.d35
      zmb2 = zmb1
      zms2 = zms1
      do j =1,nw-1
         rzmb = MAX(rmax(j),rzmb)
         rzms = MIN(rmin(j),rzms)
         zmb1 = MAX(zrmax(j),zmb1)
         zmb2 = MAX(zrmin(j),zmb2)
         zms1 = MIN(zrmin(j),zms1)
         zms2 = MIN(zrmax(j),zms2)
      enddo 
      zmb = MAX(zmb1,zmb2)
      zms = MIN(zms1,zms2)
      rsc = (rzmb-rzms)/4.
      zsc = (zmb-zms)/10.
      CALL GRAF(rzms,rzmb,rzms,rsc,zms,zmb,zms,zsc)
      CALL TITLE
      CALL CURVE(rmax,zrmax,nw-1)
      CALL CURVE(rmin,zrmin,nw-1)
      CALL ENDGRF
   
       IF(filetype =='EPS')THEN
           CALL DISFIN
           CALL SYMFIL ('PSCi', 'KEEP')
           return
       ENDIF



 50 continue
!---------------------------------------------------

!      CALL COLOR('FORE')
!      CALL DASH
!      CALL XAXGIT



      call NEWPAG


      CALL AXSPOS(450,800)
      CALL AXSLEN(1000,700)
      CALL NAME('Psi','X')
      CALL NAME('triangularity','Y')
      CALL LABDIG(4,'X')
      CALL TICKS(2,'XY')
      CALL TITLIN('  ',3)

      dtriang = (triangmax-triangmin)/5.
      dpsi = (psimax - psimin)/3.
      CALL GRAF(psimin,psimax,psimin,dpsi,triangmin, &
                triangmax,triangmin,dtriang)
      CALL TITLE
      CALL CURVE(psivlcpy,triang,npsi-1)
      CALL COLOR('RED')
      CALL CURVE(psivlcpy,triangsmth,npsi-1)
      CALL COLOR('WHITE')
      CALL ENDGRF


      CALL COLOR('WHITE')
      CALL AXSPOS(1750,800)
      CALL AXSLEN(1000,700)
      CALL NAME('Psi','X')
      CALL NAME('elongation','Y')
      CALL LABDIG(4,'X')
      CALL TICKS(2,'XY')
      CALL TITLIN('  ',3)

      delong = (elongmax-elongmin)/5.
      dpsi = (psimax - psimin)/3.
      CALL GRAF(psimin,psimax,psimin,dpsi,elongmin, &
                elongmax,elongmin,delong)
      CALL TITLE
      CALL CURVE(psivlcpy,elong,npsi-1)
      CALL COLOR('RED')
      CALL CURVE(psivlcpy,elongsmth,npsi-1)
      CALL COLOR('WHITE')
      do j =1,npsi
        print *,j,elong(j),elongsmth(j)
      enddo
      CALL ENDGRF


      CALL AXSPOS(450,1900)
      CALL AXSLEN(1000,700)
      CALL NAME('Psi','X')
      CALL NAME('indentation * 1000','Y')
      CALL LABDIG(4,'X')
      CALL LABDIG(3,'Y')
      CALL TICKS(2,'XY')
      CALL TITLIN('  ',3)

      dpindent = (pindentmax-pindentmin)/5.
      dpsi = (psimax - psimin)/3.
      CALL GRAF(psimin,psimax,psimin,dpsi,pindentmin, &
                pindentmax,pindentmin,dpindent)
      CALL TITLE
      CALL CURVE(psivlcpy,pindentnpsi,npsi-1)
      CALL COLOR('RED')
      CALL CURVE(psivlcpy,pindentsmth,npsi-1)
      CALL COLOR('WHITE')
      CALL ENDGRF




      CALL AXSPOS(1950,1900)
      CALL AXSLEN(600,800)

      CALL NAME('R,cm ','X')
      CALL NAME('Z,cm ','Y')

      CALL LABDIG(-1,'X')
      CALL LABDIG(-1,'Y')
      CALL TICKS(10,'XY')

!      CALL TITLIN('Demonstration of CURVE',1)
      CALL TITLIN('  ',3)
      xl = 100.d0 ; xr = 250.d0 ;xs = xl
      xstep = (250.-100.)/2.
      ystep = (150. -(-150.))/2.
      yb = -125.d0 ; yt = 125.d0 ; ys = yb
      CALL GRAF(xl,xr,xs,xstep,yb,yt,ys,ystep)
      CALL TITLE
      CALL COLOR('GREEN')
      do j=1,nw,6
          CALL CURVE(R(1,j),z(1,j),mpp(j))
      enddo
      CALL COLOR('WHITE')
      CALL CURVE(rzmax,zmax,nw-1)
      CALL CURVE(rzmin,zmin,nw-1)
      CALL CURVE(rmax,zrmax,nw-1)
      CALL CURVE(rmin,zrmin,nw-1)
      CALL ENDGRF
 


      CALL DISFIN
      CALL SYMFIL ('PSCi', 'KEEP')
      return
      END

 
