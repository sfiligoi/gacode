  MODULE eq_pfprim
     USE param,                 ONLY : kj
     USE rhog,                  ONLY : press,psir
     USE numbrs,                ONLY : nj
     USE soln,                  ONLY : curden
     USE constnts,              ONLY : u0
     USE geom,                  ONLY : r2cap
     USE machin,                ONLY : rmajor,btor
     USE io,                    ONLY : ncrt,nout
     USE mhdcom,                ONLY : rma

     IMPLICIT NONE
    
     CONTAINS 
       SUBROUTINE nwg_pfprim(psi_eqdsk,sf,sp,xpp,xffp,nw,psiaxis,&
                             psilim,set_grid)
!------------------------------------------------------------------
! subroutine  gets pprim,ffprim directly on equispaced psi grid
!------------------------------------------------------------------

       INTEGER i,isignn,nw,set_grid
       REAL*8 psi_eqdsk(nw),sf(nw),sp(nw),xpp(nw),xffp(nw),psiaxis,psilim
       REAL*8 curden_loc(nw),r2cap_loc(nw),dpsi,curaxis


!
! --- construct uniform, nw in length, psi grid for eqdsk
!
      IF(set_grid == 1)THEN
         psi_eqdsk(1)   =  psiaxis
         DO i=2,nw-1
            psi_eqdsk(i) = (psiaxis-(i-1)*(psiaxis-psilim)/(nw-1))
         END DO
         psi_eqdsk(nw)  =  psilim

         RETURN
       ENDIF

!
!         given press on psir  grid, get sp on psi_eqdsk grid, etc. ...
          CALL intrp (1, 1, psir,  press,  nj, psi_eqdsk, sp, nw)
          CALL intrp (1, 1, psir,  curden, nj, psi_eqdsk,curden_loc, nw)
          CALL intrp (1, 1, psir,  r2cap, nj, psi_eqdsk,r2cap_loc, nw)


          CALL difydx (psi_eqdsk,sp,xpp, nw)


          isignn = 1
          if (btor .lt. 0)  isignn = -1


          DO  i=1,nw
              xffp(i) = -u0*rmajor*(curden_loc(i)+rmajor*xpp(i)) / r2cap_loc(i)
          ENDDO
 
          sf(nw) = (rmajor*btor)**2
          DO i = nw-1,1,-1
             dpsi  = psi_eqdsk(i)-psi_eqdsk(i+1)
             sf(i) = sf(i+1)+dpsi*(xffp(i)+xffp(i+1))
          ENDDO
          DO i=1,nw
             IF (sf(i) .ge. 0.0)THEN
                sf(i) = isignn * SQRT (sf(i))
             ELSE
                IF (ncrt .ne. 0)  write (ncrt, 1)
                IF (nout .ne. 0)  write (nout, 1)
             ENDIF
          END DO
 
!          print *,'in nwg_pfprim new xffp(1:nw) = ',xffp(1:nw)
!          print *,'in nwg_pfprim xpp =',xpp(1:nw)
!          print *,'in nwg_pfprim sf =',sf(1:nw)
!          print *,'in nwg_pfprim sp =',sp(1:nw)
!          curaxis = -rma *xpp(1) - xffp(1)/(u0*rma)   ! 888888999999
!          print *,'in nwg_pfprim curaxis  =',curaxis
          RETURN

  1     format (/   ' subroutine EQ_PFPRIM reports:',/,                  &
                    '   write of EQDSK file was abandoned' ,/,           &
                    '   since code could not define a proper f(psi)', /, &
                    '   otherwise ONETWO execution will continue',    /)

       END SUBROUTINE nwg_pfprim


  END MODULE eq_pfprim
