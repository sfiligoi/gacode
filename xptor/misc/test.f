       program test
c
       implicit none
       character*132 headerline, msg
       character*2 stflg
       character*(2) namep_d(3), namei_d(3), namen_d(3)
       integer niterdb, nion_d, nprim_d, nimp_d, nneu_d,
     &         ibion_d, ishot_d
       integer i, j, nj, njj, nj_d, jmaxmt, jmaxm, nex, iflag
       real*8 aj, time_d, rgeom_d, rmag_d, rmajor_d, kappa_d,
     &       deltao_d, pindento_d, volo_d, areao_d, btor_d,
     &       tocur_d, totohm_d, totboot_d, totbeam_d, totrf_d,
     &       betap_d, beta_d, ali_d, te0_d, ti0_d
c
       parameter(niterdb=11, njj=51, nj=301, jmaxmt=300, jmaxm=80)
       real*8 blank_d(nj), x_d(nj), r_d(nj), rho_d(nj), r_new(nj)
       real*8 rho(0:80)
c
       open(unit=niterdb,status='old',access='sequential',
     &      file='iterdb.89116')
c
       read(niterdb,'(a)')headerline
       read(niterdb,'(a)') stflg
       read(niterdb,*) ishot_d
       read(niterdb,'(a)') stflg
       read(niterdb,*) nj_d
       read(niterdb,'(a)')stflg
       read(niterdb,'5(2x,i6)')nion_d
       read(niterdb,'(a)')stflg
       read(niterdb,'5(2x,i6)')nprim_d
       read(niterdb,'(a)')stflg
       read(niterdb,'5(2x,i6)')nimp_d
       read(niterdb,'(a)')stflg
       read(niterdb,'5(2x,i6)')nneu_d
       read(niterdb,'(a)')stflg
       read(niterdb,'5(2x,i6)')ibion_d
       read(niterdb,'(a)')stflg
       read(niterdb,'5(2x,a)')(namep_d(i),i=1,nprim_d)
       read(niterdb,'(a)')stflg
       read(niterdb,*)(namei_d(i),i=1,nimp_d)
       read(niterdb,'(a)')stflg
       read(niterdb,*)(namen_d(i),i=1,nneu_d)
       read(niterdb,'(a)')stflg
       read(niterdb,10)time_d

c Rgeom : major radius of geometric
         read(niterdb,'(a)')stflg
         read(niterdb,10)rgeom_d                   !meters

c Rmag :  major radius of mag axis, meters
         read(niterdb,'(a)')stflg
         read(niterdb,10)rmag_d                    !meters

c R0 : major radius of vaccuum btor ref
         read(niterdb,'(a)')stflg
         read(niterdb,10)rmajor_d                  !meters

c kappa : plasma elongation
         read(niterdb,'(a)')stflg
         read(niterdb,10)kappa_d

c delta : plasma triangularity
         read(niterdb,'(a)')stflg
         read(niterdb,10)deltao_d

c pindent : plasma indentation
         read(niterdb,'(a)')stflg
         read(niterdb,10)pindento_d

c volo : plasma volume,meters**3
         read(niterdb,'(a)')stflg
         read(niterdb,10)volo_d                    !meters**3

c cxareao :plasma cross sectional area, meters**2
         read(niterdb,'(a)')stflg
         read(niterdb,10)areao_d                   !meters**2

c Btor : vaccuum toroidal field at rmajor, tesla
         read(niterdb,'(a)')stflg
         read(niterdb,10)btor_d                    !tesla

c total,ohmic,bootstrap,beam,and rf currents, amps
         read(niterdb,'(a)')stflg
         read(niterdb,10)tocur_d,totohm_d,totboot_d,totbeam_d,totrf_d !amps

c betap : poloidal beta
         read(niterdb,'(a)')stflg
         read(niterdb,10)betap_d

c beta : toroidal beta
         read(niterdb,'(a)')stflg
         read(niterdb,10)beta_d

c ali : plasma inductance
         read(niterdb,'(a)')stflg
         read(niterdb,10)ali_d

c te0 : central electron temperature
         read(niterdb,'(a)')stflg
         read(niterdb,10)te0_d                        !kev

c ti0 : central ion temperature
         read(niterdb,'(a)')stflg
         read(niterdb,10)ti0_d                        !kev

c---psi on rho grid,volt*sec/rad

          read(niterdb,'(a)')stflg
cx          read(niterdb,10)(psir_d(j),j=1,nj_d)      !volt*se/rad
          read(niterdb,10)(blank_d(j),j=1,nj_d)

c---rho grid, meters

          read(niterdb,'(a)')stflg
          read(niterdb,10)(r_d(j),j=1,nj_d)          !meters
c
       do j=0,50
         aj=j
         rho(j)=aj/50
       enddo
       rho(0)=1.D-6
       rho(50)=rho(50)-1.D-6
c
       do j=1,51
         aj=j
         x_d(j)=aj/51
       enddo
       do j=1,51
          write(*,51) j, x_d(j), r_d(j)
       enddo
c       return
c
        nex=jmaxm+1
        do j=1,nex
          aj=j-1
          rho_d(j)=aj/jmaxm
        enddo
        rho_d(1)=1.D-6
        rho_d(nex)=rho_d(nex)-1.D-6
c
c       call w_lin_interp_r8(51,rho,r_d,nex,rho_d,r_new,iflag,msg)
c       call inter_cspl(51,rho,r_d,nex,rho_d,r_new)
        call inter_cspl(51,x_d,r_d,301,rho_d,r_new)
c        write(*,*) 'back in test ...'
c
        do j=1,jmaxm+1
          write(*,51) j, rho(j-1), r_d(j), rho_d(j), r_new(j)
        enddo
c
 10   format(5(2x,1pe14.4))! common write/read  format
 51   format(i3,2x,0p1f6.4,1p6e12.4)
c
      end
c
      SUBROUTINE W_LIN_INTERP_R8(n1,x1,y1,n2,x2,y2,iflag,message)
!***********************************************************************
!W_LIN_INTERP does a linear interpolation to obtain n2 values for the
!  target array y2 on the grid x2 from the n1 values of y1 on the grid
!  x1
!References:
!  W.A.Houlberg 3/2000
!Input:
!  n1-number of abscissas and values in the source arrays
!  x1-array of source abscissas
!  y1-array of source values
!  n2-number of target values to be found
!  x2-array of target abscissas
!Output:
!  y2-array of target values
!  iflag-error and warning flag
!       =-1 warning
!       =0 no warnings or errors
!       =1 error
!  message-warning or error message (character)
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        n1,                      n2
      REAL*8         x1(*),                   y1(*),
     &               x2(*)
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag
      REAL*8         y2(*)
!Declaration of local variables
      INTEGER        i,                       il
c      do i=1,n1
c        write(*,*) i, x1(i), y1(i), ' x1,y1'
c      enddo
      IF(n1.lt.2) THEN
        iflag=1
        message='W_LIN_INTERP/ERROR:less than 2 points in source array'
      ELSE
        il=1
        DO i=1,n2
   10     IF(x2(i).lt.x1(1)) THEN
!           Use innermost data value
            y2(i)=y1(1)
            iflag=-1
            message='W_LIN_INTERP(1)/WARNING:x<x(1), use end point'
          ELSEIF(x2(i).eq.x1(1)) THEN
            y2(i)=y1(1)
          ELSEIF(x2(i).gt.x1(il+1)) THEN
            IF(il.lt.n1-1) THEN
!             Step x1 grid forward and loop
              il=il+1
              GOTO 10
            ELSE
!             Set to last value
              y2(i)=y1(n1)
              iflag=-1
              message='W_LIN_INTERP(2)/WARNING:x>x(n1), use end point'
            ENDIF
          ELSE
!           Interpolate
            y2(i)=y1(il)
     &            +(y1(il+1)-y1(il))*(x2(i)-x1(il))/(x1(il+1)-x1(il))
          ENDIF
c          write(*,*) i, x2(i), y2(i), ' x2,y2'
        ENDDO
      ENDIF
      RETURN
      END
