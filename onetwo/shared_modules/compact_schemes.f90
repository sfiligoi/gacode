

SUBROUTINE compact_deriv(profile,grid,np,order,zctr_deriv,msp,dpl,inpt_dpl,dpr,inpt_dpr,nrhs,dpdg)
!-------------------------------------------------------------------------------------
! --- subroutine uses compact methods to determine derivative over grid
! --- simultaneously for all grid points  out to the given boundary.
! ---
! --- INPUT:
! --- profile      The profile to be differentiated (defined as profile(1)..... profile(np))
! --- grid         the (assumed 1D)   grid points over which the profile is defined
! --- np           the number of grid points
! --- zctr_deriv   derivative type,if =0 then determine face centered(eg zone edge) derivatives that extend
! ---              from j = 1  to j = np. 
! ---              IF zctr_deriv = 1 then zone centered  derivatives from j = 1 to j = np-1
! ---              are determined. 
! --- order        order of compact scheme 4th and 6th order schemes are possible.
! ---              IF order is input as 2 standard  centered finite differencing is 
! ---              implemented. If order =  6 grid spacing must be uniform!!!!
! ---              If order is 2 or 4 then non uniform meshes are allowed.
! --- msp          mesh spacing, msp =0 for uniform (order = 2,6) or msp =1 for nonuniform
! ---              grid (order = 2,4).
! --- dpl          left end derivative (optional input)
! --- dpr          right end derivative (optional input)
! --  inpt_dpl     = 0 calculate left end derivative
! --               = 1 left end derivative is input in dpl
! --  inpt_dpr     = 0 calculate right end derivative
! --               = 1 right end derivative is input in dpr
! --- nrhs         >= 1 or more rhs. The lu factorization is done onece and then used
! ---                   repeatedly for each rhs(nrhs  > 1 not operational at this time).
! ---
! ---
! ---
! --- OUTPUT:
! --- dpdg(j)      IF zctr_deriv = 0 then
! ---              derivative at grid point j = 1,..np corresponding to zone edges
! ---              with dpdg(1) = dpl and dpdg(np) = dpr if inp_dpl and/or inpt_dpr are set 
! ---              If inp_dpl and/or inpt_dpr are not set then the compact scheme is used to
! ---              determine these values correct to at least O(grid_spacing**4). 
! ---
! ---              IF zctr_deriv = 1 then drivatives on  the zone centered grid are determined
! ---              (eq half way between the given input grid).
! ---              In this case dpdg(1) = derivative at location (grid(2)+grid(1))/2
! ---              dpdg(np-1) = derivative at location (grid(np)+grid(np-1))/2
! ---              and otherwise dpdg(j) = derivative at location (grid(j+1)+grid(j))/2
! ---              if  inp_dpl and/or inpt_dpr are set then dpdg(1) and/or dpdg(np-1) are
! ---              set to The given input values (dpr and/or dpl). Note that the output ranges from
! ---              1 to np-1 for the zone centered grid.
! ---
! -- NOTE:
! --               For uniform mesh you should use the second or sixth order scheme (order =2 or 6)
! --               For non uniform mesh use the second or forth order scheme (order =2 or 4)
! --               Zone centered schemes are obtained by high order interpolation from the
! --               zone edge values.
! --
! -- REF:          S. K. Lele J. Comp. Phys. 103,16-42(1992)
! --               Gamet, Int. J. Numer.Meth. Fluids 29,159-192(1999)
! --               End points boundary  conditions -- my notes.
!----------------------------------------------------------------------HSJ------------------------


      USE nrtype,                                          ONLY : Dp,I4B

      USE error_handler,                                   ONLY : iomaxerr, lerrno,terminate

      USE io_gcnmp,                                        ONLY : nlog
 
      USE     MPI_data,                                    ONLY : myid,master

      USE common_constants,                                ONLY : zeroc,izero

      IMPLICIT  NONE 
!      INTEGER(I4B),  INTENT(IN)                   ::  np,zctr_deriv,msp,inpt_dpl,inpt_dpr
!      INTEGER(I4B),  INTENT(INOUT)                ::  order
!      REAL(DP), DIMENSION(*),INTENT(IN)           ::  grid,profile
!      REAL(DP), DIMENSION(*),INTENT(INOUT)        ::  dpdg
!     REAL(DP), INTENT(IN)                        ::  dpl,dpr
      INTEGER(I4B)                  ::  np,zctr_deriv,msp,inpt_dpl,inpt_dpr
      INTEGER(I4B)                  ::  order
      REAL(DP), DIMENSION(*)        ::  grid,profile
      REAL(DP), DIMENSION(*)        ::  dpdg
      REAL(DP), INTENT(IN)          ::  dpl,dpr
      INTEGER(I4B)   j,ndim,info,nrhs

      REAL(DP), DIMENSION(np)                :: sol,d,rhs,fctr_deriv  ! temporary arrays, lower,
                                                                      ! main and upper diagonal
      REAL(DP), DIMENSION(np-1)              :: ld,ud
      REAL(DP)  alphai,betai,betanp,d1,dsq,d2,d3,d4,AI,dx,            &
                Bi,Ci,Di,Ei,him1,hi,hip1,hip2,h2,h3,hnp,hnpm1

 
      IF(.NOT. (zctr_deriv == 0 .OR. zctr_deriv ==1))order = -1        ! force error exit 
      IF( msp == 1 .AND. order .EQ. 6)order = -2                       ! force error exit


  deriv_order:     IF(order == 2)THEN
            
            IF(zctr_deriv     == 0)THEN                                   ! face center grid
               DO j = 2,np-1
                  dpdg(j) =(profile(j+1)-profile(j-1))/(grid(j+1)-grid(j-1))
               ENDDO
               IF(inpt_dpl .GT. izero )THEN
                  dpdg(1)   = dpl
               ELSE
                  dpdg(1)  = dpdg(2)                           !   "standard"  approximation
               ENDIF
               IF(inpt_dpr .GT. izero )THEN
                  dpdg(np)  = dpr
               ELSE
                  dpdg(np)  = dpdg(np-1)                      !   "standard"  approximation
               ENDIF
            ELSE                                              !   zctr_deriv = 1, zone center grid
               DO j= 1,np-1
                  dpdg(j) = (profile(j+1)-profile(j))/(grid(j+1)-grid(j))
               ENDDO
               IF(inpt_dpl .GT.  izero ) dpdg(1)   = dpl
               IF(inpt_dpr .GT.  izero ) dpdg(np)   = dpr

            ENDIF

  ELSE IF(order == 4)THEN   deriv_order


            meshspace:       IF(msp == 0)THEN                                      ! msp =0,  uniform mesh,order 4
               lerrno = iomaxerr+130                                    
               CALL terminate(lerrno,nlog)                        ! following codeing is not complete so dont use it

               dx           = grid(2)-grid(1) 
               ndim = np
               IF(zctr_deriv     == 0)THEN                              ! face centered grid
                  DO j= 2,np-1
                     rhs(j) = 3._DP*(profile(j+1)-profile(j-1))/dx
                     ld(j)  = 1._DP
                     d(j)   = 4._DP
                     ud(j)  = 1._DP
                  ENDDO
                  IF(inpt_dpl .GT.  izero)THEN
                     dpdg(1)   = dpl
                  ELSE
                     rhs(1)   = (2.5_DP*profile(1) + 2._DP*profile(2) + 0.5_DP*profile(3))/dx
                     d(1)     = 1._DP
                     ud(1)    = 2._DP
                     ld(1)    = 0.0_DP
                  ENDIF

                  IF(inpt_dpr .GT.  izero)THEN
                     dpdg(np)  = dpr
                  ELSE
                     rhs(np)    = (2.5_DP*profile(np) + 2._DP*profile(np-1) + 0.5_DP*profile(np-2))/dx
                     d(np)     = 1._DP
                     ud(np)    = 0.0
                     ld(np)    = 2.0_DP
                  ENDIF



               ELSE       !   zctr_deriv = 1                             ! zone centered grid


               ENDIF

            ELSE meshspace                    ! msp =1, non uniform mesh,order 4

                zone_center:    IF(zctr_deriv     == 0)THEN                               !  zone edge  grid
                     hi      = grid(2) - grid(1)
                     hip1    = grid(3) - grid(2)
                     hip2    = grid(4) - grid(3)
                     alphai  = 1._DP/3._DP ; betai = alphai
                     DO j    = 3,np-2
                        him1    = hi
                        hi      = hip1
                        hip1    = hip2
                        hip2    = grid(j+2)-grid(j+1)
                        d2      = hi+hip1
                        d3      = d2+ him1
                        d4      = d3 + hip2
                        Ai      =   (  him1*hi*hip1  +hi**2*hip1 + him1*hi*hip2 + hi**2*hip2 +                            &
                             (-him1*hi**2-him1*hi*hip1-him1*hi*hip2)*alphai +                                             &
                             ( -him1*hi*hip1 - hi**2*hip1-him1*hip1**2 -2._DP*hi*hip1**2 -hip1**3                         &
                             +him1*hi*hip2 +hi**2*hip2  +2._DP*him1*hip1*hip2 + 4._DP*hi*hip1*hip2 + 3._DP*hip1**2*hip2   &
                             )*betai)/(hip1*d2*d3*hip2)
                        !print *, 'ai,corct =',ai,  7._DP/(9._DP*(grid(2)-grid(1)))
                        Bi      =   ( -(him1+hi)*(hip1**2 +hip1*hip2) +                                                   &
                             (-3._DP*him1*hi**2 -4._DP*him1*hi*hip1                                                       &
                             +hi**3 + 2._DP*hi**2*hip1 -him1*hip1**2 + hi*hip1**2 -2._DP*him1*hi*hip2 - him1*hip1*hip2    &
                             +hi**2*hip2 + hi*hip1*hip2)*alphai                                                           &
                             +( him1*hip1*hip2 + hi*hip1*hip2 + hip1**2*hip2)*betai)/                                     &
                             (him1*hi*d2*(hi+hip1+hip2))
                        !print *, 'bi,corct =',bi,  -7._DP/(9._DP*(grid(2)-grid(1)))
                        Ci      = ( -hi*hip1*(him1+ hi)                                                                   &
                             +him1*hi*(hi + hip1)*alphai                                                                  &
                             +hip1*(hi*(him1 + hi  +2._DP*hip1)  + hip1*(him1 + hip1))*betai ) /                          &
                             (hip2*(hip1+hip2)*(hi+hip1+hip2)*d4)
                        !print *, 'ci,corct =',ci,  1._DP/(36._DP*(grid(2)-grid(1)))
                        Di      = (  hi*hip1*(hip1+hip2) -hi*alphai*(hi**2+2._DP*hi*hip1+hip1**2+hi*hip2+hip1*hip2)       &
                             -hip1*hip2*betai*(hi+hip1))/(him1*(him1+hi)*d3*d4)
                        !print *, 'di,corct =',di,  -1._DP/(36._DP*(grid(2)-grid(1)))
                        Ei      = -(Ai +Bi + Ci + Di )
                        rhs(j)  = Ai*profile(j+1) + Bi*profile(j-1)+Ci*profile(j+2) + Di*profile(j-2) + Ei*profile(j)
                        d(j)    = 1._DP
                        ud(j)   = betai
                        ld(j-1) = alphai 
                     ENDDO


                     ! left side  bc at i =1,2
                     h2       = grid(2)-grid(1)
                     h3       = grid(3)-grid(2)
                     d1       = h2 + h3      ;   dsq      = d1*d1
                     alphai   = d1/h3

                     IF(inpt_dpl .GT.  izero )THEN
                        rhs(1)   = dpl
                        d(1)     = 1.0_DP
                        ud(1)    = 0.0_DP
                     ELSE
                        Ai       = -(3._DP*h2+2.*h3)/(h2*d1)
                        Bi       = d1*(2._DP*h3-h2)/(h2*h3*h3)
                        Ci       = h2*h2/(h3*h3*d1)
                        rhs(1)   = Ai*profile(1) + Bi*profile(2) + Ci*profile(3)
                        d(1)     = 1._DP
                        ud(1)    = alphai
                     ENDIF

                     Ai       = -2._DP*h3*h3*(2._dp*h2+h3)/(h2*d1*dsq)
                     Bi       = 2._DP*(h3-h2)/(h2*h3)
                     Ci       = 2._DP*h2*h2*(h2+2._DP*h3)/(h3*d1*dsq)
                     alphai   = h3*h3/dsq
                     betai    = h2*h2/dsq
                     rhs(2)   = Ai*profile(1) + Bi*profile(2) + Ci*profile(3)
                     d(2)     = 1._DP
                     ud(2)    = betai
                     ld(1)    = alphai


                     ! right side bc at i = np,np-1:
                     hnp       = grid(np)-grid(np-1) ; hnpm1 = grid(np-1)-grid(np-2)
                     betanp    = hnp + hnpm1

                     d(np-1)   = 1._DP
                     ud(np-1)  = (hnpm1/betanp)**2
                     ld(np-2)  = (hnp/betanp)**2
                     Ai        = 2._DP*(hnpm1 + 2._DP*hnp)*hnpm1**2/(hnp*betanp**3)
                     Bi        = 2._DP*(hnp**2 - hnpm1**2)/(betanp*hnp*hnpm1)
                     Ci        = -2._DP*(hnp + 2._DP*hnpm1)*hnp**2/(hnpm1*betanp**3)
                     rhs(np-1) = Ai*profile(np)+ Bi*profile(np-1) + Ci*profile(np-2)
                     IF(inpt_dpr .GT.  izero )THEN                                    ! right end deriv input 
                        rhs(np)   = dpr
                        d(np)     = 1.0
                        ld(np-1)  = 0.0_DP
                     ELSE
                        Ai        = (3._DP*hnpm1+2._Dp*hnp)/(betanp*hnp)
                        Bi        = (hnp**2 -3._DP*hnpm1**2 -2._DP*hnp*hnpm1)/(hnpm1*hnp*betanp)
                        Ci        = -hnp/(betanp*hnpm1)
                        rhs(np)   = Ai*profile(np)+ Bi*profile(np-1) + Ci*profile(np-2)
                        d(np)     = 1.0_DP
                        ld(np-1)  = 2.0_DP
                     ENDIF


                     ndim = np
                     CALL dgtsv(ndim,nrhs,ld,d,ud,rhs,ndim,info)  ! solve tridiagonal system ld = lower dagonal
                                                                  !                          ud = upper
                                                                  !                           d = diagonal
                     dpdg(1:ndim) = rhs(1:ndim)
                     IF( info .NE. 0)THEN
                        IF(MYID == MASTER)THEN
                           PRINT *,'info in dtsv =',info
                           PRINT *,'rhs =',rhs
                           PRINT *,'ld =',ld
                           PRINT *,'d =',d
                           PRINT *,'ud =',ud
                        ENDIF
                        lerrno = iomaxerr+129                                    
                        CALL terminate(lerrno,nlog)
                     ENDIF

                ELSE   zone_center           
                     IF(MYID == MASTER)THEN
                        PRINT *,'GRID =',grid(1:np)
                     ENDIF
                     !   order=4,zctr_deriv = 1, zone centered grid
                     lerrno = iomaxerr+148                              !   need to develop codeing                  
                     CALL terminate(lerrno,nlog)
                ENDIF  zone_center


         ENDIF  meshspace


  ELSE IF(order == 6)THEN            deriv_order               
          mesh_spcing:    IF(msp == 0)THEN                            ! selects uniform grid 
            dx           = grid(2)-grid(1)
            deriv_type: IF(zctr_deriv == 0 .OR. zctr_deriv == 1)THEN  ! face centered uniform grid
               alphai = 1._Dp/3._DP ; betai = alphai
               DO j = 3,np-2
                  rhs(j)   = ((7._DP/9._DP)*(profile(j+1)-profile(j-1))   &
                                +(1._DP/36._DP)*(profile(j+2)-profile(j-2)))/dx
                  ld(j-1)  = alphai                          ! LD SYMMETRY BROKEN BY DGTSV,ugly
                  d(j)     = 1._DP
                  ud(j)    = betai
               ENDDO
               IF(inpt_dpl .GT.  izero )THEN                                    ! left end deriv input 
                  rhs(1)   = dpl
                  d(1)     = 1.0_DP
                  ud(1)    = 0.0_DP


                  rhs(2)   = 0.75_DP*(profile(3)-profile(1))/dx - 0.25_DP*dpl
                  ld(1)    = 0.0_DP                                 ! ld(1) starts in row 2 for dgtsv
                  d(2)     = 1.0_DP
                  ud(2)    = 0.25_DP
               ELSE                                                !eq for  i=1,2
                  rhs(1)   = (-2.5_DP*profile(1) +2._Dp*profile(2) + 0.5_DP*profile(3))/dx
                  d(1)     = 1._DP
                  ud(1)    = 2._DP

                  rhs(2)   = 0.75_DP*(profile(3)-profile(1))/dx
                  d(2)     = 1._DP
                  ud(2)    = 0.25_DP
                  ld(1)    = 0.25_DP

               ENDIF

               IF(inpt_dpr .GT.  izero)THEN                                    ! right end deriv input 
                  rhs(np)   = dpr
                  d(np)     = 1.0
                  ld(np-1)    = 0.0_DP

                  rhs(np-1) = -0.75_DP*(profile(np-2)-profile(np))/dx -0.25*dpr
                  ld(np-2)  = 0.25_DP
                  d(np-1)   = 1.0_DP
                  ud(np-1)  = 0.0_DP
               ELSE                                                ! get deriv eq at i = np-1,np-2
                  rhs(np-1) = -0.75_DP*(profile(np-2)- profile(np))/dx
                  d(np-1)   = 1.0_DP
                  ud(np-1)  = 0.25_DP
                  ld(np-2)  = 0.25_DP

                  rhs(np)   = (2.5_DP*profile(np)-2._DP*profile(np-1) -0.5_DP*profile(np-2))/dx
                  rhs(np)   =  (3.0_DP*profile(np)-4.0_DP*profile(np-1)+ profile(np-2))/(2._DP*dx)
                  d(np)     = 1.0_DP
                  ld(np-1)  = 2.0_DP  ; ld(np-1) =0.0_DP

               ENDIF
               ndim = np
               CALL dgtsv(ndim,nrhs,ld,d,ud,rhs,ndim,info)        ! solve tridiagonal system ld = lower dagonal
                                                                  !                          ud = upper
                                                                  !                           d = diagonal
               dpdg(1:ndim) = rhs(1:ndim)
               IF( info .NE. 0)THEN
                  lerrno = iomaxerr+145                   
                  CALL terminate(lerrno,nlog)
               ENDIF

            ELSE                           deriv_type              ! error zctr_deriv must = 0 or 1 
               lerrno = iomaxerr+132                               ! not programmed
               CALL terminate(lerrno,nlog)
            ENDIF                          deriv_type               ! zctr_deriv = 0,1 option

         ELSEIF(msp == 1)THEN              mesh_spcing             ! non uniform grid not programmed for 6 th order
            lerrno = iomaxerr+131                                
            CALL terminate(lerrno,nlog)
         ENDIF                             mesh_spcing             ! msp 0,1 option

  ELSE                                 deriv_order             ! other input errors
         ! get to here if order .ne. 6,4,2
         IF(order == -1 )THEN
            lerrno = iomaxerr+127                                    
            CALL terminate(lerrno,nlog)
         ENDIF

         IF(order == -2)THEN
            lerrno = iomaxerr+128                                    
            CALL terminate(lerrno,nlog)
         ENDIF

  ENDIF                                 deriv_order 

      IF (zctr_deriv == 0  .OR. msp == 1  .OR. order .NE. 6 )RETURN
!--------------------------------------------------------------------------------------
!---  IF zone centered derivatives are requested we get them by
!---  10th order interpolation from the zone edge grid values calculated above using order=6:
!---  The implemented interpolation scheme requires that the grid be uniform(msp = 0).
!--- ( I do not have the formulation for direct determination of zone center derivatives
!---   using compact schemes worked out in uniform or non uniform grids)
!--------------------------------------------------------------------HSJ----------------
      fctr_deriv(1:np) = dpdg(1:np)                   ! dpdg face centered derivs input
      dpdg(1:np)       = zeroc
      CALL midpt_interp(fctr_deriv,np,nrhs,dpdg)      ! dpdg, zone center derivs , output

      RETURN

END SUBROUTINE compact_deriv




      SUBROUTINE midpt_interp(fin,np,nrhs,fout)
!--------------------------------------------------------------------------------------
! --- High order midpoint interpolation on UNIFORM  grid. Solves pentadiagonal system for
! --- values half way between the input grid points (which are assumed uniformly spaced).
! --- Uses dgbsv banded matrix storage mode for M. Only a single rhs is allowed (nrhs =1)
! --- INPUT:
! ---    fin(1:np) values at grid points
! ---
! --- OUTPUT:
! ---   fout(1:np-1)  values at the half way points between the grid input points.
! ---
!-----------------------------------------------------------------HSJ------------------

      USE nrtype,                                          ONLY : Dp,I4B
      USE common_constants,                                ONLY : zeroc
      USE error_handler,                                   ONLY : iomaxerr, lerrno,terminate
      USE io_gcnmp,                                        ONLY : nlog


      IMPLICIT     NONE 
      REAL(DP),    INTENT(IN)         :: fin(np)
      REAL(DP),    INTENT(INOUT)      :: fout(np-1)
      INTEGER(I4B),INTENT(IN)     :: np
      INTEGER(I4B) i,npm1,nrhs,info
      INTEGER,     PARAMETER          :: ku = 2,kl = 2           ! 2 lower , 2 upper diagonals
      INTEGER,     PARAMETER          :: ldm = 2*kl + ku +1      ! extra kl for dgbsv work space
      REAL(DP)     beta,alpha,a,b,c
      REAL(DP),    DIMENSION(np-1)    :: rhs
      REAL(DP),    DIMENSION(ldm,np-1):: M
      INTEGER,     DIMENSION(np-1)    :: ipiv


      ! for 10'th order:
      alpha     = 10._DP/21._DP ; beta = 5._Dp/126._Dp ; a = 5._DP/6._DP
      b         = 5._DP/28._DP  ; c    = 1._DP/252._DP

      npm1      = np - 1  ; nrhs = 1
      M(:,:)    = zeroc
      !banded matrix mapping is (kl+ku+1+i-j,j) <==i,j

      DO i=3,np-3
         rhs(i)            = c*(fin(i+3) +fin(i-2))+ b*(fin(i+2)+fin(i-1)) + a*(fin(i+1)+fin(i))
         M(kl+ku+1,i)      = 1._DP  !i=j,   diag element 
         M(kl+ku+2,i-1)    = alpha  !i,j-1, first lower diagonal
         M(kl+ku+3,i-2)    = beta   !i,j-2, second lower diagonal
         M(kl+ku,i+1)      = alpha  !i,j+1  first upper diagonal
         M(kl+ku-1,i+2)    = beta   !i,j+2  second upper diagonal
      ENDDO

      !boundary elements:
         M(kl+ku+1,1)        = 1._DP       !(1,1)
         M(kl+ku,2)          = zeroc       !(1,2)       
         M(kl+ku+2,1)        = zeroc       !(2,1)       
         M(kl+ku+1,2)        = 1._DP       !(2,2)
         ! fourth order accurate (eg O(h^4)) solution for f(i=3/2):
         rhs(1)              = (5._DP*fin(1) + 15._DP*fin(2) -5._DP*fin(3) +fin(4))/16._DP
         ! fourth order accurate (eg O(h^4)) solution for f(i=5/2):
         rhs(2)              = (- fin(1) + 9._DP*fin(2) + 9._DP*fin(3) -fin(4))/16._DP




         M(kl+ku+1,np-1)     = 1._DP  !(np-1,np-1)
         M(kl+ku+1,np-2)     = 1._DP  !(np-2,np-2)
         ! fourth order accurate (eg O(h^4)) solution for f(i=np-3/2):
         rhs(np-2)           =  (-fin(np) +9._DP*fin(np-1) +9._DP*fin(np-2) &
                                           - fin(np-3))/16._DP
         ! fourth order accurate (eg O(h^4)) solution for f(i=np-1/2):
          rhs(np-1)          = (5._DP*fin(np) +  15._DP*fin(np-1) -5._DP*fin(np-2) &
                                + fin(np-3))/16._DP

         CALL dgbsv(npm1,kl,ku,nrhs,M,ldm,ipiv,rhs,npm1,info)     ! solve banded system ,kl =  # lower diagonals
                                                                  !                      ku =  # upper
                                                                  ! banded storage mode

         fout(1:npm1) = rhs(1:npm1)


         IF( info .NE. 0)THEN
             lerrno = iomaxerr+146                               
             CALL terminate(lerrno,nlog)
         ENDIF
         
         RETURN
      END SUBROUTINE midpt_interp

