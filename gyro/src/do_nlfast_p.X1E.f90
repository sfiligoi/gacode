!-----------------------------------------------------------
! do_nlfast_p.X1E.f90 [caller: gyro_rhs_total]
!
! PURPOSE:
!  This routine evaluates the ExB nonlinearity with periodic 
!  boundary conditions using the (F,G)-conservative 
!  difference scheme.
!
! NOTES:
!  This version is optimized for the Cray X1E.
!-----------------------------------------------------------

subroutine do_nlfast_p

  use gyro_globals
  use gyro_pointers
  use gyro_nl_private
  use math_constants

  !--------------------------------------------
  implicit none
  !
  complex, dimension(n_x,-n_max:n_max) :: fn
  complex, dimension(n_x,-n_max:n_max) :: gn
  !
  complex, dimension(n_x,-n_max:n_max) :: fn_p
  complex, dimension(n_x,-n_max:n_max) :: gn_p
  complex, dimension(n_x,-n_max:n_max) :: fn_r
  complex, dimension(n_x,-n_max:n_max) :: gn_r
  !
  complex, dimension(n_x,0:n_max) :: fgr
  complex, dimension(n_x,0:n_max) :: fgr_p
  complex, dimension(n_x,0:n_max) :: fgp
  complex, dimension(n_x,0:n_max) :: fgp_r
  complex, dimension(n_x,0:n_max) :: fg2
  complex, dimension(n_x,0:n_max) :: nl
  !
  complex :: add1
  complex :: add2
  complex :: add3
  !--------------------------------------------


  do is=1,n_kinetic
     do i_split=1,msplit_SSUB

        do i=1,n_x
           do nn=0,n_max
              gn(i,nn) = h_M(i,i_split,i_p(nn),is)
              fn(i,nn) = gyro_u_M(i,i_split,i_p(nn),is)
           enddo ! nn
        enddo ! i

        do i=1,n_x
           do nn=1,n_max
              gn(i,-nn) = conjg(gn(i,nn))
              fn(i,-nn) = conjg(fn(i,nn))
           enddo ! nn
        enddo ! i

        !------------------------------------------------
        ! df/dp, dg/dp
        !
        do i=1,n_x
           do nn=-n_max,n_max
              fn_p(i,nn) = -i_c*n_p(nn)*fn(i,nn)
              gn_p(i,nn) = -i_c*n_p(nn)*gn(i,nn)
           enddo
        enddo
        !------------------------------------------------

        !---------------------------------------------------------------
        ! df/dr, dg/dr
        !
        fn_r = (0.0,0.0)
        gn_r = (0.0,0.0)
        !
        ! THIS IS EXPENSIVE LOOP #1
        !
        do i_diff=-m_dx,m_dx-i_dx
           do i=1,n_x
              ip = i+i_diff
              do nn=0,n_max
                 fn_r(i,nn) = fn_r(i,nn)+w_d1(i_diff)*fn(i_cyc(ip),nn)
                 gn_r(i,nn) = gn_r(i,nn)+w_d1(i_diff)*gn(i_cyc(ip),nn)
              enddo ! nn
           enddo ! i
        enddo ! i_diff

        do i=1,n_x
           do nn=1,n_max
              fn_r(i,-nn) = conjg(fn_r(i,nn))
              gn_r(i,-nn) = conjg(gn_r(i,nn))
           enddo ! nn
        enddo ! i
        !---------------------------------------------------------------

        !---------------------------------------------------------------
        ! Nonlinear convolution
        !
        do nn=0,n_max
           do i=1,n_x
              add1 = (0.0,0.0)
              add2 = (0.0,0.0)
              add3 = (0.0,0.0)
              do n1=-n_max+nn,n_max

                 ! f dg/dr - g df/dr

                 add1 = add1+&
                      fn(i,n1)*gn_r(i,nn-n1)-gn(i,n1)*fn_r(i,nn-n1)

                 ! g df/dp - f dg/dp

                 add2 = add2+&
                      gn(i,n1)*fn_p(i,nn-n1)-fn(i,n1)*gn_p(i,nn-n1)

                 ! df/dp dg/dr - df/dr dg/dp

                 add3 = add3+&
                      fn_p(i,n1)*gn_r(i,nn-n1)-fn_r(i,n1)*gn_p(i,nn-n1)

              enddo ! n1
              fgr(i,nn) = add1
              fgp(i,nn) = add2
              fg2(i,nn) = add3
           enddo ! i 
        enddo ! nn

        !---------------------------------------------------------------
        ! d/dr (g df/dp - f dg/dp)
        !
        ! THIS IS EXPENSIVE LOOP #2
        !
        fgp_r = (0.0,0.0)
        !
        do i_diff=-m_dx,m_dx-i_dx
           do i=1,n_x
              ip = i+i_diff
              do nn=0,n_max
                 fgp_r(i,nn) = fgp_r(i,nn)+w_d1(i_diff)*fgp(i_cyc(ip),nn)
              enddo ! nn 
           enddo ! i
        enddo ! i_diff
        !---------------------------------------------------------------


        !------------------------------------------------
        ! d/dp (f dg/dr - g df/dr)
        !
        do i=1,n_x
           do nn=0,n_max
              fgr_p(i,nn) = -i_c*n_p(nn)*fgr(i,nn)
           enddo ! nn
        enddo ! i
        !------------------------------------------------

        !------------------------------------------------
        ! Arakawa scheme:
        !
        ! d/dp (f dg/dr - g df/dr)
        !   + d/dr (g df/dp - f dg/dp)
        !       + df/dp dg/dr - df/dr dg/dp
        !
        do i=1,n_x
           do nn=0,n_max
              nl(i,nn) = fgr_p(i,nn)+fgp_r(i,nn)+fg2(i,nn)
           enddo ! nn
        enddo ! i
        !------------------------------------------------

        !-----------------------------------------------------
        ! Finally, update global RHS (use h_M for efficiency):
        !
        do i=1,n_x
           do nn=0,n_max
              h_M(i,i_split,i_p(nn),is) = (c_nl_i(i)/3.0)*nl(i,nn)
           enddo ! nn
        enddo ! i
        !
        !-----------------------------------------------------   

     enddo ! i_split
  enddo ! is

end subroutine do_nlfast_p
