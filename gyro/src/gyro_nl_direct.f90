!-------------------------------------------------------------------
! gyro_nl_direct.f90 
!
! PURPOSE:
!  This routine evaluates the ExB nonlinearity with both periodic
!  and nonperiodic boundary conditions using the (F,G)-conservative 
!  difference scheme.
!-------------------------------------------------------------------

subroutine gyro_nl_direct

  use gyro_globals
  use gyro_pointers
  use gyro_nl_private
  use math_constants

  !--------------------------------------------
  implicit none
  !
  complex, dimension(-n_max:n_max,i1_buffer:i2_buffer) :: fn
  complex, dimension(-n_max:n_max,i1_buffer:i2_buffer) :: gn
  complex, dimension(0:n_max,n_x) :: nl
  !
  complex, dimension(-n_max:n_max,n_x) :: fn_p
  complex, dimension(-n_max:n_max,n_x) :: gn_p
  complex, dimension(-n_max:n_max,n_x) :: fn_r
  complex, dimension(-n_max:n_max,n_x) :: gn_r
  !
  complex, dimension(0:n_max,n_x) :: fgr
  complex, dimension(0:n_max,n_x) :: fgr_p
  complex, dimension(0:n_max,i1_buffer:i2_buffer) :: fgp
  complex, dimension(0:n_max,n_x) :: fgp_r
  complex, dimension(0:n_max,n_x) :: fg2
  !
  complex :: add1
  complex :: add2
  complex :: add3
  !--------------------------------------------


  do is=1,n_kinetic
     do i_split=1,msplit_SSUB

        ! fn, gn and fgp have ip indices; 
        ! must be zeroed.

        do i=i1_buffer,0
           fn(:,i)  = (0.0,0.0)
           gn(:,i)  = (0.0,0.0)
           fgp(:,i) = (0.0,0.0)
        enddo
        do i=n_x+1,i2_buffer
           fn(:,i)  = (0.0,0.0)
           gn(:,i)  = (0.0,0.0)
           fgp(:,i) = (0.0,0.0)
        enddo

        do i=1,n_x
           do nn=0,n_max
              gn(nn,i) = h_M(i,i_split,i_p(nn),is)
              fn(nn,i) = gyro_u_M(i,i_split,i_p(nn),is)
           enddo ! nn
        enddo ! i

        do i=1,n_x
           do nn=1,n_max
              gn(-nn,i) = conjg(gn(nn,i))
              fn(-nn,i) = conjg(fn(nn,i))
           enddo ! nn
        enddo ! i

        !------------------------------------------------
        ! df/dp, dg/dp
        !
        do i=1,n_x
           do nn=-n_max,n_max
              fn_p(nn,i) = -i_c*n_p(nn)*fn(nn,i)
              gn_p(nn,i) = -i_c*n_p(nn)*gn(nn,i)
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
                 fn_r(nn,i) = fn_r(nn,i)+w_d1(i_diff)*fn(nn,i_loop(ip))
                 gn_r(nn,i) = gn_r(nn,i)+w_d1(i_diff)*gn(nn,i_loop(ip))
              enddo ! nn
           enddo ! i
        enddo ! i_diff

        do i=1,n_x
           do nn=1,n_max
              fn_r(-nn,i) = conjg(fn_r(nn,i))
              gn_r(-nn,i) = conjg(gn_r(nn,i))
           enddo ! nn
        enddo ! i
        !---------------------------------------------------------------

        !---------------------------------------------------------------
        ! Nonlinear convolution
        !
        ! ** (i,nn) order here fixes cache problem on 
        !    Opteron and Power4
        ! 
        do i=1,n_x
           do nn=0,n_max
              add1 = (0.0,0.0)
              add2 = (0.0,0.0)
              add3 = (0.0,0.0)
              do n1=-n_max+nn,n_max

                 ! f dg/dr - g df/dr

                 add1 = add1+&
                      fn(n1,i)*gn_r(nn-n1,i)-gn(n1,i)*fn_r(nn-n1,i)

                 ! g df/dp - f dg/dp

                 add2 = add2+&
                      gn(n1,i)*fn_p(nn-n1,i)-fn(n1,i)*gn_p(nn-n1,i)

                 ! df/dp dg/dr - df/dr dg/dp

                 add3 = add3+&
                      fn_p(n1,i)*gn_r(nn-n1,i)-fn_r(n1,i)*gn_p(nn-n1,i)

              enddo ! n1
              fgr(nn,i) = add1
              fgp(nn,i) = add2
              fg2(nn,i) = add3
           enddo ! nn
        enddo ! i 
        !---------------------------------------------------------------

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
                 fgp_r(nn,i) = fgp_r(nn,i)+w_d1(i_diff)*fgp(nn,i_loop(ip))
              enddo ! nn 
           enddo ! i_diff
        enddo ! i
        !---------------------------------------------------------------

        !------------------------------------------------
        ! d/dp (f dg/dr - g df/dr)
        !
        do i=1,n_x
           do nn=0,n_max
              fgr_p(nn,i) = -i_c*n_p(nn)*fgr(nn,i)
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
              nl(nn,i) = c_nl_i(i)*(fgr_p(nn,i)+fgp_r(nn,i)+fg2(nn,i))
           enddo ! nn
        enddo ! i
        !------------------------------------------------

        !-----------------------------------------------------
        ! Finally, update global RHS (use h_M for efficiency):
        !
        do nn=0,n_max
           h_M(:,i_split,i_p(nn),is) = nl(nn,:)/3.0
        enddo ! nn
        !
        !-----------------------------------------------------   

     enddo ! i_split
  enddo ! is

end subroutine gyro_nl_direct
