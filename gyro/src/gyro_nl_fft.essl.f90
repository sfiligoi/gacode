!-------------------------------------------------------------------
! gyro_nl_fft.essl.f90
!
! PURPOSE:
!  This routine evaluates the ExB nonlinearity with either periodic 
!  or nonperiodic boundary conditions using the (F,G)-conservative 
!  difference scheme using FFT in the toroidal direction.
!
! NOTES:
!  ESSL-specific version.
!-------------------------------------------------------------------

subroutine gyro_nl_fft

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
  complex, dimension(0:n_max,n_x) :: fn_p
  complex, dimension(0:n_max,n_x) :: gn_p
  complex, dimension(0:n_max,n_x) :: fn_r
  complex, dimension(0:n_max,n_x) :: gn_r
  !
  complex, dimension(0:n_max,n_x) :: fgr
  complex, dimension(0:n_max,n_x) :: fgr_p
  complex, dimension(0:n_max,i1_buffer:i2_buffer) :: fgp
  complex, dimension(0:n_max,n_x) :: fgp_r
  complex, dimension(0:n_max,n_x) :: fg2
  !
  real, dimension(0:n_fft-1,n_x) :: fg_r
  real, dimension(0:n_fft-1,n_x) :: gf_r
  real, dimension(0:n_fft-1,n_x) :: gf_p
  real, dimension(0:n_fft-1,n_x) :: fg_p
  real, dimension(0:n_fft-1,n_x) :: f_pg_r
  real, dimension(0:n_fft-1,n_x) :: f_rg_p
  !
  complex, dimension(0:n_fft/2,6*n_x) :: x_fft
  real, dimension(0:n_fft-1,6*n_x) :: y_fft
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
              gn(nn,i) = h_tran(i,i_split,i_p(nn),is)
              fn(nn,i) = gyro_u_tran(i,i_split,i_p(nn),is)
           enddo ! nn
        enddo ! i

        !------------------------------------------------
        ! df/dp, dg/dp
        !
        do i=1,n_x
           do nn=0,n_max
              fn_p(nn,i) =  -i_c*n_p(nn)*fn(nn,i)
              gn_p(nn,i) =  -i_c*n_p(nn)*gn(nn,i)
           enddo
        enddo
        !------------------------------------------------

        !---------------------------------------------------------------
        ! df/dr, dg/dr
        !
        fn_r = (0.0,0.0)
        gn_r = (0.0,0.0)
        !
        do i_diff=-m_dx,m_dx-i_dx
           do i=1,n_x
              do nn=0,n_max
                 fn_r(nn,i) = fn_r(nn,i)+w_d1(i_diff)*fn(nn,i_loop(i+i_diff))
                 gn_r(nn,i) = gn_r(nn,i)+w_d1(i_diff)*gn(nn,i_Loop(i+i_diff))
              enddo ! nn
           enddo ! i
        enddo ! i_diff
        !---------------------------------------------------------------

        !---------------------------------------------------------------
        ! Supervector loading I
        ! 
        do i=1,n_x
           do nn=0,n_max
              x_fft(nn,i) = fn(nn,i)
              x_fft(nn,n_x+i) = gn(nn,i)
              x_fft(nn,2*n_x+i) = fn_p(nn,i)
              x_fft(nn,3*n_x+i) = gn_p(nn,i)
              x_fft(nn,4*n_x+i) = fn_r(nn,i)
              x_fft(nn,5*n_x+i) = gn_r(nn,i) 
           enddo
        enddo
        !
        ! ... dealias
        !
        do nn=n_max+1,n_fft/2
           x_fft(nn,:) = (0.0,0.0)
        enddo
        !---------------------------------------------------------------

        !---------------------------------------------------
        ! Backward FFT
        !
        CALL DCRFT(0, &
             x_fft, &
             n_fft/2+1, &
             y_fft, &
             n_fft, &
             n_fft, &
             6*n_x, &
             -1, &
             1.0, &
             aux1_dcrft, &
             22000, &
             aux2_dcrft,&
             20000)
        !---------------------------------------------------

        !---------------------------------------------------------------
        ! Real space multiplications
        !
        do i=1,n_x
           do nn=0,n_fft-1
              !f dg/dr
              fg_r(nn,i) =  y_fft(nn,i)*y_fft(nn,5*n_x+i)

              !g df/dr
              gf_r(nn,i) =  y_fft(nn,n_x+i)*y_fft(nn,4*n_x+i)

              !g df/dp
              gf_p(nn,i) =  y_fft(nn,n_x+i)*y_fft(nn,2*n_x+i)

              !f dg/dp
              fg_p(nn,i) =  y_fft(nn,i)*y_fft(nn,3*n_x+i)

              ! df/dp dg/dr
              f_pg_r(nn,i) =  y_fft(nn,2*n_x+i)*y_fft(nn,5*n_x+i)

              ! df/dr dg/dp
              f_rg_p(nn,i) =  y_fft(nn,4*n_x+i)*y_fft(nn,3*n_x+i)
           enddo !nn
        enddo ! i 

        !---------------------------------------------------------------
        ! Supervector loading II
        ! 
        do i=1,n_x
           y_fft(:,i) = fg_r(:,i)
           y_fft(:,n_x+i) = gf_r(:,i)
           y_fft(:,2*n_x+i) = gf_p(:,i)
           y_fft(:,3*n_x+i) = fg_p(:,i)
           y_fft(:,4*n_x+i) = f_pg_r(:,i)
           y_fft(:,5*n_x+i) = f_rg_p(:,i)    
        enddo
        !---------------------------------------------------------------

        !---------------------------------------------------------------
        ! Forward FFT
        !
        CALL DRCFT(0, &
             y_fft, &
             n_fft, &
             x_fft, &
             n_fft/2+1, &
             n_fft, &
             6*n_x, &
             1, &
             1.0, &
             aux1_drcft, &
             22000, &
             aux2_drcft,&
             20000)
        !---------------------------------------------------------------

        !----------------------------------------------------------------
        ! f dg/dr - g df/dr, 
        !  g df/dp - f dg/dp, 
        !   df/dp dg/dr - df/dr dg/dp
        !
        do i=1,n_x
           do nn=0,n_max

              !f dg/dr - g df/dr
              fgr(nn,i) = (x_fft(nn,i)-x_fft(nn,n_x+i))/n_fft

              !g df/dp - f dg/dp
              fgp(nn,i) = (x_fft(nn,2*n_x+i)-x_fft(nn,3*n_x+i))/n_fft

              !df/dp dg/dr - df/dr dg/dp
              fg2(nn,i) = (x_fft(nn,4*n_x+i)-x_fft(nn,5*n_x+i))/n_fft

           enddo
        enddo
        !----------------------------------------------------------------

        !---------------------------------------------------------------
        ! d/dr (g df/dp - f dg/dp)
        !
        fgp_r = (0.0,0.0)
        do i_diff=-m_dx,m_dx-i_dx
           do i=1,n_x
              do nn=0,n_max
                 fgp_r(nn,i) = fgp_r(nn,i)+w_d1(i_diff)*fgp(nn,i_loop(i+i_diff))
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
              nl(nn,i) = fgr_p(nn,i)+fgp_r(nn,i)+fg2(nn,i)
           enddo ! nn
        enddo ! i
        !------------------------------------------------

        !-----------------------------------------------------
        ! Finally, update global RHS (use h_tran for efficiency):
        !
        do nn=0,n_max
           h_tran(:,i_split,i_p(nn),is) = (c_nl_i(:)/3.0)*nl(nn,:)
        enddo ! nn
        !
        !-----------------------------------------------------   

     enddo ! i_split 
  enddo ! is

end subroutine gyro_nl_fft

