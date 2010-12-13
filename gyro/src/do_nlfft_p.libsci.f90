!-----------------------------------------------------------
! do_nlfft_p.libsci.f90 [caller: gyro_rhs_total]
!
! PURPOSE:
!  This routine evaluates the ExB nonlinearity with periodic 
!  boundary conditions using the (F,G)-conservative 
!  difference scheme with FFT in the toroidal direction.
!
! NOTES:
!  LIBSCI-specific version.
!-----------------------------------------------------------

subroutine do_nlfft_p

  use gyro_globals
  use gyro_pointers
  use gyro_nl_private
  use math_constants

  !--------------------------------------------
  implicit none
  !
  complex, dimension(-n_max:n_max,n_x) :: fn
  complex, dimension(-n_max:n_max,n_x) :: gn
  complex, dimension(0:n_max,n_x) :: nl
  !
  complex :: fn_p
  complex ::  gn_p
  complex ::  fn_r
  complex ::  gn_r
  !
  complex ::  fgr
  complex  :: fgr_p
  complex, dimension(0:n_max,n_x) :: fgp
  complex :: fgp_r
  complex :: fg2
  !
  real :: fg_r
  real :: gf_r
  real :: gf_p
  real :: fg_p
  real :: f_pg_r
  real :: f_rg_p
  !
  complex, dimension(0:nx_fft-1,6*n_x) :: x_fft
  real,    dimension(0:ny_fft-1,6*n_x) :: y_fft
  !--------------------------------------------

  do is=1,n_kinetic
     do i_split=1,msplit_SSUB

!dir$ prefervector
        do nn=0,n_max
           do i=1,n_x
              gn(nn,i) = h_M(i,i_split,i_p(nn),is)
              fn(nn,i) = gyro_u_M(i,i_split,i_p(nn),is)
           enddo
        enddo
        !------------------------------------------------

        x_fft(:,:) = (0.0,0.0)

        !---------------------------------------------------------------
        ! df/dr, dg/dr
        !
        !
!dir$ concurrent,prefervector
        do nn=0,n_max
!dir$ concurrent,prefervector
           do i=1,n_x
              fn_r = (0.0,0.0)
              gn_r = (0.0,0.0)
              do i_diff=-m_dx,m_dx-i_dx
                 fn_r = fn_r+w_d1(i_diff)*fn(nn,i_cyc(i+i_diff))
                 gn_r = gn_r+w_d1(i_diff)*gn(nn,i_cyc(i+i_diff))
              enddo ! i_diff
              !------------------------------------------------
              ! df/dp, dg/dp
              !
              fn_p = -i_c*n_p(nn)*fn(nn,i)
              gn_p = -i_c*n_p(nn)*gn(nn,i)
              !---------------------------------------------------------------
              ! Supervector loading I
              ! 
              x_fft(nn,i) = fn(nn,i)
              x_fft(nn,n_x+i) = gn(nn,i)
              x_fft(nn,2*n_x+i) = fn_p
              x_fft(nn,3*n_x+i) = gn_p
              x_fft(nn,4*n_x+i) = fn_r
              x_fft(nn,5*n_x+i) = gn_r 
           enddo
        enddo
        !---------------------------------------------------------------

        !---------------------------------------------------
        ! Backward FFT
        !
        CALL ZDFFTM(-1, &
             n_fft, &
             6*n_x, &
             1.0, &
             x_fft, &
             nx_fft, &
             y_fft, &
             ny_fft, &
             table_cs, &
             work_cs, &
             0)
        !---------------------------------------------------

        if (ny_fft > n_fft) then
           y_fft(n_fft:(ny_fft-1),:) = 0.0
        endif

        !---------------------------------------------------------------
        ! Real space multiplications
        !
!dir$ concurrent,preferstream
        do i=1,n_x
!dir$ concurrent,prefervector
           do nn=0,n_fft-1
              !f dg/dr
              fg_r =  y_fft(nn,i)*y_fft(nn,5*n_x+i)

              !g df/dr
              gf_r =  y_fft(nn,n_x+i)*y_fft(nn,4*n_x+i)

              !g df/dp
              gf_p =  y_fft(nn,n_x+i)*y_fft(nn,2*n_x+i)

              !f dg/dp
              fg_p =  y_fft(nn,i)*y_fft(nn,3*n_x+i)

              ! df/dp dg/dr
              f_pg_r =  y_fft(nn,2*n_x+i)*y_fft(nn,5*n_x+i)

              ! df/dr dg/dp
              f_rg_p =  y_fft(nn,4*n_x+i)*y_fft(nn,3*n_x+i)

              !---------------------------------------------------------------
              ! Supervector loading II
              ! 
              y_fft(nn,i) = fg_r
              y_fft(nn,n_x+i) = gf_r
              y_fft(nn,2*n_x+i) = gf_p
              y_fft(nn,3*n_x+i) = fg_p
              y_fft(nn,4*n_x+i) = f_pg_r
              y_fft(nn,5*n_x+i) = f_rg_p    
           enddo
        enddo
        !---------------------------------------------------------------

        !---------------------------------------------------------------
        ! Forward FFT
        !
        CALL DZFFTM(1, &
             n_fft, &
             6*n_x, &
             1.0, &
             y_fft, &
             ny_fft, &
             x_fft, &
             nx_fft, &
             table_sc, &
             work_sc, &
             0)
        !---------------------------------------------------------------

        !----------------------------------------------------------------
        ! f dg/dr - g df/dr, 
        !  g df/dp - f dg/dp, 
        !   df/dp dg/dr - df/dr dg/dp
        !
        do i=1,n_x
           do nn=0,n_max
              !g df/dp - f dg/dp
              fgp(nn,i) = (x_fft(nn,2*n_x+i)-x_fft(nn,3*n_x+i))/n_fft
           enddo
        enddo
        !----------------------------------------------------------------

        !---------------------------------------------------------------
        ! d/dr (g df/dp - f dg/dp)
        !
!dir$ prefervector
        do nn=0,n_max
           do i=1,n_x
              fgp_r = (0.0,0.0)
              do i_diff=-m_dx,m_dx-i_dx
                 fgp_r = fgp_r+w_d1(i_diff)*fgp(nn,i_cyc(i+i_diff))
              enddo ! i_diff

              !f dg/dr - g df/dr
              fgr = (x_fft(nn,i)-x_fft(nn,n_x+i))/n_fft
              ! d/dp (f dg/dr - g df/dr)
              fgr_p = -i_c*n_p(nn)*fgr
              !df/dp dg/dr - df/dr dg/dp
              fg2 = (x_fft(nn,4*n_x+i)-x_fft(nn,5*n_x+i))/n_fft

              !------------------------------------------------
              ! Arakawa scheme:
              !
              ! d/dp (f dg/dr - g df/dr)
              !   + d/dr (g df/dp - f dg/dp)
              !       + df/dp dg/dr - df/dr dg/dp
              !
              nl(nn,i) = fgr_p+fgp_r+fg2
           enddo ! nn
        enddo ! i
        !------------------------------------------------

        !-----------------------------------------------------
        ! Finally, update global RHS (use h_M for efficiency):
        !
        do nn=0,n_max
           h_M(:,i_split,i_p(nn),is) = (c_nl_i(:)/3.0)*nl(nn,:)
        enddo ! nn
        !
        !-----------------------------------------------------   

     enddo ! i_split 
  enddo ! is

end subroutine do_nlfft_p
