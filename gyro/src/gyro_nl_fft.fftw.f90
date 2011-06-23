!----------------------------------------------------------------
! gyro_nl_fft.fftw.f90
!
! PURPOSE:
!  This routine evaluates the ExB nonlinearity with periodic or
!  nonperiodic boundary conditions using the (F,G)-conservative 
!  difference scheme with FFT in the toroidal direction.
!
! NOTES:
!  FFTW_specific version.
!----------------------------------------------------------------

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
  complex, dimension(0:n_max,i1_buffer:i2_buffer) :: fgp
  !
  complex :: fn_p
  complex :: gn_p
  complex :: fn_r
  complex :: gn_r
  !
  complex :: fgr
  complex :: fgr_p
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
  complex :: fg_r_c
  complex :: gf_r_c
  complex :: f_pg_r_c
  complex :: f_rg_p_c
  !
  real, dimension(0:n_fft-1,6*n_x) :: v_fft
  real, dimension(0:n_fft-1,6*n_x) :: vt_fft
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

        !------------------------------------------------
        !
        v_fft = 0.0
        !
        do i=1,n_x
           do nn=0,n_max
              fn_r = (0.0,0.0)
              gn_r = (0.0,0.0)
              !---------------------------------------------------------------
              ! df/dr, dg/dr
              !
              do i_diff=-m_dx,m_dx-i_dx
                 fn_r = fn_r+w_d1(i_diff)*fn(nn,i_loop(i+i_diff))
                 gn_r = gn_r+w_d1(i_diff)*gn(nn,i_loop(i+i_diff))
              enddo ! i_diff
              !------------------------------------------------
              ! df/dp, dg/dp
              !
              fn_p = -i_c*n_p(nn)*fn(nn,i)
              gn_p = -i_c*n_p(nn)*gn(nn,i)
              !---------------------------------------------------------------
              ! Dealiasing and wrap-around of 
              ! f, g, df/dp, dg/dp, df/dr, dg/dr
              !--------------------------------------------------
              ! Supervector loading I
              !
              v_fft(nn,i) = real(fn(nn,i))
              v_fft(nn,n_x+i) = real(gn(nn,i))
              v_fft(nn,2*n_x+i) = real(fn_p)
              v_fft(nn,3*n_x+i) = real(gn_p)
              v_fft(nn,4*n_x+i) = real(fn_r)
              v_fft(nn,5*n_x+i) = real(gn_r)
              if (nn /= 0) then
                 v_fft(n_fft-nn,i) = aimag(fn(nn,i))
                 v_fft(n_fft-nn,n_x+i) = aimag(gn(nn,i))
                 v_fft(n_fft-nn,2*n_x+i) = aimag(fn_p)
                 v_fft(n_fft-nn,3*n_x+i) = aimag(gn_p)
                 v_fft(n_fft-nn,4*n_x+i) = aimag(fn_r)
                 v_fft(n_fft-nn,5*n_x+i) = aimag(gn_r)
              endif

           enddo ! nn
        enddo ! i

        !---------------------------------------------------
        ! Backward FFT
        !
        call rfftw_f77(plan_b,&
             6*n_x,&
             v_fft,&
             1,&
             n_fft,&
             vt_fft,&
             1,&
             n_fft)

        !---------------------------------------------------------------
        ! Real space multiplications
        !
        do i=1,n_x
           do nn=0,n_fft-1
              !f dg/dr
              fg_r =  vt_fft(nn,i)*vt_fft(nn,5*n_x+i)

              !g df/dr
              gf_r =  vt_fft(nn,n_x+i)*vt_fft(nn,4*n_x+i)

              !g df/dp
              gf_p =  vt_fft(nn,n_x+i)*vt_fft(nn,2*n_x+i)

              !f dg/dp
              fg_p =  vt_fft(nn,i)*vt_fft(nn,3*n_x+i)

              ! df/dp dg/dr
              f_pg_r =  vt_fft(nn,2*n_x+i)*vt_fft(nn,5*n_x+i)

              ! df/dr dg/dp
              f_rg_p =  vt_fft(nn,4*n_x+i)*vt_fft(nn,3*n_x+i)

              !---------------------------------------------------------------
              ! Supervector loading II
              ! 
              vt_fft(nn,i) = fg_r
              vt_fft(nn,n_x+i) = gf_r
              vt_fft(nn,2*n_x+i) = gf_p
              vt_fft(nn,3*n_x+i) = fg_p
              vt_fft(nn,4*n_x+i) = f_pg_r
              vt_fft(nn,5*n_x+i) = f_rg_p   
           enddo !nn
        enddo
        !---------------------------------------------------------------

        !---------------------------------------------------------------
        ! Forward FFT
        !
        call rfftw_f77(plan_f,&
             6*n_x,&
             vt_fft,&
             1,&
             n_fft,&
             v_fft,&
             1,&
             n_fft)

        !--------------------------------------------------
        ! Re-Construction of complex arrays
        ! 
        !--------------------------------------------------
        !  g df/dp - f dg/dp, 
        do i=1,n_x
           fgp(0,i) = (cmplx(v_fft(0,2*n_x+i),0.0) - &
                cmplx(v_fft(0,3*n_x+i),0.0))/n_fft
           do nn=1,n_max
              fgp(nn,i) = (cmplx(v_fft(nn,2*n_x+i),v_fft(n_fft-nn,2*n_x+i))- &
                   cmplx(v_fft(nn,3*n_x+i),v_fft(n_fft-nn,3*n_x+i)))/n_fft
           enddo ! nn 
        enddo ! i
        !--------------------------------------------------

        !---------------------------------------------------------------
        do i=1,n_x
           do nn=0,n_max
              !---------------------------------------------------------------
              ! d/dr (g df/dp - f dg/dp)
              !
              fgp_r = (0.0,0.0)
              do i_diff=-m_dx,m_dx-i_dx
                 fgp_r = fgp_r+w_d1(i_diff)*fgp(nn,i_loop(i+i_diff))
              enddo ! i_diff

              if (nn == 0) then
                 fg_r_c = cmplx(v_fft(0,i),0.0)
                 gf_r_c = cmplx(v_fft(0,n_x+i),0.0)
                 f_pg_r_c = cmplx(v_fft(0,4*n_x+i),0.0)
                 f_rg_p_c = cmplx(v_fft(0,5*n_x+i),0.0)   
              else
                 fg_r_c = cmplx(v_fft(nn,i),v_fft(n_fft-nn,i))
                 gf_r_c = cmplx(v_fft(nn,n_x+i),v_fft(n_fft-nn,n_x+i))
                 f_pg_r_c = cmplx(v_fft(nn,4*n_x+i),v_fft(n_fft-nn,4*n_x+i))
                 f_rg_p_c = cmplx(v_fft(nn,5*n_x+i),v_fft(n_fft-nn,5*n_x+i)) 
              end if

              !----------------------------------------------------------------
              ! f dg/dr - g df/dr, 
              !  g df/dp - f dg/dp,   (above, before this loop)
              !   df/dp dg/dr - df/dr dg/dp
              !

              !f dg/dr - g df/dr
              fgr = (fg_r_c-gf_r_c)/n_fft

              !df/dp dg/dr - df/dr dg/dp
              fg2 = (f_pg_r_c-f_rg_p_c)/n_fft

              ! d/dp (f dg/dr - g df/dr)
              fgr_p = -i_c*n_p(nn)*fgr

              !------------------------------------------------
              ! Arakawa scheme:
              !
              ! d/dp (f dg/dr - g df/dr)
              !   + d/dr (g df/dp - f dg/dp)
              !       + df/dp dg/dr - df/dr dg/dp
              !
              nl(nn,i) = fgr_p+fgp_r+fg2
              !------------------------------------------------
           enddo
        enddo

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

end subroutine gyro_nl_fft

