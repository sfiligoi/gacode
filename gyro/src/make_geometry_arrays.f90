!-----------------------------------------------------
! make_geometry_arrays.f90 [caller BigScience]
!
! PURPOSE:
!  Definition of arrays which require calls to 
!  general geometry library.
!-----------------------------------------------------

subroutine make_geometry_arrays

  use gyro_globals
  use GEO_interface

  !--------------------------------------
  implicit none
  !--------------------------------------


  do i=1,n_x

     call gyro_to_geo(i)

     do k=1,n_lambda
        do m=1,n_stack

           call GEO_interp(theta_t(i,k,m))

           b0_t(i,k,m)      = GEO_b
           g_theta_t(i,k,m) = GEO_g_theta
           grad_r_t(i,k,m)  = GEO_grad_r
           qrat_t(i,k,m)    = GEO_gq

           if (geo_gradbcurv_flag == 1) then

              ! No explicit pressure effect.  

              ! This is also known as the grad B -> curvature 
              ! rule, where we remark that curvature has no 
              ! pressure dependence (but grad B does!).

              cos_t(i,k,m)   = GEO_gcos1
              cos_p_t(i,k,m) = 0.0

           else 

              ! Retain full grad-B and curvature effects.

              cos_t(i,k,m)   = GEO_gcos1+GEO_gcos2
              cos_p_t(i,k,m) = -GEO_gcos2

           endif

           bt_t(i,k,m)   = GEO_bt
           bp_t(i,k,m)   = GEO_bp
           bigr_t(i,k,m) = GEO_bigr

           captheta_t(i,k,m) = GEO_captheta
           sin_t(i,k,m)      = GEO_gsin

           usin_t(i,k,m) = GEO_usin
           ucos_t(i,k,m) = GEO_ucos

        enddo ! m
     enddo ! k 

     do j=1,n_theta_plot

        call GEO_interp(theta_plot(j))

        b0_plot(i,j)      = GEO_b
        g_theta_plot(i,j) = GEO_g_theta

     enddo ! j

  enddo !i

  call write_geometry_arrays(trim(path)//'geometry_arrays.out',4)

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[make_geometry_arrays done]'
  endif

end subroutine make_geometry_arrays
