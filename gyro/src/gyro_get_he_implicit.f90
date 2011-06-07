!-------------------------------------------------------
! gyro_get_he_implicit.f90
!
! PURPOSE:
!  After the ion distribution and fields are advanced,
!  we can compute the electron distribution.  
!--------------------------------------------------------

subroutine gyro_get_he_implicit

  use gyro_globals
  use gyro_pointers
  use math_constants

  !---------------------------
  implicit none
  !---------------------------

  p_nek_loc = 0

  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek)  
     k  = nek_k(p_nek)   

     ck = class(k)

     if (n_field == 1) then
        
        ! ELECTROSTATIC
        
        do i=1,n_x
           
           do m=1,n_stack
              
              m0 = m_phys(ck,m)
              
              do j=1,n_blend
                 
                 h(m,i,p_nek_loc,n_spec) = h(m,i,p_nek_loc,n_spec)+&
                      alpha_s(n_spec,i)*&
                      (c_blend(j,m0,i,p_nek_loc)-o_f(j,m,i,p_nek_loc))*&
                      field_blend(j,i,1)
                 
              enddo ! j
              
           enddo ! m
           
        enddo ! i
        
     else if (n_field == 2) then
        
        ! ELECTROMAGNETIC -- A_parallel only
        
        do i=1,n_x
           
           do m=1,n_stack
              
              m0 = m_phys(ck,m)
              
              do j=1,n_blend
                 
                 h(m,i,p_nek_loc,n_spec) = h(m,i,p_nek_loc,n_spec) &
                      +alpha_s(n_spec,i)*&
                      (c_blend(j,m0,i,p_nek_loc)-o_f(j,m,i,p_nek_loc))*field_blend(j,i,1) &
                      -alpha_s(n_spec,i)*&
                      (c_blend(j,m0,i,p_nek_loc)*v_para(m,i,p_nek_loc,n_spec)- &
                      o_fv(j,m,i,p_nek_loc))*field_blend(j,i,2)
                 
              enddo ! j
              
           enddo ! m
           
        enddo ! i
        
     else
        
        ! ELECTROMAGNETIC -- A_parallel and B_parallel
        
        do i=1,n_x
           
           do m=1,n_stack
              
              m0 = m_phys(ck,m)
              
              do j=1,n_blend
                 
                 h(m,i,p_nek_loc,n_spec) = h(m,i,p_nek_loc,n_spec) &
                      +alpha_s(n_spec,i)*&
                      (c_blend(j,m0,i,p_nek_loc)-o_f(j,m,i,p_nek_loc))*field_blend(j,i,1) &
                      -alpha_s(n_spec,i)*&
                      (c_blend(j,m0,i,p_nek_loc)*v_para(m,i,p_nek_loc,n_spec)- &
                      o_fv(j,m,i,p_nek_loc))*field_blend(j,i,2) &
                      +alpha_s(n_spec,i)*&
                      energy(ie,indx_e)*lambda(i,k)*tem_s(n_spec,i)/z(n_spec)*&
                      (c_blend(j,m0,i,p_nek_loc)-o_f(j,m,i,p_nek_loc))*field_blend(j,i,3)
                 
              enddo ! j
              
           enddo ! m
           
        enddo ! i
        
     endif

  enddo ! p_nek_loc

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_get_he_implicit done]'
  endif

end subroutine gyro_get_he_implicit
