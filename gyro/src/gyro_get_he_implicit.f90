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
  complex :: temp
  !---------------------------


  if (n_field == 1) then

     ! ELECTROSTATIC

     p_nek_loc = 0
     do p_nek=1+i_proc_1,n_nek_1,n_proc_1
        p_nek_loc = p_nek_loc+1

!$omp parallel do default(shared) private(m,j,temp)
        do i=1,n_x
           do m=1,n_stack
              temp = 0.0
              do j=1,n_blend
                 temp = temp +(c_blend(j,m,i,p_nek_loc)-o_f(j,m,i,p_nek_loc))*&
                      field_blend(j,i,1)
              enddo ! j
              h(m,i,p_nek_loc,n_spec) = h(m,i,p_nek_loc,n_spec)+&
                   alpha_s(n_spec,i)*temp
           enddo ! m
        enddo ! i
!$omp end parallel do

     enddo ! p_nek_loc

  else if (n_field == 2) then

     ! ELECTROMAGNETIC -- A_parallel only

     p_nek_loc = 0
     do p_nek=1+i_proc_1,n_nek_1,n_proc_1
        p_nek_loc = p_nek_loc+1

!$omp parallel do default(shared) private(m,j,temp)
        do i=1,n_x           
           do m=1,n_stack
              temp = 0.0
              do j=1,n_blend
                 temp = temp+(c_blend(j,m,i,p_nek_loc)-o_f(j,m,i,p_nek_loc))*&
                      field_blend(j,i,1)-(c_blend(j,m,i,p_nek_loc)*&
                      v_para(m,i,p_nek_loc,n_spec)-o_fv(j,m,i,p_nek_loc))*&
                      field_blend(j,i,2)
              enddo ! j
              h(m,i,p_nek_loc,n_spec) = h(m,i,p_nek_loc,n_spec)+alpha_s(n_spec,i)*temp 
           enddo ! m
        enddo ! i
!$omp end parallel do

     enddo ! p_nek_loc

  else

     ! ELECTROMAGNETIC -- A_parallel and B_parallel


     p_nek_loc = 0
     do p_nek=1+i_proc_1,n_nek_1,n_proc_1
        p_nek_loc = p_nek_loc+1

!$omp parallel do default(shared) private(m,j,temp)
        do i=1,n_x
           do m=1,n_stack 
              temp = 0.0
              do j=1,n_blend
                 temp = temp+ &
                      (c_blend(j,m,i,p_nek_loc)-o_f(j,m,i,p_nek_loc))*&
                      field_blend(j,i,1)-&
                      (c_blend(j,m,i,p_nek_loc)*v_para(m,i,p_nek_loc,n_spec)- &
                      o_fv(j,m,i,p_nek_loc))*field_blend(j,i,2)+&
                      energy(nek_e(p_nek),indx_e)*lambda(i,nek_k(p_nek))*&
                      tem_s(n_spec,i)/z(n_spec)*&
                      (c_blend(j,m,i,p_nek_loc)-o_f(j,m,i,p_nek_loc))*field_blend(j,i,3)
              enddo ! j  
              h(m,i,p_nek_loc,n_spec) = h(m,i,p_nek_loc,n_spec)+alpha_s(n_spec,i)*temp 
           enddo ! m
        enddo ! i
!$omp end parallel do

     enddo ! p_nek_loc

  endif


  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_get_he_implicit done]'
  endif

end subroutine gyro_get_he_implicit
