!-----------------------------------------------------
! get_delta_he.f90
!
! PURPOSE:
!  Compute the "explicit" part of h_e.
!-----------------------------------------------------

subroutine gyro_get_delta_he

  use gyro_globals
  use gyro_pointers

  !------------------------------------------
  implicit none
  !
  integer :: n_s2
  integer :: mp
  !
  complex :: h_temp(n_stack,n_x,n_nek_loc_1)
  !------------------------------------------

  n_s2 = n_stack/2

  h_temp(:,:,:)   = h(:,:,:,n_spec)
  h(:,:,:,n_spec) = (0.0,0.0)

  p_nek_loc = 0

  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek)  
     k  = nek_k(p_nek)   

     ck = class(k)

     if (ck == 1) then

!$omp parallel do default(shared) private(m,mp)
        do i=1,n_x
           do mp=1,n_s2
              do m=1,n_s2
                 h(m,i,p_nek_loc,n_spec) = h(m,i,p_nek_loc,n_spec)+&
                      o_advect(m,mp,i,p_nek_loc)*h_temp(mp,i,p_nek_loc)
              enddo ! m
           enddo ! mp
           do mp=n_s2+1,n_stack
              do m=n_s2+1,n_stack
                 h(m,i,p_nek_loc,n_spec) = h(m,i,p_nek_loc,n_spec)+&
                      o_advect(m,mp,i,p_nek_loc)*h_temp(mp,i,p_nek_loc)
              enddo ! m
           enddo ! mp
        enddo ! i
!$omp end parallel do

     else

!$omp parallel do default(shared) private(m,mp)
        do i=1,n_x
           do mp=1,n_stack
              do m=1,n_stack
                 h(m,i,p_nek_loc,n_spec) = h(m,i,p_nek_loc,n_spec)+&
                      o_advect(m,mp,i,p_nek_loc)*h_temp(mp,i,p_nek_loc)
              enddo ! m
           enddo ! mp
        enddo ! i
!$omp end parallel do

     endif

  enddo ! p_nek_loc


  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_get_delta_he done]'
  endif

end subroutine gyro_get_delta_he
