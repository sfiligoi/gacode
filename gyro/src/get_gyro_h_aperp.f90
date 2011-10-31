!-----------------------------------------------------
! get_gyro_h_aperp.f90
!
! PURPOSE:
!  Given h, compute the gyroaverage G_aperp=(1/2)(J0+J2) 
!  (gyro_h_aperp) for IONS only.  For electrons the 
!  gyroaverage is the identity operation * 0.5
!
!  This is valid for periodic OR nonperiodic boundary
!  conditions. 
!
!  NOTE:
!   This routine is expensive so some coding is inelegant 
!   for the sake of speed.
!-----------------------------------------------------

subroutine get_gyro_h_aperp

  use gyro_globals
  use gyro_pointers

  !-------------------------------------------------------
  implicit none
  !
  complex :: temp(n_stack)
  complex, dimension(n_stack,i1_buffer:i2_buffer) :: hh
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! Gyrokinetic species:
  !
  p_nek_loc = 0
  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek)  
     k  = nek_k(p_nek)   

     ck = class(k)

     do is=1,n_gk

        if (boundary_method == 1) then
           do i=1-m_gyro,0
              hh(:,i) = h(:,i+n_x,p_nek_loc,is)
           enddo
           do i=n_x+1,n_x+m_gyro-i_gyro
              hh(:,i) = h(:,i-n_x,p_nek_loc,is)
           enddo
        else
           do i=1-m_gyro,0
              hh(:,i) = 0.0
           enddo
           do i=n_x+1,n_x+m_gyro-i_gyro
              hh(:,i) = 0.0
           enddo
        endif
        do i=1,n_x
           hh(:,i) = h(:,i,p_nek_loc,is)
        enddo

        do i=1,n_x
           temp(:) = (0.0,0.0)
           do i_diff=-m_gyro,m_gyro-i_gyro
              do m=1,n_stack
                 temp(m) = temp(m) + &
                      w_gyro_aperp(m,i_diff,i,p_nek_loc,is)*hh(m,i+i_diff)

              enddo ! m
           enddo ! i_diff
           gyro_h_aperp(:,i,p_nek_loc,is) = temp(:)
        enddo ! i

     enddo ! is
  enddo ! p_nek
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! Drift-kinetic electrons:
  !
  ! (G_aperp h) = 1/2 h.
  !
  if (electron_method == 2) then
     gyro_h_aperp(:,:,:,n_spec) = 0.5*h(:,:,:,n_spec)
  endif
  !-------------------------------------------------------

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[get_gyro_h_aperp done]'
  endif

end subroutine get_gyro_h_aperp
