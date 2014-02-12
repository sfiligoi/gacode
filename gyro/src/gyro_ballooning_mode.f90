!------------------------------------------------------
! gyro_ballooning_mode.f90 [caller gyro_write_master]
!
! PURPOSE:
!  This routine decomposes real-space functions f(i,j) 
!  into ballooning mode p-harmonics, f_bar(p,j), patches 
!  together the full ballooning mode f_*(theta_*) from 
!  these p-harmonics, and then writes to disk.
!------------------------------------------------------

subroutine gyro_ballooning_mode(datafile,io,index,is_in)

  use gyro_globals
  use math_constants

  !--------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  !
  integer, intent(in) :: io
  integer, intent(in) :: index
  integer, intent(in) :: is_in
  !
  integer :: p
  integer :: pp
  integer :: l0
  integer :: np
  integer :: j_renorm(1)
  integer :: data_loop
  integer :: ampCount
  !
  complex, dimension(-n_x/2:n_x/2-1,n_theta_plot) :: f_bar
  complex, dimension(n_theta_plot,n_x) :: fplot
  complex, allocatable                 :: amp(:)
  complex :: dummy
  !--------------------------------------------------

  !------------------------------------------------------
  ! Need to define these parameters before case selection
  !
  m0 = int(box_multiplier)
  np = n_x/m0
  !------------------------------------------------------

  !----------------------------------------------------
  ! Look for acceptable in_1-value only on processor 0:
  ! 
  if (i_proc > 0) then

     return

  else

     if (n_1(in_1) == 0) return

  endif
  !
  !----------------------------------------------------

  select case (io_control)

  case(0)

     return

  case(1)

     ! Open

     open(unit=io,file=datafile,status='replace')
     close(io)

  case(2)

     ! Append

     open(unit=io,file=datafile,status='old',position='append')

     !----------------------------------------------------
     ! Select function to map to ballooning space:
     !
     select case (index)

     case (1) 
        fplot(:,:) = phi_plot(:,:,1)
     case (2) 
        fplot(:,:) = phi_plot(:,:,2)
     case (3)
        fplot(:,:) = phi_plot(:,:,3)
     case (4)
        fplot(:,:) = phi_plot(:,:,4)
     case (5) 
        fplot(:,:) = moments_plot(:,:,is_in,1)
     case (6) 
        fplot(:,:) = moments_plot(:,:,is_in,2)
     case (7) 
        fplot(:,:) = moments_plot(:,:,is_in,3)

     end select
     !-------------------------------------------------------

     !-------------------------------------------------------
     ! Find Fourier transform of radial-grid wave function
     !
     f_bar = (0.0,0.0)
     !
     if (q_s(ir_norm)*shat_s(ir_norm)*n_1(1) > 0.0) then
        ! n*q*s > 0
        do i=1,n_x 
           do ip=1,n_x  
              f_bar(i-1-n_x/2,:) = f_bar(i-1-n_x/2,:)+fplot(:,ip)*cri(i,ip)
           enddo
        enddo
     else
        ! n*q*s < 0
        do i=1,n_x 
           do ip=1,n_x  
              f_bar(i-1-n_x/2,:) = f_bar(i-1-n_x/2,:)+fplot(:,ip)*conjg(cri(i,ip))
           enddo
        enddo
     endif
     !-------------------------------------------------------

     !-------------------------------------------------------
     ! Renormalize function:
     !
     j_renorm = maxloc(abs(phi_plot(:,ir_norm,1)))
     balloon_renorm = 1.0/phi_plot(j_renorm(1),ir_norm,1)
     f_bar(:,:) = balloon_renorm*f_bar(:,:)
     !-------------------------------------------------------

     !-------------------------------------------------------
     ! Sew together mode in extended angle and write to file.
     ! 
     allocate(amp(m0*np*n_theta_plot))
     ampCount = 1
     do l0=0,m0-1
        do pp=-np/2,np/2-1
           p = m0*pp+l0
           do j_int=1,n_theta_plot

              ! Data output is ordered properly in extended angle:
              !
              !   theta_extended = theta_plot(j_int)+2*pi*pp

              amp(ampCount) = f_bar(p,j_int)*phase(in_1,ir_norm)**(real(p)/real(m0))
              write(io,fmtstr2) amp(ampCount)
              ampCount=ampCount+1
           enddo
        enddo
     enddo
     deallocate(amp)

     !-------------------------------------------------------

  case(3)

     ! Rewind

     open(unit=io,file=datafile,status='old')

     ! Mimic "reconstruct psi_star" loop

     do data_loop=0,data_step
        do l0=0,m0-1
           do pp=-np/2,np/2-1
              do j_int=1,n_theta_plot
                 read(io,fmtstr2) dummy 
              enddo
           enddo
        enddo
     enddo

     endfile(io)
     close(io)

  end select

  if (i_proc == 0 .and. debug_flag == 1) then
     print *,'[gyro_balloning_mode done]' 
  endif

end subroutine gyro_ballooning_mode
