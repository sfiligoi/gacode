subroutine read_time(dir)

  use peek_globals
  use gyro_globals

  implicit none

  real :: t_window
  real :: tmin,tmax,tmid
  character (len=*) :: dir

  open(unit=1,file=trim(dir)//'t.out',iostat=i_err)
  i_err = 0
  do while (i_err == 0)
     read(1,*,iostat=i_err) j,t_current
  enddo
  close(1)

  n_time = j

  allocate(t(0:n_time))

  open(unit=1,file=trim(dir)//'t.out')
  do i=0,n_time
     read(1,*) j,t(i)
  enddo
  close(1)

  allocate(gbflux(n_kinetic,n_field,4))
  allocate(gbflux_1(n_kinetic,n_field,4))
  allocate(gbflux_2(n_kinetic,n_field,4))
  allocate(gbflux_t(n_kinetic,n_field,4,0:n_time))

  open(unit=1,file=trim(dir)//'out.gyro.gbflux')
  do i=0,n_time
     read(1,'(100(1pe15.8,1x))') gbflux_t(:,:,:,i)
  enddo
  close(1)

  tmin = (1.0-window)*t(n_time)
  tmax = t(n_time)
  tmid = tmin+0.5*(tmax-tmin)

  !-------------------------------------------------------------------------
  ! Compute average
  t_window = 0.0
  gbflux   = 0.0
  do i=0,n_time-1
     if (t(i) > tmin) then
        gbflux = gbflux + 0.5*(gbflux_t(:,:,:,i)+gbflux_t(:,:,:,i+1)) &
             *(t(i+1)-t(i))
        t_window = t_window+(t(i+1)-t(i))
     endif
  enddo
  gbflux = gbflux/t_window
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! Compute error
  ! 
  ! interval 1: tmid,tmax
  t_window = 0.0
  gbflux_1 = 0.0
  do i=0,n_time-1
     if (t(i) > tmid) then
        gbflux_1 = gbflux_1 + 0.5*(gbflux_t(:,:,:,i)+gbflux_t(:,:,:,i+1)) &
             *(t(i+1)-t(i))
        t_window = t_window+(t(i+1)-t(i))
     endif
  enddo
  gbflux_1 = gbflux_1/t_window
  !
  ! interval 2: tmin,tmid
  t_window = 0.0
  gbflux_2 = 0.0
  do i=0,n_time-1
     if (t(i) > tmin .and. t(i) < tmid) then
        gbflux_2 = gbflux_2 + 0.5*(gbflux_t(:,:,:,i)+gbflux_t(:,:,:,i+1)) &
             *(t(i+1)-t(i))
        t_window = t_window+(t(i+1)-t(i))
     endif
  enddo
  gbflux_2 = gbflux_2/t_window
  !
  ! Error stored in gbflux_1:
  gbflux_1 = 0.5*abs(gbflux_2-gbflux_1)
  !-------------------------------------------------------------------------

  print 30,'GyroBohm fluxes','t_min=',(1.0-window)*t(n_time),'t_max=',t(n_time)
  print *
  do i_field=1,n_field
     print *,tag(i_field)
     print 20,'density','energy','momentum','exchange'
     do i_spec=1,n_kinetic
        print 10,'species: ',i_spec,gbflux(i_spec,i_field,:)
        print 10,'(error): ',i_spec,gbflux_1(i_spec,i_field,:)
     enddo
     print *
  enddo

  print *,'Total'
  print 20,'density','energy','momentum','exchange'
  do i_spec=1,n_kinetic
     print 10,'species: ',i_spec,(sum(gbflux(i_spec,:,i)),i=1,4)
     print 10,'(error): ',i_spec,(sum(gbflux_1(i_spec,:,i)),i=1,4)
  enddo
  print *

10 format(t2,a,i1,2x,4(1pe12.5,1x))
20 format(t15,a,t28,a,t41,a,t54,a)
30 format(t2,a,2x,a,f6.1,2x,a,f6.1)

end subroutine read_time
