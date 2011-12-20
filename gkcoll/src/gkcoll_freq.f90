module gkcoll_freq
  
  implicit none
  
  public :: FREQ_alloc, FREQ_do

  logical, private :: initialized = .false.
  complex, dimension(:,:), allocatable, private :: freq_loc
  real, dimension(:,:), allocatable, private :: mode_weight
  integer, parameter, private :: io_freq = 1
  character(len=80),private :: runfile_freq = 'out.gkcoll.freq'

  complex :: freq
  complex :: freq_err
  
contains
  
  subroutine FREQ_alloc(flag)
    use gkcoll_globals
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    
    if(flag == 1) then
       if(initialized) return
       
       allocate(freq_loc(n_radial, n_theta))
       allocate(mode_weight(n_radial, n_theta))
       
        if(silent_flag == 0 .and. i_proc == 0) then
           open(unit=io_freq,file=trim(path)//runfile_freq,status='replace')
           close(io_freq)
        endif

       initialized = .true.
       
    else
       if(.NOT. initialized) return
       
       deallocate(freq_loc)
       deallocate(mode_weight)
       
       initialized = .false.
    endif
    
  end subroutine FREQ_alloc

  subroutine FREQ_do
    use gkcoll_globals
    implicit none
    integer :: ir, it
    real :: df_r, df_i, total_weight
    
    ! Frequency
    total_weight = 0.0
    freq = (0.0,0.0)
    if(minval(abs(phi_old)) == 0.0) then
       freq     = 0.0
       freq_err = 1.0
    else
       do ir=1,n_radial
          do it=1,n_theta
             mode_weight(ir,it) = abs(phi(ir,it))
             total_weight = total_weight + mode_weight(ir,it)
             freq_loc(ir,it) = (i_c/delta_t)*log(phi(ir,it)/phi_old(ir,it))
             freq = freq + freq_loc(ir,it) * mode_weight(ir,it)
          enddo
       enddo
       freq = freq / total_weight
       
       ! Fractional Frequency Error
       df_r = 0.0
       df_i = 0.0
       do ir=1,n_radial
          do it=1,n_theta
             df_r = df_r + abs(real(freq_loc(ir,it)-freq)) * mode_weight(ir,it)
             df_i = df_i + abs(imag(freq_loc(ir,it)-freq)) * mode_weight(ir,it)
          enddo
       enddo
       freq_err = (df_r + i_c * df_i) / total_weight / abs(freq)
    endif

    if(silent_flag == 0 .and. i_proc == 0) then
       open(unit=io_freq,file=trim(path)//runfile_freq,status='old',&
            position='append')
       write(io_freq,20) freq, freq_err
       close(io_freq)
       print *, freq, freq_err
    endif

20  format(4(es11.4,1x))

  end subroutine FREQ_do

end module gkcoll_freq
