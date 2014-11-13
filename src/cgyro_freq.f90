module cgyro_freq
  
  implicit none
  
  public :: FREQ_alloc, FREQ_do

  logical, private :: initialized = .false.
  complex, dimension(:,:), allocatable, private :: freq_loc
  real, dimension(:,:), allocatable, private :: mode_weight
  integer, parameter, private :: io_freq = 1
  character(len=80),private :: runfile_freq = 'out.cgyro.freq'

  complex :: freq
  complex :: freq_err
  
contains
  
  subroutine FREQ_alloc(flag)
    use cgyro_globals
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

    use cgyro_globals

    implicit none
  
    real :: df_r, df_i, total_weight

    ! Use potential as gauge for frequency
    mode_weight(:,:) = abs(field(:,:,1))

    ! Define local frequencies
    freq_loc(:,:) = (i_c/delta_t)*log(field(:,:,1)/field_old(:,:,1))

    total_weight = sum(mode_weight(:,:))

    freq = sum(freq_loc(:,:)*mode_weight(:,:))/total_weight

    ! Fractional Frequency Error
    df_r = sum(abs(real(freq_loc(:,:)-freq))*mode_weight(:,:))
    df_i = sum(abs(imag(freq_loc(:,:)-freq))*mode_weight(:,:))

    freq_err = (df_r + i_c*df_i)/total_weight/abs(freq)

    if (silent_flag == 0 .and. i_proc == 0) then
       open(unit=io_freq,file=trim(path)//runfile_freq,status='old',&
            position='append')
       write(io_freq,20) freq,freq_err
       close(io_freq)
       print '(t2,1pe10.3,2x,2(1pe13.6,1x),2x,2(1pe10.3,1x))', itime*delta_t,freq,freq_err
    endif

20  format(4(es11.4,1x))

  end subroutine FREQ_do

end module cgyro_freq
