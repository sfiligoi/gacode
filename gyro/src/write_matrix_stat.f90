subroutine write_matrix_stat(nelem,nval,nindx,niter,tag)

  use gyro_globals

  !------------------------------------------------
  implicit none
  !
  integer, intent(in) :: nelem
  integer, intent(in) :: nval
  integer, intent(in) :: nindx
  integer, intent(in) :: niter
  !
  integer :: nelem_c(n_n)
  integer :: nval_c(n_n)
  integer :: nindx_c(n_n)
  integer :: niter_c(n_n)
  integer :: p
  !
  integer, intent(in) :: tag
  !
  character (len=36) :: lab1
  character (len=36) :: lab2
  character (len=36) :: lab3
  character (len=36) :: lab4
  !------------------------------------------------

  select case (output_flag)

  case (1)

     call collect_integer(nelem,nelem_c)
     call collect_integer(nval,nval_c)
     call collect_integer(nindx,nindx_c)
     call collect_integer(niter,niter_c)

     if (n(1) == 0 .and. n_n > 1) then
        p = 2
        lab1 = 'EXPLICIT POISSON:     n=0        n>0'
        lab2 = ' EXPLICIT AMPERE:     n=0        n>0'
        lab3 = '     TOTAL FIELD:     n=0        n>0'
        lab4 = '  POISSON-AMPERE:     n=0        n>0'
     endif

     if (n(1) == 0 .and. n_n == 1) then
        p = 1
        lab1 = 'EXPLICIT POISSON:     n=0'
        lab2 = ' EXPLICIT AMPERE:     n=0'
        lab3 = '     TOTAL FIELD:     n=0'
        lab4 = '  POISSON-AMPERE:     n=0'
     endif

     if (n(1) /= 0) then
        p = 1
        lab1 = 'EXPLICIT POISSON:     n>0'
        lab2 = ' EXPLICIT AMPERE:     n>0'
        lab3 = '     TOTAL FIELD:     n>0'
        lab4 = '  POISSON-AMPERE:     n>0'
     endif

     if (i_proc == 0) then
        open(unit=1,file=trim(runfile),status='old',position='append')

        select case (tag)

        case (1)

           write(1,*) '----------- SPARSE MATRIX STATS ---------------'
           write(1,*) lab1

        case (2)

           write(1,*) lab2

        case (3)

           write(1,*) lab3

        case (4)

           write(1,*) '----------- SPARSE MATRIX STATS ---------------'
           write(1,*) lab4

        end select

        write(1,40) '        nonzeros:',nelem_c(1:p)
        write(1,40) '          values:',nval_c(1:p)
        write(1,40) '         indices:',nindx_c(1:p)
        write(1,40) '      iterations:',niter_c(1:p)
        write(1,*)

        close(1)

     endif

     return

  end select

40 format(t2,a,t20,4(i8,2x))

end subroutine write_matrix_stat
