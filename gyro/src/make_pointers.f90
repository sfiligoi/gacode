!---------------------------------------------------
! make_pointers.f90
!
! PURPOSE:
!  Make pointers for domain decomposition.
!---------------------------------------------------

subroutine make_pointers

  use gyro_globals
  use gyro_pointers

  implicit none

  integer :: p


  if (linsolve_method == 3) then

     ! Keep all k (pitch-angle) on a processor

     k  = 0
     ie = 1
     i  = 0
     p  = 1
     do 
        k = k+1
        if (k > n_lambda) then
           k  = 1
           ie = ie+1
        endif
        if (ie > n_energy) exit
        nek_n(p) = in_1
        nek_e(p) = ie
        nek_k(p) = k
        p = p+n_proc
        if (p > n_nek_1) then
           i = i+1
           p = 1+i
        endif
     enddo

  else

     p = 0
     do ie=1,n_energy
        do k=1,n_lambda

           p = p+1

           nek_n(p) = in_1
           nek_e(p) = ie
           nek_k(p) = k

        enddo
     enddo

  endif


  p = 0
  do i=1,n_x
     do ie=1,n_energy

        p = p+1

        ine_i(p) = i
        ine_n(p) = in_1
        ine_e(p) = ie

     enddo
  enddo

  p = 0
  do k=1,n_lambda
     do i=1,n_x

        p = p+1

        ki_k(p) = k
        ki_i(p) = i

     enddo
  enddo


  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[make_pointers done]'
  endif

end subroutine make_pointers
