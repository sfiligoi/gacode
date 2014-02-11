! A: sin m cos n
! B: sin m sin n
! C: cos m cos n
! D: cos m sin n

subroutine le3_basis(i0,m0,n0,basis0,btype)

  use le3_globals

  implicit none
  integer, intent(in)  :: i0,m0,n0
  character(len=2) :: btype
  real, dimension(nt,np), intent(inout) :: basis0
  integer :: kt,kp

  select case (btype)

  case ('d0')

     select case (i0)

     case (1)
        do kt=1,nt
           do kp=1,np
              basis0(kt,kp) = sinm(kt,m0)*cosn(kp,n0) 
           enddo
        enddo

     case (2) 
        do kt=1,nt
           do kp=1,np
              basis0(kt,kp) = sinm(kt,m0)*sinn(kp,n0)
           enddo
        enddo

     case (3)
        do kt=1,nt
           do kp=1,np
              basis0(kt,kp) = cosm(kt,m0)*cosn(kp,n0)
           enddo
        enddo

     case (4)
        do kt=1,nt
           do kp=1,np
              basis0(kt,kp) = cosm(kt,m0)*sinn(kp,n0)
           enddo
        enddo

     end select

  case ('dt')

     select case (i0)

     case (1)
        do kt=1,nt
           do kp=1,np
              basis0(kt,kp) = m0 * cosm(kt,m0)*cosn(kp,n0) 
           enddo
        enddo

     case (2)
        do kt=1,nt
           do kp=1,np
              basis0(kt,kp) = m0 * cosm(kt,m0)*sinn(kp,n0)
           enddo
        enddo

     case (3)
        do kt=1,nt
           do kp=1,np
              basis0(kt,kp) = -m0 * sinm(kt,m0)*cosn(kp,n0)
           enddo
        enddo

     case (4)
        do kt=1,nt
           do kp=1,np
              basis0(kt,kp) = -m0 * sinm(kt,m0)*sinn(kp,n0)
           enddo
        enddo

     end select

  case ('dp')

     select case (i0)

     case (1)

        do kt=1,nt
           do kp=1,np
              basis0(kt,kp) = -n0 * sinm(kt,m0)*sinn(kp,n0) 
           enddo
        enddo

     case (2)
        do kt=1,nt
           do kp=1,np
              basis0(kt,kp) = n0 * sinm(kt,m0)*cosn(kp,n0)
           enddo
        enddo

     case (3)
        do kt=1,nt
           do kp=1,np
              basis0(kt,kp) = -n0 * cosm(kt,m0)*sinn(kp,n0)
           enddo
        enddo

     case (4)
        do kt=1,nt
           do kp=1,np
              basis0(kt,kp) = n0 * cosm(kt,m0)*cosn(kp,n0)
           enddo
        enddo

     end select

  end select

end subroutine le3_basis
