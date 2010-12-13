subroutine write_theta_operators

  use gyro_globals

  implicit none

  do k=1,n_lambda

     ck = class(k)

     print *,' tau    theta_t'
     do m=1,n_tau(ck)
        print '(2(f8.5,2x))',tau(1,k,m),theta_t(1,k,m)
     enddo
     print *

     print *,' j    m_map(j,1)  m_map(j,2) '
     do j=1,n_theta(ck)
        print '(3(i3,2x))',j,m_map(ck,j,1),m_map(ck,j,2)
     enddo
     print *

     print *,' m   m_cyc(m,1)  m_cyc(m,2) '
     do m=-n_stack+1,2*n_stack
        print '(3(i3,2x),4(f6.3,1x))',m,m_cyc(ck,m,:),p_cyc(ck,2,m,:)
     enddo
     print *

     print *,' m  m_phys  p_phys'
     do m=1,n_stack
        print '(3(i3,2x))',m,m_phys(ck,m),p_phys(ck,m)
     enddo
     print *

  enddo

end subroutine write_theta_operators
