subroutine cleanup_mumps

  use gyro_mumps_private

  integer :: i

  do i=1,3
     if (associated(m_mumps(i)%A)) deallocate(m_mumps(i)%A)
     if (associated(m_mumps(i)%IRN)) deallocate(m_mumps(i)%IRN)
     if (associated(m_mumps(i)%JCN)) deallocate(m_mumps(i)%JCN)
     if (associated(m_mumps(i)%rhs)) deallocate(m_mumps(i)%rhs)
  enddo

end subroutine cleanup_mumps
