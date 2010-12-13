subroutine BLEND_cleanup

  use BLEND_private
  
  implicit none

  deallocate(i_piv)
  deallocate(cs)
  deallocate(c0)

end subroutine BLEND_cleanup
