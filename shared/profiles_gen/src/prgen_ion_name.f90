subroutine prgen_ion_name(m,z,name)

  integer, intent(in) :: m
  integer, intent(in) :: z
  character (len=3), intent(inout) :: name

  select case (m)

  case (1)
     name = 'H'

  case (2)
     name = 'D'

  case (3)
     if (z == 1) then
        name = 'T'
     else if (z == 2) then
        name = 'He3'
     else
        name = '?'
     endif

  case (4)
     name = 'He4'

  case (7)
     name = 'Li'

  case (9)
     name = 'Be'

  case (12)
     name = 'C'

  case (40)
     name = 'Ar'

  case default
     name ='?'

  end select

end subroutine prgen_ion_name
