subroutine prgen_ion_name(m,z,name)

  integer, intent(in) :: m
  integer, intent(in) :: z
  character (len=*), intent(inout) :: name

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
        name = 'LUMPED'
     endif

  case (4)
     ! 2021/02/15: Changed alpha name to He from He4
     name = 'He'

  case (7)
     name = 'Li'

  case (9)
     name = 'Be'

  case (11)
     ! 2021/04/12: Added Boron
     name = 'B'

  case (12)
     name = 'C'

  case (14)
     name = 'N'

  case (16)
     name = 'O'

  case (20)
     name = 'Ne'

  case (40)
     name = 'Ar'

  case (58:64)
     name = 'Ni'

  case (180:188)
     name = 'W'

  case default
     name = 'LUMPED'

  end select
  
end subroutine prgen_ion_name
