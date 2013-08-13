 
module periodic_table_mod
  implicit none
 
  !
  ! --- a structure which describes an element ---
  !
  type, private :: ELM
    character*2 :: symbol   ! symbol for the element
    integer     :: number   ! atomic number of the element
    real        :: weight   ! atomic weight of the element
  end type ELM
 
  !
  ! --- parameters ---
  !
  integer, parameter,private  :: max_elements = 109 ! number of elements
 
  type (elm), dimension(max_elements), parameter, private :: elex = (/ &
       elm('Ac', 89,227.0000) , elm('Ag', 47,107.8680) , elm('Al', 13, 26.9815) ,  &
       elm('Am', 95,243.0000) , elm('Ar', 18, 39.9480) , elm('As', 33, 74.9216) ,  &
       elm('At', 85,210.0000) , elm('Au', 79,196.9665) , elm('B ',  5, 10.8100) ,  &
       elm('Ba', 56,137.3400) , elm('Be',  4,  9.0122) , elm('Bh',107,262.0000) ,  &
       elm('Bi', 83,208.9804) , elm('Bk', 97,247.0000) , elm('Br', 35, 79.9040) ,  &
       elm('C ',  6, 12.0110) , elm('Ca', 20, 40.0800) , elm('Cd', 48,112.4000) ,  &
       elm('Ce', 58,140.1200) , elm('Cf', 98,251.0000) , elm('Cl', 17, 35.4530) ,  &
       elm('Cm', 96,247.0000) , elm('Co', 27, 58.9332) , elm('Cr', 24, 51.9960) ,  &
       elm('Cs', 55,132.9054) , elm('Cu', 29, 63.5460) , elm('Db',105,262.0000) ,  &
       elm('Dy', 66,162.5000) , elm('Er', 68,167.2600) , elm('Es', 99,252.0000) ,  &
       elm('Eu', 63,151.9600) , elm('F ',  9, 18.9984) , elm('Fe', 26, 55.8470) ,  &
       elm('Fm',100,257.0000) , elm('Fr', 87,223.0000) , elm('Ga', 31, 69.7200) ,  &
       elm('Gd', 64,157.2500) , elm('Ge', 32, 72.5900) , elm('H ',  1,  1.0079) ,  &
       elm('He',  2,  4.0026) , elm('Hf', 72,178.4900) , elm('Hg', 80,200.5900) ,  &
       elm('Ho', 67,164.9304) , elm('Hs',108,265.0000) , elm('I ', 53,126.9045) ,  &
       elm('In', 49,114.8200) , elm('Ir', 77,192.2200) , elm('K ', 19, 39.0980) ,  &
       elm('Kr', 36, 83.8000) , elm('La', 57,138.9055) , elm('Li',  3,  6.9410) ,  &
       elm('Lr',103,262.0000) , elm('Lu', 71,174.9700) , elm('Md',101,258.0000) ,  &
       elm('Mg', 12, 24.3050) , elm('Mn', 25, 54.9380) , elm('Mo', 42, 95.9400) ,  &
       elm('Mt',109,265.0000) , elm('N ',  7, 14.0067) , elm('Na', 11, 22.9898) ,  &
       elm('Nb', 41, 92.9064) , elm('Nd', 60,144.2400) , elm('Ne', 10, 20.1790) ,  &
       elm('Ni', 28, 58.7000) , elm('No',102,259.0000) , elm('Np', 93,237.0482) ,  &
       elm('O ',  8, 15.9994) , elm('Os', 76,190.2000) , elm('P ', 15, 30.9738) ,  &
       elm('Pa', 91,231.0359) , elm('Pb', 82,207.2000) , elm('Pd', 46,106.4000) ,  &
       elm('Pm', 61,145.0000) , elm('Po', 84,209.0000) , elm('Pr', 59,140.9077) ,  &
       elm('Pt', 78,195.0900) , elm('Pu', 94,244.0000) , elm('Ra', 88,226.0254) ,  &
       elm('Rb', 37, 85.4678) , elm('Re', 75,186.2070) , elm('Rf',104,261.0000) ,  &
       elm('Rh', 45,102.9055) , elm('Rn', 86,222.0000) , elm('Ru', 44,101.0700) ,  &
       elm('S ', 16, 32.0600) , elm('Sb', 51,121.7500) , elm('Sc', 21, 44.9559) ,  &
       elm('Se', 34, 78.9600) , elm('Sg',106,263.0000) , elm('Si', 14, 28.0860) ,  &
       elm('Sm', 62,150.4000) , elm('Sn', 50,118.6900) , elm('Sr', 38, 87.6200) ,  &
       elm('Ta', 73,180.9479) , elm('Tb', 65,158.9254) , elm('Tc', 43, 97.0000) ,  &
       elm('Te', 52,127.6000) , elm('Th', 90,232.0381) , elm('Ti', 22, 47.9000) ,  &
       elm('Tl', 81,204.3700) , elm('Tm', 69,168.9342) , elm('U ', 92,238.0290) ,  &
       elm('V ', 23, 50.9414) , elm('W ', 74,183.5000) , elm('Xe', 54,131.3000) ,  &
       elm('Y ', 39, 88.9059) , elm('Yb', 70,173.0400) , elm('Zn', 30, 65.3800) ,  &
       elm('Zr', 40, 91.2200) /)
 
  character*26, parameter, private :: lower  = 'abcdefghijklmnopqrstuvwxyz'
  character*26, parameter, private :: upper  = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  character*10, parameter, private :: number = '0123456789'
 
  character*2, dimension(max_elements), parameter, private :: element = (/ &
         'H ', 'He', &
         'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
         'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', &
         'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', &
         'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', &
         'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', &
         'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', &
         'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', &
         'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
         'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', &
         'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', &
         'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U ', 'Np', 'Pu', &
         'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', &
         'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt' /)
 
  character*13, dimension(max_elements), parameter, private :: element_name = (/ &
         'Hydrogen     ', 'Helium       ', &
         'Lithium      ', 'Beryllium    ', 'Boron        ', 'Carbon       ', &
         'Nitrogen     ', &
         'Oxygen       ', 'Fluorine     ', 'Neon         ', &
         'Sodium       ', 'Magnesium    ', 'Aluminum     ', 'Silicon      ', &
         'Phosphorus   ', 'Sulfur       ', 'Chlorine     ', 'Argon        ', &
         'Potassium    ', 'Calcium      ', 'Scandium     ', 'Titanium     ', &
         'Vanadium     ', 'Chromium     ', 'Manganese    ', 'Iron         ', &
         'Cobalt       ', 'Nickel       ', 'Copper       ', 'Zinc         ', &
         'Gallium      ', 'Germanium    ', 'Arsenic      ', 'Selenium     ', &
         'Bromine      ', 'Krypton      ', &
         'Rubidium     ', 'Strontium    ', 'Yttrium      ', 'Zirconium    ', &
         'Niobium      ', 'Molybdenum   ', 'Technetium   ', 'Ruthenium    ', &
         'Rhodium      ', 'Palladium    ', 'Silver       ', 'Cadmium      ', &
         'Indium       ', 'Tin          ', 'Antimony     ', 'Tellurium    ', &
         'Iodine       ', 'Xenon        ', &
         'Cesium       ', 'Barium       ', 'Lanthanum    ', 'Cerium       ', &
         'Praseodymium ', 'Neodymium    ', 'Promethium   ', 'Samarium     ', &
         'Europium     ', 'Gadolinium   ', 'Terbium      ', 'Dysprosium   ', &
         'Holmium      ', 'Erbium       ', 'Thulium      ', 'Ytterbium    ', &
         'Lutetium     ', 'Hafnium      ', 'Tantalum     ', 'Tungsten     ', &
         'Rhenium      ', 'Osmium       ', 'Iridium      ', 'Platinum     ', &
         'Gold         ', 'Mercury      ', &
         'Thallium     ', 'Lead         ', 'Bismuth      ', 'Polonium     ', &
         'Astatine     ', 'Radon        ', &
         'Francium     ', 'Radium       ', 'Actinium     ', 'Thorium      ', &
         'Protactinium ', 'Uranium      ', 'Neptunium    ', 'Plutonium    ', &
         'Americium    ', 'Curium       ', 'Berkelium    ', 'Californium  ', &
         'Einsteinium  ', 'Fermium      ', 'Mendeleev    ', 'Nobelium     ', &
         'Lawrencium   ', 'Rutherfordium', 'Dubnium      ', 'Seaborgium   ', &
         'Bohrium      ', 'Hassium      ', 'Meitnerium   ' /)
 
contains
  !
  !
  ! -------------------------- to_periodic_table ---------------------------
  !
  ! this function returns the periodic table symbol for some element,
  ! if there is an error there will be a ? in the first position of
  ! the returned string followed by a short message
  !
  ! examples:    to_periodic_table(6,0,-1,1)       ->  'C           '
  !              to_periodic_table(6,12.0107,-1,1) ->  'C12         '
  !              to_periodic_table(6,0,4,1)        ->  'C+4         '
  !              to_periodic_table(6,12.0107,4,1)  ->  'C12+4       '
  !              to_periodic_table(6,12.0107,14,1) ->  '?Charge>Z   '
  !
  !   z = atomic number -- must be greater than 0
  !   a = atomic weight or 0
  !   charge = charge on element or -1
  !
  function to_periodic_table(z, a, charge, isymbol) result(out)
    implicit none
 
    integer, intent(in) :: z       ! the atomic number of the element
    real,    intent(in) :: a       ! the atomic weight of the element
    integer, intent(in) :: charge  ! the charge on the element -- ignore if <0
    integer, intent(in) :: isymbol ! symbol to use to separate element name and charge
                                   ! 0 -> _  this is an underscore
                                   ! 1 -> +
    character*12        :: out     ! the symbol returned
    character*12        :: swork   ! working character string
    character*1         :: csymbol ! charge symbol
 
    if (z<=0) then
       out = '? Z<=0      '
       return
    end if
 
    if (z>max_elements) then
       out = '?Z too large'
       return
    end if
 
    out = element(z)   ! form the symbol
 
    ! add in atomic weight
    if (a>0) then
       write(swork, '(I6)') nint(a)
       out = trim(out) // trim(adjustl(swork))
    end if
 
    ! add in charge state
    if (charge>=0) then
       if (charge>z) then
          out = '?Charge>Z   '
          return
       end if
       write(swork, '(I6)') charge
       if (isymbol > 0) then
          csymbol = '+'
       else
          csymbol = '_'   ! this is an underscore
       end if
       out = trim(out) // csymbol // trim(adjustl(swork))
    end if
  end function to_periodic_table
 
  !
  ! ------------------------- inv_periodic_table -------------------
  !
  ! this subroutine maps a character string returned by the function
  ! to_periodic_table to the atomic number, atomic weight and charge
  ! or zero if the data is unavailable.
  !
  ! examples:   inv_periodic_table('C',    .true.,  z, a, c) -> z=6,  a=12.0110, c=-1
  !             inv_periodic_table('C',    .false., z, a, c) -> z=6,  a=0.,      c=-1
  !             inv_periodic_table('Mg',   .true.,  z, a, c) -> z=12, a=24.3050, c=-1
  !             inv_periodic_table('C12',  .true.,  z, a, c) -> z=6,  a=12.,     c=-1
  !             inv_periodic_table('C+4',  .true.,  z, a, c) -> z=6,  a=12.0110, c=4
  !             inv_periodic_table('C14+4',.true.,  z, a, c) -> z=6,  a=14,      c=4
  !             inv_periodic_table('C14_4',.true.,  z, a, c) -> z=6,  a=14,      c=4
  !             inv_periodic_table('Q+7',  .true.,  z, a, c) -> z=0,  a=0,       c=-1
  !             inv_periodic_table('',     .true.,  z, a, c) -> z=0,  a=0,       c=-1
  !
  subroutine inv_periodic_table(symbol, default, z, a, charge)
    character*(*), intent(in) :: symbol   ! string containing element (e.g. 'C12+4 ')
                                          ! must be left justified
    logical, intent(in)       :: default  ! .true. to return an atomic weight even
                                          ! if data is not contained in the string
 
    integer, intent(out) :: z             ! atomic number of the element
    real,    intent(out) :: a             ! atomic weight of the element
    integer, intent(out) :: charge        ! charge state of the element
 
    character*2 :: as     ! atomic symbol in string
    character*10 zread    ! integer read buffer
 
    integer :: il         ! length of input string without trailing whitespace
    integer :: is         ! length of atomic symbol part in string
    integer :: j1,j2,j    ! temporary indices into element array
 
    il = len_trim(symbol)     ! length of string
    if (il == 0) goto 10      ! huh?
 
    ! --- get 2 character symbol id ---
    is = min(il,2)                                  ! length of symbol part
    as = symbol(1:is)                               ! atomic symbol
    if (index(upper,symbol(1:1)) == 0) goto 10      ! first letter is not upper case
    if (il<2 .or. index(lower,symbol(is:is)) == 0) then  ! lowercase second character?
       is=1           ! symbol has only one character
       as(2:2) = ' '  ! make sure padded with space
    end if
 
    ! --- search for element ---
    j1 = 1
    j2 = max_elements
    do while ((j2-j1)>1)            ! stop when j2=j1+1
       j = (j1+j2)/2                ! binary section search
       if (as<elex(j)%symbol) then  ! alphabetic ordering
          j2=j
       else
          j1=j
          if (as == elex(j1)%symbol) exit  ! exact match
       end if
    end do
    if (as == elex(j1)%symbol) then        ! matches j1,j2 or nothing
       j = j1
    else if (as == elex(j2)%symbol) then
       j = j2
    else
       goto 10     ! not found in list of elements
    end if
 
    ! --- found element, set defaults ---
    z = elex(j)%number
    charge = -1
    if (default) then
       a = elex(j)%weight
    else
       a = 0.
    end if
 
    ! --- get atomic weight ---
    j1 = is+1        ! first character of atomic weight
    j2 = j1          ! first character after atomic weight
    do while(j2<=il)
       if (index(number,symbol(j2:j2)) == 0) exit  ! character is not a number
       j2 = j2+1
    end do
    if (j1<j2) then
       zread=' '
       zread(10-(j2-j1)+1:10)=symbol(j1:j2-1)
       read(zread,'(i10)') j               ! convert string to atomic weight
       a = j
    end if
 
    ! --- get charge ---
    if (j2<il) then
       if (symbol(j2:j2) == '+' .or. symbol(j2:j2) == '_') then
          j1 = j2+1  ! first character of charge
          j2 = j1    ! first character after charge
          do while(j2<=il)
             if (index(number,symbol(j2:j2)) == 0) exit ! not a number
             j2 = j2+1
          end do
          if (j1<j2) then
             zread=' '
             zread(10-(j2-j1)+1:10)=symbol(j1:j2-1)
             read(zread,'(i10)') charge
          end if
       end if
    end if
 
    return  ! ok
 
    ! ----- error exit ------
10  continue
    z = 0
    a = 0.
    charge = -1
 
  end subroutine inv_periodic_table
 
  !
  !
  ! -------------------------- to_pseudo_periodic_table ---------------------------
  ! return the periodic table symbol except the isotopes of hydrogen
  ! are returned as H,D,T if ipseudo/=0
  !
  function to_pseudo_periodic_table(ipseudo, z, a, charge, isymbol) result(out)
    implicit none
 
    integer, intent(in) :: ipseudo ! 0 -> standard periodic table
                                   ! 1  -> return H,D,T for hydrogen isotopes
    integer, intent(in) :: z       ! the atomic number of the element
    real,    intent(in) :: a       ! the atomic weight of the element
    integer, intent(in) :: charge  ! the charge on the element
    integer, intent(in) :: isymbol ! symbol to use to separate element name and charge
                                   ! 0 -> _  this is an underscore
                                   ! 1 -> +
    character*12        :: cout    ! the symbol returned from to_periodic_table
    character*12        :: out     ! the symbol returned
    character*12        :: swork   ! working character string
 
    cout = to_periodic_table(z, a, charge, isymbol)
    if (ipseudo/=0 .and. cout(1:1) == 'H') then
       if (cout(2:2) == '1') then
          out = 'H'//cout(3:12)//' '
       else if (cout(2:2) == '2') then
          out = 'D'//cout(3:12)//' '
       else if (cout(2:2) == '3') then
          out = 'T'//cout(3:12)//' '
       else
          out = cout
       end if
    else
       out = cout
    end if
  end function to_pseudo_periodic_table
 
 
  !
  ! ---------------------- inv_pseudo_periodic_table --------------------
  ! return the inv_periodic_table data where the hydrogen isotopes
  ! are given by H,D,T if ipseudo/=0
  !
  subroutine inv_pseudo_periodic_table(ipseudo, symbol, default, z, a, charge)
    implicit none
 
    integer, intent(in) :: ipseudo        ! 0 -> standard periodic table
                                          ! 1 -> use H,D,T for hydrogen isotopes
    character*(*), intent(in) :: symbol   ! string containing element (e.g. 'C12+4 ')
                                          ! must be left justified
    logical, intent(in)       :: default  ! .true. to return an atomic weight even
                                          ! if data is not contained in the string
 
    integer, intent(out) :: z             ! atomic number of the element
    real,    intent(out) :: a             ! atomic weight of the element
    integer, intent(out) :: charge        ! charge state of the element
    integer              :: k             ! isotope number of hydrogen
 
    character*12          :: out,lsymbol     ! construct a new symbol string
    character*3,parameter :: first  = 'HDT'  ! allowed hydrogen first character
    character*3,parameter :: second = ' +_'  ! allowed hydrogen second character
 
    out = symbol              ! copy, length is at least 12
    lsymbol = adjustl(out)    ! left justify
    k = index(first,lsymbol(1:1))
    if (ipseudo/=0 .and. k/=0 .and. index(second,lsymbol(2:2))/=0) then
       if (k==1) then
          out = 'H1'//lsymbol(2:11)
       else if (k==2) then
          out = 'H2'//lsymbol(2:11)
       else
          out = 'H3'//lsymbol(2:11)
       end if
    end if
    call inv_periodic_table(out, default, z, a, charge)
 
  end subroutine inv_pseudo_periodic_table
 
  !
  ! ----------------------- name_periodic_table --------------
  ! return the name of the element with atomic number 'z' or
  ! a '?' if it is out of bounds -- add atomic weight to end
  ! if a>0
  ! example:
  !    name_periodic_table(20,0)     -> Calcium
  !    name_periodic_table(20,40.08) -> Calcium40
  !
  function name_periodic_table(z,a) result(out)
    integer, intent(in) :: z       ! the atomic number of the element
    real,    intent(in) :: a       ! the atomic weight of the element
    character*16        :: out     ! name of the element
    integer             :: i       ! length of trimmed name
    character*3         :: as      ! a as a string
 
    if (z>0 .and. z<=max_elements) then
       out(1:13) = element_name(z)
       out(14:16) = '   '          ! len_trim fails without this
       if (a>0 .and. a<1000) then
          i = len_trim(out)
          write(as,'(I3)') int(nint(a))
          out(i+1:i+3) = adjustl(as)
       end if
    else
       out = "?"
    end if
  end function name_periodic_table
 
 
  !
  ! ----------------------- name_pseudo_periodic_table --------------
  ! return the name of the element with atomic number 'z' or
  ! a '?' if it is out of bounds but return Deuterium, Tritium
  ! as appropriate for hydrogen
  ! example:
  !    name_periodic_table(1,20,0)     -> Calcium
  !    name_periodic_table(1,20,40.08) -> Calcium40
  !    name_periodic_table(1,1,1)      -> Hydrogen
  !    name_periodic_table(1,1,2)      -> Deuterium
  !    name_periodic_table(1,1,3)      -> Tritium
  !    name_periodic_table(0,1,3)      -> Hydrogen3
  !
  function name_pseudo_periodic_table(ipseudo,z,a) result(out)
    integer, intent(in) :: ipseudo ! 0 -> standard periodic table
                                   ! 1 -> use H,D,T for hydrogen isotopes
    integer, intent(in) :: z       ! the atomic number of the element
    real,    intent(in) :: a       ! the atomic weight of the element
    character*16        :: out     ! name of the element
    integer             :: ia      ! a as an integer
 
    if (z==1 .and. ipseudo>0) then
       ia = nint(a)
       if (ia==1) then
          out = 'Hydrogen'
       else if (ia==2) then
          out = 'Deuterium'
       else
          out = 'Tritium'
       end if
    else
       out = name_periodic_table(z,a)
    end if
  end function name_pseudo_periodic_table
 
 
  !
  ! ------------------------- standard_amu -------------------
  !
  ! this subroutine maps a charge number (z) to the corresponding
  !  atom's standard mass in amu, A; A=0.0 is returned if z is out
  !  of range.
  !
  ! examples:   call standard_amu(6,A)  --> A=12.011  (carbon)
  !             call standard_amu(26,A) --> A=55.847  (iron)
  !
 
  subroutine standard_amu(z,a)
 
    integer, intent(in) :: z     ! atomic number (in)
    real, intent(out) :: a       ! standard atomic mass, amu (out)
 
  ! -----------------
 
    character*12 tmp
    integer iz,icharge
 
  ! -----------------
 
    if( (z.le.0).or.(z.gt.max_elements) ) then
       a = 0.0
    else
       a = 0.0
       tmp = to_periodic_table(z,a,-1,1)
       call inv_periodic_table(tmp,.true.,iz, a , icharge)
    endif
 
  end subroutine standard_amu
 
end module periodic_table_mod
 
