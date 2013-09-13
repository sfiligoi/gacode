MODULE Vector_class          ! Original Taken from Atkins book  filename: class_Vector.f90
  ! Modified for use in transport code
  USE nrtype,   ONLY : DP,I4B

  IMPLICIT NONE
  ! public, everything by default, but can specify any
  TYPE Vector
     !PRIVATE                          !means only code in this module can 
                                       !access the vector type
     INTEGER(I4B)                     :: size   ! vector length
     REAL(DP), POINTER, DIMENSION(:)  :: DATA   ! component values
  END TYPE Vector

  !              Overload common operators
  ! ---------------------------------------------------------------------
  INTERFACE OPERATOR (+)                   ! add others later 
     MODULE PROCEDURE add_Vector, add_Real_to_Vector 
  END INTERFACE
  INTERFACE OPERATOR (-)                   ! add unary versions later
     MODULE PROCEDURE subtract_Vector, subtract_Real 
  END INTERFACE
  INTERFACE OPERATOR (*)                   ! overload *
     MODULE PROCEDURE dot_Vector, real_mult_Vector, Vector_mult_real
  END INTERFACE
  INTERFACE ASSIGNMENT (=)                 ! overload =
     MODULE PROCEDURE equal_Real 
  END INTERFACE
  INTERFACE OPERATOR (==)                  ! overload ==
     MODULE PROCEDURE is_equal_to 
  END INTERFACE
! -------------------------------------------------------------------------


CONTAINS  ! functions & operators

  FUNCTION assign_Vector (values) RESULT (name)   ! array to vector constructor
    REAL(DP), INTENT(in) :: values(:)             ! given rank 1 array
    INTEGER(I4B)         :: length                ! array size
    TYPE (Vector)        :: name                  ! Vector to create
    length = SIZE(values) ; ALLOCATE ( name%data(length) ) 
    name % size = length  ; name % DATA = values 
  END FUNCTION assign_Vector

  FUNCTION assignrt_Vector (values,rt) RESULT (name)  ! array to rt*vector constructor
    REAL(DP), INTENT(in) :: values(:)             ! given rank 1 array
    REAL(DP), INTENT(in) :: rt                    ! multiplier
    INTEGER(I4B)         :: length                ! array size
    TYPE (Vector)        :: name                  ! Vector to create
    length = SIZE(values) ; ALLOCATE ( name%data(length) ) 
    name % size = length  ; name % DATA = values 
    name%data   = rt * name%data
  END FUNCTION assignrt_Vector

  FUNCTION LOAD_VECTOR (values) RESULT (name)     ! array to vector constructor
    REAL(DP), INTENT(in) :: values(:)             ! given rank 1 array
    INTEGER(I4B)         :: length                ! array size
    TYPE (Vector)        :: name                  ! Vector to load
    length = SIZE(values) 
    name % size = length  ; name % DATA = values 
  END FUNCTION LOAD_VECTOR

  FUNCTION ZERO_VECTOR (n) RESULT (name)     ! zero vector constructor
    INTEGER(I4B), INTENT(in) :: n            ! array size
    TYPE (Vector)            :: name         ! Vector to create
    ALLOCATE ( name%data(n) )
    name % size = n  ; name % DATA = 0.0_DP 
  END FUNCTION ZERO_VECTOR

  FUNCTION make_Vector (len, values) RESULT(v)      ! Optional Constructor
    INTEGER(I4B), OPTIONAL, INTENT(in) :: len       ! number of values
    REAL(DP),    OPTIONAL, INTENT(in) :: values(:)  ! given values
    TYPE (Vector)                 :: v
    IF ( PRESENT (len) ) THEN                       ! create vector data
       v%size = len ; ALLOCATE ( v%data(len) )
       IF ( PRESENT (values)) THEN ; v%data = values      ! vector
       ELSE                      ; v%data = 0.d0          ! null vector
       END IF ! values present
    ELSE                                      ! scalar constant
       v%size = 1_I4B                  ; ALLOCATE ( v%data(1) ) ! default
       IF ( PRESENT (values)) THEN ; v%data(1) = values(1)  ! scalar
       ELSE                      ; v%data(1) = 0._DP       ! null 
       END IF ! value present
    END IF ! len present
  END FUNCTION make_Vector

  FUNCTION make_uniform_0_1_Vector (len) RESULT(v) ! Optional Constructor
    INTEGER(I4B), INTENT(in) :: len       ! number of values
    TYPE (Vector)                 :: v
    INTEGER(I4B) j
    REAL(DP)  dv
    v = make_Vector(len)
    dv = 1._DP/(len-1._DP)
    v%data(1) = 0.0
    v%data(len) =1._DP
    DO j=2,len-1
       v%data(j) = v%data(j-1)+ dv
    ENDDO
  END FUNCTION make_uniform_0_1_Vector

  FUNCTION add_Real_to_Vector (v, r) RESULT (NEW)  ! overload +
    TYPE (Vector), INTENT(in) :: v
    REAL(DP),          INTENT(in) :: r
    TYPE (Vector)                 :: NEW         ! new = v + r
    IF ( v%size < 1 ) STOP "No sizes in add_Real_to_Vector"
    new%size = v%size
    ALLOCATE ( new%data(v%size) ) ; 
     new%data(:) = v%data(:) + r       ! as array operation, or use implied loop
  END FUNCTION add_Real_to_Vector

  FUNCTION add_Vector (a, b) RESULT (NEW)    ! vector + vector
    TYPE (Vector), INTENT(in) :: a, b
    TYPE (Vector)             :: NEW         ! new = a + b
    IF ( a%size /= b%size ) STOP "Sizes differ in add_Vector"
    ALLOCATE ( new%data(a%size) ) ; new%size = a%size
    new%data = a%data + b%data    
  END FUNCTION add_Vector

  FUNCTION copy_Vector (name) RESULT (NEW)
    TYPE (Vector), INTENT(in) :: name
    TYPE (Vector)             :: NEW 
    ALLOCATE ( new%data(name%size) ) ; new%size = name%size
    new%data = name%data             
  END FUNCTION copy_Vector
     
  SUBROUTINE  delete_Vector (name)                ! deallocate allocated items
    TYPE (Vector), INTENT(inout) :: name
    INTEGER(I4B)                      :: ok       ! check deallocate status 
    DEALLOCATE (name%data, stat = ok )
    !IF (ok /=0) &
     !   STOP "Vector not allocated in delete_Vector"
     name%size = 0 
  END SUBROUTINE delete_Vector

 
 FUNCTION delete_Vector_nf (name) RESULT(ok)   ! deallocate possibly allocated items
    TYPE (Vector), INTENT(inout) :: name
    INTEGER(I4B)                      :: ok     ! non fatal
    DEALLOCATE (name%data, stat = ok ) 
    name%size = 0 
  END FUNCTION delete_Vector_nf


  FUNCTION dot_Vector (a, b) RESULT (c)      ! overload *
    TYPE (Vector), INTENT(in) :: a, b
    REAL(DP)                      :: c     
    IF ( a%size /= b%size ) STOP "Sizes differ in dot_Vector"
    c = DOT_PRODUCT (a%data, b%data) 
  END FUNCTION dot_Vector

  FUNCTION crossprod_Vector (a, b,scalar) RESULT (c)      ! a transpose times b 
    TYPE (Vector), INTENT(in) :: a, b
    TYPE (Vector)             :: c
    INTEGER(I4B) j
    REAL(DP),OPTIONAL,INTENT(IN) ::  scalar
    IF ( a%size /= b%size ) STOP "Sizes differ in dot_Vector"
    ALLOCATE ( c%data(a%size) ) ; c%size = a%size
    DO j=1, c%size
       c%data(j) = a%data(j)*b%data(j)
       IF(PRESENT(scalar))c%data(j) = c%data(j)*scalar
    ENDDO
  END FUNCTION crossprod_Vector

  SUBROUTINE equal_Real (NEW, R)        ! overload =, real to vector
    TYPE (Vector), INTENT(inout) :: NEW 
    REAL(DP),          INTENT(in)    :: R   
    IF ( ASSOCIATED (new%data) ) DEALLOCATE (new%data)
    ALLOCATE ( new%data(1) ); new%size = 1
    new%data = R            
  END SUBROUTINE equal_Real

  LOGICAL FUNCTION is_equal_to (a, b) RESULT (t_f)   ! overload ==
    TYPE (Vector), INTENT(in) :: a, b                ! left & right of ==
    t_f = .FALSE.                                    ! initialize
    IF ( a%size /= b%size ) RETURN                   ! same size ?
    t_f = ALL ( a%data == b%data )                   ! and all values match
  END FUNCTION is_equal_to

  FUNCTION length_Vector (name) RESULT (n)          ! accessor member
    TYPE (Vector), INTENT(in) :: name
    INTEGER(I4B)                   :: n
    n = name % size 
  END FUNCTION length_Vector

  SUBROUTINE list (name)                     ! accessor member
    TYPE (Vector), INTENT(in) :: name
    PRINT *, "[", name % DATA(1:name%size), "]" 
  END SUBROUTINE list

  FUNCTION normalize_Vector (name)  RESULT (NEW)
    TYPE (Vector), INTENT(in) :: name
    TYPE (Vector)             :: NEW 
    REAL(DP)                  :: total, nil = EPSILON(1.0)     ! tolerance
    ALLOCATE ( new%data(name%size) ) ; new%size = name%size
    total = SQRT ( SUM ( name%data**2 ) )           ! intrinsic functions
    IF ( total < nil ) THEN ; new%data = 0.d0       ! avoid division by 0
    ELSE                  ; new%data = name%data/total
    END IF
  END FUNCTION normalize_Vector

  SUBROUTINE read_Vector (name)              ! read array, assign 
    TYPE (Vector), INTENT(inout) :: name
    INTEGER(I4B), PARAMETER           :: max = 999
    INTEGER(I4B)                      :: length
    READ (*,'(i1)', advance = 'no') length
    IF ( length <= 0 )   STOP "Invalid length in read_Vector"
    IF ( length >= max ) STOP "Maximum length in read_Vector"
    ALLOCATE ( name % DATA(length) ) ; name % size = length 
    READ *, name % DATA(1:length)    
  END SUBROUTINE read_Vector

  FUNCTION real_mult_Vector (r, v) RESULT (NEW) ! overload *
    REAL(DP),          INTENT(in) :: r
    TYPE (Vector), INTENT(in) :: v
    TYPE (Vector)             :: NEW         ! new = r * v
    IF ( v%size < 1 ) STOP "Zero size in real_mult_Vector"
    ALLOCATE ( new%data(v%size) ) ; new%size = v%size
    new%data = r * v%data         
  END FUNCTION real_mult_Vector

  FUNCTION size_Vector (name) RESULT (n)     ! accessor member
    TYPE (Vector), INTENT(in) :: name
    INTEGER(I4B)                   :: n
    n = name % size 
  END FUNCTION size_Vector

  FUNCTION subtract_Real (v, r) RESULT (NEW) ! vector - real, overload -
    TYPE (Vector), INTENT(in) :: v
    REAL(DP),          INTENT(in) :: r
    TYPE (Vector)             :: NEW         ! new = v + r
    IF ( v%size < 1 ) STOP "Zero length in subtract_Real"
    ALLOCATE ( new%data(v%size) ) ; new%size = v%size
    new%data = v%data - r         
  END FUNCTION subtract_Real

  FUNCTION subtract_Vector (a, b) RESULT (NEW) ! overload -
    TYPE (Vector), INTENT(in) :: a, b
    TYPE (Vector)             :: NEW  
    IF ( a%size /= b%size ) STOP "Sizes differ in subtract_Vector"
    ALLOCATE ( new%data(a%size) ) ; new%size = a%size
     new%data = a%data - b%data    
  END FUNCTION subtract_Vector

  FUNCTION get_values (name) RESULT (array)        ! accessor member
    TYPE (Vector), INTENT(in) :: name
    REAL(DP)                      :: array(name%size)
    array = name % DATA 
  END FUNCTION get_values

  FUNCTION get_element(name,j) RESULT (elmt)        ! accessor member
    TYPE (Vector), INTENT(in) :: name
    INTEGER,INTENT(IN) :: j
    REAL(DP)                      :: elmt
    elmt = name % DATA (j)
  END FUNCTION get_element


  FUNCTION new_Vector (length, values) RESULT(name)  ! Public constructor
    INTEGER(I4B),      INTENT(in) :: length          ! array size
    REAL(DP), TARGET, INTENT(in)  :: values(length)  ! given array
    REAL(DP), POINTER             :: pt_to_val(:)    ! pointer to array   
    TYPE (Vector)                 :: name            ! Vector to create     
    INTEGER(I4B)                  :: get_a           ! allocate flag
    ALLOCATE ( pt_to_val (length), stat = get_a )    ! allocate
    IF ( get_a /= 0 ) STOP 'allocate error'          ! check
    pt_to_val = values                         ! dereference values
    name      = Vector(length, pt_to_val)      ! intrinsic constructor
  END FUNCTION new_Vector

  FUNCTION Vector_max_value (a) RESULT (v)       ! accessor member
    TYPE (Vector), INTENT(in) :: a
    REAL(DP)                      :: v
    v = MAXVAL ( a%data(1:a%size) ) 
  END FUNCTION Vector_max_value

  FUNCTION Vector_min_value (a) RESULT (v)       ! accessor member
    TYPE (Vector), INTENT(in) :: a
    REAL(DP)                      :: v
    v = MINVAL ( a%data(1:a%size) ) 
  END FUNCTION Vector_min_value

  SUBROUTINE Vector_loc_min (a,jmin,minval)      

    TYPE (Vector), INTENT(in)   :: a
    INTEGER(I4B),INTENT(OUT)    :: jmin
    REAL(DP),INTENT(OUT)        :: minval
    INTEGER(I4b) j
    REAL(DP) val

    minval = HUGE(1._DP)
    DO j  = 1, a%size
       val = a%data(j)
       minval = MIN(minval,val)
       IF(ABS(minval-val) .lt. 2._DP*SPACING(val)) jmin = j
    ENDDO
    minval = a%data(jmin)
    RETURN

  END SUBROUTINE Vector_loc_min



  SUBROUTINE Vector_loc_max (a,jmax,maxval)       ! Portland problem with minloc,maxloc
    TYPE (Vector), INTENT(in)   :: a
    INTEGER(I4B),INTENT(OUT)    :: jmax
    REAL(DP),INTENT(OUT)        :: maxval
    INTEGER(I4b) j
    REAL(DP) val
    maxval = -HUGE(1._DP)
    DO j  = 1, a%size
       val = a%data(j)
       maxval = MAX(maxval,val)
       IF(ABS(maxval-val) .lt. 2._DP*SPACING(val)) jmax = j
    ENDDO
    maxval = a%data(jmax)
    RETURN
  END SUBROUTINE Vector_loc_max

  FUNCTION Vector_mult_real (v, r) RESULT (NEW)  ! vector*real, overload *
    TYPE (Vector), INTENT(in) :: v
    REAL(DP),          INTENT(in) :: r
    TYPE (Vector)             :: NEW             ! new = v * r
    IF ( v%size < 1 ) STOP "Zero size in Vector_mult_real"    
    NEW = Real_mult_Vector (r, v) 
  END FUNCTION Vector_mult_real

  SUBROUTINE Vector_mult (a,mult)                ! in place multiplication
    TYPE (Vector), INTENT(inout)   :: a
    REAL(DP),INTENT(in)            :: mult
    INTEGER(I4b) j
    DO j  = 1, a%size
       a%data(j) =  mult*a%data(j)
    ENDDO
    RETURN
  END SUBROUTINE  Vector_mult

END MODULE Vector_class

