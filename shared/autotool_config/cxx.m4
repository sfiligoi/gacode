dnl ######################################################################
dnl
dnl Determine C++ compiler characteristics
dnl
dnl ######################################################################

dnl ######################################################################
dnl
dnl Check whether c++ compiler supports exception handling
dnl
dnl ######################################################################

AC_CACHE_CHECK(whether c++ compiler supports exception handling,
ac_cv_exception_handling,
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE_GLOBAL(
,
void f() { throw "abc"; }
void g() { try { f(); } catch(char*) { } }
,
ac_cv_exception_handling=yes, ac_cv_exception_handling=no)
AC_LANG_RESTORE)

if test $ac_cv_exception_handling = no
then
    AC_MSG_WARN(Compiler does not support exception handling!)
fi

dnl ######################################################################
dnl
dnl Check whether c++ compiler supports typename
dnl
dnl ######################################################################

AC_CACHE_CHECK(whether c++ compiler supports typename,
ac_cv_typename,
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE_GLOBAL(
,
template<class T> class X { typedef typename T::Y Y; };
class Z { public: typedef int Y; };
X<Z> x;
,
ac_cv_typename=yes, ac_cv_typename=no)
AC_LANG_RESTORE)

if test $ac_cv_typename = no
then
    AC_MSG_WARN(Compiler does not support typename!)
fi

dnl ######################################################################
dnl
dnl Check whether c++ compiler recognizes size_t is the same as uint
dnl
dnl ######################################################################

AC_CACHE_CHECK(whether c++ compiler recognizes size_t is the same as uint,
ac_cv_uintissize_t,
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE_GLOBAL(
,
#include <cstdio>
struct Holder {
  void operator<<(size_t i) { }
  void operator<<(unsigned int i) { }
};
,
ac_cv_uintissize_t=no, ac_cv_uintissize_t=yes)
AC_LANG_RESTORE)

if test $ac_cv_uintissize_t = yes; then
  AC_MSG_WARN(Compiler recognizes size_t is the same as uint)
  AC_DEFINE(UINT_IS_SIZE_T, , Define if C++ compiler recognizes uint as size_t)
fi

dnl ######################################################################
dnl
dnl Check whether c++ compiler recognizes ssize_t is the same as int
dnl
dnl ######################################################################

AC_CACHE_CHECK(whether c++ compiler recognizes ssize_t is the same as int,
ac_cv_intisssize_t,
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE_GLOBAL(
,
#include <iostream>
struct Holder {
  void operator<<(ssize_t i) { }
  void operator<<(int i) { }
};
,
ac_cv_intisssize_t=no, ac_cv_intisssize_t=yes)
AC_LANG_RESTORE)

if test $ac_cv_intisssize_t = yes; then
  AC_MSG_WARN(Compiler recognizes ssize_t is the same as int)
  AC_DEFINE(INT_IS_SSIZE_T, , Define if C++ compiler recognizes int as ssize_t)
fi

dnl ######################################################################
dnl
dnl Check whether c++ compiler can explicitly instantiate templates
dnl
dnl ######################################################################

AC_CACHE_CHECK(whether c++ compiler can explicitly instantiate templates,
ac_cv_explicit_templates,
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE_GLOBAL(
,
template<class T> class X { T t; }; template class X<int>; X<int> x;
,
ac_cv_explicit_templates=yes, ac_cv_explicit_templates=no)
AC_LANG_RESTORE)

if test $ac_cv_explicit_templates = no
then
    AC_MSG_WARN(Compiler does not support explicit instantiation of templates!)
fi

dnl ######################################################################
dnl
dnl Check whether c++ compiler supports RTTI
dnl
dnl ######################################################################

AC_CACHE_CHECK(whether c++ compiler supports RTTI,
ac_cv_rtti,
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE_GLOBAL(
,
struct X { virtual ~X() { } };
struct Y : virtual public X { virtual ~Y() { } };
void f() { Y y; X* xp = &y; Y* yp = dynamic_cast<Y*>(xp); }
,
ac_cv_rtti=yes, ac_cv_rtti=no)
AC_LANG_RESTORE)

if test $ac_cv_rtti = no
then
    AC_DEFINE(HAVE_NO_RTTI, , Define if C++ compiler does not support RTTI)
fi

dnl ######################################################################
dnl
dnl Check whether c++ compiler supports namespaces
dnl
dnl ######################################################################

AC_CACHE_CHECK(whether c++ compiler supports namespaces,
ac_cv_namespaces,
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE_GLOBAL(
,
namespace N { typedef int I; } N::I i;
,
ac_cv_namespaces=yes, ac_cv_namespaces=no)
AC_LANG_RESTORE)

if test $ac_cv_namespaces = yes
then
    AC_DEFINE(HAVE_NAMESPACES, , Define if C++ compiler supports namespaces)
fi

dnl ######################################################################
dnl
dnl Check whether c++ compiler has complex in the namespace std
dnl
dnl ######################################################################

if test $ac_cv_namespaces = yes
then
    AC_CACHE_CHECK(whether c++ compiler has complex in the namespace std,
    ac_cv_std_complex,
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
    AC_TRY_COMPILE_GLOBAL(
    ,
    #include <complex>
    std::complex<double> cd;
    ,
    ac_cv_std_complex=yes, ac_cv_std_complex=no)
    AC_LANG_RESTORE)

    if test $ac_cv_std_complex = yes
    then
	AC_DEFINE(HAVE_STD_COMPLEX, 1, Define if your C++ compiler has complex in the namespace std)
    fi
fi

dnl ######################################################################
dnl
dnl Check whether c++ compiler has streams in the namespace std
dnl
dnl ######################################################################

if test $ac_cv_namespaces = yes
then
    AC_CACHE_CHECK(whether c++ compiler has streams in the namespace std,
    ac_cv_std_streams,
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
    # if test $ac_cv_header_iostream = yes
    # then
    AC_TRY_COMPILE_GLOBAL(
    ,
    #include <iostream>
    using namespace std;
    ,
    ac_cv_std_streams=yes, ac_cv_std_streams=no)
    # else
    # AC_TRY_COMPILE_GLOBAL(
    # ,
    # #include <iostream.h>
    # using namespace std;
    # ,
    # ac_cv_std_streams=yes, ac_cv_std_streams=no)
    # fi
    AC_LANG_RESTORE)

    if test $ac_cv_std_streams = yes
    then
	AC_DEFINE(HAVE_STD_STREAMS, 1, Define if your C++ compiler has streams in the namespace std )
    fi
fi

dnl ######################################################################
dnl
dnl Check whether c++ compiler can overload const type conversions
dnl
dnl ######################################################################

AC_CACHE_CHECK(whether c++ compiler can overload const type conversions,
ac_cv_const_type_conversion_overload,
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE_GLOBAL(
,
struct X { int* p_; X(int* p) : p_(p) { }
operator int&() { return *p_; } operator const int&() const { return *p_; }
operator int*() { return p_; } operator const int*() const { return p_; } };
void f() { int i = 1; int j = 2;
X x(&i); const X cx(&j); int k = x; int l = cx; k = *x; l = *cx; }
,
ac_cv_const_type_conversion_overload=yes,
ac_cv_const_type_conversion_overload=no)
if test $ac_cv_const_type_conversion_overload = yes
then
    AC_TRY_COMPILE_GLOBAL(
    ,
    struct X { int* p_; X(int* p) : p_(p) { }
    operator int&() { return *p_; } operator int() const { return *p_; } };
    void f() { int i = 1; int j = 2;
    X x(&i); const X cx(&j); int k = x; int l = cx; }
    ,
    ac_cv_const_type_conversion_overload=yes,
    ac_cv_const_type_conversion_overload=no)
fi
AC_LANG_RESTORE)

if test $ac_cv_const_type_conversion_overload = no
then
    AC_DEFINE(HAVE_NO_CONST_TYPE_CONVERSION_OVERLOAD, ,
	Define if C++ compiler cannot overload const type conversions)
fi

dnl ######################################################################
dnl
dnl Check whether c++ compiler knows mutable
dnl
dnl ######################################################################

AC_CACHE_CHECK(whether c++ compiler knows mutable,
ac_cv_mutable,
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE_GLOBAL(
,
class X { mutable int i; }; X x;
,
ac_cv_mutable=yes, ac_cv_mutable=no)
AC_LANG_RESTORE)

if test $ac_cv_mutable = no
then
    AC_DEFINE(HAVE_NO_MUTABLE, , Define if your C++ compiler doesn't know mutable.)
fi

dnl ######################################################################
dnl
dnl Check whether template friends need brackets
dnl
dnl ######################################################################

AC_CACHE_CHECK(whether template friends need brackets,
ac_cv_template_friends_need_brackets,
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE_GLOBAL(
,
template <class T2> int foo(T2 t2);
template <class T1> class X { T1 t; friend int foo<> (T1); };
X<double> x;
,
ac_cv_template_friends_need_brackets=yes,
ac_cv_template_friends_need_brackets=no)
AC_LANG_RESTORE)

if test $ac_cv_template_friends_need_brackets = yes
then
    AC_DEFINE(TEMPLATE_FRIENDS_NEED_BRACKETS, 1,
	Define if your C++ compiler requires brackets in declarations of template friends.)
fi

dnl ######################################################################
dnl
dnl Check whether template friend operators allowed.  Works only for ansi
dnl declaration with <>.
dnl
dnl ######################################################################

AC_CACHE_CHECK(whether nontype template operators are allowed,
ac_cv_nontype_template_operators,
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE_GLOBAL(
,
template <int d> class X {
  friend void operator<< <> (int, X<d>&);
};
template <int d> void operator<< (int i, X<d>& xx){}
X<1> x;
,
ac_cv_nontype_template_operators=yes,
ac_cv_nontype_template_operators=no)
AC_LANG_RESTORE)

if test $ac_cv_nontype_template_operators = yes
then
    AC_DEFINE(HAVE_NONTYPE_TEMPLATE_OPERATORS, 1,
	Define if your C++ compiler can handle nontyped templated operators)
fi

dnl ######################################################################
dnl
dnl Whether the c++ compiler knows std::abs(double)
dnl
dnl ######################################################################

AC_CACHE_CHECK(whether the c++ compiler knows std::abs(double),
ac_cv_stdabsdbl,
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE_GLOBAL(
,
#include <cmath>
double a = 0;
double b = std::abs(a);
,
ac_cv_stdabsdbl=yes, ac_cv_stdabsdbl=no)
AC_LANG_RESTORE)
if test $ac_cv_stdabsdbl = no; then
  AC_DEFINE(NOT_HAVE_STD_ABS_DOUBLE, 1, Define if your C++ compiler does not know std::abs for double arg)
fi
AM_CONDITIONAL(NOT_HAVE_STD_ABS_DOUBLE, test $ac_cv_stdabsdbl = no)

dnl ######################################################################
dnl
dnl Check whether static variables can be declared generally
dnl
dnl ######################################################################

AC_CACHE_CHECK(whether static variables can be declared generally,
ac_cv_generally_declared_statics,
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_LINK_IFELSE([AC_LANG_SOURCE([
AC_LANG_PROGRAM(
template <class TYPE> class X {
  public:
    static int r;
};
template <class TYPE> int X<TYPE>::r = 0;
,
X<double> x;
int rr = x.r + X<float>::r;
)])],
ac_cv_generally_declared_statics=yes,
ac_cv_generally_declared_statics=no)
AC_LANG_RESTORE)

case $SERIALCXX in
  *-pgCC | pgCC)
# Detection does not work on pgCC
    :
    ;;
  *)
    if test $ac_cv_generally_declared_statics  = yes; then
      AC_DEFINE(HAVE_GENERALLY_DECLARED_STATICS, 1,
	[Define if your C++ compiler allows static variables to be declared generally])
    fi
    ;;
esac

dnl ######################################################################
dnl
dnl Check whether compiler handles templated member functions
dnl
dnl ######################################################################

AC_CACHE_CHECK(whether compiler handles templated member functions,
ac_cv_have_templated_member_functions,
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_LINK_IFELSE([AC_LANG_SOURCE([
AC_LANG_PROGRAM(
class foobar {
public:
  template <typename T>
  void getAttribNames() {
    T bar;
    return ;
  }
};
template <typename T>
void foo() {
  foobar ted;
  ted.getAttribNames< T >();
}
,
  foo<int>();
  return 0;
)])],
ac_cv_have_templated_member_functions=yes,
ac_cv_have_templated_member_functions=no)
AC_LANG_RESTORE)

if test $ac_cv_have_templated_member_functions = yes; then
    AC_DEFINE(HAVE_TEMPLATED_MEMBER_FUNCTIONS, 1,
	[Define if your C++ compiler can compile templated member functions])
fi



