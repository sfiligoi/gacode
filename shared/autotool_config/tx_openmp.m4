#
# SYNOPSIS
#
#   TX_OPENMP
#
# DESCRIPTION
#
#  Check whether compiler supports openMP directives
#
#  $Id: tx_openmp.m4 3596 2010-07-14 16:44:21Z pletzer $
#
AC_DEFUN([TX_C_OPENMP],[
AC_CACHE_CHECK([whether C compiler supports OpenMP directives],
ac_cv_have_c_openmp,
[AC_LANG_PUSH(C)
AC_OPENMP
CFLAGS="$OPENMP_CFLAGS $CFLAGS"
CPPFLAGS="$OPENMP_CFLAGS $CPPFLAGS"
AC_LANG_POP(C)
])
if test x"$ac_cv_have_c_openmp" = xyes; then
   AC_DEFINE(HAVE_C_OPENMP, [],[Defined if C compiler supports openMP directives])
   echo "C compiler supports openMP directives"  >> $config_summary_file
fi
])

AC_DEFUN([TX_CXX_OPENMP],[
AC_CACHE_CHECK([whether C++ compiler supports OpenMP directives],
ac_cv_have_cxx_openmp,
[AC_LANG_PUSH(C++)
AC_OPENMP
CXXFLAGS="$OPENMP_CXXFLAGS $CXXFLAGS"
CPPFLAGS="$OPENMP_CXXFLAGS $CPPFLAGS"
AC_LANG_POP(C++)
])
if test x"$ac_cv_have_cxx_openmp" = xyes; then
   AC_DEFINE(HAVE_CXX_OPENMP, [],[Defined if C++ compiler supports openMP directives])
   echo "C++ compiler supports openMP directives"  >> $config_summary_file
fi
])

AC_DEFUN([TX_FC_OPENMP],[
AC_CACHE_CHECK([whether Fortran compiler supports OpenMP directives],
ac_cv_have_fc_openmp,
[AC_LANG_PUSH(Fortran)
AC_OPENMP
FCFLAGS="$OPENMP_FCFLAGS $FCFLAGS"
AC_LANG_POP(Fortran)
])
if test x"$ac_cv_have_fc_openmp" = xyes; then
   AC_DEFINE(HAVE_FC_OPENMP, [],[Defined if Fortran compiler supports openMP directives])
   echo "Fortran compiler supports openMP directives"  >> $config_summary_file
fi
])

AC_DEFUN([TX_F77_OPENMP],[
AC_CACHE_CHECK([whether F77 compiler supports OpenMP directives],
ac_cv_have_f77_openmp,
[AC_LANG_PUSH(F77)
AC_OPENMP
FFLAGS="$OPENMP_FFLAGS $FFLAGS"
AC_LANG_POP(F77)
])
if test x"$ac_cv_have_f77_openmp" = xyes; then
   AC_DEFINE(HAVE_F77_OPENMP, [],[Defined if F77 compiler supports openMP directives])
   echo "F77 compiler supports openMP directives"  >> $config_summary_file
fi
])


