dnl
dnl @synopsis LLNL_ISO_C_BINDING
dnl
dnl
dnl @author muszala
dnl
dnl Note:  Define a configuraion macro to see if bindC is requested and then
dnl if the Fortran compiler actually has the iso_c_binding module.
  AC_MSG_CHECKING([if BindC is requested])
    ac_arg_with_bindc=no
    AC_ARG_WITH([bindc],
      AS_HELP_STRING(--with-bindc@<:@=bindc@:>@,F2003 iso_c_binding support @<:@default=no@:>@),
      [ case $withval in
          no) ac_arg_with_bindc=no ;; 
          yes)
           ac_arg_with_bindc=yes 
          ;; 
          *) ac_arg_with_bindc=no;
        esac]) # end AC_ARG_WITH
  AC_MSG_RESULT([$ac_arg_with_bindc])
  if test $ac_arg_with_bindc = yes; then
    AC_CACHE_CHECK([if this Fortran compiler actually has iso_c_binding],
    ac_cv_enable_bindc,
    [AC_LANG_PUSH(Fortran)
    AC_COMPILE_IFELSE([
      program main  
          use iso_c_binding
      end
      ],
      [ac_cv_enable_bindc="yes"],
      [ac_cv_enable_bindc="no"])
    AC_LANG_POP(Fortran)
    ])
    if test "$ac_cv_enable_bindc" != no; then
      AC_DEFINE(SIDL_HAS_ISO_C_BINDING,1,
	    [Define to 1 if the Fortran compiler supports iso_c_binding.])
      msgs="$msgs
          BindC requested and available."
    else 
      msgs="$msgs
          BindC request disabled against user request--iso_c_binding not supported."
    fi
  fi
  AM_CONDITIONAL(HAS_BINDC,test x$ac_cv_enable_bindc = xyes)
