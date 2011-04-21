dnl ----------------------------------------------------------------------
dnl
dnl Find rst2html.py for documentation
dnl
dnl ----------------------------------------------------------------------
AC_PATH_PROGS(RST2HTML, rst2html.py, "")
if test -z "$RST2HTML"; then
  AC_MSG_WARN(The rst2html.py utility was not found in your path.  Documentation will not be made.)
  enable_rst2html=true 
  if test -n "$config_summary_file"; then
    echo "Will not make html documentation"         >> $config_summary_file
    echo                                            >> $config_summary_file
  fi
else
  enable_rst2html=false 
  if test -n "$config_summary_file"; then
    echo "Will make documentation"                  >> $config_summary_file
    echo                                            >> $config_summary_file
  fi
fi
AC_SUBST(RST2HTML)
AM_CONDITIONAL(enable_rst2html, test x$enable_rst2html = xtrue)
