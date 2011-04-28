dnl ######################################################################
dnl
dnl File:	tx_latex.m4
dnl
dnl Purpose:	Determine the location of pdflatex and associated
dnl
dnl Version: $Id: tx_latex.m4 3184 2009-08-09 17:02:35Z cary $
dnl
dnl ######################################################################

AC_PATH_PROGS(PDFLATEX, pdflatex, "", $PATH)
AM_CONDITIONAL(HAVE_PDFLATEX, test -n "$PDFLATEX")
if test -z "$PDFLATEX"; then
  AC_MSG_WARN(pdflatex was not found in your path.  Manuals will not be made.)
else
  AC_PATH_PROGS(MAKEINDEX, makeindex, "", $PATH)
  AC_PATH_PROGS(BIBTEX, bibtex, "", $PATH)
fi
AM_CONDITIONAL(HAVE_MAKEINDEX, test -n "$MAKEINDEX")
AM_CONDITIONAL(HAVE_BIBTEX, test -n "$BIBTEX")
AC_SUBST(PDFLATEX)
AC_SUBST(MAKEINDEX)
AC_SUBST(BIBTEX)
# echo "  PDFLATEX:  $PDFLATEX"   >> $config_summary_file
# echo "  MAKEINDEX: $MAKEINDEX"  >> $config_summary_file
# echo "  BIBTEX:    $BIBTEX"     >> $config_summary_file
TX_PRINT_VAR(PDFLATEX)
TX_PRINT_VAR(MAKEINDEX)
TX_PRINT_VAR(BIBTEX)

