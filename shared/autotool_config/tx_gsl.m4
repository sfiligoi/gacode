dnl ######################################################################
dnl
dnl File:	    tx_gsl.m4
dnl
dnl Purpose:	Determine where the GSL files are.
dnl
dnl Version:	$Id: tx_gsl.m4 3645 2010-09-02 21:29:15Z paulm $
dnl
dnl Copyright 2001-2010, Tech-X Corporation.  Redistribution allowed provided
dnl this copyright statement remains intact.
dnl
dnl ######################################################################

builtin(include, config/txsearch.m4)

dnl Locate all the babel related includes and libraries
TX_LOCATE_PKG(
  [gsl],
  [$HOME/software/gsl:/contrib/gsl:/usr/local:/contrib/gsl-1.11],
  [gsl],
  [gsl,gslcblas],
  [],
  [])


