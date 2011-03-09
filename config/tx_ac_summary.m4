# ######################################################################
#
# File:		tx_ac_summary.m4
#
# Purpose:	To initialize the configuaration summary file
#
# Version:	$Id: tx_ac_summary.m4 3366 2010-01-15 18:43:10Z dws $
#
# Copyright 2007-2010, Tech-X Corporation. Redistribution allowed provided
# this copyright statement remains intact.
#
# ######################################################################

dnl ######################################################################
dnl
dnl Function:   TX_AC_SUMMARY_INIT(summary-file-name, calling-program)
dnl
dnl Purpose:    Initialize the configuration summary.
dnl
dnl ######################################################################

AC_DEFUN([TX_AC_SUMMARY_INIT], [

  config_summary_file=$1
  echo Configuration summary file is $config_summary_file
  top_config_summary_file=$abs_top_builddir/$config_summary_file
  AC_SUBST(top_config_summary_file)
  echo ""                                           >  $config_summary_file
  echo "Config line:"                               >> $config_summary_file
  echo $2 $ac_configure_args                        >> $config_summary_file
  AC_DEFINE_UNQUOTED(CONFIG_LINE, "$2 $ac_configure_args",
    [The command line used for configuration])
  echo ""                                           >> $config_summary_file
  echo "Directory: $PWD"                            >> $config_summary_file
  if test -z "$host"; then
    AC_MSG_WARN([You must define host before setting up config_summary_file])
  else
    echo "Host: $host"                                >> $config_summary_file
  fi
  config_date=`date`
  echo "Date: $config_date"                           >> $config_summary_file
  AC_DEFINE_UNQUOTED(CONFIG_DATE, "$config_date",
    [The date when configured])

])

dnl ######################################################################
dnl
dnl Function:   TX_AC_SUMMARY_FINALIZE
dnl
dnl Purpose:    Finalize the configuration summary.
dnl
dnl ######################################################################

AC_DEFUN([TX_AC_SUMMARY_FINALIZE], [

  if test -n "$config_summary_file"; then
    echo >>$config_summary_file
# The following line is needed to detect success
    echo "Successfully configured." >>$config_summary_file
    echo
    echo "-----------------------------------"
    echo "Configure summary (from file $config_summary_file):"
    cat $config_summary_file
    echo "-----------------------------------"
    echo
  else
    echo config_summary_file not defined.
  fi

])

