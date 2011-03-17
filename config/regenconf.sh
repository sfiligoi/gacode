#!/bin/sh
######################################################################
#
# File:         regenconf.sh
#
# Purpose:      to regenerate autotools files
#
# Version:      $Id: regenconf.sh 3562 2010-06-07 01:53:37Z cary $
#
# Copyright 2003-2010, Tech-X Corporation
#
######################################################################

# Get rid of old files
echo Remaking all configuration files
mydir=`dirname $0`

# See if libtoolize should be run
if test -f configure.in; then
  res=`grep LIBTOOL configure.in | grep -v dnl`
else
  res=`grep LIBTOOL configure.ac | grep -v dnl`
fi
if test -n "$res"; then
  echo Redoing libtoolize
  ltcmd="libtoolize --force --copy"
  if test "`uname`" == "Darwin"; then

# NDS: glibtoolize works on OS X 10.4 (Tiger)
# JRC: glibtoolize does NOT work on my 10.4.
# JRC: Should look for (g)libtoolize in this order:
#
# 1. .../libtoolize (i.e., NOT /usr/bin)
# 2. glibtoolize
# 3. /usr/bin/libtoolize

    lt=/usr/bin/libtoolize
    tmp=`which glibtoolize | sed 's/ /_/g'`
    if test -x "$tmp"; then
      lt=$tmp
    fi
    tmp=`which libtoolize | sed 's/ /_/g'`
    if test -x "$tmp" -a "$tmp" != /usr/bin/libtoolize; then
      lt=$tmp
    fi
    ltcmd="$lt --force --copy"
  fi
else
  ltcmd=":"
fi

# Creating new files for this platform
for i in \
    "rm -rf aclocal.m4 autom4te.cache autom4te-*.cache configure" \
    "$ltcmd" \
    "env AUTOMAKE=config/automake.sh autoreconf" \
    aclocal \
    autoheader \
    "autoconf --force"; do
  echo $i
  if test "$i" == autoheader; then
    if test -e config.h; then
      $i
    fi
  else
  $i
  fi
  res=$?
  if test $res != 0; then
    echo Failure of $i.
    exit $res
  fi
done

# Do automake
cmd="config/automake.sh"
echo $cmd
$cmd
res=$?
if test $res != 0; then
  echo config/automake.sh failed.
  exit $res
fi

