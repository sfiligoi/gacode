dnl ######################################################################
dnl
dnl File:       txsearch.m4
dnl
dnl Purpose:    Defines a series of macros to be used to search for common
dnl		files, directores, libraries, and includes.
dnl
dnl Functions defined:
dnl   TX_CLEAN_LIBS
dnl   TX_PRINT_VAR
dnl   TX_PATH_FILE
dnl   TX_PATH_FILES
dnl
dnl Version:    $Id: txsearch.m4 3695 2010-10-03 21:14:30Z cary $
dnl
dnl Copyright 2001-2010, Tech-X Corporation.  Redistribution allowed
dnl provided this copyright statement remains intact.
dnl
dnl ######################################################################

dnl ######################################################################
dnl
dnl Function: TX_CLEAN_LIBS(variable): Clean up the libs variable
dnl by removing system libdirs.
dnl for VARLIBS define
dnl   VARRPLIBS to contain same with RPATH <libdir> for each libdir
dnl   VARLTLIBS to contain same with -rpath <libdir> for each libdir
dnl   VARALIBS to contain collection of fully qualified static libs
dnl
dnl ######################################################################

AC_DEFUN([TX_CLEAN_LIBS], [
  unset tx_clean_libs
  unset tx_clean_libdirs
  unset tx_clean_other
  rm -f tx_clean_libdirs$$
  touch tx_clean_libdirs$$
  eval tx_clean_val="\$$1"
  # echo For $1: tx_clean_val = $tx_clean_val
  grabnext=false
  for i in $tx_clean_val; do
# get a second arg
    if $grabnext; then
      tx_clean_other="$tx_clean_other $i"
      grabnext=false
      continue
    fi
    case $i in
      -framework)
        tx_clean_other="$tx_clean_other $i"
        grabnext=true
        ;;
      -L*)
        tx_clean_libdir=`echo $i | sed 's/^-L//'`
# Clean up directory if possible, may not be possible due to cross compilation
        if test -d $tx_clean_libdir; then
          tx_clean_libdir=`(cd $tx_clean_libdir; pwd -P)`
        fi
        case $tx_clean_libdir in
          /lib | /lib64 | /usr/lib | /usr/lib64 | /usr/libexec | /usr/lib/gcc/* )
# Skip these, since system has them
            ;;
          *)
            hasit=`grep ^$tx_clean_libdir\$ tx_clean_libdirs$$`
            if test -z "$hasit"; then
              echo $tx_clean_libdir >>tx_clean_libdirs$$
              tx_clean_libdirs="$tx_clean_libdirs $tx_clean_libdir"
            fi
            ;;
        esac
        ;;
      -l*)
        tx_clean_lib=`echo $i | sed 's/^-l//'`
        tx_clean_libs="$tx_clean_libs $tx_clean_lib"
        ;;
      /*/lib*.a)
        tx_clean_lib=`echo $i | sed -e 's?^.*/lib??' -e 's?\.a$??'`
        tx_clean_libdir=`echo $i | sed -e "s?/lib$tx_clean_lib.a??"`
        tx_clean_libs="$tx_clean_libs $tx_clean_lib"
        tx_clean_libdirs="$tx_clean_libdirs $tx_clean_libdir"
        ;;
      /*/*.a)
# Not starting with lib
        tx_clean_libs="$tx_clean_libs $i"
        ;;
      *)
        tx_clean_other="$tx_clean_other $i"
        ;;
    esac
  done
  # echo For $1: tx_clean_other = $tx_clean_other
dnl make clean library var with libdirs
  unset $1
  for i in $tx_clean_libdirs; do
    eval $1=\"\$$1 -L$i\"
  done
dnl Determine whether any dynamic libs present
  hasdynlibs=false
  for i in $tx_clean_libs; do
    for j in $tx_clean_libdirs; do
      if test -f $j/lib${i}.so -o -f $j/lib${i}.dylib; then
        # echo Found shared lib, lib${i}.so, in $j
        hasdynlibs=true
        break
      else
        : # echo Did not find shared lib, lib${i}.so, in $j
      fi
    done
    if $hasdynlibs; then
      break
    fi
  done
dnl Above not accurate when no library in the
  dnl echo tx_clean_libdirs = $tx_clean_libdirs
  dnl echo tx_clean_libs = $tx_clean_libs
  dnl hasdynlibs=true
dnl Make auxiliary variables with rpath vars
  baselibsname=`echo $1 | sed 's/LIBS\$//'`
  rplibsname=${baselibsname}RPLIBS
  ltlibsname=${baselibsname}LTLIBS
  alibsname=${baselibsname}ALIBS
  unset $rplibsname
  unset $ltlibsname
  unset $alibsname
  if $hasdynlibs; then
    for i in $tx_clean_libdirs; do
      if test -n "$RPATH_FLAG"; then
        if test "$RPATH_FLAG" != "-L"; then
          eval $rplibsname=\"\$$rplibsname -L$i ${RPATH_FLAG}$i\"
        else
          eval $rplibsname=\"\$$rplibsname -L$i\"
        fi
      else
        eval $rplibsname=\"\$$rplibsname -L$i \\$\(RPATH_FLAG\)$i\"
      fi
      eval $ltlibsname=\"\$$ltlibsname -L$i -rpath $i\"
    done
  else
    eval $rplibsname=\"\$$1\"
    eval $ltlibsname=\"\$$1\"
  fi
dnl Add in libraries
  for i in $tx_clean_libs; do
    eval $1=\"\$$1 -l$i\"
    eval $rplibsname=\"\$$rplibsname -l$i\"
    eval $ltlibsname=\"\$$ltlibsname -l$i\"
    unset alibname
    for j in $tx_clean_libdirs; do
      if test -f $j/lib${i}.a; then
        alibname=$j/lib${i}.a
        break
      fi
    done
dnl Add in .a file if found, otherwise just -l<libname>
    if test -n "$alibname"; then
      eval $alibsname=\"\$${alibsname} $alibname\"
    else
      eval $alibsname=\"\$${alibsname} -l$i\"
    fi
  done
  rm -f tx_clean_libdirs$$
dnl Add in other
  eval $1=\"\$$1 $tx_clean_other\"
  eval $rplibsname=\"\$$rplibsname $tx_clean_other\"
  eval $ltlibsname=\"\$$ltlibsname $tx_clean_other\"
  eval $alibsname=\"\$$alibsname $tx_clean_other\"
dnl clean leading and trailing spaces
  dnl for i in $1 $rplibsname $ltlibsname $alibsname; do
  for i in $1 $ltlibsname $alibsname; do
    eval tx_pkg_val=\"\$$i\"
    tx_pkg_val=`echo $tx_pkg_val | sed -e 's/^  *//' -e 's/  *$//'`
    eval $i=\"$tx_pkg_val\"
  done
dnl substs occur outside
  dnl AC_SUBST($rplibsname)
  dnl AC_SUBST($ltlibsname)
  dnl AC_SUBST($alibsname)
])

dnl ######################################################################
dnl
dnl Print an environment variable according to standard format
dnl
dnl ######################################################################

AC_DEFUN([TX_PRINT_VAR], [
  if test -n "$2"; then
    tx_print_val="$2"
  else
    eval tx_print_val="\$$1"
  fi
  printf "    %-27s: %s\n" $1 "$tx_print_val" >>$config_summary_file
])

dnl ######################################################################
dnl
dnl Function:	TX_PATH_FILE( variable , file-to-check-for ,
dnl	[value-if-not-found] , [path] )
dnl
dnl Purpose:	Like AC_PATH_PROG but searches for a file instead of a
dnl	program.  Calls AC_SUBST for "variable".
dnl
dnl ######################################################################

AC_DEFUN([TX_PATH_FILE],[

  if test -z "$5"; then
    AC_MSG_CHECKING([for $2])
  fi

  dnl Check for a path argument.  If none is given, set to current
  dnl directory, otherwise expand the path argument.
  if test -z "$4"; then
    tx_search_path="."
  else
    tx_search_path=`echo "$4" | tr ':' ' '`
  fi

  tx_return=""
  for tx_path in $tx_search_path; do
    ## if test -z `echo "$tx_path" | grep "/$"`; then
    if test -a "$tx_path"/$2; then
dnl JRC: This is not correct on octet.  Cannot trim like this.
dnl      tx_path=`(cd $tx_path; pwd -P | \
dnl			 sed -e 's@/usr/local/contrib@/contrib@' \
dnl			     -e 's@/usr/local/internal@/internal@')`
      tx_path=`(cd $tx_path; pwd -P)`
      tx_return="$tx_path"/$2
      break
    fi
  done

  if test -z "$tx_return"; then
    $1="$3"
    if test -z "$5"; then
      AC_MSG_RESULT([no])
    fi
  else
    $1="$tx_return"
    if test -z "$5"; then
      AC_MSG_RESULT([$tx_return])
    fi
  fi

  AC_SUBST([$1])

])

dnl ######################################################################
dnl
dnl Function:	TX_PATH_FILES( variable , files-to-check-for ,
dnl	[value-if-not-found] , [path] )
dnl
dnl Purpose:	Like AC_PATH_PROGS but searches for files instead of
dnl	programs.  Calls AC_SUBST for "variable".
dnl
dnl ######################################################################

AC_DEFUN([TX_PATH_FILES],[

  srchfile=`echo $2 | tr ' ' ','`
  AC_MSG_CHECKING([for $srchfile])

  for tx_dirs in `echo "$4" | tr ':' ' '`; do
    for tx_files in $2; do
      TX_PATH_FILE([$1],[$tx_files],[],[$tx_dirs],[noprint])
      if test -n "$$1"; then break; fi
    done
    if test -n "$$1"; then break; fi
  done

  if test -z "$$1"; then
    $1="$3"
    AC_MSG_RESULT([no])
  else
    AC_MSG_RESULT([$$1])
  fi

  AC_SUBST([$1])

])

dnl ######################################################################
dnl
dnl Function:	TX_LOCATE_PKG( pkgname , default-path , [inc-files] ,
dnl	[lib-files] , [default-inc-path] , [default-lib-path] )
dnl
dnl Purpose:	All inclusive search for package files in the case of a
dnl	simple package.
dnl
dnl Inputs:
dnl     pkgname:
dnl		    Name of the package to search for.
dnl     default-path:
dnl         Default path in which to search for the package,
dnl	        colon delimited for multiple possiblities.
dnl     inc-files: (optional)
dnl	        Comma separated list of critical includes required for
dnl	        the package to function.  Explicitly passing [-] to
dnl         this argument will disable include file search
dnl	        defaults to pkgname.h
dnl	    lib-files: (optional)
dnl	        Comma separated list of critical libs required for the
dnl	        package to function.  Explicitly passing [-] to
dnl         this argument will disable library file search
dnl	        defaults to libpkgname.so(a,dylib)
dnl	    default-inc-path: (optional)
dnl	        Location of the includes for this package,
dnl	        colon delimited for multiple possibilities.
dnl	        defaults to default-path/include/
dnl     default-lib-path: (optional)
dnl	        location of the libraries for this package,
dnl	        colon delimited for multiple possibilities.
dnl	        defaults to default-path/lib
dnl
dnl Variables read during input:
dnl     pkgname_DOSEARCH:   if set to 'no', then the package finding
dnl                         system will automatically disable the
dnl                         searching of this package.  This is
dnl                         overridden by command line input.
dnl
dnl Outputs (via AC_SUBST):
dnl     pkgname_INCDIR:     directory of the last include
dnl     pkgname_INC:        complete, flagged, includes variable
dnl     pkgname_INC_file:   absolute path to a specific header
dnl                         (one for each file)
dnl     pkgname_LIBDIR:     directory of the first library
dnl     pkgname_LIBS:       complete, flagged, libraries variable
dnl     pkgname_LIB_lib:    absolute path to a specific library
dnl                         (one for each library)
dnl
dnl Outputs (via AC_DEFINE):
dnl	    HAVE_pkgname:       defined if includes and libs are found
dnl
dnl Outputs (via AM_CONDITIONAL):
dnl     HAVE_pkgname:       defined if includes and libs are found
dnl
dnl Variables set for later use:
dnl     pkgname_DOSEARCH:   set to 'no' if the package was disabled by
dnl                         user request
dnl     FOUND_pkgname:      set to "yes" if package was found, or "no"
dnl                         otherwise
dnl
dnl Also defines these user override flags
dnl     --disable-pkgname:          Disables the package
dnl     --with-pkgname-dir:         Prepends default
dnl     --with-pkgname-incdir:      Prepends default
dnl     --with-pkgname-inc-file:    Overrides default (one for each file)
dnl     --with-pkgname-libdir:      Prepends default
dnl     --with-pkgname-lib-lib:     Overrides default (one for each file)
dnl
dnl ######################################################################

dnl ----------------------------------------------------------------------
dnl Definition of TX_LOCATE_PKG
dnl ----------------------------------------------------------------------

AC_DEFUN([TX_LOCATE_PKG],[

  tx_failure="no"

  m4_ifval([$1],[],[AC_DIAGNOSE([syntax],[No package name specified])])

  dnl Create the package's unique variable prepender
  PKGNM=`echo $1 | tr '[a-z./-]' '[A-Z___]'`

dnl Set up the disable flag to allow for the disabling of certain packages
  eval tmpvar=$`echo "$PKGNM"_DOSEARCH`
  dnl if test -n "$tmpvar"; then
  if test -z "$tmpvar"; then
    AC_ARG_ENABLE(
      translit($1,'A-Z./','a-z__'),
      [AC_HELP_STRING([--disable-]translit($1,'A-Z./','a-z__'), [Disables the use of the $1 package])],
      [eval "$PKGNM"_DOSEARCH=$enableval])
  fi

dnl --------------------------------------------------------------------
dnl Set the default package installation directory
dnl --------------------------------------------------------------------

  m4_ifval([$2],[],[AC_DIAGNOSE([syntax],[No default directory specified for $1 package])])

  dnl Additional user input
  AC_ARG_WITH(
    translit($1,'A-Z./','a-z__')[-dir],
    [AC_HELP_STRING([--with-]translit($1,'A-Z./','a-z__')[-dir=<dir>],[Installation directory of $1 package])],
    [eval "$PKGNM"_PATH="$withval:$2"],
    [eval "$PKGNM"_PATH="$2"])

  eval tmpvar=$`echo "$PKGNM"_DOSEARCH`
  if test x$tmpvar = xno; then
    pkgdir=""
  else
    eval pkgdir=$`echo "$PKGNM"_PATH`
  fi
  eval "$PKGNM"_PATH=\"\"
  for pdir in `echo $pkgdir | tr ':' ' '`; do
    pkgdirsubstr=`echo "$pdir" | grep "^/"`
    if test -z "$pkgdirsubstr"; then
dnl Relative path specified
      pdir=$PWD/$pdir
    fi
    dnl Trailing / missing
    pkgdirsubstr=`echo "$pdir" | grep "/$"`
    if test -z "$pkgdirsubstr"; then
dnl Relative path specified
      eval "$PKGNM"_PATH=$`echo "$PKGNM"_PATH`:$pdir/
    else
      eval "$PKGNM"_PATH=$`echo "$PKGNM"_PATH`:$pdir
    fi
  done
  eval pkgdir=$`echo "$PKGNM"_PATH`
  eval tmpvar=$`echo "$PKGNM"_DOSEARCH`
  if test -z "$pkgdir" -a x$tmpvar != xno; then
    AC_MSG_ERROR([Package $1 has no default directories specified, use --with-]translit($1,'A-Z./','a-z__')[-dir to specify one])
  fi

dnl Write to summary file if defined
  if test -n "$config_summary_file"; then
    eval tmpvar=$`echo "$PKGNM"_DOSEARCH`
    if test x$tmpvar = xno; then
      echo                                         >> $config_summary_file
      eval printf \"NOT using %s\\n\" \"$1\"       >> $config_summary_file
    else
      echo                                         >> $config_summary_file
      eval printf \"Using %s with\\n\" \"$1\"      >> $config_summary_file
      eval printf \""  "Directories searched:\\n\" >> $config_summary_file
      for pdir in `echo $pkgdir | tr ':' ' '`; do
        eval printf \""    "%s\\n\" \"$pdir\"	   >> $config_summary_file
      done
    fi
  fi

dnl --------------------------------------------------------------------
dnl Set the default include path search directories
dnl --------------------------------------------------------------------

  eval "$PKGNM"_INCDIR_PATH=\"\"
  tx_inclist_glob=""

dnl Additional user input
  AC_ARG_WITH(
    translit($1,'A-Z./','a-z__')[-incdir],
    [AC_HELP_STRING([--with-]translit($1,'A-Z./','a-z__')[-incdir=<dir>],[Includes directory of $1 package])],
    [tempvar=$withval],
    [unset tempvar])
  tx_inclist_glob="$tempvar":"$tx_inclist_glob"

dnl Additional specified locations
  tx_inclist_glob="$tx_inclist_glob":"$5"

dnl Default search location
  tx_inclist_glob="$tx_inclist_glob":include/

dnl Compile the list
  tx_inclist=""
  for pdir in `echo $pkgdir | tr ':' ' '`; do
    for itr in `echo "$tx_inclist_glob" | tr ':' ' '`; do
      tx_isabs=`echo "$itr" | grep "^/"`
      if test -z "$tx_isabs"; then
dnl Relative path specified
        incpath=${pdir}${itr}
        hastrailing=`echo "$incpath" | grep "/$"`
        if test -z "$hastrailing"; then
          eval incpath="$incpath"/
        fi
        TX_LOCATE_PKG_UNIQUE_CAT([$incpath],[tx_inclist],[:])
      else
dnl Absolute path specified
        incpath=${itr}
        hastrailing=`echo "$incpath" | grep "/$"`
        if test -z "$hastrailing"; then
          eval incpath="$incpath"/
        fi
        TX_LOCATE_PKG_UNIQUE_CAT([$incpath],[tx_inclist],[:])
      fi
    done
  done

dnl Push the list into the proper container
  eval "$PKGNM"_INCDIR_PATH="$tx_inclist"

dnl --------------------------------------------------------------------
dnl Locate and set the actual requested include files
dnl --------------------------------------------------------------------

  eval tx_incpath=$`echo "$PKGNM"_INCDIR_PATH`
  # echo Looking for includes in $tx_incpath

dnl First off, check and see if we even need to look for headers
  tx_heads=""
  m4_ifval(
    [$3],
    [if test "$3" = "-"; then tx_heads="-- none --"; else tx_heads="$3"; fi],
    [tx_heads="translit($1,'A-Z./','a-z__').h"])

dnl Write to summary file if defined
  eval tmpvar=$`echo "$PKGNM"_DOSEARCH`
  if test -n "$config_summary_file" -a x$tmpvar != xno; then
    eval printf \""  "Includes sought: %s\\n\" \"$tx_heads\" >> $config_summary_file
  fi

dnl Skip the header check if we're not looking for any headers
  eval tmpvar=$`echo "$PKGNM"_DOSEARCH`
dnl jrc11oct08: but must still unset the suspect dir.
  tx_suspect_dir=""
  if test "$tx_heads" != "-- none --" -a x$tmpvar != xno; then

dnl Catch and save variables that are previously defined at the command prompt
    test_save=""; eval test_save=\"$`echo "$PKGNM"_INC`\"
    if test -n "$test_save"; then eval SAVE_"$PKGNM"_INC=\"$test_save\"; fi
    test_save=""; eval test_save=\"$`echo "$PKGNM"_INCDIR`\"
    if test -n "$test_save"; then eval SAVE_"$PKGNM"_INCDIR=\"$test_save\"; fi

dnl Begin the actual tests
    eval "$PKGNM"_INC=""

dnl Loop over all includes specified
    m4_ifval(
      [$3],
      [TX_LOCATE_PKG_HEADLOOP([$1],[$3])],
      [TX_LOCATE_PKG_HEADLOOP([$1],translit($1,'A-Z./','a-z__')[.h])])

dnl Push the data into the correct container.
    # echo tx_incs = $tx_incs
    for itr in $tx_incs; do
      eval "$PKGNM"_INC=\"$`echo "$PKGNM"_INC` -I$itr\"
    done
    # eval val="\$${PKGNM}_INC"
    # echo ${PKGNM}_INC = $val

dnl Replace variables with command prompt overrides if supplied.
dnl If this override is used, reset the failure flag
    eval test_save=\"$`echo SAVE_"$PKGNM"_INC`\"
dnl The modules system can leave a bogus value in here, so must
dnl start with -I to use
    case "$test_save" in
      -I*)
        eval "$PKGNM"_INC=\"$test_save\"
        tx_failure="no"
        ;;
    esac
    eval test_save=\"$`echo SAVE_"$PKGNM"_INCDIR`\"
    if test -n "$test_save"; then
      eval "$PKGNM"_INCDIR=\"$test_save\"
      tx_failure="no"
    fi

dnl Clean up and set for substitution
    eval tx_pkg_inc=\"$`echo "$PKGNM"_INC`\"
    tx_pkg_inc=`echo $tx_pkg_inc | sed -e 's/^  *//' -e 's/  *$//'`
    eval "$PKGNM"_INC=\"$tx_pkg_inc\"
    AC_SUBST(translit($1,'a-z./-','A-Z____')[_INC])
    AC_SUBST(translit($1,'a-z./-','A-Z____')[_INCDIR])

dnl Write to summary file if defined
    if test -n "$config_summary_file"; then
      TX_PRINT_VAR(${PKGNM}_INCDIR)
      TX_PRINT_VAR(${PKGNM}_INC)
    fi

  fi

dnl --------------------------------------------------------------------
dnl Set the default library path search directories
dnl --------------------------------------------------------------------

  eval "$PKGNM"_LIBDIR_PATH=\"\"
  tx_liblist_glob=""

dnl Additional user input
  AC_ARG_WITH(
    translit($1,'A-Z./','a-z__')[-libdir],
    [AC_HELP_STRING([--with-]translit($1,'A-Z./','a-z__')[-libdir=<dir>],[Libraries directory of $1 package])],
    [tempvar=$withval],
    [unset tempvar])
  tx_liblist_glob="$tempvar":"$tx_liblist_glob"

dnl Additional specified locations
  tx_liblist_glob="$tx_liblist_glob":"$6"

dnl Default search location
  tx_liblist_glob="$tx_liblist_glob":lib/

dnl Compile the list
  tx_liblist=""
  for pdir in `echo $pkgdir | tr ':' ' '`; do
    for itr in `echo "$tx_liblist_glob" | tr ':' ' '`; do
      tx_isabs=`echo "$itr" | grep "^/"`
      if test -z "$tx_isabs"; then
dnl Relative path specified
        libpath=${pdir}${itr}
        hastrailing=`echo "$libpath" | grep "/$"`
        if test -z "$hastrailing"; then
          eval libpath="$libpath"/
        fi
        TX_LOCATE_PKG_UNIQUE_CAT([$libpath],[tx_liblist],[:])
      else
dnl Absolute path specified
        libpath=${itr}
        hastrailing=`echo "$libpath" | grep "/$"`
        if test -z "$hastrailing"; then
          eval libpath="$libpath"/
        fi
        TX_LOCATE_PKG_UNIQUE_CAT([$libpath],[tx_liblist],[:])
      fi
    done
  done

dnl echo "*************************"
dnl echo "$tx_liblist" | tr ':' '\n'
dnl echo "*************************"

dnl Push the list into the proper container
  eval "$PKGNM"_LIBDIR_PATH="$tx_liblist"

dnl --------------------------------------------------------------------
dnl Locate and set the actual required library files
dnl --------------------------------------------------------------------

dnl First off, check and see if we even need to look for libraries
  tx_libs=""
  m4_ifval(
    [$4],
    [if test "$4" = "-"; then tx_libs="-- none --"; else tx_libs="$4"; fi],
    [tx_libs="translit($1,'A-Z./','a-z__')"])

dnl Write to summary file if defined
  eval tmpvar=$`echo "$PKGNM"_DOSEARCH`
  if test -n "$config_summary_file" -a x$tmpvar != xno; then
    eval printf \""  "Libraries sought: %s\\n\" \"$tx_libs\"		        				>> $config_summary_file
  fi

dnl Dont bother if already failed, skip the library check if we're not
dnl looking for any libraries
  dnl if test "$tx_failure" = "no" -a "$tx_libs" != "-- none --"; then
dnl JRC: continue looking for libs even if continue failed
  eval tmpvar=$`echo "$PKGNM"_DOSEARCH`
  if test "$tx_libs" != "-- none --" -a x$tmpvar != xno; then

dnl Catch and save variables that are previously defined at the command prompt
    test_save=""; eval test_save=\"$`echo "$PKGNM"_LIBS`\"
    if test -n "$test_save"; then eval SAVE_"$PKGNM"_LIBS=\"$test_save\"; fi
    test_save=""; eval test_save=\"$`echo "$PKGNM"_LIBDIR`\"
    if test -n "$test_save"; then eval SAVE_"$PKGNM"_LIBDIR=\"$test_save\"; fi

    dnl Do the actual search
    eval "$PKGNM"_LIBS=""

    tx_dynamic=no
    m4_ifval(
      [$4],
      [TX_LOCATE_PKG_LIBLOOP([$1],[$4])],
      [TX_LOCATE_PKG_LIBLOOP([$1],translit($1,'A-Z./','a-z__'))])

dnl Add in regular dirs and libraries
    for itr in $tx_lib_paths; do
      eval "$PKGNM"_LIBS=\"$`echo "$PKGNM"_LIBS` -L$itr\"
    done
    for itr in $tx_libs; do
      eval "$PKGNM"_LIBS=\"$`echo "$PKGNM"_LIBS` -l$itr\"
    done

dnl Replace variables with command prompt overrides if supplied
dnl If this override is used, reset the failure flag
    test_save=""
    eval test_save=\"$`echo SAVE_"$PKGNM"_LIBS`\"
    if test -n "$test_save"; then
      eval "$PKGNM"_LIBS=\"$test_save\"
      tx_failure="no"
    fi
    test_save=""
    eval test_save=\"$`echo SAVE_"$PKGNM"_LIBDIR`\"
    if test -n "$test_save"; then
      eval "$PKGNM"_LIBDIR=\"$test_save\"
      tx_failure="no"
    fi

dnl Export
    TX_CLEAN_LIBS([translit($1,'a-z./-','A-Z____')[_LIBS]])
    AC_SUBST(translit($1,'a-z./-','A-Z____')[_LIBS])
    AC_SUBST(translit($1,'a-z./-','A-Z____')[_RPLIBS])
    AC_SUBST(translit($1,'a-z./-','A-Z____')[_LTLIBS])
    AC_SUBST(translit($1,'a-z./-','A-Z____')[_ALIBS])
    AC_SUBST(translit($1,'a-z./-','A-Z____')[_LIBDIR])

dnl Write to summary file if defined
    if test -n "$config_summary_file"; then
      TX_PRINT_VAR(${PKGNM}_LIBDIR)
      TX_PRINT_VAR(${PKGNM}_LIBS)
      TX_PRINT_VAR(${PKGNM}_RPLIBS)
      TX_PRINT_VAR(${PKGNM}_LTLIBS)
      TX_PRINT_VAR(${PKGNM}_ALIBS)
    fi

  fi

dnl --------------------------------------------------------------------
dnl Confirm that the package was found
dnl --------------------------------------------------------------------
  eval tmpvar=$`echo "$PKGNM"_DOSEARCH`
dnl  echo ************"$PKGNM"_DOSEARCH="$tmpvar" 
  if test x$tx_failure = xno -a x$tmpvar != xno; then
    AC_DEFINE([HAVE_]translit($1,'a-z./-','A-Z____'),[],[Defined if $1 found])
    eval FOUND_"$PKGNM"="yes"
  else
    eval FOUND_"$PKGNM"="no"
  fi
  AM_CONDITIONAL([HAVE_]translit($1,'a-z./-','A-Z____'),[test x$tx_failure = xno -a x$tmpvar != xno])

])

dnl ----------------------------------------------------------------------
dnl Helper functions TX_LOCATE_PKG_FILELOOP used by TX_LOCATE_PKG
dnl ----------------------------------------------------------------------

dnl This is a generic function used by TX_LOCATE_PKG that on the m4 level
dnl loops through the different possible header files, adds and sets
dnl options for each file, and returns outputs accordingly
AC_DEFUN([TX_LOCATE_PKG_HEADLOOP],[

  dnl Initialize local container variables
  tx_incs=""

  m4_foreach([incfile],m4_dquote($2),[
    dnl Create this file's unique variable postpender
    FILENM=`echo "incfile" | tr '[a-z./-]' '[A-Z___]'`

    dnl Additional user input
    AC_ARG_WITH(
      translit($1,'A-Z./','a-z__')[-inc-]translit(incfile,'A-Z./','a-z__'),
      [  --with-]translit($1,'A-Z./','a-z__')[-inc-]translit(incfile,'A-Z./','a-z__')[=<file>
                          Full path location of the ]incfile[ header file],
      [eval searchdir=`dirname $withval`],
      [eval searchdir=$`echo "$PKGNM"_INCDIR_PATH`])

dnl Find the file
    echo Looking for includes in $searchdir
    TX_PATH_FILE(translit($1,'a-z./-','A-Z____')[_INC_]translit(incfile,'a-z./-','A-Z____'),incfile,[],[$searchdir])
    eval tempvar="$`echo "$PKGNM"_INC_"$FILENM"`"

    tmpfile=incfile
    while true; do
      tmpfile=`dirname tmpfile`
      if test "$tmpfile" = "."; then
        break
      fi
      tempvar=`dirname tempvar`
    done

    dnl Confirm the integrity of the file
    test_integ=""
    if test -n "$tempvar"; then
      for tx_itr in `echo $tx_inclist_glob | tr ':' ' '`; do
        abspath=`echo $tx_itr | grep '^/'`
        if test -n "$abspath"; then
          test_integ=""
          break
        fi
        TX_LOCATE_PKG_EXTRACT_ROOT([`dirname $tempvar`],[$tx_itr])
        if test -n "$return"; then
          test_integ=$test_integ" "$return
        fi
      done
    fi
    dnl We're going to iteratively grep each string in test_integ, and then pick the one that is the
    dnl greatest common factor
    cand=""
    if test -n "$test_integ"; then
      for tx_iitr in $test_integ; do
        loopcontrol=true
        for tx_jitr in $test_integ; do
          tmp=`echo $tx_jitr | grep $tx_iitr`
          if test -z "$tmp"; then
            loopcontrol=false
            break
          fi
        done
        if $loopcontrol; then
          cand="$tx_iitr"
          break
        fi
      done
    fi
dnl Now that we have a candidate, we're going to compare it to our
dnl already acquired tx_suspect_dir.  If it does not match tx_suspect_dir,
dnl we will zero out $tempvar because it's from a different
dnl build installation.  If tx_suspect_dir is empty, we'll fill it with
dnl the candidate.  If the candidate is empty, we're going to let it
dnl pass (pathological case).
    if test -z "$cand"; then
      tempvar="$tempvar"
    elif test -n "$cand" -a -z "$tx_suspect_dir"; then
      tx_suspect_dir=$cand
    elif test x"$cand" != x"$tx_suspect_dir"; then
      tempvar="-- conflict --"
    fi

dnl If the file is found, append the flag to tx_incs and write
dnl the directory to _INCDIR, otherwise warn
    if test -n "$tempvar" -a x"$tempvar" != x"-- conflict --"; then
      TX_LOCATE_PKG_UNIQUE_CAT([`dirname $tempvar`],[tx_incs])
      echo Adding `dirname $tempvar` to tx_incs.
      eval "$PKGNM"_INCDIR=`dirname $tempvar`
    elif test x"$tempvar" = x"-- conflict --"; then
      AC_MSG_WARN([Conflicting ]incfile[ header file found in $searchdir])
      tx_failure="yes"
    else
      AC_MSG_WARN([Could not find ]incfile[ header file in $searchdir])
      tx_failure="yes"
    fi

dnl Write to summary file if defined
    if test -n "$config_summary_file"; then
      if test -n "$tempvar" -a x"$tempvar" != x"-- conflict --"; then
        TX_PRINT_VAR(${PKGNM}_INC_$FILENM)
      elif test x"$tempvar" = x"-- conflict --"; then
        TX_PRINT_VAR(${PKGNM}_INC_$FILENM, [-- conflict --])
      else
        TX_PRINT_VAR(${PKGNM}_INC_$FILENM, [-- not found --])
      fi
    fi
  ])

])

dnl This is a generic function used by TX_LOCATE_PKG that on the m4 level
dnl loops through the different possible library files, adds and sets
dnl options for each file, and returns outputs accordingly
AC_DEFUN([TX_LOCATE_PKG_LIBLOOP],[

  dnl Initialize local container variables
  tx_libs=""
  tx_lib_paths=""

  m4_foreach([libfile],m4_dquote($2),[
    dnl Create this file's unique variable postpender
    FILENM=`echo "libfile" | tr '[a-z./-]' '[A-Z___]'`

    searchfile=""
    dnl Additional user input
    AC_ARG_WITH(
      translit($1,'A-Z./','a-z__')[-lib-]translit(libfile,'A-Z./','a-z__'),
      [  --with-]translit($1,'A-Z./','a-z__')[-lib-]translit(libfile,'A-Z./','a-z__')[=<file>
                          Full path location of the lib]libfile[.so(a) library file],
      [eval searchdir=`dirname $withval`; eval searchfile=`basename $withval`],
      [eval searchdir=$`echo "$PKGNM"_LIBDIR_PATH`])

dnl Find the library file (searches dynamic, then static)
    # echo Looking for libraries in $searchdir
    if test -z "$searchfile"; then
      TX_PATH_FILES(translit($1,'a-z./-','A-Z____')[_LIB_]translit(libfile,'a-z./-','A-Z____'),[lib]libfile[.so] [lib]libfile[.dylib] [lib]libfile[.a] libfile[.lib],[],[$searchdir])
    else
      TX_PATH_FILE(translit($1,'a-z./-','A-Z____')[_LIB_]translit(libfile,'a-z./-','A-Z____'),[$searchfile],[],[$searchdir])
    fi
    eval tempvar="$`echo "$PKGNM"_LIB_"$FILENM"`"

    dnl Confirm the integrity of the file, skip absolute path checking
    test_integ=""
    if test -n "$tempvar"; then
      for tx_itr in `echo $tx_liblist_glob | tr ':' ' '`; do
        abspath=`echo $tx_itr | grep '^/'`
        if test -n "$abspath"; then
          test_integ=""
          break
        fi
        TX_LOCATE_PKG_EXTRACT_ROOT([`dirname $tempvar`],[$tx_itr])
        if test -n "$return"; then
          test_integ=$test_integ" "$return
        fi
      done
    fi
    dnl We're going to iteratively grep each string in test_integ, and then pick the one that is the
    dnl greatest common factor
    cand=""
    if test -n "$test_integ"; then
      for tx_iitr in $test_integ; do
        loopcontrol=true
        for tx_jitr in $test_integ; do
          tmp=`echo $tx_jitr | grep $tx_iitr`
          if test -z "$tmp"; then
            loopcontrol=false
            break
          fi
        done
        if $loopcontrol; then
          cand="$tx_iitr"
          break
        fi
      done
    fi
dnl Now that we have a candidate, we're going to compare it to our
dnl already acquired tx_suspect_dir.  If it does not match
dnl tx_suspect_dir, we will zero out $tempvar because it's from a different
dnl build installation.  If tx_suspect_dir is empty, we'll fill it with
dnl the candidate.  If the candidate is empty, we're going to let it
dnl pass (pathological case).
    if test -z "$cand"; then
      tempvar="$tempvar"
    elif test -n "$cand" -a -z "$tx_suspect_dir"; then
      tx_suspect_dir=$cand
dnl jrc 11oct08: Sometimes on is looking for only a library
    elif test x"$cand" != x"$tx_suspect_dir"; then
      tempvar="-- conflict --"
      AC_MSG_WARN([Candidate, $cand, not same as suspect, $tx_suspect_dir.])
    fi

dnl If the file is found, append the flag to _LIB and write the directory to _LIBDIR, otherwise warn
    if test -n "$tempvar" -a x"$tempvar" != x"-- conflict --"; then
      TX_LOCATE_PKG_UNIQUE_CAT([`dirname $tempvar`],[tx_lib_paths])
      TX_LOCATE_PKG_UNIQUE_CAT([libfile],[tx_libs])
      eval "$PKGNM"_LIBDIR=`dirname $tempvar`
    elif test x"$tempvar" = x"-- conflict --"; then
      if test -n "$searchfile"; then
        AC_MSG_WARN([Conflicting $searchfile library file found])
      else
        AC_MSG_WARN([Conflicting lib]libfile[.so(a) library file found])
      fi
      tx_failure="yes"
    else
      if test -n "$searchfile"; then
        AC_MSG_WARN([Could not find $searchfile library file])
      else
        AC_MSG_WARN([Could not find lib]libfile[.so(a) library file])
      fi
      tx_failure="yes"
    fi

dnl Finally, check to see if if any dynamic libs were found
    tmp=`echo $tempvar | grep '\.so$'`
    if test -n "$tmp"; then
      tx_dynamic=yes
    fi

    dnl Write to summary file if defined
    if test -n "$config_summary_file"; then
      if test -n "$tempvar" -a x"$tempvar" != x"-- conflict --"; then
        TX_PRINT_VAR(${PKGNM}_LIB_$FILENM)
      elif test x"$tempvar" = x"-- conflict --"; then
        TX_PRINT_VAR(${PKGNM}_LIB_$FILENM, [-- conflict --])
      else
        TX_PRINT_VAR(${PKGNM}_LIB_$FILENM, [-- not found --])
      fi
    fi
  ])

])

dnl This is a generic function used by TX_LOCATE_PKG that keeps tracks
dnl of pseudo-lists and only adds to a list an item if it's a unique
dnl entity.  Effectively it is the python
dnl 	if not item in list: list.append(item)
dnl Argument 1 is "item"; argument 2 is a container containing "list".  List should be space separated.
AC_DEFUN([TX_LOCATE_PKG_UNIQUE_CAT],[

  testbool="1"
  if test -z "$3"; then
    tx_dlm=" "
  else
    tx_dlm=$3
  fi
  # echo 1 = $1
  # echo 2 = $2
  # echo tx_dlm = $tx_dlm
  tx_temp=`eval echo $$2 | tr "'$tx_dlm'" "' '"`
  for itr in $tx_temp; do
    if test "$1" = "$itr"; then
      testbool="0"
      break
    fi
  done
  if test "$testbool" = "1"; then
    $2="$$2${tx_dlm}$1"
  fi

])

dnl This is another generic function used by TX_LOCATE_PKG that is
dnl used to determine a root directory given a full directory and
dnl a relative path.  Basically, it confirms that the relative path is
dnl indeed the suffix and then sets up a container variable
dnl $return containing the root of the path.
dnl If the suffix is not found, then $return will be empty
AC_DEFUN([TX_LOCATE_PKG_EXTRACT_ROOT],[

  return=""
  arg1="$1"
  arg2="$2"

  while true; do
    pop1=`basename $arg1`
    pop2=`basename $arg2`
    if test x"$pop1" = x"$pop2" -a x"$pop1" != x"/"; then
      arg1=`dirname $arg1`
      arg2=`dirname $arg2`
      continue
    else
      break
    fi
  done
  if test x"$arg2" = x"."; then
    return="$arg1"
  fi

])
