AC_MSG_CHECKING(command-line arguments)
if test ! -f $abs_top_srcdir/.config.opts -o \
        ! -f $abs_top_srcdir/.config.with; then
    cat <<_
setup files not found.
(Run config/cleanconf.sh if you want arg-checking functionality in this script.)
_
else
    for arg in `set | grep '^with_' | \
			 sed -e 's/^with_//' -e 's/=.*$//'`; do
        if ! grep '^'$arg'$' $abs_top_srcdir/.config.with > /dev/null; then
	    echo ' '
	    AC_MSG_ERROR(Unknown arg --with-$arg)
  	fi
    done
    for arg in `set | grep '^enable_' | \
			 sed -e 's/^enable_//' -e 's/=.*$//'`; do
        if ! grep '^'$arg'$' $abs_top_srcdir/.config.opts > /dev/null; then
	    echo ' '
	    AC_MSG_ERROR(Unknown arg --enable-$arg (or --disable-$arg))
  	fi
    done
    echo "success"
fi
