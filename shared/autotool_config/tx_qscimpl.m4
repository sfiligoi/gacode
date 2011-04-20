dnl ######################################################################
dnl
dnl File:	tx_qt4.m4
dnl
dnl Purpose:	Determine the locations of QT includes and libraries
dnl
dnl Version:	$Id: tx_qscimpl.m4 3282 2009-11-22 16:00:39Z cary $
dnl
dnl Tech-X configure system
dnl
dnl ######################################################################

dnl Required for TX_PATH_FILE(S) to work
builtin(include, config/txsearch.m4)

dnl ######################################################################
dnl
dnl find the location of the qt installation
dnl
dnl ######################################################################

if test -n "$SUPRA_SEARCH_PATH"; then
  for i in `echo $SUPRA_SEARCH_PATH | tr ':' ' '`; do
    QSCIMPL_PATH="$QSCIMPL_PATH:$i/qscimpl"
  done
else
  echo "SUPRA_SEARCH_PATH not set!"
fi

TX_LOCATE_PKG(
  [qscimpl],
  [$QSCIMPL_PATH],
  [qscimpl.h],
  [txqtrol, txqmovie, txqattrib, txqeditor, txq3d, txq2d, txqbase, txmodel],
  [qt/include/Qt],
  [lib/$COMPDIR])

