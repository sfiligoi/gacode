dnl ######################################################################
dnl
dnl File:	tx_mysql.m4
dnl
dnl Purpose:	Determine where the mysql files are.
dnl
dnl Version:	$Id: tx_mysql.m4 3366 2010-01-15 18:43:10Z dws $
dnl
dnl Copyright 2001-2010, Tech-X Corporation.  Redistribution allowed provided
dnl this copyright statement remains intact.
dnl
dnl ######################################################################

builtin(include, config/txsearch.m4)

dnl Default search paths for mysql
mysql_searchpath="$HOME"/mysql:/contrib/mysql:/usr/local/mysql:/opt/mysql

dnl Locate activemq-cpp package
TX_LOCATE_PKG(
	[mysql],
	[$mysql_searchpath],
	[mysql.h],
	[mysqlclient],
	[include/mysql:include],
	[lib/mysql:lib]
)

