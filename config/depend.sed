##########
#
# File:		depend.sed
#
# Purpose:	sed script for cleaning up garbage put out by some compilers
#
# $Id: reform.sed 2178 2007-02-05 12:43:46Z cary $
#
##########
# IBM puts a space before semicolon, which it then does not like
s/\.o :/.o:/g
# PGI put these overflow messages in
s/^overflow.*$//
s/overflow error//
# Sometimes the depends are accumulated on a line
s/ \([a-zA-Z][a-zA-Z0-9_]*.o:\)/\
\1/g
