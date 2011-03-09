## ##########################################################################
##
## File:	Makefile.am
##
## Purpose:     Run automake in top directory to make Makefile.in
##	      in all subdirectories
##
## Version:     $Id: Makefile.am 1201 2009-08-18 00:22:45Z mmiah $
##
## Copyright Tech-X Corporation (2002-2008)
## See LICENSE file for conditions of use.
##
## ##########################################################################

# Above not used as must decide on case-by-case basis what is installed
SOURCE_FILES = ga_transport.csh ga_transport.sh

install-data-local:
	$(INSTALL) -d $(prefix)/share
	for i in $(SOURCE_FILES); do \
	  cp -f $$i $$i.bak; \
	  sed -e "s%$(abs_top_builddir)%$(prefix)%g" <$$i.bak > $$i; \
	  $(INSTALL) $$i $(prefix)/share; \
	  mv -f $$i.bak $$i; \
	done

uninstall-local:
	for i in $(SOURCE_FILES); do \
	  rm -f $(prefix)/share/$$i; \
	done


# Test -e needed to check for file or link

clean-local:
	:

## #######
##
## Other targets
##
## #######

include ../config/filelists.am

