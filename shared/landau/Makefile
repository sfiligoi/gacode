include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

OBJECTS = landau.o \
          gyrotransformation.o

.PHONY: clean all

.SUFFIXES:
.SUFFIXES: .o .f90 .f .F90

all: $(OBJECTS)

.f90.o .F90.o:
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -rf landau.o gyrotransformation.o
	rm -f ${GACODE_ROOT}/modules/landau.mod
	rm -f ${GACODE_ROOT}/modules/gyrotransformation.mod

# module landau in landau.F90
gyrotransformation.o: landau.o
