c
c --- NOTE that common block /SHAY/ is now properly data-aligned.
c               PLEASE maintain this when you make modifications. ... JF
c
      parameter    (maxtimes = 50)  ! must match kinvrt in plot programs
c
      logical       suserho     , scsmult
      integer       irinshay    , iroutshay, irinshaycev,
     .              iroutshaycev, ishayform, ktimes
c
      common /shay/ smult, snexp, sbpexp, srexp, sbigrexp, stexp,
     .              sdtdrexp, srin, srout, skimult, schie(kj),
     .              schii(kj), rmbpav(kj), xkeshay(kj), xkishay(kj),
     .              wshay, timeshay, relaxshay, schiesv(kj),
     .              schiisv(kj), scheff(kj), dshay(kj), srincev,
     .              sroutcev,
     .              smultstder(maxtimes), skimultstder(maxtimes),
     .              smulta    (maxtimes), skimulta    (maxtimes),
     .              slim95    (maxtimes), skilim95    (maxtimes),
     .              sdenscale, xkangsch(kj), xkangschsv(kj), skimass,
     .              irinshay    , iroutshay, irinshaycev, ! all INTEGERs
     .              iroutshaycev, ishayform, ktimes     , ! all INTEGERs
     .              suserho     , scsmult                 ! all LOGICALs
c
c --- set wshay > 0 if Hsieh model of ke is to be used
c --- scsmult = true if constant multiplier in Hsieh model is
c ---           to be found (rather than input by user)
c --- suserho = true if flux surface average of Hsieh ke is to be used
c --- maxtimes is the maximum number of time points at which
c --- output for plotting will be collected.
c
