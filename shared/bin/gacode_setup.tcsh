#--------------------------------------------
# Environment variable setup for gacode
#--------------------------------------------
#!/bin/tcsh

setenv PATH ${PATH}:$GACODE_ROOT/tgyro/bin
setenv PATH ${PATH}:$GACODE_ROOT/gyro/bin
setenv PATH ${PATH}:$GACODE_ROOT/cgyro/bin
setenv PATH ${PATH}:$GACODE_ROOT/neo/bin
setenv PATH ${PATH}:$GACODE_ROOT/tglf/bin
setenv PATH ${PATH}:$GACODE_ROOT/glf23/bin
setenv PATH ${PATH}:$GACODE_ROOT/le3/bin
setenv PATH ${PATH}:$GACODE_ROOT/shared/bin

if ( $?PYTHONPATH ) then
 setenv PYTHONPATH ${PYTHONPATH}:$GACODE_ROOT/python
else
 setenv PYTHONPATH $GACODE_ROOT/python
endif

if ( $?IDL_PATH ) then
 setenv IDL_PATH ${IDL_PATH}:$GACODE_ROOT/gyro/vugyro
else
 setenv IDL_PATH $GACODE_ROOT/gyro/vugyro
endif

if ( $?GACODE_PLATFORM ) then
  set GA = (BANACH DROP LOHAN SATURN VENUS)
  foreach server ($GA)
    if ( $GACODE_PLATFORM == $server ) then
      setenv HARVEST_HOST venus.gat.com
      setenv HARVEST_PORT 32000
    endif
  end
endif