#--------------------------------------------
# Environment variable setup for gacode
#--------------------------------------------
#!/bin/tcsh

setenv PATH ${PATH}:$GACODE_ROOT/tgyro/bin
setenv PATH ${PATH}:$GACODE_ROOT/gyro/bin
setenv PATH ${PATH}:$GACODE_ROOT/cgyro/bin
setenv PATH ${PATH}:$GACODE_ROOT/neo/bin
setenv PATH ${PATH}:$GACODE_ROOT/vgen/bin
setenv PATH ${PATH}:$GACODE_ROOT/tglf/bin
setenv PATH ${PATH}:$GACODE_ROOT/glf23/bin
setenv PATH ${PATH}:$GACODE_ROOT/le3/bin
setenv PATH ${PATH}:$GACODE_ROOT/shared/bin

if ( $?GACODE_ADD_ROOT ) then       
   setenv PATH=${PATH}:$GACODE_ADD_ROOT/freya/bin
endif

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
