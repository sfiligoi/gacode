#--------------------------------------------
# Environment variable setup for gacode
#--------------------------------------------
#!/bin/tcsh

setenv PATH ${PATH}:$GACODE_ROOT/tgyro/bin
setenv PATH ${PATH}:$GACODE_ROOT/gyro/bin
setenv PATH ${PATH}:$GACODE_ROOT/neo/bin
setenv PATH ${PATH}:$GACODE_ROOT/tglf/bin
setenv PATH ${PATH}:$GACODE_ROOT/shared/bin

setenv PYTHONPATH $GACODE_ROOT/shared/python

