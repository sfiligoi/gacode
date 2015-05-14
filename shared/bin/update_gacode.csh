date
echo '     Setting up the environment'
env | grep -i gacode
source $GACODE_ROOT/shared/env/env.${GACODE_PLATFORM}.tcsh
source $GACODE_ROOT/shared/bin/gacode_setup.tcsh
echo 'Done setting up the environment'
echo 'Changing to directory /tmp/gacode'
(cd /tmp/gacode && git pull) || (cd /tmp && rm -rf gacode && git clone git@github.com:gafusion/gacode.git )
cd /tmp/gacode

echo '     Cleaning GACODE; log in /tmp/gacode.make_clean.log'
make clean > /tmp/gacode.make_clean.log
echo 'Done cleaning GACODE'

echo '     Making GACODE; log in /tmp/gacode.make.log'
make > /tmp/gacode.make.log
echo 'Done making GACODE'

tail -n5 /tmp/gacode.make.log

set tmp_dir=`dirname $GACODE_ROOT`

cp -ru /tmp/gacode ${tmp_dir} 
python -m compileall ${GACODE_ROOT}/python
chgrp -R atom $GACODE_ROOT
chmod -R g+rwX $GACODE_ROOT
date
