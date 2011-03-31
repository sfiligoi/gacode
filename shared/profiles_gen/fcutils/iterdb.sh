#!/bin/sh
cmd="getIterDbData.py -i iterdb.125730 -o temp.pre"
echo "Executing: "$cmd
eval $cmd
echo "Done executing "

sleep 6 ; echo; echo
echo "Output in temp.pre."
cmd="cat temp.pre"
echo "Executing: "$cmd
eval $cmd

sleep 7 ; echo; echo
echo "Print available fields and redirect"
cmd="getIterDbData.py -p  iterdb.125730 > out; cat out"
echo "Executing: "$cmd
sleep 16 
eval $cmd

sleep 7 ; echo; echo
echo "Rescale temperatures to eV and rename"
cmd='getIterDbData.py -f "Te,Ti" -c "1000,1000" -r"electron_temperation,ion_temperature" -o temp2.pre iterdb.125730'
echo "Executing: "$cmd
eval $cmd
echo "Done executing "
echo "Output in temp2.pre."

sleep 7 ; echo; echo
cmd="cat temp2.pre"
echo "Executing: "$cmd
eval $cmd

