#!/bin/sh

echo "Print all data to standard out"
cmd="geqdsk.py -a -f g118898.03400"
echo "Executing: "$cmd
sleep 10 
eval $cmd

sleep 7 ; echo; echo
echo "Print the outer boundary data"
cmd="geqdsk.py -v rbbbs,zbbbs -f g118898.03400"
echo "Executing: "$cmd
sleep 8 
eval $cmd



