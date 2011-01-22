#!/bin/bash
#
# Extracts a given column from input.profiles
# Argument is the name of the input parameter. 
# Example:
# ./extractcol.sh er
# will copy the er variable into the file extract_er 


param1="$1(" 
param2="$1[[:space:]](" 
filename=$2
: ${filename:="input.profiles"}
         
nblocks=8                # number of data blocks
colwidth=16              # number of characters (including spaces) in a data column
# number of non-commented lines which are not data column lines
nextracomm=`cut -c 1 $filename | grep -v "#" | grep -v [[:space:]] | grep -v "-" | wc -l`     
nlines=`cat $filename | wc -l`            #total number of lines in the file   
ncomm=`cat $filename | grep "#" | wc -l`  # number of comment lines
npure=`expr $nlines - $ncomm`             
npure=`expr $npure - $nextracomm`         # number of data column lines
nline=`expr $npure / $nblocks`            # number of lines in a data column    
grep -A $nline $param1 $filename > temp1  
isspace=`cat temp1 | wc -c`
if [ "$isspace" = "0" ] ; then 
   grep -A $nline $param2 $filename > temp1
   grep $param2 $filename | sed "s/$param2/\n$param2/" | head -n 1 > temp2
else
   grep $param1 $filename | sed "s/$param1/\n$param1/" | head -n 1 > temp2
fi
nchar=`cat temp2 | wc -c`
colnumber=`expr $nchar / $colwidth`
colnumber=`expr $colnumber + 1`
cat temp1 | awk "(FNR>1){ print \" \"\$$colnumber}" | sed s/" -"/"-"/ > "extract_$1"

rm temp1
rm temp2
exit 0
