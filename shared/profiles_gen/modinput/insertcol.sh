#!/bin/bash
#
# inserts data into a given column of input.profiles
# 1st rgument is the name of the input parameter which will be replaced 
# 2nd argument is the name of a file containing the new data column  
# Example:
# ./insertcol.sh er extract_er
# will insert the data in extract_er into the column of the er variable in 
#   input.profiles 


param1="$1(" 
param2="$1[[:space:]](" 
inp="$2"
filename="$3"
: ${filename:="input.profiles"}
outputname="$4"
: ${outputname:="input.profiles_mod"}

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
   grep -n $param2 $filename | sed 's/:/\n/' | head -1 > temp3
else
   grep $param1 $filename | sed "s/$param1/\n$param1/" | head -n 1 > temp2
   grep -n $param1 $filename | sed 's/:/\n/' | head -1 > temp3
fi
npos=`cat temp3 | awk '{ print $1 }'`

head -$npos $filename  > $outputname
ntail=`expr $nlines - $npos`
ntail=`expr $ntail - $nline`

tail $filename -n $ntail > temptail 

nchar=`cat temp2 | wc -c`
colnumber=`expr $nchar / $colwidth`
colnumber=`expr $colnumber + 1`         # which column is the data in?
tail temp1 -n $nline > body
case $colnumber in
      1) pr -m -J -t body $inp  | gawk '{print " " $6 "   " $2 "   " $3 "   " $4 "   " $5}' | sed s/" -"/"-"/g >> $outputname ;;
      2) pr -m -J -t body $inp  | gawk '{print " " $1 "   " $6 "   " $3 "   " $4 "   " $5}' | sed s/" -"/"-"/g >> $outputname ;;
      3) pr -m -J -t body $inp  | gawk '{print " " $1 "   " $2 "   " $6 "   " $4 "   " $5}' | sed s/" -"/"-"/g >> $outputname ;;
      4) pr -m -J -t body $inp  | gawk '{print " " $1 "   " $2 "   " $3 "   " $6 "   " $5}' | sed s/" -"/"-"/g >> $outputname ;;
      5) pr -m -J -t body $inp  | gawk '{print " " $1 "   " $2 "   " $3 "   " $4 "   " $6}' | sed s/" -"/"-"/g >> $outputname ;;
esac

cat temptail >> $outputname

rm body
rm temptail
rm temp1
rm temp2
rm temp3
exit 0
