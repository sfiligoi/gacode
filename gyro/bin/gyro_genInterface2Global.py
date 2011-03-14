#!/usr/bin/python
# This is a simple little file to generate the gyro_dump_interface.f90 file
# based on gyro_interface.f90. 
#It produces 1 filewhich shows 
# interface variable name = interface value, global variale name = global variable value

import sys, os
import re
import math

def print_header(fw):
      fw.write("!--------------------------------------------------------------\n")
      fw.write("! gyro_dump_interface.f90\n")
      fw.write("!\n")
      fw.write("! PURPOSE:\n")
      fw.write("!  Complete dump of interface parameters matched with.\n")
      fw.write("!\n")
      fw.write("! NOTES:\n")
      fw.write("!  Also necessary for debugging exactly all interface variables used in a simulation.\n")
      fw.write("!\n")
      fw.write("! \n")
      fw.write("!---------------------------------------------------------------\n")
      fw.write("\n")
      fw.write(" subroutine gyro_dump_interface\n")
      fw.write("\n")
      fw.write(" use gyro_globals\n")
      fw.write(" use gyro_interface\n")
      fw.write(" implicit none\n")
      fw.write("!\n")
      fw.write(" integer :: funit=21\n")
      fw.write(" character(100) :: fname\n")
      fw.write("!-------------------------\n")
      fw.write(" if (i_proc == 0) THEN\n")
      fw.write('   fname=trim(path)//"dump_interface_inputs.txt"\n')
      fw.write("   open(unit=funit,file=trim(fname),status='replace')\n")
      fw.write("\n")
      fw.write("\n")

def print_footer(fw):
      fw.write("  close(funit)\n")
      fw.write(" end if\n")
      fw.write("!\n")
      fw.write(" end subroutine gyro_dump_interface\n")
      fw.write("!--------------------------------------------------------\n")
      fw.write("\n")
      fw.write(" subroutine dumpIntInterface(unit1, inInterfaceVarName, inInterfaceVar, globalVarName, globalVar)\n")
      fw.write("\n")
      fw.write(" character(*), intent(in) :: inInterfaceVarName,globalVarName  \n")
      fw.write(" integer, intent(in) :: inInterfaceVar, globalVar, unit1\n")
      fw.write(" character(*), parameter :: FMT1 = \"(A30,A3,I5.2,A3,A30,A3,I5.2)\"\n")
      #fw.write(" character(*), parameter :: FMT1 = \"(4A30)\"\n")
      fw.write("\n")
      fw.write(' write(unit1,FMT1) inInterfaceVarName,  " = ",inInterfaceVar, " , ", globalVarName, " = ",  globalVar\n')
      fw.write(' return\n')
      fw.write(' end subroutine dumpIntInterface\n')
      fw.write("!\n")
      fw.write("!--------------------------------------------------------\n")
      fw.write("\n")
      fw.write(" subroutine dumpRealInterface(unit1, inInterfaceVarName, inInterfaceVar, globalVarName, globalVar)\n")
      fw.write("\n")
      fw.write(" character(*), intent(in) :: inInterfaceVarName,globalVarName  \n")
      fw.write(" real, intent(in) :: inInterfaceVar, globalVar\n")
      fw.write(" integer, intent(in) ::  unit1\n")
      fw.write(" character(*), parameter :: FMT1 = \"(A30,A3,F10.6,A3,A30,A3,F10.6)\"\n")
      #fw.write(" character(*), parameter :: FMT1 = \"(4A30)\"\n")
      fw.write("\n")
      fw.write(' write(unit1,FMT1) inInterfaceVarName,  " = ",inInterfaceVar, " , ", globalVarName, " = ",  globalVar\n')
      fw.write(' return\n')
      fw.write(' end subroutine dumpRealInterface\n')

def buildGyroDict(filename,gyroDict):
    fr=open(filename,"r")
    lines = fr.readlines()
    inSegment = False
    for line in lines:
      if "subroutine map_global2interface" in line: 
        inSegment = True
      if "end subroutine map_global2interface" in line:
        break
      if inSegment and "=" in line and ":" not in line and "==" not in line:
        sline = line.split("=")
        gyroDict[sline[1].strip()]=sline[0]
      

varPairs = dict()
buildGyroDict("gyro_interface.f90",varPairs)
fr=open("gyro_read_input.f90","r")
fw=open("gyro_dump_interface.f90","w")
print_header(fw)
lines=fr.readlines()
for line in lines:
  if "readbc" in line:
    ttype = line.split("readbc_")[1].split("(")[0].strip()
    if ";" not in line:
      globalVar = line.split("(")[1].split(")")[0].lstrip().strip()
    else:
      globalVar = line.split(";")[1].split("=")[0].lstrip().strip()
    if ttype  == "int":
      subName="dumpIntInterface"
    elif ttype=="real":
      subName="dumpRealInterface"
    else:
      print "Error:  Type not recognized: "+ttype
    wrline1 = "      call "+subName+ "(21,\""+varPairs[globalVar]+"\","+varPairs[globalVar]+", &\n"
    wrline2 = "            \""+globalVar+"\","+globalVar+")\n"
    fw.write(wrline1)
    fw.write(wrline2)
print_footer(fw)





