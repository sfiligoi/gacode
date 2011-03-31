#!/usr/bin/python
# This is a simple little file to generate the dump_input.f90 file
# based on gyro_read_input.f90.  It's nothing too fancy

def print_header(fw):
      fw.write("!--------------------------------------------------------------\n")
      fw.write("! gyro_dump_input.f90\n")
      fw.write("!\n")
      fw.write("! PURPOSE:\n")
      fw.write("!  Complete dump of input parameters.\n")
      fw.write("!\n")
      fw.write("! NOTES:\n")
      fw.write("!  Necessary for debugging exactly all input variables used in a simulation.\n")
      fw.write("!\n")
      fw.write("! \n")
      fw.write("!---------------------------------------------------------------\n")
      fw.write("\n")
      fw.write(" subroutine gyro_dump_input\n")
      fw.write("\n")
      fw.write(" use gyro_globals\n")
      fw.write(" implicit none\n")
      fw.write("!\n")
      fw.write(" integer :: funit=21\n")
      fw.write(" integer :: gunit=22\n")
      fw.write(" character(100) :: fname\n")
      fw.write(" character(100) :: gname\n")
      fw.write("!-------------------------\n")
      fw.write(" if (i_proc == 0) THEN\n")
      fw.write('   fname=trim(path)//"dump_inputs.txt"\n')
      fw.write("   open(unit=funit,file=trim(fname),status='replace')\n")
      fw.write('   gname=trim(path)//"gen_input.dat"\n')
      fw.write("   open(unit=gunit,file=trim(gname),status='replace')\n")
      fw.write("\n")
      fw.write("\n")

def print_footer(fw):
      fw.write("  close(funit)\n")
      fw.write(" end if\n")
      fw.write("!\n")
      fw.write(" end subroutine gyro_dump_input\n")
      fw.write("!--------------------------------------------------------\n")
      fw.write("\n")
      fw.write(" subroutine dumpIntVar(unit1, unit2, inVarName, inVar)\n")
      fw.write("\n")
      fw.write(" character(*), intent(in) :: inVarName\n")
      fw.write(" integer, intent(in) :: inVar, unit1, unit2\n")
      fw.write("\n")
      fw.write(' write(unit1,*) inVarName,  " = ",inVar\n')
      fw.write(' write(unit2,*) inVar\n')
      fw.write(' return\n')
      fw.write(' end subroutine dumpIntVar\n')
      fw.write("!\n")
      fw.write("!--------------------------------------------------------\n")
      fw.write("\n")
      fw.write(" subroutine dumpRealVar(unit1, unit2, inVarName, inVar)\n")
      fw.write("\n")
      fw.write(" character(*), intent(in) :: inVarName\n")
      fw.write(" real, intent(in) :: inVar\n")
      fw.write(" integer, intent(in) :: unit1, unit2\n")
      fw.write("\n")
      fw.write(' write(unit1,*) inVarName,  " = ",inVar\n')
      fw.write(' write(unit2,*) inVar\n')
      fw.write(' return\n')
      fw.write(' end subroutine dumpRealVar\n')



fr=open("gyro_read_input.f90","r")
fw=open("gyro_dump_input.f90","w")
print_header(fw)
lines=fr.readlines()
for line in lines:
  if "readbc" in line:
    if "subroutine" not in line:
      #print line
      ttype = line.split("readbc_")[1].split("(")[0].strip()
      if ";" not in line:
        globalVar = line.split("(")[1].split(")")[0]
#        print "no ; ", ttype, " ", globalVar
#        print line
      else:
        globalVar = line.split(";")[1].split("=")[0]
#        print "with ; ", ttype, " ", globalVar
#        print line
      if ttype  == "int":
        subName="dumpIntVar"
      elif ttype=="real":
        subName="dumpRealVar"
      else:
        print "Error:  Type not recognized: "+ttype
      wrline = "      call "+subName+ "(21,22,\""+globalVar+"\","+globalVar+")\n"
      fw.write(wrline)
print_footer(fw)
