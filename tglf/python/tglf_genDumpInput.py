#!/usr/bin/python
# This is a simple little file to generate the dump_input.f90 file
# based on tglf_read_input.f90.  It's nothing too fancy

def print_header(fw):
      fw.write("!--------------------------------------------------------------\n")
      fw.write("! tglf_dump_input.f90\n")
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
      fw.write(" subroutine tglf_dump_input\n")
      fw.write("\n")
      fw.write(" use tglf_interface \n")
      fw.write(" implicit none\n")
      fw.write("!\n")
      fw.write(" integer :: funit=21\n")
      fw.write(" integer :: gunit=22\n")
      fw.write(" character(100) :: fname\n")
      fw.write(" character(100) :: gname\n")
      fw.write("!-------------------------\n")
      fw.write('   fname=trim(tglf_dump_path_in)//"dump_inputs.txt"\n')
      fw.write("   open(unit=funit,file=trim(fname),status='replace')\n")
      fw.write('   gname=trim(tglf_dump_path_in)//"gen_input.dat"\n')
      fw.write("   open(unit=gunit,file=trim(gname),status='replace')\n")
      fw.write("\n")
      fw.write("\n")

def print_footer(fw):
      fw.write("  close(funit)\n")
      fw.write("!\n")
      fw.write(" end subroutine tglf_dump_input\n")
      fw.write("!--------------------------------------------------------\n")
      fw.write("\n")
      


fr=open("tglf_read_input.f90","r")
fw=open("tglf_dump_input.f90","w")
print_header(fw)
lines=fr.readlines()
for line in lines:
  if "read(" in line:
    globalVar = line.split("read(1,*)")[1].strip()
    wrline1 = " write(21,*) \""+globalVar+" = \","+globalVar+"\n"   
    wrline2 = " write(22,*) "+globalVar+"\n"   
    fw.write(wrline1)
    fw.write(wrline2)
print_footer(fw)
