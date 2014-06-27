#----------------------------------------------------------------------
# gacodeinput.py
#
# PURPOSE:
#  Collection of classes for parsing of GACODE free-format input files.
#
# NOTES: 
#  SimpleInput  : input.gyro, input.neo, input.tglf 
#  ProfileInput : input.profiles
#  ManagerInput : input.tgyro [see tgyro/bin/tgyro_parse.py]
#
# AUTHOR(S): 
#  Jeff Candy
#----------------------------------------------------------------------

import string
import os

class SimpleInput:
    """Input parser for simple input.* files"""
    def __init__(self):
        self.data_dict = {}
        self.data_orderlist = []
        self.dep_dict = {}
        self.dep_orderlist = []
        self.user_dict = {}
        self.error = 0
        self.error_msg = ""
        self.extension = ".gen"

    def add(self,param,default):
        self.data_dict[param]=default
        self.data_orderlist.append(param)

    def dep(self,param,default):
        self.dep_dict[param]=default
        self.dep_orderlist.append(param)

    def printdebug(self):
        print self.data_dict
        print self.data_orderlist
        print self.dep_dict
        print self.dep_orderlist
        print self.user_dict
        if self.error == 1:
            print self.error_msg

    def printmsg(self):
        if self.error == 1:
            print self.error_msg
        
    def set_extension(self,text):
        self.extension = text

    def read_input(self,inputfile):
        # 1. read user input file
        for line in open(inputfile,'r').readlines():

            # Remove leading and trailing whitespace from line
            line = string.strip(line)

            # Skip blank lines
            if len(line) > 0 and line[0] != '#':
                x = string.splitfields(line,'=')
                y = string.splitfields(x[1],'#')
                #y = string.splitfields(y[0],' ')
                arg = string.strip(x[0])
                val = string.strip(y[0])

                self.user_dict[arg] = val

        # 2. build complete input file, looking for errors
        for x in self.user_dict.keys():
            if self.data_dict.has_key(x) == 1:
                self.data_dict[x] = self.user_dict[x]
            elif self.dep_dict.has_key(x) == 1:
                self.error=1
                self.error_msg=self.error_msg+'ERROR: (gacodeinput) Deprecated parameter '+x+'\n'
                self.error_msg=self.error_msg+'       '+self.dep_dict[x]+'\n'
            else:
                self.error=1
                self.error_msg=self.error_msg+"ERROR: (gacodeinput) Bogus parameter "+x+'\n'

        if self.error == 0:
            f=open(inputfile+self.extension,'w')
            for x in self.data_orderlist:
                f.write(self.data_dict[x]+'  '+x+'\n')
        

class ProfileInput:
    """Input parser for input.profiles"""
    def __init__(self):
        self.data_dict = {}
        self.data_orderlist = []
        self.user_dict = {}
        self.error = 0
        self.error_msg = ""
        self.extension = ".gen"

    def add(self,param,default):
        self.data_dict[param]=default
        self.data_orderlist.append(param)

    def printmsg(self):
        if self.error == 1:
            print self.error_msg

    def set_extension(self,text):
        self.extension = text

    def read_input(self,inputfile):
        # 1. read user input file
        for line in open(inputfile,'r').readlines():

            # Remove leading and trailing whitespace from line
            line = string.strip(line)

            # Skip blank lines
            if len(line) > 0 and line[0] != '#':
                x = string.splitfields(line,'=')
                y = string.splitfields(x[1],'#')
                arg = string.strip(x[0])
                val = string.strip(y[0])

                self.user_dict[arg] = val

        # 2. build complete input file, looking for errors
        for x in self.user_dict.keys():
            if self.data_dict.has_key(x) == 1:
                self.data_dict[x] = self.user_dict[x]
            else:
                self.error=1
                self.error_msg=self.error_msg+"ERROR: Bogus parameter "+x+'\n'


    def read_profile(self,inputfile):

        outputfile = inputfile+self.extension

        file_out = open(outputfile,'w')
        file_temp = open('parse_temp','w')

        profile_data = []

        for line in open(inputfile,'r').readlines():

            line = string.strip(line)

            # Split inputfile (input.profiles) into scalar and 
            # vector data:

            if (len(line) > 0) and (line[0] != '#'):
                if string.find(line,'=') > -1:
                    # Write scalar data into temp file
                    file_temp.write(line+'\n')
                else:
                    # Originally '  ' was the pattern:       
                    #  data = string.splitfields(line,'   ')
                    data = string.split(line)
                    ncol = len(data)
                    for j in range (0,ncol):
                        # Save vector data in variable v.
                        profile_data.append(data[j])

        file_temp.close()

        # Parse the scalar data in file_temp
        self.read_input('parse_temp')

        # Compute number of rows for profile data
        nrow = int(self.data_dict['N_EXP'])
        nblock = len(profile_data)/(nrow*ncol)

        file_out.write(str(ncol)+' ncol\n')
        file_out.write(str(nblock)+' nblock\n')

        for x in self.data_orderlist:
            file_out.write(self.data_dict[x]+'  '+x+'\n')

        # Write vector data
        for k in range (0,nblock):
            for j in range (0,ncol):
                for i in range(0,nrow):
                    indx = ncol*i+j+k*nrow*ncol
                    file_out.write(profile_data[indx]+'\n')
    
        # Clean up temporary file
        os.system('rm parse_temp')

class ManagerInput:
    """Input parser for input.tgyro"""
    def __init__(self):
        self.data_dict = {}
        self.data_orderlist = []
        self.dep_dict = {}
        self.dep_orderlist = []
        self.user_dict = {}
        self.error = 0
        self.error_msg = ""
        self.extension = ".gen"
        self.tag = 'DIR'
        self.slavepath = []
        self.slaveproc = []
        self.overlayfile = []
        self.sum_proc = 0

    def add(self,param,default):
        self.data_dict[param]=default
        self.data_orderlist.append(param)

    def dep(self,param,default):
        self.dep_dict[param]=default
        self.dep_orderlist.append(param)

    def printmsg(self):
        if self.error == 1:
            print self.error_msg

    def write_proc(self,datafile):
        f = open(datafile,'w')
        f.write(str(len(self.slavepath))+'\n')
        f.write(str(self.sum_proc))
        f.close()

    def set_extension(self,text):
        self.extension = text

    def getcpu(self,path):

        for line in open(path+'/input.gyro.gen','r').readlines():
            line = string.strip(line)
            data = string.splitfields(line,'  ')
            if (data[1] == 'ENERGY_GRID'):
                ne = int(data[0])
            if (data[1] == 'PASS_GRID'):
                n_pass = int(data[0])
            if (data[1] == 'TRAP_GRID'):
                n_trap = int(data[0])
            if (data[1] == 'TOROIDAL_GRID'):
                n_n = int(data[0])

        return ne*(n_pass+n_trap)*n_n

    def strip_tag(self,datafile):

        file_input = open(datafile+'.input','w')

        n = 0
        for line in open(datafile,'r').readlines():
            line_s = string.strip(line)

            # Look for occurence of tag and put item in list.
            if (line_s[0:3] == self.tag):   
                n = n+1
                data = string.splitfields(line_s,' ')
                self.slavepath.append(data[1])
                self.slaveproc.append(data[2]) 
                self.overlayfile.append('overlay.'+str(n))
                file_overlay = open('overlay.'+str(n),'w')

                #----------------------------------------------------------
                # This loop writes each overlay parameter list to overlay.*
                for j in range(len(data)-3):
                    file_overlay.write(data[j+3]+'\n')
                file_overlay.close()
                #----------------------------------------------------------
 
            else:

                file_input.write(line)

        file_input.close()


    def read_input(self,datafile,gyro_start):

        # Eventual output file
        outfile = datafile+self.extension

        # Split datafile into datafile.input (pure input) 
        # and overlay.n, where n is the number of tags.
        self.strip_tag(datafile)

        # Parse simple input part
        z = SimpleInput()
        z.data_dict = self.data_dict
        z.data_orderlist = self.data_orderlist
        z.dep_dict = self.dep_dict
        z.read_input(datafile+'.input')
        self.error = z.error
        self.error_msg = z.error_msg

        # Something went wrong in SimpleInput (bogus parameter),
        # so let's just return
        if self.error == 1:
            return
        
        # File generated above
        genfile = datafile+'.input'+z.extension 

        n_path = len(self.slavepath)

        os.system('mv '+genfile+' '+outfile)
        os.system('rm '+datafile+'.input')
        file_outfile = open(outfile,'a')
        file_outfile.write(str(n_path)+'\n')

        # Logging
        print 'INFO: (gacodeinput) Number of code instances: '+str(n_path)

        # Add special entries (DIR) to output file, and overlay
        # extra parameters onto GYRO or TGLF input files.
        #
        # NOTE:
        #  - GYRO directories must be of the form GYRO*
        #  - IFS directories must be of the form IFS*
        #  - TGLF directories must be of the form TGLF*
        #  - FUN directories must be of the form FUN*

        for p in range(len(self.slavepath)):
            self.sum_proc = self.sum_proc + int(self.slaveproc[p])
            basedir = self.slavepath[p]
            file_outfile.write(basedir+' '+self.slaveproc[p]+'\n') 

            if basedir[0:3] == 'IFS':
                print 'INFO: (gacodeinput) Detected '+basedir+'; CPU_max=1'

            elif basedir[0:3] == 'FUN':
                 print 'INFO: (gacodeinput) Detected '+basedir+'; CPU_max=1'
               
            elif basedir[0:4] == 'TGLF':
                basefile = basedir+'/input.tglf' 
                tempfile = basedir+'/input.tglf.temp' 
                file_base = open(basefile,'r')
                file_temp = open(tempfile,'w')

                for line in file_base.readlines():
                    if line[0:18] <> "# -- Begin overlay":
                        file_temp.write(line)
                    else:
                        break

                file_base.close()
                file_temp.close()

                os.system('echo "# -- Begin overlay [Add parameters above this line]" >> '+tempfile)
                os.system('cat '+self.overlayfile[p]+' >> '+tempfile)
                os.system('mv '+tempfile+' '+basefile)

                os.system('tglf -i '+basedir+' -p $PWD')
                print 'INFO: (gacodeinput) Processed input.* in '+basedir+'; CPU_max=1'
                
            else:
                basefile = basedir+'/input.gyro' 
                tempfile = basedir+'/input.gyro.temp' 
                file_base = open(basefile,'r')
                file_temp = open(tempfile,'w')

                for line in file_base.readlines():
                    if line[0:18] <> "# -- Begin overlay":
                        file_temp.write(line)
                    else:
                        break

                file_base.close()
                file_temp.close()

                os.system('echo "# -- Begin overlay" >> '+tempfile)
                os.system('cat '+self.overlayfile[p]+' >> '+tempfile)
                os.system('mv '+tempfile+' '+basefile)

                os.system('gyro -i '+basedir+' -n 1 -nomp 1 -p $PWD '+gyro_start)
                cpu = self.getcpu(basedir)
                print 'INFO: (gacodeinput) Processed input.* in '+basedir+'; CPU_max='+str(cpu)

            os.system('rm '+self.overlayfile[p])

        print 'INFO: (gacodeinput) Required MPI tasks in TGYRO: '+str(self.sum_proc)
        
