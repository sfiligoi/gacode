#!/usr/bin/env python
#xmlrpc can take the place of a more complicated
#scheme such as CORBA.



from SimpleXMLRPCServer import SimpleXMLRPCServer
import sys ,string, os,time,re,shutil
from socket import *
#s_port = raw_input('Enter port number ( 10000 to 60000): ')
#port = int(s_port)
#filename   = raw_input('Enter fully qualified name   of rpc config file ')
filename = 'P_Nfreya_rpc_config'

#   \s+        one or more whitespace characters
#   \d+        one or more digits (i.e. an integer)
#  (\d+)      can extract the integer later
#  \W      any non-word character

real_in = r'(\d+)'
ipattern = r'\s*P_Nfreya_port_no\s*\W\s*' + real_in

# read the file and pick ot port number:
try :
    file = open(filename,'r')
except :
    print "Did not find ",filename
    sys.exit(1)
while 1:
    line = file.readline()
    if not line: break
    matchi = re.search(ipattern, line)
    if matchi :
        port = int(matchi.group(1))
        break
print "Selected port is port = ",port



local_host = gethostname()
if string.find(local_host,'lohan4')  >=0   :
                 server = SimpleXMLRPCServer(("192.5.166.128",port))
elif  string.find(local_host,'lohan5')  >=0  :
                  server = SimpleXMLRPCServer(("192.5.166.142",port))   #lohan5
elif  string.find(local_host,'lohan6')  >=0  :
                  server = SimpleXMLRPCServer(("192.5.166.186",port))   #lohan6
elif  string.find(local_host,'head')  >=0  :
                  server = SimpleXMLRPCServer(("192.5.166.186",port))   #lohan6
elif  string.find(local_host,'lohan7')  >=0  :
                  server = SimpleXMLRPCServer(("192.5.166.195",port))   #lohan7
elif  string.find(local_host,'node01')  >=0  :
                  server = SimpleXMLRPCServer(("192.5.166.195",port))   #lohan7
elif  string.find(local_host,'node06.cluster')  >=0  :
                  server = SimpleXMLRPCServer(("192.5.166.28",port))   #lohan12
elif  string.find(local_host,'lohan1')  >=0  :
                  server = SimpleXMLRPCServer(("192.5.166.181",port))   #lohan1
elif  string.find(local_host,'venusa')  >=0  :
                  server = SimpleXMLRPCServer(("192.5.166.160",port))  
elif  string.find(local_host,'star7.gat.com')  >=0  :
                  server = SimpleXMLRPCServer(("192.5.166.248",port))  
elif  string.find(local_host,'star8.gat.com')  >=0  :
                  server = SimpleXMLRPCServer(("192.5.166.249",port))
elif  string.find(local_host,'star9.gat.com')  >=0  :
                  server = SimpleXMLRPCServer(("192.5.166.250",port))  
elif  string.find(local_host,'star10.gat.com')  >=0  :
                  server = SimpleXMLRPCServer(("192.5.166.251",port))
elif  string.find(local_host,'star11.gat.com')  >=0  :
                  server = SimpleXMLRPCServer(("192.5.166.252",port))
elif  string.find(local_host,'star12.gat.com')  >=0  :
                  server = SimpleXMLRPCServer(("192.5.166.252",port))  
                                    
server.register_introspection_functions()

#my reminder of how to use:
# Register a function under a  name
#def adder_function(x,y):
#    return x + y
#server.register_function(adder_function, 'add')

def run_P_Nfreya_com(nprocs,statefile,outpt_statefile,run_directives_file,\
              local_dir,remote_dir,P_Nfreya_dir,hostname,client,mpi_hostfile):
    #local dir is dir on server side in which to run P_Nfreya
    #P_Nfreya_dir  is the dir on server side where
    #executable P_Nfreya will be found
    #remote dir is dir on client side to which output will be shipped
    error_code = 0
    #os.mkdir(	path[, mode])  mode =permissions
    #os.rmdir(local_path) !remove directory
    #os.remove(file)
    cm =' '
    os.chdir(local_dir)
    logfile = 'rpc_P_Nfreya_run_log'
    if local_host == 'lohan5.gat.com' :
       mpirun = '/c/source/PGI/pgi/linux86/9.0/mpi/mpich/bin/mpirun'
       error_code = 1 # no valid mpi copilation on loahn4,5
       return error_code
    if local_host == 'head' or local_host == 'node01' or local_host == 'node02'\
       or local_host == 'node03' or local_host == 'node04' \
       or local_host == 'node05' or local_host == 'node06' :
        #mpirun = '/export/mvapich2-1.4rc1_open-fabric/pgi-9.0-1/bin/mpirun_rsh'
        #mpirun_intel '/export/mvapich2-1.4rc1_open-fabric/intel-11.0/bin/mpirun_rsh
        #mpirun = 'mpirun_pgf '
        mpirun ='/export/mvapich2-1.4rc1_open-fabric/pgi-9.0-1/bin/mpirun_rsh'
        #mpirun = 'mpirun_intel'
        cm =' -hostfile ' + mpi_hostfile +' '
    if local_host == 'lohan1.gat.com' :
       mpirun = '/task/imd/apps/mpi/mpich2/pgi/64/bin/mpirun'
       cm =' '
    if local_host == 'venusa' :
       mpirun = '/u/stjohn/mpich2/gfortran/bin/mpirun -launcher ssh'
       cm =' -f ' + mpi_hostfile +' '
    if local_host == 'star7.gat.com' :
       mpirun = '/u/stjohn/mpich2/gfortran/bin/mpirun -launcher ssh'
       cm =' -f ' + mpi_hostfile +' '
    if local_host == 'star8.gat.com' :
       mpirun = '/u/stjohn/mpich2/gfortran/bin/mpirun -launcher ssh'
       cm =' -f ' + mpi_hostfile +' '
    if local_host == 'star9.gat.com' :
       mpirun = '/u/stjohn/mpich2/gfortran/bin/mpirun -launcher ssh'
       cm =' -f ' + mpi_hostfile +' '
    if local_host == 'star10.gat.com' :
       mpirun = '/u/stjohn/mpich2/gfortran/bin/mpirun -launcher ssh'
       cm =' -f ' + mpi_hostfile +' '
    if local_host == 'star11.gat.com' :
       mpirun = '/u/stjohn/mpich2/gfortran/bin/mpirun -launcher ssh'
       cm =' -f ' + mpi_hostfile +' '
    if local_host == 'star12.gat.com' :
       mpirun = '/u/stjohn/mpich2/gfortran/bin/mpirun -launcher ssh'
       cm =' -f ' + mpi_hostfile +' '
    # run P_Nfreya 
    command = mpirun +'  -np ' + str(nprocs) + cm + P_Nfreya_dir +' '\
              + statefile +' ' + run_directives_file + ' ' + logfile
    print 'command = ',command
    os.system( command )


    
    #ship P_Nfreya output statefile back to client:
    shipclient = client + ':' 
    if shipclient.find('venusa') >=0  :
        shipclient =" "
    command = ' scp ' + outpt_statefile  +' ' + shipclient  + remote_dir
    print 'command = ',command
    os.system(command)
    
    return error_code
    #return command
server.register_function(run_P_Nfreya_com, 'run_P_Nfreya')


def set_work_dir_com(local_dir ) :
    ignore_errors = 5 > 3
    shutil.rmtree(local_dir,ignore_errors)
    os.mkdir(local_dir)
    os.chdir(local_dir)
    filename ="leave_footprint"
    FILE = open(filename,"w")
    linetwrt =[]
    linetwrt.append( "testing removal/creation of files")
    linetwrt.append("\n next line")
    linetwrt.append("\n")
    FILE.writelines(linetwrt)
    FILE.close()
    return
server.register_function(set_work_dir_com, 'set_work_dir')


#reminder of use:
# Register an instance; all the methods of the instance are 
# published as XML-RPC methods (in this case, just 'div').
class MyFuncs:
    def div(self, x, y): 
        return x // y
    
server.register_instance(MyFuncs())

# Run the server's main loop
server.serve_forever()

