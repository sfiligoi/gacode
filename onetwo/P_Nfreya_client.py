#!/usr/bin/env python
import xmlrpclib,sys,os


# this client (and server) code assumes that we are behind an
# appropriate firewall, hence no additional security features
# are installed !!!!!!!


#INPUT ARGUMENTS TO THIS SCRIPT:
#-------------------------------------------------------------------------------------------------
# nprocs                is number  of processors to run P_nfreya with on the server
#                       The system command in pgf FOrtran does not pass more
#                       than 15 arguments correctly we append either T or F to nprocs
#                       to indicate that the ufile is to be sent or not
#statefile_input        is the statefile file that P_Nfreya will use for startup
#statefile_output       is the statefile file that P_Nfreya procudes as output
#run_dirf               is name of run directives file (namelist file) that P_Nfreya  will read
#remote_dir             is directory on server in which P_nfreya will be run
#local_dir              is directory on local machine to which output will be shipped
#P_Nfreya_dir           is directory/executable  where  P_Nfreya  will be found
# server_name:          is the name of machine on which P_Nfreya will execute the master process
#                       valid server names are: gvine in server list below
# client  name         is name of host on which this client is running
#                      (machine to which p_Nfreya output statefile will be shipped)
#nubeam_nml            name of nubeam namelist to be shipped to server
#nubeam_ufile          name of ufile to be shipped to server
#mpi_hostfile          is the name of the mpi file that contains the list of hosts on which
#                      parallel execution of P_Nfreya will take place
#-------------------------------------------------------------------------------------------------

n=len(sys.argv)
print 'nargs =',n
print sys.argv
if  n != 15 :
    print 'P_Nfreya_client.py reports that number of arguments it was given is wrong:'
    print sys.argv
    print "Usage:", sys.argv[0], "nprocs port_no P_Nfreya_remote_dir P_Nfreya_remote_host MPI_hostfile \
    P_Nfreya_inpt_statefile P_Nfreya_outpt_statefile P_Nfreya_nml_filename nubeam_namelist nubeam_ufile \
    local_dir P_Nfreya_dir local_host dir_dispatch"
    sys.exit(1)


snprocs               = sys.argv[1] # exmaple snprocs = 14F
nprocs               = int(snprocs[0:-1]) # yields nprocs = 14
sent_ufile = snprocs[-1]  # yields F

port                 = sys.argv[2]
remote_dir           = sys.argv[3]
server_name          = sys.argv[4]    # machine on which P_Nreyawill be run
mpi_hostfile         = sys.argv[5]
statefile_input      = sys.argv[6]
statefile_output     = sys.argv[7]
run_dirf             = sys.argv[8]
nubeam_nml           = sys.argv[9]
nubeam_ufile         = sys.argv[10]
local_dir            = sys.argv[11]    # ship output from P_Nreya to this directory
                                       # on machine client_name
P_Nfreya_dir         = sys.argv[12]
wipe_remote_dir      = sys.argv[13]
client_name          = sys.argv[14]    # we are calling P_Nreyafrom this machine
#send_ufile           = sys.argv[15]
print 'nprocs =',nprocs
print 'port =',port
print 'remote_dir=', remote_dir
print 'server_name=', server_name
print 'mpi_hostfile =',mpi_hostfile
print 'statefile_input',statefile_input
print 'run_dirf =',run_dirf 
print 'nubeam_nml  =',nubeam_nml  
print 'nubeam_ufile =',nubeam_ufile
print 'local_dir =',local_dir
print 'P_Nfreya_dir =', P_Nfreya_dir
print 'client_name=',client_name
print 'wipe_remote_dir =',wipe_remote_dir
#print ' send ufile =', send_ufile

#server_list  = ['lohan6','lohan7','lohan8','benten','scyld']
server_list  = ['lohan5','lohan6','lohan7','lohan8','lohan9','lohan10','lohan11','lohan12','node01','lohan1','venusa','star7','star8','star9','star10','star11','star12']
if server_name not in server_list :
    print ' Error, host:',server_name ,' is not recognized'
    print 'Possible host names (depending on 32/64 bit rpc code)  are :'
    print server_list
    sys.exit(1)



    
#s = xmlrpclib.Server('http://localhost:9000')
#port ='9005'
site  ='http://' + server_name + '.gat.com'
server = site +  ':' + port
#s = xmlrpclib.Server('http://lohan1.gat.com:9000')
s = xmlrpclib.Server(server)
if wipe_remote_dir == 'WIPE' :
    try :
        print 'calling s.set_work_dir'
        print s.set_work_dir(remote_dir)
        sys.exit(0)
    except :
        sys.exit(2)

#ship namelist and iterdb to remote host:
#command = ' scp ' + statefile_input + ' ' + server_name +'.gat.com:' + remote_dir
#print 'command =',command
#os.system(command)
#command = ' scp ' + run_dirf +' ' + server_name +'.gat.com:' + remote_dir
#print 'command =',command
#os.system(command)
#command = ' scp ' + ' ' +  nubeam_nml + ' '  + server_name +'.gat.com:' + remote_dir
#print 'command =',command

command = ' scp ' + statefile_input +' ' + run_dirf +' ' +  nubeam_nml + ' ' + server_name +'.gat.com:' + remote_dir
print 'command =',command
os.system(command)

if sent_ufile == 'F'  :
    command = ' scp ' + nubeam_ufile +' ' + server_name +'.gat.com:' + remote_dir
    print 'sending ufile with command =',command
    os.system(command)
else :
    print 'skipping ufile send'
    




print 'P_nfreya_client reports:  making rpc call'
try :
    print s.run_P_Nfreya(nprocs,statefile_input,statefile_output,run_dirf,remote_dir,local_dir,\
                        P_Nfreya_dir,server_name,client_name,mpi_hostfile)
except :
    print 'error,connection on ',server_name,'refused'
    sys.exit(2)
# Print list of available methods
# print s.system.listMethods()
  
