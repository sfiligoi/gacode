class profiles_genData:
    """profiles_gen output data class.

     """

    #---------------------------------------------------------------------------#
    # Methods

    def __init__(self,infile):
        """
        Constructor reads in data from directory and creates new object.
        """

        import string
        import numpy as np

        self.infile=infile

        # Initialize data
        self.data  = {}
        self.n_exp = 0
        self.geo   = {}

        # Read data
        row  = 0
        for line in open(infile,'r').readlines():
            row = row+1
            if line[0:5] == 'N_EXP':
                self.n_exp = int(string.splitfields(line,'=')[1])
            if line[0:8] == 'ARHO_EXP':
                self.max_rho = float(line.split('=')[1])
            if line[0:7] == '#rho(-)':
                break

        try:
            data = np.loadtxt(infile,skiprows=row)
            print '(INFO): (profiles_gen_plot) '+infile+' found.'
        except:
            raise IOError('(ERROR): (profiles_gen_plot) '+infile+' not found.')

        n = self.n_exp

        self.data['rho']     = data[0:n,0]
        self.data['rmin']    = data[0:n,1]
        self.data['polflux'] = data[0*n,2]
        self.data['q']       = data[0:n,3]
        self.data['omega0']  = data[0:n,4]
        
        self.data['rmaj']  = data[n:2*n,0]
        self.data['zmag']  = data[n:2*n,1]
        self.data['kappa'] = data[n:2*n,2]
        self.data['delta'] = data[n:2*n,3]
        self.data['zeta']  = data[n:2*n,4]
 
        self.data['ne']     = data[2*n:3*n,0]
        self.data['Te']     = data[2*n:3*n,1]
        self.data['ptot']   = data[2*n:3*n,2]
        self.data['z_eff']  = data[2*n:3*n,3]
        # NULL
        
        self.data['ni_1'] = data[3*n:4*n,0]
        self.data['ni_2'] = data[3*n:4*n,1]
        self.data['ni_3'] = data[3*n:4*n,2]
        self.data['ni_4'] = data[3*n:4*n,3]
        self.data['ni_5'] = data[3*n:4*n,4]
        
        self.data['ni_6'] = data[4*n:5*n,0]
        self.data['ni_7'] = data[4*n:5*n,1]
        self.data['ni_8'] = data[4*n:5*n,2]
        self.data['ni_9'] = data[4*n:5*n,3]
        self.data['ni_10'] = data[4*n:5*n,4]

        self.data['Ti_1'] = data[5*n:6*n,0]
        self.data['Ti_2'] = data[5*n:6*n,1]
        self.data['Ti_3'] = data[5*n:6*n,2]
        self.data['Ti_4'] = data[5*n:6*n,3]
        self.data['Ti_5'] = data[5*n:6*n,4]

        self.data['Ti_6'] = data[6*n:7*n,0]
        self.data['Ti_7'] = data[6*n:7*n,1]
        self.data['Ti_8'] = data[6*n:7*n,2]
        self.data['Ti_9'] = data[6*n:7*n,3]
        self.data['Ti_10'] = data[6*n:7*n,4]

        self.data['vtor_1'] = data[7*n:8*n,0]
        self.data['vtor_2'] = data[7*n:8*n,1]
        self.data['vtor_3'] = data[7*n:8*n,2]
        self.data['vtor_4'] = data[7*n:8*n,3]
        self.data['vtor_5'] = data[7*n:8*n,4]

        self.data['vtor_6'] = data[8*n:9*n,0]
        self.data['vtor_7'] = data[8*n:9*n,1]
        self.data['vtor_8'] = data[8*n:9*n,2]
        self.data['vtor_9'] = data[8*n:9*n,3]
        self.data['vtor_10'] = data[8*n:9*n,4]

        self.data['vpol_1'] = data[10*n:11*n,0]
        self.data['vpol_2'] = data[10*n:11*n,1]
        self.data['vpol_3'] = data[10*n:11*n,2]
        self.data['vpol_4'] = data[10*n:11*n,3]
        self.data['vpol_5'] = data[10*n:11*n,4]

        self.data['vpol_6'] = data[12*n:13*n,0]
        self.data['vpol_7'] = data[12*n:13*n,1]
        self.data['vpol_8'] = data[12*n:13*n,2]
        self.data['vpol_9'] = data[12*n:13*n,3]
        self.data['vpol_10'] = data[12*n:13*n,4]
 
        self.data['flow_beam'] = data[13*n:14*n,0]
        self.data['flow_wall'] = data[13*n:14*n,1]
        self.data['flow_mom']  = data[13*n:14*n,2]
        # NULL
        # NULL
        
        self.data['pow_e']     = data[14*n:15*n,0]
        self.data['pow_i']     = data[14*n:15*n,1]
        self.data['pow_ei']    = data[14*n:15*n,2]
        self.data['pow_e_aux'] = data[14*n:15*n,3]
        self.data['pow_i_aux'] = data[14*n:15*n,4]

        self.data['pow_e_fus']  = data[15*n:16*n,0]
        self.data['pow_i_fus']  = data[15*n:16*n,1]
        self.data['pow_e_sync'] = data[15*n:16*n,2]
        self.data['pow_e_brem'] = data[15*n:16*n,3]
        self.data['pow_e_line'] = data[15*n:16*n,4]

        # OPTIONAL: Read input.profiles.extra if it exists
        try:
            data = np.loadtxt(infile+'.extra',comments='#')
            print '(INFO): (profiles_gen_plot) '+infile+'.extra found.'
            x = data.reshape((n,45),order='F')

            self.data['bunit']     = x[0:n,0]
            self.data['s']         = x[0:n,1]
            self.data['drmaj']     = x[0:n,2]
            self.data['dzmag']     = x[0:n,3]
            self.data['sdelta']    = x[0:n,4]
            self.data['skappa']    = x[0:n,5]
            self.data['szeta']     = x[0:n,6]
            self.data['dlnnedr']   = x[0:n,7]
            self.data['dlntedr']   = x[0:n,8]
            self.data['dlnnidr_1'] = x[0:n,9]
            self.data['dlnnidr_2'] = x[0:n,10]
            self.data['dlnnidr_3'] = x[0:n,11]
            self.data['dlnnidr_4'] = x[0:n,12]
            self.data['dlnnidr_5'] = x[0:n,13]
            self.data['dlnnidr_6'] = x[0:n,14]
            self.data['dlnnidr_7'] = x[0:n,15]
            self.data['dlnnidr_8'] = x[0:n,16]
            self.data['dlnnidr_9'] = x[0:n,17]
            self.data['dlnnidr_10'] = x[0:n,18]
            self.data['dlntidr_1'] = x[0:n,19]
            self.data['dlntidr_2'] = x[0:n,20]
            self.data['dlntidr_3'] = x[0:n,21]
            self.data['dlntidr_4'] = x[0:n,22]
            self.data['dlntidr_5'] = x[0:n,23]
            self.data['dlntidr_6'] = x[0:n,24]
            self.data['dlntidr_7'] = x[0:n,25]
            self.data['dlntidr_8'] = x[0:n,26]
            self.data['dlntidr_9'] = x[0:n,27]
            self.data['dlntidr_10'] = x[0:n,28]
            self.data['dlnptotdr'] = x[0:n,29]
            self.data['drdrho']    = x[0:n,30]
            self.data['w0p']       = x[0:n,31]
            self.data['vol']       = x[0:n,32]
            self.data['volp']      = x[0:n,33]
            self.data['cs']        = x[0:n,34]
            self.data['rhos']      = x[0:n,35]
            self.data['ni_new']    = x[0:n,36]
            self.data['dlnnidr_new'] = x[0:n,37]
            self.data['grad_r0']   = x[0:n,38]
            self.data['ave_grad_r']= x[0:n,39]
            self.data['bp0']       = x[0:n,40]
            self.data['bt0']       = x[0:n,41]
            self.data['gamma_e']   = x[0:n,42]
            self.data['gamma_p']   = x[0:n,43]
            self.data['mach']      = x[0:n,44]
        except:
            print '(INFO): (profiles_gen_plot) '+infile+'.extra NOT found.'
          
        # OPTIONAL: Read input.profiles.geo if it exists
        try:
            # First, get number of Fourier modes
            fp = open(infile+'.geo')
            for i, line in enumerate(fp):
                if i == 11:
                    self.nfourier=int(line)
                    break
            fp.close()

            data = np.loadtxt(infile+'.geo',skiprows=12)
            # Dimension 9 assumes nfourier=8
            print '(INFO): (profiles_gen_plot) '+infile+'.geo found.'
            x = data.reshape((4,self.nfourier+1,n),order='F')
            self.geo['ar']=x[0,:,:]
            self.geo['br']=x[1,:,:]
            self.geo['az']=x[2,:,:]
            self.geo['bz']=x[3,:,:]
        except:
            print '(INFO): (profiles_gen_plot) '+infile+'.geo NOT found.'

        # OPTIONAL: Read input.profiles.jbs if it exists
        try:
            data = np.loadtxt(infile+'.jbs',comments='#')
            print '(INFO): (profiles_gen_plot) '+infile+'.jbs found.'
            x = data.reshape((n,6),order='F')

            self.data['jbs_err']    = x[0:n,1]
            self.data['jbs_neo']    = x[0:n,2]
            self.data['jbs_sauter'] = x[0:n,3]
            self.data['jbs_nclass'] = x[0:n,4]
            self.data['jbs_koh']    = x[0:n,5]
        except:
            print '(INFO): (profiles_gen_plot) '+infile+'.jbs NOT found.'
