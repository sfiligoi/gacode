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

        data = np.loadtxt(infile,skiprows=row)

        n = self.n_exp

        self.data['rho']   = data[0:n,0]
        self.data['rmin']  = data[0:n,1]
        self.data['rmaj']  = data[0:n,2]
        self.data['q']     = data[0:n,3]
        self.data['kappa'] = data[0:n,4]

        self.data['delta']  = data[n:2*n,0]
        self.data['Te']     = data[n:2*n,1]
        self.data['ne']     = data[n:2*n,2]
        self.data['z_eff']  = data[n:2*n,3]
        self.data['omega0'] = data[n:2*n,4]

        self.data['flow_mom'] = data[2*n:3*n,0]
        self.data['pow_e']    = data[2*n:3*n,1]
        self.data['pow_i']    = data[2*n:3*n,2]
        self.data['pow_ei']   = data[2*n:3*n,3]
        self.data['zeta']     = data[2*n:3*n,4]

        self.data['flow_beam'] = data[3*n:4*n,0]
        self.data['flow_wall'] = data[3*n:4*n,1]
        self.data['zmag']      = data[3*n:4*n,2]
        self.data['ptot']      = data[3*n:4*n,3]
        self.data['polflux']   = data[3*n:4*n,4]

        self.data['ni_1'] = data[4*n:5*n,0]
        self.data['ni_2'] = data[4*n:5*n,1]
        self.data['ni_3'] = data[4*n:5*n,2]
        self.data['ni_4'] = data[4*n:5*n,3]
        self.data['ni_5'] = data[4*n:5*n,4]

        self.data['Ti_1'] = data[5*n:6*n,0]
        self.data['Ti_2'] = data[5*n:6*n,1]
        self.data['Ti_3'] = data[5*n:6*n,2]
        self.data['Ti_4'] = data[5*n:6*n,3]
        self.data['Ti_5'] = data[5*n:6*n,4]

        self.data['vtor_1'] = data[6*n:7*n,0]
        self.data['vtor_2'] = data[6*n:7*n,1]
        self.data['vtor_3'] = data[6*n:7*n,2]
        self.data['vtor_4'] = data[6*n:7*n,3]
        self.data['vtor_5'] = data[6*n:7*n,4]

        self.data['vpol_1'] = data[7*n:8*n,0]
        self.data['vpol_2'] = data[7*n:8*n,1]
        self.data['vpol_3'] = data[7*n:8*n,2]
        self.data['vpol_4'] = data[7*n:8*n,3]
        self.data['vpol_5'] = data[7*n:8*n,4]

        data = np.loadtxt(infile+'.extra')
 
        try:
            data = np.loadtxt(infile+'.extra',comments='#')
            x = data.reshape((n,35),order='F')

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
            self.data['dlntidr_5'] = x[0:n,13]
            self.data['dlntidr_1'] = x[0:n,14]
            self.data['dlntidr_2'] = x[0:n,15]
            self.data['dlntidr_3'] = x[0:n,16]
            self.data['dlntidr_4'] = x[0:n,17]
            self.data['dlntidr_5'] = x[0:n,18]
            self.data['dlnptotdr'] = x[0:n,19]
            self.data['drdrho']    = x[0:n,20]
            self.data['w0p']       = x[0:n,21]
            self.data['vol']       = x[0:n,22]
            self.data['volp']      = x[0:n,23]
            self.data['cs']        = x[0:n,24]
            self.data['rhos']      = x[0:n,25]
            self.data['ni_new']    = x[0:n,26]
            self.data['dlnnidr_new'] = x[0:n,27]
            self.data['grad_r0']   = x[0:n,28]
            self.data['ave_grad_r']= x[0:n,29]
            self.data['bp0']       = x[0:n,30]
            self.data['bt0']       = x[0:n,31]
            self.data['gamma_e']   = x[0:n,32]
            self.data['gamma_p']   = x[0:n,33]
            self.data['mach']      = x[0:n,34]
        except:
            print infile+'.extra not available.'


        try:
            data = np.loadtxt(infile+'.geo',skiprows=12)
            # Dimension 9 assumes nfourier=8
            x = data.reshape((4,9,n),order='F')
            self.geo['ar']=x[0,:,:]
            self.geo['br']=x[1,:,:]
            self.geo['az']=x[2,:,:]
            self.geo['bz']=x[3,:,:]
        except:
            print infile+'.geo not available.'
