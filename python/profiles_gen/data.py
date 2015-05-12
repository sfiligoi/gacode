vars_input_profiles = [
    ['rho', 'rmin', 'polflux', 'q', 'omega0'],
    ['rmaj', 'zmag', 'kappa', 'delta', 'zeta'],
    ['ne', 'Te', 'ptot', 'z_eff', 'NULL'],
    ['ni_1', 'ni_2', 'ni_3', 'ni_4', 'ni_5'],
    ['ni_6', 'ni_7', 'ni_8', 'ni_9', 'ni_10'],
    ['Ti_1', 'Ti_2', 'Ti_3', 'Ti_4', 'Ti_5'],
    ['Ti_6', 'Ti_7', 'Ti_8', 'Ti_9', 'Ti_10'],
    ['vtor_1', 'vtor_2', 'vtor_3', 'vtor_4', 'vtor_5'],
    ['vtor_6', 'vtor_7', 'vtor_8', 'vtor_9', 'vtor_10'],
    ['vpol_1', 'vpol_2', 'vpol_3', 'vpol_4', 'vpol_5'],
    ['vpol_6', 'vpol_7', 'vpol_8', 'vpol_9', 'vpol_10'],
    ['flow_beam', 'flow_wall', 'flow_mom', 'NULL', 'NULL'],
    ['pow_e', 'pow_i', 'pow_ei', 'pow_e_aux', 'pow_i_aux'],
    ['pow_e_fus', 'pow_i_fus', 'pow_e_sync', 'pow_e_brem', 'pow_e_line']
]

vars_input_profiles_extra = [
    'bunit', 's', 'drmaj', 'dzmag', 'sdelta', 'skappa', 'szeta', 'dlnnedr', 'dlntedr',
    'dlnnidr_1', 'dlnnidr_2', 'dlnnidr_3', 'dlnnidr_4', 'dlnnidr_5',
    'dlnnidr_6', 'dlnnidr_7', 'dlnnidr_8', 'dlnnidr_9', 'dlnnidr_10',
    'dlntidr_1', 'dlntidr_2', 'dlntidr_3', 'dlntidr_4', 'dlntidr_5',
    'dlntidr_6', 'dlntidr_7', 'dlntidr_8', 'dlntidr_9', 'dlntidr_10',
    'dlnptotdr', 'drdrho', 'w0p', 'vol', 'volp', 'cs', 'rhos', 'ni_new', 'dlnnidr_new',
    'grad_r0', 'ave_grad_r', 'bp0', 'bt0', 'gamma_e', 'gamma_p', 'mach'
]

vars_input_profiles_jbs = ['expro_rho', 'jbs_err', 'jbs_neo', 'jbs_sauter', 'jbs_nclass', 'jbs_koh']

#tuple: 1) latex 2) units 3) profiles column name
fancyNames = {'NULL'       : (''                ,''             ,'[null]'),
              'rho'        : ('\\rho'           ,''             ,'rho(-)'),
              'rmin'       : ('r_{min}'         ,'m'            ,'rmin(m)'),
              'rmaj'       : ('r_{maj}'         ,''             ,'rmax(m)'),
              'q'          : ('q'               ,''             ,'q(-)'),
              'kappa'      : ('\\kappa'         ,''             ,'kappa(-)'),
              'delta'      : ('\\delta'         ,''             ,'delta(-)'),
              'Te'         : ('T_e'             ,'keV'          ,'Te(keV)'),
              'ne'         : ('n_e'             ,'10^{19}/m^3'  ,'ne(10^19/m^3)'),
              'z_eff'      : ('z_{eff}'         ,''             ,'zeff(-)'),
              'omega0'     : ('\\omega_0'       ,'1/s'          ,'omega0(1/s)'),
              'flow_mom'   : ('flow_{mom}'      ,'Nm'           ,'flow_mom(Nm)'),
              'pow_e'      : ('p_e'             ,'MW'           ,'pow_e(MW)'),
              'pow_i'      : ('p_i'             ,'MW'           ,'pow_i(MW)'),
              'pow_ei'     : ('p_{ei}'          ,'MW'           ,'pow_ei(MW)'),
              'zeta'       : ('\\zeta'          ,''             ,'zeta(-)'),
              'flow_beam'  : ('flow_{beam}'     ,'kW/eV'        ,'flow_beam(kW/eV)'),
              'flow_wall'  : ('flow_{wall}'     ,'kW/eV'        ,'flow_wall(kW/eV)'),
              'zmag'       : ('z_{mag}'         ,'m'            ,'zmag(m)'),
              'ptot'       : ('p_{tot}'         ,'Pa'           ,'ptot(Pa)'),
              'polflux'    : ('\\psi'           ,'Wb/rad'       ,'polflux(Wb/rad)'),
              'pow_e_aux'  : ('p_{e,aux}'       ,'MW'           ,'pow_e_aux(MW)'),
              'pow_i_aux'  : ('p_{i,aux}'       ,'MW'           ,'pow_i_aux(MW)'),
              'pow_e_fus'  : ('p_{e,fus}'       ,'MW'           ,'pow_e_fus(MW)'),
              'pow_i_fus'  : ('p_{i,fus}'       ,'MW'           ,'pow_i_fus(MW)'),
              'pow_e_sync' : ('p_{e,sync}'      ,'MW'           ,'pow_e_sync(MW)'),
              'pow_e_brem' : ('p_{e,brem}'      ,'MW'           ,'pow_e_brem(MW)'),
              'pow_e_line' : ('p_{e,line}'      ,'MW'           ,'pow_e_line(MW)'),
              # extra (needs units)
              'bunit'      : ('bunit'           ,'T'            ,''),
              's'          : ('s'               ,''             ,''),
              'drmaj'      : ('drmaj'           ,''             ,''),
              'dzmag'      : ('dzmag'           ,''             ,''),
              'sdelta'     : ('sdelta'          ,''             ,''),
              'skappa'     : ('skappa'          ,''             ,''),
              'szeta'      : ('szeta'           ,''             ,''),
              'dlnnedr'    : ('dlnnedr'         ,'1/m'          ,''),
              'dlntedr'    : ('dlntedr'         ,'1/m'          ,''),
              'dlnptotdr'  : ('dlnptotdr'       ,'1/m'          ,''),
              'drdrho'     : ('drdrho'          ,'1/m'          ,''),
              'w0p'        : ('w0p'             ,'1/s/m'        ,''),
              'vol'        : ('vol'             ,'m^3'          ,''),
              'volp'       : ('volp'            ,'m^2'          ,''),
              'cs'         : ('cs'              ,'m/s'          ,''),
              'rhos'       : ('rhos'            ,'m'            ,''),
              'ni_new'     : ('ni_new'          ,'10^{19}/m^3'  ,''),
              'dlnnidr_new': ('dlnnidr_new'     ,'1/m'          ,''),
              'grad_r0'    : ('grad_r0'         ,''             ,''),
              'ave_grad_r' : ('ave_grad_r'      ,''             ,''),
              'bp0'        : ('bp0'             ,'T'            ,''),
              'bt0'        : ('bt0'             ,'T'            ,''),
              'gamma_e'    : ('gamma_e'         ,'1/s'          ,''),
              'gamma_p'    : ('gamma_p'         ,'1/s'          ,''),
              'mach'       : ('mach'            ,''             ,''),
              #jbs
              'expro_rho'  : ('rho'             ,''             ,''),
              'jbs_err'    : ('jbs_{err}'       ,'MA/m^2'       ,''),
              'jbs_neo'    : ('jbs_{neo}'       ,'MA/m^2'       ,''),
              'jbs_sauter' : ('jbs_{sauter}'    ,'MA/m^2'       ,''),
              'jbs_nclass' : ('jbs_{nclass}'    ,'MA/m^2'       ,''),
              'jbs_koh'    : ('jbs_{koh}'       ,'MA/m^2'       ,''),
              }
for _k in range(1, 11):
    fancyNames['ni_%d'%_k]     =('ni_%d'%_k     ,'10^19/m^3'    ,'ni_%d(10^19/m^3)'%_k)
    fancyNames['Ti_%d'%_k]     =('Ti_%d'%_k     ,'keV'          ,'Ti_%d(keV)'%_k)
    fancyNames['vtor_%d'%_k]   =('vtor_%d'%_k   ,'m/s'          ,'vtor_%d(m/s)')
    fancyNames['vpol_%d'%_k]   =('vpol_%d'%_k   ,'m/s'          ,'vpol_%d(m/s)')
    fancyNames['dlntedr_%d'%_k]=('dlntedr_%s'%_k,''             ,'')  # needs to be fixed
    fancyNames['dlntidr_%d'%_k]=('dlntidr_%s'%_k,''             ,'')  # needs to be fixed
    fancyNames['dlnnidr_%d'%_k]=('dlnnidr_%s'%_k,''             ,'')  # needs to be fixed

class profiles_genData:
    """
        profiles_gen output data class.
        """

    # ---------------------------------------------------------------------------#
    # Methods

    def __init__(self, infile):
        """
        Constructor reads in data from directory and creates new object.
        """

        import numpy as np

        # Initialize data
        self.infile = infile
        self.data = {}
        self.n_exp = 0
        self.geo = {}

        # Read input.profiles
        try:
            tmp = profiles_gen(infile)
            self.data.update(tmp.data)
            self.n_exp = tmp.n_exp
            print('(INFO): (profiles_genData) ' + infile + ' found.')
        except Exception as E:
            raise (IOError('(ERROR): (profiles_genData) ' + infile + ' not found: ' + str(E)))

        # OPTIONAL: Read input.profiles.extra if it exists
        try:
            self.data.update(profiles_gen_extra(infile + '.extra').data)
            print('(INFO): (profiles_genData) ' + infile + '.extra found.')
        except Exception as E:
            print('(INFO): (profiles_genData) ' + infile + '.extra NOT loaded: ' + str(E))

        # OPTIONAL: Read input.profiles.geo if it exists
        try:
            self.geo.update(profiles_gen_geo(infile + '.geo').geo)
            print('(INFO): (profiles_genData) ' + infile + '.geo found.')
        except Exception as E:
            print('(INFO): (profiles_genData) ' + infile + '.geo NOT loaded: ' + str(E))

        # OPTIONAL: Read input.profiles.jbs if it exists
        try:
            self.data.update(profiles_gen_jbs(infile + '.jbs').data)
            print('(INFO): (profiles_genData) ' + infile + '.jbs found.')
        except Exception as E:
            print('(INFO): (profiles_genData) ' + infile + '.jbs NOT loaded: ' + str(E))

class profiles_gen:
    def __init__(self, infile):
        import numpy as np

        self.data = {}
        self.n_exp = 0

        row = 0
        for line in open(infile, 'r').readlines():
            row = row + 1
            if line[0:5] == 'N_EXP':
                self.n_exp = int(line.split('=')[1])
            if line[0:8] == 'ARHO_EXP':
                self.max_rho = float(line.split('=')[1])
            if line[0:7] == '#rho(-)':
                break

        data = np.loadtxt(infile, skiprows=row)

        for k1, line in enumerate(vars_input_profiles):
            if self.n_exp * (k1 + 1) > data.shape[0]:
                break
            for k2, var in enumerate(line):
                self.data[var] = data[self.n_exp * k1:self.n_exp * (k1 + 1), k2]
        if 'NULL' in self.data:
            del self.data['NULL']

class profiles_gen_extra:
    def __init__(self, infile):
        import numpy as np

        data = np.loadtxt(infile, comments='#')
        data = data.reshape((-1, len(vars_input_profiles_extra)), order='F')

        self.data = {}
        for k, var in enumerate(vars_input_profiles_extra):
            self.data[var] = data[:, k]

class profiles_gen_geo:
    def __init__(self, infile):
        import numpy as np

        # First, get number of Fourier modes
        fp = open(infile)
        for i, line in enumerate(fp):
            if i == 11:
                self.nfourier = int(line)
                break
        fp.close()

        data = np.loadtxt(infile, skiprows=12)

        # Dimension 9 assumes nfourier=8
        self.geo = {}
        data = data.reshape((4, self.nfourier + 1, -1), order='F')
        self.geo['ar'] = data[0, :, :]
        self.geo['br'] = data[1, :, :]
        self.geo['az'] = data[2, :, :]
        self.geo['bz'] = data[3, :, :]

class profiles_gen_jbs:
    def __init__(self, infile):
        import numpy as np

        data = np.loadtxt(infile, comments='#')
        data = data.reshape((-1, len(vars_input_profiles_jbs)), order='F')

        self.data = {}
        for k, var in enumerate(vars_input_profiles_jbs):
            self.data[var] = data[:, k]

if __name__ == '__main__':
    # check that all fancy names have been assigned
    import numpy

    tmp = set(numpy.array(vars_input_profiles).flatten())
    tmp = tmp.union(set(vars_input_profiles_jbs))
    tmp = tmp.union(set(vars_input_profiles_extra))
    print tmp.difference(set(fancyNames.keys()))
    
    profiles_genData('/Users/meneghini/Dropbox/pycode/OMFIT-source/samples/input.profiles')
