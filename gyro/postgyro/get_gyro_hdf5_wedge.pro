FUNCTION get_gyro_HDF5_wedge, data, N_TIME=n_time, TIME_SKIP=tskip
;
; C. Holland, UCSD
; v1.0: 3.26.2012
;
; routine to load in wedge HDF5 files in directory specified by
; data.simdir
;

  DEFAULT, n_time, 41
  DEFAULT, tskip, 50
  dirpath  =GETENV('GYRO_DIR') + '/sim/' + data.simdir

  FOR it = 0, n_time-1 DO BEGIN
      tlabel = STRCOMPRESS(STRING(it*tskip),/REMOVE_ALL)
      WHILE (STRLEN(tlabel) LT 5) DO tlabel = '0' + tlabel 
      print, tlabel
      wedge_label= 'gyrowedge'+tlabel+'.h5'
      wedge = H5_PARSE(dirpath+'/'+wedge_label, /READ)
      IF (it EQ 0) THEN BEGIN
          R_phys = wedge.r._data ;in m
          R_gyro = wedge.rgyro._data ;R/a
          Z_phys = wedge.z._data ;in m
          Z_gyro = wedge.zgyro._data ;z/a
          r_min = wedge.r_min._data
          theta = wedge
          alpha = wedge.alpha._data ;[n_tor_frac, n_rmin, n_y)
          n_tor_frac = wedge.alpha._dimensions[0]
          n_r = wedge.alpha._dimensions[1]
          n_y = wedge.alpha._dimensions[2]

          phi_RZ = 0
          n_RZ = 0
          E_RZ = 0
          T_RZ = 0
          V_RZ = 0

          IF (data.exists_n) THEN n_RZ = FLTARR(data.n_kinetic, n_tor_frac, n_r, n_y, n_time)
;          IF (data.exists_E) THEN E_RZ = FLTARR(data.n_kinetic, n_tor_frac, n_r, n_y, n_time)
;          IF (data.exists_T) THEN T_RZ = FLTARR(data.n_kinetic, n_tor_frac, n_r, n_y, n_time)
;          IF (data.exists_V) THEN V_RZ = FLTARR(data.n_kinetic, n_tor_frac, n_r, n_y, n_time)
      ENDIF

      ;load ion fluctuations
      IY = FLTARR(n_y) + 1
      FOR i_kin = 0, data.n_ion-1 DO BEGIN
          i_kin_str = NUMTOSTRING(i_kin)
 
         com = 'n_RZ['+i_kin_str+',*,*,*,'+NUMTOSTRING(it) + $
                '] = wedge.density_ion'+ i_kin_str + '_toroidal' + $
                '.density_ion' + i_kin_str+ '._data'
          print, com
          IF (data.exists_n) THEN r = EXECUTE(com)

         com = 'E_RZ['+i_kin_str+',*,*,*,'+NUMTOSTRING(it) + $
                '] = wedge.energy_ion'+ i_kin_str + '_toroidal' + $
                '.energy_ion' + i_kin_str+ '._data'
          print, com
          IF (data.exists_E) THEN r = EXECUTE(com)

         com = 'V_RZ['+i_kin_str+',*,*,*,'+NUMTOSTRING(it) + $
                '] = wedge.v_par_ion'+ i_kin_str + '_toroidal' + $
                '.v_par_ion' + i_kin_str+ '._data'
          print, com
          IF (data.exists_V) THEN r = EXECUTE(com)

          IF (data.exists_T) THEN FOR i_tor_frac = 0, n_tor_frac-1 DO $
            T_RZ[i_kin, i_tor_frac,*,*,it] = $
              ((2./3)*REFORM(E_RZ[i_kin,i_tor_frac,*,*,it]) - $
               (REFORM(data.T_eq[i_kin,*]) # IY)*REFORM(n_RZ[i_kin,i_tor_frac,*,*,it]))/$
              (REFORM(data.n_eq[i_kin,*]) # IY)
      ENDFOR

      IF (data.n_kinetic GT data.n_ion) THEN BEGIN ;load electron fluctuations
          i_kin = data.n_kinetic-1
          i_kin_str = NUMTOSTRING(i_kin)
 
         com = 'n_RZ['+i_kin_str+',*,*,*,'+NUMTOSTRING(it) + $
                '] = wedge.density_electron_toroidal' + $
                '.density_electron._data'
          print, com
          IF (data.exists_n) THEN r = EXECUTE(com)

         com = 'E_RZ['+i_kin_str+',*,*,*,'+NUMTOSTRING(it) + $
                '] = wedge.energy_electron_toroidal' + $
                '.energy_electron._data'
          print, com
          IF (data.exists_E) THEN r = EXECUTE(com)

         com = 'V_RZ['+i_kin_str+',*,*,*,'+NUMTOSTRING(it) + $
                '] = wedge.v_par_electron_toroidal' + $
                '.v_par_electron._data'
          print, com
          IF (data.exists_V) THEN r = EXECUTE(com)
 
          IF (data.exists_T) THEN FOR i_tor_frac = 0, n_tor_frac-1 DO $
            T_RZ[i_kin, i_tor_frac,*,*,it] = $
              ((2./3)*REFORM(E_RZ[i_kin,i_tor_frac,*,*,it]) - $
               (REFORM(data.T_eq[i_kin,*]) # IY)*REFORM(n_RZ[i_kin,i_tor_frac,*,*,it]))/$
              (REFORM(data.n_eq[i_kin,*]) # IY)
     ENDIF

  ENDFOR

  result = {simdir: data.simdir, $ ;#simulation directory
            n_time: n_time, $ ;# time points
            time_skip: tskip, $ ;# dt
            n_tor_frac: n_tor_frac, $ ;#torodial grid points for (R,Z) outputs
            n_r: n_r, $ ;#radial grid points
            n_y: n_y, $ ;#polodial grid points
            R_GYRO: R_GYRO, $ ;R/a
            Z_GYRO: Z_GYRO, $ ;Z/a
            R_phys: R_phys, $ ;m
            Z_phys: Z_phys, $, ;m
            r_min: r_min, $ ;r_min/a
            theta: theta, $ ;radians
            alpha: alpha, $ ;[n_tor_frac x n_r x_ny]
            exists_phi: 0, $ ;data.exists_phi, $
            exists_n: data.exists_n, $
            exists_E: data.exists_E, $
            exists_T: data.exists_T, $
            exists_V: data.exists_V, $
            phi_RZ: phi_RZ, $
            n_RZ: n_RZ, $
            E_RZ: E_RZ, $
            T_RZ: T_RZ, $
            V_RZ: V_RZ}
  RETURN, result
END ;get_gyro_HDF5_wedge
