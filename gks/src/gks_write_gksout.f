      subroutine write_gksout
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c     This routine writes out GKS results to a NetCDF 
c     file 'gksout.nc'
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      use gks_var
      implicit none
c
      include 'netcdf.inc'
      include 'input.m'
      include 'data_exp.m'
      include 'glf.m'
      include 'gks_out.m'
c
      character*(*) time_units, temp_units, leng_units, den_units, 
     &              field_units, cur_units, velocity_units
      integer txlen
      parameter (time_units='secs',temp_units='keV',leng_units='m')
      parameter (velocity_units='km/sec')
      parameter (den_units='1.e19 m^-3')
      parameter (field_units='T', cur_units='MA')
      parameter (txlen=15)
      integer i,j, k, val, lenid, ncid, rcode
      integer rhidd, rhid (200), rhid1 (200), chid(4)
      integer strval(4),scalval(200), oneval(200)
      integer stat, rdim
      integer timeid, txid
      integer txdims, cstart, ccount
      real*8 relval
      character*40 txval
c
      txval(txlen:txlen) = char(0)  ! null terminate
c open netcdf output file
      val = nf_create('gksout.nc', nf_clobber, ncid)
c      if (val .ne. nf_noerr) call handle_err(val)
c
c Define scalars and attributes
c
      i=0
      j=0
      k=0
c  model and run selection
      i=i+1
      strval(i) = nf_def_dim(ncid,'l_version',len(version),txdims)
      strval(i) = nf_def_var(ncid,'version',nf_char,1,
     >txdims,chid(i))
      strval(i) = nf_put_att_text(ncid,chid(i),'definition',
     >len('TGLF version number'),'TGLF version number')
      j=j+1
      scalval(j) = nf_def_var(ncid,'igks_model',nf_int,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('GKS=0, TGLF=1'),'GKS=0, TGLF=1')
      j=j+1
      scalval(j) = nf_def_var(ncid,'igks',nf_int,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('run type'),'run type')
      j=j+1
      scalval(j) = nf_def_var(ncid,'jin_m',nf_int,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('inner gridpoint'),'inner gridpoint')
      j=j+1
      scalval(j) = nf_def_var(ncid,'jout_m',nf_int,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('outer gridpoint'),'outer gridpoint')
c general settings shared by both models
      j=j+1
      scalval(j) = nf_def_var(ncid,'gks_defaults',nf_int,0,0,
     >rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('default selection'),'default selection')
      j=j+1
      scalval(j) = nf_def_var(ncid,'kys0',nf_double,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('reference wavenumber in TGLF units'),
     >  'reference wavenumber in TGLF units')
      j=j+1
      scalval(j) = nf_def_var(ncid,'aky1',nf_double,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('reference wavenumber in GKS units'),
     >  'reference wavenumber in GKS units')
      j=j+1
      scalval(j) = nf_def_var(ncid,'igeo',nf_int,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('s-alpha=0,miller=1'),'s-alpha=0,miller=1')
      j=j+1
      scalval(j) = nf_def_var(ncid,'igyro_fix',nf_int,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('B in FLR: no=0, yes=1'),'B in FLR: no=0, yes=1')
      j=j+1
      scalval(j) = nf_def_var(ncid,'i_bpar',nf_int,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('use B_par: no=0, yes=1'),'use B_par: no=0, yes=1')
      j=j+1
      scalval(j) = nf_def_var(ncid,'cdebye',nf_double,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('debye term multiplier'),'debye term multiplier')
      j=j+1
      scalval(j) = nf_def_var(ncid,'cbetae',nf_double,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('betae multiplier'),'betae multiplier')
      j=j+1
      scalval(j) = nf_def_var(ncid,'xalpha',nf_double,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('alpha multiplier'),'alpha multiplier')
      j=j+1
      scalval(j) = nf_def_var(ncid,'cnewk3',nf_double,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('vnewk3 multiplier'),'vnewk3 multiplier')
      j=j+1
      scalval(j) = nf_def_var(ncid,'cnewk2',nf_double,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('vnewk2 multiplier'),'vnewk2 multiplier')
      j=j+1
      scalval(j) = nf_def_var(ncid,'alpha_mach',nf_double,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('parallel ion velocity multiplier'),
     >  'parallel ion velocity multiplier')
      j=j+1
      scalval(j) = nf_def_var(ncid,'alpha_cur',nf_double,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('parallel current velocity multiplier'),
     >  'parallel current velocity multiplier')
      j=j+1
      scalval(j) = nf_def_var(ncid,'alpha_p',nf_double,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('parallel ion velocity shear multiplier'),
     >  'parallel ion velocity shear multiplier')
      j=j+1
      scalval(j) = nf_def_var(ncid,'alpha_p_cur',nf_double,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('parallel electron velocity shear multiplier'),
     >  'parallel electron velocity shear multiplier')
      j=j+1
      scalval(j) = nf_def_var(ncid,'alpha_e',nf_double,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('exb velocity shear multiplier'),
     >  'exb velocity shear multiplier')
      j=j+1
      scalval(j) = nf_def_var(ncid,'z2',nf_double,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('impurity ion charge'),'impurity ion charge')
      j=j+1
      scalval(j) = nf_def_var(ncid,'z3',nf_double,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('electron charge'),'electron charge')
      j=j+1
      scalval(j) = nf_def_var(ncid,'amass2',nf_double,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('impurity ion mass'),'impurity ion mass')
      j=j+1
      scalval(j) = nf_def_var(ncid,'amass3',nf_double,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('electron mass'),'electron mass')
      j=j+1
      scalval(j) = nf_def_var(ncid,'ncspec1',nf_int,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('use ion 1: no=0, yes=1'),'use ion1: no=0, yes=1')
      j=j+1
      scalval(j) = nf_def_var(ncid,'ncspec2',nf_int,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('use ion 2: no=0, yes=1'),'use ion 2: no=0, yes=1')
      j=j+1
      scalval(j) = nf_def_var(ncid,'ncspec3',nf_int,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('use electrons: no=0, yes=1'),
     >  'use electrons: no=0, yes=1')
      j=j+1
      scalval(j) = nf_def_var(ncid,'ncspec4',nf_int,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('use ion 4: no=0, yes=1'),'use ion 4: no=0, yes=1')
      j=j+1
      scalval(j) = nf_def_var(ncid,'ncspec5',nf_int,0,0,rhid(j))
      scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('use ion 5: no=0, yes=1'),'use ion 5: no=0, yes=1')
      if(igks.ne.0)then
c data will be read echo file and data handling settings
        j=j+1
        scalval(j) = nf_def_var(ncid,'idata',nf_int,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('data type: TRANSP=0, ONETWO=1'),
     >  'data type: TRANSP=0, ONETWO=1')
        j=j+1
        scalval(j) = nf_def_var(ncid,'xp_time',nf_double,0,
     >  0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('analysis time'),'analysis time')
        i=i+1
        strval(i) = nf_def_dim(ncid,'l_tok',len(tok),txdims)
        strval(i) = nf_def_var(ncid,'tok',nf_char,1,
     >  txdims,chid(i))
        strval(i) = nf_put_att_text(ncid,chid(i),'definition',
     >  len('tokamak id'),'tokamak id')
        i=i+1
        strval(i) = nf_def_dim(ncid,'l_shot',len(shot),txdims)
        strval(j) = nf_def_var(ncid,'shot',nf_char,1,
     >  txdims,chid(i))
        strval(i) = nf_put_att_text(ncid,chid(i),'definition',
     >  len('discharge number'),'discharge number')
        i=i+1
        strval(i) = nf_def_dim(ncid,'l_cudir',len(cudir),txdims)
        strval(i) = nf_def_var(ncid,'cudir',nf_char,1,
     >  txdims,chid(i))
        strval(i) = nf_put_att_text(ncid,chid(i),'definition',
     >  len('path to data'),'path to data')
        j=j+1
        scalval(j) = nf_def_var(ncid,'iptot',nf_int,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('read total pressure: no=0, yes=1'),
     >  'read total pressure: no=0, yes=1')
        j=j+1
        scalval(j) = nf_def_var(ncid,'ismooth_data',nf_int,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('smooth data: no=0, yes=1'),'smooth data: no=0, yes=1')
        j=j+1
        scalval(j) = nf_def_var(ncid,'corot',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('toroidal rotation multiplier'),
     >  'toroidal rotation multiplier')
        j=j+1
        scalval(j) = nf_def_var(ncid,'irot1',nf_int,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('toroidal rotation option 1'),'toroidal rotation option 1')
        j=j+1
        scalval(j) = nf_def_var(ncid,'irot2',nf_int,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('toroidal rotation option 2'),'toroidal rotation option 2')
        j=j+1
        scalval(j) = nf_def_var(ncid,'arho_exp',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('unit lenght'),'unit length')
        j=j+1
        scalval(j) = nf_def_var(ncid,'Bt_exp',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('magnetic field in definition of rho'),
     >  'magnetic field in definition of rho')
      endif
      if(igks_model.eq.0)then
c GKS specific settings
        j=j+1
        scalval(j) = nf_def_var(ncid,'nstep',nf_int,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('number of GKS timesteps'),'number of GKS timesteps')
        j=j+1
        scalval(j) = nf_def_var(ncid,'ntheta',nf_int,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('number of GKS grid pts/period'),
     >  'number of GKS grid pts/period')
        j=j+1
        scalval(j) = nf_def_var(ncid,'nperiod',nf_int,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('number of GKS periods'),'number of GKS periods')
        j=j+1
        scalval(j) = nf_def_var(ncid,'ngauss',nf_int,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('maximum number of GKS energy pts'),
     >  'maximum of GKS energy pts')
        j=j+1
        scalval(j) = nf_def_var(ncid,'negrid',nf_int,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('number of GKS energy pts'),'number of GKS energy pts')
        j=j+1
        scalval(j) = nf_def_var(ncid,'ipar',nf_int,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('symmetric GKS wavefunction: no=0, yes=1'),
     >  'symmetric GKS wavefunction: no=0, yes=1')
        j=j+1
        scalval(j) = nf_def_var(ncid,'delt',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('GKS timestep scale'),'GKS timestep scale')
        j=j+1
        scalval(j) = nf_def_var(ncid,'tol',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('GKS error tolerance'),'GKS error tolerance')
      elseif(igks_model.eq.1)then
c TGLF specific settings
        j=j+1
        scalval(j) = nf_def_var(ncid,'nbasis_min',nf_int,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('TGLF nbasis_min'),'TGLF nbasis_min')
        j=j+1
        scalval(j) = nf_def_var(ncid,'nbasis_max',nf_int,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('TGLF nbasis_max'),'TGLF nbasis_max')
      endif  
c output for specific run types igks
      if(igks.lt.3.or.igks.eq.7.or.(jin_m.eq.jout_m))then
        j=j+1
        scalval(j) = nf_def_var(ncid,'rho_k',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('rho position'),'rho position')
        j=j+1
        scalval(j) = nf_def_var(ncid,'tprim1',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('tprim1'),'tprim1')
        j=j+1
        scalval(j) = nf_def_var(ncid,'tprim2',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('tprim2'),'tprim2')
        j=j+1
        scalval(j) = nf_def_var(ncid,'tprim3',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('tprim3'),'tprim3')
        j=j+1
        scalval(j) = nf_def_var(ncid,'fprim1',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('fprim1'),'fprim1')
        j=j+1
        scalval(j) = nf_def_var(ncid,'fprim2',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('fprim2'),'fprim2')
        j=j+1
        scalval(j) = nf_def_var(ncid,'fprim3',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('fprim3'),'fprim3')
        j=j+1
        scalval(j) = nf_def_var(ncid,'mach1',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('mach1'),'mach1')
        j=j+1
        scalval(j) = nf_def_var(ncid,'mach2',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('mach2'),'mach2')
        j=j+1
        scalval(j) = nf_def_var(ncid,'mach3',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('mach3'),'mach3')
        j=j+1
        scalval(j) = nf_def_var(ncid,'uprim1',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('uprim1'),'uprim1')
        j=j+1
        scalval(j) = nf_def_var(ncid,'uprim2',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('uprim2'),'uprim2')
        j=j+1
        scalval(j) = nf_def_var(ncid,'uprim3',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('uprim3'),'uprim3')
        j=j+1
        scalval(j) = nf_def_var(ncid,'beta',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('beta'),'beta')
        j=j+1
        scalval(j) = nf_def_var(ncid,'zeff',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('zeff'),'zeff')
        j=j+1
        scalval(j) = nf_def_var(ncid,'temp3',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('temp3'),'temp3')
        j=j+1
        scalval(j) = nf_def_var(ncid,'vnewk1',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('vnewk1'),'vnewk1')
        j=j+1
        scalval(j) = nf_def_var(ncid,'vnewk2',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('vnewk2'),'vnewk2')
        j=j+1
        scalval(j) = nf_def_var(ncid,'vnewk3',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('vnewk3'),'vnewk3')
        if(igeo.eq.0)then
          j=j+1
          scalval(j) = nf_def_var(ncid,'dbeam',nf_double,0,0,rhid(j))
          scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >    len('dbeam'),'dbeam')
          j=j+1
          scalval(j) = nf_def_var(ncid,'teti',nf_double,0,0,rhid(j))
          scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >    len('teti'),'teti')
          j=j+1
          scalval(j) = nf_def_var(ncid,'pk',nf_double,0,0,rhid(j))
          scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >    len('pk'),'pk')
          j=j+1
          scalval(j) = nf_def_var(ncid,'eps',nf_double,0,0,rhid(j))
          scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >    len('eps'),'eps')
          j=j+1
          scalval(j) = nf_def_var(ncid,'epsa',nf_double,0,0,rhid(j))
          scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >    len('epsa'),'epsa')
          j=j+1
          scalval(j) = nf_def_var(ncid,'epsl',nf_double,0,0,rhid(j))
          scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >    len('epsl'),'epsl')
          j=j+1
          scalval(j) = nf_def_var(ncid,'shat',nf_double,0,0,rhid(j))
          scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >    len('shat'),'shat')
          j=j+1
          scalval(j) = nf_def_var(ncid,'shift',nf_double,0,0,rhid(j))
          scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >    len('shift'),'shift')
        elseif(igeo.eq.1)then
          j=j+1
          scalval(j) = nf_def_var(ncid,'nione_loc',
     >    nf_double,0,0,rhid(j))
          scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >    len('nione_loc'),'nione_loc')
          j=j+1
          scalval(j) = nf_def_var(ncid,'fastionfrac_loc',
     >    nf_double,0,0,rhid(j))
          scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >    len('fastionfrac_loc'),'fastionfrac_loc')
          j=j+1
          scalval(j) = nf_def_var(ncid,'tiote_loc',
     >    nf_double,0,0,rhid(j))
          scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >    len('tiote_loc'),'tiote_loc')
          j=j+1
          scalval(j) = nf_def_var(ncid,'rmin_loc',
     >    nf_double,0,0,rhid(j))
          scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >    len('rmin_loc'),'rmin_loc')
          j=j+1
          scalval(j) = nf_def_var(ncid,'rmaj0_loc',
     >    nf_double,0,0,rhid(j))
          scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >    len('rmaj0_loc'),'rmaj0_loc')
          j=j+1
          scalval(j) = nf_def_var(ncid,'q_loc',
     >    nf_double,0,0,rhid(j))
          scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >    len('q_loc'),'q_loc')
          j=j+1
          scalval(j) = nf_def_var(ncid,'shat_mhd_loc',
     >    nf_double,0,0,rhid(j))
          scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >    len('shat_mhd_loc'),'shat_mhd_loc')
          j=j+1
          scalval(j) = nf_def_var(ncid,'alpha_mhd_loc',
     >    nf_double,0,0,rhid(j))
          scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >    len('alpha_mhd_loc'),'alpha_mhd_loc')
          j=j+1
          scalval(j) = nf_def_var(ncid,'shift_loc',
     >    nf_double,0,0,rhid(j))
          scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >    len('shift_loc'),'shift_loc')
          j=j+1
          scalval(j) = nf_def_var(ncid,'kappa_loc',
     >    nf_double,0,0,rhid(j))
          scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >    len('kappa_loc'),'kappa_loc')
          j=j+1
          scalval(j) = nf_def_var(ncid,'delta_loc',
     >    nf_double,0,0,rhid(j))
          scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >    len('delta_loc'),'delta_loc')
          j=j+1
          scalval(j) = nf_def_var(ncid,'s_kappa_loc',
     >    nf_double,0,0,rhid(j))
          scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >    len('s_kappa_loc'),'s_kappa_loc')
          j=j+1
          scalval(j) = nf_def_var(ncid,'s_delta_loc',
     >    nf_double,0,0,rhid(j))
          scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >    len('s_delta_loc'),'s_delta_loc')
        endif
      endif
      if(igks.eq.0.or.(igks.ge.3.and.(jin_m.eq.jout_m)))then
        j=j+1
        scalval(j) = nf_def_var(ncid,'agammas',
     >  nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('agammas'),'agammas')
        j=j+1
        scalval(j) = nf_def_var(ncid,'dgammas',
     >  nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('dgammas'),'dgammas')
        j=j+1
        scalval(j) = nf_def_var(ncid,'afreqs',
     >  nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('afreqs'),'afreqs')
        j=j+1
        scalval(j) = nf_def_var(ncid,'kys',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('kys'),'kys')
        j=j+1
        scalval(j) = nf_def_var(ncid,'peflxa',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('peflxa'),'peflxa')
        j=j+1
        scalval(j) = nf_def_var(ncid,'eeflxa',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('eeflxa'),'eeflxa')
        j=j+1
        scalval(j) = nf_def_var(ncid,'eiflxa',nf_double,0,0,rhid(j))
        scalval(j) = nf_put_att_text(ncid,rhid(j),'definition',
     >  len('eiflxa'),'eiflxa')
c
c define rdim
        stat = nf_def_dim(ncid, 'l_theta',ntgridl+ntgrid+1,rdim)
c          if (stat .ne. nf_noerr) call handle_err(stat)
        k=k+1
        oneval(k) = nf_def_var (ncid,'rho',
     >  nf_double,1,rdim,rhid1(k))
        oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >  len('theta'),'theta')
        k=k+1
        oneval(k) = nf_def_var (ncid,'rephi',
     >  nf_double,1,rdim,rhid1(k))
        oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >  len('dimensionless'),'dimensionless')
        k=k+1
        oneval(k) = nf_def_var (ncid,'imphi',
     >  nf_double,1,rdim,rhid1(k))
        oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >  len('dimensionless'),'dimensionless')
        if(cbetae.gt.0.0)then
           k=k+1
           oneval(k) = nf_def_var (ncid,'reapar',
     >     nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >      len('dimensionless'),'dimensionless')
            k=k+1
            oneval(k) = nf_def_var (ncid,'imapar',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >      len('dimensionless'),'dimensionless')
            if(i_bpar.ne.0)then
              k=k+1
              oneval(k) = nf_def_var (ncid,'rebpar',
     >        nf_double,1,rdim,rhid1(k))
              oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >        len('dimensionless'),'dimensionless')
              k=k+1
              oneval(k) = nf_def_var (ncid,'imbpar',
     >        nf_double,1,rdim,rhid1(k))
              oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >        len('dimensionless'),'dimensionless')
            endif
          endif
          if(igeo.eq.1)then
            k=k+1
            oneval(k) = nf_def_var (ncid,'sintheta_geo',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >      len('dimensionless'),'dimensionless')
            k=k+1
            oneval(k) = nf_def_var (ncid,'costheta_geo',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >      len('dimensionless'),'dimensionless')
            k=k+1
            oneval(k) = nf_def_var (ncid,'b_geo',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >      len('dimensionless'),'dimensionless')
          endif
        elseif(igks.eq.1.or.igks.eq.2.or.igks.eq.7)then  ! igks = 1,2 or 7 only
c define rdim
          stat = nf_def_dim(ncid,'l_ky',21,rdim)
c          if (stat .ne. nf_noerr) call handle_err(stat)
          k=k+1
          oneval(k) = nf_def_var (ncid,'rho',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >    len('aky'),'aky')
          k=k+1
          oneval(k) = nf_def_var (ncid,'kys_k',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >    len('kys'),'kys')
          k=k+1
          oneval(k) = nf_def_var (ncid,'gamma_k',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >    len('cs/a'),'cs/a')
          k=k+1
          oneval(k) = nf_def_var (ncid,'dgamma_k',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >    len('cs/a'),'cs/a')
          k=k+1
          oneval(k) = nf_def_var (ncid,'freq_k',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >    len('cs/a'),'cs/a')
          k=k+1
          oneval(k) = nf_def_var (ncid,'ky_mks',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >    len('1/m'),'1/m')
          k=k+1
          oneval(k) = nf_def_var (ncid,'gamma_mks',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >    len('rad/s'),'rad/s')
          k=k+1
          oneval(k) = nf_def_var (ncid,'dgamma_mks',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >    len('rad/s'),'rad/s')
          k=k+1
          oneval(k) = nf_def_var (ncid,'freq_mks',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >    len('rad/s'),'rad/s')
          k=k+1
          oneval(k) = nf_def_var (ncid,'peflx_m',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >    len('gyrobohm'),'gyrobohm')
          k=k+1
          oneval(k) = nf_def_var (ncid,'qeflx_m',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >    len('gyrobohm'),'gyrobohm')
          k=k+1
          oneval(k) = nf_def_var (ncid,'qiflx_m',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >    len('gyrobohm'),'gyrobohm')
          if(igks_model.eq.1)then ! TGLF
            k=k+1
            oneval(k) = nf_def_var (ncid,'gamma_ion',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >      len('cs/a'),'cs/a')
            k=k+1
            oneval(k) = nf_def_var (ncid,'freq_ion',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >      len('cs/a'),'cs/a')
            k=k+1
            oneval(k) = nf_def_var (ncid,'gamma_electron',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >      len('cs/a'),'cs/a')
            k=k+1
            oneval(k) = nf_def_var (ncid,'freq_electron',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >      len('cs/a'),'cs/a')
            k=k+1
            oneval(k) = nf_def_var (ncid,'phi_bar_m',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >      len('dimensionless'),'dimensionless')
            k=k+1
            oneval(k) = nf_def_var (ncid,'ne_bar_m',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >      len('dimensionless'),'dimensionless')
            k=k+1
            oneval(k) = nf_def_var (ncid,'te_bar_m',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >      len('dimensionless'),'dimensionless')
            k=k+1
            oneval(k) = nf_def_var (ncid,'ti_bar_m',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >      len('dimensionless'),'dimensionless')
            k=k+1
            oneval(k) = nf_def_var (ncid,'ne_te_phase_ion',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >      len('dimensionless'),'dimensionless')
            k=k+1
            oneval(k) = nf_def_var (ncid,
     >      'ne_te_phase_electron',nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >      len('dimensionless'),'dimensionless')
          endif
        endif
        if(igks.ge.3.and.(jin_m.ne.jout_m))then
c define rdim
          stat = nf_def_dim(ncid,'l_rho',jmaxm+1,rdim)
c          if (stat .ne. nf_noerr) call handle_err(stat)
          k=k+1
          oneval(k) = nf_def_var (ncid,'rho',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          k=k+1
          oneval(k) = nf_def_var (ncid,'ne_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len(den_units),den_units)
          k=k+1
          oneval(k) = nf_def_var (ncid,'ni_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len(den_units),den_units)
          k=k+1
          oneval(k) = nf_def_var (ncid,'nz_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len(den_units),den_units)
          k=k+1
          oneval(k) = nf_def_var (ncid,'nfst_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len(den_units),den_units)
          k=k+1
          oneval(k) = nf_def_var (ncid,'te_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len(temp_units),temp_units)
          k=k+1
          oneval(k) = nf_def_var (ncid,'ti_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len(temp_units),temp_units)
          k=k+1
          oneval(k) = nf_def_var (ncid,'ptot_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('total pressure'),'total pressure')
          k=k+1
          oneval(k) = nf_def_var (ncid,'pfast_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('fast ion pressure'),'fast ion pressure')
          k=k+1
          oneval(k) = nf_def_var (ncid,'vphi_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len(velocity_units),velocity_units)
          k=k+1
          oneval(k) = nf_def_var (ncid,'vecur_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len(velocity_units),velocity_units)
          k=k+1
          oneval(k) = nf_def_var (ncid,'zeff_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          k=k+1
          oneval(k) = nf_def_var (ncid,'betae_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          k=k+1
          oneval(k) = nf_def_var (ncid,'q_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          k=k+1
          oneval(k) = nf_def_var (ncid,'shat_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          k=k+1
          oneval(k) = nf_def_var (ncid,'alpha_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          k=k+1
          oneval(k) = nf_def_var (ncid,'delta_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          k=k+1
          oneval(k) = nf_def_var (ncid,'elong_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          k=k+1
          oneval(k) = nf_def_var (ncid,'drhodr',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          k=k+1
          oneval(k) = nf_def_var (ncid,'zpne_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          k=k+1
          oneval(k) = nf_def_var (ncid,'zpni_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          k=k+1
          oneval(k) = nf_def_var (ncid,'zpnz_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          k=k+1
          oneval(k) = nf_def_var (ncid,'zpte_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          k=k+1
          oneval(k) = nf_def_var (ncid,'zpti_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          k=k+1
          oneval(k) = nf_def_var (ncid,'mach_i_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          k=k+1
          oneval(k) = nf_def_var (ncid,'mach_e_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          k=k+1
          oneval(k) = nf_def_var (ncid,'gamma_p_i_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          k=k+1
          oneval(k) = nf_def_var (ncid,'gamma_p_e_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          k=k+1
          oneval(k) = nf_def_var (ncid,'egamma_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          k=k+1
          oneval(k) = nf_def_var (ncid,'csda_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          k=k+1
          oneval(k) = nf_def_var (ncid,'rhosda_exp',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          k=k+1
          oneval(k) = nf_def_var (ncid,'ky_m',
     >    nf_double,1,rdim,rhid1(k))
          oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
          if(igks.le.4)then
            k=k+1
            oneval(k) = nf_def_var (ncid,'anrate_m',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('cs/a'),'cs/a')
            k=k+1
            oneval(k) = nf_def_var (ncid,'dnrate_m',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('cs/a'),'cs/a')
            k=k+1
            oneval(k) = nf_def_var (ncid,'anfreq_m',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('cs/a'),'cs/a')
            k=k+1
            oneval(k) = nf_def_var (ncid,'ky_mks',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('1/m'),'1/m')
            k=k+1
            oneval(k) = nf_def_var (ncid,'gamma_mks',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('rad/s'),'rad/s')
            k=k+1
            oneval(k) = nf_def_var (ncid,'dgamma_mks',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('rad/s'),'rad/s')
            k=k+1
            oneval(k) = nf_def_var (ncid,'freq_mks',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('rad/s'),'rad/s')
            if(igks.eq.3.and.igks_model.eq.1)then
              k=k+1
              oneval(k) = nf_def_var (ncid,'gamma_ion',
     >        nf_double,1,rdim,rhid1(k))
              oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('cs/a'),'cs/a')
              k=k+1
              oneval(k) = nf_def_var (ncid,'freq_ion',
     >        nf_double,1,rdim,rhid1(k))
              oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('cs/a'),'cs/a')
              k=k+1
              oneval(k) = nf_def_var (ncid,'gamma_electron',
     >        nf_double,1,rdim,rhid1(k))
              oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('cs/a'),'cs/a')
              k=k+1
              oneval(k) = nf_def_var (ncid,'freq_electron',
     >        nf_double,1,rdim,rhid1(k))
              oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('cs/a'),'cs/a')
              k=k+1
              oneval(k) = nf_def_var (ncid,'phi_bar_m',
     >        nf_double,1,rdim,rhid1(k))
              oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
              k=k+1
              oneval(k) = nf_def_var (ncid,'ne_bar_m',
     >        nf_double,1,rdim,rhid1(k))
              oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
              k=k+1
              oneval(k) = nf_def_var (ncid,'te_bar_m',
     >        nf_double,1,rdim,rhid1(k))
              oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
              k=k+1
              oneval(k) = nf_def_var (ncid,'ti_bar_m',
     >        nf_double,1,rdim,rhid1(k))
              oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
              k=k+1
              oneval(k) = nf_def_var (ncid,'ne_te_phase_ion',
     >        nf_double,1,rdim,rhid1(k))
              oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
              k=k+1
              oneval(k) = nf_def_var (ncid,
     >        'ne_te_phase_electron',nf_double,1,rdim,rhid1(k))
              oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
              k=k+1
              oneval(k) = nf_def_var (ncid,'mhd_DR_m',
     >        nf_double,1,rdim,rhid1(k))
              oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
            endif
          elseif(igks.eq.5)then
            k=k+1
            oneval(k) = nf_def_var (ncid,'zpte_crit_m',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
            k=k+1
            oneval(k) = nf_def_var (ncid,'ky_mks',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('1/m'),'1/m')
          elseif(igks.eq.6)then
            k=k+1
            oneval(k) = nf_def_var (ncid,'zpti_crit_m',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('dimensionless'),'dimensionless')
            k=k+1
            oneval(k) = nf_def_var (ncid,'anfreq_m',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('cs/a'),'cs/a')
            k=k+1
            oneval(k) = nf_def_var (ncid,'ky_mks',
     >      nf_double,1,rdim,rhid1(k))
            oneval(k)=nf_put_att_text(ncid,rhid1(k),'units',
     >          len('1/m'),'1/m')
          endif
        endif
c
c end definitions
c
        val = nf_enddef(ncid)
c
c
c write data
c
      i=0
      j=0
      k=0
c
      j=j+1
c  model and run selection
      i=i+1
      strval(i) = nf_put_var_text(ncid, chid(i),version )
      scalval(j) = nf_put_var_int(ncid, rhid(j),igks_model )
      j=j+1
      scalval(j) = nf_put_var_int(ncid, rhid(j),igks )
      j=j+1
      scalval(j) = nf_put_var_int(ncid, rhid(j),jin_m )
      j=j+1
      scalval(j) = nf_put_var_int(ncid, rhid(j),jout_m )
c general settings shared by both models
      j=j+1
      scalval(j) = nf_put_var_int(ncid, rhid(j),gks_defaults )
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j),kys0 )
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j),aky1 )
      j=j+1
      scalval(j) = nf_put_var_int(ncid, rhid(j),igeo )
      j=j+1
      scalval(j) = nf_put_var_int(ncid, rhid(j),igyro_fix)
      j=j+1
      scalval(j) = nf_put_var_int(ncid, rhid(j),i_bpar )
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j),cdebye )
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j),cbetae )
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j),xalpha )
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j),cnewk3 )
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j),cnewk2 )
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j),alpha_mach)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j),alpha_cur)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j),alpha_p)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j),alpha_p_cur)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j),alpha_e)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j),z2)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j),z3)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j),amass2)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j),amass3)
      j=j+1
      scalval(j) = nf_put_var_int(ncid, rhid(j),ncspec1)
      j=j+1
      scalval(j) = nf_put_var_int(ncid, rhid(j),ncspec2)
      j=j+1
      scalval(j) = nf_put_var_int(ncid, rhid(j),ncspec3)
      j=j+1
      scalval(j) = nf_put_var_int(ncid, rhid(j),ncspec4)
      j=j+1
      scalval(j) = nf_put_var_int(ncid, rhid(j),ncspec5)
      if(igks.ne.0)then
c data will be read echo file and data handling settings
        j=j+1
        scalval(j) = nf_put_var_int(ncid, rhid(j),idata )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),xp_time)
        i=i+1
        strval(i) = nf_put_var_text(ncid, chid(i),tok )
        i=i+1
        strval(i) = nf_put_var_text(ncid, chid(i),shot )
        i=i+1
        strval(i) = nf_put_var_text(ncid, chid(i),cudir)
        j=j+1
        scalval(j) = nf_put_var_int(ncid, rhid(j),iptot )
        j=j+1
        scalval(j) = nf_put_var_int(ncid, rhid(j),ismooth_data)
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),corot )
        j=j+1
        scalval(j) = nf_put_var_int(ncid, rhid(j),irot1 )
        j=j+1
        scalval(j) = nf_put_var_int(ncid, rhid(j),irot2 )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),arho_exp )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),Bt_exp )
      endif
      if(igks_model.eq.0)then
c GKS specific settings
        j=j+1
        scalval(j) = nf_put_var_int(ncid, rhid(j),nstep )
        j=j+1
        scalval(j) = nf_put_var_int(ncid, rhid(j),ntheta )
        j=j+1
        scalval(j) = nf_put_var_int(ncid, rhid(j),nperiod )
        j=j+1
        scalval(j) = nf_put_var_int(ncid, rhid(j),ngauss )
        j=j+1
        scalval(j) = nf_put_var_int(ncid, rhid(j),negrid )
        j=j+1
        scalval(j) = nf_put_var_int(ncid, rhid(j),ipar )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),delt )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),tol )
      elseif(igks_model.eq.1)then
c TGLF specific settings
        j=j+1
        scalval(j) = nf_put_var_int(ncid, rhid(j),nbasis_min )
        j=j+1
        scalval(j) = nf_put_var_int(ncid, rhid(j),nbasis_max )
      endif  
c output for specific run types igks
      if(igks.lt.3.or.igks.eq.7.or.(jin_m.eq.jout_m))then
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),rho(jin_m) )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),tprim1 )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),tprim2 )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),tprim3 )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),fprim1 )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),fprim2 )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),fprim3 )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),mach1 )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),mach2 )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),mach3 )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),uprim1 )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),uprim2 )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),uprim3 )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),beta )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),zeff )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),temp3 )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),vnewk1 )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),vnewk2 )
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),vnewk3 )
        if(igeo.eq.0)then
          j=j+1
          scalval(j) = nf_put_var_double(ncid, rhid(j),dbeam )
          j=j+1
          scalval(j) = nf_put_var_double(ncid, rhid(j),teti )
          j=j+1
          scalval(j) = nf_put_var_double(ncid, rhid(j),pk )
          j=j+1
          scalval(j) = nf_put_var_double(ncid, rhid(j),eps )
          j=j+1
          scalval(j) = nf_put_var_double(ncid, rhid(j),epsa )
          j=j+1
          scalval(j) = nf_put_var_double(ncid, rhid(j),epsl )
          j=j+1
          scalval(j) = nf_put_var_double(ncid, rhid(j),shat )
          j=j+1
          scalval(j) = nf_put_var_double(ncid, rhid(j),shift )
        elseif(igeo.eq.1)then
          j=j+1
          scalval(j) = nf_put_var_double(ncid, rhid(j),nione_loc )
          j=j+1
          scalval(j) = nf_put_var_double(ncid, rhid(j),fastionfrac_loc )
          j=j+1
          scalval(j) = nf_put_var_double(ncid, rhid(j),tiote_loc )
          j=j+1
          scalval(j) = nf_put_var_double(ncid, rhid(j),rmin_loc )
          j=j+1
          scalval(j) = nf_put_var_double(ncid, rhid(j),rmaj0_loc )
          j=j+1
          scalval(j) = nf_put_var_double(ncid, rhid(j),q_loc )
          j=j+1
          scalval(j) = nf_put_var_double(ncid, rhid(j),shat_mhd_loc )
          j=j+1
          scalval(j) = nf_put_var_double(ncid, rhid(j),alpha_mhd_loc )
          j=j+1
          scalval(j) = nf_put_var_double(ncid, rhid(j),shift_loc )
          j=j+1
          scalval(j) = nf_put_var_double(ncid, rhid(j),kappa_loc )
          j=j+1
          scalval(j) = nf_put_var_double(ncid, rhid(j),delta_loc )
          j=j+1
          scalval(j) = nf_put_var_double(ncid, rhid(j),s_kappa_loc )
          j=j+1
          scalval(j) = nf_put_var_double(ncid, rhid(j),s_delta_loc )
        endif
      endif
      if(igks.eq.0.or.(igks.ge.3.and.(jin_m.eq.jout_m)))then
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),agammas(1))
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),dgammas(1))
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),afreqs(1))
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),kys(1))
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),peflxa)
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),eeflxa)
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j),eiflxa)
c
c define rdim
c          stat = nf_def_dim(ncid, 'one',ntgridl+ntgrid+1,rdim)
c          if (stat .ne. nf_noerr) call handle_err(stat)
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    theta(-ntgridl:ntgrid))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    real(phinorm(-ntgridl:ntgrid,1)))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    aimag(phinorm(-ntgridl:ntgrid,1)))
          if(cbetae.gt.0.0)then
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      real(aparnorm(-ntgridl:ntgrid,1)))
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      aimag(aparnorm(-ntgridl:ntgrid,1)))
            if(i_bpar.ne.0)then
              k=k+1
              oneval(k) = nf_put_var_double(ncid, 
     >        rhid1(k),real(bparnorm(-ntgridl:ntgrid,1)))
              k=k+1
              oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >        aimag(bparnorm(-ntgridl:ntgrid,1)))
            endif
          endif
          if(igeo.eq.1)then
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      sintheta_geo(-ntgridl:ntgrid))
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      costheta_geo(-ntgridl:ntgrid))
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      b_geo(-ntgridl:ntgrid))
          endif
        elseif(igks.eq.1.or.igks.eq.2.or.igks.eq.7)then  ! igks = 1,2 or 7 only
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    aky_m(0:20))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    ky_m(0:20))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    anrate_m(0:20))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    dnrate_m(0:20))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    anfreq_m(0:20))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    ky_mks(0:20))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    gamma_mks(0:20))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    dgamma_mks(0:20))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    freq_mks(0:20))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    peflx_m(0:20))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    qeflx_m(0:20))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    qiflx_m(0:20))
          if(igks_model.eq.1)then ! TGLF
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      gamma_ion(0:20))
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      freq_ion(0:20))
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      gamma_electron(0:20))
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      freq_electron(0:20))
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      phi_bar_m(0:20))
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      ne_bar_m(0:20))
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      te_bar_m(0:20))
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      ti_bar_m(0:20))
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      ne_te_phase_ion(0:20))
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      ne_te_phase_electron(0:20))
          endif
        endif
        if(igks.ge.3.and.(jin_m.ne.jout_m))then
c define rdim
c          stat = nf_def_dim(ncid, 'one',jmaxm+1,rdim)
c          if (stat .ne. nf_noerr) call handle_err(stat)
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    rho(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    ne_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    ni_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    nz_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    nfst_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    te_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    ti_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    ptot_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    pfast_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    vphi_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    vecur_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    zeff_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    betae_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    q_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    shat_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    alpha_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    delta_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    elong_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    drhodr(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    zpne_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    zpni_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    zpnz_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    zpte_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    zpti_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    mach_i_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    mach_e_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    gamma_p_i_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    gamma_p_e_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    egamma_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    csda_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    rhosda_exp(0:jmaxm))
          k=k+1
          oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >    ky_m(0:jmaxm))
          if(igks.le.4)then
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      anrate_m(0:jmaxm))
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      dnrate_m(0:jmaxm))
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      anfreq_m(0:jmaxm))
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      ky_mks(0:jmaxm))
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      gamma_mks(0:jmaxm))
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      dgamma_mks(0:jmaxm))
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      freq_mks(0:jmaxm))
            if(igks.eq.3.and.igks_model.eq.1)then
              k=k+1
              oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >        gamma_ion(0:jmaxm))
              k=k+1
              oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >        freq_ion(0:jmaxm))
              k=k+1
              oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >        gamma_electron(0:jmaxm))
              k=k+1
              oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >        freq_electron(0:jmaxm))
              k=k+1
              oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >        phi_bar_m(0:jmaxm))
              k=k+1
              oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >        ne_bar_m(0:jmaxm))
              k=k+1
              oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >        te_bar_m(0:jmaxm))
              k=k+1
              oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >        ti_bar_m(0:jmaxm))
              k=k+1
              oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >        ne_te_phase_ion(0:jmaxm))
              k=k+1
              oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >        ne_te_phase_electron(0:jmaxm))
              k=k+1
              oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >        mhd_DR_m(0:jmaxm))
            endif
          elseif(igks.eq.5)then
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      zpte_crit_m(0:jmaxm))
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      ky_mks(0:jmaxm))
          elseif(igks.eq.6)then
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      zpti_crit_m(0:jmaxm))
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      anfreq_m(0:jmaxm))
            k=k+1
            oneval(k) = nf_put_var_double(ncid, rhid1(k),
     >      ky_mks(0:jmaxm))
          endif
        endif
c
      call ncclos(ncid, rcode)

c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
  50  format('write',2x,i2,2x,0p1f9.6,1p6e14.6)
 100  format(/,' Writing output to gksout.nc ...')
 150  format(/,3x,' Writing output to ... ',a30,/)
c
      return
      end
