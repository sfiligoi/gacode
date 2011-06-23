      subroutine write_netcdf
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c     This routine writes out xptor results to a NetCDF 
c     file 'xptor_out.nc'
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      implicit none
c
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/ptor.m'
      include '../inc/glf.m'
      include '../inc/vtrans.m'
      include 'netcdf.inc'
      include 'mpif.h'
c
      character*(*) time_units, temp_units, leng_units, den_units, 
     &              field_units, cur_units, pow_units, mom_units
      integer txlen
      parameter (time_units='secs',temp_units='keV',leng_units='m')
      parameter (den_units='1.e19 m^-3', mom_units='NT-M^2')
      parameter (field_units='T', cur_units='MA', pow_units='MW')
      parameter (txlen=15)
      integer i,j,k,val, lenid, ncid, rcode
      integer rhidd, rhid (500), rhid1 (500)
      integer scalval(500), oneval(500)
      integer stat, rdim
      integer CHID, timeid, txid
      integer txdims, cstart, ccount
      real*8 relval
      character*40 txval
c      data txval /"example string"/
c
      txval(txlen:txlen) = char(0)  ! null terminate
c
      if(i_proc.eq.0) then
c
      val = nf_create('xptor_out.nc', nf_clobber, ncid)
      if (val .ne. nf_noerr) call handle_err(val)
c
c Define scalars and attributes
c
      val = nf_def_var (ncid, 'nsteps', nf_short, 0,0,rhidd)
      j=1
      scalval(j)=nf_def_dim(ncid,"chid", 40, CHID)
      txdims = CHID
      scalval(j) = nf_def_var (ncid,"discharge",nf_char,
     &             1,txdims,txid)
c      if(scalval(j) .ne. nf_noerr) call handle_err(scalval(j))
      scalval(j)=nf_put_att_text(ncid,txid,'definition',
     &           len('shot number'), 'shot number')
c
      j=j+1
      scalval(j) = nf_def_var (ncid, 'nrho', nf_int, 0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'definition',15,
     &           'no. of grid pts')
      j=j+1
      scalval(j) = nf_def_var (ncid, 'nrho_bc', nf_int, 0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'definition',30,
     &           'grid pt for boundary condition')
      j=j+1
      scalval(j) = nf_def_var (ncid, 'xp_time', nf_double, 0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',len(time_units),
     &           time_units)
      j=j+1
      scalval(j) = nf_def_var(ncid,'rmajor_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',len(leng_units),
     &           leng_units)
      j=j+1
      scalval(j) = nf_def_var(ncid,'rgeo_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',len(leng_units),
     &           leng_units)
      j=j+1
      scalval(j) = nf_def_var(ncid,'rminor_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',len(leng_units),
     &           leng_units)
      j=j+1
      scalval(j) = nf_def_var(ncid,'bt_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',len(field_units),
     &           field_units)
      j=j+1
      scalval(j) = nf_def_var(ncid,'ip_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',len(cur_units),
     &           cur_units)
      j=j+1
      scalval(j) = nf_def_var(ncid,'ohmcur_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',len(cur_units),
     &           cur_units)
      j=j+1
      scalval(j) = nf_def_var(ncid,'bootcur_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',len(cur_units),
     &           cur_units)
      j=j+1
      scalval(j) = nf_def_var(ncid,'nbcur_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',len(cur_units),
     &           cur_units)
      j=j+1
      scalval(j) = nf_def_var(ncid,'rfcur_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',len(cur_units),
     &           cur_units)
      j=j+1
      scalval(j) = nf_def_var(ncid,'neline_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',
     &           len('line avg elec density 10^19/m^3'),
     &           'line avg elec density 10^19/m^3')
      j=j+1
      scalval(j) = nf_def_var(ncid,'nevol_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',
     &           len('vol avg elec density 10^19/m^3'),
     &           'vol avg elec density 10^19/m^3')
      j=j+1
      scalval(j) = nf_def_var(ncid,'niline_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',
     &           len('line avg ion density 10^19/m^3'),
     &           'line avg ion density 10^19/m^3')
      j=j+1
      scalval(j)=nf_def_var(ncid,'kappaa_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',15,
     &           'kappa at rho(a)')
      j=j+1
      scalval(j)=nf_def_var(ncid,'deltaa_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',15,
     &           'delta at rho(a)')
      j=j+1
      scalval(j)=nf_def_var(ncid,'wstr_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',25,
     &           'total stored energy in MJ')
      j=j+1
      scalval(j)=nf_def_var(ncid,'wstr_ped_exp',nf_double,
     &            0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',35,
     &           'total pedestal stored energy in MJ')
      j=j+1
      scalval(j)=nf_def_var(ncid,'wstr_inc_exp',nf_double,
     &            0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',37,
     &           'total incremental stored energy in MJ')
      j=j+1
      scalval(j)=nf_def_var(ncid,'wstr_core_exp',nf_double,
     &            0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',24,
     &           'core stored energy in MJ')
      j=j+1
      scalval(j)=nf_def_var(ncid,'taue_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'def',34,
     &           'exp energy confinement time (secs)')
      j=j+1
      scalval(j)=nf_def_var(ncid,'taue_89p',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'def',32,
     &           'ITER89P energy confinement time')
      j=j+1
      scalval(j)=nf_def_var(ncid,'taue_97L',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'def',31,
     &           'ITER97L energy confinement time')
      j=j+1
      scalval(j)=nf_def_var(ncid,'taue_98y2',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'def',32,
     &           'ITER98y2 energy confinement time')
      j=j+1
      scalval(j)=nf_def_var(ncid,'powetot_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'def',25,
     &           'Total exp elec power (MW)')
      j=j+1
      scalval(j)=nf_def_var(ncid,'powitot_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'def',24,
     &           'Total exp ion power (MW)')
      j=j+1
      scalval(j)=nf_def_var(ncid,'powe_aux_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'def',23,
     &           'Aux exp elec power (MW)')
      j=j+1
      scalval(j)=nf_def_var(ncid,'powi_aux_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'def',22,
     &           'Aux exp ion power (MW)')
      j=j+1
      scalval(j)=nf_def_var(ncid,'pohm_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'def',20,
     &           'Exp Ohmic power (MW)')
      j=j+1
      scalval(j)=nf_def_var(ncid,'prad_exp',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'def',23,
     &           'Exp Radiated power (MW)')
c
      j=j+1
      scalval(j)=nf_def_var(ncid,'wstr_m',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',29,
     &           'predicted stored energy in MJ')
      j=j+1
      scalval(j)=nf_def_var(ncid,'wstr_inc_m',nf_double,
     &            0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',41,
     &           'predicted incremental stored energy in MJ')
      j=j+1
      scalval(j)=nf_def_var(ncid,'wstr_core_m',nf_double,
     &            0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'units',34,
     &           'predicted core stored energy in MJ')
      j=j+1
      scalval(j)=nf_def_var(ncid,'taue_m',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'def',33,
     &           'predicted energy confinement time')
      j=j+1
      scalval(j)=nf_def_var(ncid,'nt_taue_m',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'def',17,
     &           'predicted nT*taue')
      j=j+1
      scalval(j)=nf_def_var(ncid,'volavete_m',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'def',25,
     &           'vol avg predicted Te(keV)')
      j=j+1
      scalval(j)=nf_def_var(ncid,'volaveti_m',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'def',25,
     &           'vol avg predicted Ti(keV)')
      j=j+1
      scalval(j)=nf_def_var(ncid,'volaveniti_m',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'def',22,
     &           'vol avg predicted ni*Ti')
      j=j+1
      scalval(j)=nf_def_var(ncid,'neped_m',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'def',
     &           len('predicted elec density at BC'),
     &           'predicted elec density at BC')
      j=j+1
      scalval(j)=nf_def_var(ncid,'teped_m',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'def',
     &           len('Te at BC (kev)'), 'Te at BC (kev)')
      j=j+1
      scalval(j)=nf_def_var(ncid,'tiped_m',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'def',
     &           len('Ti at BC (kev)'), 'Ti at BC (kev)')
      if(abs(irad).gt.0) then
       j=j+1
       scalval(j)=nf_def_var(ncid,'prad_brem_m',nf_double,0,0,rhid(j))
       scalval(j)=nf_put_att_text(ncid,rhid(j),'def',
     &    len('predicted Brem. radiation (MW)'), 
     &    'predicted Brem. radiation (MW)')
       j=j+1
       scalval(j)=nf_def_var(ncid,'prad_sync_m',nf_double,0,0,rhid(j))
       scalval(j)=nf_put_att_text(ncid,rhid(j),'def',
     &    len('predicted Sync. radiation (MW)'), 
     &    'predicted Sync. radiation (MW)')
      endif
      if(ialpha.eq.1) then
       j=j+1
       scalval(j)=nf_def_var(ncid,'pfus_m',nf_double,0,0,rhid(j))
       scalval(j)=nf_put_att_text(ncid,rhid(j),'def',
     &    len('predicted P_fusion (MW)'), 'predicted P_fusion (MW)')
       j=j+1
       scalval(j)=nf_def_var(ncid,'pnet_m',nf_double,0,0,rhid(j))
       scalval(j)=nf_put_att_text(ncid,rhid(j),'def',
     &    len('predicted P_net (MW)'), 'predicted P_net (MW)')
      endif
      if(abs(iohm).gt.0) then
       j=j+1
       scalval(j)=nf_def_var(ncid,'pohm_m',nf_double,0,0,rhid(j))
       scalval(j)=nf_put_att_text(ncid,rhid(j),'def',
     &    len('predicted P_ohmic (MW)'), 'predicted P_ohmic (MW)')
      endif
      j=j+1
      scalval(j)=nf_def_var(ncid,'normS',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'def',
     &           len('normS'), 'normS')
      j=j+1
      scalval(j)=nf_def_var(ncid,'normDT',nf_double,0,0,rhid(j))
      scalval(j)=nf_put_att_text(ncid,rhid(j),'def',
     &           len('normDT'), 'normDT')

c
c Define 1d vars and attributes
c
      stat = nf_def_dim(ncid, 'one',jmaxm+1,rdim)
      if (stat .ne. nf_noerr) call handle_err(stat)
      j=1
      oneval(j) = nf_def_var (ncid,'rho',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &          len('dimensionless'),'dimensionless')
      j=j+1
      oneval(j) = nf_def_var (ncid,'te_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',len(temp_units),
     &           temp_units)
      j=j+1
      oneval(j) = nf_def_var (ncid,'ti_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',len(temp_units),
     &           temp_units)
      j=j+1
      oneval(j)=nf_def_var(ncid,'tfast_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',len(temp_units),
     &           temp_units)
      j=j+1
      oneval(j) = nf_def_var (ncid,'ne_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. elec density 10^19 m^-3'),
     & 'exp. elec density 10^19 m^-3')
      j=j+1
      oneval(j) = nf_def_var (ncid,'ni_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. main ion density 10^19 m^-3'),
     & 'exp. main ion density 10^19 m^-3')
      j=j+1
      oneval(j) = nf_def_var (ncid,'nz_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. impurity density 10^19 m^-3'),
     & 'exp. impurity density 10^19 m^-3')
      j=j+1
      oneval(j)=nf_def_var(ncid,'nfast_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. fast ion density 10^19 m^-3'),
     & 'exp. fast ion density 10^19 m^-3')
      j=j+1
      oneval(j) = nf_def_var (ncid,'zeff_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. Zeff profile'),'exp. Zeff profile')
      j=j+1
      oneval(j) = nf_def_var (ncid,'curden_exp',
     &  nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &  len('total current density (A/m^2)'),
     &  'total current density (A/m^2)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'curboot_exp',
     &  nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &  len('bootstrap current density (A/m^2)'),
     &  'bootstrap current density (A/m^2)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'curbeam_exp',
     &  nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &  len('beam current density (A/m^2)'),
     &  'beam current density (A/m^2)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'currf_exp',
     &  nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &  len('RF current density (A/m^2)'),
     &  'RF current density (A/m^2)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'curohm_exp',
     &  nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &  len('Ohmic current density (A/m^2)'),
     &  'Ohmic current density (A/m^2)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'q_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. safety factor'),'exp. safety factor')
      j=j+1
      oneval(j) = nf_def_var (ncid,'shat_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. magnetic shear'),'exp. magnetic shear')
      j=j+1
      oneval(j) = nf_def_var(ncid,'alpha_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. MHD alpha'),'exp. MHD alpha')
      j=j+1
      oneval(j) = nf_def_var(ncid,'betae_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. betae'),'exp. betae')
      j=j+1
      oneval(j) = nf_def_var(ncid,'ptot_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. total pressure (keV/m^3)'),
     & 'exp. total pressure (keV/m^3)')
      j=j+1
      oneval(j)=nf_def_var(ncid,'torque_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. torque density (nt-m/m^3)'),
     & 'exp. torque density (nt-m/m^3)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vphi_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. ion toroidal velocity 10^5 m/s'),
     & 'exp. ion toroidal velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var(ncid,'vphie_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. electron plasma toroidal velocity 10^5 m/s'),
     & 'exp. electron plasma toroidal velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var(ncid,'vphiz_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. impurity plasma toroidal velocity 10^5 m/s'),
     & 'exp. impurity plasma toroidal velocity 10^5 m/s')
      j=j+1
      oneval(j)=nf_def_var (ncid,'angrot_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &      len('exp. impurity angular rotation velocity 10^5 rad/s'),
     &      'exp. impurity angular rotation velocity 10^5 rad/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vpar_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. parallel velocity 10^5 m/s'),
     & 'exp. parallel velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var(ncid,'vpare_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. electron parallel velocity 10^5 m/s'),
     & 'exp. electron parallel velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var(ncid,'vparz_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp.impurity  parallel velocity 10^5 m/s'),
     & 'exp. impurity parallel velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vexb_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. ExB velocity 10^5 m/s'),
     & 'exp. ExB velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vdia_exp(1)'
     & ,nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. electrondiamagnetic velocity 10^5 m/s'),
     & 'exp. electron diamagnetic velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vdia_exp(2)'
     & ,nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. ion diamagnetic velocity 10^5 m/s'),
     & 'exp. ion diamagnetic velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vdia_exp(3)'
     & ,nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. impurity diamagnetic velocity 10^5 m/s'),
     & 'exp. impurity diamagnetic velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vneo_exp(1)'
     & ,nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. electron neo. pol. velocity 10^5 m/s'),
     & 'exp. electron neo. pol. velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vneo_exp(2)'
     & ,nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. ion neo. pol. velocity 10^5 m/s'),
     & 'exp. ion neo. pol. velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vneo_exp(3)'
     & ,nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. impurity neo. pol. velocity 10^5 m/s'),
     & 'exp. impurity neo. pol. velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'mach_exp(1)'
     & ,nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. electron parallel mach number'),
     & 'exp. electron parallel mach number')
      j=j+1
      oneval(j) = nf_def_var (ncid,'mach_exp(2)'
     & ,nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. ion parallel mach number'),
     & 'exp. ion parallel mach number')
      j=j+1
      oneval(j) = nf_def_var (ncid,'mach_exp(3)'
     & ,nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp. impurity parallel mach number'),
     & 'exp. impurity parallel mach number')
      j=j+1
      oneval(j) = nf_def_var (ncid,'zpne_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &          len('exp a/Lne (dimensionless)'),
     &          'exp a/Lne (dimensionless)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'zpni_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &          len('exp a/Lni (dimensionless)'),
     &          'exp a/Lni (dimensionless)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'zpte_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &          len('exp a/Lte (dimensionless)'),
     &          'exp a/Lte (dimensionless)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'zpti_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &          len('exp a/Lti (dimensionless)'),
     &          'exp a/Lti (dimensionless)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'zpne2_exp',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('exp der. of a/Lne'),'exp der. of a/Lne')
      j=j+1
      oneval(j) = nf_def_var (ncid,'zpte2_exp',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('exp der. of a/Lte'),'exp der. of a/Lte')
      j=j+1
      oneval(j) = nf_def_var (ncid,'zpti2_exp',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('exp der. of a/Lti'),'exp der. of a/Lti')
      j=j+1
      oneval(j)=nf_def_var(ncid,'flow_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('total exp particle flow KA'),
     &           'total exp particle flow KA')
      j=j+1
      oneval(j) = nf_def_var (ncid,'powe_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('exp integrated elec power (MW)'),
     &           'exp integrated elec power (MW)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'powi_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('exp integrated ion power (MW)'),
     &          'exp integrated ion power (MW)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'stress_tor_exp',nf_double,1,
     &           rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('integrated toroidal momentum flux NT-M'),
     &          'integrated toroidal momentum flux NT-M')
      j=j+1
      oneval(j)=nf_def_var(ncid,'pow_ei_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('exp exchange power MW'),
     &           'exp exchange power MW')
      j=j+1
      oneval(j)=nf_def_var(ncid,'powe_beam_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('exp elec beam power (MW)'),
     &           'exp elec beam power (MW)')
      j=j+1
      oneval(j)=nf_def_var(ncid,'powi_beam_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('exp ion beam power (MW)'),
     &           'exp ion beam power (MW)')
      j=j+1
      oneval(j)=nf_def_var(ncid,'powe_rf_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('exp elec RF power (MW)'),
     &           'exp elec RF power (MW)')
      j=j+1
      oneval(j)=nf_def_var(ncid,'powi_rf_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('exp ion RF power (MW)'),
     &           'exp ion RF power (MW)')
      j=j+1
      oneval(j)=nf_def_var(ncid,'powe_oh_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('exp Ohmic heating power (MW)'),
     &           'exp Ohmic heating power (MW)')
      j=j+1
      oneval(j)=nf_def_var(ncid,'powe_rad_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('exp elec radiated power (MW)'),
     &           'exp elec radiated power (MW)')
      j=j+1
      oneval(j)=nf_def_var(ncid,'powi_ion_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('exp ion ionization power (MW)'),
     &           'exp ion ionization power (MW)')
      j=j+1
      oneval(j)=nf_def_var(ncid,'powi_cx_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('exp ion CX power (MW)'),
     &           'exp ion CX power (MW)')
      j=j+1
      oneval(j)=nf_def_var(ncid,'powe_wdot_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('exp elec wdot power (MW)'),
     &           'exp elec wdot power (MW)')
      j=j+1
      oneval(j)=nf_def_var(ncid,'powi_wdot_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('exp ion wdot power (MW)'),
     &           'exp ion wdot power (MW)')
      j=j+1
      oneval(j)=nf_def_var(ncid,'flow_wall_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('exp wall particle flow KA'),
     &           'exp wall particle flow KA')
      j=j+1
      oneval(j)=nf_def_var(ncid,'flow_beam_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('exp beam particle flow KA'),
     &           'exp beam particle flow KA')
      j=j+1
      oneval(j)=nf_def_var(ncid,'sion_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &      len('exp source due to ionization (#/(m**3*sec)), d'),
     &          'exp source due to ionization (#/(m**3*sec)), d')
      j=j+1
      oneval(j)=nf_def_var(ncid,'srecom_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &      len('exp source due to recombination (#/(m**3*sec)), d'),
     &          'exp source due to recombination (#/(m**3*sec)), d')
      j=j+1
      oneval(j)=nf_def_var(ncid,'scx_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp source due to cx with thermal neut (#/(m**3*sec)), d'),
     &     'exp source due to cx with thermal neut (#/(m**3*sec)), d')
      j=j+1
      oneval(j)=nf_def_var(ncid,'sbcx_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &   len('exp source due to cx with beam neut (#/(m**3*sec)), d'),
     &       'exp source due to cx with beam neut (#/(m**3*sec)), d')
      j=j+1
      oneval(j)=nf_def_var(ncid,'s_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &   len('exp total source rate (#/(m**3*sec)), d'),
     &       'exp total source rate (#/(m**3*sec)), d')
      j=j+1
      oneval(j)=nf_def_var(ncid,'sdot_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &   len('exp time deriv of total source rate (#/(m**3*sec)), d'),
     &       'exp time deriv of total source rate (#/(m**3*sec)), d')
      j=j+1
      oneval(j)=nf_def_var(ncid,'enn_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &          len('exp neut density (#/m**3), d'),
     &              'exp neut density (#/m**3), d')
      j=j+1
      oneval(j)=nf_def_var(ncid,'ennw_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &   len('exp neut density due to wall source (#/m**3), d'),
     &       'exp neut density due to wall source (#/m**3), d')
      j=j+1
      oneval(j)=nf_def_var(ncid,'ennv_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &   len('exp neut density due to vol source (#/m**3), d'),
     &       'exp neut density due to vol source (#/m**3), d')
      j=j+1
      oneval(j)=nf_def_var(ncid,'volsn_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &   len('exp vol source of neutrals (#/(m**3*sec)), d'),
     &       'exp vol source of neutrals (#/(m**3*sec)), d')
      j=j+1
      oneval(j)=nf_def_var(ncid,'sbeame_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &   len('exp neut density due to beam elec source (#/(m**3*sec))'),
     &       'exp neutl density due to beam elec source (#/(m**3*sec))')
      j=j+1
      oneval(j)=nf_def_var(ncid,'sbeam_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &len('exp neut density due to beam ion source (#/(m**3*sec))'),
     &    'exp neut density due to beam ion source (#/(m**3*sec))')
      j=j+1
      oneval(j)=nf_def_var(ncid,'qbeame_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('exp elec beam power density (W/m^3)'),
     &           'exp elec beam power density (W/m^3)')
      j=j+1
      oneval(j)=nf_def_var(ncid,'qbeami_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('exp ion beam power density (W/m^3)'),
     &           'exp ion beam power density (W/m^3)')
      j=j+1
      oneval(j)=nf_def_var(ncid,'qrfe_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('exp elec rf power density (W/m^3)'),
     &           'exp elec rf power density (W/m^3)')
      j=j+1
      oneval(j)=nf_def_var(ncid,'qrfi_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('exp ion rf power density (W/m^3)'),
     &           'exp ion rf power density (W/m^3)')
      j=j+1
      oneval(j)=nf_def_var(ncid,'qohm_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('exp ohmic power density (W/m^3)'),
     &           'exp ohmic power density (W/m^3)')
      j=j+1
      oneval(j)=nf_def_var(ncid,'qrad_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('exp radiated power density (W/m^3)'),
     &           'exp radiated power density (W/m^3)')
      j=j+1
      oneval(j)=nf_def_var(ncid,'qione_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &   len('exp elec recomb and ionization power density (W/m^3)'),
     &   'exp elec recomb and ionization power density (W/m^3)')
      j=j+1
      oneval(j)=nf_def_var(ncid,'qioni_exp',nf_double,1,
     &          rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &   len('exp ion recomb and ionization power density (W/m^3)'),
     &   'exp ion recomb and ionization power density (W/m^3)')
      j=j+1
      oneval(j) = nf_def_var(ncid,'rmin_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp rminor (m)'),'exp rminor (m)')
      j=j+1
      oneval(j) = nf_def_var(ncid,'rmaj_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp rmajor (m)'),'exp rmajor (m)')
      j=j+1
      oneval(j) = nf_def_var(ncid,'elong_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp kappa'),'exp kappa')
      j=j+1
      oneval(j) = nf_def_var(ncid,'delta_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp avg delta'),'exp avg delta')
      j=j+1
      oneval(j) = nf_def_var(ncid,'grho1_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp grad rho'),'exp grad rho')
      j=j+1
      oneval(j) = nf_def_var(ncid,'grho2_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp grad rho squared'),'exp grad rho squared')
      j=j+1
      oneval(j) = nf_def_var(ncid,'fcap_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp f(psilim)/f(psi)'),'exp f(psilim)/f(psi)')
      j=j+1
      oneval(j) = nf_def_var(ncid,'gcap_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp <(grad rho)**2*(R0/R)**2>'),'
     & exp <(grad rho)**2*(R0/R)**2>')
      j=j+1
      oneval(j) = nf_def_var(ncid,'hcap_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp (dvolume/drho)/(4*pi*pi*R0*rho)'),'
     & exp (dvolume/drho)/(4*pi*pi*R0*rho)')
      j=j+1
      oneval(j) = nf_def_var(ncid,'vol_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp volume of each flux surface (m^3)'),
     &    'exp volume of each flux surface (m^3)')
      j=j+1
      oneval(j)=nf_def_var(ncid,'drhodr_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp drho/dr'),'exp drho/dr')
      j=j+1
      oneval(j) = nf_def_var (ncid,'theta_exp',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('exp theta'), 'exp theta')
      j=j+1
      oneval(j) = nf_def_var (ncid,'zptheta_exp',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('exp zptheta'), 'exp zptheta')
      j=j+1
      oneval(j)=nf_def_var(ncid,'a_pol',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('a_pol'),'a_pol')
      j=j+1
      oneval(j)=nf_def_var(ncid,'a_tor',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('a_tor'),'a_tor')
      j=j+1
      oneval(j)=nf_def_var(ncid,'c_par',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('c_par'),'c_par')
      j=j+1
      oneval(j)=nf_def_var(ncid,'c_per',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('c_per'),'c_per')
      j=j+1
      oneval(j)=nf_def_var(ncid,'c_tor',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('c_tor'),'c_tor')
      j=j+1
      oneval(j)=nf_def_var(ncid,'vprime1',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp vprime1'),'exp vprime1')
      j=j+1
      oneval(j)=nf_def_var(ncid,'vprime2',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp vprime2'),'exp vprime2')
      j=j+1
      oneval(j)=nf_def_var(ncid,'sfarea_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp surface area of each flux surface (m^2)'),
     & 'exp surface area of each flux surface (m^2)')
      j=j+1
      oneval(j)=nf_def_var(ncid,'cxarea_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp cross sect area of each flux surface (m^2)'),
     & 'exp cross sect area of each flux surface (m^2)')
      j=j+1
      oneval(j)=nf_def_var(ncid,'geofac',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp geometry factor'),'exp geometry factor')
      j=j+1
      oneval(j) = nf_def_var(ncid,'csda_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp c_s/a'),'exp c_s/a')
      j=j+1
      oneval(j)=nf_def_var(ncid,'rhosda_exp',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('exp rho_s/a'),'exp rho_s/a')
c
      j=j+1
      oneval(j) = nf_def_var (ncid,'ti_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',len(temp_units),
     &           temp_units)
      j=j+1
      oneval(j) = nf_def_var (ncid,'te_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',len(temp_units),
     &           temp_units)
      j=j+1
      oneval(j) = nf_def_var (ncid,'ne_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('elec density 10^19 m^-3'),'elec density 10^19 m^-3')
      j=j+1
      oneval(j) = nf_def_var (ncid,'ni_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('ion density 10^19 m^-3'),'ion density 10^19 m^-3')
      j=j+1
      oneval(j) = nf_def_var (ncid,'nz_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('impurity density 10^19 m^-3'),'impurity density 10^19 m^-3')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vphi_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('ion toroidal velocity 10^5 m/s'),
     & 'ion toroidal velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vphie_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('electron toroidal velocity 10^5 m/s'),
     & 'electron toroidal velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vphiz_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('impurity toroidal velocity 10^5 m/s'),
     & 'impurity toroidal velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vpar_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('ion parallel velocity 10^5 m/s'),
     & 'ion parallel velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vpare_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('electron parallel velocity 10^5 m/s'),
     & 'electron parallel velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vparz_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('impurity parallel velocity 10^5 m/s'),
     & 'impurity parallel velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vpol_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('poloidal velocity 10^5 m/s'),'poloidal velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'xnu_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('nuei / cs/a'),'nuei / cs/a')
      j=j+1
      oneval(j) = nf_def_var (ncid,'nuei_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('nuei rad/sec '),
     & 'collisional electron-ion energy exchange rate')
      j=j+1
      oneval(j) = nf_def_var (ncid,'nu_pol_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('nu_pol rad/sec '),'collisional poloidal flow damping rate')
      j=j+1
      oneval(j) = nf_def_var (ncid,'kpol_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     & len('kpol dimensionless '),
     & 'neoclassical poloidal flow coefficient')
      j=j+1
      oneval(j) = nf_def_var (ncid,'alpha_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('MHD alpha'),'MHD alpha')
      j=j+1
      oneval(j) = nf_def_var (ncid,'betae_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('electron beta'),'electron beta')
      j=j+1
      oneval(j) = nf_def_var (ncid,'zpne_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('a/Lne (dimensionless)'),'a/Lne (dimensionless)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'zpni_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('a/Lni (dimensionless)'),'a/Lni (dimensionless)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'zpte_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('a/LTe (dimensionless)'),'a/LTe (dimensionless)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'zpti_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('a/LTi (dimensionless)'),'a/LTi (dimensionless)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'interchange_DR_m',nf_double,1, 
     & rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('DR (dimensionless)'),'DR (dimensionless)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'anrate_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('gamma / (cs/a)'),'gamma / (cs/a)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'anfreq_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('omega / (cs/a)'),'omega / (cs/a)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'egamma_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('gamma_E / (cs/a)'),'gamma_E / (cs/a)')
      j=j+1
      oneval(j)=nf_def_var (ncid,'gamma_p_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('gamma_p / (cs/a)'),'gamma_p / (cs/a)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vexb_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('v_ExB, ExB velocity 10^5 m/s'),
     & 'v_ExB, ExB velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vdia_m(1)'
     & ,nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('v_dia, electron diamagnetic velocity 10^5 m/s'),
     & 'v_dia, electron diamagnetic velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vdia_m(2)'
     & ,nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('v_dia, ion diamagnetic velocity 10^5 m/s'),
     & 'v_dia, ion diamagnetic velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vdia_m(3)'
     & ,nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('v_dia, impurity diamagnetic velocity 10^5 m/s'),
     & 'v_dia, impurity diamagnetic velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vneo_m(1)'
     & ,nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('v_neo, electron neo. pol. velocity 10^5 m/s'),
     & 'v_neo, electron neo. pol. velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vneo_m(2)'
     & ,nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('v_neo, ion neo. pol. velocity 10^5 m/s'),
     & 'v_neo, ion neo. pol. velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vneo_m(3)'
     & ,nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('v_neo, impurity neo. pol. velocity 10^5 m/s'),
     & 'v_neo, impurity neo. pol. velocity 10^5 m/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'mach_m(1)'
     & ,nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('model electron parallel mach number'),
     & 'model electron parallel mach number')
      j=j+1
      oneval(j) = nf_def_var (ncid,'mach_m(2)'
     & ,nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('model ion parallel mach number'),
     & 'model ion parallel mach number')
      j=j+1
      oneval(j) = nf_def_var (ncid,'mach_m(3)'
     & ,nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('model impurity parallel mach number'),
     & 'model impurity parallel mach number')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vstar_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('v_star'),'v_star')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vepol_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('vepol'),'vepol')
      j=j+1
      oneval(j) = nf_def_var (ncid,'vetor_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def.',
     & len('vetor'),'vetor')
      j=j+1
      oneval(j) = nf_def_var (ncid,'diff_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &          len('particle diffusivity m^2/s'),
     &          'particle diffusivity m^2/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'chii_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &          len('ion thermal diffusivity m^2/s'),
     &          'ion thermal diffusivity m^2/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'chie_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &          len('electron thermal diffusivity m^2/s'),
     &          'electron thermal diffusivity m^2/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'chie_e_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &          len('electron high-k thermal diffusivity m^2/s'),
     &          'electron high-k thermal diffusivity m^2/s')
      j=j+1
      oneval(j)=nf_def_var(ncid,'chii_neo_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &          len('ion neo thermal diffusivity m^2/s'),
     &          'ion neo thermal diffusivity m^2/s')
      j=j+1
      oneval(j) = nf_def_var(ncid,'etaphi_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &          len('tor. momentum diffusivity m^2/s'),
     &          'tor. momentum diffusivity m^2/s')
      j=j+1
      oneval(j) = nf_def_var (ncid,'cgyrobohm_m',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &          len('gb diffusivity m^2/s'),
     &          'gb diffusivity m^2/s')
      j=j+1
      oneval(j) = nf_def_var(ncid,'flowe_glf',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('glf electron particle flux KA'),
     &          'glf electron particle flux KA')
      j=j+1
      oneval(j) = nf_def_var(ncid,'flowi_glf',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('glf ion particle flux KA'),
     &          'glf ion particle flux KA')
      j=j+1
      oneval(j) = nf_def_var(ncid,'flowz_glf',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('glf impurity particle flux KA'),
     &          'glf impurity particle flux KA')
      j=j+1
      oneval(j) = nf_def_var (ncid,'powe_glf',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('glf electron power flux MW'),
     &          'glf electron power flux MW')
      j=j+1
      oneval(j) = nf_def_var (ncid,'powi_glf',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('glf ion power flux MW'),
     &          'glf ion power flux MW')
      j=j+1
      oneval(j) = nf_def_var (ncid,'powz_glf',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('glf impurity power flux MW'),
     &          'glf impurity power flux MW')
      j=j+1
      oneval(j) = nf_def_var (ncid,'stress_tor_i_glf',nf_double,1,
     &           rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('glf toroidal ion momentum flux NT-M'),
     &          'glf toroidal ion momentum flux NT-M')
      j=j+1
      oneval(j) = nf_def_var (ncid,'stress_tor_z_glf',nf_double,1,
     &           rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('glf toroidal impurity momentum flux NT-M'),
     &          'glf toroidal impurity momentum flux NT-M')
      j=j+1
      oneval(j) = nf_def_var (ncid,'stress_par_i_glf',nf_double,1,
     &           rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('glf parallel ion momentum flux NT-M'),
     &          'glf parallel ion momentum flux NT-M')
      j=j+1
      oneval(j) = nf_def_var (ncid,'stress_par_z_glf',nf_double,1,
     &           rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('glf parallel impurity momentum flux NT-M'),
     &          'glf parallel impurity momentum flux NT-M')
      j=j+1
      oneval(j) = nf_def_var(ncid,'flowe_neo',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('neo electron particle flux KA'),
     &          'neo electron particle flux KA')
      j=j+1
      oneval(j) = nf_def_var(ncid,'flowi_neo',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('neo ion particle flux KA'),
     &          'neo ion particle flux KA')
      j=j+1
      oneval(j) = nf_def_var(ncid,'flowz_neo',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('neo impurity particle flux KA'),
     &          'neo impurity particle flux KA')
      j=j+1
      oneval(j) = nf_def_var (ncid,'powe_neo',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('neo electron power flux MW'),
     &          'neo electron power flux MW')
      j=j+1
      oneval(j) = nf_def_var (ncid,'powi_neo',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('neo ion power flux MW'),
     &          'neo ion power flux MW')
      j=j+1
      oneval(j) = nf_def_var (ncid,'powz_neo',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('neo impurity power flux MW'),
     &          'neo impurity power flux MW')
      j=j+1
      oneval(j) = nf_def_var (ncid,'stress_tor_i_neo',nf_double,1,
     &           rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('neo ion toroidal momentum flux NT-M'),
     &          'neo ion toroidal momentum flux NT-M')
      j=j+1
      oneval(j) = nf_def_var (ncid,'stress_tor_z_neo',nf_double,1,
     &           rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('neo impurity toroidal momentum flux NT-M'),
     &          'neo impurity toroidal momentum flux NT-M')
      j=j+1
      oneval(j) = nf_def_var (ncid,'stress_par_i_neo',nf_double,1,
     &           rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('neo ion parallel momentum flux NT-M'),
     &          'neo ion parallel momentum flux NT-M')
      j=j+1
      oneval(j) = nf_def_var (ncid,'stress_par_z_neo',nf_double,1,
     &           rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('neo impurity parallel momentum flux NT-M'),
     &          'neo impurity parallel momentum flux NT-M')
      j=j+1
      oneval(j) = nf_def_var (ncid,'flowe_adhoc',nf_double
     &  ,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('adhoc electron particle flux KA'),
     &          'adhoc electron particle flux KA')
      j=j+1
      oneval(j) = nf_def_var (ncid,'flowi_adhoc',nf_double
     &  ,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('adhoc ion particle flux KA'),
     &          'adhoc ion particle flux KA')
      j=j+1
      oneval(j) = nf_def_var (ncid,'flowz_adhoc',nf_double
     &  ,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('adhoc impurity particle flux KA'),
     &          'adhoc impurity particle flux KA')
      j=j+1
      oneval(j) = nf_def_var (ncid,'powe_adhoc',nf_double
     &  ,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('adhoc electron power flux MW'),
     &          'adhoc electron power flux MW')
      j=j+1
      oneval(j) = nf_def_var (ncid,'powi_adhoc',nf_double
     &  ,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('adhoc ion power flux MW'),
     &          'adhoc ion power flux MW')
      j=j+1
      oneval(j) = nf_def_var (ncid,'powz_adhoc',nf_double
     &  ,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('adhoc impurity power flux MW'),
     &          'adhoc impurity power flux MW')
      j=j+1
      oneval(j) = nf_def_var (ncid,'stress_tor_i_adhoc',nf_double
     &  ,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('adhoc ion toroidal momentum flux NT-M'),
     &          'adhoc ion toroidal momentum flux NT-M')
      j=j+1
      oneval(j) = nf_def_var (ncid,'stress_tor_z_adhoc',nf_double
     &  ,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('adhoc impurity toroidal momentum flux NT-M'),
     &          'adhoc impurity toroidal momentum flux NT-M')
      j=j+1
      oneval(j) = nf_def_var (ncid,'stress_par_i_adhoc',nf_double
     &  ,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('adhoc ion parallel momentum flux NT-M'),
     &          'adhoc ion parallel momentum flux NT-M')
      j=j+1
      oneval(j) = nf_def_var (ncid,'stress_par_z_adhoc',nf_double
     &  ,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('adhoc impurity parallel momentum flux NT-M'),
     &          'adhoc impurity parallel momentum flux NT-M')
      j=j+1
      oneval(j) = nf_def_var(ncid,'flowe_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('predicted electron particle flux KA'),
     &          'predictedelectron  particle flux KA')
      j=j+1
      oneval(j) = nf_def_var(ncid,'flowi_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('predicted ion particle flux KA'),
     &          'predicted ion particle flux KA')
      j=j+1
      oneval(j) = nf_def_var(ncid,'flowz_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('predicted impurity particle flux KA'),
     &          'predicted impurity particle flux KA')
      j=j+1
      oneval(j) = nf_def_var (ncid,'powe_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('predicted electron power flux MW'),
     &          'predicted electron power flux MW')
      j=j+1
      oneval(j) = nf_def_var (ncid,'powi_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('predicted ion power flux MW'),
     &          'predicted ion power flux MW')
      j=j+1
      oneval(j) = nf_def_var (ncid,'powz_m',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('predicted impurity power flux MW'),
     &          'predicted impurity power flux MW')
      j=j+1
      oneval(j) = nf_def_var (ncid,'stress_tor_i_m',nf_double,1,
     &           rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('predicted ion toroidal momentum flux NT-M'),
     &          'predicted ion toroidal momentum flux NT-M')
      j=j+1
      oneval(j) = nf_def_var (ncid,'stress_tor_z_m',nf_double,1,
     &           rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('predicted impurity toroidal momentum flux NT-M'),
     &          'predicted impurity toroidal momentum flux NT-M')
      j=j+1
      oneval(j) = nf_def_var (ncid,'stress_par_i_m',nf_double,1,
     &           rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('predicted ion parallel momentum flux NT-M'),
     &          'predicted ion parallel momentum flux NT-M')
      j=j+1
      oneval(j) = nf_def_var (ncid,'stress_par_z_m',nf_double,1,
     &           rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('predicted impurity parallel momentum flux NT-M'),
     &          'predicted impurity parallel momentum flux NT-M')
      j=j+1
      oneval(j) = nf_def_var (ncid,'stress_par_cor_m',nf_double,1,
     &           rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('predicted parallel damping flux NT-M'),
     &          'predicted parallel damping flux NT-M')
      j=j+1
      oneval(j) = nf_def_var (ncid,'Pi_alpha',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('predicted ion alpha heating power MW'),
     &          'predicted ion alpha heating power MW')
      j=j+1
      oneval(j) = nf_def_var (ncid,'Pe_alpha',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('predicted elec alpha heating power MW'),
     &          'predicted elec alpha heating power MW')
      j=j+1
      oneval(j) = nf_def_var (ncid,'Pi_aux',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &        len('predicted ion aux heating power density (MW/m^3)'),
     &        'predicted ion aux heating power density (MW/m^3)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'Pe_aux',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &        len('predicted ion aux heating power density (MW/m^3)'),
     &        'predicted ion aux heating power density (MW/m^3)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'Pohpro',nf_double,1,rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &        len('predicted Ohmic heating power density (MW/m^3)'),
     &        'predicted Ohmic heating power density (MW/m^3)')
      j=j+1
      oneval(j) = nf_def_var (ncid,'pow_ei_cor_m',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('predicted exchange power'),
     &          'predicted exchange power')
      j=j+1
      oneval(j) = nf_def_var (ncid,'pow_ei_glf',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('glf exchange power MW'),
     &          'glf exchange power MW')
      j=j+1
      oneval(j) = nf_def_var (ncid,'exch_glf',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'units',
     &           len('glf exchange heating MW/m^2'),
     &          'glf exchang heating (MW/m^2')
      j=j+1
      oneval(j) = nf_def_var (ncid,'S1',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('S1'), 'S1')
      j=j+1
      oneval(j) = nf_def_var (ncid,'S2',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('S2'), 'S2')
      j=j+1
      oneval(j) = nf_def_var (ncid,'S3',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('S3'), 'S3')
      j=j+1
      oneval(j) = nf_def_var (ncid,'S4',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('S4'), 'S4')
      j=j+1
      oneval(j) = nf_def_var (ncid,'S5',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('S5'), 'S5')
      j=j+1
      oneval(j) = nf_def_var (ncid,'INTEGRAL_LHS1',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('INTEGRAL_LHS1'), 'INTEGRAL_LHS1')
      j=j+1
      oneval(j) = nf_def_var (ncid,'INTEGRAL_LHS2',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('INTEGRAL_LHS2'), 'INTEGRAL_LHS2')
      j=j+1
      oneval(j) = nf_def_var (ncid,'INTEGRAL_LHS3',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('INTEGRAL_LHS3'), 'INTEGRAL_LHS3')
      j=j+1
      oneval(j) = nf_def_var (ncid,'INTEGRAL_LHS4',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('INTEGRAL_LHS4'), 'INTEGRAL_LHS4')
      j=j+1
      oneval(j) = nf_def_var (ncid,'INTEGRAL_LHS5',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('INTEGRAL_LHS5'), 'INTEGRAL_LHS5')
      j=j+1
      oneval(j) = nf_def_var (ncid,'INTEGRAL_RHS1',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('INTEGRAL_RHS1'), 'INTEGRAL_RHS1')
      j=j+1
      oneval(j) = nf_def_var (ncid,'INTEGRAL_RHS2',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('INTEGRAL_RHS2'), 'INTEGRAL_RHS2')
      j=j+1
      oneval(j) = nf_def_var (ncid,'INTEGRAL_RHS3',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('INTEGRAL_RHS3'), 'INTEGRAL_RHS3')
      j=j+1
      oneval(j) = nf_def_var (ncid,'INTEGRAL_RHS4',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('INTEGRAL_RHS4'), 'INTEGRAL_RHS4')
      j=j+1
      oneval(j) = nf_def_var (ncid,'INTEGRAL_RHS5',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('INTEGRAL_RHS5'), 'INTEGRAL_RHS5')
      j=j+1
      oneval(j) = nf_def_var (ncid,'DIFF',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('DIFF'), 'DIFF')
      j=j+1
      oneval(j) = nf_def_var (ncid,'CONV',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('CONV'), 'CONV')
      j=j+1
      oneval(j) = nf_def_var (ncid,'NU',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('NU'), 'NU')
      j=j+1
      oneval(j) = nf_def_var (ncid,'STIFFee',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('STIFFee'), 'STIFFee')
      j=j+1
      oneval(j) = nf_def_var (ncid,'STIFFie',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('STIFFie'), 'STIFFie')
      j=j+1
      oneval(j) = nf_def_var (ncid,'STIFFei',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('STIFFei'), 'STIFFei')
      j=j+1
      oneval(j) = nf_def_var (ncid,'STIFFii',nf_double,1,
     &            rdim,rhid1(j))
      oneval(j)=nf_put_att_text(ncid,rhid1(j),'def',
     &           len('STIFFii'), 'STIFFii')
c
c End definitions
c
      val = nf_enddef(ncid)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c Write scalars
c
      val = nf_put_var_int(ncid, rhidd, 101)
      j=1
      cstart = 1
c      ccount = txlen
      ccount = len(shot)
c      write(*,*) 'txval = ',txval,txid
      scalval(j) = nf_put_vara_text(ncid,txid,cstart,ccount,shot)
      j=j+1
      scalval(j) = nf_put_var_int(ncid, rhid(j), jmaxm)
      j=j+1
      scalval(j) = nf_put_var_int(ncid, rhid(j), jout_m)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), xp_time)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), rmajor_exp)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), rgeom_d)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), amin_d)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), bt_exp)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), tocur_d/1.D6)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), totohm_d/1.D6)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), totboot_d/1.D6)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), totbeam_d/1.D6)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), totrf_d/1.D6)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), neline_exp)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), volavene_exp)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), niline_exp)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), elonga_exp)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), deltaa_exp)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), wstr_exp)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), wstr_ped_exp)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), wstr_inc_exp)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), wstr_core_exp)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), gtaut_exp)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), tauescale_exp(1))
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), tauescale_exp(2))
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), tauescale_exp(3))
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), powe_exp(jmaxm))
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), powi_exp(jmaxm))
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), 
     &             (powe_beam_exp(jmaxm) + powe_rf_exp(jmaxm)) )
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), 
     &             (powi_beam_exp(jmaxm) + powi_rf_exp(jmaxm)) )
      j=j+1
      scalval(j) = nf_put_var_double(ncid,rhid(j),
     &             powe_oh_exp(jmaxm))
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), 
     &             powe_rad_exp(jmaxm))
c
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), wstr_m)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), wstr_inc_m)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), wstr_core_m)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), gtaut_m)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), n_t_tau_m)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), volavete_m)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), volaveti_m)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), volaveniti_m)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), ne_m(jout_m))
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), te_m(jout_m))
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), ti_m(jout_m))
      if(abs(irad).gt.0) then
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j), -Pradb_tot)
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j), -Prads_tot)
      endif
      if(ialpha.gt.0) then
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j), Pfusion)
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j), Pfusion/5.0 -
     &            Pradb_tot-Prads_tot+
     &            pbescale*powe_beam_exp(jmaxm)+
     &            prfscale*powe_rf_exp(jmaxm)+
     &            pbiscale*powi_beam_exp(jmaxm)+
     &            prfscale*powi_rf_exp(jmaxm)+powe_oh_exp(jmaxm) )
      endif
      if(abs(iohm).gt.0) then
        j=j+1
        scalval(j) = nf_put_var_double(ncid, rhid(j), Pohmic_tot)
      endif
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), norm_s)
      j=j+1
      scalval(j) = nf_put_var_double(ncid, rhid(j), norm_dt)
c
c Write 1d vars
c
      j=1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), rho)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), te_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), ti_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), tfast_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), ne_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), ni_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), nz_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), nfast_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), zeff_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), curden_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), curboot_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), curbeam_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), currf_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), curohm_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), q_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), shat_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), alpha_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), betae_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), ptot_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), torque_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),1.D-2*vphi_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),1.D-2*vphie_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),1.D-2*vphiz_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), 1.D-5*angrot_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),1.D-2*vpar_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),1.D-2*vpare_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),1.D-2*vparz_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), 1.D-2*vexb_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), 
     >   1.D-2*vdia_exp(1,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), 
     >   1.D-2*vdia_exp(2,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), 
     >   1.D-2*vdia_exp(3,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),
     >    1.D-2*vneo_exp(1,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),
     >    1.D-2*vneo_exp(2,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),
     >    1.D-2*vneo_exp(3,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),
     >    mach_exp(1,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),
     >    mach_exp(2,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),
     >    mach_exp(3,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), zpne_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), zpni_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), zpte_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), zpti_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), zpne2_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), zpte2_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), zpti2_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), flow_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powe_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powi_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), stress_tor_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), pow_ei_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powe_beam_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powi_beam_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powe_rf_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powi_rf_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powe_oh_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powe_rad_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powi_ion_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powi_cx_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powe_wdot_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powi_wdot_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), flow_wall_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), flow_beam_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), sion_exp(:,1))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), srecom_exp(:,1))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), scx_exp(:,1))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), sbcx_exp(:,1))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), s_exp(:,1))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), sdot_exp(:,1))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), enn_exp(:,1))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), ennw_exp(:,1))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), ennv_exp(:,1))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), volsn_exp(:,1))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), sbeame_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), sbeam_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), qbeame_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), qbeami_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), qrfe_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), qrfi_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), qohm_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), qrad_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), qione_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), qioni_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), rmin_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), rmaj_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), elong_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), delta_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), gradrho_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), gradrhosq_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), f_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), g_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), h_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), volume_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), drhodr)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), theta_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), zptheta_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), a_pol)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), a_tor)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), c_par)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), c_per)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), c_tor)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), 
     &            vprime(1:mxgrid+1,1))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), 
     &            vprime(1:mxgrid+1,2))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), sfarea_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), cxarea_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), geofac)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), csda_exp)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), rhosda_exp)
c
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), ti_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), te_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), ne_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), ni_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), nz_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), 1.D-2*vphi_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), 1.D-2*vphie_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), 1.D-2*vphiz_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), 1.D-2*vpar_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), 1.D-2*vpare_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), 1.D-2*vparz_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),1.D-2*vpol_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), xnu_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), nuei_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), nu_pol_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), kpol_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), alpha_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), betae_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), zpne_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), zpni_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), zpte_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), zpti_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), interchange_DR_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), anrate_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), anfreq_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), egamma_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), gamma_p_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), 1.D-2*vexb_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),
     >   1.D-2*vdia_m(1,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),
     >   1.D-2*vdia_m(2,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),
     >   1.D-2*vdia_m(3,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), 
     >  1.D-2*vneo_m(1,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), 
     >  1.D-2*vneo_m(2,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), 
     >  1.D-2*vneo_m(3,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), 
     >  mach_m(1,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), 
     >  mach_m(2,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), 
     >  mach_m(3,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), vstar_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), vepol_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), vetor_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid,rhid1(j),diffgb_m*cgyrobohm_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid,rhid1(j),chiigb_m*cgyrobohm_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid,rhid1(j),chiegb_m*cgyrobohm_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid,rhid1(j),
     &            chie_e_gb_m*cgyrobohm_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid,rhid1(j),chiineo_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid,rhid1(j),
     &            etagb_phi_m*cgyrobohm_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid,rhid1(j),cgyrobohm_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), flowe_glf)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), flowi_glf)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), flowz_glf)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powe_glf)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powi_glf)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powz_glf)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), stress_tor_i_glf)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), stress_tor_z_glf)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), stress_par_i_glf)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), stress_par_z_glf)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), flowe_neo)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), flowi_neo)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), flowz_neo)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powe_neo)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powi_neo)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powz_neo)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), stress_tor_i_neo)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), stress_tor_z_neo)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), stress_par_i_neo)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), stress_par_z_neo)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), flowe_adhoc)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), flowi_adhoc)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), flowz_adhoc)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powe_adhoc)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powi_adhoc)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powz_adhoc)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), stress_tor_i_adhoc)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), stress_tor_z_adhoc)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), stress_par_i_adhoc)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), stress_par_z_adhoc)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), flowe_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), flowi_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), flowz_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powe_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powi_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), powz_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), stress_tor_i_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), stress_tor_z_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), stress_par_i_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), stress_par_z_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), stress_par_cor_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), Pi_alpha)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), Pe_alpha)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), Piaux)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), Peaux)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), Pohpro)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), pow_ei_cor_m)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), pow_ei_glf)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), exch_glf)
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), S(1,1:mxgrid+1))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), S(2,1:mxgrid+1))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), S(3,1:mxgrid+1))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), S(4,1:mxgrid+1))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j), S(5,1:mxgrid+1))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),
     >  INTEGRAL_LHS(1,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),
     >  INTEGRAL_LHS(2,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),
     >  INTEGRAL_LHS(3,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),
     >  INTEGRAL_LHS(4,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),
     >  INTEGRAL_LHS(5,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),
     >  INTEGRAL_RHS(1,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),
     >  INTEGRAL_RHS(2,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),
     >  INTEGRAL_RHS(3,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),
     >  INTEGRAL_RHS(4,0:mxgrid))
      j=j+1
      oneval(j) = nf_put_var_double(ncid, rhid1(j),
     >  INTEGRAL_RHS(5,0:mxgrid))
      j=j+1
      oneval(j)=nf_put_var_double(ncid, rhid1(j),diff(1,1,1:mxgrid+1))
      j=j+1
      oneval(j)=nf_put_var_double(ncid, rhid1(j),conv3(1,1,1:mxgrid+1))
      j=j+1
      oneval(j)=nf_put_var_double(ncid, rhid1(j),nu(1,1,1:mxgrid+1))
      j=j+1
      oneval(j)=nf_put_var_double(ncid, rhid1(j),stiff(1,1,1:mxgrid+1))
      j=j+1
      oneval(j)=nf_put_var_double(ncid, rhid1(j),stiff(2,1,1:mxgrid+1))
      j=j+1
      oneval(j)=nf_put_var_double(ncid, rhid1(j),stiff(1,2,1:mxgrid+1))
      j=j+1
      oneval(j)=nf_put_var_double(ncid, rhid1(j),stiff(2,2,1:mxgrid+1))
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c

c
c      dum = nf_def_dim(ncid, 'len1', 6, lenid)
c      if (dum .ne. nf_noerr) call handle_err(dum)
c
      call ncclos(ncid, rcode)
c
      endif
c
      end
