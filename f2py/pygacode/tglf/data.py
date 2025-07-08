
import os
import numpy as np
import sys
import time

class tglfdata:
    """A class for management of TGLF output data.

    Data:

    Example Usage:
        >>> from pygacode.tglf.data import tglfdata
        >>> sim1 = tglfdata('example_directory')
    """

    #-------------------------------------------------------------------------#
    # Methods

    def __init__(self,sim_directory,silent=False,fast=False):

      self.silent = silent
      self.dir = sim_directory

      
      
      self.getdata()
      self.get_eigenvalue_spectrum()
      self.get_flux_spectrum()
      self.get_total_flux()

    # standard routine to read binary or ASCII data
    def extract(self,f, cmplx=False):

      self.BYTE = 'float64'
      self.CBYTE = 'complex64'
      
      if os.path.isfile(self.dir+'bin'+f):
         
         fmt = 'bin'
         if cmplx:
            dtype = self.CBYTE
         else:
            dtype = self.BYTE   
         
         data = np.fromfile(self.dir+'bin'+f,dtype=dtype)
      elif os.path.isfile(self.dir+'out'+f):
         filename = self.dir+'out'+f
         with open(filename, 'r') as f:
            data = f.readlines()

         fmt = 'out'
         #data = np.fromfile(self.dir+'out'+f,dtype='float')
      else:
         print("File not found")
         # File not found
         fmt = 'null'
         data = []

     

      return fmt,data
        
    #-------------------------------------------------------------------------#

    
    def getdata(self):
      fmt, data = self.extract('.tglf.QL_flux_spectrum')
      
      
      
      ntype, nspecies, nfield, nky, nmodes = list(map(int, data[3].strip().split()))
      self.nmodes = nmodes
      self.nspecies=nspecies
      self.nfield=nfield
      self.nky=nky
      self.ntype = ntype
    #-------------------------------------------------------------------------#
    def get_eigenvalue_spectrum(self):
        nmodes = self.nmodes
      
        fmt,eigen=self.extract('.tglf.eigenvalue_spectrum')
        fmt, ky = self.extract('.tglf.ky_spectrum')

        eigen = ''.join(eigen[2:]).split()
        ky = ''.join(ky[2:]).split()
        
        gamma = []
        freq = []
        for k in range(self.nmodes):
            gamma.append(np.array(eigen[2 * k :: nmodes * 2], dtype=float))
            freq.append(np.array(eigen[2 * k + 1 :: nmodes * 2], dtype=float))
        
        self.gamma=gamma      
        self.freq=freq
        self.ky=ky
       
    def get_flux_spectrum(self):
        nmodes = self.nmodes
        
        fmt,flux=self.extract('.tglf.sum_flux_spectrum')
        
        
        flux_spectrum_particle = np.zeros((self.nky, self.nmodes))
        flux_spectrum_energy = np.zeros((self.nky, self.nmodes))
        flux_spectrum_toroidal_stress = np.zeros((self.nky, self.nmodes))
        flux_spectrum_parallel_stress =np.zeros((self.nky, self.nmodes))
        flux_spectrum_exchange =np.zeros((self.nky, self.nmodes))
       
        for n in range(self.nmodes):
            data = flux[self.nky*n+2:self.nky*(n+1)+2] 
            
            if any([x.startswith((" s")) for x in data]):
               
               
                data = flux[self.nky*n+4:self.nky*(n+1)+4] 
                
           
            for k in range(0, self.nky*5, 5):
           
                
                split_data = [element for line in data for element in line.split()]
                
                 
                index = int(k/5)
                
                flux_spectrum_particle[index,n]= split_data[k]
                flux_spectrum_energy[index,n]= split_data[k+1]
                flux_spectrum_toroidal_stress[index,n]= split_data[k+2]
                flux_spectrum_parallel_stress[index,n]= split_data[k+3]
                flux_spectrum_exchange[index,n]= split_data[k+4]
               


       # types = ['particle', 'energy', 'toroidal stress', 'parallel stress', 'exchange']
        self.flux_spectrum_particle = flux_spectrum_particle
        self.flux_spectrum_energy = flux_spectrum_energy
        self.flux_spectrum_toroidal_stress = flux_spectrum_toroidal_stress
        self.flux_spectrum_parallel_stress = flux_spectrum_parallel_stress
        self.flux_spectrum_exchange = flux_spectrum_exchange
    
    
    def get_total_flux(self):
       print("")