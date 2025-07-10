
import os
import numpy as np
import sys
import time

class tglfdata:
    """A class for management of TGLF output data.

    This class provides methods to extract and manage data from TGLF output files.
    It reads various spectra and total flux information related to the simulation.

    Attributes:
        silent (bool): If True, suppresses output messages.
        dir (str): The directory where the TGLF output files are located.
        nmodes (int): The number of modes from the simulation.
        nspecies (int): The number of species in the simulation.
        nfield (int): The number of fields in the simulation.
        nky (int): The number of ky values.
        ntype (int): The type of simulation.
        gamma (list): List of eigenvalue growth rates for each mode.
        freq (list): List of eigenvalue frequencies for each mode.
        ky (numpy.ndarray): Array of ky spectrum values.
        flux_spectrum_particle (numpy.ndarray): Particle flux spectrum.
        flux_spectrum_energy (numpy.ndarray): Energy flux spectrum.
        flux_spectrum_toroidal_stress (numpy.ndarray): Toroidal stress flux spectrum.
        flux_spectrum_parallel_stress (numpy.ndarray): Parallel stress flux spectrum.
        flux_spectrum_exchange (numpy.ndarray): Exchange flux spectrum.
        particle_GB (list): Total particle flux for each species.
        energy_GB (list): Total energy flux for each species.
        momentum_GB (list): Total momentum flux for each species.

    Example Usage:
        >>> from pygacode.tglf.data import tglfdata
        >>> sim1 = tglfdata('example_directory')
    """
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
    def extract(self,f):

      
        
      if os.path.isfile(self.dir+'out'+f):
         filename = self.dir+'out'+f
         with open(filename, 'r') as f:
            data = f.readlines()

         
      else:
         print("File not found")
         # File not found
        
         data = []

     

      return data
        
    #-------------------------------------------------------------------------#

    
    def getdata(self):
      data = self.extract('.tglf.QL_flux_spectrum')
      
      
      
      ntype, nspecies, nfield, nky, nmodes = list(map(int, data[3].strip().split()))
      self.nmodes = nmodes
      self.nspecies=nspecies
      self.nfield=nfield
      self.nky=nky
      self.ntype = ntype
    #-------------------------------------------------------------------------#
    def get_eigenvalue_spectrum(self):
        nmodes = self.nmodes
      
        eigen=self.extract('.tglf.eigenvalue_spectrum')
        ky = self.extract('.tglf.ky_spectrum')

        eigen = ''.join(eigen[2:]).split()
        ky = ''.join(ky[2:]).split()
        
        gamma = []
        freq = []
        for k in range(self.nmodes):
            gamma.append(np.array(eigen[2 * k :: nmodes * 2], dtype=float))
            freq.append(np.array(eigen[2 * k + 1 :: nmodes * 2], dtype=float))
        
        self.gamma=gamma      
        self.freq=freq
        self.ky=np.array(ky, dtype=float)
       
    def get_flux_spectrum(self):
        
        
        flux=self.extract('.tglf.sum_flux_spectrum')
        
        
        flux_spectrum_particle = np.zeros((self.nky, self.nspecies))
        flux_spectrum_energy = np.zeros((self.nky, self.nspecies))
        flux_spectrum_toroidal_stress = np.zeros((self.nky, self.nspecies))
        flux_spectrum_parallel_stress =np.zeros((self.nky, self.nspecies))
        flux_spectrum_exchange =np.zeros((self.nky, self.nspecies))
       
        for n in range(self.nspecies):
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
       
       flux_GB=self.extract('.tglf.run')
       flux_GB = ''.join(flux_GB[5:]).split()
       particle_GB=[]
       energy_GB=[]
       momentum_GB=[]
   
       for ns in range(self.nspecies):
          particle_GB.append(float(flux_GB[ns*6+1]))
          energy_GB.append(float(flux_GB[ns*6+2]))
          momentum_GB.append(float(flux_GB[ns*6+4]))
       self.particle_GB = particle_GB
       self.energy_GB = energy_GB
       self.momentum_GB = momentum_GB
       