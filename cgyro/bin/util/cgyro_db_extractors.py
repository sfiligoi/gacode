#!/usr/bin/env python
"""
CGYRO AI Database Extraction Utilities

Provides functions to extract metadata, species data, flux data, and provenance
information from CGYRO simulations for storage in the AI database.

Based on patterns from existing cgyro_db tool.
"""

import os
import json
import hashlib
import numpy as np
from gacodefuncs import time_average, time_index


def extract_metadata(sim):
    """
    Extract grid and equilibrium metadata from cgyrodata object

    Parameters:
    -----------
    sim : cgyrodata
        CGYRO simulation data object

    Returns:
    --------
    dict : Metadata dictionary with grid and equilibrium parameters
    """
    metadata = {
        # Grid parameters
        'n_n': int(sim.n_n),
        'n_radial': int(sim.n_radial),
        'n_theta': int(sim.n_theta),
        'n_species': int(sim.n_species),
        'n_field': int(sim.n_field),
        'n_energy': int(sim.n_energy),
        'n_xi': int(sim.n_xi),

        # Equilibrium parameters - geometry
        'rmin': float(sim.rmin),
        'rmaj': float(sim.rmaj),
        'q': float(sim.q),
        'shear': float(sim.shear),
        'shift': float(sim.shift),
        'kappa': float(sim.kappa),
        's_kappa': float(sim.s_kappa),
        'delta': float(sim.delta),
        's_delta': float(sim.s_delta),
        'zeta': float(sim.zeta),
        's_zeta': float(sim.s_zeta),
        'zmag': float(sim.zmag),
        'dzmag': float(sim.dzmag),

        # Resolution parameters
        'ky0': float(sim.ky0),

        # Plasma parameters
        'betae_unit': float(sim.betae_unit),
        'beta_star': float(sim.beta_star),
        'gamma_e': float(sim.gamma_e),
        'gamma_p': float(sim.gamma_p),
        'mach': float(sim.mach),
    }

    return metadata


def extract_species(sim):
    """
    Extract species arrays from cgyrodata object

    Parameters:
    -----------
    sim : cgyrodata
        CGYRO simulation data object

    Returns:
    --------
    dict : Species data with arrays converted to lists for JSON serialization
    """
    species_data = {
        'z': sim.z.tolist(),
        'mass': sim.mass.tolist(),
        'dens': sim.dens.tolist(),
        'temp': sim.temp.tolist(),
        'dlnndr': sim.dlnndr.tolist(),
        'dlntdr': sim.dlntdr.tolist(),
        'nu': sim.nu.tolist(),
    }

    return species_data


def extract_flux(sim, time_window='0.5'):
    """
    Extract and time-average flux data from cgyrodata object

    Uses the same method as cgyro_db tool:
    - Sum over fields (axis 2) and radial modes (axis 3)
    - Time average over specified window (default: last 50%)

    Parameters:
    -----------
    sim : cgyrodata
        CGYRO simulation data object
    time_window : str, optional
        Time window for averaging. Default '0.5' means last 50% of simulation.
        Can also be 't_min,t_max' for explicit range.

    Returns:
    --------
    list : List of flux dictionaries, one per species
    """
    # Load flux data
    sim.getflux(total=True)

    # Determine time averaging window
    imin,imax = time_index(sim.t,time_window)

    fluxes = []
    for ispec in range(sim.n_species):
        # Time-average fluxes
        fluxn = time_average(sim.fluxtot[ispec,0,:],sim.t,imin,imax)
        fluxe = time_average(sim.fluxtot[ispec,1,:],sim.t,imin,imax)
        fluxv = time_average(sim.fluxtot[ispec,2,:],sim.t,imin,imax)

        flux_dict = {
            'species': int(ispec),
            'fluxn': fluxn.tolist(),
            'fluxe': fluxe.tolist(),
            'fluxv': fluxv.tolist(),
            'time_window': time_window,
            'time_max': float(sim.t[-1])

        }
        fluxes.append(flux_dict)

    return fluxes


def compute_file_hash(filepath, algorithm='md5'):
    """
    Compute hash of a file for change detection

    Parameters:
    -----------
    filepath : str
        Path to file
    algorithm : str, optional
        Hash algorithm ('md5' or 'sha256'). Default 'md5'.

    Returns:
    --------
    str : Hexadecimal hash string, or None if file doesn't exist
    """
    if not os.path.exists(filepath):
        return None

    if algorithm == 'md5':
        h = hashlib.md5()
    elif algorithm == 'sha256':
        h = hashlib.sha256()
    else:
        raise ValueError(f"Unsupported hash algorithm: {algorithm}")

    try:
        with open(filepath, 'rb') as f:
            # Read in chunks for large files
            for chunk in iter(lambda: f.read(4096), b""):
                h.update(chunk)
        return h.hexdigest()
    except Exception as e:
        print(f"Warning: Could not hash {filepath}: {e}")
        return None


def compute_input_hash(case_path):
    """
    Compute unique hash of input configuration

    This creates a unique identifier for the input parameters.
    Uses input.cgyro file, which contains all simulation parameters.

    Parameters:
    -----------
    case_path : str
        Path to simulation directory

    Returns:
    --------
    str : MD5 hash of input.cgyro, or None if file doesn't exist
    """
    input_file = os.path.join(case_path, 'input.cgyro')
    return compute_file_hash(input_file, algorithm='md5')


def extract_provenance(case_path, t_final, i_final):
    """
    Extract simulation provenance information

    Parameters:
    -----------
    case_path : str
        Path to simulation directory
    t_final : float
        Final simulation time (sim.t[-1])
    i_final : int
        Number of time steps (len(sim.t))

    Returns:
    --------
    dict : Provenance data including timestamps, simulation state, file hashes
    """
    time_file = os.path.join(case_path, 'out.cgyro.time')

    # Get file modification times
    created_timestamp = None
    modified_timestamp = None

    if os.path.exists(time_file):
        modified_timestamp = os.path.getmtime(time_file)
        # Use creation time if available, else use modified time
        try:
            created_timestamp = os.path.getctime(time_file)
        except:
            created_timestamp = modified_timestamp

    # Count restart files
    restart_count = 0
    try:
        restart_files = [f for f in os.listdir(case_path) if f.startswith('bin.cgyro.restart')]
        restart_count = len(restart_files)
    except:
        pass

    # Compute hash of input configuration (unique identifier)
    input_hash = compute_input_hash(case_path)

    # Compute hashes of key output files for change detection
    file_hashes = {}
    key_files = ['bin.cgyro.freq', 'out.cgyro.time', 'out.cgyro.prec']
    for fname in key_files:
        fpath = os.path.join(case_path, fname)
        h = compute_file_hash(fpath)
        if h is not None:
            file_hashes[fname] = h

    provenance = {
        't_final': t_final,
        'i_final': i_final,
        'created_timestamp': created_timestamp,
        'modified_timestamp': modified_timestamp,
        'n_restarts': restart_count,
        'input_hash': input_hash,  # Unique identifier for input configuration
        'file_hashes': json.dumps(file_hashes),
    }

    return provenance


def get_git_version(directory):
    """
    Get git version/commit hash for a directory

    Parameters:
    -----------
    directory : str
        Directory path

    Returns:
    --------
    str : Git commit hash, or None if not in git repo
    """
    import subprocess
    try:
        result = subprocess.run(
            ['git', 'rev-parse', 'HEAD'],
            cwd=directory,
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except:
        pass

    return None


def is_cgyro_case(directory):
    """
    Check if a directory is a valid CGYRO simulation case

    Parameters:
    -----------
    directory : str
        Directory path to check

    Returns:
    --------
    bool : True if directory contains CGYRO output files
    """
    required_files = ['out.cgyro.grids']
    for fname in required_files:
        if not os.path.exists(os.path.join(directory, fname)):
            return False
    return True
