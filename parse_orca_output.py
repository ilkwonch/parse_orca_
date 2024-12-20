import sys
import os
import numpy as np  # Import NumPy for numerical operations
from ase import units

# https://github.com/psi4/psi4/blob/103e1d6b5bf8cb12005fbfc9d5b695d6866830ea/psi4/include/psi4/physconst.h#L402
bohr2angstroms = 0.52917721067
hartree2eV = units.Hartree

from parser import (
    total_energy,
    gradient,
    hirshfeld_chg,
    dipole,
    ParseError,
    dftd_energy,
    dftd_gradient,
)

def parse_orca_output(filename):
    with open(filename, 'r') as f:
        results = {}

        # Parse total energy
        try:
            f.seek(0)
            results.update(total_energy(f))
        except ParseError as e:
            results['total_energy'] = None

        # Parse forces (gradients)
        try:
            f.seek(0)
            results.update(gradient(f))
        except ParseError as e:
            results['gradients'] = None

        # Parse Hirshfeld charges
        try:
            f.seek(0)
            results.update(hirshfeld_chg(f))
        except ParseError as e:
            results['hirshfeld_charges'] = None

        # Parse Dipole
        try:
            f.seek(0)
            results.update(dipole(f))
        except ParseError as e:
            results['dipole'] = None

        # Parse DFT-D energy
        try:
            f.seek(0)
            results.update(dftd_energy(f))
        except ParseError as e:
            results['dftd_energy'] = None

        # Parse DFT-D gradient
        try:
            f.seek(0)
            results.update(dftd_gradient(f))
        except ParseError as e:
            results['dftd_gradient'] = None

    return results

def process_parsed_data(data):
    # Ensure all entries are converted to NumPy arrays of dtype float64
    parsed_data = {key: np.array(val, dtype=np.float64) if val is not None else None for key, val in data.items()}

    # Compute derived properties
    energy = None
    forces = None
    charge = None
    dipole = None

    # Compute energy difference
    if parsed_data['total_energy'] is not None and parsed_data['dftd_energy'] is not None:
        energy = (parsed_data['total_energy'] - parsed_data['dftd_energy']) * hartree2eV

    # Compute forces difference and reshape to (n, 3)
    if parsed_data['gradients'] is not None and parsed_data['dftd_gradients'] is not None:
        forces = parsed_data['gradients'] - parsed_data['dftd_gradients']
        forces = - (forces.reshape((-1, 3)) * hartree2eV / bohr2angstroms)

    # Assign charge and dipole
    if parsed_data['hirshfeld_charges'] is not None:
        charge = parsed_data['hirshfeld_charges']
    if parsed_data['dipole'] is not None:
        dipole = (parsed_data['dipole'] * bohr2angstroms)
        

    return {
        'energy': energy,
        'forces': forces,
        'charge': charge,
        'dipole': dipole,
    }
