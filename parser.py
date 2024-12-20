from util import goto, iter_until, ParseError
import numpy as np


# -------------------------   --------------------
# FINAL SINGLE POINT ENERGY      -455.766609105682
# -------------------------   --------------------
def total_energy(f):
    try:
        line = goto(f, 'FINAL SINGLE POINT ENERGY')
        result = {'total_energy': float(line.split()[-1])}
        return result
    except:
        raise ParseError(f'{line}')


# ----------------
# TOTAL SCF ENERGY
# ----------------
#
# Total Energy       :         -455.76660911 Eh          -12402.03994 eV
#
def scf_energy(f):
    goto(f, 'TOTAL SCF ENERGY')
    try:
        line = goto(f, 'Total Energy')
        result = {'scf_energy': float(line.split()[3])}
        return result
    except:
        raise ParseError(f'{line}')


# -------------------------
# THE POLARIZABILITY TENSOR
# -------------------------
#
# The raw cartesian tensor (atomic units):
#   38.34080      0.62734     -1.56520
#    0.62734     34.89230     -2.20299
#   -1.56520     -2.20299     28.09135
# diagonalized tensor:
#   27.27992     35.21723     38.82730
#
#    0.12018      0.29760     -0.94710
#    0.26682     -0.92859     -0.25792
#    0.95622      0.22171      0.19100
#
# Isotropic polarizability :  33.77482
def polar(f):
    goto(f, 'THE POLARIZABILITY TENSOR')
    try:
        line = goto(f, 'Isotropic polarizability')
        result = {'alpha': float(line.split()[3])}
        return result
    except:
        raise ParseError(f'{line}')


# -------------
# DIPOLE MOMENT
# -------------
#                                X             Y             Z
# Electronic contribution:     -0.01701      -2.20693       1.57967
# Nuclear contribution   :      0.06820       1.14335      -0.02857
#                        -----------------------------------------
# Total Dipole Moment    :      0.05119      -1.06357       1.55110
#                        -----------------------------------------
# Magnitude (a.u.)       :      1.88141
# Magnitude (Debye)      :      4.78217
def dipole(f):
    goto(f, 'DIPOLE MOMENT')
    try:
        line = goto(f, 'Total Dipole Moment')
        result = {'dipole': [float(i) for i in line.split()[4:7]]}
        return result
    except:
        raise ParseError(f'{line}')


def mp2_dipole(f):
    goto(f, 'MP2 RELAXED DENSITY')
    d = dipole(f)
    return {'mp2_dipole': d['dipole']}


def mp2_quadrupole(f):
    goto(f, 'MP2 RELAXED DENSITY')
    d = quadrupole(f)
    return {'mp2_quadrupole': d['quadrupole']}

# ------------------
# CARTESIAN GRADIENT
# ------------------
#
#   1   C   :   -0.080651102    0.000000348   -0.000000057
#   2   C   :    0.080646539    0.000002113   -0.000000084
#
def gradient(f):
    line = goto(f, 'CARTESIAN GRADIENT')
    try:
        next(f)
        next(f)
        g = list()
        for line in iter_until(f, '^\n'):
            g.extend(float(x) for x in line.split()[-3:])
        result = dict(gradients=g)
        return result
    except:
        raise ParseError(f'{line}')


#DFT components:
#N(Alpha)           :        9.000000193095 electrons
#N(Beta)            :        9.000000193095 electrons
#N(Total)           :       18.000000386190 electrons
#E(X)               :      -19.484899019128 Eh
#E(C)               :       -0.731059133957 Eh
#NL Energy, E(C,NL) :        0.075791465926 Eh
#E(XC)              :      -20.140166687159 Eh
#DFET-embed. en.    :        0.000000000000 Eh
def dft_nl(f):
    line = goto(f, '^NL Energy, E\(C,NL\)')
    try:
        v = float(line.split()[-2])
        return {'nl_energy': v}
    except:
        try:
            f.seek(0)
            line = goto(f, '^NL    Energy:')
            v = float(line.split()[-1])
            return {'nl_energy': v}
        except:
            raise ParseError(f'{line}')


#------------------------
#QUADRUPOLE MOMENT (A.U.)
#------------------------
#
#                XX           YY           ZZ           XY           XZ           YZ
#  NUC       218.21144    377.01684    168.02767     88.55736     78.04989     60.34841
#  EL       -256.01416   -418.80024   -206.36632    -83.14691    -77.99073    -61.51749
#  TOT       -37.80272    -41.78339    -38.33865      5.41045      0.05916     -1.16908 (a.u.)
#            -50.84595    -56.20009    -51.56679      7.27724      0.07958     -1.57246 (Buckingham)
def quadrupole(f):
    goto(f, '^QUADRUPOLE MOMENT')
    try:
        line = goto(f, '^\s+TOT')
        result = {'quadrupole': [float(i) for i in line.split()[1:7]]}
        return result
    except:
        raise ParseError(f'{line}')


#-------------------------   ----------------
#Dispersion correction           -0.042136596
#-------------------------   ----------------
def dftd_energy(f):
    line = goto(f, '^Dispersion correction')
    try:
        result = {'dftd_energy': float(line.split()[-1])}
        return result
    except:
        raise ParseError(f'{line}')


#-------------------
#DISPERSION GRADIENT
#-------------------
#   1  :    0.000036082    0.001141919    0.000151236
# ...
#  18  :    0.000487689    0.000075138   -0.000609330
#
def dftd_gradient(f):
    line = goto(f, '^DISPERSION GRADIENT')
    try:
        next(f)
        next(f)
        g = list()
        for line in iter_until(f, '^\n'):
            g.extend(float(x) for x in line.split()[-3:])
        result = dict(dftd_gradients=g)
        return result
    except:
        raise ParseError(f'{line}')







def engrad(f):
    goto(f, 'The current gradient in Eh/bohr')
    line = None
    try:
        g = list()
        next(f)
        line = next(f)
        while not line.startswith('#'):
            g.append(float(line))
            line = next(f)
        return {'engrad': g}
    except:
        raise ParseError(f'{line}')


def fod(f):
    line = goto(f, '^N_FOD =')
    try:
        v = float(line.split()[-1])
        return {'fod': v}
    except:
        raise ParseError(f'{line}')


#SMD CDS free energy correction energy :                 3.57207     Kcal/mol
def smd_cds(f):
    line = goto(f, '^SMD CDS free energy correction energy')
    try:
        v = float(line.split()[-2])
        return {'smd_cds': v}
    except:
        raise ParseError(f'{line}')


def _readomat(f, mat):
    while True:
        line = f.readline()
        if line.startswith('$'+mat):
            dim = [int(i) for i in f.readline().split()]
            # guess right array dimensions
            if len(dim) == 1:
                x = f.tell()
                line1 = f.readline()
                line2 = f.readline()
                if len(line1.split()) < len(line2.split()):
                    dim = [dim[0], dim[0]]
                else:
                    dim = [dim[0], 1]
                f.seek(x)
            arr = np.empty(dim, dtype=np.float64)
            x = f.tell()
            _n = len(next(f).split())
            f.seek(x)
            for _ in range(int(np.ceil(dim[1]/_n))):
                if dim[1] == 1:
                    j = [0]
                else:
                    j = [int(x) for x in f.readline().split()]
                for i in range(dim[0]):
                    arr[i, j] = [float(x) for x in f.readline().split()[1:]]
            break
    return arr


def normal_modes(f):
    """ Parse .hess file to get normal modes
    """
    return {'normal_modes': _readomat(f, 'normal_modes').flatten().tolist()}

def vib_greq(f):
    """ Parse .hess file to get vibrational frequencies
    """
    return {'vib_greq': _readomat(f, 'vibrational_frequencies').flatten().tolist()}

def coord(f):
    """ return coord from .xyz file
    """
    lines = f.readlines()
    n = int(lines[0].rstrip().lstrip())
    coord = list()
    for i in range(n):
        coord.extend([float(x) for x in lines[i+2].split()[1:]])
    return {'coord': coord}


def hirshfeld_chg(f):
    goto(f, '^HIRSHFELD ANALYSIS')
    goto(f, '^  ATOM     CHARGE')
    chg = list()
    spin_chg = list()
    for line in f:
        p = line.split()
        if len(p) < 2:
            break
        if len(p) >= 3:
            chg.append(float(p[2]))
        if len(p) == 4:
            spin_chg.append(float(p[3]))
    ret = {'hirshfeld_charges': chg}
    if spin_chg:
        ret['hirshfeld_spin_charges'] = spin_chg
    return ret


def orb_occ_e(f):
    goto(f, '^ORBITAL ENERGIES')
    goto(f, '  NO   OCC')
    occ, ene = [], []
    for line in f:
        ll = line.split()
        if len(ll) < 4:
            break
        occ.append(float(ll[1]))
        ene.append(float(ll[3]))
    return({'orbital_occupation': occ, 'orbital_energy': ene})

def gap(f):
    d = orb_occ_e(f)
    ene = d['orbital_energy']
    occ = d['orbital_occupation']
    homo, lumo = None, None
    for e, c in zip(ene, occ):
        if c == 2.0:
            homo = e
        elif c == 0.0:
            lumo = e
            break
    assert homo is not None and lumo is not None
    gap = lumo - homo
    return {'gap': gap}

