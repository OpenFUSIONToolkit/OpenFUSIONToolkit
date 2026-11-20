#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
"""Plasma initiation (burn-through) model for the Breakdown_Startup machine and circuit models, after
DYON: H.-T. Kim, PhD thesis, Imperial College London (2013); Kim et al., Nucl. Fusion 60 126049
(2020), 62 126012 (2022).

    AtomicSpecies           ADAS iso-nuclear rates and powers of one element, per charge state
    fractional_abundance    charge-state balance with a transport loss (coronal for tau -> inf)

Atomic data: OPEN-ADAS adf11 files (H.P. Summers, The ADAS User Manual), downloaded on first use
into ADAS_DIR ($OFT_ADAS_DIR, default ~/.cache/oft_breakdown_startup/adas). They are not shipped:
the OPEN-ADAS terms forbid redistributing the files with a modelling code.

SI units (m, s, m^-3, m^3/s, W m^3) except temperatures [eV].
"""
import os
import urllib.request

import numpy as np
from scipy.interpolate import RegularGridInterpolator

ADAS_DIR = os.environ.get('OFT_ADAS_DIR', os.path.join(os.path.expanduser('~'), '.cache', 'oft_breakdown_startup', 'adas'))
ADAS_URL = 'https://open.adas.ac.uk/download/adf11/{cls}/{name}.dat'

# adf11 files per element: scd ionization, acd recombination, ccd charge-exchange recombination by
# neutral hydrogen, plt line power, prb recombination + bremsstrahlung power. 96 (the GCR project)
# everywhere, as DYON (K20 sec. 2.3); there is no ccd96 for O and Be, so ccd89 there.
ADAS_FILES = {
    'h': {'scd': 'scd96_h', 'acd': 'acd96_h', 'ccd': 'ccd96_h', 'plt': 'plt96_h', 'prb': 'prb96_h'},
    'd': {'scd': 'scd96_h', 'acd': 'acd96_h', 'ccd': 'ccd96_d', 'plt': 'plt96_h', 'prb': 'prb96_h'},
    'be': {'scd': 'scd96_be', 'acd': 'acd96_be', 'ccd': 'ccd89_be', 'plt': 'plt96_be', 'prb': 'prb96_be'},
    'c': {'scd': 'scd96_c', 'acd': 'acd96_c', 'ccd': 'ccd96_c', 'plt': 'plt96_c', 'prb': 'prb96_c'},
    'n': {'scd': 'scd96_n', 'acd': 'acd96_n', 'ccd': 'ccd96_n', 'plt': 'plt96_n', 'prb': 'prb96_n'},
    'o': {'scd': 'scd96_o', 'acd': 'acd96_o', 'ccd': 'ccd89_o', 'plt': 'plt96_o', 'prb': 'prb96_o'},
}
# Ionization energies of charge state z -> z+1 [eV]: NIST Atomic Spectra Database (retrieved 2026-11-19)
IONIZATION_EV = {
    'h': [13.598434599702], 'd': [13.598434599702],
    'be': [9.322699, 18.21115, 153.896205, 217.71858459],
    'c': [11.2602880, 24.383143, 47.88778, 64.49352, 392.09056, 489.99320779],
    'n': [14.53413, 29.60125, 47.4453, 77.4735, 97.8901, 552.06741, 667.0461377],
    'o': [13.618055, 35.12112, 54.93554, 77.41350, 113.8990, 138.1189, 739.32697, 871.4099138],
}
MASS_AMU = {'h': 1.00794, 'd': 2.01410, 'be': 9.01218, 'c': 12.0107, 'n': 14.0067, 'o': 15.9994}


def adas_file(name):
    """Local path of OPEN-ADAS adf11 file `name` (eg. 'scd96_c'), downloaded into ADAS_DIR if missing."""
    path = os.path.join(ADAS_DIR, name + '.dat')
    if not os.path.isfile(path):
        os.makedirs(ADAS_DIR, exist_ok=True)
        url = ADAS_URL.format(cls=name.split('_')[0], name=name)
        with urllib.request.urlopen(url, timeout=60) as r:
            text = r.read().decode('latin-1')
        if 'DOCTYPE' in text[:200] or '/' not in text.split('\n')[0]:     # OPEN-ADAS answers a missing file with a web page
            raise FileNotFoundError(f'{name}: not an adf11 file at {url}')
        with open(path, 'w') as f:
            f.write(text)
    return path


def read_adf11(path):
    """An adf11 (iso-nuclear) file: {'z_max', 'logne' [log10 m^-3], 'logT' [log10 eV], 'blocks'
    [n_block, n_T, n_ne] log10 SI (m^3/s or W m^3), 'z1' [n_block] (the file's Z1 per block)}."""
    lines = open(path, encoding='latin-1').read().split('\n')
    head = lines[0].split()
    z_max, n_ne, n_T = int(head[0]), int(head[1]), int(head[2])
    nums, blocks, z1, k = [], [], [], 2                 # line 1 is a separator
    while len(nums) < n_ne + n_T:
        nums += [float(x) for x in lines[k].split()]
        k += 1
    logne, logT = np.array(nums[:n_ne]) + 6.0, np.array(nums[n_ne:n_ne + n_T])     # cm^-3 -> m^-3
    while k < len(lines):
        line = lines[k]
        if line.startswith('C') or line.strip() == '':
            break
        if 'Z1=' in line:
            z1.append(int(line.split('Z1=')[1].split('/')[0]))
            vals = []
            k += 1
            while len(vals) < n_ne * n_T:
                vals += [float(x) for x in lines[k].split()]
                k += 1
            blocks.append(np.array(vals).reshape(n_T, n_ne) - 6.0)   # cm^3/s, W cm^3 -> m^3/s, W m^3
            continue
        k += 1
    return {'z_max': z_max, 'logne': logne, 'logT': logT, 'blocks': np.array(blocks), 'z1': np.array(z1)}


class AtomicSpecies:
    """ADAS rates and powers of one element ('h', 'd', 'be', 'c', 'n', 'o'), per charge state z.

    Attributes: name, Z (nuclear charge), mass_amu, W [Z] ionization energy z -> z+1 [eV].
    rates(Te, ne) evaluates, at electron temperature Te [eV] and density ne [m^-3] (arrays of one shape):
        'S'    [Z, ...]  ionization z -> z+1, z = 0..Z-1                     [m^3/s]
        'alpha' [Z, ...] recombination z -> z-1, entry z-1 for z = 1..Z       [m^3/s]
        'cx'   [Z, ...]  charge-exchange recombination with neutral H/D, z -> z-1, entry z-1   [m^3/s]
        'plt'  [Z, ...]  line power of charge state z, z = 0..Z-1             [W m^3]
        'prb'  [Z, ...]  recombination + bremsstrahlung power of z, entry z-1 [W m^3]
    Interpolation is linear in (log T, log ne) and clamped at the table edges. ccd tables are
    evaluated at the temperature given (which temperature belongs there is a model choice)."""

    def __init__(self, name):
        self.name = name.lower()
        files = ADAS_FILES[self.name]
        self.W = np.array(IONIZATION_EV[self.name])
        self.Z, self.mass_amu = len(self.W), MASS_AMU[self.name]
        self.tables, self._interp = {}, {}
        for kind, fname in files.items():
            t = self.tables[kind] = read_adf11(adas_file(fname))
            if t['z_max'] != self.Z or len(t['blocks']) != self.Z:
                raise ValueError(f'{fname}: {len(t["blocks"])} blocks for Z = {self.Z}')
            self._interp[kind] = [RegularGridInterpolator((t['logT'], t['logne']), b) for b in t['blocks']]

    def _eval(self, kind, Te, ne):
        t = self.tables[kind]
        lT = np.clip(np.log10(np.maximum(Te, 1e-30)), t['logT'][0], t['logT'][-1])
        ln = np.clip(np.log10(np.maximum(ne, 1e-30)), t['logne'][0], t['logne'][-1])
        lT, ln = np.broadcast_arrays(lT, ln)
        pts = np.stack([lT.ravel(), ln.ravel()], axis=-1)
        return np.array([10.0 ** f(pts) for f in self._interp[kind]]).reshape((-1,) + lT.shape)

    def rates(self, Te, ne):
        Te, ne = np.asarray(Te, float), np.asarray(ne, float)
        return {'S': self._eval('scd', Te, ne), 'alpha': self._eval('acd', Te, ne), 'cx': self._eval('ccd', Te, ne),
                'plt': self._eval('plt', Te, ne), 'prb': self._eval('prb', Te, ne)}


def fractional_abundance(sp, Te, ne, tau=np.inf):
    """Steady-state charge-state fractions [Z+1] of species sp at Te [eV], ne [m^-3] with ion
    confinement time tau [s], every lost ion returning as a neutral (thesis eq. 5.2; coronal for
    tau = inf). Scalar Te, ne."""
    r = sp.rates(Te, ne)
    S, al = np.ravel(r['S']), np.ravel(r['alpha'])      # [Z] at the one (Te, ne)
    Z = sp.Z
    A = np.zeros((Z + 1, Z + 1))                         # dn/dt / ne = A n
    loss = 0.0 if np.isinf(tau) else 1.0 / (ne * tau)
    for z in range(Z + 1):
        if z < Z:                                        # ionization z -> z+1
            A[z, z] -= S[z]
            A[z + 1, z] += S[z]
        if z > 0:                                        # recombination z -> z-1, transport loss back to neutral
            A[z, z] -= al[z - 1] + loss
            A[z - 1, z] += al[z - 1]
            A[0, z] += loss
    A[-1, :] = 1.0                                       # replace one equation by sum(n) = 1
    rhs = np.zeros(Z + 1)
    rhs[-1] = 1.0
    return np.linalg.solve(A, rhs)
