"""
Script for computing F12 integrals using Libint2
"""

__authors__    = "Daniel G. A. Smith"
__credits__   = ["Daniel G. A. Smith", "Dominic A. Sirianni", "Rob Parrish"]

__copyright__ = "(c) 2014-2018, The Psi4NumPy Developers"
__license__   = "BSD-3-Clause"
__date__      = "2017-05-23"

import time
import numpy as np
import sys
from stggtg_6810 import *
np.set_printoptions(precision=10, linewidth=200, suppress=True)
import psi4
from numpy import linalg as LA

# Memory for Psi4 in GB
psi4.set_memory('2 GB')
psi4.core.set_output_file('output.dat', False)

# Memory for numpy in GB
numpy_memory = 2

mol = psi4.geometry("""
O
H 1 1.1
H 1 1.1 2 104
symmetry c1
""")

psi4.set_options({'basis': 'cc-pvdz',
                  'scf_type': 'pk',
                  'mp2_type': 'conv',
                  'e_convergence': 1e-8,
                  'd_convergence': 1e-8})

# Check energy against psi4?
check_energy = False

print('\nStarting SCF and integral build...')
t = time.time()

# First compute SCF energy using Psi4
scf_e, wfn = psi4.energy('SCF', return_wfn=True)

# Grab data from wavfunction class 
nmo = wfn.nmo()

# Compute size of ERI tensor in GB
ERI_Size = (nmo ** 4) * 8e-9
print('Size of the ERI/MO tensor will be %4.2f GB.' % ERI_Size)
memory_footprint = ERI_Size * 2.5
if memory_footprint > numpy_memory:
    clean()
    raise Exception("Estimated memory utilization (%4.2f GB) exceeds numpy_memory \
                    limit of %4.2f GB." % (memory_footprint, numpy_memory))

print('Building ERI integrals.')
# Integral generation from Psi4's MintsHelper
wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option('BASIS'))
conv = psi4.core.BasisSet.build(mol,'BASIS',psi4.core.get_global_option('BASIS'))
mints = psi4.core.MintsHelper(wfn.basisset())
I = np.asarray(mints.ao_eri(conv,conv,conv,conv))

# Compute F12/STG integrals 
# ---------------------------------------
f12stg_ints = np.asarray(mints.ao_f12_stg(1.0))

# Compute F12/cGTG integrals
# -----------------------------------------------
cgtg_params = stggtg(1.0, cgtg = 6)
f12cgtg_ints = np.asarray(mints.ao_f12(cgtg_params))
print("\nF12/cGTG = 6: ")
print(round(LA.norm(f12stg_ints-f12cgtg_ints),8))

cgtg_params = stggtg(1.0, cgtg = 8)
f12cgtg_ints = np.asarray(mints.ao_f12(cgtg_params))
print("\nF12/cGTG = 8: ")
print(round(LA.norm(f12stg_ints-f12cgtg_ints),8))

cgtg_params = stggtg(1.0, cgtg = 10)
f12cgtg_ints = np.asarray(mints.ao_f12(cgtg_params))
print("\nF12/cGTG = 10: ")
print(round(LA.norm(f12stg_ints-f12cgtg_ints),8))

cgtg_params = stggtg(1.0, cgtg = 6, eventempered=True)
f12cgtg_ints = np.asarray(mints.ao_f12(cgtg_params))
print("\nF12/cGTG = 6, even-tempered: ")
print(round(LA.norm(f12stg_ints-f12cgtg_ints),8))

cgtg_params = stggtg(1.0, cgtg = 8, eventempered=True)
f12cgtg_ints = np.asarray(mints.ao_f12(cgtg_params))
print("\nF12/cGTG = 8, even-tempered: ")
print(round(LA.norm(f12stg_ints-f12cgtg_ints),8))

cgtg_params = stggtg(1.0, cgtg = 10, eventempered=True)
f12cgtg_ints = np.asarray(mints.ao_f12(cgtg_params))
print("\nF12/cGTG = 10, even-tempered: ")
