"""
Tests for F12DOUBLE_COMMUTATOR integrals - libint2
"""
import psi4
import pytest
from data_molecules.molecules import *
from stggtg import *
from data_libint1 import *

def test_f12double_commutator_he_sto3g():
    """He STO-3G"""
    # Psi4 Setup
    psi4.set_memory('1 GiB')
    psi4.core.set_output_file('output.dat', False)
    psi4.set_options({'basis': 'STO-3G',
                      'scf_type': 'pk',
                      'mp2_type': 'conv',
                      'e_convergence': 1e-8,
                      'd_convergence': 1e-8})
    mol = psi4.geometry(moldict["He"])
    rhf_e, rhf_wfn = psi4.energy('SCF', return_wfn=True)

    # Integral generation from Psi4's MintsHelper
    wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option('BASIS'))
    conv = psi4.core.BasisSet.build(mol,'BASIS',psi4.core.get_global_option('BASIS'))
    mints = psi4.core.MintsHelper(wfn.basisset())

    cgtg_params = stggtg(1.0)
    f12double_commutator_ints = np.asarray(mints.ao_f12_double_commutator(cgtg_params))
    f12double_commutator_ints2 = np.asarray(mints.ao_f12_double_commutator(cgtg_params,conv,conv,conv,conv))

    assert np.allclose(f12double_commutator_ints,f12double_commutator_ints2)

    # Test aginst libint1       
    f12double_commutator_ints_libint1 = np.load("data_libint1/he_sto3g_f12double_commutator.npy")

    assert np.allclose(f12double_commutator_ints,f12double_commutator_ints_libint1)

def test_f12double_commutator_he_cc_pvdz():
    """He cc-pvdz"""
    # Psi4 Setup
    psi4.set_memory('1 GiB')
    psi4.core.set_output_file('output.dat', False)
    psi4.set_options({'basis': 'cc-pvdz',
                      'scf_type': 'pk',
                      'mp2_type': 'conv',
                      'e_convergence': 1e-8,
                      'd_convergence': 1e-8})
    mol = psi4.geometry(moldict["He"])
    rhf_e, rhf_wfn = psi4.energy('SCF', return_wfn=True)

    # Integral generation from Psi4's MintsHelper
    wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option('BASIS'))
    conv = psi4.core.BasisSet.build(mol,'BASIS',psi4.core.get_global_option('BASIS'))
    mints = psi4.core.MintsHelper(wfn.basisset())

    cgtg_params = stggtg(1.0)
    f12double_commutator_ints = np.asarray(mints.ao_f12_double_commutator(cgtg_params))
    f12double_commutator_ints2 = np.asarray(mints.ao_f12_double_commutator(cgtg_params,conv,conv,conv,conv))

    assert np.allclose(f12double_commutator_ints,f12double_commutator_ints2)

    # Test aginst libint1       
    f12double_commutator_ints_libint1 = np.load("data_libint1/he_cc-pvdz_f12double_commutator.npy")

    assert np.allclose(f12double_commutator_ints,f12double_commutator_ints_libint1)

def test_f12double_commutator_he_aug_cc_pvdz():
    """He cc-pvdz"""
    # Psi4 Setup
    psi4.set_memory('1 GiB')
    psi4.core.set_output_file('output.dat', False)
    psi4.set_options({'basis': 'aug-cc-pvdz',
                      'scf_type': 'pk',
                      'mp2_type': 'conv',
                      'e_convergence': 1e-8,
                      'd_convergence': 1e-8})
    mol = psi4.geometry(moldict["He"])
    rhf_e, rhf_wfn = psi4.energy('SCF', return_wfn=True)

    # Integral generation from Psi4's MintsHelper
    wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option('BASIS'))
    conv = psi4.core.BasisSet.build(mol,'BASIS',psi4.core.get_global_option('BASIS'))
    mints = psi4.core.MintsHelper(wfn.basisset())

    cgtg_params = stggtg(1.0)
    f12double_commutator_ints = np.asarray(mints.ao_f12_double_commutator(cgtg_params))
    f12double_commutator_ints2 = np.asarray(mints.ao_f12_double_commutator(cgtg_params,conv,conv,conv,conv))

    assert np.allclose(f12double_commutator_ints,f12double_commutator_ints2)

    # Test aginst libint1       
    f12double_commutator_ints_libint1 = np.load("data_libint1/he_aug-cc-pvdz_f12double_commutator.npy")

    assert np.allclose(f12double_commutator_ints,f12double_commutator_ints_libint1)

def test_f12double_commutator_h2o_sto3g():
    """H2O STO-3G"""
    # Psi4 Setup
    psi4.set_memory('1 GiB')
    psi4.core.set_output_file('output.dat', False)
    psi4.set_options({'basis': 'STO-3G',
                      'scf_type': 'pk',
		      'mp2_type': 'conv',
                      'e_convergence': 1e-8,
                      'd_convergence': 1e-8})    
    mol = psi4.geometry(moldict["H2O"])
    rhf_e, rhf_wfn = psi4.energy('SCF', return_wfn=True)

    # Integral generation from Psi4's MintsHelper
    wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option('BASIS'))
    conv = psi4.core.BasisSet.build(mol,'BASIS',psi4.core.get_global_option('BASIS'))
    mints = psi4.core.MintsHelper(wfn.basisset())

    cgtg_params = stggtg(1.0)
    f12double_commutator_ints = np.asarray(mints.ao_f12_double_commutator(cgtg_params))
    f12double_commutator_ints2 = np.asarray(mints.ao_f12_double_commutator(cgtg_params,conv,conv,conv,conv))

    assert np.allclose(f12double_commutator_ints,f12double_commutator_ints2)	

    # Test aginst libint1	
    f12double_commutator_ints_libint1 = np.load("data_libint1/h2o_sto3g_f12double_commutator.npy")

    assert np.allclose(f12double_commutator_ints,f12double_commutator_ints_libint1)	

def test_f12double_commutator_h2o_cc_pvdz():
    """H2O cc-pVDZ"""
    # Psi4 Setup
    psi4.set_memory('1 GiB')
    psi4.core.set_output_file('output.dat', False)
    psi4.set_options({'basis': 'cc-pVDZ',
                      'scf_type': 'pk',
                      'mp2_type': 'conv',
                      'e_convergence': 1e-8,
                      'd_convergence': 1e-8})
    mol = psi4.geometry(moldict["H2O"])
    rhf_e, rhf_wfn = psi4.energy('SCF', return_wfn=True)

    # Integral generation from Psi4's MintsHelper
    wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option('BASIS'))
    conv = psi4.core.BasisSet.build(mol,'BASIS',psi4.core.get_global_option('BASIS'))
    mints = psi4.core.MintsHelper(wfn.basisset())

    cgtg_params = stggtg(1.0)
    f12double_commutator_ints = np.asarray(mints.ao_f12_double_commutator(cgtg_params))
    f12double_commutator_ints2 = np.asarray(mints.ao_f12_double_commutator(cgtg_params,conv,conv,conv,conv))

    assert np.allclose(f12double_commutator_ints,f12double_commutator_ints2)

    # Test aginst libint1       
    f12double_commutator_ints_libint1 = np.load("data_libint1/h2o_cc-pvdz_f12double_commutator.npy")

    assert np.allclose(f12double_commutator_ints,f12double_commutator_ints_libint1)

def test_f12double_commutator_h2o_aug_cc_pvdz():
    """H2O aug-cc-pVDZ"""
    # Psi4 Setup
    psi4.set_memory('1 GiB')
    psi4.core.set_output_file('output.dat', False)
    psi4.set_options({'basis': 'aug-cc-pVDZ',
                      'scf_type': 'pk',
                      'mp2_type': 'conv',
                      'e_convergence': 1e-8,
                      'd_convergence': 1e-8})
    mol = psi4.geometry(moldict["H2O"])
    rhf_e, rhf_wfn = psi4.energy('SCF', return_wfn=True)

    # Integral generation from Psi4's MintsHelper
    wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option('BASIS'))
    conv = psi4.core.BasisSet.build(mol,'BASIS',psi4.core.get_global_option('BASIS'))
    mints = psi4.core.MintsHelper(wfn.basisset())

    cgtg_params = stggtg(1.0)
    f12double_commutator_ints = np.asarray(mints.ao_f12_double_commutator(cgtg_params))
    f12double_commutator_ints2 = np.asarray(mints.ao_f12_double_commutator(cgtg_params,conv,conv,conv,conv))

    assert np.allclose(f12double_commutator_ints,f12double_commutator_ints2)

    # Test aginst libint1       
    f12double_commutator_ints_libint1 = np.load("data_libint1/h2o_aug-cc-pvdz_f12double_commutator.npy")

    assert np.allclose(f12double_commutator_ints,f12double_commutator_ints_libint1)
	
