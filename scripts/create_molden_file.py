#!/usr/bin/python3
import argparse
import pyscf as pyscf
from pyscf import cc, gto, scf, mp, mcscf, fci
from pyscf.tools import molden, fcidump
import os
import numpy as np

# Build parser.
parser = argparse.ArgumentParser()
parser.add_argument('--mol','-m', help =
   'Molecule definition ("atom" argument to pyscf, e.g. "H 0 0 0; F 0 0 1").',
   type = str, default="Be 0 0 0")
parser.add_argument('--basis', '-b', help = 'basis', type = str,
   default = 'sto-3g')
parser.add_argument('--charge', '-c', help = 'charge', type = int, default = 0)
parser.add_argument('--spin', '-s', help = 'spin multiplicity', type = int,
   default = 0)
parser.add_argument('--coupled-cluster', help='Do a Coupled Cluster calculation',
                    action='store_true')
parser.add_argument('--no-fcidump', help='turn off FCIDUMP printing',
                    action='store_true')
parser.add_argument('--print', '-p', help='print out more info', action='store_true')
parser.add_argument('--frozen', '-f', help='frozen core for CC', type=int, default=0)
parser.add_argument('--ecp', help='Effective core potential', type=str, default='')
parser.add_argument('--fci', help='Do FCI calculation', action='store_true')

# Parse arguments and report.
args = parser.parse_args()
print('System:')
print('* Molecule: "' + args.mol + '"')
print('* Basis: "' + args.basis + '"')
print('* Charge: ', args.charge)
print('* Spin: ', args.spin)

# Build molecule and run RHF.
print("Running SCF.")
mol = pyscf.gto.M(atom = args.mol, basis = args.basis,
   spin = args.spin, charge = args.charge, ecp=args.ecp)
mol_rhf = pyscf.scf.RHF(mol).run()

print("Number of spatial-orbitals: ", len(mol_rhf.mo_coeff))
print("Number of spin-orbitals: ", 2 * len(mol_rhf.mo_coeff))
print("Creating molden file.")
molden.from_scf(mol_rhf, 'input.molden')
if not(args.no_fcidump):
    print("Creating reference FCIDUMP file.")
    fcidump.from_scf(mol_rhf, 'FCIDUMP')

