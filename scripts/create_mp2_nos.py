#!/usr/bin/python3
import argparse
import pyscf as pyscf
from pyscf import cc, gto, scf, mp, mcscf
from pyscf.tools import molden, fcidump
import os
import scipy
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--mol','-m', help =
   'Molecule definition ("atom" argument to pyscf, e.g. "H 0 0 0; F 0 0 1").',
   type=str, default="Be 0 0 0")
parser.add_argument('--basis', '-b', help='basis', type=str, default='sto-3g')
parser.add_argument('--charge', '-c', help='charge', type=int, default=0)
parser.add_argument('--spin', '-s', help='spin multiplicity', type=int, default=0)
parser.add_argument('--coupled-cluster', help='Do a Coupled Cluster calculation', action='store_true')
parser.add_argument('--no-fcidump', help='turn off FCIDUMP printing', action='store_true')
parser.add_argument('--frozen', '-f', help='frozen core for CC', type=int, default=0)
parser.add_argument('--print', '-p', help='print out more info', action='store_true')

args = parser.parse_args()


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
   spin = args.spin, charge = args.charge)

if args.spin > 0:
    mol_rhf = pyscf.scf.HF(mol).run()
else:
    mol_rhf = pyscf.scf.RHF(mol).run()

print("Number of spin-orbitals: ", 2 * len(mol_rhf.mo_coeff))

print("Running MP2")
pt = mp.MP2(mol_rhf, frozen=args.frozen).run()

rdm1 = pt.make_rdm1()

if (args.spin > 0):
    occ, no = scipy.linalg.eigh(rdm1[0])
else:
    occ, no = scipy.linalg.eigh(rdm1)

occ = occ[::-1]
no = no[:,::-1]

np.savetxt('mp-no-occupation', occ)

if args.spin > 0:
    mp2_coeffs = mol_rhf.mo_coeff[0].dot(no)
else:
    mp2_coeffs = mol_rhf.mo_coeff.dot(no)


if not(args.no_fcidump):
    print("Creating reference FCIDUMP file.")
    fcidump.from_mo(mol, 'FCIDUMP.bare', mp2_coeffs)

print("Creating molden file.")
if args.spin > 0:
    molden.from_mo(mol, 'input.molden', mp2_coeffs, occ=(pt.mo_occ[0]+pt.mo_occ[1]))
else:
    molden.from_mo(mol, 'input.molden', mp2_coeffs, occ=pt.mo_occ)

