#!/usr/bin/env python3
from math import factorial, sqrt, comb
from datetime import datetime
import numpy as np
import contextlib
import warnings
import shutil
import copy
import math
import time
import sys
import os
import re

BOHR2ANGST = 0.52917726
OCC_TOL = 1e-4

DD_BANNER = (
    "\n"
    " ****************************************************************\n"
    " ***                                                          ***\n"
    " ***             GIMIC – Density Decomposition                ***\n"
    " ***                                                          ***\n"
    " ****************************************************************\n"
    f" {datetime.now():%a %b %d %H:%M:%S %Y}\n"
    "\n"
)

cart_nbf_per_l = {'s': 1, 'p': 3, 'd': 6, 'f': 10, 'g': 15, 'h': 21}
spher_nbf_per_l = {'s': 1, 'p': 3, 'd': 5, 'f': 7, 'g': 9, 'h': 11}
l_value = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5}
l_label = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g', 5: 'h'}
"""
Order of Cartesian GTO in XDENS file
"""
cartesian_orders = {
    's': [(0, 0, 0)],
    'p': [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
    'd': [(2, 0, 0), (0, 2, 0), (0, 0, 2),
          (1, 1, 0), (1, 0, 1), (0, 1, 1)],
    'f': [(3, 0, 0), (0, 3, 0), (0, 0, 3),
          (1, 2, 0), (2, 1, 0), (2, 0, 1),
          (1, 0, 2), (0, 1, 2), (0, 2, 1), (1, 1, 1)],
    'g': [(4, 0, 0), (0, 4, 0), (0, 0, 4),
          (3, 1, 0), (3, 0, 1), (1, 3, 0),
          (0, 3, 1), (1, 0, 3), (0, 1, 3),
          (2, 2, 0), (2, 0, 2), (0, 2, 2),
          (2, 1, 1), (1, 2, 1), (1, 1, 2)],
    # CHECK h
    'h': [(5, 0, 0), (0, 5, 0), (0, 0, 5),
          (4, 1, 0), (4, 0, 1), (1, 4, 0),
          (0, 4, 1), (1, 0, 4), (0, 1, 4),
          (3, 1, 1), (1, 3, 1), (1, 1, 3),
          (3, 2, 0), (3, 0, 2), (2, 3, 0),
          (0, 3, 2), (2, 0, 3), (0, 2, 3),
          (2, 2, 1), (2, 1, 2), (1, 2, 2)]}

spherical_orders = {'s': [0],
                    'p': [1, -1, 0],
                    'd': [2, 1, -1, -2, 0],
                    'f': [-3, -2, -1, 0, 1, 2, 3],
                    'g': [-4, -3, -2, -1, 0, 1, 2, 3, 4],
                    'h': [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]}


def factorial2(n):
    if n == -1 or n == 0:
        return 1
    if n < -1:
        raise ValueError("not defined for n < -1")

    result = 1
    while n > 1:
        result *= n
        n -= 2
    return result


class ContractedGTO:
    def __init__(self, exponents, coeffs, l_val):
        """
        Contraction coefficients will be multiplied by CGTO normalization.
        Parameters:
            exponents (list of float): List of exponents for the primitives.
            coeffs (list of float): List of contraction coefficients.
            l_val (str): "s", "p", ...
        """
        self.exponents = np.array(exponents, dtype='float')
        self.coeffs = np.array(coeffs, dtype='float')
        self.cart_numbers = cartesian_orders[l_val]
        self.norm_primitives = self.primitive_norm()
        self.norm_cgto = self.cgto_norm()
        self.coeffs *= self.norm_cgto[tuple(self.cart_numbers[0])]

    def primitive_norm(self):
        """
        N_{lmn}(a) = (2a/π)^(3/4) * sqrt((4a)^(l+m+n) / \\Prod_{k in order} ((2k-1)!!))
        """
        norm_list = {}
        prefactor = (2 * self.exponents / math.pi) ** 0.75
        for cn in self.cart_numbers:
            tot_exp = sum(cn)
            power_term = np.power(4 * self.exponents, tot_exp)
            denom = np.prod([factorial2(2 * k - 1) if k > 0 else 1 for k in cn])
            norm_list[tuple(cn)] = prefactor * np.sqrt(power_term / denom)
        return norm_list

    def cgto_norm(self):
        """
        Calculates the normalization coefficient for the contracted GTO for each
        Cartesian combination specified in self.cart_numbers.
        Returns a dictionary where the key is a tuple of angular orders and the value
        is the normalization coefficient.
        """
        norm_dict = {}
        for cn in self.cart_numbers:
            S = 0.0
            N_prim = self.primitive_norm()[tuple(cn)]
            n = len(self.exponents)
            for i in range(n):
                for j in range(n):
                    gamma_ij = self.exponents[i] + self.exponents[j]
                    n_x, n_y, n_z = cn
                    s_x = math.sqrt(math.pi) * (factorial2(2 * n_x - 1) if n_x > 0 else 1) / (
                            2 ** n_x * gamma_ij ** ((2 * n_x + 1) / 2))
                    s_y = math.sqrt(math.pi) * (factorial2(2 * n_y - 1) if n_y > 0 else 1) / (
                            2 ** n_y * gamma_ij ** ((2 * n_y + 1) / 2))
                    s_z = math.sqrt(math.pi) * (factorial2(2 * n_z - 1) if n_z > 0 else 1) / (
                            2 ** n_z * gamma_ij ** ((2 * n_z + 1) / 2))
                    s_tot = s_x * s_y * s_z
                    S += self.coeffs[i] * self.coeffs[j] * N_prim[i] * N_prim[j] * s_tot
            norm_dict[tuple(cn)] = 1.0 / np.sqrt(S)
        return norm_dict


class Shell:
    """
    List of contracted Gaussian functions with the same l for an atom.

    Parameters:
      l (str): The orbital type ('s', 'p', 'd', etc.).
      primitives (list of ContractedGTO): List of primitive GTOs that form the contracted function.
    """

    def __init__(self, l, contracted_primitives):
        self.l = l.lower()
        self.contracted_primitives = contracted_primitives


class Atom:
    def __init__(self, element, number, xyz, basis):
        """
        Parameters:
          element (str): Atomic symbol.
          xyz (list or np.array): Cartesian coordinates.
          number (int): Sequential atom number.
          basis (list of Shell): List of shells associated with the atom.
        """
        self.element = element
        self.position = np.array(xyz, dtype=float)
        self.number = number
        self.basis = basis
        self.cart_nbf = sum([cart_nbf_per_l[s.l] * len(s.contracted_primitives) for s in basis])
        self.spher_nbf = sum([spher_nbf_per_l[s.l] * len(s.contracted_primitives) for s in basis])

    def __repr__(self):
        return f"Atom({self.element}, {self.position}, {self.cart_nbf})"

    def distance(self, other_atom):
        return np.linalg.norm(self.position - other_atom.position)


class Calculator:
    def __init__(self, atoms):
        """
        Parameters:
            atoms (list of Atom): List of atoms in the molecule.
        """
        self.atoms = atoms
        self.basis = self.build_basis()
        self.Ncao = sum([a.cart_nbf for a in atoms])
        self.Nsao = sum([a.spher_nbf for a in atoms])

        # self.max_l = max(l_value[cg[1]] for cg in self.basis)
        # self.cao2sao_blocks = self.calculate_sao2cao_blocks()
        # self.sao2cao_blocks = { m: np.linalg.pinv(v) for m, v in self.cao2sao_blocks.items() }
        # self.c2s = self.cartesian_to_spherical_matrix()
        self.cache = {}

    def sort_shells_by_l(self):
        """ Sort shells first by l, then by atom"""
        shell_per_atom = np.array([cart_nbf_per_l[s.l]
                                   for a in self.atoms
                                   for s in a.basis
                                   for cg in s.contracted_primitives
                                   for _ in range(cart_nbf_per_l[s.l])])
        atom_per_shell = np.array([a.number
                                   for a in self.atoms
                                   for s in a.basis
                                   for cg in s.contracted_primitives
                                   for _ in range(cart_nbf_per_l[s.l])])
        idces1 = np.argsort(atom_per_shell, kind="mergesort")
        shell_per_atom = shell_per_atom[idces1]
        idces2 = np.argsort(shell_per_atom, kind="mergesort")
        final_idx = idces1[idces2]
        return final_idx

    def sort_shells_by_atom(self):
        """ Sort shells first by atoms, then by l"""
        shell_per_atom = np.array([cart_nbf_per_l[s.l] for a in atom_list for s in a.basis])
        atom_per_shell = np.array([a.number for a in atom_list for s in a.basis])
        idces1 = np.argsort(shell_per_atom, kind="mergesort")
        atom_per_shell = atom_per_shell[idces1]
        idces2 = np.argsort(atom_per_shell, kind="mergesort")
        final_idx = idces1[idces2]
        return final_idx

    def build_basis(self):
        """
        Extracts all basis functions from the atoms.
        For each shell, the Cartesian components are expanded.

        Returns:
            list of tuples: (center, l, ContractedGTO, order)
        """
        basis_functions = []
        for atom in self.atoms:
            for shell in atom.basis:
                orders = cartesian_orders[shell.l]
                for gto in shell.contracted_primitives:
                    for order in orders:
                        basis_functions.append((atom.position, shell.l, gto, order))
        return basis_functions

    def calculate_sao2cao_blocks(self):
        def N_S_lm(l, m):
            """
            (9.1.10)
            """
            res = sqrt(2 * factorial(l + abs(m)) * factorial(l - abs(m))) / (factorial(l) * 2 ** abs(m))
            if m == 0:
                res /= sqrt(2)
            return res

        def C_lm_tuv(l, m, t, u, v, vm):
            """
            (9.1.11)
            """
            res = (-1) ** (t + v - vm)
            res *= 0.25 ** t
            res *= comb(l, t) * comb(l - t, abs(m) + t) * comb(t, u) * comb(abs(m), int(2 * v))

            return res

        blocks = {}
        for l_int in range(self.max_l + 1):
            l_str = l_label[l_int]
            n_spher = spher_nbf_per_l[l_str]
            n_cart = cart_nbf_per_l[l_str]
            T_block = np.zeros((n_spher, n_cart))
            for ii, m in enumerate(spherical_orders[l_str]):
                N = N_S_lm(l_int, m)
                vm = 0 if m >= 0 else 0.5
                max_t = math.floor((l_int - abs(m)) / 2)
                max_v = math.floor(abs(m) / 2 - vm) + vm
                for t in range(max_t + 1):
                    for u in range(t + 1):
                        v_range = np.arange(vm, max_v, 1)
                        if v_range.size == 0 or not np.isclose(v_range[-1], max_v):
                            v_range = np.append(v_range, max_v)
                        for v in v_range:
                            i = 2 * t + abs(m) - 2 * (u + v)
                            j = 2 * (u + v)
                            k = l_int - 2 * t - abs(m)
                            j = cartesian_orders[l_str].index((i, j, k))
                            T_block[ii, j] += N * C_lm_tuv(l_int, m, t, u, v, vm)
            if not np.allclose(T_block.T @ T_block, np.eye(n_cart), atol=1e-8):
                print(f"Warning: T_block for l={l_str} is not orthogonal!")
            blocks[l_str] = T_block
        return blocks

    def cartesian_to_spherical_matrix(self):
        """
        Collects the full block-diagonal transformation matrix for all atoms.
        For each atom and its shells the corresponding block from self.sao2cao_blocks is taken.
        """
        cao2sao = np.zeros((self.Nsao, self.Ncao))
        cart_index = 0
        spher_index = 0
        for atom in self.atoms:
            for shell in atom.basis:
                l = shell.l.lower()
                n_cart = len(cartesian_orders[l])
                n_spher = spher_nbf_per_l[l]
                T_block = self.cao2sao_blocks[l]
                for _ in range(len(shell.contracted_primitives)):
                    cao2sao[spher_index:spher_index + n_spher, cart_index:cart_index + n_cart] = T_block
                    cart_index += n_cart
                    spher_index += n_spher
        return cao2sao

    def os_overlap_1d(self, i, j, a, b, A_coord, B_coord):
        """
        Calculates the one-dimensional overlap between primitive Gaussians
        with exponents a and b and angular momentum quantum numbers i and j
        aaccording to the Obara-Saika recurrence relations (9.3.8)-(9.3.9):

          Base case:
            S(0,0) = sqrt(pi/(a+b)) * exp[- (a*b/(a+b)) * (A_coord - B_coord)^2 ]

          For i = 0, j >= 1:
            S(0,j) = (P - B_coord)*S(0,j-1) + ((j-1)/(2*(a+b)))*S(0,j-2)

          For j = 0, i >= 1:
            S(i,0) = (P - A_coord)*S(i-1,0) + ((i-1)/(2*(a+b)))*S(i-2,0)

          For i >= 1, j >= 1:
            S(i,j) = (P - A_coord)*S(i-1,j) + (1/(2*(a+b)))*[ (i-1)*S(i-2,j) + j*S(i-1,j-1) ]

        where:
          p = a + b,
          P = (a*A_coord + b*B_coord)/p.

        The result is cached.
        """
        if i < 0 or j < 0:
            return 0.0
        p = a + b
        P = (a * A_coord + b * B_coord) / p
        R_AB = A_coord - B_coord
        key = (i, j, a, b, A_coord, B_coord)
        if key in self.cache:
            return self.cache[key]

        if i == 0 and j == 0:
            result = math.sqrt(math.pi / p) * math.exp(- (a * b / p) * (R_AB ** 2))
        elif i == 0:
            # For i = 0, j >= 1:
            if j == 1:
                result = (P - B_coord) * self.os_overlap_1d(0, 0, a, b, A_coord, B_coord)
            else:
                result = (P - B_coord) * self.os_overlap_1d(0, j - 1, a, b, A_coord, B_coord) \
                         + ((j - 1) / (2 * p)) * self.os_overlap_1d(0, j - 2, a, b, A_coord, B_coord)
        elif j == 0:
            # For j = 0, i >= 1:
            if i == 1:
                result = (P - A_coord) * self.os_overlap_1d(0, 0, a, b, A_coord, B_coord)
            else:
                result = (P - A_coord) * self.os_overlap_1d(i - 1, 0, a, b, A_coord, B_coord) \
                         + ((i - 1) / (2 * p)) * self.os_overlap_1d(i - 2, 0, a, b, A_coord, B_coord)
        else:
            # For i >= 1 and j >= 1, use the recurrence relation (9.3.8):
            result = (P - A_coord) * self.os_overlap_1d(i - 1, j, a, b, A_coord, B_coord) \
                     + (1 / (2 * p)) * ((i - 1) * self.os_overlap_1d(i - 2, j, a, b, A_coord, B_coord) +
                                        j * self.os_overlap_1d(i - 1, j - 1, a, b, A_coord, B_coord))
        self.cache[key] = result
        return result

    def compute_overlap_matrix(self):
        """
        Computes the complete overlap matrix for the entire molecule.
        Returns:
            np.array: A square numpy overlap matrix of size Ncao x Ncao.
        """
        overlap_matrix = np.zeros((self.Ncao, self.Ncao))
        for i in range(self.Ncao):
            center_i, l_i, gto_i, order_i = self.basis[i]
            norm_coeffs_i = gto_i.coeffs
            prim_norm_i = gto_i.norm_primitives[order_i]
            for j in range(i, self.Ncao):
                center_j, l_j, gto_j, order_j = self.basis[j]
                norm_coeffs_j = gto_j.coeffs
                prim_norm_j = gto_j.norm_primitives[order_j]
                Sij = 0.0
                for a, ca, Npi in zip(gto_i.exponents, norm_coeffs_i, prim_norm_i):
                    for b, cb, Npj in zip(gto_j.exponents, norm_coeffs_j, prim_norm_j):
                        Sx = self.os_overlap_1d(order_i[0], order_j[0], a, b, center_i[0], center_j[0])
                        Sy = self.os_overlap_1d(order_i[1], order_j[1], a, b, center_i[1], center_j[1])
                        Sz = self.os_overlap_1d(order_i[2], order_j[2], a, b, center_i[2], center_j[2])
                        Sij += ca * cb * Npi * Npj * (Sx * Sy * Sz)
                overlap_matrix[i, j] = Sij
                overlap_matrix[j, i] = Sij
        return overlap_matrix


def index_str(str, word, number=0):
    for i in range(number, len(str)):
        if word in str[i]:
            return i


def read_MOL(fname):
    """
    Reads MOL file and returns a list of Atom objects.
    """
    atom_list = []
    Nelec = 0
    atom_number = 1
    n = 0

    try:
        with open(fname) as f:
            lines = f.readlines()[5:]

        while n < len(lines):
            line = lines[n].split()
            if line[-1].isdigit():
                Nelec += float(line[0])
                # Read number of contracted GTO per shell
                nshell = list(map(int, line[3:]))
                n += 1

                atom_label, dublicat_number = lines[n][:2].strip(), int(lines[n][3:5])
                if dublicat_number != 1:
                    print(f"dublicat_number != 1 n: {n}")
                atom_coord = np.array(lines[n].split()[-3:], dtype=float)
                n += 1
                shell_list = []
                for l, o in enumerate(nshell):
                    shell_temp = []
                    for _ in range(o):
                        ngauss = int(lines[n].split()[0])
                        exp_temp, coef_temp = zip(*[lines[n + i].split()[:2] for i in range(1, ngauss + 1)])
                        shell_temp.append(ContractedGTO(np.array(exp_temp), np.array(coef_temp), l_label[l]))
                        n += ngauss + 1
                    shell_list.append(Shell(l_label[l], shell_temp))

                atom_list.append(Atom(atom_label, atom_number, atom_coord, shell_list))
                atom_number += 1
            else:
                n += 1
    except FileNotFoundError:
        sys.stderr.write(f"File with basis info {fname} not found.\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"Failed to read line n={n} from {fname}: {e}\n")
        sys.exit(1)

    return atom_list


class Turbomole:

    def __init__(self, workdir, aomix_fname='aomix.in', xdens_fname='XDENS'):
        self.workdir = workdir
        self.xdens_fname = xdens_fname
        self.aomix_fname = aomix_fname
        if self.aomix_fname is not None:
            self.read_AOMix()
        norm = self.C.T @ self.S @ self.C
        if not (np.isclose(np.trace(norm), self.n_sao, atol=1e-8)):
            warnings.warn(f"Tr(C.T @ S @ C) = {np.trace(norm)} is not fulfilled with accuracy 1e-8 "
                          f"for matrices from {self.aomix_fname}")
        self.read_density()

    def read_AOMix(self):
        mos = []
        irrep_list = []
        energy_list = []
        mo_names = []
        n = 0
        try:
            with open(self.aomix_fname) as f:
                lines = f.readlines()

                n = index_str(lines, "[MO]") + 1
                while all(s not in lines[n] for s in ['[', ']']):
                    irrep_list.append(lines[n].split()[1])
                    energy_list.append(float(lines[n + 1].split()[1]))
                    n += 4
                    temp = []
                    while all(s not in lines[n] for s in ['[', ']', 'Sym=']):
                        temp.append(float(lines[n].split()[5]))
                        if len(mos) == 0:
                            mo_names.append(lines[n].split()[4])
                        n += 1
                    mos.append(temp)
                self.C = np.array(mos).T
                self.n_cao = self.C.shape[0]
                self.n_sao = self.C.shape[1]

                n = index_str(lines, "overlap matrix", n)

                self.S = np.zeros((self.n_cao, self.n_cao))
                for i in range(self.n_cao):
                    for j in range(i + 1):
                        n += 1
                        self.S[i, j] = float(lines[n].split()[6])
                        self.S[j, i] = float(lines[n].split()[6])

            self.irrep_list = irrep_list
            self.energy_list = energy_list
            R = np.eye(self.n_cao)

            # R - renormalization matrix to connect C and S with densities from XDENS
            for i, m in enumerate(mo_names):
                if m in ["dx2", "dy2", "dz2"]:
                    R[i, i] *= 2 / np.sqrt(12)

                elif m in ["fx3", "fy3", "fz3"]:
                    R[i, i] *= 2 / np.sqrt(60)

                elif m == "fy2x":
                    R[i, i] = 0.0
                    R[i + 2, i] = 2 / np.sqrt(12)

                elif m == "fx2y":
                    R[i, i] = 0.0
                    R[i - 1, i] = 2 / np.sqrt(12)

                elif m == "fx2z":
                    R[i, i] = 0.0
                    R[i - 1, i] = 2 / np.sqrt(12)

                elif m == "fz2x":
                    R[i, i] = 0.0
                    R[i + 1, i] = 2 / np.sqrt(12)

                elif m == "fz2y":
                    R[i, i] = 0.0
                    R[i + 1, i] = 2 / np.sqrt(12)

                elif m == "fy2z":
                    R[i, i] = 0.0
                    R[i - 2, i] = 2 / np.sqrt(12)

                elif m in ["gx4", "gy4", "gz4"]:
                    R[i, i] *= 8 / np.sqrt(6720)
                elif m in ["gx2yz", "gy2xz", "gz2xy"]:
                    R[i, i] *= 2 / np.sqrt(12)
                elif m in ["gx3y", "gy3x", "gx3z", "gy3z", "gz3x", "gz3y"]:
                    R[i, i] *= 1 / np.sqrt(15)
                elif m in ["gx2y2", "gx2z2", "gy2z2"]:
                    R[i, i] *= 1 / 3

            self.C = R.dot(self.C)
            self.S = np.linalg.inv(R).T.dot(self.S).dot(np.linalg.inv(R))

        except FileNotFoundError:
            sys.stderr.write(f"Missing file with MO CAO coefficients: {self.aomix_fname}.\n")
            sys.exit(1)
        except Exception as e:
            sys.stderr.write(f"Failed to read line n={n} from {self.aomix_fname}: {e}.\n")
            sys.exit(1)


    def read_density(self):
        """
        Get the density matrix
        TM 7.6 ['CAODENS', 'XCAODENS"]
        TM 7.8 ['XDENS']
        """
        try:
            D = np.loadtxt(self.xdens_fname)
            if self.n_cao is None:
                self.n_cao = int(np.sqrt(len(D) / 4))

            self.densities = D.reshape(4, self.n_cao, self.n_cao)
            self.densities = np.array([D.T for D in self.densities])

        except FileNotFoundError:
            sys.stderr.write(f"File XDENS {self.xdens_fname} is not found.\n")
            sys.exit(1)
        except Exception as e:
            sys.stderr.write(f"Failed to read {self.xdens_fname}: {e}.\n")
            sys.exit(1)


def decompose_and_writeXDENSes(workdir, mo_groups_dict, n_occ, n_sao, dD_MO, tmat, C,
                              dirname='XDENSes', basis_fname='MOL'):
    """
    Decompose density matrices and write XDENS files per group.

    workdir: base output directory
    mo_groups_dict: {group_name: list of orbital indices order from 0}
    n_occ: number of occupied orbitals
    n_sao: number of basis functions
    dD_MO: list of 3 perturbation density matrices (Dx, Dy, Dz) in MO basis
    tmat: reordering matrix
    C: MO coefficients
    """
    base_dir = os.path.join(workdir, dirname)
    os.makedirs(base_dir, exist_ok=True)

    for group_name, mo_indices in mo_groups_dict.items():
        print(f" Decomposing and writing XDENS for group: {group_name}", flush=True)
        start_time = time.perf_counter()

        # Build density matrices in MO
        D = np.zeros((4, n_sao, n_sao))
        for orb in mo_indices:
            D[0, orb, orb] = 2.0
            for d in range(1, 4):
                DD = np.zeros((n_sao, n_sao))
                DD[orb, :] = dD_MO[d - 1][orb, :]
                DD[:, orb] = dD_MO[d - 1][:, orb]
                DD[:n_occ, :n_occ] /= 2
                D[d] += DD

        # Transform to AO and orthogonal AO basis
        D_AO = [C @ mat @ C.T for mat in D]
        D_AO_ortho = [tmat @ mat @ tmat.T for mat in D_AO]

        # Prepare output directory
        group_dir = os.path.join(base_dir, group_name)
        os.makedirs(group_dir, exist_ok=True)

        # Copy MOL file
        try:
            shutil.copy(os.path.join(workdir, basis_fname), os.path.join(group_dir, "MOL"))
        except Exception as e:
            sys.stderr.write(f"[ERROR] Failed to copy MOL to {group_dir}: {e}\n")
            sys.exit(1)

        # Write XDENS file
        try:
            with open(os.path.join(group_dir, "XDENS"), 'w') as f:
                for mat in D_AO_ortho:
                    for e in mat.T.flatten():
                        f.write(f"{e:0.14E}\n")
                    f.write('\n')
        except Exception as e:
            sys.stderr.write(f"[ERROR] Failed to write XDENS for group '{group_name}': {e}\n")
            sys.exit(1)

        # Log time
        elapsed = time.perf_counter() - start_time
        h, rem = divmod(elapsed, 3600)
        m, s = divmod(rem, 60)
        print(f" Finished {group_name} in {int(h):02}:{int(m):02}:{int(s):02}", flush=True)


def parse_groups_from_inkeys(dd_sect, tm):
    """
    Parse dd_sect.getkw('groups') and return dict of groups
    """
    raw = dd_sect.getkw('groups')
    groups = {}
    rng = re.compile(r'(\d+)\s*-\s*(\d+)')
    for line in raw:
        try:
            line = line.split("#")[0]
            if not line.strip():
                continue

            if line == 'all':
                for i in tm.occ_mo_list:
                    name = f"{i + 1}_" + tm.irrep_list[i]
                    if tm.irrep_list[i] in groups.keys():
                        name += '_1'
                    groups[name] = np.asarray([i], int)
            else:
                key, rhs = line.split('=', 1)
                name = key.strip()
                rhs = rhs.strip().strip('[]')
                orbs = []
                for tok in rhs.split(','):
                    tok = tok.strip()
                    m = rng.fullmatch(tok)
                    if m:
                        a, b = map(int, m.groups())
                        orbs.extend(range(a, b + 1))
                    else:
                        orbs.append(int(tok))
                # Make count of MO from 0
                groups[name] = np.asarray(orbs, int) - 1
        except:
            sys.stderr.write(f"Failed to read from input line: {line}.\n")
            sys.exit(1)
    return groups


def resolve_paths(dd_sect, args, inkeys):
    """Return absolute paths for XDENS, AOMIX, MOL"""
    workdir = os.path.dirname(os.path.abspath(args.infile))
    basis_fname = inkeys.getkw("basis")[0]
    aomix = os.path.join(workdir, dd_sect.getkw('mofile')[0])
    xdens = os.path.join(workdir, inkeys.getkw("xdens")[0])
    mol = os.path.join(workdir, basis_fname)
    return workdir, aomix, xdens, mol, basis_fname


def sanity(tm, groups):
    print(" Sanity check.", flush=True)

    # Check normalization
    norm = tm.C.T @ tm.S @ tm.C
    print(f" Tr(C.T S C) = {np.trace(norm): 7.5e}, n_sao = {tm.n_sao}.", flush=True)
    if not np.isclose(np.trace(norm), tm.n_sao, atol=1e-8):
        print("       Tr(C.T S C) differs from n_sao by >1e-8 ")

    # Check if virtual orbitals are occupied by user input
    fail_flag = False
    for gdir, v in groups.items():
        if not all(i in tm.occ_mo_list for i in v):
            sys.stderr.write(f"Group '{gdir}' contains not occupied orbitals: {v + 1}.\n")
            fail_flag = True
    if fail_flag:
        sys.exit(1)

    # Check if smearing is used in calculation
    correct_occupation = np.zeros(tm.n_sao)
    for i in tm.occ_mo_list:
        correct_occupation[i] = 2.0
    error = np.diag(tm.densities[0]) - correct_occupation
    bad_indices = np.where(np.abs(error) > OCC_TOL)[0]

    if bad_indices.size > 0:
        print("[WARN] Fractional MO occupations detected (possible Fermi smearing):", flush=True)
        for i in bad_indices:
            print(f"MO {i + 1}: occ={tm.densities[0][i, i]:.6f}, error = {error[i]:.6f}", flush=True)


def cdens_calculation(workdir, groups, args, inkeys):
    """Optionally run cdens calculation for each <group>."""

    base = os.path.join(workdir, "XDENSes")
    cwd = os.getcwd()  # Save current working directory

    for g in groups:
        print(f" Running cdens in: {g}")
        gdir = os.path.join(base, g)
        os.chdir(gdir)  # Change to group directory

        # Make a local copy of inkeys and override filenames
        local_log_path = os.path.join(gdir, "gimic.out")
        local_inkeys = copy.deepcopy(inkeys)
        local_inkeys.setkw('xdens', 'XDENS')
        local_inkeys.setkw('basis', 'MOL')

        if inkeys.getkw('backend')[0] == 'fgimic':
            from fgimic.gimic import GimicDriver
        else:
            from pygimic.pygimic import GimicDriver

        driver = GimicDriver(args, local_inkeys)

        try:
            with open(local_log_path, "w") as log_fh, \
                 contextlib.redirect_stdout(log_fh), \
                 contextlib.redirect_stderr(log_fh):
                # Redirect stdout and stderr to log file
                log_fd = log_fh.fileno()
                save_out, save_err = os.dup(1), os.dup(2)
                os.dup2(log_fd, 1)
                os.dup2(log_fd, 2)
                try:
                    driver.run()
                finally:
                    # Restore original stdout and stderr
                    os.dup2(save_out, 1)
                    os.dup2(save_err, 2)
                    os.close(save_out)
                    os.close(save_err)

        except:
            sys.stderr.write(f"Failed to run gimic in {gdir}")
        finally:
            os.chdir(cwd)  # Always return to the original directory

def parse_index_list(lst):
    """
    Converts a list like [1, "2-4", 7] into [1, 2, 3, 4, 7]
    Assumes 1-based indexing.
    """
    result = []
    for item in lst:
        if isinstance(item, int):
            result.append(item)
        elif isinstance(item, str) and "-" in item:
            start, end = map(int, item.split("-"))
            result.extend(range(start, end + 1))
        else:
            raise ValueError(f"Invalid group entry: {item}")
    return np.array(result)

# main entry for gimic
def run(dd_sect, args, inkeys):
    """
    Decompose the density into user-defined MO groups and, if requested,
    run a cdens calculation for every group.
    """
    print(DD_BANNER, end="")

    # 1. Parse user input
    workdir, aomix, xdens, mol, basis_fname = resolve_paths(dd_sect, args, inkeys)
    # print(f" Workdir: {workdir}", flush=True)
    # print(f" Input files:\n  AOMIX: {aomix}\n  XDENS: {xdens}\n  MOL: {mol}", flush=True)

    do_cdens = (dd_sect.getkw('cdens_calc')[0].lower() in ('1', 'yes', 'true', 'on')
                and inkeys.getkw('calc')[0] == 'cdens')
    if do_cdens:
        print(f" cdens calculations enabled for each group.", flush=True)
    else:
        print(f" cdens calculations for each group disabled.", flush=True)

    # 2. Read Turbomole files
    tm = Turbomole(workdir, aomix, xdens)
    atom_list = read_MOL(mol)
    mol_integrator = Calculator(atom_list)
    if getattr(tm, 'S', None) is None:
        tm.S = mol_integrator.compute_overlap_matrix()

    # 3. Form matrix to make reorder
    l_order = mol_integrator.sort_shells_by_l()
    tmat = np.eye(tm.n_cao)[l_order, :]

    tm.densities = np.array([tm.C.T @ tm.S @ (tmat.T @ D @ tmat) @ tm.S @ tm.C
                             for D in tm.densities])
    tm.occ_mo_list = np.where(np.abs(np.diag(tm.densities[0]) - 2.0) < OCC_TOL)[0]
    tm.nocc = len(tm.occ_mo_list)

    # 4. Get groups of orbitals
    groups = parse_groups_from_inkeys(dd_sect, tm)
    if not groups:
        sys.stderr.write("DensityDecomposition: no groups specified.")
        return
    print(" Groups:", flush=True)
    for g, orbs in groups.items():
        print(f"  {g:<15} = {orbs}", flush=True)

    # 5. Sanity
    sanity(tm, groups)

    # 6. Prepare and write XDENSes
    decompose_and_writeXDENSes(workdir, groups, tm.nocc, tm.n_sao, tm.densities[1:], tmat, tm.C,
                               basis_fname=basis_fname)
    # 7. optional cdens runs
    if do_cdens:
        # Run the user file gimic.inp in each group directory, changing only 'basis' and 'xdens'.
        cdens_calculation(workdir, groups, args, inkeys)
    print("Density decomposition finished.", flush=True)

if __name__ == "__main__":
    import yaml

    workdir = os.getcwd()
    aomix_fname = "aomix.in"
    if len(sys.argv) > 1:
        aomix_fname = sys.argv[1]

    # 1. read Turbomole data
    tm = Turbomole(workdir, aomix_fname)
    atoms = read_MOL('MOL')
    integ = Calculator(atoms)
    if getattr(tm, "S", None) is None:
        tm.S = integ.compute_overlap_matrix()

    # 2. Build transformation
    l_order = integ.sort_shells_by_l()
    tmat = np.eye(tm.n_cao)[l_order, :]
    tm.densities = np.array([tm.C.T @ tm.S @ (tmat.T @ D @ tmat) @ tm.S @ tm.C
                             for D in tm.densities])
    tm.occ_mo_list = np.where(np.abs(np.diag(tm.densities[0]) - 2.0) < OCC_TOL)[0]
    tm.nocc = len(tm.occ_mo_list)

    # 3. read groups.yaml
    grp_path = os.path.join(workdir, "groups.yaml")
    if os.path.isfile(grp_path):
        try:
            with open("groups.yaml") as f:
                groups = yaml.safe_load(f)

                if isinstance(groups, dict):
                    groups = {
                        k: parse_index_list(v) - 1 for k, v in groups.items()
                    }
        except:
            print("Failed to read groups.yaml.")
            sys.exit(1)
    else:
        print("groups.yaml not found, using all occupied orbitals.")
        groups = {f"{i+1}_{tm.irrep_list[i]}": [i+1] for i in tm.occ_mo_list}

    # 4. sanity-check
    sanity(tm, groups)

    # 6. Prepare and write XDENSes
    decompose_and_writeXDENSes(workdir, groups, tm.nocc, tm.n_sao, tm.densities[1:], tmat, tm.C,
                               basis_fname=os.path.join(workdir, 'MOL'))