#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" song2005.py: Runs song2005 three-level V-type setup."""

# mbsolve: An open-source solver tool for the Maxwell-Bloch equations.
#
# Copyright (c) 2016, Computational Photonics Group, Technical University
# of Munich.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

import mbsolve.lib as mb
import mbsolve.solvercpu
import mbsolve.writerhdf5

import math
import time

# Hamiltonian
energies = [ 0, 2.3717e15 * mb.HBAR, 2.4165e15 * mb.HBAR ]
H = mb.qm_operator(energies)

# dipole moment operator
dipoles = [ 9.2374e-11 * mb.E0, 9.2374e-11 * math.sqrt(2) * mb.E0, 0]
u = mb.qm_operator([ 0, 0, 0 ], dipoles)

# relaxation superoperator
rate = 1e10
rates = [ [ 0, rate, rate ], [ rate, 0, rate ], [ rate, rate, 0 ] ]
relax_sop = mb.qm_lindblad_relaxation(rates)

# initial density matrix
rho_init = mb.qm_operator([ 1, 0, 0 ])

# quantum mechanical description
qm = mb.qm_description(6e24, H, u, relax_sop)
mat_ar = mb.material("AR_Song", qm)
mb.material.add_to_library(mat_ar)

# Song setup
dev = mb.device("Song")
dev.add_region(mb.region("Active region (single point)", mat_ar, 0, 0))

# scenario
ic_d = mb.ic_density_const(rho_init)
ic_e = mb.ic_field_const(0.0)
ic_m = mb.ic_field_const(0.0)
sce = mb.scenario("Basic", 1, 80e-15, ic_d, ic_e, ic_m, 10000)
sce.add_record(mb.record("e", 0.0, 0.0))
sce.add_record(mb.record("d11", mb.record.density, 1, 1, 0.0, 0.0))
sce.add_record(mb.record("d22", mb.record.density, 2, 2, 0.0, 0.0))
sce.add_record(mb.record("d33", mb.record.density, 3, 3, 0.0, 0.0))

# add source
sce.add_source(mb.sech_pulse("sech", 0.0, mb.source.hard_source, 3.5471e9,
                             3.8118e14, 17.248, 1.76/5e-15, -math.pi/2))

# run solver
sol = mb.solver.create_instance("cpu-fdtd-3lvl-cvr-rodr", dev, sce)
print('Solver ' + sol.get_name() + ' started')
tic = time.time()
sol.run()
toc = time.time()
print('Solver ' + sol.get_name() + ' finished in ' + str(toc - tic) + ' sec')

# write results
wri = mb.writer.create_instance("hdf5")
outfile = dev.get_name() + "_" + sce.get_name() + "." + wri.get_extension()
results = sol.get_results()
wri.write(outfile, sol.get_results(), dev, sce)
