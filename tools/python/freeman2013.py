import sys
sys.path.append('./mbsolve-lib/')
sys.path.append('./solver-openmp/')
sys.path.append('./writer-matlab/')

import pymbsolvelib as mb
import pysolveropenmp
import pywritermatlab

import numpy as np
import math
import time

# absorber material
# rel permittivity 1
# overlap factor 1
# losses 3e4/cm -> 3e6/m
mat_abs = mb.material("Absorber", None, 1, 1, 3e6)
mb.material.add_to_library(mat_abs)

# Freeman setup varies cavity length (1.5mm/3mm)
if True:
    # cavity length 1.5mm
    L = 1.5e-3
    # total loss 11/cm
    loss = 1100.0
else:
    # cavity length 3mm
    L = 3.0e-3
    # total loss 7/cm
    loss = 700


print(L)
print(loss)
# Freeman 2013 active region
# varies gain recovery time T1
T1 = 20e-12
qm = mb.qm_desc_2lvl(3.7e20, 2 * math.pi * 2.45e12, 6.2e-9, 1/T1, 1/2.35e-12, 1.0)
# background rel permittivity 12.9
# overlap factor 1
mat_ar = mb.material("AR_Freeman", qm, 12.9, 1, loss)
mb.material.add_to_library(mat_ar)

L_abs = 0.25e-3

dev = mb.device("Freeman")
dev.add_region(mb.region("Vacuum left", mat_abs, 0, L_abs))
dev.add_region(mb.region("Active region", mat_ar, L_abs, L_abs + L))
dev.add_region(mb.region("Vacuum right", mat_abs, L_abs + L, 2 * L_abs + L))

# scenario
# approx 14000 grid points -> choose 16k?
# rather set d_x directly
# courant number?
sce = mb.scenario("Basic", 16384, 230e-12)
sce.add_record(mb.record("inv12", 2e-12))
sce.add_record(mb.record("e", 1e-13, L_abs + L))

#TODO check whether position > device length

# add source
sce.add_source(mb.single_cycle_pulse("seed-pulse", L_abs,
                                     mb.source.soft_source, 1e5,
                                     2.45e12, 1/(0.35e-12/2.634), 1e-12))
# TODO: add phase to carrier sine wave ?

# initialization: perfect inversion
sce.set_dm_init_type(mb.scenario.upper_full)

sol = mb.solver("openmp-2lvl-pc", dev, sce)

print('Solver ' + sol.get_name() + ' started')

tic = time.time()
sol.run()
toc = time.time()

print('Solver ' + sol.get_name() + ' finished in ' + str(toc - tic) + ' sec')

wri = mb.writer("matlab")

outfile = dev.get_name() + "_" + sce.get_name() + "." + wri.get_extension()

results = sol.get_results()

#print("I have " + str(len(results)) + " result(s)")

#d11 = np.matrix(results[0].get_data_complex()).reshape(results[0].get_rows(),
#                                                       results[0].get_cols())

#print(d11)



wri.write(outfile, sol.get_results(), dev, sce)
