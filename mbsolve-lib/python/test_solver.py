import sys
sys.path.append('./mbsolve-lib/')
sys.path.append('./solver-generic/')
sys.path.append('./writer-matlab/')

import pymbsolvelib as mb
import pysolvergeneric
import pywritermatlab

import numpy as np

dev = mb.device("test")

sce = mb.scenario("basic", 32768, 200e-15)

sce.add_record(mb.record("d11"))

sol = mb.solver("generic", dev, sce)

print('Solver ' + sol.get_name() + ' started')

sol.run()

print('Solver ' + sol.get_name() + ' finished')

wri = mb.writer("matlab")

outfile = dev.get_name() + "_" + sce.get_name() + "." + wri.get_extension()

results = sol.get_results()

print("I have " + str(len(results)) + " result(s)")

d11 = np.matrix(results[0].get_data_complex()).reshape(results[0].get_rows(),
                                                       results[0].get_cols())

print(d11)



wri.write(outfile, sol.get_results(), dev, sce)
