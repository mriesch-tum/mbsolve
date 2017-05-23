import sys
sys.path.append('./mbsolve-lib/')
sys.path.append('./solver-generic/')
import pymbsolvelib as mb
#import pysolvergeneric as gen

dev = mb.device("test")

sce = mb.scenario("basic", 32768, 200e-15)

sol = mb.solver("generic", dev, sce)

print('Solver ' + sol.get_name() + ' started')

sol.run()

print('Solver ' + sol.get_name() + ' finished')

#except ValueError:
#    print('Error: ')
#    raise
