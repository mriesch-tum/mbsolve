import pymbsolvelib as mb


vacuum = mb.material('vacuum')

reg = mb.region('Vacuum right', vacuum, 0, 7.5e-6)

dev = mb.Device('QCL')

dev.add_region(reg)
