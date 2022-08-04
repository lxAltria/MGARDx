from netCDF4 import Dataset
import sys
import numpy as np

r1 = int(sys.argv[2])
r2 = int(sys.argv[3])
r3 = int(sys.argv[4])
dataset = Dataset(sys.argv[1], 'w', format='NETCDF4_CLASSIC')
dataset.createDimension('x', r1)
dataset.createDimension('y', r2)
dataset.createDimension('z', r3)
dataset.createVariable('velocity_x', np.float32, ('x', 'y', 'z'))
velocity_x = np.fromfile(sys.argv[5], dtype=np.float32).reshape([r1, r2, r3])
print(dataset.variables['velocity_x'].shape)
dataset.variables['velocity_x'][:] = velocity_x
print(dataset.variables['velocity_x'].shape)
print(dataset.variables['velocity_x'].dtype)
print(np.max(dataset.variables['velocity_x']))
dataset.close()

dataset2 = Dataset(sys.argv[1], 'r', format='NETCDF4_CLASSIC')
print(np.max(dataset2.variables['velocity_x']))
dataset2.close()