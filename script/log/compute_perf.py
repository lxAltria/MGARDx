import numpy as np
import sys

dims = {}
dims['hurricane'] = 100 * 500 * 500
dims['nyx'] = 512 * 512 * 512
dims['scale'] = 98 * 1200 * 1200
dims['qmcpack'] = 33120 * 69 * 69

time = np.loadtxt('tmp')
size = dims[sys.argv[1]] * 4.0 / 1024 / 1024
perf = size / np.mean(time)
print("Performance = {} MB/s".format(perf))
