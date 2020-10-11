import numpy as np
import sys

dims = {}
dims['hurricane'] = 100 * 500 * 500
dims['nyx'] = 512 * 512 * 512
dims['scale'] = 98 * 1200 * 1200
dims['qmcpack'] = 33120 * 69 * 69

perf = np.zeros([4, 4])
time = np.loadtxt('tmp')
start = 0
end = 13 * 4
time_field = time[start:end].reshape([13, 4])
perf[0, :] = dims['hurricane'] * 4.0 / 1024 / 1024 / np.mean(time_field, axis = 0)
print(perf[0, :])
start = end
end = end + 6 * 4
time_field = time[start:end].reshape([6, 4])
print(time_field)
perf[1, :] = dims['nyx'] * 4.0 / 1024 / 1024 / np.mean(time_field, axis = 0)
print(perf[1, :])
start = end
end = end + 12 * 4
time_field = time[start:end].reshape([12, 4])
perf[2, :] = dims['scale'] * 4.0 / 1024 / 1024 / np.mean(time_field, axis = 0)
print(perf[2, :])
start = end
end = end + 1 * 4
time_field = time[start:end].reshape([1, 4])
perf[3, :] = dims['qmcpack'] * 4.0 / 1024 / 1024 / np.mean(time_field, axis = 0)
print(perf[3, :])
np.savetxt("perf.txt", perf)

