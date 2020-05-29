import numpy as np
import sys
from matplotlib import pyplot as plt

SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 20
LEGEND_SIZE = 16

plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=LEGEND_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
data = np.fromfile(sys.argv[1], dtype=np.float32).reshape([512, 512, 512])
data = np.abs(data)
n = 256
max_v = np.max(data)
print(max_v)
scale = float(sys.argv[2])
hist_range = max_v * scale
for i in range(9):
    n_total = 2*n
    tmp = np.copy(data[:n_total, :n_total, :n_total])
    if n>1:
        tmp[:n, :n, :n] = -1
    tmp2 = tmp[tmp > -1]
    bin_height,bin_boundary = np.histogram(tmp2, range=(0, hist_range), bins=100)
    bin_boundary[:-1] = bin_boundary[:-1] / (hist_range) * scale
    width = bin_boundary[1]-bin_boundary[0]
    print(tmp2.size)
    if np.sum(bin_height) > 0:
        bin_height = bin_height/float(np.sum(bin_height))
    axs.bar(bin_boundary[:-1], bin_height,width = width, label='level{}'.format(i + 1) )
    n = n//2

axs.legend()
# plt.show()
axs.set_xlabel('Error normalized to value range')
axs.set_ylabel('Distribution normalized within level')
axs.legend(loc='upper right', bbox_to_anchor= (1.02, 1.02))
plt.tight_layout()
plt.savefig('err_dist.pdf', bbox_inches='tight')
