import numpy as np
from matplotlib import pyplot as plt

data = np.fromfile("decomposed.dat", dtype=np.float32).reshape([512, 512, 512])
data = np.abs(data)
n = 256
max_v = np.max(data)
print(max_v)
scale = 0.001
hist_range = max_v * scale
for i in range(9):
    n_total = 2*n
    tmp = np.copy(data[:n_total, :n_total, :n_total])
    if n>1:
        tmp[:n, :n, :n] = -1
    tmp2 = tmp[tmp > -1]
    print(tmp2.size)
    bin_height,bin_boundary = np.histogram(tmp2, range=(0, hist_range), bins=100)
    bin_boundary[:-1] = bin_boundary[:-1] / (hist_range) * scale
    width = bin_boundary[1]-bin_boundary[0]
    bin_height = bin_height/float(np.sum(bin_height))
    plt.bar(bin_boundary[:-1], bin_height,width = width, label='level{}'.format(i + 1) )
    n = n//2

plt.legend()
plt.show()