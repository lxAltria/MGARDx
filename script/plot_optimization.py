import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys 

TINY_SIZE = 8
SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

Hurricane_size = 100 * 500 * 500 * 4.0 / 1024 / 1024
NYX_size = 512 * 512 * 512 * 4.0 / 1024 / 1024
SCALE_size = 98 * 1200 * 1200 * 4.0 / 1024 / 1024
QMCPACK_size = 288 * 115 * 69 * 69 * 4.0 / 1024 / 1024

sizes = {}
sizes['Hurricane'] = Hurricane_size
sizes['NYX'] = NYX_size
sizes['SCALE'] = SCALE_size
sizes['QMCPACK'] = QMCPACK_size

num_opt = 4
def read_time_opt(filename, dataset):
    time = np.loadtxt(filename).reshape([num_opt, -1])
    # print(time.shape)
    _, num_var = time.shape
    total_time = np.sum(time, axis=1)
    perf = sizes[dataset] * num_var / total_time
    return perf

ops=['decomposition', 'recomposition']
name_map = {}
name_map[ops[0]] = 'Decomposition'
name_map[ops[1]] = 'Recomposition'

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8,5))
width = 0.2  # the width of the bars
labels = ['Hurricane', 'NYX', 'SCALE-LETKF', 'QMCPACK']
x = np.arange(len(labels))  # the label locations

ax=[axs[0], axs[1]]

i = 0
for op in ops:
    perfs = np.zeros([4, 4])
    perfs[:, 0] = read_time_opt('./result/optimizations/hurricane_{}_time.txt'.format(op), 'Hurricane')
    perfs[:, 1] = read_time_opt('./result/optimizations/nyx_{}_time.txt'.format(op), 'NYX')
    perfs[:, 2] = read_time_opt('./result/optimizations/scale_{}_time.txt'.format(op), 'SCALE')
    perfs[:, 3] = read_time_opt('./result/optimizations/qmc_{}_time.txt'.format(op), 'QMCPACK')
    rects1 = ax[i].bar(x - 1.5 * width, perfs[0], width, label='MGARD')
    rects2 = ax[i].bar(x - 0.5 * width, perfs[1], width, label='OurMethod w/ opt. ABC')
    rects3 = ax[i].bar(x + 0.5 * width, perfs[2], width, label='OurMethod w/ opt. ABCD')
    rects3 = ax[i].bar(x + 1.5 * width, perfs[3], width, label='OurMethod w/ opt. ABCDE')
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax[i].set_ylabel('Performance (MB/s)')
    ax[i].set_ylim(0, 350)
    ax[i].set_title(name_map[op])
    ax[i].set_xticks(x)
    ax[i].set_xticklabels(labels)
    ax[i].tick_params(axis='x', rotation=45)
    # ax[i].set_xlabel('Level of Decomposition')
    ax[i].legend(loc='upper right', ncol=1, bbox_to_anchor= (1, 1.02))
    i += 1

# def autolabel(rects):
#     """Attach a text label above each bar in *rects*, displaying its height."""
#     for rect in rects:
#         height = rect.get_height()
#         ax.annotate('{}'.format(height),
#                     xy=(rect.get_x() + rect.get_width() / 2, height),
#                     xytext=(0, 3),  # 3 points vertical offset
#                     textcoords="offset points",
#                     ha='center', va='bottom')
# autolabel(rects1)
# autolabel(rects2)
plt.tight_layout()
plt.savefig('optimizations.pdf')
    #plt.show()