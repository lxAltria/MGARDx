import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys 

SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 14

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
sizes['hurricane'] = Hurricane_size * 13
sizes['nyx'] = NYX_size * 6
sizes['scale'] = SCALE_size * 12
sizes['qmc'] = QMCPACK_size * 1

num_opt = 4


datasets=['hurricane', 'nyx', 'scale', 'qmc']
name_map = {}
name_map[datasets[0]] = 'Hurricane'
name_map[datasets[1]] = 'NYX'
name_map[datasets[2]] = 'SCALE-LETKF'
name_map[datasets[3]] = 'QMCPACK'

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,6))
width = 0.2  # the width of the bars
labels = ['0.1', '0.01', '0.001', '0.0001']
x = np.arange(len(labels))  # the label locations

ax=[axs[0,0], axs[0,1], axs[1,0], axs[1,1]]

i = 0
for dataset in datasets:
    time = np.loadtxt('./result/rate/{}_decomp_time.txt'.format(dataset))
    perfs = sizes[dataset] / time
    print(perfs)
    
    rects1 = ax[i].bar(x - 1.5 * width, perfs[0], width, label='SZ')
    rects2 = ax[i].bar(x - 0.5 * width, perfs[1], width, label='ZFP')
    rects3 = ax[i].bar(x + 0.5 * width, perfs[2], width, label='MGARD')
    rects4 = ax[i].bar(x + 1.5 * width, perfs[3], width, label='OurMethod')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax[i].set_ylabel('Decompression Perf. (MB/s)')
    ax[i].set_ylim(0, 800)
    ax[i].set_title(name_map[dataset])
    ax[i].set_xticks(x)
    ax[i].set_xticklabels(labels)
    ax[i].set_xlabel('Error Bound')
    ax[i].legend()
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
plt.savefig('decomp_perf.pdf')
    #plt.show()