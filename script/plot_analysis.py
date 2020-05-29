import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys 

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


fields=['vx', 'temperature']
name_map = {}
name_map[fields[0]] = 'velocity_x'
name_map[fields[1]] = 'temperature'
mgard_decompose_time = {}
mgard_decompose_time[fields[0]] = [0, 57.1874, 70.4244, 74.2488]
mgard_decompose_time[fields[1]] = [0, 58.8711, 72.7401, 76.6835]
ourmethod_decompose_time = {}
ourmethod_decompose_time[fields[0]] = [0, 2.25921, 2.5132, 2.52631]
ourmethod_decompose_time[fields[1]] = [0, 2.26155, 2.50248, 2.52636]

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8,4))
width = 0.2  # the width of the bars
labels = ['0', '1', '2', '3']
x = np.arange(len(labels))  # the label locations

ax=[axs[0], axs[1]]

i = 0
for field in fields:
    data = np.loadtxt("./result/parallel/level_scalability_{}.txt".format(field))[0]
    t_mgard = data + mgard_decompose_time[field]
    t_om = data + ourmethod_decompose_time[field]
    scale = np.loadtxt("./result/parallel/resource_scalability_{}.txt".format(field))[0]
    

    rects1 = ax[i].bar(x - 0.5 * width, t_mgard, width, label='MGARD')
    rects2 = ax[i].bar(x + 0.5 * width, t_om, width, label='OurMethod')

    l0 = ax[i].axhline(y=scale[0], color='k', linestyle='--')
    l1 = ax[i].axhline(y=scale[1], color='g', linestyle='--')
    l2 = ax[i].axhline(y=scale[2], color='r', linestyle='--')
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax[i].set_ylabel('Analysis Time (s)')
    ax[i].set_ylim(0, 150)
    ax[i].set_title(name_map[field])
    ax[i].set_xticks(x)
    ax[i].set_xticklabels(labels)
    ax[i].set_xlabel('Level of Decomposition')
    ax[i].legend(loc='upper right', ncol=1, bbox_to_anchor= (1, 0.95))
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
plt.savefig('analysis.pdf')
    #plt.show()