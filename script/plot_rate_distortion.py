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




# input dimension: num_fields x num_cases
def get_total_rate_distortion(nrmse, cr):
    if len(nrmse.shape) == 1:
        nrmse = nrmse.reshape([1, -1])
        cr = cr.reshape([1, -1])
    total_nrmse = np.sqrt(np.mean(nrmse**2, axis=0))
    total_psnr = -20 * np.log10(total_nrmse)
    total_br = np.mean(32.0 / cr, axis=0)
    return total_br, total_psnr

def load_rate_distortion(dataset, compressor):
    nrmse = np.loadtxt("result/{}_{}_nrmse.txt".format(dataset, compressor))
    ratio = np.loadtxt("result/{}_{}_ratio.txt".format(dataset, compressor))
    if compressor == 'old_mgard':
        nrmse2 = np.loadtxt("result/{}_{}_nrmse.txt.cont".format(dataset, compressor))
        ratio2 = np.loadtxt("result/{}_{}_ratio.txt.cont".format(dataset, compressor))
        # print(ratio.shape)
        # print(ratio2.shape)
        if len(ratio.shape) == 1:
            num_field = 1
            nrmse2 = nrmse2[-2:]
            ratio2 = ratio2[-2:]
            num_sample = ratio.size
            num_sample2 = ratio2.size
        else:
            num_field, num_sample = ratio.shape
            nrmse2 = nrmse2[:, -2:]
            ratio2 = ratio2[:, -2:]
            _, num_sample2 = ratio2.shape
        tmp_nrmse = np.zeros([num_field, num_sample + num_sample2])
        tmp_nrmse[:, :num_sample] = nrmse
        tmp_nrmse[:, num_sample:] = nrmse2
        tmp_ratio = np.zeros([num_field, num_sample + num_sample2])
        tmp_ratio[:, :num_sample] = ratio
        tmp_ratio[:, num_sample:] = ratio2
        nrmse = tmp_nrmse
        ratio = tmp_ratio        
    return get_total_rate_distortion(nrmse, ratio)


datasets=['Hurricane', 'NYX', 'SCALE', 'QMCPACK']
name_map = {}
name_map[datasets[0]] = 'Hurricane'
name_map[datasets[1]] = 'NYX'
name_map[datasets[2]] = 'SCALE-LETKF'
name_map[datasets[3]] = 'QMCPACK'

compressors=['sz', 'zfp', 'old_mgard', 'mgard_level1']
compressor_name_map=['SZ', 'ZFP', 'MGARD', 'OurMethod']
styles=['b-^', 'g--s','c-.+', 'r:x']
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,6))


ax=[axs[0,0], axs[0,1], axs[1,0], axs[1,1]]

i = 0
for dataset in datasets:
    j = 0
    for compressor in compressors:
        br, psnr = load_rate_distortion(dataset, compressor)

        p, = ax[i].plot(br, psnr, styles[j], label='{}'.format(compressor_name_map[j]))
        j += 1
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax[i].set_ylabel('PSNR')
    ax[i].set_ylim(top=90)
    ax[i].set_xlim(0, 1.6)
    ax[i].grid(which='major', axis='y')
    ax[i].set_title(name_map[dataset])

    ax[i].set_xlabel('Bit-rate')
    ax[i].legend(loc='lower right')
    i += 1


plt.tight_layout()
plt.savefig('rate_distortion.pdf')
    #plt.show()