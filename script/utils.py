import numpy as np 

# input dimension: num_fields x num_cases
def get_total_rate_distortion(nrmse, cr):
    if len(nrmse.shape) == 1:
        nrmse = nrmse.reshape([1, -1])
        cr = cr.reshape([1, -1])
    total_nrmse = np.sqrt(np.mean(nrmse**2, axis=0))
    total_psnr = -20 * np.log10(total_nrmse)
    total_br = np.mean(32.0 / cr, axis=0)
    return total_br, total_psnr

def PSNR(data, dec_data):
	value_range = np.max(data) - np.min(data)
	diff = data - dec_data
	nrmse = np.sqrt(np.mean(diff**2)) / value_range
	psnr = 20 * np.log10(1.0 / nrmse)
	return psnr, nrmse

def normalize(data):
	min_data = np.min(data)
	max_data = np.max(data)
	value_range = max_data - min_data
	data = (data - min_data) / value_range
	return data, min_data, value_range

def denormalize(data, min_data, value_range):
	data = data * value_range + min_data
	return data

def load_rate_distorition(dataset, compressor):
    nrmse = np.loadtxt("result/{}_{}_nrmse.txt".format(dataset, compressor))
    ratio = np.loadtxt("result/{}_{}_ratio.txt".format(dataset, compressor))
    return get_total_rate_distortion(nrmse, ratio)

def load_rate_distorition_given_field(dataset, compressor, field):
    psnr = np.loadtxt("result/{}_{}_psnr.txt".format(dataset, compressor))[field]
    ratio = np.loadtxt("result/{}_{}_ratio.txt".format(dataset, compressor))[field]
    return 32.0/ratio, psnr

from matplotlib import pyplot as plt
def plot_rate_distortion(dataset, compressor):
    br, psnr = load_rate_distorition(dataset, compressor)
    plt.plot(br, psnr, label='{}'.format(compressor))

def plot_mgard_all_levels(dataset, levels):
    for i in range(levels):
        br, psnr = load_rate_distorition(dataset, 'mgard_level{}'.format(i+1))
        plt.plot(br, psnr, label='level{}'.format(i+1))

def plot_mgard_given_level(dataset, level):
    br, psnr = load_rate_distorition(dataset, 'mgard_level{}'.format(level))
    plt.plot(br, psnr, label='level{}'.format(level))

def plot_mgard_given_level_given_field(dataset, level, field):
    br, psnr = load_rate_distorition_given_field(dataset, 'mgard_level{}'.format(level), field)
    plt.plot(br, psnr, label='level{}'.format(level))

from os import system
from os import listdir
from os.path import getsize
def run_mgard(folder, dims, level=20):
    dataset_name = folder[folder.rfind('/') + 1:]
    if dataset_name == 'step48':
        dataset_name = 'Hurricane'
    comp_exec='/Users/xin/github/MGARD/build/test/test_compress'
    decomp_exec='/Users/xin/github/MGARD/build/test/test_decompress'
    data_files = sorted([f for f in listdir(folder) if f.endswith(".dat")])
    print(data_files)
    ebs = np.array([1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001])
    num_fields = len(data_files)
    num_eb = len(ebs)
    psnr = np.zeros([num_fields, num_eb])
    nrmse = np.zeros([num_fields, num_eb])
    ratio = np.zeros([num_fields, num_eb])
    for i in range(num_fields):
        file = data_files[i]
        filename = "{}/{}".format(folder, file)
        filename_comp = "{}/{}.mgard".format(folder, file)
        filename_decomp = "{}/{}.mgard.out".format(folder, file)
        data = np.fromfile(filename, dtype=np.float32)
        value_range = np.max(data) - np.min(data)
        for j in range(num_eb):
            eb = ebs[j]
            system("{} {} 0 {} {} 1 3 {} {} {}".format(comp_exec, filename, eb * value_range, level, dims[0], dims[1], dims[2]))
            system("{} {} {}.mgard 0 3 {} {} {}".format(decomp_exec, filename, filename, dims[0], dims[1], dims[2]))
            dec_data = np.fromfile(filename_decomp, dtype=np.float32)
            psnr[i, j], nrmse[i, j] = PSNR(data, dec_data)
            ratio[i, j] = getsize(filename) * 1.0 / getsize(filename_comp)
    np.savetxt("result/{}_mgard_level{}_psnr.txt".format(dataset_name, level), psnr)
    np.savetxt("result/{}_mgard_level{}_nrmse.txt".format(dataset_name, level), nrmse)
    np.savetxt("result/{}_mgard_level{}_ratio.txt".format(dataset_name, level), ratio)
    br_overall, psnr_overall = get_total_rate_distortion(nrmse, ratio)
    print(br_overall)
    print(psnr_overall)

def run_zfp(folder, dims):
    dataset_name = folder[folder.rfind('/') + 1:]
    if dataset_name == 'step48':
        dataset_name = 'Hurricane'
    comp_exec='/Users/xin/github/zfp/bin/zfp'
    data_files = sorted([f for f in listdir(folder) if f.endswith(".dat")])
    print(data_files)
    ebs = np.array([2, 1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005])
    num_fields = len(data_files)
    num_eb = len(ebs)
    psnr = np.zeros([num_fields, num_eb])
    nrmse = np.zeros([num_fields, num_eb])
    ratio = np.zeros([num_fields, num_eb])
    for i in range(num_fields):
        file = data_files[i]
        filename = "{}/{}".format(folder, file)
        filename_comp = "{}/{}.zfp".format(folder, file)
        filename_decomp = "{}/{}.zfp.out".format(folder, file)
        data = np.fromfile(filename, dtype=np.float32)
        value_range = np.max(data) - np.min(data)
        for j in range(num_eb):
            eb = ebs[j]
            system("{} -f -i {} -z {} -o {} -a {} -3 {} {} {} -s".format(comp_exec, filename, filename_comp, filename_decomp, eb * value_range, dims[2], dims[1], dims[0]))
            dec_data = np.fromfile(filename_decomp, dtype=np.float32)
            psnr[i, j], nrmse[i, j] = PSNR(data, dec_data)
            ratio[i, j] = getsize(filename) * 1.0 / getsize(filename_comp)
    np.savetxt("result/{}_zfp_psnr.txt".format(dataset_name), psnr)
    np.savetxt("result/{}_zfp_nrmse.txt".format(dataset_name), nrmse)
    np.savetxt("result/{}_zfp_ratio.txt".format(dataset_name), ratio)
    br_overall, psnr_overall = get_total_rate_distortion(nrmse, ratio)
    print(br_overall)
    print(psnr_overall)

def run_sz(folder, dims):
    dataset_name = folder[folder.rfind('/') + 1:]
    if dataset_name == 'step48':
        dataset_name = 'Hurricane'
    comp_exec='/Users/xin/utils/sz_master/bin/sz'
    data_files = sorted([f for f in listdir(folder) if f.endswith(".dat")])
    print(data_files)
    ebs = np.array([0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005])
    num_fields = len(data_files)
    num_eb = len(ebs)
    psnr = np.zeros([num_fields, num_eb])
    nrmse = np.zeros([num_fields, num_eb])
    ratio = np.zeros([num_fields, num_eb])
    for i in range(num_fields):
        file = data_files[i]
        filename = "{}/{}".format(folder, file)
        filename_comp = "{}/{}.sz".format(folder, file)
        filename_decomp = "{}/{}.sz.out".format(folder, file)
        data = np.fromfile(filename, dtype=np.float32)
        value_range = np.max(data) - np.min(data)
        for j in range(num_eb):
            eb = ebs[j]
            system("{} -z -f -i {} -M ABS -A {} -3 {} {} {}".format(comp_exec, filename, eb * value_range, dims[2], dims[1], dims[0]))
            system("{} -x -f -i {} -s {} -a -3 {} {} {}".format(comp_exec, filename, filename_comp, dims[2], dims[1], dims[0]))
            dec_data = np.fromfile(filename_decomp, dtype=np.float32)
            psnr[i, j], nrmse[i, j] = PSNR(data, dec_data)
            ratio[i, j] = getsize(filename) * 1.0 / getsize(filename_comp)
    np.savetxt("result/{}_sz_psnr.txt".format(dataset_name), psnr)
    np.savetxt("result/{}_sz_nrmse.txt".format(dataset_name), nrmse)
    np.savetxt("result/{}_sz_ratio.txt".format(dataset_name), ratio)
    br_overall, psnr_overall = get_total_rate_distortion(nrmse, ratio)
    print(br_overall)
    print(psnr_overall)
