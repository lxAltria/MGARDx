import os
import numpy as np
import sys

target_level = int(sys.argv[1])
block_size = int(sys.argv[2])
opt = int(sys.argv[3])
cmd = "./test/test_convergence Uf48.bin.dat {} 3 100 500 500 {} {} > tmp.log".format(target_level, block_size, opt)
# cmd = "./test/test_convergence velocity_x.dat {} 3 512 512 512 {} {} > tmp.log".format(target_level, block_size, opt)
# cmd_grep = "grep \"rmse\" tmp.log | grep -o \"[0-9].*\" > tmp.result"
cmd_grep = "grep \"iso_change_count\" tmp.log | grep -o \"[0-9].*\" > tmp.result"

os.system(cmd)
os.system(cmd_grep)

if opt == 0:
	difference = np.loadtxt("tmp.result").reshape([target_level + 1, -1])
	np.save("errors.npy", difference)
if opt == 1:
	difference = np.loadtxt("tmp.result").reshape([target_level, -1])
	np.save("difference.npy", difference)

# errors = np.load("errors.npy")
# errors.shape
# difference = np.load("difference.npy")
# difference.shape
# convergence = difference[1:]/difference[:-1]

# num_result, num_sample = difference.shape
# convergence = np.zeros([num_result - 1, num_sample])
# for i in range(num_result - 1):
# 	for j in range(num_sample):
# 		if difference[i, j] == 0:
# 			if difference[i+1, j] != 0:
# 				convergence[i, j] = 1
# 			else:
# 				convergence[i, j] = 0
# 		else:
# 			convergence[i, j] = difference[i+1, j] / difference[i, j]

# inds = np.zeros([num_sample], dtype=np.int32)
# for i in range(num_sample):
# 	inds_i = np.where(difference[:, i] == 0)
# 	if inds_i[0].size == 0:
# 		inds[i] = num_result - 1
# 	else:
# 		inds[i] = inds_i[0][0]

# for i in range(target_level):
# 	convergence[i] = errors[i+1] / errors[i]
# n1, n2 = 10, 20
# def plot(n1, n2):
# 	num_levels, num_samples = errors.shape
# 	print("actual errors:")
# 	print(errors[:, n1])
# 	print(errors[:, n2])
# 	print("adjacent difference:")
# 	print(difference[:, n1])
# 	print(difference[:, n2])
# 	print("convergence rate:")
# 	print(convergence[:, n1])
# 	print(convergence[:, n2])
# 	fig,ax = plt.subplots()
# 	ax.plot(np.arange(num_levels), errors[:, n1], 'r--')
# 	ax.plot(np.arange(num_levels), errors[:, n2], 'g--')
# 	ax2=ax.twinx()
# 	ax2.plot(np.arange(num_levels - 2) + 2, convergence[:, n1], 'r-')
# 	ax2.plot(np.arange(num_levels - 2) + 2, convergence[:, n2], 'g-')
# 	plt.show()

# def plot(n1):
# 	num_levels, num_samples = errors.shape
# 	print("actual errors:")
# 	print(errors[:, n1])
# 	print("adjacent difference:")
# 	print(difference[:, n1])
# 	print("convergence rate:")
# 	print(convergence[:, n1])
# 	fig,ax = plt.subplots()
# 	l1 = ax.plot(np.arange(num_levels), errors[:, n1], 'r-', label="actual difference")
# 	l2 = ax.plot(np.arange(num_levels - 1) + 1, difference[:, n1], 'g--', label="adjacent difference")
# 	ax.set_xlabel("Level")
# 	ax.set_ylabel("Error")
# 	ax2 = ax.twinx()
# 	l3 = ax2.plot(np.arange(num_levels - 2) + 2, convergence[:, n1], 'b:', label="convergence")
# 	ax2.set_ylabel("Convergence")
# 	ax2.set_ylim(0, 1)
# 	lines = l1 + l2 + l3
# 	labels = [l.get_label() for l in lines]
# 	ax.legend(lines, labels)
# 	plt.title('block_{}'.format(n1))
# 	# plt.show()
# 	plt.savefig('block_{}.png'.format(n1), bbox_inches='tight')
# 	plt.close()

# num_levels, num_samples = convergence.shape
# for i in range(num_samples):
# 	plt.plot(np.arange(num_levels), convergence[:, i])

