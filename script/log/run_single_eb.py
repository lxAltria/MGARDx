import numpy as np
import sys
import utils

dataset = sys.argv[1]
eb = float(sys.argv[2])
if dataset == "Hurricane":
	dims = np.array([100, 500, 500])
elif dataset == "NYX":
	dims = np.array([512, 512, 512])
elif dataset == "SCALE":
	dims = np.array([98, 1200, 1200])
elif dataset == "QMCPACK":
	dims = np.array([33120, 69, 69])
else:
	print("no such dataset")
	exit(0)
ebs = np.array([eb])
method = sys.argv[3]
if method == "MGARD":
	utils.run_mgard(dataset, dims, ebs)
elif method == "SZ":
	utils.run_sz(dataset, dims, ebs)
elif method == "ZFP":
	utils.run_zfp(dataset, dims, ebs)
elif method == "HYBRID":
	utils.run_hybrid_model(dataset, dims, ebs)
else:
	print("no such option")
	exit(0)

