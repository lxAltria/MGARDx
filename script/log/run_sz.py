import numpy as np
import utils

eb = np.array([0.5, 0.1, 0.05, 0.01, 0.005, 0.003, 0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000005, 0.000001])

dims = np.array([100, 500, 500])
utils.run_sz("Hurricane/step48", dims, eb)

dims = np.array([512, 512, 512])
utils.run_sz("NYX", dims, eb)

dims = np.array([98, 1200, 1200])
utils.run_sz("SCALE", dims, eb)

dims = np.array([33120, 69, 69])
utils.run_sz("QMCPACK", dims, eb)

