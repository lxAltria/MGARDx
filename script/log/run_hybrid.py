import numpy as np
import utils

#eb = np.array([0.5, 0.1, 0.05, 0.01, 0.005, 0.003, 0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000005, 0.000001])
eb = np.array([0.1, 0.01, 0.001, 0.0001])

dims = np.array([100, 500, 500])
utils.run_hybrid_model("Hurricane", dims, eb)

dims = np.array([512, 512, 512])
utils.run_hybrid_model("NYX", dims, eb)

dims = np.array([98, 1200, 1200])
utils.run_hybrid_model("SCALE", dims, eb)

dims = np.array([33120, 69, 69])
utils.run_hybrid_model("QMCPACK", dims, eb)

