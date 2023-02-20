# SNMF_of_HIC_data
# This code inputs a .hic file and an positive integer n. It calculates the optimal symmetric non-negative rank r approximation of HiC contact 
# matrix for all r less than or equal to n. The results are displayed as a r(x-axis)-by-"error of optimal-rank-r approximation"(y-axis) graph. 

#the low rank approximation is based off code from https://github.com/hiroyuki-kasai/NMFLibrary
#hic-straw is based off code from https://github.com/aidenlab/straw


#installation requirements. need to run these install commands in the terminal before running the program
#python3 -m pip install hic-straw
#python3 -m pip install numpy
#python3 -m pip install scipy
#python3 -m pip install matlabengine //need to have matlab installed on your computer first

import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import hicstraw
import seaborn as sns
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap






# hard-coded .hic file to practice with
hic = hicstraw.HiCFile("https://www.encodeproject.org/files/ENCFF718AWL/@@download/ENCFF718AWL.hic")

#matrix zoom data object
matrix_object_chr4 = hic.getMatrixZoomData('4', '4', "observed", "KR", "BP", 5000)


#get numpy matrix
numpy_matrix = mzd.getRecordsAsMatrix(10000000, 12000000, 10000000, 12000000)




# generate random matrix to test procedure on
m = 500
n = 100
V = np.random.random_sample((m,n))

# define rank to be factorized
rank = 5

# set options for solvers
options = dict()
options['verbose'] = '1'
options['max_epoch'] = '100'

# import MATLAB module
import matlab.engine
eng = matlab.engine.start_matlab()

# add MATLAB path
eng.run_me_first(0, nargout=0)

# convert numpy array to matlab.double
V_m = matlab.double(V.tolist())

# convert dictioary to matlab.struct
options_m = eng.struct(options);

# perform solvers in MATLAB 
# Fro-MU
[w_mu, infos_mu] = eng.fro_mu_nmf(V_m, rank, options_m, nargout=2)
# HALS
[w_hals, infos_hals] = eng.als_nmf(V_m, rank, options_m, nargout=2)
# ACC-HALS
options['alg'] = 'acc_hals'
options_m = eng.struct(options);
[w_acc_hals, infos_acc_hals] = eng.als_nmf(V_m, rank, options_m, nargout=2)

# convert matlab.struct to list
iter_mu = list(infos_mu['iter'][0])
time_mu = list(infos_mu['time'][0])
cost_mu = list(infos_mu['cost'][0])
iter_hals = list(infos_hals['iter'][0])
time_hals = list(infos_hals['time'][0])
cost_hals = list(infos_hals['cost'][0])
iter_acc_hals = list(infos_acc_hals['iter'][0])
time_acc_hals = list(infos_acc_hals['time'][0])
cost_acc_hals = list(infos_acc_hals['cost'][0])

# plotting
plt.figure()
plt.xlabel('Epoch')
plt.ylabel('Cost')
plt.plot(iter_mu, cost_mu, label ="Fro-MU")
plt.plot(iter_hals, cost_hals, label ="HALS")
plt.plot(iter_acc_hals, cost_acc_hals, label ="ACC-HALS")
plt.legend() 
plt.title('Epoch vs. Cost')

plt.figure()
plt.xlabel('Time [sec]')
plt.ylabel('Cost')
plt.plot(time_mu, cost_mu, label ="Fro-MU")
plt.plot(time_hals, cost_hals, label ="HALS")
plt.plot(time_acc_hals, cost_acc_hals, label ="ACC-HALS")
plt.legend() 
plt.title('Time vs. Cost')

plt.show()

eng.quit()
