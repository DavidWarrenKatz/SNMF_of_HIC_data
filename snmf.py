# SNMF_of_HIC_data
# This code inputs a .hic file and an positive integer n. It calculates the optimal symmetric non-negative rank r approximation of HiC contact
# matrix for all r less than or equal to n. The results are displayed as a r(x-axis)-by-"error of optimal-rank-r approximation"(y-axis) graph.

import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import hicstraw
import seaborn as sns
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap

# hard-coded .hic file to practice with. This command returns a hicstraw.HiCFile object. This file contains hic data with many resolutions.
hic = hicstraw.HiCFile("https://www.encodeproject.org/files/ENCFF718AWL/@@download/ENCFF718AWL.hic")

#get the input n from user
n = int(input("Input an integer n: "))

#Get the matrix object for chromosome 4, no normalization, at 5kb base pair resolution. hicstraw.MatrixZoomData returns an object of type hicstraw.MatrixZoomData, a class defined in hic-straw
matrix_object_chr4 = hic.getMatrixZoomData('4', '4', "observed", "NONE", "BP", 5000)

#Get a numpy matrix for the loci between 10MB and 12MB. This is of type numpy.ndarray
#This is a 401 by 401 matrix, represneting two million base paris at 5000 bp resolution
V = matrix_object_chr4.getRecordsAsMatrix(10000000, 12000000, 10000000, 12000000)

# import MATLAB module
import matlab.engine
eng = matlab.engine.start_matlab()

# add MATLAB path
eng.run_me_first(0, nargout=0)

# convert numpy array to matlab.double
V_m = matlab.double(V.tolist())

# set options for symmetric solver
options = dict()
options['verbose'] = '1'
options['max_epoch'] = '500'
#options.calc_symmetry = true;    #this line might only work in the matlab code


# convert dictioary to matlab.struct
#options_m = eng.struct(options);

#create array of ranks
ranks = np.arange(0., n+1., 1)

#create an empy array of costs to be filled in the following loop
errors_mu = np.zeros(n+1)
errors_hals = np.zeros(n+1)
errors_symm_newton = np.zeros(n+1)
errors_acc_hals = np.zeros(n+1)

for i in range(n):
    # define rank to be factorized
    rank = i+1
    
    # perform solvers in MATLAB
    
    # Fro-MU
    options['alg'] = 'mu'
    options_m = eng.struct(options);
    [w_mu, infos_mu] = eng.fro_mu_nmf(V_m, rank, options_m, nargout=2)
    cost_mu = list(infos_mu['cost'][0])
    errors_mu[i] = cost_mu[0]
    
    # HALS
    options['alg'] = 'hals'
    options_m = eng.struct(options);
    [w_hals, infos_hals] = eng.als_nmf(V_m, rank, options_m, nargout=2)
    cost_hals = list(infos_hals['cost'][0])
    errors_hals[i] = cost_hals[0]
    
#    # Symm-Newton
#    options['calc_symmetry'] = True;
#    options['alg'] = 'symm_newton';
#    options['alpha'] = 0.6;
#    options_m = eng.struct(options);
#    [w_symm_newton, infos_symm_newton] = eng.symm_newton(V_m, rank, options_m);
#    cost_symm_newton = list(infos_symm_newton['cost'][0])
#    errors_symm_neton[i] = cost_symm_newton[0]
    
    # ACC-HALS
    options['alg'] = 'acc_hals'
    options_m = eng.struct(options);
    [w_acc_hals, infos_acc_hals] = eng.als_nmf(V_m, rank, options_m, nargout=2)
    cost_acc_hals = list(infos_acc_hals['cost'][0])
    errors_acc_hals[i] = cost_acc_hals[0]
    
    
    
# plotting
plt.figure()
plt.xlabel('Rank')
plt.ylabel('Cost')
plt.plot(ranks, errors_mu, label ="Fro-MU")
plt.plot(ranks, errors_hals, label ="HALS")
#plt.plot(ranks, errors_symm_newton, label ="Symm-Newton")
plt.plot(ranks, errors_hals, label ="ACC-HALS")
plt.legend()
plt.title('Rank vs. Cost')
plt.show()
    


eng.quit()
