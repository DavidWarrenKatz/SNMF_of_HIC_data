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

'''
On the bmi cluster, I received the following error when trying to install hic-straw and numpy:

python3 -m pip install hic-straw
Collecting hic-straw
  Retrying (Retry(total=4, connect=None, read=None, redirect=None, status=None)) after connection broken by 'ProxyError('Cannot connect to proxy.', OSError('Tunnel connection failed: 407 Proxy Authentication Required',))': /simple/hic-straw/
  Retrying (Retry(total=3, connect=None, read=None, redirect=None, status=None)) after connection broken by 'ProxyError('Cannot connect to proxy.', OSError('Tunnel connection failed: 407 Proxy Authentication Required',))': /simple/hic-straw/
  Retrying (Retry(total=2, connect=None, read=None, redirect=None, status=None)) after connection broken by 'ProxyError('Cannot connect to proxy.', OSError('Tunnel connection failed: 407 Proxy Authentication Required',))': /simple/hic-straw/
  Retrying (Retry(total=1, connect=None, read=None, redirect=None, status=None)) after connection broken by 'ProxyError('Cannot connect to proxy.', OSError('Tunnel connection failed: 407 Proxy Authentication Required',))': /simple/hic-straw/
  Retrying (Retry(total=0, connect=None, read=None, redirect=None, status=None)) after connection broken by 'ProxyError('Cannot connect to proxy.', OSError('Tunnel connection failed: 407 Proxy Authentication Required',))': /simple/hic-straw/
  Could not find a version that satisfies the requirement hic-straw (from versions: )
No matching distribution found for hic-straw
'''





# hard-coded .hic file to practice with
hic = hicstraw.HiCFile("https://www.encodeproject.org/files/ENCFF718AWL/@@download/ENCFF718AWL.hic")

'''
This is where the program is crashing locally. I received the following comment:

Python(10695,0x11ac0c600) malloc: Heap corruption detected, free list is damaged at 0x600001899d20
*** Incorrect guard value: 278612453268736
Python(10695,0x11ac0c600) malloc: *** set a breakpoint in malloc_error_break to debug
zsh: abort      python3

#matrix zoom data object
matrix_object_chr4 = hic.getMatrixZoomData('4', '4', "observed", "KR", "BP", 5000)
'''
#get the input n from user
n = int(input("Input an integer n: "))

#get numpy matrix
numpy_matrix = mzd.getRecordsAsMatrix(10000000, 12000000, 10000000, 12000000)

# import MATLAB module
import matlab.engine
eng = matlab.engine.start_matlab()

# add MATLAB path
eng.run_me_first(0, nargout=0)

for i in range(n):
    # define rank to be factorized
    rank = i

    # set options for symmetric solver
    options = dict()
    options['verbose'] = '1'
    options['max_epoch'] = '100'
    options.calc_symmetry = true;    #this line might only work in the matlab code


    # convert numpy array to matlab.double
    V_m = matlab.double(V.tolist())

    # convert dictioary to matlab.struct
    options_m = eng.struct(options);

    # perform solvers in MATLAB
    # Fro-MU
    [w_mu, infos_mu] = eng.fro_mu_nmf(V_m, rank, options_m, nargout=2)


    
    ## perform symmetric factroization
    # Symm-ANLS
    options.alpha = 0.6;
    [w_symm_anls, infos_symm_anls] = eng.symm_anls(V, rank, options);
    # Symm-Newton
    [w_symm_newton, infos_symm_newton] = eng.symm_newton(V, rank, options);
    # Symm-Hals
    options.lambda = 0.6;
    [w_symm_halsacc, infos_symm_halsacc] = eng.symm_halsacc(V, rank, options);
    
    # convert matlab.struct to list for symmetric alg
    iter_anls = list(infos_symm_anls['iter'][0])
    cost_anls = list(infos_symm_anls['cost'][0])
    iter_newton = list(infos_symm_newton['iter'][0])
    cost_newton = list(infos_symm_newton['cost'][0])
    iter_hals = list(infos_symm_halsacc['iter'][0])
    cost_hals = list(infos_symm_halsacc['cost'][0])
    

    # plotting
    plt.figure()
    plt.xlabel('Epoch')
    plt.ylabel('Cost')
    plt.plot(iter_anls, cost_anls, label ="norm(W = Wt)")
    plt.legend()
    plt.title(f'Epoch vs. Cost for {i}')

    plt.show()




eng.quit()
