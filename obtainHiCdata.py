import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import hicstraw
import seaborn as sns
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap
from scipy.io import savemat


# hard-coded .hic file to practice with. This command returns a hicstraw.HiCFile object. This file contains hic data with many resolutions.
hic = hicstraw.HiCFile("https://www.encodeproject.org/files/ENCFF718AWL/@@download/ENCFF718AWL.hic")

#Get the matrix object for chromosome 4, no normalization, at 5kb base pair resolution. hicstraw.MatrixZoomData returns an object of type hicstraw.MatrixZoomData, a class defined in hic-straw
matrix_object_chr4 = hic.getMatrixZoomData('4', '4', "observed", "NONE", "BP", 5000)

#Get a numpy matrix for the loci between 10MB and 12MB. This is of type numpy.ndarray
#This is a 401 by 401 matrix, represneting two million base paris at 5000 bp resolution
V = matrix_object_chr4.getRecordsAsMatrix(10000000, 12000000, 10000000, 12000000)

#save the numpy array into a dictionary to use in the savemat function
mdic = {"V": V}

savemat("matlabfile.mat", mdic)
