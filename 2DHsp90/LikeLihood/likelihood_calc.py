### CODE TO CALCULATE THE LIKELIHOOD GIVEN AN INPUT Gmax and PATH
### Run:
### python3.6 likelihood_calc.py [FILE_PATH] [FILE_GMAX] > Log-L-image 


import numpy as np
import random as r
import scipy.special as ss
from numba import jit ,njit
import sys

@njit
## Alex's routine
def w_logsumexp(x, weight):  ##To calculate the weighted logsumexp (Alex's idea). In our case, x is the
                              ##not-normalized BioEM values of the path, and the weights G that max the likelihod from cryo-BIFE.    
    C = np.zeros(x.shape[0])
    D = np.zeros(x.shape[0])
    norm = np.sum(np.exp(-weight))
    rho = np.exp(-weight)/norm
    
    for i in range(x.shape[0]):
        
        C[i] = np.max(x[i][:])
        
    for j in range(x.shape[0]):
        
        D[j] = np.log(np.sum(rho*np.exp(x[j][:] - C[j]))) + C[j]
        ## Printing the likelihood for each image, and separated into terms.
        print('likeliI', j,D[j],np.log(np.sum(rho*np.exp(x[j][:] - C[j]))),C[j])
        
    return (D)

@njit
def initial_val(path_ini,all_images): ## Get the initial path with BioEM values

    nframes = all_images.shape[0]

    c = 0

    ini_path = np.zeros((nframes,len(path_ini)))

    for i in path_ini:

        for j in range(nframes):

            ini_path[j][c] = all_images[j][i-1]

        c+=1

    return(ini_path)

#########################################################
#### MAIN

### Reading input

## Name of file with path node indexes (following Julian's notation)
path_ini = np.loadtxt(sys.argv[1])
path_ini = path_ini.astype(int)
## Name of file with Gs that maximize the likelihood of that path 
Gtest = np.loadtxt(sys.argv[2])

## Non-normalized BioEM probabilities for all images and the 225 nodes 
all_images = np.loadtxt('all_images') 

### Getting BioEM values of that path
ini_pathb = initial_val(path_ini,all_images)

### Logsumexp; returns a vector with the likelihood of each image
log_Likelihood_alphab = w_logsumexp(ini_pathb, Gtest)

## Summing up the likelihoods
log_c_alpha = np.sum(log_Likelihood_alphab)

print('\n Likelihood path = ',log_c_alpha)
