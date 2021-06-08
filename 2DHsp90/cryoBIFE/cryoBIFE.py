#### CRYOBIFE Routine ###
### Run
### python3.6 cryoBIFE.py [FILE_PATH] > CryoBIFE-MCMC
 

import time
import numpy as np
import random as r
import scipy.special as ss
import matplotlib.pyplot as plt
from numba import config,jit ,njit, threading_layer
import sys

start_time = time.time()

################ Part 1: Defining functions #########################
config.THREADING_LAYER = 'threadsafe'
@njit
## Alex's routine
def w_logsumexp(x, weight, maxalpha):  ##To calculate the weighted logsumexp (Alex's idea). In our case, x is the
                              ##not-normalized BioEM values of the path, and the weights G that max the likelihod from cryo-BIFE.    
    C = np.zeros(x.shape[0])
    D = np.zeros(x.shape[0])
    norm = np.sum(np.exp(-weight))
    rho = np.exp(-weight)/norm

### No need for doing this at every MC step::: values sorted in maxalpha
    #for i in range(x.shape[0]):

    #    C[i] = np.max(x[i][:])

    for j in range(x.shape[0]):

        D[j] = np.log(np.sum(rho*np.exp(x[j][:] - maxalpha[j]))) + maxalpha[j]

    return (D)

config.THREADING_LAYER = 'threadsafe'
@njit(parallel=True)
#@jit
### Routine for initial path: 
def initial_val(path_ini, all_images):
    
    nframes = all_images.shape[0]
    
    process_path = np.zeros((nframes,len(path_ini)))
    G_old = np.zeros(len(path_ini))
    
    c = 0

    ini_path = np.zeros((nframes,len(path_ini)))

## Getting BioEM probs for the specificed path
    for i in path_ini:

        for j in range(nframes):

            ini_path[j][c] = all_images[j][i-1]

        c+=1
        
    Max_alpha = np.zeros((nframes))
    
## Getting Maximum of each image (to be used in logsumexp)
    for m in range(nframes):
        
        Max_alpha[m] = np.max(ini_path[m][:])        

### logP-BioEM values for path without normalizing
    for k in range(nframes):

        process_path[k] = ini_path[k,:] #- np.max(ini_path[k,:]))

## Setting random-uniform initial G distribution
        
    for x in range(0,len(path_ini)):

        G_old[x] = 4.*r.random()-2.0
        
    return(process_path, G_old, Max_alpha)


config.THREADING_LAYER = 'threadsafe'

@njit(parallel=True)
#@jit
def Gibbs(trials, nmodels, nframes, mcsteps, BioEM, G_old, Max_alpha):
    
    
    G = np.zeros((mcsteps,nmodels))   #Matriz para guardar los valores aceptados de energia
                                      # en cada mc-step

    Acc = np.zeros((trials,mcsteps,nmodels)) #Guarda las matrices de G aceptadas
    
    Gstart = np.zeros((trials,nmodels)) #Matriz para guardar los valores iniciales de G
                                          #en cada uno de los procesos.
    Prob = np.zeros((trials,mcsteps))
    
    exProb = np.zeros((trials,mcsteps))
    
    
    for i in range(trials):
        
        pold = -9999999999

        Gold = G_old
        Gnew = np.zeros(nmodels)
        Gold_1 = np.zeros(nmodels)

        Gstart[i] = Gold        
        
        for mc in range(mcsteps): 

## Likelihood
            A = np.zeros(nframes)
## Log-L
            lnA = 0

            rnum2 = r.randint(0,nmodels)   #Numero aleatorio en (0,20) 

            dis = -0.5 + r.random() #Valor del desplazamiento

            for j in range (nmodels):

                if j==rnum2:

                    Gnew[j] = Gold[j] + dis

                else:

                    Gnew[j] = Gold[j]

### Normalizing to mean G = 0
            Gmean = np.sum(Gnew)/nmodels
            Gnew = Gnew - Gmean
      
            norm = np.sum(np.exp(-Gnew))   ##Calcula la normalizacion

            prior = 0
            for l in range (nmodels-1):
                prior += (Gnew[l+1]-Gnew[l])**2

### Alex's routine
            A = w_logsumexp(BioEM, Gnew, Max_alpha)

            lnA = np.sum(A)+np.log(1/prior**2)

### CHECK OLD 
#            D2 = np.zeros((nframes,nmodels))

#            for k in range(nframes):
#                 D2[k,:] = np.exp(BioEM[k,:] - Max_alpha[k])
 
#            A2=np.dot(D2,np.exp(-Gnew)/norm)
            
#            lnA2 = np.sum(Max_alpha+np.log(A2))+np.log(1/prior**2)

#            print('check',lnA,lnA2,lnA-lnA2)
#################################
#            lnA = lnA2                      

            rr=np.log(r.random())

            if (lnA > pold):

                Prob[i,mc] = lnA    #criterio de aceptacion de metropolis
                
                exProb[i,mc] = np.exp(lnA)

                pold = lnA

                Gold = np.copy(Gnew)

            elif (rr < -(pold-lnA)):

                Prob[i,mc] = lnA  #criterio de aceptacion de metropolis
                
                exProb[i,mc] = np.exp(lnA)

                pold = lnA

                Gold = np.copy(Gnew)

            else:

                Prob[i,mc] = pold   #Rechazo del valor nuevo
                
                exProb[i,mc] = np.exp(pold)

            G[mc] = Gold 
            print('Steps',i,mc,Prob[i,mc],lnA,np.log(1/prior**2))

        Acc[i] = G
        
    return(Acc, Prob)

################# MAIN ######################

print('Starting cryo-BIFE')

## Name of file with path node indexes (following Julian's notation)
path_ini = np.loadtxt(sys.argv[1])
path_ini = path_ini.astype(int)

#Non-normalized BioEM probabilities for all images and the 225 nodes 
all_images = np.loadtxt('all_images')

### Intializing routine
BioEM_ini, G_old, Max_alpha = initial_val(path_ini, all_images)

######### MC PARAMS ####
#MC runs
trials = 1
#number of nodes along the path
nmodels = BioEM_ini.shape[1]
#number of cryo-EM particles
nframes = BioEM_ini.shape[0]
#MC steps
mcsteps = 500000

Gb_ini, prob_ini = Gibbs(trials, nmodels, nframes, mcsteps, BioEM_ini,G_old,Max_alpha)

#### Printing Gs for postprocesing #####

f = open('Gs','w')

#for k in range(np.array(Gb_ini).shape[0]):
for k in range(0,mcsteps):
    f.write(str(k))

    for j in range(nmodels):

        f.write(' ' + str(Gb_ini[0][k][j]))

    f.write('\n')

f.close()


