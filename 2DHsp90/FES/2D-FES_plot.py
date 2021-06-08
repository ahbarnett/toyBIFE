### CODE TO PLOT 2D FES WITH PATHS 
### RUN WITH python3.6 
### Num_Images_grid has to be in the same folde


import numpy as np
import matplotlib.pyplot as plt

z = np.loadtxt('Num_Images_grid') 

models_grid = np.zeros((15,15)) #Matrix where each value is an index of particular node over the 2D-Free Energy Surface (2D-FES)
kappa=10

for i in range(1,16):
    
    for j in range(1,16):
        
        models_grid[i-1][j-1] = (15-i)*15 + j
        
points_model = np.array([i for i in range (0,20)])

mdl = np.zeros((len(points_model),len(points_model)))

for i in range(len(points_model)):
    for j in range(len(points_model)):
        
        mdl[i,j] = points_model[i]

Path1_a = [6,7,8,9,10,11,12,13,12,11,10,9,8,7] - np.ones(14)  ##Black path
Path1_b = [6,7,7,8,9,9,10,11,12,13,13,14,14,15] - np.ones(14) 

Path2_a = [6,8,9,11,12,13,14,13,14,13,12,9,8,7] - np.ones(14) ##Orange path
Path2_b = [6,5,3,3,4,6,9,11,11,15,17,16,16,15] - np.ones(14)

x = [i for i in range(20)]
y = [i for i in range(20)]
z1 = -np.log(z)

fig, ax1 = plt.subplots(figsize=(15,13))
ax1.contour(x, y, z1-np.min(z1), levels=14, linewidths=0.5, colors='k')
cntr1 = ax1.contourf(x, y, z1-np.min(z1), levels=14, cmap="RdBu_r")

for i in range(len(points_model)):

    plt.plot(points_model,mdl[i,:],'.', color = 'darkgreen')    


plt.plot(Path1_b,Path1_a,'-o' ,color = 'black', label = 'Path 1',linewidth=3)
plt.plot(Path2_b,Path2_a,'-o' ,color = 'orange', label = 'Path 2',linewidth=3)

plt.xlabel('x')
plt.grid()

plt.plot(x[5], mdl[5,5],'*',color='k')
plt.plot(x[14], mdl[6,14],'*',color='k')
plt.plot(x[10], mdl[12,10],'*',color='k')

#fig.colorbar(im, orientation="horizontal", pad=0.2)
fig.colorbar(cntr1, ax=ax1)
plt.legend(bbox_to_anchor=(0.5, 1.11), loc='upper center',ncol=3,borderaxespad=0., borderpad=2,prop={'size':15})
plt.savefig('3-wells-paths.png', bbox_inches='tight', dpi=500)
plt.close()
