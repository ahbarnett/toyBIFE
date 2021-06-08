# toyBIFE
simple experiments in Bayesian inference of free energy and path selection for cryo-EM, for collaboration with Pilar Cossio and group.

Alex Barnett, June 2021.

See ``notes`` and ``code``.

Added a new folder **2DHsp90** with 3 subfolders:

## FES: Data for vizualizing the 2D distribution or free energy surface.
Files in folder:
 * Num\_Images\_grid: number of images at each x,y point (in matrix form 20 by 20). This is the 'scatter plot'.
 * 2D-FES\_plot.py: python code to graph FES using Num\_Images\_grid and the two reference paths (black and orange). 
 Run: python3.6 2D-FES\_plot.py 
 Output: 3-wells-paths.png

## LikeLihood: Data for calculating the log-likelihood given a predetermined path and G that maximizes it.
Files in folder:
* likelihood\_calc.py: calculates the log-L (total and for each image) given a path and Gmax. 
Run: python3.6 likelihood\_calc.py [FILE\_PATH] [FILE\_GMAX] > Log-L-image
e.g., FILE\_PATH=Path-black and FILE\_GMAX=Gmax-black
Output: log-L of individual images and final line total log-L. 
* Path-black and Path-orange: Node indexes (following Julian's notation) of the black and orange paths
* Gmax-black and Gmax-orange: Gs that maximize cryoBIFE likelihood. 
* all\_images: LogP BioEM for all images (rows) and the 225 models (columns) within the green square (slide 8 alexdiscussion.pdf; typo there not 200 but 225)

## cryoBIFE: Data for running cryoBIFE given a predetermined path. 
Files in folder:
* cryoBIFE.py: MCMC from paper using Alex's logsumexp function. 
Run:  python3.6 cryoBIFE.py [FILE\_PATH] > LogL-[FILE\_PATH]
Output: LogL-[FILE\_PATH] logL at each MCMC step, an additional output file "Gs" has the Gs at each MCMC step.
* get-max.sh: bash script to get Gs that maximize Log-L. 
Run: ./get-max.sh LogL-[FILE\_PATH] 
* Check-logsumexp: Checking the results from using Alex's routine and the old (correct) one that we had (see commented-out line 148 - 160)
* all\_images: same as above.
* Path-black and Path-orange: same as above.
