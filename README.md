# toyBIFE
simple experiments in Bayesian inference of free energy and path selection for cryo-EM, for collaboration with Pilar Cossio and group.

Alex Barnett, June 2021.

See ``notes`` and ``code``.

Added a new folder 2DHsp90 with 3 subfolders:

* Item 1: FES: Data for vizualizing the 2D distribution or free energy surface.
 * Sub Item 1: Num\_Images\_grid: number of images at each x,y point (in matrix form 20 by 20). This is the 'scatter plot'.
 * Sub Item 2: 2D-FES\_plot.py: python code to graph FES using Num\_Images\_grid and the two reference paths (black and orange), output: 3-wells-paths.png

* LikeLihood: Data for calculating the log-likelihood given a predetermined path and G that maximizes the cryoBIFE likelilood.
** likelihood\_calc.py: calculates the total log-L and for each image, given a path and Gmax. Run: python3.6 likelihood\_calc.py [FILE\_PATH] [FILE\_GMAX] > Log-L-image (individual log-L and final line total log-L). e.g., FILE\_PATH=Path-black and FILE\_GMAX=Gmax-black
** Path-black and Path-orange: Node indexes (following Julian's notation) of the black and orange paths
** Gmax-black and Gmax-orange: Gs that maximize cryoBIFE likelihood. 
** all\_images: LogP BioEM for all images (rows) and the 225 models (columns) within the green square (slide 8 alexdiscussion.pdf; typo there not 200 but 225)

* cryoBIFE: Data for running cryoBIFE given a predetermined path. 
** cryoBIFE.py: MCMC from paper using Alex's logsumexp function. Run:  python3.6 cryoBIFE.py [FILE\_PATH] > LogL-[FILE\_PATH]; an additional output file Gs (has the Gs for each MCMC step)
** get-max.sh: bash script to get Gs that maximize Log-L. Run: ./get-max.sh LogL-[FILE\_PATH] 
** all\_images: LogP BioEM for all images (rows) and the 225 models (columns) within the green square (slide 8 alexdiscussion.pdf; typo there not 200 but 225)
** Path-black and Path-orange: Node indexes (following Julian's notation) of the black and orange paths
