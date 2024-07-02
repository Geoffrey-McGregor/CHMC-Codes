# CHMC-Codes

Welcome to our CHMC Matlab repository

You will find two folders containing the codes used to run the examples found in the CHMC paper

**Compact Chi** contains codes for sampling the p-Generalized Chi distribution 

**Compact PGauss** contains codes for sampling the p-Generalized Gaussian distribution.

# Compact Chi: 
Here you find four executable files associated with Figures from the paper: CompactConvergence, CompactHeatmapD, CompactHeatmapP and Compact Violin.

Running CompactViolin will reproduce Figure 1 from the paper
Running CompactHeatD will reproduce the top two panels of Figure 2 from the paper
Running CompactHeatP will reproduce the bottom two panels of Figure 2 from the paper
Running CompactConverence will reproduce Figure 3 from the paper

If you wish to sample the p-Generalized Chi distribution with your own set of parameters, call ChiSampler(d,p,N,Chains,T,dt), where

d and p are parameters associated with the distribution
N is the number of iterations
Chains is the number of chains
T is the integration time
dt is the time step

The output of ChiSampler is two matrices containing two Chains X N matrices, the first being the CHMC samples and the second being the HMC Leapfrog samples


# PCompact PGauss:
Here you find one executable file, PGaussSampler.

Running PGaussSampler with a specified dimension d on line 9 will reproduce a convergence plot as seen in Figure 4 of the paper.

You may choose the value of PGauss on line 5 to be 2,4 or 6. If you wish to choose a different value for p, you will need to reference our paper's supplemental for how to adjust the Newton solve.
