# CHMC-Codes

Welcome to our CHMC Matlab repository

You will find two directories containing the codes used to run the examples found in our CHMC paper

1. `Compact Chi contains codes for sampling the p-Generalized Chi distribution.` 

2. `Compact PGauss contains codes for sampling the p-Generalized Gaussian distribution.`

# Compact Chi 
**Executible files to reproduce figures** _(No modification to the code or parameters is required)_

1. `CompactViolin_ produces Figure 1`

2. `CompactHeatD produces the top two panels of Figure 2`

3. `CompactHeatP_ will produces the bottom two panels of **Figure 2`

4. `CompactConverence_ produces Figure 3`


If you wish to sample the p-Generalized Chi distribution with your own set of parameters, call _ChiSampler(d,p,N,Chains,T,dt)_, where

**d** and **p** are parameters associated with the distribution

**N** is the number of samples in each chain

**Chains** is the number of chains

**T** is the integration time

**dt** is the time step

The output of ChiSampler is two matrices containing two Chains by N matrices, the first being the CHMC samples and the second being the HMC-Leapfrog samples


# PCompact PGauss:
Here you find one executable file, PGaussSampler.

Running _PGaussSampler_ with a specified dimension d (1024, 2560 ,5120, 10240, 20480 or 40960) on line 9 will reproduce a convergence plot as seen in **Figure 4** of the paper.

You may choose the value of PGauss on line 5 to be 2,4 or 6. If you wish to choose a different value for p, you will need to reference our paper's supplemental for how to adjust the Newton solve.
