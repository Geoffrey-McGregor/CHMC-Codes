# CHMC-Codes

Welcome to our CHMC Matlab repository

You will find two directories containing the codes used to run the examples found in our CHMC paper

1. `Compact Chi contains codes for sampling the p-Generalized Chi distribution.` 

2. `Compact PGauss contains codes for sampling the p-Generalized Gaussian distribution.`

# Compact Chi 
**Executible files to reproduce figures** _(No modification to the code or parameters is required)_

1. `CompactViolin_ produces Figure 1`

2. `CompactHeatD produces the top two panels of Figure 2`

3. `CompactHeatP_ will produces the bottom two panels of Figure 2`

4. `CompactConverence_ produces Figure 3`


**Additional executible file**

1. `ChiSampler(d,p,N,Chains,T,dt) samples the p-Generalized Chi distribution with your choice of parameters`  
 - See preamble in .m file for parameter descriptions
 - The output of ChiSampler is two matrices, the first being the CHMC samples and the second being the HMC-Leapfrog samples. The number of rows corresponds to the number of chains and the number of columns corresponds to the number of samples used.


# PCompact PGauss:
**Executible files**

1. `PGaussSampler produces one of the convergence plots from Figure 4.`
- Choosing the dimension **d** on line 9 from 1024, 2560 ,5120, 10240, 20480 or 40960 will reproduce the corresponding panels from Figure 4.
- You may choose the value of **PGauss** on line 5 to be 2,4 or 6. Other choices will require alterations of the solvers in _CHMCVectorSolver.m_ and _HMCSolver.m_ See the supplemental for implementation details.

# Matlab requirements

Executing these files requires MATLAB version R2023b or later and the [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html)