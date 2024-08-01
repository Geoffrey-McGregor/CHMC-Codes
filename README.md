# CHMC-Codes

Welcome to our CHMC Matlab repository

You will find two directories containing the codes used to run the examples found in our CHMC paper

1. `Compact Chi contains codes for sampling the p-Generalized Chi distribution.` 

2. `Compact PGauss contains codes for sampling the p-Generalized Gaussian distribution.`

# Compact Chi 
**Executible files to reproduce figures** _(No modification to the code or parameters is required)_

1. `CompactViolin.m produces Figure 1.`

2. `CompactHeatD.m produces the top two panels of Figure 2.`

3. `CompactHeatP.m will produces the bottom two panels of Figure 2.`

4. `CompactConvergence.m produces Figure 3.`
- Each of these four files calls ChiSampler.m with specific parameter values. ChiSampler.m runs the HMC and CHMC algorithms to sample from the p-generalized Chi distribution.


**Additional executible file**

1. `ChiSampler(d,p,N,Chains,T,dt) samples the p-Generalized Chi distribution with your choice of parameters.`  
 - **Parameters:** See preamble in .m file for parameter descriptions.
 - **Purpose:** ChiSampler.m runs the HMC and CHMC algorithms to sample from the p-generalized Chi distribution.
 - **Output:** ChiSampler.m produces two matrices, the first being the CHMC samples and the second being the HMC-Leapfrog samples. The number of rows corresponds to the number of chains and the number of columns corresponds to the number of samples used.

**Additional .m files**
1. `AGSamChi.m is an exact sampler built for sampling the p-generalized Chi distribution. We use this to compute the CHMC and HMC sampling errors.` 

2. `HeatmapGen.m is called by CompactHeatD.m and CompactHeatP.m to generate the heatmap from Figure 2.`

3. `Violin.m and Violinplot.m are called by CompactViolin.m to genereate the violin plots from Figure 1.`

4. `ws_distance.m is used to compute the Wasserstein 1 distance between the discrete samples generated by ChiSampler.m and AGSamChi.m`

# Compact PGauss:
**Executible files**

1. `PGaussSampler.m produces one of the convergence plots from Figure 4.`
- Choosing the dimension **d** on line 9 from 1024, 2560 ,5120, 10240, 20480 or 40960 will reproduce the corresponding panels from Figure 4.
- You may choose the value of **PGauss** on line 5 to be 2,4 or 6. Other choices will require alterations of the solvers in _CHMCVectorSolver.m_ and _HMCSolver.m_ See the paper's supplemental materials for implementation details.
- PGaussSampler.m calls **HMCSolver.m** to sample using the Leapfrog integrator and calls **CHMCVectorSolver.m** to sample using the Conservative integrator.

**Additional .m files**

1. `AGSam.m is an exact sampler built for sampling one-dimensional slices of the p-generalized Chi distribution. We use this to compute the CHMC and HMC sampling errors.`
2. `CHMCVectorSolver.m obtains samples using the conservative integrator.`
3. `HMCSolver.m obtains samples using the Leapfrog integrator.`
4. `PGenCDF.m evalutes the CDF of a one-dimensional slice of the p-generalized Gaussian Distribution. We use this for computing CHMC and HMC errors.`
5. `VarComp.m is used for compute errors in the covariance for CHMC and HMC samplers.` 


# Matlab requirements

Executing these files requires MATLAB version R2023b or later and the [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html)