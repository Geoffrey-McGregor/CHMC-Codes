# CHMC-Codes

Welcome to our CHMC Matlab repository

You will find two directories containing the codes used to run the examples found in our CHMC paper

1. `Compact Chi contains codes for sampling the p-Generalized Chi distribution.` 

2. `Compact PGauss contains codes for sampling the p-Generalized Gaussian distribution.`

3. `Compact PNorm contains codes for sampling the p-Generalized Gaussian distribution then transforming the samples to the p-generalized Chi distribution.`

# Compact Chi 
**Executible files to reproduce figures** _(No modification to the code or parameters is required)_

1. `CompactHistogramViolinConvergence.m produces Figure 1.`

2. `CompactHeat.m produces the two panels of Figure 2.`

- These files calls **ChiSampler.m** with specific parameter values, and runs the HMC and CHMC algorithms to sample from the p-generalized Chi distribution.


**Additional executible file**

1. `ChiSampler(d,p,N,Chains,T,dt) samples the p-Generalized Chi distribution with your choice of parameters.`  
 - **Parameters:** See preamble in .m file for parameter descriptions.
 - **Purpose:** ChiSampler.m runs the HMC and CHMC algorithms to sample from the p-generalized Chi distribution.
 - **Output:** ChiSampler.m produces two matrices, the first being the CHMC samples and the second being the HMC-Leapfrog samples. The number of rows corresponds to the number of chains and the number of columns corresponds to the number of samples used.

**Additional .m files**
1. `AGSamChi.m is an exact sampler built for sampling the p-generalized Chi distribution. We use this to compute the CHMC and HMC sampling errors.` 

2. `Violin.m and Violinplot.m are called by CompactViolin.m to genereate the violin plots from Figure 1.`

3. `ws_distance.m is used to compute the Wasserstein 1 distance between the discrete samples generated by ChiSampler.m and AGSamChi.m`

# Compact PGauss:
**Executible files to reproduce figures** _(No modification to the code or parameters is required)_

1.`CompactPGaussComparison.m produces Figure 3.`

**Additional .m files**

1. `PGaussSampler.m is called by **CompactPGaussComparison.m** to generate the convergence data for the specified parameters.`
  - PGaussSampler.m calls **HMCSolver.m** to sample using the Leapfrog integrator and calls **CHMCVectorSolver.m** to sample using the conservative integrator.
  - dimArray contains the dimensions **d** and the value of **p** corresponds to the p-generalized Gaussian distributions being sampled.
  - **N** is the number of samples and **Chains** is the number of chains,**T** **dt** are the integration time and step size respectively.
  - **numPoints** specifies the number of iterations between each metric computation.
2. `AGSam.m is an exact sampler built for sampling one-dimensional slices of the p-generalized Gaussian distribution. We use this to compute the CHMC and HMC sampling errors.`
3. `CHMCSolver.m obtains samples using the conservative integrator.`
4. `HMCSolver.m obtains samples using the Leapfrog integrator.`
5. `PGenCDF.m evalutes the CDF of a one-dimensional slice of the p-generalized Gaussian Distribution. We use this for computing CHMC and HMC errors.`
6. `VarComp.m is used for compute errors in the covariance for CHMC and HMC samplers.`


# Compact PNorm
**Executible files to reproduce figures** _(No modification to the code or parameters is required)_

1. `CompactPNormComparison.m produces Figure 4.`
- Calls **PGaussSampler.m** using **d**=2560, 10240 and 40960 and **dt**=0.1, 0.05 and 0.025 to obtain samples from the p-generalized Gaussian.

**Additional .m files**

1. `CHMCVectorSolverJ0.m obtains samples using the conservative integrator with approximate Jacobian 1.`
2. `CHMCVectorSolverFullJ.m obtains samples using the conservative integrator using the exact Jacobian.`
3. `HMCSolver.m obtains samples using the Leapfrog integrator.`
4. `Violin.m and Violinplot.m are called by CompactPNormComparison.m to genereate the violin plots from Figure 5.`



# Matlab requirements

Executing these files requires MATLAB version R2023b or later and the [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html)
