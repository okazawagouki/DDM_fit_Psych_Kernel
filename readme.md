# Fitting drift diffusion model to behavioral data (shared version)

---

## Overview

This code fits the drift diffusion model (DDM) to choice and reaction time data using maximum likelihood estimation. On each iteration of parameter search, it solves the Fokker-Planck equation with the finite difference method to compute a joint probability distribution of choice and reaction time. The code also includes the comparison of psychophysical kernels between model and data.
 

## Requirements
The code is tested using MATLAB 2019b (with Statistics Toolbox).  
It uses a mex file solving the Fokker-Planck equation developed in Kiani lab: (https://github.com/KianiLab/FokkerPlanck)[https://github.com/KianiLab/FokkerPlanck]. The mex file is included in this repository.


## Package

[tutorial\_DDM.m](./tutorial_DDM.m) - code to run DDM fit.

[tutorial\_DDM\_PsychKernel.m](./tutorial_DDM_PsychKernel.m) - code to run DDM fit with fluctuation in stimulus strength. With the fluctuation, psychophysical kernel can be calculated.

`example_data.mat` - example behavioral data of one subject from Okazawa (2018).

`example_PK_data.mat` - example behavioral data from Luo (2025).

## Important Notes

1. The fittings in the tutorial could take more than half an hour, although they use a speeded option (it runs fits with relaxed convergence criteria). If you turn this option off, the fits would take even longer.
2. The package contains example data and example results of model fits. If you run the tutorial, you will see the results of the model fits.
3. This tutorial runs only one fit using one initial starting point, but ideally you want to use ten or more starting points and choose the one that converges with the highest log likelihood.


## Reference

**Fokker-Planck equation**:

Karlin S. & Taylor H. E. A (1981) *Second Course in Stochastic Processes* (Elsevier, Amsterdam).

**DDM fit**:

Kiani R., and Shadlen M.N. (2009). Representation of confidence associated with a decision by neurons in the parietal cortex. *Science* 324, 759–764.


**Psychophysical kernel**:

Okazawa G., Sha L., Purcell B.A., and Kiani R. (2018). Psychophysical reverse correlation reflects both sensory and decision-making processes. *Nat. Commun.* 9, 3479.

Luo T., Zheng Z., Xu M., Okazawa G. (2025).  Limitation of switching sensory information flow in flexible perceptual decision making.  *Nat. Commun.*   16: 172


---
Okazawa lab @ ION

[http://english.cebsit.cas.cn/lab/okazawagoki/research/](http://english.cebsit.cas.cn/lab/okazawagoki/research/)


[okazawa@ion.ac.cn](mailto:okazawa@ion.ac.cn)

Package developed by Luo Tianlin and Gouki Okazawa.




