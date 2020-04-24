# SoS_LIP_Model

**Purpose:** This repository reproduces the balance and step recovery capturability results using the Sums-of-Squares (SOS) Optimization procedure described in the follwing paper: 
* *Balancing and Step Recovery Capturability via Sums-of-Squares Optimization*. Posa, M. Koolen, T. Tedrake, R. (2017)
* The paper can be found [here](http://rss2017.lids.mit.edu/static/papers/09.pdf).

## Outer Approximation
To solve the outer approximation of the 0-step viable-capture region for the LIP Model run [Code/LIPM_Outer_0step.m](https://github.com/grantgib/SoS_LIP_Model/blob/master/Code/LIPM_Outer_0step.m). To solve the outer approximation of the 1-step viable-capture region run [Code/LIPM_Outer_1step.m](https://github.com/grantgib/SoS_LIP_Model/blob/master/Code/LIPM_Outer_1step.m). Note that the solution of the 0-step region is used to compute the 1-step solution.
