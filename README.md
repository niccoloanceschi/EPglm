# EPGLM
R codes for efficient expectation propagation for generalized linear models under Gaussian prior.

Authors: Niccolò Anceschi, Augusto Fasano, Beatrice Franzolini, and Giovanni Rebaudo.

Overview

This repository is associated with the article Anceschi, Fasano, Franzolini, and Rebaudo (2024+) Scalable expectation propagation for generalized linear models. 

The key contribution of the paper is outlined below.
 
> [...] We [...] deriv[e] a novel efficient formulation of \ep for \textsc{glm}s, whose cost scales linearly in the number of covariates p, and reducing the state-of-the-art O(p<sup>2</sup>n) per-iteration computational cost of the \ep routine for GLMs to O(p n  min{p,n}), with n being the sample size.
We also show that, for binary models and log-linear GLMs approximate predictive means can be obtained at no additional cost.
To preserve efficient moment matching for count data, we propose employing a combination of log-normal Laplace transform approximations, avoiding numerical integration.
These novel results open the possibility of employing EP in settings that were believed to be practically impossible.

More precisely, we provide the R code to implement Algorithms 1 and 2 in Anceschi, Fasano, Franzolini, and Rebaudo (2024+) and replicate the results of their Illustration.

The repository contains the following:

1. `EPGLM_illustrationProbit.R` code to reproduce the results in Section 5.1;
2. `EPGLM_illustrationPoisson.R` code to reproduce the results in Section 5.2;
3. `EPGLM_fcts.R` functions needed to run the main code;
4. `Data-and-Results` folder with data and results of the analyses.


