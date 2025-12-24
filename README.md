# Blocked, Partially Collapsed, and Nested Gibbs Sampling for the Bayesian Lasso

This repository contains code accompanying the paper **“Blocked, Partially Collapsed, and Nested Gibbs Sampling for the Bayesian Lasso.”**  

In this work, we develop Gibbs samplers using blocking, partial collapsing, and nesting strategies to improve sampling efficiency. These strategies lead to more efficient exploration of the posterior distribution, particularly improving the mixing of the residual variance and penalty tuning parameters, which are often slow to converge in standard implementations. We show that partial collapsing typically results in moderate improvements, while nested Gibbs samplers offer significant gains over traditional two-block Gibbs samplers. The scalability of our approach allows fitting high-dimensional Bayesian Lasso models with tens of thousands of predictors in seconds to minutes. We further extend the methodology to the horseshoe penalty and demonstrate that similar conclusions hold. 

The purpose of this repository is to **reproduce the results presented in the research paper**.

---

## 📂 Repository Structure

- `cpp/` – C++ scripts for:
  - Partially Collapsed (PC), Hans, and slice sampling algorithms  
  - Functions for generating inverse Gaussian random numbers
  - Code for sampling from the Lasso distribution

- `R/` – R scripts implementing:
  - PC, Hans, two-block, three-block, and four-block Gibbs sampling algorithms  
  - Functions for sampling from `bayeslm`, `bayesreg`, `monomvn`, and `rstan` packages  
  - Functions for generating and normalizing data  
  - Effective sample size computation and MCMC diagnostics  
  - Horseshoe penalty and nested algorithms  

- `data/` – Datasets used in the paper  

---

## 📑 R Files

- **benchmark_horseshoe_ngtp.Rmd**  
  Benchmarks Bayesian Lasso when *n > p* using the horseshoe penalty. Runs the modified Park and Casella Gibbs sampler, nested Gibbs, two-block Gibbs, partially collapsed Gibbs, `bayeslm`, `bayesreg`, `monomvn`, and `rstan` algorithms.  

- **benchmark_horseshoe_pgtn.Rmd**  
  Benchmarks Bayesian Lasso when *p > n* using the horseshoe penalty. Includes Park and Casella Gibbs sampler, benchmarking functions for two-block, nested, and partially collapsed Gibbs samplers, and exact/approximate horseshoe (EHS/AHS).  

- **benchmark_lasso_ngtp.Rmd**  
  Benchmarks Bayesian Lasso when *n > p* using the Lasso penalty. Includes Park and Casella, nested Gibbs, two-block Gibbs, modified Hans Gibbs, partially collapsed Gibbs, `bayeslm`, `bayesreg`, `monomvn`, and `rstan`.  

- **benchmark_lasso_pgtn.Rmd**  
  Benchmarks Bayesian Lasso when *p > n* using the Lasso penalty. Similar to above, with additional algorithms included.  

- **create_tables.Rmd**  
  Generates the tables provided in the paper.  

- **benchmark_blasso_bayeslm.R**  
  Samples from `bayeslm` (Lasso penalty).  

- **benchmark_blasso_bayesreg.R**  
  Samples from `bayesreg` (Lasso penalty).  

- **benchmark_blasso_hans.R**  
  Samples from the modified Hans Gibbs sampler.  

- **benchmark_blasso_monomvn.R**  
  Samples from `monomvn` (Lasso penalty).  

- **benchmark_blasso_rstan.R**  
  Samples from `rstan`.  

- **benchmark_horseshoe_bayeslm.R**  
  Samples from `bayeslm` (horseshoe penalty).  

- **benchmark_horseshoe_bayesreg.R**  
  Samples from `bayesreg` (horseshoe penalty).  

- **benchmark_horseshoe_brms.R**  
  Samples from `brm` in the `brms` package (Stan backend) with the horseshoe penalty.  

- **benchmark_nested_gibbs.R**  
  Implements the nested Gibbs sampler.  

- **benchmark_2bg.R**  
  Two-block Gibbs sampler:  
  - `benchmark_blasso_2BG_beta_sigma2`: β and σ² in one block  
  - `benchmark_blasso_2BG_beta_lambda2`: β and λ² in one block  

- **benchmark_3bg.R**  
  Three-block Gibbs sampler.  

- **benchmark_4bg.R**  
  Four-block Gibbs sampler.  

- **benchmark_pcg.R**  
  Partially Collapsed Gibbs samplers, where some conditional parameters are marginalized from the parent Gibbs sampler.  

