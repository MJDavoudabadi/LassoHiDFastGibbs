pkgname <- "LassoHiDFastGibbs"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "LassoHiDFastGibbs-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('LassoHiDFastGibbs')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("blasso_gibbs_2block_bl")
### * blasso_gibbs_2block_bl

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: blasso_gibbs_2block_bl
### Title: Bayesian lasso Gibbs sampler: 2-block (beta–lambda2) variant
### Aliases: blasso_gibbs_2block_bl

### ** Examples

set.seed(1)
n <- 30; p <- 6
X <- matrix(rnorm(n * p), n, p)
y <- rnorm(n)
out <- blasso_gibbs_2block_bl(
  vy = y, mX = X,
  a = 1, b = 1, u = 1, v = 1,
  nsamples = 200,
  lambda_init = 1, sigma2_init = 1,
  verbose = 0
)
str(out)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("blasso_gibbs_2block_bl", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("blasso_gibbs_2block_bs")
### * blasso_gibbs_2block_bs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: blasso_gibbs_2block_bs
### Title: Bayesian lasso Gibbs sampler: 2-block (beta–sigma2) variant
### Aliases: blasso_gibbs_2block_bs

### ** Examples

set.seed(1)
n <- 30; p <- 6
X <- matrix(rnorm(n * p), n, p)
y <- rnorm(n)
out <- blasso_gibbs_2block_bs(
  vy = y, mX = X,
  a = 1, b = 1, u = 1, v = 1,
  nsamples = 200,
  lambda_init = 1, sigma2_init = 1,
  verbose = 0
)
str(out)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("blasso_gibbs_2block_bs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("blasso_pcg_lambda2_va")
### * blasso_pcg_lambda2_va

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: blasso_pcg_lambda2_va
### Title: Bayesian lasso PCG sampler: lambda2 collapsed over local scales
### Aliases: blasso_pcg_lambda2_va

### ** Examples

set.seed(1)
n <- 40; p <- 6
X <- matrix(rnorm(n * p), n, p)
beta <- c(1.2, 2, -1, 0.5, 0.75, 2.5)
y <- X %*% beta + rnorm(n)

out <- blasso_pcg_lambda2_va(
  vy = y, mX = X,
  a = 1, b = 1, u = 1, v = 1,
  nsamples = 200,
  lambda_init = 1, sigma2_init = 1,
  verbose = 0
)

summary(out$vlambda2)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("blasso_pcg_lambda2_va", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("blasso_pcg_sigma2_va")
### * blasso_pcg_sigma2_va

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: blasso_pcg_sigma2_va
### Title: Bayesian lasso PCG sampler: sigma2 collapsed over local scales
### Aliases: blasso_pcg_sigma2_va

### ** Examples

set.seed(1)
n <- 40; p <- 6
X <- matrix(rnorm(n * p), n, p)
beta <- c(1.2, 2, -1, 0.5, 0.75, 2.5)
y <- X %*% beta + rnorm(n)

out <- blasso_pcg_sigma2_va(
  vy = y, mX = X,
  a = 1, b = 1, u = 1, v = 1,
  nsamples = 200,
  lambda_init = 1, sigma2_init = 1,
  va_init = rep(1, p),
  verbose = 0
)

summary(out$vsigma2)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("blasso_pcg_sigma2_va", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("normalize")
### * normalize

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: normalize
### Title: Normalize Response and Covariates
### Aliases: normalize

### ** Examples

set.seed(1)
X <- matrix(rnorm(100 * 10), 100, 10)
beta <- c(2, -3, rep(0, 8))
y <- as.vector(X %*% beta + rnorm(100))
norm_result <- normalize(y, X)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("normalize", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("penalized_nested_Gibbs")
### * penalized_nested_Gibbs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: penalized_nested_Gibbs
### Title: Penalized nested Gibbs sampler for Bayesian linear regression
### Aliases: penalized_nested_Gibbs

### ** Examples

set.seed(1)
n <- 50; p <- 10
X <- matrix(rnorm(n * p), n, p)
y <- rnorm(n)

out <- penalized_nested_Gibbs(
  vy = y, mX = X,
  penalty_type = "lasso",
  a = 1, b = 1, u = 1, v = 1,
  nsamples = 200,
  lambda_init = 1,
  va_init = NULL,
  verbose = 0,
  lower = 1e-12,
  upper = 5000,
  s_beta = 1,
  s_siglam = 1
)
str(out)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("penalized_nested_Gibbs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("penalized_pcg_beta_sigma2")
### * penalized_pcg_beta_sigma2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: penalized_pcg_beta_sigma2
### Title: Penalized PCG sampler: beta block, lambda2 collapsed over sigma2
### Aliases: penalized_pcg_beta_sigma2

### ** Examples

set.seed(1)
n <- 40; p <- 6
X <- matrix(rnorm(n * p), n, p)
beta <- c(1.2, 2, -1, 0.5, 0.75, 2.5)
y <- X %*% beta + rnorm(n)

out <- penalized_pcg_beta_sigma2(
  vy = y, mX = X, penalty_type = "horseshoe",
  a = 1, b = 1, u = 1, v = 1,
  nsamples = 200,
  lambda_init = 1, sigma2_init = 1,
  verbose = 0
)

summary(out$mBeta)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("penalized_pcg_beta_sigma2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("penalized_pcg_lambda2_sigma2")
### * penalized_pcg_lambda2_sigma2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: penalized_pcg_lambda2_sigma2
### Title: Penalized PCG sampler: lambda2 collapsed over sigma2
### Aliases: penalized_pcg_lambda2_sigma2

### ** Examples

set.seed(1)
n <- 40; p <- 6
X <- matrix(rnorm(n * p), n, p)
beta <- c(1.2, 2, -1, 0.5, 0.75, 2.5)
y <- X %*% beta + rnorm(n)
out <- penalized_pcg_lambda2_sigma2(
  vy = y, mX = X, penalty_type = "lasso",
  a = 1, b = 1, u = 1, v = 1,
  nsamples = 200, lambda_init = 1, sigma2_init = 1,
  verbose = 0
)
str(out)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("penalized_pcg_lambda2_sigma2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("penalized_pcg_sigma2_beta")
### * penalized_pcg_sigma2_beta

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: penalized_pcg_sigma2_beta
### Title: Penalized PCG sampler: sigma2 collapsed over beta
### Aliases: penalized_pcg_sigma2_beta

### ** Examples

set.seed(1)
n <- 40; p <- 6
X <- matrix(rnorm(n * p), n, p)
beta <- c(1.2, 2, -1, 0.5, 0.75, 2.5)
y <- X %*% beta + rnorm(n)

out <- penalized_pcg_sigma2_beta(
  vy = y, mX = X, penalty_type = "horseshoe",
  a = 1, b = 1, u = 1, v = 1,
  nsamples = 200,
  lambda_init = 1, sigma2_init = 1,
  va_init = rep(1, p),
  verbose = 0
)

summary(out$vsigma2)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("penalized_pcg_sigma2_beta", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("penalized_pcg_sigma2_lambda2")
### * penalized_pcg_sigma2_lambda2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: penalized_pcg_sigma2_lambda2
### Title: Penalized PCG sampler: sigma2 collapsed over lambda2
### Aliases: penalized_pcg_sigma2_lambda2

### ** Examples

set.seed(1)
n <- 40; p <- 6
X <- matrix(rnorm(n * p), n, p)
beta <- c(1.2, 2, -1, 0.5, 0.75, 2.5)
y <- X %*% beta + rnorm(n)

out <- penalized_pcg_sigma2_lambda2(
  vy = y, mX = X, penalty_type = "lasso",
  a = 1, b = 1, u = 1, v = 1,
  nsamples = 200,
  lambda_init = 1, sigma2_init = 1,
  va_init = rep(1, p),
  verbose = 0
)

str(out)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("penalized_pcg_sigma2_lambda2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
