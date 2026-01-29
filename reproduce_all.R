pkg <- "LassoHiDFastGibbs"



dir.create("reproduction/helpers", recursive = TRUE, showWarnings = FALSE)

file.copy(
  from = c("inst/benchmarks/helpers/benchmark_3bg.R",
           "inst/benchmarks/helpers/benchmark_4bg.R",
           "inst/benchmarks/helpers/benchmark_blasso_bayeslm.R",
           "inst/benchmarks/helpers/benchmark_blasso_bayesreg.R",
           "inst/benchmarks/helpers/benchmark_blasso_hans.R",
           "inst/benchmarks/helpers/benchmark_blasso_monomvn.R",
           "inst/benchmarks/helpers/benchmark_blasso_rstan.R",
           "inst/benchmarks/helpers/benchmark_EHS.R",
           "inst/benchmarks/helpers/benchmark_horseshoe.R",
           "inst/benchmarks/helpers/benchmark_horseshoe_bayeslm.R",
           "inst/benchmarks/helpers/benchmark_horseshoe_bayesreg.R",
           "inst/benchmarks/helpers/benchmark_horseshoe_brms.R",
           "inst/benchmarks/helpers/generate_data.R",
           "inst/benchmarks/helpers/plot_densities.R",
           # these two are in R/
           "R/effective_sample_size.R",
           "R/mcmc_stats.R"),
  to = "reproduction/helpers",
  overwrite = TRUE
)



if (!requireNamespace(pkg, quietly = TRUE)) {
  stop("Please install LassoHiDFastGibbs first.")
}

bench_dir <- system.file("benchmarks", package = pkg)
if (bench_dir == "") stop("Benchmarks directory not found in installed package.")

out_dir <- file.path(getwd(), "reproduction")
dir.create(out_dir, showWarnings = FALSE)

# 1) copy all Rmd files
rmds <- list.files(bench_dir, pattern = "\\.Rmd$", full.names = TRUE)
file.copy(rmds, out_dir, overwrite = TRUE)

# 2) copy the helpers directory (recursively)
helpers_src <- file.path(bench_dir, "helpers")
helpers_dst <- file.path(out_dir, "helpers")
if (dir.exists(helpers_src)) {
  dir.create(helpers_dst, showWarnings = FALSE)
  helper_files <- list.files(helpers_src, full.names = TRUE, recursive = TRUE)
  # preserve subfolders
  rel_paths <- sub(paste0("^", gsub("\\\\", "/", helpers_src), "/?"), "", gsub("\\\\", "/", helper_files))
  dst_paths <- file.path(helpers_dst, rel_paths)
  dir.create(unique(dirname(dst_paths)), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(helper_files, dst_paths, overwrite = TRUE)
  if (!all(ok)) warning("Some helper files failed to copy.")
}

message("Benchmark materials copied to: ", normalizePath(out_dir))
message("Setwd('reproduction') then knit any Rmd in that folder.")
