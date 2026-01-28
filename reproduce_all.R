pkg <- "LassoHiDFastGibbs"

if (!requireNamespace(pkg, quietly = TRUE)) {
  stop("Please install LassoHiDFastGibbs first.")
}

bench_dir <- system.file("benchmarks", package = pkg)
if (bench_dir == "") stop("Benchmarks directory not found in installed package.")

out_dir <- file.path(getwd(), "reproduction")
dir.create(out_dir, showWarnings = FALSE)

files <- list.files(bench_dir, pattern = "\\.Rmd$", full.names = TRUE)
file.copy(files, out_dir, overwrite = TRUE)

message("Benchmark Rmd files copied to: ", normalizePath(out_dir))
message("Open any Rmd in this folder and knit it to reproduce results.")
