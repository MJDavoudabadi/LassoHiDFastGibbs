pkg <- "LassoHiDFastGibbs"

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
