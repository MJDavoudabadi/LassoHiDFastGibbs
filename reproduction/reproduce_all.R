#!/usr/bin/env Rscript

# JSS replication entry point
# From the repo root, run:
#   Rscript reproduction/reproduce_all.R
#
# This runs the 4 benchmark Rmd files sequentially, then runs create_tables.Rmd.
# It assumes all Rmds are in reproduction/ and helper scripts are in reproduction/helpers/.

options(stringsAsFactors = FALSE)

pkg <- "LassoHiDFastGibbs"

# --------- locate repo root from THIS script's location (no hard-coded paths) ---------
get_script_path <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", ca, value = TRUE)
  if (length(file_arg) == 1) return(normalizePath(sub("^--file=", "", file_arg), winslash = "/"))
  # If someone runs via source() in RStudio, we fall back to current wd
  return(normalizePath(getwd(), winslash = "/"))
}

find_repo_root <- function(start_dir) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in 1:60) {
    if (file.exists(file.path(cur, "DESCRIPTION"))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  stop("Could not find repo root (DESCRIPTION not found). Start: ", start_dir)
}

script_path <- get_script_path()
script_dir  <- dirname(script_path)

# script is in <root>/reproduction, so root is its parent
root <- find_repo_root(script_dir)

repro_dir <- file.path(root, "reproduction")
helpers_dir <- file.path(repro_dir, "helpers")
out_dir <- file.path(repro_dir, "output")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Repo root: ", root)
message("Reproduction dir: ", repro_dir)
message("Helpers dir: ", helpers_dir)
message("Output dir: ", out_dir)

# --------- install deps + (re)install package from this repo ---------
install_if_missing <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    message("Installing missing packages: ", paste(miss, collapse = ", "))
    install.packages(miss, repos = "https://cloud.r-project.org")
  }
}

install_if_missing(c("remotes", "rmarkdown", "knitr"))

# Ensure compiled code matches this checkout
message("Installing ", pkg, " from local repo ...")
remotes::install_local(root, dependencies = TRUE, upgrade = "never", force = TRUE)

suppressPackageStartupMessages(library(pkg, character.only = TRUE))

# --------- Rmd paths (exactly your 4 scenarios + tables) ---------
bench_rmds <- c(
  file.path(repro_dir, "benchmark_lasso_ngtp.Rmd"),
  file.path(repro_dir, "benchmark_lasso_pgtn.Rmd"),
  file.path(repro_dir, "benchmark_horseshoe_ngtp.Rmd"),
  file.path(repro_dir, "benchmark_horseshoe_pgtn.Rmd")
)
tables_rmd <- file.path(repro_dir, "create_tables.Rmd")

missing <- c(bench_rmds, tables_rmd)[!file.exists(c(bench_rmds, tables_rmd))]
if (length(missing)) stop("Missing required file(s):\n- ", paste(missing, collapse = "\n- "))

if (!dir.exists(helpers_dir)) stop("helpers folder not found: ", helpers_dir)

# --------- IMPORTANT: render in reproduction/ so relative 'helpers/...' works ---------
# Many of your Rmds use relative paths like 'helpers/...', so we set knit_root_dir.
knitr::opts_knit$set(root.dir = repro_dir)

# One shared environment so create_tables sees objects created by benchmark runs
env <- new.env(parent = globalenv())

set.seed(1)

# Patch eval=interactive() WITHOUT moving the file (write a temp patched file in reproduction/)
patch_in_place <- function(rmd_path) {
  patched <- file.path(dirname(rmd_path), paste0(".patched_", basename(rmd_path)))
  txt <- readLines(rmd_path, warn = FALSE)
  txt <- gsub("eval\\s*=\\s*interactive\\(\\)", "eval=TRUE", txt)
  writeLines(txt, patched)
  patched
}

render_one <- function(rmd_path) {
  patched <- patch_in_place(rmd_path)
  on.exit(unlink(patched), add = TRUE)

  message("\n=== Rendering: ", basename(rmd_path), " ===")
  rmarkdown::render(
    input      = patched,
    output_dir = out_dir,
    envir      = env,
    quiet      = FALSE
  )
}

# --------- run 4 benchmark Rmds sequentially ---------
for (f in bench_rmds) render_one(f)

# --------- then build paper tables ---------
render_one(tables_rmd)

# Save session info
sink(file.path(out_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("\nDone. Outputs are in: ", out_dir)
