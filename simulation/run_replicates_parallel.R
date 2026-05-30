#!/usr/bin/env Rscript

allowed_settings <- c(
  "Setting_I.r",
  "Setting_II_df=3.r",
  "Setting_II_df=5.r",
  "Setting_III.r",
  "Setting_IV.r",
  "Setting_SUPP.R"
)

usage <- paste(
  "Usage:",
  "Rscript simulation/run_replicates_parallel.R --setting Setting_I.r --replicates 20 --jobs 4 --bootstrap-runs 10 --run-bootstrap --max-attempts 5 --bootstrap-max-attempts 3",
  "",
  "Options:",
  "  --setting <file>                  One of Setting_I.r, Setting_II_df=3.r, Setting_II_df=5.r, Setting_III.r, Setting_IV.r, Setting_SUPP.R",
  "  --replicates <int>                Number of independent simulation replicates",
  "  --jobs <int>                      Number of concurrent Rscript workers",
  "  --bootstrap-runs <int>            Number of bootstrap runs per replicate",
  "  --run-bootstrap                  Run bootstrap; omit to skip bootstrap",
  "  --max-attempts <int>              Main simulation convergence retry limit",
  "  --bootstrap-max-attempts <int>    Bootstrap convergence retry limit",
  "  --bias-cutoff <number>            rBias removal cutoff used in summary code",
  "  --base-seed <int>                 Seed for first replicate; replicate i uses base_seed + i - 1",
  "  --output-root <path>              Output directory; defaults to simulation_output_parallel/<setting>_<timestamp>",
  sep="\n"
)

parse_args <- function(args){
  opts <- list(
    setting="Setting_I.r",
    replicates=10L,
    jobs=4L,
    bootstrap_runs=0L,
    run_bootstrap=FALSE,
    max_attempts=5L,
    bootstrap_max_attempts=5L,
    bias_cutoff=1000,
    base_seed=1L,
    output_root=NULL
  )

  i <- 1L
  while(i <= length(args)){
    arg <- args[[i]]
    if(arg == "--run-bootstrap"){
      opts$run_bootstrap <- TRUE
      i <- i + 1L
      next
    }

    needs_value <- c(
      "--setting", "--replicates", "--jobs", "--bootstrap-runs",
      "--max-attempts", "--bootstrap-max-attempts", "--bias-cutoff",
      "--base-seed", "--output-root"
    )
    if(!(arg %in% needs_value)){
      stop("Unknown argument: ", arg, "\n\n", usage)
    }
    if(i == length(args)){
      stop("Missing value for: ", arg, "\n\n", usage)
    }
    value <- args[[i + 1L]]
    if(arg == "--setting") opts$setting <- value
    if(arg == "--replicates") opts$replicates <- as.integer(value)
    if(arg == "--jobs") opts$jobs <- as.integer(value)
    if(arg == "--bootstrap-runs") opts$bootstrap_runs <- as.integer(value)
    if(arg == "--max-attempts") opts$max_attempts <- as.integer(value)
    if(arg == "--bootstrap-max-attempts") opts$bootstrap_max_attempts <- as.integer(value)
    if(arg == "--bias-cutoff") opts$bias_cutoff <- as.numeric(value)
    if(arg == "--base-seed") opts$base_seed <- as.integer(value)
    if(arg == "--output-root") opts$output_root <- value
    i <- i + 2L
  }
  opts
}

`%||%` <- function(x, y){
  if(is.null(x)) y else x
}

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork=FALSE), error=function(e) NA_character_)
if(is.na(script_path)){
  cmd <- commandArgs(trailingOnly=FALSE)
  file_arg <- grep("^--file=", cmd, value=TRUE)
  if(length(file_arg) > 0){
    script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork=FALSE)
  } else {
    script_path <- file.path(getwd(), "simulation", "run_replicates_parallel.R")
  }
}
script_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(script_dir, ".."), mustWork=TRUE)

opts <- parse_args(commandArgs(trailingOnly=TRUE))

if(!(opts$setting %in% allowed_settings)){
  stop("Invalid --setting: ", opts$setting, "\n\n", usage)
}
if(!is.finite(opts$replicates) || opts$replicates < 0) stop("--replicates must be >= 0")
if(!is.finite(opts$jobs) || opts$jobs < 1) stop("--jobs must be >= 1")
if(!is.finite(opts$bootstrap_runs) || opts$bootstrap_runs < 0) stop("--bootstrap-runs must be >= 0")

if(is.null(opts$output_root)){
  setting_name <- tools::file_path_sans_ext(opts$setting)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  opts$output_root <- file.path(repo_root, "simulation_output_parallel", paste0(setting_name, "_", timestamp))
}
dir.create(opts$output_root, recursive=TRUE, showWarnings=FALSE)

setting_path <- file.path(repo_root, "simulation", opts$setting)
if(!file.exists(setting_path)) stop("Setting file not found: ", setting_path)

run_one <- function(i, opts, repo_root, setting_path){
  rep_output <- file.path(opts$output_root, sprintf("rep_%03d", i))
  dir.create(rep_output, recursive=TRUE, showWarnings=FALSE)
  seed <- opts$base_seed + i - 1L
  setting_name <- tools::file_path_sans_ext(opts$setting)
  log_file <- file.path(rep_output, paste0(setting_name, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))

  env <- c(
    SIM_REP="1",
    SIM_SEED=as.character(seed),
    SIM_K_RUNS=as.character(opts$bootstrap_runs),
    SIM_SKIP_BOOTSTRAP=if(opts$run_bootstrap) "false" else "true",
    SIM_MAX_ATTEMPTS=as.character(opts$max_attempts),
    SIM_BOOTSTRAP_MAX_ATTEMPTS=as.character(opts$bootstrap_max_attempts),
    SIM_RBIAS_REMOVE_CUTOFF=as.character(opts$bias_cutoff),
    SIM_OUTPUT_DIR=normalizePath(rep_output, mustWork=FALSE)
  )

  header <- c(
    paste0("Repository: ", repo_root),
    paste0("Setting: ", opts$setting),
    paste0("Replicate: ", i),
    paste0("SIM_SEED=", seed),
    paste0("SIM_K_RUNS=", opts$bootstrap_runs),
    paste0("SIM_SKIP_BOOTSTRAP=", env[["SIM_SKIP_BOOTSTRAP"]]),
    paste0("SIM_MAX_ATTEMPTS=", opts$max_attempts),
    paste0("SIM_BOOTSTRAP_MAX_ATTEMPTS=", opts$bootstrap_max_attempts),
    paste0("SIM_RBIAS_REMOVE_CUTOFF=", opts$bias_cutoff),
    paste0("SIM_OUTPUT_DIR=", env[["SIM_OUTPUT_DIR"]]),
    paste0("LogFile=", log_file),
    ""
  )
  writeLines(header, log_file, useBytes=TRUE)

  message("Starting replicate ", i, " seed=", seed, " output=", rep_output)
  start_time <- proc.time()[["elapsed"]]
  old_wd <- getwd()
  on.exit(setwd(old_wd), add=TRUE)
  setwd(repo_root)

  env_strings <- paste(names(env), env, sep="=")
  output <- system2("Rscript", c(setting_path), env=env_strings, stdout=TRUE, stderr=TRUE)
  status <- attr(output, "status")
  if(is.null(status)) status <- 0L
  write(output, file=log_file, append=TRUE)
  elapsed <- proc.time()[["elapsed"]] - start_time
  write(paste0("ElapsedSeconds=", round(elapsed, 1)), file=log_file, append=TRUE)

  if(status != 0L){
    stop("Replicate ", i, " failed with exit code ", status, ". Log: ", log_file)
  }

  list(replicate=i, seed=seed, output=rep_output, log=log_file, elapsed=elapsed)
}

cat("Repository:", repo_root, "\n")
cat("Setting:", opts$setting, "\n")
cat("Replicates:", opts$replicates, "\n")
cat("Jobs:", opts$jobs, "\n")
cat("BootstrapRuns:", opts$bootstrap_runs, "\n")
cat("RunBootstrap:", opts$run_bootstrap, "\n")
cat("OutputRoot:", normalizePath(opts$output_root, mustWork=FALSE), "\n")

if(opts$replicates == 0L){
  cat("Parallel replicates complete. OutputRoot=", normalizePath(opts$output_root, mustWork=FALSE), "\n", sep="")
  quit(status=0)
}

jobs <- min(opts$jobs, opts$replicates)
cl <- parallel::makeCluster(jobs, type="PSOCK")
on.exit(parallel::stopCluster(cl), add=TRUE)

parallel::clusterExport(cl, c("run_one", "opts", "repo_root", "setting_path"), envir=environment())
results <- parallel::parLapplyLB(cl, seq_len(opts$replicates), function(i){
  run_one(i, opts, repo_root, setting_path)
})

summary_file <- file.path(opts$output_root, "parallel_run_summary.rds")
saveRDS(results, summary_file)
cat("Parallel replicates complete. OutputRoot=", normalizePath(opts$output_root, mustWork=FALSE), "\n", sep="")
