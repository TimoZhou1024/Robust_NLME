args <- commandArgs(trailingOnly=TRUE)

usage <- paste(
  "Usage:",
  "Rscript simulation/evaluate_se_stability.R <combined.rds_or_output_dir> [more paths] [--out <output_dir>]",
  sep="\n"
)

if(length(args) == 0) stop(usage)

out.idx <- which(args == "--out")
if(length(out.idx) > 0){
  if(out.idx == length(args)) stop("--out requires an output directory")
  output.dir <- args[[out.idx + 1]]
  input.paths <- args[-c(out.idx, out.idx + 1)]
} else {
  output.dir <- file.path(getwd(), "se_evaluation")
  input.paths <- args
}
if(length(input.paths) == 0) stop(usage)

dir.create(output.dir, recursive=TRUE, showWarnings=FALSE)

as_num <- function(x) suppressWarnings(as.numeric(x))

find_combined <- function(path){
  if(file.exists(path) && grepl("\\.rds$", path, ignore.case=TRUE)) return(normalizePath(path))
  candidate <- file.path(path, "combined.rds")
  if(file.exists(candidate)) return(normalizePath(candidate))
  rds <- list.files(path, pattern="combined\\.rds$", recursive=TRUE, full.names=TRUE)
  if(length(rds) == 0) stop("No combined.rds found for: ", path)
  normalizePath(rds[[1]])
}

parse_logs <- function(root){
  logs <- list.files(root, pattern="\\.log$", recursive=TRUE, full.names=TRUE)
  if(length(logs) == 0){
    return(data.frame(
      log_count=0, sim_k_runs=NA_real_, bt_elapsed_count=0,
      bt_elapsed_mean=NA_real_, bt_elapsed_min=NA_real_, bt_elapsed_max=NA_real_,
      bt1_valid_sum=NA_real_, bt2_valid_sum=NA_real_, error_lines=NA_integer_
    ))
  }

  read_log <- function(path){
    con <- file(path, open="r", encoding="UTF-16LE")
    on.exit(close(con), add=TRUE)
    x <- try(readLines(con, warn=FALSE), silent=TRUE)
    if(inherits(x, "try-error") || length(x) == 0){
      x <- readLines(path, warn=FALSE)
    }
    x
  }

  log.lines <- lapply(logs, read_log)
  lines <- unlist(log.lines, use.names=FALSE)
  get_nums <- function(pattern){
    hits <- grep(pattern, lines, value=TRUE)
    if(length(hits) == 0) return(numeric(0))
    as_num(sub(pattern, "\\1", hits))
  }

  sim.k <- get_nums(".*SIM_K_RUNS=([0-9.]+).*")
  bt.run.counts <- vapply(log.lines, function(x) length(grep("--- BT run", x)), integer(1))
  bt.elapsed <- get_nums(".*BT Rnlme elapsed= *([0-9.]+) s.*")
  bt1 <- get_nums(".*Average runs for BT1 is *([0-9.]+).*")
  bt2 <- get_nums(".*Average runs for BT2 is *([0-9.]+).*")
  errors <- grep("Error in|Exceeded|Execution halted", lines, value=TRUE)

  data.frame(
    log_count=length(logs),
    sim_k_runs=if(length(sim.k)) sim.k[[1]] else if(any(bt.run.counts > 0)) median(bt.run.counts[bt.run.counts > 0]) else NA_real_,
    bt_elapsed_count=length(bt.elapsed),
    bt_elapsed_mean=if(length(bt.elapsed)) mean(bt.elapsed, na.rm=TRUE) else NA_real_,
    bt_elapsed_min=if(length(bt.elapsed)) min(bt.elapsed, na.rm=TRUE) else NA_real_,
    bt_elapsed_max=if(length(bt.elapsed)) max(bt.elapsed, na.rm=TRUE) else NA_real_,
    bt1_valid_sum=if(length(bt1)) sum(bt1, na.rm=TRUE) else NA_real_,
    bt2_valid_sum=if(length(bt2)) sum(bt2, na.rm=TRUE) else NA_real_,
    error_lines=length(errors)
  )
}

truth_for <- function(out, param.names, p){
  truth <- out$True
  if(is.null(truth)) return(rep(NA_real_, p))
  if(!is.null(names(truth)) && all(param.names %in% names(truth))){
    return(as_num(truth[param.names]))
  }
  as_num(truth[seq_len(min(length(truth), p))])
}

eval_se_matrix <- function(est, se, truth, se.type){
  if(is.null(est) || is.null(se)) return(NULL)
  est <- as.matrix(est)
  se <- as.matrix(se)
  p <- min(ncol(est), ncol(se), length(truth))
  if(p == 0) return(NULL)

  rows <- vector("list", p)
  for(j in seq_len(p)){
    e <- as_num(est[, j])
    s <- as_num(se[, j])
    t0 <- truth[[j]]
    ok <- is.finite(e) & is.finite(s) & is.finite(t0)
    n <- sum(ok)
    emp.se <- if(n > 1) sd(e[ok]) else NA_real_
    avg.se <- if(n > 0) sqrt(mean(s[ok]^2, na.rm=TRUE)) else NA_real_
    ratio <- avg.se / emp.se
    cover <- if(n > 0) mean((e[ok] - 1.96*s[ok] <= t0) & (t0 <= e[ok] + 1.96*s[ok])) else NA_real_

    rows[[j]] <- data.frame(
      parameter=colnames(est)[[j]],
      se_type=se.type,
      n_complete=n,
      empirical_se=emp.se,
      avg_reported_se=avg.se,
      se_ratio=ratio,
      coverage=cover,
      rel_mcse_empirical_se=if(n > 1) 1 / sqrt(2 * (n - 1)) else NA_real_,
      coverage_mcse=if(n > 0 && is.finite(cover)) sqrt(cover * (1 - cover) / n) else NA_real_
    )
  }
  do.call(rbind, rows)
}

eval_est <- function(out){
  est <- as.matrix(out$Est)
  p <- ncol(est)
  param.names <- colnames(est)
  truth <- truth_for(out, param.names, p)

  rows <- vector("list", p)
  for(j in seq_len(p)){
    e <- as_num(est[, j])
    t0 <- truth[[j]]
    ok <- is.finite(e) & is.finite(t0)
    n <- sum(ok)
    mean.est <- if(n > 0) mean(e[ok]) else NA_real_
    emp.se <- if(n > 1) sd(e[ok]) else NA_real_
    bias <- mean.est - t0
    rows[[j]] <- data.frame(
      parameter=param.names[[j]],
      n_est=n,
      truth=t0,
      mean_est=mean.est,
      bias=bias,
      rel_bias_pct=if(is.finite(t0) && abs(t0) > 0) abs(bias) / abs(t0) * 100 else NA_real_,
      empirical_se=emp.se,
      rel_mcse_mean=if(n > 0 && is.finite(emp.se) && is.finite(t0) && abs(t0) > 0) emp.se / sqrt(n) / abs(t0) * 100 else NA_real_,
      rel_mcse_empirical_se=if(n > 1) 1 / sqrt(2 * (n - 1)) else NA_real_
    )
  }
  do.call(rbind, rows)
}

combined.files <- vapply(input.paths, find_combined, character(1))
all.est <- list()
all.se <- list()
all.run <- list()

for(file in combined.files){
  root <- dirname(file)
  setting <- basename(root)
  x <- readRDS(file)
  run.diag <- parse_logs(root)
  run.diag$setting <- setting
  run.diag$combined_file <- file
  all.run[[setting]] <- run.diag

  out.names <- grep("\\.out$", names(x), value=TRUE)
  for(out.name in out.names){
    out <- x[[out.name]]
    if(!is.list(out) || is.null(out$Est)) next
    model <- sub("\\.out$", "", out.name)
    est.metrics <- eval_est(out)
    est.metrics$setting <- setting
    est.metrics$model <- model
    all.est[[paste(setting, model, sep=":")]] <- est.metrics

    truth <- truth_for(out, colnames(as.matrix(out$Est)), ncol(as.matrix(out$Est)))
    se.types <- c("SD", "SD.BT", "SD.BT1", "SD.BT2")
    for(se.type in se.types){
      se.metrics <- eval_se_matrix(out$Est, out[[se.type]], truth, se.type)
      if(!is.null(se.metrics)){
        se.metrics$setting <- setting
        se.metrics$model <- model
        all.se[[paste(setting, model, se.type, sep=":")]] <- se.metrics
      }
    }
  }
}

est.table <- if(length(all.est)) do.call(rbind, all.est) else data.frame()
se.table <- if(length(all.se)) do.call(rbind, all.se) else data.frame()
run.table <- do.call(rbind, all.run)

if(nrow(run.table) > 0){
  run.table$bootstrap_se_rel_mcse_nominal <- ifelse(
    is.finite(run.table$sim_k_runs) & run.table$sim_k_runs > 1,
    1 / sqrt(2 * (run.table$sim_k_runs - 1)),
    NA_real_
  )
}

flag <- function(x, good, warn){
  if(!is.finite(x)) return("insufficient")
  if(good(x)) "pass" else if(warn(x)) "warn" else "fail"
}

if(nrow(se.table) > 0){
  se.table$se_ratio_flag <- vapply(se.table$se_ratio, flag, character(1),
                                   good=function(z) z >= 0.8 && z <= 1.25,
                                   warn=function(z) z >= 0.67 && z <= 1.5)
  se.table$empirical_se_mc_flag <- vapply(se.table$rel_mcse_empirical_se, flag, character(1),
                                          good=function(z) z <= 0.15,
                                          warn=function(z) z <= 0.25)
}

write.csv(est.table, file.path(output.dir, "estimate_metrics.csv"), row.names=FALSE)
write.csv(se.table, file.path(output.dir, "se_metrics.csv"), row.names=FALSE)
write.csv(run.table, file.path(output.dir, "run_diagnostics.csv"), row.names=FALSE)

summary.lines <- c(
  "# SE Stability Evaluation",
  "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Files",
  paste("- `", combined.files, "`", sep=""),
  "",
  "## Run Diagnostics",
  paste(capture.output(print(run.table)), collapse="\n"),
  "",
  "## Interpretation",
  "- `rel_mcse_empirical_se` estimates Monte Carlo uncertainty of the empirical SE from simulation repetitions: approximately `1 / sqrt(2 * (R - 1))`.",
  "- `bootstrap_se_rel_mcse_nominal` estimates Monte Carlo uncertainty of a bootstrap SE from `B` bootstrap samples: approximately `1 / sqrt(2 * (B - 1))`.",
  "- A reported SE is treated as practically aligned with empirical SE when `se_ratio = avg_reported_se / empirical_se` is roughly between `0.8` and `1.25` for pilot work.",
  "- Coverage is reported but should not be used for final decisions unless the number of simulation repetitions is large enough; with small R it is too noisy.",
  "",
  "## Outputs",
  "- `estimate_metrics.csv`",
  "- `se_metrics.csv`",
  "- `run_diagnostics.csv`"
)

writeLines(summary.lines, file.path(output.dir, "SE_EVALUATION_SUMMARY.md"))

cat("Wrote evaluation outputs to:", normalizePath(output.dir), "\n")
