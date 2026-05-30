args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 1){
  stop("Usage: Rscript simulation/combine_parallel_results.R <parallel_output_root> [output_rds]")
}

root <- normalizePath(args[[1]], mustWork=TRUE)
output <- if(length(args) >= 2) args[[2]] else file.path(root, "combined.rds")

rds.files <- list.files(root, pattern="\\.rds$", recursive=TRUE, full.names=TRUE)
rds.files <- rds.files[basename(rds.files)!="combined.rds"]
if(length(rds.files) == 0){
  stop("No replicate .rds files found under: ", root)
}

combine_field <- function(values){
  values <- values[!vapply(values, is.null, logical(1))]
  if(length(values) == 0) return(NULL)
  if(all(vapply(values, function(x) is.matrix(x) || is.data.frame(x), logical(1)))){
    return(do.call(rbind, values))
  }
  if(all(vapply(values, is.atomic, logical(1)))){
    return(do.call(rbind, lapply(values, function(x) matrix(x, nrow=1))))
  }
  values
}

combine_out <- function(outputs){
  output.names <- unique(unlist(lapply(outputs, names)))
  combined <- list()
  for(name in output.names){
    values <- lapply(outputs, function(x) x[[name]])
    if(name == "True"){
      combined[[name]] <- values[[which(!vapply(values, is.null, logical(1)))[1]]]
    } else {
      combined[[name]] <- combine_field(values)
    }
  }
  combined
}

replicates <- lapply(rds.files, readRDS)
top.names <- unique(unlist(lapply(replicates, names)))

combined <- list()
for(name in top.names){
  values <- lapply(replicates, function(x) x[[name]])
  if(all(vapply(values, is.list, logical(1)))){
    combined[[name]] <- combine_out(values)
  } else {
    combined[[name]] <- combine_field(values)
  }
}

combined$source.files <- rds.files
saveRDS(combined, output)

cat("Combined", length(rds.files), "files\n")
cat("Output:", output, "\n")
