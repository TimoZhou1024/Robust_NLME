.robust_nlme_deriv_cache <- new.env(parent=emptyenv())

get_Deriv_cached <- function(expr, vars){
  key <- paste(paste(expr, collapse="\n"), paste(vars, collapse=","), sep="\n---\n")
  if(!exists(key, envir=.robust_nlme_deriv_cache, inherits=FALSE)){
    assign(key, Deriv(expr, vars), envir=.robust_nlme_deriv_cache)
  }
  get(key, envir=.robust_nlme_deriv_cache, inherits=FALSE)
}
