# 1) find_outlier_snps()
#   - x: numeric vector of loadings (e.g. load.rda[,1])
#   - allele_freqs: data.frame or matrix of allele freqs (columns = SNP names)
#   - z: how many SDs from the mean to call an outlier
# Returns a data.frame with SNP_name, loading, and SNP_pos
find_outlier_snps <- function(x, allele_freqs, z = 2.5) {
  # compute limits
  lims <- mean(x) + c(-1, 1) * z * sd(x)
  # pick the loadings beyond those limits
  out_vals <- x[x < lims[1] | x > lims[2]]
  
  # assemble
  df <- data.frame(
    SNP_name = names(out_vals),
    loading  = as.numeric(out_vals),
    SNP_pos  = match(names(out_vals), colnames(allele_freqs)),
    stringsAsFactors = FALSE
  )
  return(df)
}


# 2) annotate_snps_with_env()
#   - cand: data.frame produced by find_outlier_snps()
#   - allele_freqs: data.frame/matrix of allele freqs
#   - env: data.frame of environmental variables (cols are predictors)
# Returns `cand` augmented with one column per env var (the correlations),
# plus “predictor” (the env var most strongly correlated) and “correlation”
annotate_snps_with_env <- function(cand, allele_freqs, env) {
  # pre-allocate correlation matrix
  n   <- nrow(cand)
  m   <- ncol(env)
  cors <- matrix(NA, nrow = n, ncol = m,
                 dimnames = list(NULL, colnames(env)))
  
  # compute all correlations
  for (i in seq_len(n)) {
    snp_gen  <- allele_freqs[, cand$SNP_name[i]]
    cors[i,] <- apply(env, 2, function(v) cor(v, snp_gen, use = "pairwise.complete.obs"))
  }
  
  # bind correlations onto cand
  out <- cbind(cand, cors)
  
  # for each SNP, find the predictor with the max absolute correlation
  abs_cors <- abs(cors)
  best_idx <- apply(abs_cors, 1, which.max)
  
  out$predictor   <- colnames(env)[best_idx]
  out$correlation <- apply(abs_cors, 1, max)
  
  return(out)
}
