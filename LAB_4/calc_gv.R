library(dplyr)
library(tidyr)
library(stringr)

calculate_genomic_vulnerability <- function(freq_df, meta_df, pop_names, model_name = "RCP26") {
  # Pivot to wider format
  freq_wide <- pivot_wider(freq_df, names_from = Allele, 
                           values_from = c(Allele.freq, Freq.e2), 
                           id_cols = c(Pop))
  
  # Extract column names and SNP names
  cols <- colnames(freq_wide)
  snps <- unique(gsub("^(Allele.freq_|Freq.e2_)", "", cols))
  
  # Order columns: Pop, then alternating Allele.freq_ and Freq.e2_
  new_order <- c("Pop", 
                 unlist(lapply(snps, function(snp) {
                   c(paste0("Allele.freq_", snp), paste0("Freq.e2_", snp))
                 })))
  
  # Remove duplicate or unwanted entries (your example removed 2nd and 3rd)
  new_order <- new_order[new_order %in% cols]  # Safely keep existing columns only
  
  # Reorder the dataframe
  freq_ordered <- freq_wide %>% dplyr::select(all_of(new_order))
  
  # Calculate absolute differences
  df_diff <- freq_ordered %>% dplyr::select(Pop)
  for (snp in snps) {
    allele_col <- paste0("Allele.freq_", snp)
    freq_col <- paste0("Freq.e2_", snp)
    if (allele_col %in% cols & freq_col %in% cols) {
      df_diff[[paste0("Abs.diff_", snp)]] <- abs(freq_ordered[[allele_col]] - freq_ordered[[freq_col]])
    }
  }
  
  # Compute row means as the genetic variation value
  gv_df <- data.frame(
    pops = pop_names,
    gv = rowMeans(df_diff[-1], na.rm = TRUE)
  )
  
  # Add factors and metadata
  gv_df$pops_order <- factor(gv_df$pops, levels = rev(meta_df$pop))
  gv_df$model <- model_name
  gv_df$region <- meta_df$region
  
  return(gv_df)
}

#gv_df_RCP26 <- calculate_genomic_vulnerability(
 # freq_df = Poptri.freq.future1,
 #meta_df = meta_pops,
 # pop_names = meta_pops$pop,
 #model_name = "RCP26"
 #)
