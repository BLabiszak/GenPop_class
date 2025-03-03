library(ape)
library(pegas)
library(rehh)

compute_fay_wu_H <- function(fasta_file, outgroup_file = NULL, use_outgroup = TRUE, seq_range = NULL) {
  # Load the DNA sequence data
  dna <- read.dna(fasta_file, format = "fasta")
  
  # Subset sequences if a range is provided
  if (!is.null(seq_range)) {
    dna <- dna[seq_range, ]
  }
  
  # Identify segregating sites
  seg_sites <- seg.sites(dna)
  if (length(seg_sites) == 0) {
    stop("No segregating sites found in the dataset.")
  }
  
  # Convert DNAbin to character matrices
  dna_matrix <- as.character(as.matrix(dna))
  dna_snp_matrix <- dna_matrix[, seg_sites, drop = FALSE]
  
  # Initialize lists to store ancestral and derived alleles
  ancestral_alleles <- character(ncol(dna_snp_matrix))
  derived_alleles <- character(ncol(dna_snp_matrix))
  
  if (use_outgroup && !is.null(outgroup_file)) {
    # Load outgroup sequences
    outgroup <- read.dna(outgroup_file, format = "fasta")
    outgroup_matrix <- as.character(as.matrix(outgroup))
    outgroup_snp_matrix <- outgroup_matrix[, seg_sites, drop = FALSE]
    
    # Assign ancestral and derived alleles using outgroup
    for (i in seq_len(ncol(dna_snp_matrix))) {
      outgroup_allele <- unique(outgroup_snp_matrix[, i])
      outgroup_allele <- outgroup_allele[outgroup_allele %in% c("a", "c", "g", "t")]  # Remove gaps
      
      if (length(outgroup_allele) == 1) {  # Ensure a single clear ancestral allele
        ancestral_alleles[i] <- outgroup_allele  # Assign as ancestral allele
        
        # Get allele counts in ingroup
        allele_counts <- table(dna_snp_matrix[, i])
        allele_counts <- allele_counts[names(allele_counts) %in% c("a", "c", "g", "t")]  # Use lowercase
        
        # Identify the derived allele (most frequent non-ancestral)
        non_ancestral_alleles <- names(allele_counts)[names(allele_counts) != outgroup_allele]
        if (length(non_ancestral_alleles) > 0) {
          derived_alleles[i] <- names(sort(allele_counts[non_ancestral_alleles], decreasing = TRUE))[1]
        } else {
          derived_alleles[i] <- NA
        }
      } else {
        ancestral_alleles[i] <- NA
        derived_alleles[i] <- NA
      }
    }
  } else {
    # Assign ancestral and derived alleles using major allele frequency
    for (i in seq_len(ncol(dna_snp_matrix))) {
      allele_counts <- table(dna_snp_matrix[, i])
      allele_counts <- allele_counts[names(allele_counts) %in% c("a", "c", "g", "t")]  # Use lowercase
      
      if (length(allele_counts) >= 2) {
        sorted_alleles <- sort(allele_counts, decreasing = TRUE)
        ancestral_alleles[i] <- names(sorted_alleles)[1]  # Major allele as ancestral
        derived_alleles[i] <- names(sorted_alleles)[2]  # Second most frequent allele as derived
      } else {
        ancestral_alleles[i] <- NA
        derived_alleles[i] <- NA
      }
    }
  }
  
  # Create a map file
  snp_positions <- 1:ncol(dna_snp_matrix)
  map_file <- data.frame(
    SNP_ID = paste0("SNP", snp_positions),
    CHR = 1,
    POSITION = snp_positions,
    ANCESTRAL = ancestral_alleles,
    DERIVED = derived_alleles
  )
  
  # Save the map file
  write.table(map_file, "map_file.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Create a haplotype file
  hap_file <- as.data.frame(dna_snp_matrix)
  hap_file <- cbind(ID = rownames(hap_file), hap_file)  # Add an ID column
  colnames(hap_file)[-1] <- paste0("SNP", snp_positions)
  
  # Save the haplotype file
  write.table(hap_file, "hap_file.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Convert to haplohh object
  haplohh_obj <- data2haplohh(hap_file = "hap_file.txt", map_file = "map_file.txt", allele_coding = "map")
  
  # Compute Fay and Wu’s H statistic
  sfs_results <- calc_sfs_tests(haplohh_obj, polarized = TRUE)
  
  # Extract Fay and Wu's H statistic
  fay_wu_H <- sfs_results$FAY_WU_H
  
  # Return the Fay and Wu's H statistic
  return(fay_wu_H)
}

# compute_fay_wu_H(fasta_file = "Gen_1_fay.fasta",use_outgroup = T,outgroup_file = "Gen_1_outgrup.fasta",seq_range = 1:40)
# compute_fay_wu_H(fasta_file = "Gen_2_fay.fas",use_outgroup = F,outgroup_file = "Gen_2_outgrup.fasta",seq_range = 1:40)
# compute_fay_wu_H(fasta_file = "Gen_3_fay.fas",use_outgroup = T,outgroup_file = "Gen_3_outgrup.fasta",seq_range = 1:40)


