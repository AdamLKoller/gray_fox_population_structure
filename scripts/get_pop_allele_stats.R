#!/usr/bin/env Rscript

# Install packages
#install.packages("argparse")
#install.packages("hierfstat")

# libraries
library(argparse)
library(hierfstat)
#\library(dplyr)

# parse data
parser <- ArgumentParser(description= 'This progrom calculates pairwise Fst with hierfstat ')
parser$add_argument('--input', '-i', help= 'input file')
parser$add_argument('--output','-o', help='path to write output file')

xargs = parser$parse_args()

geno = read.csv(xargs$input, row.names=1, header=TRUE)



# Calculate basic statistics
basic_stats <- basic.stats(geno, diploid = TRUE)

# Calculate allelic richness
allelic_richness_result <- allelic.richness(geno)

# Calculate private alleles

find_private_allele <- function(table) {
  # Initialize the result list with zeros
  result <- rep(0, ncol(table))
  
  # Iterate over each row
  for (i in 1:nrow(table)) {
    # Check if the max value in the row equals the sum of the row
    if (max(table[i, ]) == sum(table[i, ])) {
      # Find the column index of the max value
      max_col <- which.max(table[i, ])
        
      cat(table, max(table[i, ]), sum(table[i, ]), sep = "     ", fill = TRUE)
      # Set the corresponding index in the result list to 1
      result[max_col] <- 1
      # Return the result list
      return(result)
    }
  }
  
  # If no row meets the condition, return the result list with zeros
  return(result)
}


#print(basic_stats$pop.freq)

private_alleles = lapply(basic_stats$pop.freq, find_private_allele)
private_alleles = as.data.frame(do.call(rbind, private_alleles))

# Create a summary table
summary_table <- data.frame(
  Allelic_Richness = colMeans(allelic_richness_result$Ar),
  Observed_Heterozygosity = colMeans(basic_stats$Ho),
  Expected_Heterozygosity = colMeans(basic_stats$Hs),
  N = colMeans(basic_stats$n.ind.samp),
  Fis = colMeans(basic_stats$Fis, na.rm = TRUE),
  Private_Alleles = colSums(private_alleles)
)

# Print the summary table
print(summary_table)
write.table(summary_table,file=xargs$output, sep = "\t", row.names = FALSE, col.names = TRUE )



