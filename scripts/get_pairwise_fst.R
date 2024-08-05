#!/usr/bin/env Rscript

# Install packages
#install.packages("argparse")
#install.packages("hierfstat")

# libraries
library(argparse)
library(hierfstat)
library(dplyr)

# parse data
parser <- ArgumentParser(description= 'This progrom calculates pairwise Fst with hierfstat ')
parser$add_argument('--input', '-i', help= 'input file')
parser$add_argument('--output','-o', help='path to write output file')

xargs = parser$parse_args()

geno = read.csv(xargs$input, row.names=1, header=TRUE)


# bootstrap pairwise Fst calculation
#fst_table_boot = boot.ppfst(geno)
#write.table(fst_table_boot$ul,file=paste(xargs$output,'_ul'), sep = ',')
#write.table(fst_table_boot$ll,file=paste(xargs$output,'_ll'), sep = ',')

# Actual calculation (Nei's)
fst_table_nei = pairwise.neifst(geno)
write.table(fst_table_nei,file=paste(xargs$output,sep=''), sep = ',')

# Actual calculation (Weir and Cockerham)
#fst_table_WC = pairwise.WCfst(geno)
#write.table(fst_table_WC,file=paste(xargs$output,'_wc',sep=''), sep = ',')

#
