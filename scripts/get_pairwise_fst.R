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

# Perform pairwise Nei's Fst calculation
fst_table = pairwise.neifst(geno)
write.table(fst_table,file=xargs$output, sep = ',', row.names = FALSE, col.names = FALSE )



