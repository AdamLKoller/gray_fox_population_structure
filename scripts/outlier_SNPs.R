# Load necessary libraries
library(pcadapt)
library(argparse)

# Parse data
parser <- ArgumentParser(description = 'This program calculates pairwise Fst with hierfstat')
parser$add_argument('--input', '-i', help = 'input file')
parser$add_argument('--output', '-o', help = 'path to write output file')

xargs = parser$parse_args()

# Read and process data
filename = read.pcadapt(xargs$input, type = 'lfmm')
x <- pcadapt(input = filename, K = 20)

# Open a PNG device
png(xargs$output, width = 1600, height = 1200)

# Set up the layout for multiple plots
par(mfrow = c(2, 2))  # 3 rows and 2 columns

# Plot the screeplot
plot(x$singular.values, type = "b", main = "Screeplot", xlab = "Components", ylab = "Proportion of explained variance")

# Recompute pcadapt with K=3
x <- pcadapt(filename, K = 3)

# Plot the Manhattan plot
plot(-log10(x$pvalues), main = "Manhattan Plot", xlab = "SNP", ylab = "-log10(p-values)", pch = 20, col = "blue")

# Plot the histogram of p-values
hist(x$pvalues, xlab = "p-values", main = "Histogram of p-values", breaks = 50, col = "orange")


# Close the PNG device
dev.off()
