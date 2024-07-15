# Install necessary packages if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("tess3r", quietly = TRUE)) devtools::install_github("bcm-uga/TESS3_encho_sen")
if (!requireNamespace("sf", quietly = TRUE)) install.packages("sf")

# Loading libraries
library(argparse)
library(tess3r)
library(dplyr)
library(ggplot2)
library(sf)
library(maps)
library(mapdata)
library(raster)
library(fields)
library(maptools)

# Parse data
parser <- ArgumentParser(description = 'This program calculates pairwise Fst with hierfstat')
parser$add_argument('--input_lfmm', '-l', help = 'path to lfmm file')
parser$add_argument('--meta', '-m', help = 'path to meta file')
parser$add_argument('--output', '-o', help = 'basename for output')

xargs = parser$parse_args()

meta = read.csv(xargs$meta)
coordinates = data.matrix(meta %>% filter(to_exclude == "False") %>% dplyr::select(c('longitude','latitude')))

lfmm = read.table(xargs$input_lfmm, sep=" ")
lfmm[lfmm == 9] = NA

tess3.obj = tess3(X = lfmm, coord = coordinates, K = 1:10,
                 method = "projected.ls", ploidy = 2, rep=10, keep='best')

my.colors <- c('orangered', 'mediumturquoise', 'gold', 'magenta', 'orange', 'red', 'aquamarine1', 'hotpink', 'limegreen', 'blue')


par(mfrow = c(1, 1))
png(paste0(xargs$output, '_cross_validation_plot.png'))
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")
dev.off()



for (k in 1:10){
my.palette <- CreatePalette(my.colors[0:k], palette.length=10)
q.matrix <- qmatrix(tess3.obj, K = k)

us_map <- st_as_sf(map("state", plot = FALSE, fill = TRUE))

coordinates_sf <- st_as_sf(as.data.frame(coordinates), coords = c("longitude", "latitude"), crs = 4326)
us_map <- st_transform(us_map, st_crs(coordinates_sf))

lat_limits <- c(30, 48.5)
lon_limits <- c(-98, -78)

spatial_plot <- ggtess3Q(q.matrix, coordinates, col.palette = my.palette) +
  geom_sf(data = us_map, fill = NA, color = "black") +   
  xlim(lon_limits) + 
  ylim(lat_limits) + 
  coord_sf() + 
  geom_sf(data = coordinates_sf, size = 0.2) + 
  labs(x = "Longitude", y = "Latitude") + 
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold", size = 5, family="Arial"),
        axis.title.y = element_text(face = "bold", size = 5, family="Arial"),
        axis.text.x = element_text(size = 4, family = "Arial", color = 'black'),
        axis.text.y = element_text(size = 4, family = "Arial", color='black'),
       panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(paste0(xargs$output, '_barplot_map_', k,'.png'), plot = spatial_plot, width = 3, height = 3, dpi = 300)

}

