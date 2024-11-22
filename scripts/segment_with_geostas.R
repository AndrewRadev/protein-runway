library(bio3d)
library(jsonlite)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("USAGE: Rscript scripts/segment_with_bio3d_geostas.R <input.dcd> <amsm_path> <output_dir>", call.=FALSE)
}

traj       <- read.dcd(args[1])
amsm_path  <- args[2]
output_dir <- args[3]

dir.create(output_dir, recursive=T, showWarnings=F)

write_clustering <- function(clustering, filename) {
  clusters <- lapply(clustering$inds, function(ind) ind$atom)
  json     <- toJSON(clusters, pretty=T)

  f <- file(filename)
  writeLines(json, f)
  close(f)
}

amsm <- as.matrix(read.table(amsm_path, sep=',', header=F))

# K-means clustering:
for (k in 2:9) {
  clustering <- geostas(traj, amsm=amsm, k=k)
  write_clustering(clustering, paste0(output_dir, '/clustering_kmeans_0', k, '.json'))
}
for (k in 10:10) {
  clustering <- geostas(traj, amsm=amsm, k=k)
  write_clustering(clustering, paste0(output_dir, '/clustering_kmeans_', k, '.json'))
}

# Hierarchical clustering (experimental, so only 03 for now):
for (k in 2:9) {
  clustering <- geostas(traj, amsm=amsm, k=k, clustalg='hclust')
  write_clustering(clustering, paste0(output_dir, '/clustering_hier_0', k, '.json'))
}
for (k in 10:10) {
  clustering <- geostas(traj, amsm=amsm, k=k, clustalg='hclust')
  write_clustering(clustering, paste0(output_dir, '/clustering_hier_', k, '.json'))
}
