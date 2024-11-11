library(bio3d)
library(jsonlite)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("USAGE: Rscript scripts/segment_with_bio3d_geostas.R <input.dcd> <output_dir>", call.=FALSE)
}

traj       <- read.dcd(args[1])
output_dir <- args[2]

dir.create(output_dir, recursive=T, showWarnings=F)

write_clustering <- function(clustering, filename) {
  clusters <- lapply(clustering$inds, function(ind) ind$atom)
  json     <- toJSON(clusters, pretty=T)

  f <- file(filename)
  writeLines(json, f)
  close(f)
}

clustering <- geostas(traj, k=2)

# Reuse computed atomic movement similarity matrix (AMSM):
amsm <- clustering$amsm

write_clustering(clustering, paste0(output_dir, '/clustering02.json'))

for (k in 3:9) {
  clustering <- geostas(traj, amsm=amsm, k=k)
  write_clustering(clustering, paste0(output_dir, '/clustering0', k, '.json'))
}

for (k in 10:10) {
  clustering <- geostas(traj, amsm=amsm, k=k)
  write_clustering(clustering, paste0(output_dir, '/clustering', k, '.json'))
}
