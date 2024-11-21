library(bio3d)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("USAGE: Rscript scripts/generate_amsm.R <input.dcd> <amsm_path>", call.=FALSE)
}

traj      <- read.dcd(args[1])
amsm_path <- args[2]

clustering <- geostas(traj, fit=T, pairwise=T)

write.table(clustering$amsm, sep=',', amsm_path, row.names=F, col.names=F)
