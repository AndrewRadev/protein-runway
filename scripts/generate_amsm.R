library(bio3d)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("USAGE: Rscript scripts/generate_amsm.R <input.dcd> <amsm_path>", call.=FALSE)
}

traj      <- read.dcd(args[1])
amsm_path <- args[2]

# We could just use amsm.xyz here, but it doesn't fit the atom positions before
# measuring movement similarity, so let's run the full function and extract the
# matrix. The actual clustering is fast enough that the wasted work doesn't
# matter.
clustering <- geostas(traj, fit=T, pairwise=T)

write.table(clustering$amsm, sep=',', amsm_path, row.names=F, col.names=F)
