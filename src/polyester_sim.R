#######################################
#                                     #
#       Simulate with Polyester       #
#             SQANTI-SIM              #
#                                     #
#######################################

#Author: Jorge Mestre Tomas (jormart2@alumni.uv.es)

suppressMessages(library(polyester))
suppressMessages(library(Biostrings))

args <- commandArgs(trailingOnly = TRUE)
ref.trans <- args[1]
count.mat <- args[2]
out.dir <- args[3]
in.seed <- args[4]
READ_LENGTH <- 100

count.mat <- read.table(count.mat, header=F)
count.mat <- as.matrix(count.mat)

simulate_experiment_countmat(ref.trans, readmat=count.mat, outdir=out.dir,
                             paired = TRUE, readlen=READ_LENGTH, seed=in.seed) 