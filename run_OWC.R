# example script for running OWC
source("OWC.R")
load("owc_example_input.RData")
rho.file <- read.table("rho.txt")

Zs <- data$Z
R <- data$R
rho <- as.matrix(rho.file)
n.perm <- 1e3
W <- 1/sqrt(data$W*(1-data$W))
pval <- OWC(Zs, R, n.perm, rho, W)
print(pval)
