rm(list=ls())



#library(ape)
#library(MCMCglmm)
#phylo<-read.nexus("/home/jamie/Documents/2020_gene_expression_study/dat/cladogram.nex")
#print (phylo)

dat <- read.csv(file = "/home/jamie/Documents/2020_gene_expression_study/cali_test/CaliWD/out/CaliOut_2.csv", header = FALSE)

#names(dat) <- as.matrix(dat[1, ])
#dat <- dat[-1, ]
#dat[] <- lapply(dat, function(x) type.convert(as.character(x)))


head (dat)

quit()

dat.pca <- prcomp(dat, center = TRUE, scale. =TRUE)

summary(dat.pca)

