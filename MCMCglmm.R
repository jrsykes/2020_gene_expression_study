#library(ape)
#library(MCMCglmm)
#phylo<-read.nexus("/home/jamie/Documents/2020_gene_expression_study/dat/cladogram.nex")
#print (phylo)


dat <- read.csv(file = "/home/jamie/Documents/2020_gene_expression_study/cali_test/CaliWD/out/cali_out.csv", header = FALSE)

dat1 <- dat[-1]
row.names(dat1) <- dat$V1

names(dat1) <- as.matrix(dat1[1, ])
dat1 <- dat1[-1, ]
dat1[] <- lapply(dat1, function(x) type.convert(as.character(x)))


dat <- t(dat1)

dat.pca <- prcomp(dat, center = TRUE, scale. =TRUE)

summary(dat.pca)

