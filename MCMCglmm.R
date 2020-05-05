#library(ape)
#library(MCMCglmm)
#phylo<-read.nexus("/home/jamie/Documents/2020_gene_expression_study/dat/cladogram.nex")
#print (phylo)

dat <- read.csv(file = "/home/jamie/Documents/2020_gene_expression_study/cali_test/CaliWD/out/cali_out.csv", header = FALSE)
#dat <- t(dat)

#names <- read.csv(file = "/home/jamie/Documents/2020_gene_expression_study/cali_test/CaliWD/out/caliOut_2.csv", header = FALSE)[,1][-1]


#column_to_rownames(dat, var='X')

head (dat)




quit()

dat.pca <- prcomp(dat, center = TRUE, scale. =TRUE)

summary(dat.pca)

