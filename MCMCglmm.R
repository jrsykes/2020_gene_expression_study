#library(ape)
#library(MCMCglmm)
#phylo<-read.nexus("/home/jamie/Documents/2020_gene_expression_study/dat/cladogram.nex")
#print (phylo)


dat <- read.csv("/home/jamie/Documents/2020_gene_expression_study/cali_test/CaliWD/out/CaliOut.csv")

dat <- subset(dat,select=-c(X0))

gene_dat <- tail(dat, -3)


head(mtcars)

quit()

dat.pca <- prcomp(gene_dat, center = TRUE, scale. =TRUE)

summary(dat.pca)