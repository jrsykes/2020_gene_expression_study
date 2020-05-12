library(ape)
library(MCMCglmm)
library(ggbiplot)
library(dplyr)
library(data.table)
library(tidyverse)

########### Load data
# Load phylogentic tree
phylo<-read.nexus("/home/jamie/Documents/2020_gene_expression_study/dat/cladogram.nex")

# Load data from cali.py program
dat <- read.csv(file = "/home/jamie/Documents/2020_gene_expression_study/dat/cali_out.csv", header = FALSE, 
                stringsAsFactors=FALSE)

############# Tidy data set and split gene data from non-gene data
dat = setNames(data.frame(t(dat[,-1])), dat[,1])
rownames(dat) <- NULL

non_gene_dat <- subset(dat, select = c(SRR, species, SexDeterm, sex))
gene_dat <- subset(dat, select = -c(SRR, species, SexDeterm, sex))

# Transform gene tmp values from factor to numeric
gene_dat[] <- lapply(gene_dat, function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else x
})
sapply(gene_dat, class  )

# Remove genes with var = 0
gene_dat <- gene_dat[,apply(gene_dat, 2, var, na.rm=TRUE) != 0]

# Histogram of gene data
hist(gene_dat)

################ Reduce dimentionality of gene transcriptiion data via principal component analysis
# Run PCA
prin_comp <- prcomp(t(gene_dat), scale. = T)

# Asses first n principal components
principal.components <- prin_comp$rotation[,1:10]
prin_comp$rotation[,1:10]

# Standard deviation and variance of PCs
std_dev <- prin_comp$sdev
pr_var <- std_dev^2

#Total var explained/PC
pr_var[1:10]

# % var explained by first n PCs
prop_varex <- pr_var/sum(pr_var)
prop_varex[1:10]
sum(prop_varex[1:7])

# Plot first two principal comonenets
ggbiplot(prin_comp, scale = 0, pc.biplot = TRUE, labels = row.names(t(clean_gene_dat)), varname.size = 0)

ggbiplot(prin_comp, scale = 0, pc.biplot = TRUE, labels = final.df['SexD.sex'] , varname.size = 0)


# Scree plot
plot(prop_varex, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")

######## Preparing final data frame for MCMCglmm

# Convert principal component values to numeric
pc.matrix <- data.matrix(principal.components)
principal.components <- as.data.frame(pc.matrix)

# Set row names of PC data to SRR numbers
rownames(principal.components) <- non_gene_dat[,1]
principal.components <- as.data.frame(principal.components)

# Add SRR numbers as a column to allow merging with non-genetic data frame
principal.components <- tibble::rownames_to_column(principal.components, "VALUE")
colnames(principal.components)[1] <- 'SRR'

# Merge principal components with non-genetic data to create final data frame
final.df <- merge(non_gene_dat, principal.components, by='SRR')
colnames(final.df)[2] <- 'phylo'

# Set sex and sex determination system columns to factor and then merge into one column. 
# e.g. haplo male, haplo female...
final.df$sex <- as.factor(final.df$sex)
final.df$SexDeterm <- as.factor(final.df$SexDeterm)
final.df$SexD.sex <- paste(final.df$SexDeterm,final.df$sex)

# Create inverse matrix of phylogenetic correlation data from phylo object
inv.phylo<-inverseA(phylo,nodes="TIPS",scale=TRUE)  

########### MCMCglmm

### Full model
# Runs
prior2.2 <- list(G = list(G1 = list(V = diag(10), n = 9.002)), 
                 R = list(V = diag(10), n = 9.002))

model2.2<-MCMCglmm(cbind(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)~trait-1 + trait:SexD.sex,
                                  random=~us(trait):phylo,
                                  rcov=~us(trait):units,family=c("gaussian", "gaussian", "gaussian", "gaussian", "gaussian", "gaussian", "gaussian", "gaussian", "gaussian", "gaussian"),
                                  ginverse=list(phylo=inv.phylo$Ainv),
                                  data=final.df, 
                                  prior = prior2.2,
                                  nitt=2000000,burnin=1000,thin=500)

 # Output assesment

par(mar=c(1,1,1,1))
plot(model2.2$Sol)
plot(model2.2$VCV)

summary(model2.2)
autocorr(model2.2$VCV)

# Test effect of sex*sex determination system on transcriptome composition
posterior.mode(model2.2$Sol)
HPDinterval(model2.2$Sol)

############################################################################3
# Sex determ + sex

prior2.3 <- list(G = list(G1 = list(V = diag(10), n = 9.002)), 
                 R = list(V = diag(10), n = 9.002))

model2.3<-MCMCglmm(cbind(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)~
                   trait-1 + trait:SexDeterm + trait:sex,
                   random=~us(trait):phylo,
                   rcov=~us(trait):units,family=c("gaussian", "gaussian", "gaussian", "gaussian", "gaussian", "gaussian", "gaussian", "gaussian", "gaussian", "gaussian"),
                   ginverse=list(phylo=inv.phylo$Ainv),
                   data=final.df, 
                   prior = prior2.3
                   ,nitt=2000000,burnin=1000,thin=500)

# Output assesment

par(mar=c(1,1,1,1))
plot(model2.3$Sol)
plot(model2.3$VCV)

summary(model2.3)
autocorr(model2.3$VCV)

# Test effect of sex*sex determination system on transcriptome composition
posterior.mode(model2.3$Sol)
HPDinterval(model2.3$Sol)

#########################################################################
#Runs
# Testing phylogenetic effect via deviance information criterion

model2.2.no.phy <- MCMCglmm(cbind(PC1, PC2, PC3, PC4, PC5, PC6, PC7)~trait-1 + trait:SexD.sex,
                   data=final.df, 
                   nitt=500000,burnin=1000,thin=500)


summary(model2.2.no.phy)



