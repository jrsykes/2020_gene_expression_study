library(ape)
library(MCMCglmm)
library(ggbiplot)
library(dplyr)
#library(rlist)
library(tidyverse)

phylo<-read.nexus("/home/jamie/Documents/2020_gene_expression_study/dat/cladogram.nex")

#dat <- read.csv(file = "/home/jamie/Documents/2020_gene_expression_study/cali_test/CaliWD/out/cali_out.csv", header = FALSE)
dat <- read.csv(file = "/home/jamie/Documents/2020_gene_expression_study/dat/cali_out.csv", header = FALSE, 
                stringsAsFactors=FALSE)

##### Tidy data set and split gene data from non-gene data
dat = setNames(data.frame(t(dat[,-1])), dat[,1])
rownames(dat) <- NULL

non_gene_dat <- subset(dat, select = c(SRR, species, SexDeterm, sex))
gene_dat <- subset(dat, select = -c(SRR, species, SexDeterm, sex))

####### Transform all values from factor to numeric
gene_dat[] <- lapply(gene_dat, function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else x
})
#sapply(gene_dat, class  )


#### Merege gene and non-gene data
cleaned_dat <- bind_cols(non_gene_dat, gene_dat)

### Set row names as SRR
rownames(cleaned_dat) <- cleaned_dat[,1]
cleaned_dat[,1] <- NULL

#rownames(non_gene_dat) <- non_gene_dat[,1]
#non_gene_dat[,1] <- NULL

clean_gene_dat <- subset(cleaned_dat, select = -c(species, SexDeterm, sex))

# Remove genes with var = 0
clean_gene_dat <- clean_gene_dat[,apply(clean_gene_dat, 2, var, na.rm=TRUE) != 0]


# PCA
prin_comp <- prcomp(t(clean_gene_dat), scale. = T)

principal.components <- prin_comp$rotation
prin_comp$rotation[,1:10]

ggbiplot(prin_comp, scale = 0, pc.biplot = TRUE, labels = row.names(t(clean_gene_dat)), varname.size = 0)

std_dev <- prin_comp$sdev
pr_var <- std_dev^2

#Total var explained/PC
pr_var[1:14]

# % var explained by first n PCs
prop_varex <- pr_var/sum(pr_var)
prop_varex[1:14]
sum(prop_varex[1:14])

# Scree plot
plot(prop_varex, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")

SRR = non_gene_dat[,1]
principal.components <- cbind(principal.components, SRR)
colnames(principal.components) <- c(1:ncol(principal.components))

1:ncol(principal.components)

non_gene_dat[,1]
principal.components[1:2,1:2]
non_gene_dat[1:5,1:4]

final.df <- merge(non_gene_dat, principal.components, by=rownames)

########## Phylo
inv.phylo<-inverseA(phylo,nodes="TIPS",scale=TRUE)  

###### MCMCglmm
prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))

model<-MCMCglmm(SexDeterm~random=~phylo*sex*,
  family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),prior=prior,
  data=dat,nitt=5000000,burnin=1000,thin=500)





