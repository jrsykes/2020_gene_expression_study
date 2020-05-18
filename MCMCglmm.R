library(ape)
library(MCMCglmm)
library(ggbiplot)
library(dplyr)
library(data.table)
library(tidyverse)
library(grid)
library(gridExtra)
library(ggplot2)

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
#gene_dat$SexD.sex <- paste(dat$SexDeterm,dat$sex)

# Transform gene tmp values from factor to numeric
gene_dat[] <- lapply(gene_dat, function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else x
})
sapply(gene_dat, class  )

# Remove genes with var = 0
gene_dat <- gene_dat[,apply(gene_dat, 2, var, na.rm=TRUE) != 0]


################ Reduce dimentionality of gene transcriptiion data via principal component analysis
# Run PCA
prin_comp <- prcomp(t(gene_dat), scale. = T)

# Asses first n principal components
principal.components <- prin_comp$rotation
prin_comp$rotation[,1:7]

# Standard deviation and variance of PCs
std_dev <- prin_comp$sdev
pr_var <- std_dev^2

# Scree plot
plot(prop_varex, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")

#Total var explained/PC
pr_var[1:10]

# % var explained by first n PCs
prop_varex <- pr_var/sum(pr_var)
prop_varex[1:6]
sum(prop_varex[1:6])

principal.components <- prin_comp$rotation[,1:6]

# Plot first two principal comonenets
#ggbiplot(prin_comp, scale = 0, pc.biplot = TRUE, labels = row.names(t(clean_gene_dat)), varname.size = 0)

#ggbiplot(prin_comp, scale = 0, pc.biplot = TRUE, labels = final.df ['SexD.sex'] , varname.size = 0)




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


############################################
# Asses distributions of Principal Components
hist(final.df$PC1)
shapiro.test(final.df$PC1)

hist(final.df$PC2)
shapiro.test(final.df$PC2)

hist(final.df$PC3)
shapiro.test(final.df$PC3)


#################################################################################################
# Sex determ * sex

p1.1.2 <- ggplot(final.df, aes(x=(PC1), y=(PC2), color = SexD.sex)) + xlab("PC1") + ylab("PC2") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p1.1.3 <- ggplot(final.df, aes(x=(PC1), y=(PC3), color = SexD.sex)) + xlab("PC1") + ylab("PC3") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p1.1.4 <- ggplot(final.df, aes(x=(PC1), y=(PC4), color = SexD.sex)) + xlab("PC1") + ylab("PC4") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p1.1.5 <- ggplot(final.df, aes(x=(PC1), y=(PC5), color = SexD.sex)) + xlab("PC1") + ylab("PC5") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p1.2.3 <- ggplot(final.df, aes(x=(PC2), y=(PC3), color = SexD.sex)) + xlab("PC2") + ylab("PC3") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p1.2.4 <- ggplot(final.df, aes(x=(PC2), y=(PC4), color = SexD.sex)) + xlab("PC2") + ylab("PC4") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p1.2.5 <- ggplot(final.df, aes(x=(PC2), y=(PC5), color = SexD.sex)) + xlab("PC2") + ylab("PC5") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p1.3.4 <- ggplot(final.df, aes(x=(PC3), y=(PC4), color = SexD.sex)) + xlab("PC3") + ylab("PC4") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
  

p1.3.5 <- ggplot(final.df, aes(x=(PC3), y=(PC1), color = SexD.sex)) + xlab("PC3") + ylab("PC5") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") +
  theme(legend.text=element_text(size=17)) + theme(legend.title = element_blank())

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

mylegend <- g_legend(p1.3.5)

plot1 <- grid.arrange(p1.1.2,p1.1.3,p1.1.4,p1.1.5,p1.2.3,p1.2.4,p1.2.5,p1.3.4,mylegend, ncol = 3, nrow = 3)

plot1
#################################################################################################
# Sex determ
p2.1.2 <- ggplot(final.df, aes(x=(PC1), y=(PC2), color = SexDeterm)) + xlab("PC1") + ylab("PC2") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p2.1.3 <- ggplot(final.df, aes(x=(PC1), y=(PC3), color = SexDeterm)) + xlab("PC1") + ylab("PC3") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p2.1.4 <- ggplot(final.df, aes(x=(PC1), y=(PC4), color = SexDeterm)) + xlab("PC1") + ylab("PC4") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p2.1.5 <- ggplot(final.df, aes(x=(PC1), y=(PC5), color = SexDeterm)) + xlab("PC1") + ylab("PC5") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p2.2.3 <- ggplot(final.df, aes(x=(PC2), y=(PC3), color = SexDeterm)) + xlab("PC2") + ylab("PC3") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p2.2.4 <- ggplot(final.df, aes(x=(PC2), y=(PC4), color = SexDeterm)) + xlab("PC2") + ylab("PC4") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p2.2.5 <- ggplot(final.df, aes(x=(PC2), y=(PC5), color = SexDeterm)) + xlab("PC2") + ylab("PC5") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p2.3.4 <- ggplot(final.df, aes(x=(PC3), y=(PC4), color = SexDeterm)) + xlab("PC3") + ylab("PC4") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p2.3.5 <- ggplot(final.df, aes(x=(PC3), y=(PC5), color = SexDeterm)) + xlab("PC3") + ylab("PC5") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm")

#################################################################################################
# sex

p3.1.2 <- ggplot(final.df, aes(x=(PC1), y=(PC2), color = sex)) + xlab("PC1") + ylab("PC2") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p3.1.3 <- ggplot(final.df, aes(x=(PC1), y=(PC3), color = sex)) + xlab("PC1") + ylab("PC3") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p3.1.4 <- ggplot(final.df, aes(x=(PC1), y=(PC4), color = sex)) + xlab("PC1") + ylab("PC4") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p3.1.5 <- ggplot(final.df, aes(x=(PC1), y=(PC5), color = sex)) + xlab("PC1") + ylab("PC5") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p3.2.3 <- ggplot(final.df, aes(x=(PC2), y=(PC3), color = sex)) + xlab("PC2") + ylab("PC3") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p3.2.4 <- ggplot(final.df, aes(x=(PC2), y=(PC4), color = sex)) + xlab("PC2") + ylab("PC4") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p3.2.5 <- ggplot(final.df, aes(x=(PC2), y=(PC5), color = sex)) + xlab("PC2") + ylab("PC5") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p3.3.4 <- ggplot(final.df, aes(x=(PC3), y=(PC4), color = sex)) + xlab("PC3") + ylab("PC4") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.position = 'none')
p3.3.5 <- ggplot(final.df, aes(x=(PC3), y=(PC5), color = sex)) + xlab("PC3") + ylab("PC5") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm")

plot2 <- grid.arrange(p2.1.2,p2.1.3,p2.1.4,p2.1.5,p2.2.3,p2.2.4,p2.2.5,p2.3.4,p2.3.5, ncol = 3, nrow = 3)
plot3 <- grid.arrange(p3.1.2,p3.1.3,p3.1.4,p3.1.5,p3.2.3,p3.2.4,p3.2.5,p3.3.4,p3.3.5, ncol = 3, nrow = 3)

#################################################################################################

# Create inverse matrix of phylogenetic correlation data from phylo object
inv.phylo<-inverseA(phylo,nodes="TIPS",scale=TRUE)  
inv.phylo
    ########### MCMCglmm

############################################################################3
# Sex determ + sex

prior2.1 <- list(G = list(G1 = list(V = diag(6), n = 5.002)), 
                 R = list(V = diag(6), n = 5.002))

model2.1<-MCMCglmm(cbind(PC1, PC2, PC3, PC4, PC5, PC6)
                   ~trait-1 + trait:SexDeterm + trait:sex,
                   random=~us(trait):phylo,
                   rcov=~us(trait):units,family=c("gaussian", "gaussian", "gaussian", "gaussian", "gaussian", "gaussian"),
                   ginverse=list(phylo=inv.phylo$Ainv),
                   data=final.df, 
                   prior = prior2.1
                   ,nitt=2000000,burnin=1000,thin=500)

# Output assesment
options(max.print = 99999)
autocorr(model2.1$VCV)

par(mar=c(1,1,1,1))
plot(model2.1$Sol)
plot(model2.1$VCV)

summary(model2.1)

# Addative genetic and residual variance
posterior.mode(model2.1$VCV)
HPDinterval(model2.1$VCV)

# Test effect of sex*sex determination system on transcriptome composition
posterior.mode(model2.1$Sol)
HPDinterval(model2.1$Sol)

###################################################################33

### Full model
# Runs
prior2.2 <- list(G = list(G1 = list(V = diag(6), n = 5.002)), 
                 R = list(V = diag(6), n = 5.002))

model2.2<-MCMCglmm(cbind(PC1, PC2, PC3, PC4, PC5, PC6)
                                  ~trait-1 + trait:SexD.sex,
                                  random=~us(trait):phylo,
                                        rcov=~us(trait):units,family=c("gaussian", "gaussian", "gaussian", "gaussian", "gaussian", "gaussian"),
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



#########################################################################
#Runs
# Testing phylogenetic effect via deviance information criterion

model2.2.no.phy <- MCMCglmm(cbind(PC1, PC2, PC3, PC4, PC5, PC6, PC7)~trait-1 + trait:SexD.sex,
                   data=final.df, 
                   nitt=500000,burnin=1000,thin=500)


summary(model2.2.no.phy)



