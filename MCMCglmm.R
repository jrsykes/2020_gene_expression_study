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
prin_comp$rotation[,1:10]

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
sum(prop_varex[1:50])

principal.components <- prin_comp$rotation[,1:50]

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

#write.csv(final.df, '/home/jamie/Documents/2020_gene_expression_study/dat/cali_out_PCA.csv', row.names = FALSE)

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
inv.phylo<-inverseA(phylo)  

########### MCMCglmm PCA

############################################################################3
# Sex determ + sex

prior2.1 <- list(G = list(G1 = list(V = diag(6), n = 5.002)), 
                 R = list(V = diag(6), n = 5.002))
model2.1 <- 
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

###################################################################33

### Full model 10 PC
# Runs
prior2.2.2 <- list(G = list(G1 = list(V = diag(10), n = 9.002)), 
                 R = list(V = diag(10), n = 9.002))

model2.2.2<-MCMCglmm(cbind(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)
                   ~trait-1 + trait:SexD.sex,
                   random=~us(trait):phylo,
                   rcov=~us(trait):units,family=c("gaussian", "gaussian", "gaussian", "gaussian","gaussian", "gaussian", "gaussian", "gaussian", "gaussian", "gaussian"),
                   ginverse=list(phylo=inv.phylo$Ainv),
                   data=final.df, 
                   prior = prior2.2.2,
                   nitt=2000000,burnin=1000,thin=500)

# Output assesment

par(mar=c(1,1,1,1))
plot(model2.2.2$Sol)
plot(model2.2.2$VCV)

summary(model2.2.2)
autocorr(model2.2.2$VCV)

# Test effect of sex*sex determination system on transcriptome composition
posterior.mode(model2.2.2$Sol)
HPDinterval(model2.2.2$Sol)

#########################################################################
#Runs
# Testing phylogenetic effect via deviance information criterion

model2.2.no.phy <- MCMCglmm(cbind(PC1, PC2, PC3, PC4, PC5, PC6)
                            ~trait-1 + trait:SexD.sex,
                            data=final.df, rcov=~us(trait):units,
                            family=c("gaussian", "gaussian", "gaussian", "gaussian", "gaussian", "gaussian"),
                            nitt=2000000,burnin=1000,thin=500)


summary(model2.2.no.phy)

###################################################################33

### Full model with informative prior
# Runs

p.var1 <- var(final.df$PC1)
p.var2 <- var(final.df$PC2)
p.var3 <- var(final.df$PC3)
p.var4 <- var(final.df$PC4)
p.var5 <- var(final.df$PC5)
p.var6 <- var(final.df$PC6)

p.var <- diag(c(p.var1,p.var2,p.var3,p.var4,p.var5,p.var6))

prior2.3.st_pri <- list(G = list(G1 = list(V = p.var*0.05, n = 6)), 
                 R = list(V = p.var*0.95, n = 6))

model2.3.st_pri<-MCMCglmm(cbind(PC1, PC2, PC3, PC4, PC5, PC6)
                   ~trait-1 + trait:SexD.sex,
                   random=~us(trait):phylo,
                   rcov=~us(trait):units,family=c("gaussian", "gaussian", "gaussian", "gaussian", "gaussian", "gaussian"),
                   ginverse=list(phylo=inv.phylo$Ainv),
                   data=final.df, 
                   prior = prior2.3.st_pri,
                   nitt=2000000,burnin=1000,thin=500)

posterior.mode(model2.2$VCV)
posterior.mode(model2.3.st_pri$VCV)
summary(model2.3.st_pri)


###########################
#UMAP
#labels <- (final.df$SexD.sex)
library(umap)
raw.umap = umap(gene_dat, random_state = 123)

##################
#testing UMAP
colnames(gene_dat)
gene_dat.wnoise <- gene_dat + matrix(rnorm(150*40, 0, 0.1), ncol = ncol(gene_dat))

raw.umap.wnoise <- as.data.frame(predict(raw.umap, gene_dat.wnoise))
raw.umap.wnoise$Sex.SexD <- umap.df$Sex.SexD

raw.umap.nonoise <- umap.df[c('V1', 'V2', 'Sex.SexD')]

raw.umap.wnoise$condition <- 'with noise'
raw.umap.nonoise$condition <- 'no noise'

dim(raw.umap.wnoise)
head(raw.umap.nonoise, 15)

umap.projection <- rbind(raw.umap.wnoise,raw.umap.nonoise)

plt.projection <- ggplot(umap.projection, aes(x=(V1), y=(V2), color = Sex.SexD, shape=condition)) + xlab("UMAP1") + ylab("UMAP2") + 
   geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.title = element_blank())

plt.projection

######################################3

species.list <- read.csv(file = '/home/jamie/Documents/2020_gene_expression_study/dat/species_list.csv', header = TRUE)
order <- as.data.frame(species.list[c('order', 'species')])
names(order)[names(order) == "species"] <- "phylo"

umap.df <- as.data.frame(raw.umap$layout)
umap.df$Sex.SexD <- final.df$SexD.sex
umap.df$sex <- final.df$sex
umap.df$SexDeterm <- final.df$SexDeterm
umap.df$phylo <- final.df$phylo

umap.df$Sex.SexD[umap.df$Sex.SexD == "diplo male"] <-("Diplo male")
umap.df$Sex.SexD[umap.df$Sex.SexD == "diplo female"] <-("Diplo female")
umap.df$Sex.SexD[umap.df$Sex.SexD == "haplo male"] <-("Haplo male")
umap.df$Sex.SexD[umap.df$Sex.SexD == "haplo female"] <-("Haplo female")

umap.df.order <- merge(umap.df, order, by.x="phylo")


plt.umap.SexD_sex <- ggplot(umap.df, aes(x=(V1), y=(V2), color = Sex.SexD)) + xlab("UMAP1") + ylab("UMAP2") + 
  geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.title = element_blank())

plt.umap.SexD_sex

plt.umap.order <- ggplot(umap.df.order, aes(x=(V1), y=(V2), color = order)) + xlab("UMAP1") + ylab("UMAP2") + 
  geom_point()# + stat_ellipse(type = "norm")

plt.umap.order


plt.umap.sex <- ggplot(umap.df, aes(x=(V1), y=(V2), color = sex)) + xlab("UMAP1") + ylab("UMAP2") + 
  geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm")

plt.umap.sex  

plt.umap.sexD <- ggplot(umap.df, aes(x=(V1), y=(V2), color = SexDeterm)) + xlab("UMAP1") + ylab("UMAP2") + 
  geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm")

plt.umap.sexD

hist(1/(umap.df$V1))
hist(1/umap.df$V2)

shapiro.test(umap.df$V1)
shapiro.test(umap.df$V2)

########### MCMCglmm UMAP

############################################################################3
# Sex determ + sex + interaction

prior3.1 <- list(G = list(G1 = list(V = diag(2), n = 1.002)), 
                 R = list(V = diag(2), n = 1.002))

model3.1<-MCMCglmm(cbind(V1, V2)
                   ~trait-1 + trait:SexDeterm + trait:sex + trait:SexDeterm:sex,
                   random=~us(trait):phylo,
                   rcov=~us(trait):units,family=c("gaussian", "gaussian"),
                   ginverse=list(phylo=inv.phylo$Ainv),
                   data=umap.df, 
                   prior = prior3.1,
                   nitt=10000000,burnin=10000,thin=500)

# Output assesment
options(max.print = 99999)
autocorr(model3.1$VCV)

par(mar=c(1,1,1,1))
plot(model3.1$Sol)
plot(model3.1$VCV)

#######################################################################
# Testing effect of reciprical transformation

model3.1.2<-MCMCglmm(cbind(1/V1, 1/V2)
                   ~trait-1 + trait:SexDeterm + trait:sex + trait:SexDeterm:sex,
                   random=~us(trait):phylo,
                   rcov=~us(trait):units,family=c("gaussian", "gaussian"),
                   ginverse=list(phylo=inv.phylo$Ainv),
                   data=umap.df, 
                   prior = prior3.1
                   ,nitt=10000000,burnin=10000,thin=500)

posterior.mode(model3.1$VCV)
posterior.mode(model3.1.2$VCV)

# Output assesment
options(max.print = 99999)
autocorr(model3.1.2$VCV)

par(mar=c(1,1,1,1))
plot(model3.1.2$sol)
plot(model3.1.2$VCV)
###################################################################33

### Model with informative prior. Strong genetic effect.

p.var1 <- var(umap.df$V1)
p.var2 <- var(umap.df$V2)
p.var <- diag(c(p.var1,p.var2))

prior3.2.str_pri <- list(G = list(G1 = list(V = p.var*0.05, n = 2)), 
                        R = list(V = p.var*0.95, n = 2))

model3.2.str_pri <- MCMCglmm(cbind(V1, V2)
                   ~trait-1 + trait:SexDeterm + trait:sex + trait:SexDeterm:sex,
                   random=~us(trait):phylo,
                   rcov=~us(trait):units,family=c("gaussian", "gaussian"),
                   ginverse=list(phylo=inv.phylo$Ainv),
                   data=umap.df, 
                   prior = prior3.2.str_pri,
                   nitt=10000000,burnin=10000,thin=500)

posterior.mode(model3.1$VCV)
posterior.mode(model3.2.str_pri$VCV)



#################################################
#Test fixed effects

summary(model3.1)

# Addative genetic and residual variance
posterior.mode(model3.1$VCV)
HPDinterval(model3.1$VCV)

# Test effect of sex*sex determination system on transcriptome composition
posterior.mode(model3.1$Sol)
HPDinterval(model3.1$Sol)


###################################################################33
# Testing phylogenetic effect via deviance information criterion

model3.3.no.phy <- MCMCglmm(cbind(V1, V2)
                            ~trait-1 + trait:SexDeterm + trait:sex + trait:SexDeterm:sex,
                            data=umap.df, rcov=~us(trait):units,
                            family=c("gaussian", "gaussian"),
                            nitt=10000000,burnin=10000,thin=500)

model3.1$DIC
model3.3.no.phy$DIC
