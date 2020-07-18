library(ape)
library(MCMCglmm)
library(ggbiplot)
library(dplyr)
library(data.table)
library(tidyverse)
library(grid)
library(gridExtra)
library(ggplot2)
library(umap)

load(file='/home/jamie/Documents/2020_gene_expression_study/dat/R.workspace.RData')

########### Load data
# Load phylogentic tree
phylo<-read.nexus("/home/jamie/Documents/2020_gene_expression_study/dat/cladogram_reduced.nex")

# Load data from cali.py program
dat <- read.csv(file = "/home/jamie/Documents/2020_gene_expression_study/dat/cali_out.csv", header = FALSE, 
                stringsAsFactors=FALSE)

# Load species list data
species.list <- read.csv(file = '/home/jamie/Documents/2020_gene_expression_study/dat/species_list.csv', header = TRUE)

############# Tidy data set and split gene data from non-gene data
dat = setNames(data.frame(t(dat[,-1])), dat[,1])
rownames(dat) <- NULL

##### Sperating out leptinotarsa_decemlineata gene data for UMAp on tissue type
dat.L.decemlineata <- subset(dat, species=='leptinotarsa_decemlineata')
dat.L.decemlineata_gene <- subset(dat.L.decemlineata, select = -c(SRR, species, SexDeterm, sex))

dat.L.decemlineata_gene[] <- lapply(dat.L.decemlineata_gene, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
sapply(gene_dat, class  )



### Sepertating out gene data from non-gene data
non_gene_dat <- subset(dat, select = c(SRR, species, SexDeterm, sex))
gene_dat <- subset(dat, select = -c(SRR, species, SexDeterm, sex))
non_gene_dat$SexD.sex <- paste(dat$SexDeterm,dat$sex)



# Transform gene tmp values from factor to numeric
gene_dat[] <- lapply(gene_dat, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
sapply(gene_dat, class  )

# Remove genes with var = 0
gene_dat <- gene_dat[,apply(gene_dat, 2, var, na.rm=TRUE) != 0]

###########################
#UMAP

raw.umap = umap(gene_dat, random_state = 123)


#order <- as.data.frame(species.list[c('order', 'species')])
#names(order)[names(order) == "species"] <- "phylo"

umap.df <- as.data.frame(raw.umap$layout)
umap.df$SexD.sex <- non_gene_dat$SexD.sex
umap.df$sex <- non_gene_dat$sex
umap.df$SexDeterm <- non_gene_dat$SexDeterm
umap.df$phylo <- non_gene_dat$species


#umap.df.order <- merge(umap.df, order, by.x="phylo")

umap.df$sex <- as.factor(umap.df$sex)
umap.df$SexDeterm <- as.factor(umap.df$SexDeterm)
umap.df$SexD.sex <- paste(umap.df$SexDeterm,umap.df$sex)

umap.df$SexD.sex[umap.df$SexD.sex == "diplo male"] <-("XY male")
umap.df$SexD.sex[umap.df$SexD.sex == "diplo female"] <-("XY female")
umap.df$SexD.sex[umap.df$SexD.sex == "haplo male"] <-("Haplodiploid male")
umap.df$SexD.sex[umap.df$SexD.sex == "haplo female"] <-("Haplodiploid female")

##################
#testing UMAP
gene_dat.wnoise <- gene_dat + matrix(rnorm(150*40, 0, 0.1), ncol = ncol(gene_dat))

raw.umap.wnoise <- as.data.frame(predict(raw.umap, gene_dat.wnoise))
raw.umap.wnoise$SexD.sex <- umap.df$SexD.sex

raw.umap.nonoise <- umap.df[c('V1', 'V2', 'SexD.sex')]

raw.umap.wnoise$condition <- 'with noise'
raw.umap.nonoise$condition <- 'no noise'

umap.projection <- rbind(raw.umap.wnoise,raw.umap.nonoise)

umap.projection$SexD.sex[umap.projection$SexD.sex == "Diplo male"] <-("XY male")
umap.projection$SexD.sex[umap.projection$SexD.sex == "Diplo female"] <-("XY female")
umap.projection$SexD.sex[umap.projection$SexD.sex == "Haplo male"] <-("Haplodiploid male")
umap.projection$SexD.sex[umap.projection$SexD.sex == "Haplo female"] <-("Haplodiploid female")

plt.projection <- ggplot(umap.projection, aes(x=(V1), y=(V2), color = SexD.sex, shape=condition)) + 
  xlab("UMAP 1") + ylab("UMAP 2") + geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + 
  theme(legend.title = element_blank()) +
  theme_bw() + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black")) +
  theme(legend.title=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plt.projection

####################################
#UMAP & MCMCglmm for leptinotarsa_decemlineata tissue type

leptinotarsa_decemlineata.umap = umap(dat.L.decemlineata_gene, random_state = 123)

l.d.umap.df <- as.data.frame(leptinotarsa_decemlineata.umap$layout)

tissue <- c('Head', 'Head', 'Head', 'Midgut', 'Midgut', 'Midgut', 'Head', 'Head', 'Head', 'Head', 'Head', 'Head', 'Midgut', 'Midgut', 'Midgut')
l.d.umap.df$tissue <- tissue

plt.l.d.tissue <- ggplot(l.d.umap.df, aes(x=(V1), y=(V2), color = tissue)) + xlab("UMAP 1") + ylab("UMAP 2") + 
  geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.title = element_blank()) +
  theme_bw() + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black")) +
  theme(legend.title=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  
plt.l.d.tissue

model3.4.l.d.tissue <- MCMCglmm(cbind(V1, V2)
                            ~trait-1 + trait:tissue,
                            data=l.d.umap.df, rcov=~us(trait):units,
                            family=c("gaussian", "gaussian"),
                            nitt=60010000,burnin=10000,thin=40000)

posterior.mode(model3.4.l.d.tissue$Sol)
HPDinterval(model3.4.l.d.tissue$Sol, 0.95)

######################################

plt.umap.SexD_sex <- ggplot(umap.df, aes(x=(V1), y=(V2), color = SexD.sex)) + xlab("UMAP 1") + ylab("UMAP 2") + 
  geom_point(, show.legend = FALSE) + stat_ellipse(type = "norm") + theme(legend.title = element_blank()) +
  theme_bw() + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black")) +
  theme(legend.title=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotate(geom="text", x=9, y=7.5, label="L. decemlineata", color="black") +
  annotate(geom="text", x=4, y=-1, label="L. decemlineata", color="black") +
  annotate(geom="text", x=10, y=-4, label="R. robini", color="black") +
  annotate(geom="text", x=1, y=-6, label="D. ponderosae", color="black")

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

hist(umap.df$V1, main='', xlab='UMAP 1' )
hist(umap.df$V2, main='', xlab='UMAP 2')

shapiro.test(umap.df$V1)
shapiro.test(umap.df$V2)

############################################################################3
###### UMAP on L. decemlineata by tissue type


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
                   nitt=60010000,burnin=10000,thin=40000)


# Output assesment
options(max.print = 99999)
autocorr(model3.1$Sol)
autocorr(model3.1$VCV)

par(mar=c(1,1,1,1))
plot(model3.1$Sol)
plot(model3.1$VCV)

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
                             nitt=60010000,burnin=10000,thin=40000)

posterior.mode(model3.1$VCV)
HPDinterval(model3.1$VCV)
posterior.mode(model3.2.str_pri$VCV)
HPDinterval(model3.2.str_pri$VCV)


#################################################
#Test fixed effects

summary(model3.1)

# Addative genetic and residual variance
posterior.mode(model3.1$VCV)
HPDinterval(model3.1$VCV)


# Test effect of sex*sex determination system on transcriptome composition
posterior.mode(model3.1$Sol)
HPDinterval(model3.1$Sol, 0.1)


###################################################################33
# Testing phylogenetic effect via deviance information criterion

model3.3.no.phy <- MCMCglmm(cbind(V1, V2)
                            ~trait-1 + trait:SexDeterm + trait:sex + trait:SexDeterm:sex,
                            data=umap.df, rcov=~us(trait):units,
                            family=c("gaussian", "gaussian"),
                            nitt=60010000,burnin=10000,thin=40000)

model3.1$DIC
model3.3.no.phy$DIC
