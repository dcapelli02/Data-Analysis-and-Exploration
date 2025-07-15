### MULTIOMICS ANALYSIS ###

setwd("C:/Users/Daniele/Desktop/2024 - 25 Secondo Semestre/Data Analysis and Exploration/Lung Adercinoma")

## Useful libraries

library("GEOquery")
library("useful")

## IMPORT THE DATA

gse<- getGEO(file = "GSE32863_series_matrix.txt.gz", getGPL = FALSE)

ex <- exprs(gse)
dim(ex)

# The idea now is to simplify a bit the dataset
# As we have data on the same patient both for the tumor tissue and the 
# non-tumor tissue, we could pick for both of them the tumor and
# for the other half the non-tumor

set.seed(1234)

subjects <- unique(gse$`subject:ch1`)

tumor_subjects_index <- sample(1:length(subjects), 30)
tumor_subjects <- subjects[tumor_subjects_index]

prova <- 1:length(gse$`subject:ch1`)
indici_tumore <- grep("Lung adenocarcinoma", gse$`type:ch1`)
individui_tumore <- prova[gse$`subject:ch1` %in% tumor_subjects]
prova_tumori <- prova[(gse$`subject:ch1` %in% tumor_subjects)]
prova_tumori <- prova_tumori[prova_tumori %in% indici_tumore]

healthy_subjects <- setdiff(subjects, tumor_subjects)
indici_sani <- setdiff(prova, indici_tumore)
individui_sani <- prova[gse$`subject:ch1` %in% healthy_subjects]
prova_sani <- prova[(gse$`subject:ch1` %in% healthy_subjects)]
prova_sani <- prova_sani[prova_sani %in% indici_sani]

index_dataset <- c(prova_tumori, prova_sani)

gse <- gse[, index_dataset]

## Work on the data ##

ex <- exprs(gse)
dim(ex)

boxplot(ex, outline = FALSE)

# The result is ok
indici_tumore <- grep("Lung adenocarcinoma", gse$`type:ch1`)  # cerca stringa in quel vettore


# Now we can create the groups
gruppi <- c(rep("C", length(gse$`type:ch1`)))
gruppi[indici_tumore] = "A"
gruppi <- as.factor(gruppi)


### METHYLATION DATA ###

library(mixOmics)

## Read the data
gse1 <- getGEO(file = "GSE32861_series_matrix.txt.gz", getGPL = FALSE)
gse1 <- gse1[, index_dataset]

ex1 <- exprs(gse1)
dim(ex1)

## Feature selection
library("genefilter")

tt40 <- rowttests(ex,gruppi)
#keepers <- which(p.adjust(tt40$p.value)<0.5)
keepers <- which(tt40$p.value<0.000001)
ex3 <- ex[keepers,]
tex3 <- t(ex3)

# Do the same on the other dataset
tt40_1 <- rowttests(ex1,gruppi)
#keepers <- which(p.adjust(tt40$p.value)<0.5)
keepers_1 <- which(tt40_1$p.value<0.01)
ex4 <- ex1[keepers_1,]
tex4 <- t(ex4)

## prepare the first argument of mixOmics functions:

MyResult1.spls <- spls(tex3, tex4, keepX = c(10, 10), keepY
                       = c(7,7))
# basic plots
plotIndiv(MyResult1.spls)
plotVar(MyResult1.spls, cex=c(3,2), legend = TRUE)

## Save the results ##
png("Omics 1.png", width = 800, height = 600)
plotVar(MyResult1.spls, cex = c(2, 2))
dev.off()

# advanced plot: use average of components from the Xand the Y-space
# and group ID ('Genotype') for labeling
plotIndiv(MyResult1.spls, group = gruppi,
          rep.space = "XY-variate", legend = TRUE)

# advanced sample plot: arrow plot to overlap X- and Yspace plots

## Arrow Plot ##

png("Frecce omics.png", width = 800, height = 600)
plotArrow(MyResult1.spls,group=gruppi, legend =
            TRUE, legend.title = 'Group',
          X.label = 'PLS comp 1', Y.label = 'PLS comp 2')
dev.off()



### VARIABLES ###

# Variables Loadings
load.X <- loadings(MyResult1.spls)$X
load.Y <- loadings(MyResult1.spls)$Y

# take the first two components
load.X12 <- as.matrix(load.X[, 1:2])
load.Y12 <- as.matrix(load.Y[, 1:2])

# Normalize
norm.X <- load.X12 / sqrt(rowSums(load.X12^2))
norm.Y <- load.Y12 / sqrt(rowSums(load.Y12^2))

# Scalar product
dot.products <- norm.X %*% t(norm.Y)

# Most opposite couples
soglia <- -0.95
opp_coppie <- which(dot.products < soglia, arr.ind = TRUE)

cat("Opposite variables:\n")
for (i in seq_len(nrow(opp_coppie))) {
  xname <- rownames(norm.X)[opp_coppie[i, 1]]
  yname <- rownames(norm.Y)[opp_coppie[i, 2]]
  score <- dot.products[opp_coppie[i, 1], opp_coppie[i, 2]]
  cat(sprintf(" - %s ↔ %s  (dot product = %.3f)\n", xname, yname, score))
}

cat("Uncorrelated Variables:\n")
for (i in seq_len(nrow(opp_coppie))) {
  xname <- rownames(norm.X)[opp_coppie[i, 1]]
  yname <- rownames(norm.Y)[opp_coppie[i, 2]]
  score <- dot.products[opp_coppie[i, 1], opp_coppie[i, 2]]
  cat(sprintf(" - %s ↔ %s\n", xname, yname))
}


# Loadings of X 
sel.vars.X <- rownames(loadings(MyResult1.spls)$X)[
  apply(loadings(MyResult1.spls)$X[, 1:2], 1, function(x) any(x != 0))
]

# Loadings of Y
sel.vars.Y <- rownames(loadings(MyResult1.spls)$Y)[
  apply(loadings(MyResult1.spls)$Y[, 1:2], 1, function(x) any(x != 0))
]

# Stampa
cat("VX variables:\n")
print(sel.vars.X)

cat("Y variables:\n")
print(sel.vars.Y)


### CCA ###

# Repeat with less genes

tt40 <- rowttests(ex,gruppi)
#keepers <- which(p.adjust(tt40$p.value)<0.5)
keepers <- which(tt40$p.value<0.0000000000000000001)
ex3 <- ex[keepers,]
tex3 <- t(ex3)

tt40_1 <- rowttests(ex1,gruppi)
#keepers <- which(p.adjust(tt40$p.value)<0.5)
keepers_1 <- which(tt40_1$p.value<0.0000000001)
ex4 <- ex1[keepers_1,]
tex4 <- t(ex4)

result.cca <- rcc(tex4, tex3, method = "ridge", 
                  lambda1 = 0.5, lambda2 = 0.05) # using the ridge method

# plot projection into canonical variate subspace
plotIndiv(result.cca)
# plot original variables' correlation with canonical variates
plotVar(result.cca)

plotArrow(result.cca, group=gruppi, legend =
            TRUE, legend.title = 'Groups',
          X.label = 'PLS comp 1', Y.label = 'PLS comp 2')

## Save figures ##

png("Omics CCA.png", width = 800, height = 600)
plotVar(result.cca, cex = c(2,2))
dev.off()

png("Arrows CCA.png", width = 800, height = 600)
plotArrow(result.cca, group=gruppi, legend =
            TRUE, legend.title = 'Groups',
          X.label = 'PLS comp 1', Y.label = 'PLS comp 2')
dev.off()

