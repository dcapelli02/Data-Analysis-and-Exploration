###  Genomic profiling in Adenocarcinoma ###
###  GSE32863                            ###

setwd("C:/Users/Daniele/Desktop/2024 - 25 Secondo Semestre/Data Analysis and Exploration/Lung Adercinoma")

library("GEOquery")
library("useful")

gse<- getGEO(file = "GSE32863_series_matrix.txt.gz", getGPL = FALSE)

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

## TGO BACK TO THE METHODS ##

ex <- exprs(gse)
dim(ex)

png("Boxplot.png", width = 800, height = 600)
boxplot(ex, outline = FALSE)
dev.off()

# Now we create the groups
indici_tumore <- grep("Lung adenocarcinoma", gse$`type:ch1`)

gruppi <- c(rep("C", length(gse$`type:ch1`)))
gruppi[indici_tumore] = "A"
gruppi <- as.factor(gruppi)


### PCA ###

pca <- prcomp(t(ex))

summary(pca)
png("Screeplot PCA.png")
screeplot(pca)
dev.off()

# draw PCA plot
png("PCA component 1-2.png")
grpcol <- c(rep("blue",29), rep("red",29))
plot(pca$x[,1], pca$x[,2], xlab="PCA1", ylab="PCA2", main="PCA for components 1&2", type="p", pch=10, col=grpcol)
text(pca$x[,1], pca$x[,2], 
     #rownames(pca$x), 
     gruppi,
     cex=0.75) 
dev.off()


# Other components
png("PCA component 2-3.png")
plot(pca$x[,2], pca$x[,3], xlab="PCA1", ylab="PCA3", main="PCA for components 1&2", type="p", pch=10, col=grpcol)
text(pca$x[,2], pca$x[,3], 
     #rownames(pca$x),
     gruppi,
     cex=0.75) 
dev.off()

## K-MEANS ##

# Two groups

k <- 2
kmeans_result <- kmeans(t(ex), k)
table(kmeans_result$cluster)

png("K-means.png")
plot(kmeans_result, data=t(ex)) + geom_text(aes(label=gruppi),
                                            hjust=0,vjust=0)
dev.off()

# Repeat with 3

k <- 3
kmeans_result <- kmeans(t(ex), k)
table(kmeans_result$cluster)

png("K-means 3 gruppi.png")
plot(kmeans_result, data=t(ex)) + geom_text(aes(label=gruppi),
                                            hjust=0,vjust=0)
dev.off()

## DENDROGRAM

dist_matrix <- dist(t(ex))
hc_result <- hclust(dist_matrix, method = "ave")
k <- 3
groups <- cutree(hc_result, k=k)
table(groups)
plot(hc_result, hang <- -1, labels=gruppi)
rect.hclust(hc_result, k = k, which = NULL, x = NULL, h = NULL, border = 2, cluster = NULL) # red boxes to show groups

# Another distance choice

hc_result_single <- hclust(dist_matrix, method = "single")
k <- 2
groups <- cutree(hc_result_single, k=k)
table(groups)
plot(hc_result_single, hang <- -1, labels=gruppi)
rect.hclust(hc_result_single, k = k, which = NULL, x = NULL, h = NULL, border = 2, cluster = NULL) # red boxes to show groups

## Complete linkage

hc_result_comp <- hclust(dist_matrix, method = "complete")
k <- 2
groups <- cutree(hc_result_comp, k=k)
table(groups)
png("Cluster Dendrogram.png")
plot(hc_result_comp, hang <- -1, labels=gruppi)
rect.hclust(hc_result_comp, k = k, which = NULL, x = NULL, h = NULL, border = 2, cluster = NULL) # red boxes to show groups
dev.off()

# With 3 groups
hc_result_comp <- hclust(dist_matrix, method = "complete")
k <- 3
groups <- cutree(hc_result_comp, k=k)
table(groups)
png("Cluster Dendrogram 3 gruppi.png")
plot(hc_result_comp, hang <- -1, labels=gruppi)
rect.hclust(hc_result_comp, k = k, which = NULL, x = NULL, h = NULL, border = 2, cluster = NULL) # red boxes to show groups
dev.off()


### RANDOM FOREST ###

library(randomForest)
set.seed(1234)
rf <- randomForest(x=t(ex), y=as.factor(gruppi), ntree=1000)
png("Rf.png", width = 800, height = 600)
plot(rf)
dev.off()

# graph of sorted importance values
png("RF importance.png", width = 800, height = 600)
plot(sort(rf$importance, decreasing=TRUE))
dev.off()

# Zoom
png("Head Importan RF.png", width = 800, height = 600)
plot(head(sort(rf$importance, decreasing=TRUE), 200))
dev.off()

# can also use: 
varImpPlot(rf)

#extract the most 'important' genes
probe.names <- rownames(rf$importance)
top200 <- probe.names[order(rf$importance, decreasing=TRUE)[1:200]]
write.csv(top200, file = "probes-top200-rifatta.txt", quote=FALSE, 
          row.names = FALSE, col.names=FALSE)

top100 <- probe.names[order(rf$importance, decreasing=TRUE)[1:100]]
write.csv(top100, file = "probes-top100-rifatta.txt", quote=FALSE, 
          row.names = FALSE, col.names=FALSE)

top400 <- probe.names[order(rf$importance, decreasing=TRUE)[1:400]]
write.csv(top400, file = "probes-top400-rifatta.txt", quote=FALSE, 
          row.names = FALSE, col.names=FALSE)

### LDA ###

library("genefilter")

tt40 <- rowttests(ex,gruppi)
#keepers <- which(p.adjust(tt40$p.value)<0.5)
keepers <- which(tt40$p.value<0.0001)
ex3 <- ex[keepers,]
tex3 <- t(ex3)

dat <- cbind(as.data.frame(tex3),gruppi)
colnames(dat)[ncol(dat)] <- "Adeno"


library("caret")
#
# Run algorithms using 10-fold cross validation
#
control <- trainControl(method="cv", number=10,
                        classProbs = TRUE,
                        savePredictions = "final")
metric <- "Accuracy"

# Long time period
fit.lda <- train(Adeno~., data=dat, method="lda",
                 metric=metric, trControl=control)

# ROC Curve

library("pROC")
pred <- fit.lda$pred

best <- fit.lda$bestTune
pred <- pred[pred$parameter == best$parameter, ]

# ROC
roc_obj <- roc(response = pred$obs, predictor = pred$A)
plot(roc_obj, main = "ROC Curve LDA", col = "blue")
auc(roc_obj)

png("ROC LDA.png", width = 800, height = 600)
plot(roc_obj, main = "ROC Curve LDA", col = "blue")
dev.off()

# Projection for LDA
proj1 <- predict(fit.lda$finalModel, dat[, 1:5652])

# Plot
plot(proj1$x[,1], 
     ylab = "LDA Axis", 
     xlab = "PCA1",
     main = "Projection on the LDA axis")
text(proj1$x[,1], labels = dat$A, col = as.numeric(dat$A) + 1)

png("LDA Axis.png", width = 800, height = 600)
plot(proj1$x[,1], 
     ylab = "LDA Axis", 
     xlab = "PCA1",
     main = "Projection on the LDA axis")
text(proj1$x[,1], labels = dat$A, col = as.numeric(dat$A) + 1)
dev.off()


# Compare with Random Forest
fit.rf <- train(Adeno~., data=dat, method="rf",
                metric=metric, trControl=control)

# ROC RF

pred_rf <- fit.rf$pred
roc_obj_rf <- roc(response = pred_rf$obs, predictor = pred_rf$A)
plot(roc_obj_rf, main = "ROC Curve RF", col = "blue")
auc(roc_obj_rf)

## Save the figure ##

png("ROC RF.png", width = 800, height = 600)
plot(roc_obj_rf, main = "ROC Curve RF", col = "blue")
dev.off()

# Comparison

results <- resamples(list(LDA=fit.lda, RF=fit.rf))
summary(results)
ggplot(results) + labs(y = "Accuracy")



### LASSO ###

fit.lasso <- train(Adeno~., data=dat, method="glmnet",
                   family = "binomial",
                   tuneGrid = expand.grid(alpha = 1,
                                          lambda = seq(0,1,by=0.05)),
                   trControl = control,
                   metric = metric)

png("Lasso.png", width = 800, height = 600)
plot(fit.lasso)
dev.off()

# Comparison with other classification methods

results <- resamples(list(RF=fit.rf, LDA=fit.lda, Lasso=fit.lasso))

summary(results)
png("Comparison.png", width = 800, height = 600)
ggplot(results) + labs(y = "Accuracy")
dev.off()


### SCUDO ###

dat_scudo <- ex[]

set.seed(123)
inTrain <- createDataPartition(gruppi, list = FALSE)
trainData <- dat_scudo[, inTrain]
testData <- dat_scudo[, -inTrain]

# analyze training set

library("rScudo")
trainRes <- scudoTrain(trainData, groups = gruppi[inTrain],
                       nTop = 250, nBottom = 250, alpha = 0.05)

trainRes

# inspect signatures
upSignatures(trainRes)[1:5,1:5]
consensusUpSignatures(trainRes)[1:5, ]

# generate and plot map of training samples
trainNet <- scudoNetwork(trainRes, N = 0.3)


png("Scudo plo.png", width = 800, height = 600)
scudoPlot(trainNet, vertex.label = NA, main = "Train")
dev.off()

# perform validation using testing samples
testRes <- scudoTest(trainRes, testData, gruppi[-inTrain],
                     nTop = 100, nBottom = 100)

testNet <- scudoNetwork(testRes, N = 0.5)
png("Scudo Test Net.png", width = 800, height = 600)
scudoPlot(testNet, vertex.label = NA, main = "Test")
dev.off()

# identify clusters on map

library("igraph")

testClust <- igraph::cluster_spinglass(testNet, spins = 2)

png("Scudo testClust.png", width = 800, height = 600)
plot(testClust, testNet, vertex.label = NA, main = "Cluster")
dev.off()

# perform classification
classRes <- scudoClassify(trainData, testData, N = 0.25,
                          nTop = 250, nBottom = 250,
                          trainGroups = gruppi[inTrain], alpha = 0.5)
caret::confusionMatrix(classRes$predicted, gruppi[-inTrain])


## FINAL COMPARISON ##
scudo_res <- caret::confusionMatrix(classRes$predicted, gruppi[-inTrain])


# 1) SCUDO accuracy
scudo_acc <- as.numeric(scudo_res$overall["Accuracy"])

# 2) data frame
scudo_df <- data.frame(
  Model    = "SCUDO",
  Accuracy = scudo_acc
)

ggplot(results) + labs(y = "Accuracy") +
  geom_point(
    data = scudo_df,
    aes(x = Model, y = Accuracy),
    colour = "black"
  ) 

# 3) plot
png("Comparison_with_SCUDO.png", width = 800, height = 600)
library(ggplot2)

ggplot(results, metric = "Accuracy") +
  labs(y = "Accuracy") +
  geom_point(
    data = scudo_df,
    aes(x = Model, y = Accuracy),
  )

dev.off()
