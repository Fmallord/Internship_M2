####################################################
##### for linguistic variants (Maio & Santiago)
#install.packages("WriteXLS")  # for large genetic tables, with thousands of rows bcz of SNPs
library("WriteXLS")
#install.packages('ggplot2')
library('ggplot2')
#install.packages("abcrf")
library("abcrf")
#install.packages("readxl")
#install.packages("tidyr")
library("tidyr")
library("readxl")
library("stringr")
#install.packages("xlsx") # to export in an excel sheet
library("xlsx")

param <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/C++/Maio1_better/param.txt", header=TRUE)
stat_L_Maio <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/raw_data/for_C++/stat_L_Maio.csv", header=TRUE)
stat_L_pairwise <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/raw_data/for_C++/stat_L_pairwise.csv", header=TRUE)
stat_L_Santiago <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/raw_data/for_C++/stat_L_Santiago.csv", header=TRUE)
stat_L_out_0_1__M1 <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/C++/Maio1_better/Maio1_better_stock/stat_L_out_0_1", header=TRUE)
stat_L_out_0__M1 <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/C++/Maio1_better/Maio1_better_stock/stat_L_out_0", header=TRUE)
stat_L_out_1__M1 <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/C++/Maio1_better/Maio1_better_stock/stat_L_out_1", header=TRUE)
stat_L_out_0_1__M2 <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/C++/Maio2_better/Maio2_better_stock/stat_L_out_0_1", header=TRUE)
stat_L_out_0__M2 <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/C++/Maio2_better/Maio2_better_stock/stat_L_out_0", header=TRUE)
stat_L_out_1__M2 <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/C++/Maio2_better/Maio2_better_stock/stat_L_out_1", header=TRUE)

real_data = cbind(stat_L_Santiago, stat_L_Maio, stat_L_pairwise)

################################# ORDONNER LES DONNÉES POUR LES COMPARER AU MIEUX ##################
# permet de comparer par simul les 2 scénarios, en classant les fichiers par n_sim et en conservant que des n_sim communes aux 2 fichiers
a=(sort(stat_L_out_0_1__M1$n_sim)==sort(stat_L_out_0_1__M2$n_sim)) # (pour checker les indices n_sim en commun)

stat_L_out_0__M1 <- stat_L_out_0__M1[(stat_L_out_0__M1$n_sim %in% stat_L_out_0_1__M1$n_sim),] #on garde que les indices communs à tous les fichiers
stat_L_out_0__M2 <- stat_L_out_0__M2[(stat_L_out_0__M2$n_sim %in% stat_L_out_0_1__M1$n_sim),] 
stat_L_out_1__M1 <- stat_L_out_1__M1[(stat_L_out_1__M1$n_sim %in% stat_L_out_0_1__M1$n_sim),] 
stat_L_out_1__M2 <- stat_L_out_1__M2[(stat_L_out_1__M2$n_sim %in% stat_L_out_0_1__M1$n_sim),] 
stat_L_out_0_1__M2 <- stat_L_out_0_1__M2[(stat_L_out_0_1__M2$n_sim %in% stat_L_out_0_1__M1$n_sim),] 

stat_L_out_0__M1 <- stat_L_out_0__M1[(stat_L_out_0__M1$n_sim != 0),] #on garde que les indices communs à tous les fichiers
stat_L_out_0__M2 <- stat_L_out_0__M2[(stat_L_out_0__M2$n_sim != 0),] 
stat_L_out_1__M1 <- stat_L_out_1__M1[(stat_L_out_1__M1$n_sim != 0),] 
stat_L_out_1__M2 <- stat_L_out_1__M2[(stat_L_out_1__M2$n_sim != 0),] 
stat_L_out_0_1__M2 <- stat_L_out_0_1__M2[(stat_L_out_0_1__M2$n_sim != 0),] 
stat_L_out_0_1__M1 <- stat_L_out_0_1__M1[(stat_L_out_0_1__M1$n_sim != 0),] 

stat_L_out_0_1__M1 <- stat_L_out_0_1__M1[order(stat_L_out_0_1__M1$n_sim),]
stat_L_out_0_1__M2 <- stat_L_out_0_1__M2[order(stat_L_out_0_1__M2$n_sim),]
stat_L_out_0__M1 <- stat_L_out_0__M1[order(stat_L_out_0__M1$n_sim),] 
stat_L_out_0__M2 <- stat_L_out_0__M2[order(stat_L_out_0__M2$n_sim),]
stat_L_out_1__M1 <- stat_L_out_1__M1[order(stat_L_out_1__M1$n_sim),]
stat_L_out_1__M2 <- stat_L_out_1__M2[order(stat_L_out_1__M2$n_sim),]

################################ ANALYSE GRAPHIQUE DES DONNEES ######################
######################

stat_L_M1 <- cbind(stat_L_out_0__M1[,-1], stat_L_out_1__M1[,-1], stat_L_out_0_1__M1[,-1])
stat_L_M2 <- cbind(stat_L_out_0__M2[,-1], stat_L_out_1__M2[,-1], stat_L_out_0_1__M2[,-1])
stat_L_M1[stat_L_M1 =="NaN"] <- -1
stat_L_M2[stat_L_M2 =="NaN"] <- -1

colnames(stat_L_M1)[23] = "Fst"
colnames(stat_L_M2)[23] = "Fst"

stat_final <- rbind(stat_L_M1, stat_L_M2)

nulcol = c()
for (i in 1:length(stat_L_M1)) {
  if (length(unique(stat_L_M1[,i]))==1) {
    nulcol = c(nulcol, i)
  }
}



for (i in 1:length(stat_L_M2)) {
  if (length(unique(stat_L_M2[,i]))==1) {
    nulcol = c(nulcol, i)
  }
}

nulcol = unique(nulcol)

stat_final = stat_final[,-nulcol]
stat_L_M1 = stat_L_M1[,-nulcol]

real_data = real_data[,-nulcol]

ind1 = matrix(1,dim(stat_L_M1)[1],1)
ind2 = matrix(2,dim(stat_L_M2)[1],1)

ind = as.factor(c(ind1, ind2))

ind1 = matrix("blue",dim(stat_L_M1)[1],1)
ind2 = matrix("red",dim(stat_L_M2)[1],1)


ind = as.factor(c(ind1, ind2))

modindex <- ind
sumsta <- stat_final
data1 <- data.frame(modindex, sumsta)
model.rf1 <- abcrf(modindex~., data = data1, ntree=500, paral = TRUE, ncore = 40, lda=FALSE)
model.rf1

pred = predict(model.rf1, real_data, data1, ntree = 500, paral = TRUE)

param2 = param[stat_L_out_0__M1[,1],]
# pour un paramètre (temps t_MA_b), ou on regarde les quantiles
t <- param2$t_MA_b
sumsta <- stat_L_M1
data2 <- data.frame(t, sumsta)

model.rf.t <- regAbcrf(t~., data2, ntree=500, paral = TRUE)
model.rf.t

pred2 = predict(model.rf.t, real_data, data2)

# same_M1.pca <- prcomp(~ ., data=stat_same_M1[2:17], center = TRUE, scale = TRUE)  # là ou on fait les data PCA
# summary(same_M1.pca)
# same_M2.pca <- prcomp(~ ., data=stat_same_M2[2:17], center = TRUE, scale = TRUE)  # là ou on fait les data PCA
# summary(same_M2.pca)
# diff_M1.pca <- prcomp(~ ., data=stat_diff_M1[2:17], center = TRUE, scale = TRUE)  # là ou on fait les data PCA
# summary(diff_M1.pca)
# diff_M2.pca <- prcomp(~ ., data=stat_diff_M2[2:17], center = TRUE, scale = TRUE)  # là ou on fait les data PCA
# summary(diff_M2.pca)

stat_final.pca <- prcomp(~ ., data=stat_final, center = TRUE, scale = TRUE)  # là ou on fait les data PCA
summary(stat_final.pca)
plot(stat_final.pca$x[,1], stat_final.pca$x[,2], col=ind, pch=19, asp=1, xlab="PC1", ylab="PC2")    

temp1 <- predict(stat_final.pca,real_data)
points(temp1, col="yellow", pch=8, cex=2, lwd=3)

################################## Goodness of Fit 19 stats maio -> Prior Error Check

#install.packages("abc")
library(abc)

realtemp <- real_data[,-c(8,13,17)]
realtemp

stattemp <- stat_final[,-c(8,13,17)]
head(stattemp)

GoodFit <- gfit(real_data, stat_final, nb.replicate= 1000, tol=0.0001)
summary(GoodFit)
plot(GoodFit)
################################# Comparison; prior vs posterior data
#install.packages("fields")
library("fields")
t_col <- function(color, percent = 50, name = NULL) {
  
  #	  color = color name
  #	percent = % transparency
  #	   name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
  
  ## Save the color
  invisible(t.col)
  
}

pred.t = predict(model.rf.t, real_data, data2)
densityPlot(model.rf.t, kernel = c("epanechnikov"), from=241, to=815, real_data, data2, main = "Posterior density for t_MA_b", xlab="t_MA_b", ylab="distribution of t_MA_b") # to see what the posterior prob look like
xline(x = pred.t$med[1], col="steelblue2", lwd=2)
mycol <- t_col("lightblue", perc = 85, name = "lt.blue")
polygon(border="lightskyblue", x=c(pred.t$quantiles[1], pred.t$quantiles[1], pred.t$quantiles[2], pred.t$quantiles[2]), y=c(-10e10, 10e10, 10e10, -10e10), col=mycol)
xline(x = pred.t$expectation, col="red", lwd=2)


n0 <- param2$N0
sumsta <- stat_L_M1
data2 <- data.frame(n0, sumsta)
model.rf.n0 <- regAbcrf(n0~., data2, ntree=500, paral = TRUE)
model.rf.n0
densityPlot(model.rf.n0, real_data, data2, main = "Posterior density for N0", xlab="N0", ylab="distribution of N0") # to see what the posterior prob look like
pred.t = predict(model.rf.n0, real_data, data2)
xline(x = pred.t$med[1], col="steelblue2", lwd=2)
mycol <- t_col("lightblue", perc = 85, name = "lt.blue")
polygon(border="lightskyblue", x=c(pred.t$quantiles[1], pred.t$quantiles[1], pred.t$quantiles[2], pred.t$quantiles[2]), y=c(-10e10, 10e10, 10e10, -10e10), col=mycol)
xline(x = pred.t$expectation, col="red", lwd=2)


n_rl <- param2$n_Rl
sumsta <- stat_L_M1
data2 <- data.frame(n_rl, sumsta)
model.rf.n_rl <- regAbcrf(n_rl~., data2, ntree=500, paral = TRUE)
model.rf.n_rl
densityPlot(model.rf.n_rl, real_data, data2, main = "Posterior density for n_Rl", xlab="n_Rl", ylab="distribution of n_Rl") # to see what the posterior prob look like
pred.t = predict(model.rf.n_rl, real_data, data2)
xline(x = pred.t$med[1], col="steelblue2", lwd=2)
mycol <- t_col("lightblue", perc = 85, name = "lt.blue")
polygon(border="lightskyblue", x=c(pred.t$quantiles[1], pred.t$quantiles[1], pred.t$quantiles[2], pred.t$quantiles[2]), y=c(-10e10, 10e10, 10e10, -10e10), col=mycol)
xline(x = pred.t$expectation, col="red", lwd=2)


n1_st <- param2$N1_ST
sumsta <- stat_L_M1
data2 <- data.frame(n1_st, sumsta)
model.rf.n1_st <- regAbcrf(n1_st~., data2, ntree=500, paral = TRUE)
model.rf.n1_st
densityPlot(model.rf.n1_st, real_data, data2, main = "Posterior density for N1_ST", xlab="N1_ST", ylab="distribution of N1_ST") # to see what the posterior prob look like
pred.t = predict(model.rf.n1_st, real_data, data2)
xline(x = pred.t$med[1], col="steelblue2", lwd=2)
mycol <- t_col("lightblue", perc = 85, name = "lt.blue")
polygon(border="lightskyblue", x=c(pred.t$quantiles[1], pred.t$quantiles[1], pred.t$quantiles[2], pred.t$quantiles[2]), y=c(-10e10, 10e10, 10e10, -10e10), col=mycol)
xline(x = pred.t$expectation, col="red", lwd=2)


n2_st <- param2$N2_ST
sumsta <- stat_L_M1
data2 <- data.frame(n2_st, sumsta)
model.rf.n2_st <- regAbcrf(n2_st~., data2, ntree=500, paral = TRUE)
model.rf.n2_st
densityPlot(model.rf.n2_st, real_data, data2, main = "Posterior density for N2_ST", xlab="N2_ST", ylab="distribution of N2_ST") # to see what the posterior prob look like
pred.t = predict(model.rf.n2_st, real_data, data2)
xline(x = pred.t$med[1], col="steelblue2", lwd=2)
mycol <- t_col("lightblue", perc = 85, name = "lt.blue")
polygon(border="lightskyblue", x=c(pred.t$quantiles[1], pred.t$quantiles[1], pred.t$quantiles[2], pred.t$quantiles[2]), y=c(-10e10, 10e10, 10e10, -10e10), col=mycol)
xline(x = pred.t$expectation, col="red", lwd=2)


n1_ma <- param2$N1_MA
sumsta <- stat_L_M1
data2 <- data.frame(n1_ma, sumsta)
model.rf.n1_ma <- regAbcrf(n1_ma~., data2, ntree=500, paral = TRUE)
model.rf.n1_ma
densityPlot(model.rf.n1_ma, real_data, data2, main = "Posterior density for N1_MA", xlab="N1_MA", ylab="distribution of N1_MA") # to see what the posterior prob look like
pred.t = predict(model.rf.n1_ma, real_data, data2)
xline(x = pred.t$med[1], col="steelblue2", lwd=2)
mycol <- t_col("lightblue", perc = 85, name = "lt.blue")
polygon(border="lightskyblue", x=c(pred.t$quantiles[1], pred.t$quantiles[1], pred.t$quantiles[2], pred.t$quantiles[2]), y=c(-10e10, 10e10, 10e10, -10e10), col=mycol)
xline(x = pred.t$expectation, col="red", lwd=2)


n2_ma <- param2$N2_MA
sumsta <- stat_L_M1
data2 <- data.frame(n2_ma, sumsta)
model.rf.n2_ma <- regAbcrf(n2_ma~., data2, ntree=500, paral = TRUE)
model.rf.n2_ma
densityPlot(model.rf.n2_ma, real_data, data2, main = "Posterior density for N2_MA", xlab="N2_MA", ylab="distribution of N2_MA") # to see what the posterior prob look like
pred.t = predict(model.rf.n2_ma, real_data, data2)
xline(x = pred.t$med[1], col="steelblue2", lwd=2)
mycol <- t_col("lightblue", perc = 85, name = "lt.blue")
polygon(border="lightskyblue", x=c(pred.t$quantiles[1], pred.t$quantiles[1], pred.t$quantiles[2], pred.t$quantiles[2]), y=c(-10e10, 10e10, 10e10, -10e10), col=mycol)
xline(x = pred.t$expectation, col="red", lwd=2)


nb_ma <- param2$Nb_MA
sumsta <- stat_L_M1
data2 <- data.frame(nb_ma, sumsta)
model.rf.nb_ma <- regAbcrf(nb_ma~., data2, ntree=500, paral = TRUE)
model.rf.nb_ma
densityPlot(model.rf.nb_ma, real_data, data2, main = "Posterior density for Nb_MA", xlab="Nb_MA", ylab="distribution of Nb_MA") # to see what the posterior prob look like
pred.t = predict(model.rf.nb_ma, real_data, data2)
xline(x = pred.t$med[1], col="steelblue2", lwd=2)
mycol <- t_col("lightblue", perc = 85, name = "lt.blue")
polygon(border="lightskyblue", x=c(pred.t$quantiles[1], pred.t$quantiles[1], pred.t$quantiles[2], pred.t$quantiles[2]), y=c(-10e10, 10e10, 10e10, -10e10), col=mycol)
xline(x = pred.t$expectation, col="red", lwd=2)


n1_st_ma <- param2$N1_ST_MA
sumsta <- stat_L_M1
data2 <- data.frame(n1_st_ma, sumsta)
model.rf.n1_st_ma <- regAbcrf(n1_st_ma~., data2, ntree=500, paral = TRUE)
model.rf.n1_st_ma
densityPlot(model.rf.n1_st_ma, real_data, data2, main = "Posterior density for N1_ST_MA", xlab="N1_ST_MA", ylab="distribution of N1_ST_MA") # to see what the posterior prob look like
pred.t = predict(model.rf.n1_st_ma, real_data, data2)
xline(x = pred.t$med[1], col="steelblue2", lwd=2)
mycol <- t_col("lightblue", perc = 85, name = "lt.blue")
polygon(border="lightskyblue", x=c(pred.t$quantiles[1], pred.t$quantiles[1], pred.t$quantiles[2], pred.t$quantiles[2]), y=c(-10e10, 10e10, 10e10, -10e10), col=mycol)
xline(x = pred.t$expectation, col="red", lwd=2)


n2_st_ma <- param2$N2_ST_MA
sumsta <- stat_L_M1
data2 <- data.frame(n2_st_ma, sumsta)
model.rf.n2_st_ma <- regAbcrf(n2_st_ma~., data2, ntree=500, paral = TRUE)
model.rf.n2_st_ma
densityPlot(model.rf.n2_st_ma, real_data, data2, main = "Posterior density for N2_ST_MA", xlab="N2_ST_MA", ylab="distribution of N2_ST_MA") # to see what the posterior prob look like
pred.t = predict(model.rf.n2_st_ma, real_data, data2)
xline(x = pred.t$med[1], col="steelblue2", lwd=2)
mycol <- t_col("lightblue", perc = 85, name = "lt.blue")
polygon(border="lightskyblue", x=c(pred.t$quantiles[1], pred.t$quantiles[1], pred.t$quantiles[2], pred.t$quantiles[2]), y=c(-10e10, 10e10, 10e10, -10e10), col=mycol)
xline(x = pred.t$expectation, col="red", lwd=2)


t_nexp <- param2$t_Nexp
sumsta <- stat_L_M1
data2 <- data.frame(t_nexp, sumsta)
model.rf.t_nexp <- regAbcrf(t_nexp~., data2, ntree=500, paral = TRUE)
model.rf.t_nexp
densityPlot(model.rf.t_nexp, real_data, data2, main = "Posterior density for t_Nexp", xlab="t_Nexp", ylab="distribution of t_Nexp") # to see what the posterior prob look like
pred.t = predict(model.rf.t_nexp, real_data, data2)
xline(x = pred.t$med[1], col="steelblue2", lwd=2)
mycol <- t_col("lightblue", perc = 85, name = "lt.blue")
polygon(border="lightskyblue", x=c(pred.t$quantiles[1], pred.t$quantiles[1], pred.t$quantiles[2], pred.t$quantiles[2]), y=c(-10e10, 10e10, 10e10, -10e10), col=mycol)
xline(x = pred.t$expectation, col="red", lwd=2)


t_ma_e <- param2$t_MA_e
sumsta <- stat_L_M1
data2 <- data.frame(t_ma_e, sumsta)
model.rf.t_ma_e <- regAbcrf(t_ma_e~., data2, ntree=500, paral = TRUE)
model.rf.t_ma_e
densityPlot(model.rf.t_ma_e, real_data, data2, main = "Posterior density for t_MA_e", xlab="t_MA_e", ylab="distribution of t_MA_e") # to see what the posterior prob look like
pred.t = predict(model.rf.t_ma_e, real_data, data2)
xline(x = pred.t$med[1], col="steelblue2", lwd=2)
mycol <- t_col("lightblue", perc = 85, name = "lt.blue")
polygon(border="lightskyblue", x=c(pred.t$quantiles[1], pred.t$quantiles[1], pred.t$quantiles[2], pred.t$quantiles[2]), y=c(-10e10, 10e10, 10e10, -10e10), col=mycol)
xline(x = pred.t$expectation, col="red", lwd=2)


l_a <- param2$L_alpha
sumsta <- stat_L_M1
data2 <- data.frame(l_a, sumsta)
model.rf.l_a <- regAbcrf(l_a~., data2, ntree=500, paral = TRUE)
model.rf.l_a
densityPlot(model.rf.l_a, real_data, data2, main = "Posterior density for L_alpha", xlab="L_alpha", ylab="distribution of L_alpha") # to see what the posterior prob look like
pred.t = predict(model.rf.l_a, real_data, data2)
xline(x = pred.t$med[1], col="steelblue2", lwd=2)
mycol <- t_col("lightblue", perc = 85, name = "lt.blue")
polygon(border="lightskyblue", x=c(pred.t$quantiles[1], pred.t$quantiles[1], pred.t$quantiles[2], pred.t$quantiles[2]), y=c(-10e10, 10e10, 10e10, -10e10), col=mycol)
xline(x = pred.t$expectation, col="red", lwd=2)


l_mut <- param2$L_mut
sumsta <- stat_L_M1
data2 <- data.frame(l_mut, sumsta)
model.rf.l_mut <- regAbcrf(l_mut~., data2, ntree=500, paral = TRUE)
model.rf.l_mut
densityPlot(model.rf.l_mut, real_data, data2, main = "Posterior density for L_mut", xlab="L_mut", ylab="distribution of L_mut") # to see what the posterior prob look like
pred.t = predict(model.rf.l_mut, real_data, data2)
xline(x = pred.t$med[1], col="steelblue2", lwd=2)
mycol <- t_col("lightblue", perc = 85, name = "lt.blue")
polygon(border="lightskyblue", x=c(pred.t$quantiles[1], pred.t$quantiles[1], pred.t$quantiles[2], pred.t$quantiles[2]), y=c(-10e10, 10e10, 10e10, -10e10), col=mycol)
xline(x = pred.t$expectation, col="red", lwd=2)

g_mut <- param2$G_mut
sumsta <- stat_L_M1
data2 <- data.frame(g_mut, sumsta)
model.rf.g_mut <- regAbcrf(g_mut~., data2, ntree=500, paral = TRUE)
model.rf.g_mut
densityPlot(model.rf.g_mut, real_data, data2, main = "Posterior density for G_mut", xlab="G_mut", ylab="distribution of G_mut") # to see what the posterior prob look like
#densityPlot(model.rf.g_mut, real_data, data2, xlim=c(1e-10, 1e-7), log="x", main = "Posterior density for G_mut", xlab="G_mut", ylab="distribution of G_mut") # to see what the posterior prob look like
pred.t = predict(model.rf.g_mut, real_data, data2)
xline(x = pred.t$med[1], col="steelblue2", lwd=2)
mycol <- t_col("lightblue", perc = 85, name = "lt.blue")
polygon(border="lightskyblue", x=c(pred.t$quantiles[1], pred.t$quantiles[1], pred.t$quantiles[2], pred.t$quantiles[2]), y=c(-10e10, 10e10, 10e10, -10e10), col=mycol)
xline(x = pred.t$expectation, col="red", lwd=2)
