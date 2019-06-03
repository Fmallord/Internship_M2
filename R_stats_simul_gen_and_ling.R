####################################################
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
abl=c(1:20000)
param=cbind(abl, param)

stat_L_Maio <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/raw_data/for_C++/stat_L_Maio.csv", header=TRUE)
stat_L_pairwise <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/raw_data/for_C++/stat_L_pairwise.csv", header=TRUE)
stat_L_Santiago <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/raw_data/for_C++/stat_L_Santiago.csv", header=TRUE)
stat_G_Maio <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/raw_data/for_C++/stat_G_Maio.csv", header=TRUE)
stat_G_pairwise <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/raw_data/for_C++/stat_G_pairwise.csv", header=TRUE)
stat_G_Santiago <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/raw_data/for_C++/stat_G_Santiago.csv", header=TRUE)
stat_G_out_0_1__M1 <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/C++/Maio1_better/Maio1_better_stock/stat_G_out_0_1", header=TRUE)
stat_G_out_0__M1 <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/C++/Maio1_better/Maio1_better_stock/stat_G_out_0", header=TRUE)
stat_G_out_1__M1 <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/C++/Maio1_better/Maio1_better_stock/stat_G_out_1", header=TRUE)
stat_L_out_0_1__M1 <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/C++/Maio1_better/Maio1_better_stock/stat_L_out_0_1", header=TRUE)
stat_L_out_0__M1 <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/C++/Maio1_better/Maio1_better_stock/stat_L_out_0", header=TRUE)
stat_L_out_1__M1 <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/C++/Maio1_better/Maio1_better_stock/stat_L_out_1", header=TRUE)
stat_G_out_0_1__M2 <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/C++/Maio2_better/Maio2_better_stock/stat_G_out_0_1", header=TRUE)
stat_G_out_0__M2 <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/C++/Maio2_better/Maio2_better_stock/stat_G_out_0", header=TRUE)
stat_G_out_1__M2 <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/C++/Maio2_better/Maio2_better_stock/stat_G_out_1", header=TRUE)
stat_L_out_0_1__M2 <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/C++/Maio2_better/Maio2_better_stock/stat_L_out_0_1", header=TRUE)
stat_L_out_0__M2 <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/C++/Maio2_better/Maio2_better_stock/stat_L_out_0", header=TRUE)
stat_L_out_1__M2 <- read.table("/media/fmallordy/DATA1/FMallordy/files_CV_2019/C++/Maio2_better/Maio2_better_stock/stat_L_out_1", header=TRUE)
param=param[(param$abl %in% stat_L_out_0_1__M1$n_sim),]

################################# ORDONNER LES DONNÉES POUR LES COMPARER AU MIEUX ##################
# permet de comparer par simul les 2 scénarios, en classant les fichiers par n_sim et en conservant que des n_sim communes aux 2 fichiers
a=(sort(stat_G_out_0_1__M1$n_sim)==sort(stat_L_out_0_1__M1$n_sim)) # (pour checker les indices n_sim en commun)

stat_L_out_0__M1 <- stat_L_out_0__M1[(stat_L_out_0__M1$n_sim %in% stat_L_out_0_1__M1$n_sim),] #on garde que les indices communs à tous les fichiers
stat_L_out_0__M2 <- stat_L_out_0__M2[(stat_L_out_0__M2$n_sim %in% stat_L_out_0_1__M1$n_sim),] 
stat_L_out_1__M1 <- stat_L_out_1__M1[(stat_L_out_1__M1$n_sim %in% stat_L_out_0_1__M1$n_sim),] 
stat_L_out_1__M2 <- stat_L_out_1__M2[(stat_L_out_1__M2$n_sim %in% stat_L_out_0_1__M1$n_sim),] 
stat_G_out_0__M1 <- stat_G_out_0__M1[(stat_G_out_0__M1$n_sim %in% stat_L_out_0_1__M1$n_sim),] #on garde que les indices communs à tous les fichiers
stat_G_out_0__M2 <- stat_G_out_0__M2[(stat_G_out_0__M2$n_sim %in% stat_L_out_0_1__M1$n_sim),] 
stat_G_out_1__M1 <- stat_G_out_1__M1[(stat_G_out_1__M1$n_sim %in% stat_L_out_0_1__M1$n_sim),] 
stat_G_out_1__M2 <- stat_G_out_1__M2[(stat_G_out_1__M2$n_sim %in% stat_L_out_0_1__M1$n_sim),] 
stat_G_out_0_1__M2 <- stat_G_out_0_1__M2[(stat_G_out_0_1__M2$n_sim %in% stat_L_out_0_1__M1$n_sim),] 
stat_L_out_0_1__M2 <- stat_L_out_0_1__M2[(stat_L_out_0_1__M2$n_sim %in% stat_L_out_0_1__M1$n_sim),] 

stat_L_out_0_1__M1 <- stat_L_out_0_1__M1[order(stat_L_out_0_1__M1$n_sim),]
stat_L_out_0_1__M2 <- stat_L_out_0_1__M2[order(stat_L_out_0_1__M2$n_sim),]
stat_L_out_0__M1 <- stat_L_out_0__M1[order(stat_L_out_0__M1$n_sim),] 
stat_L_out_0__M2 <- stat_L_out_0__M2[order(stat_L_out_0__M2$n_sim),]
stat_L_out_1__M1 <- stat_L_out_1__M1[order(stat_L_out_1__M1$n_sim),]
stat_L_out_1__M2 <- stat_L_out_1__M2[order(stat_L_out_1__M2$n_sim),]
stat_G_out_0_1__M1 <- stat_G_out_0_1__M1[order(stat_G_out_0_1__M1$n_sim),]
stat_G_out_0_1__M2 <- stat_G_out_0_1__M2[order(stat_G_out_0_1__M2$n_sim),]
stat_G_out_0__M1 <- stat_G_out_0__M1[order(stat_G_out_0__M1$n_sim),] 
stat_G_out_0__M2 <- stat_G_out_0__M2[order(stat_G_out_0__M2$n_sim),]
stat_G_out_1__M1 <- stat_G_out_1__M1[order(stat_G_out_1__M1$n_sim),]
stat_G_out_1__M2 <- stat_G_out_1__M2[order(stat_G_out_1__M2$n_sim),]

cuta=order(stat_L_out_0_1__M1$n_sim)
cutb=order(stat_L_out_0_1__M2$n_sim)
finall=cbind(cuta, cutb)
################################ ANALYSE GRAPHIQUE DES DONNEES ######################
######################
stat_L_M1 <- cbind(stat_L_out_0__M1[,-1], stat_L_out_1__M1[,-1], stat_L_out_0_1__M1[,-1])
stat_L_M2 <- cbind(stat_L_out_0__M2[,-1], stat_L_out_1__M2[,-1], stat_L_out_0_1__M2[,-1])
stat_G_M1 <- cbind(stat_G_out_0__M1[,-1], stat_G_out_1__M1[,-1], stat_G_out_0_1__M1[,-1])
stat_G_M2 <- cbind(stat_G_out_0__M2[,-1], stat_G_out_1__M2[,-1], stat_G_out_0_1__M2[,-1])
colnames(stat_L_M1)[23] = "Fst"
colnames(stat_L_M2)[23] = "Fst"
colnames(stat_G_M1)[13] = "Fst"
colnames(stat_G_M2)[13] = "Fst"

stat_same_M1 <- cbind(stat_G_M1, stat_L_M1)
stat_same_M1 <- stat_same_M1[1:7404,]
stat_same_M2 <- cbind(stat_G_M2, stat_L_M2)
stat_same_M2 <- stat_same_M2[1:7404,]
stat_diff_M1 <- cbind(stat_G_M1, stat_L_M2)
stat_diff_M1 <- stat_diff_M1[7405:14808,]
stat_diff_M2 <- cbind(stat_G_M2, stat_L_M1)
stat_diff_M2 <- stat_diff_M2[7405:14808,]

stat_same_M1b <- stat_same_M1
stat_same_M2b <- stat_same_M2
stat_diff_M1b <- stat_diff_M1
stat_diff_M2b <- stat_diff_M2

stat_final <- rbind(stat_same_M1, stat_same_M2, stat_diff_M1, stat_diff_M2)

stat_final[stat_final=="NaN"] <- -1

nulcol = c()
for (i in 1:length(stat_same_M1)) {
  if (length(unique(stat_same_M1[,i]))==1) {
    nulcol = c(nulcol, i)
  }
}

for (i in 1:length(stat_same_M2)) {
  if (length(unique(stat_same_M2[,i]))==1) {
    nulcol = c(nulcol, i)
  }
}

for (i in 1:length(stat_diff_M1)) {
    if (length(unique(stat_diff_M1[,i]))==1) {
      nulcol = c(nulcol, i)
    }
}

for (i in 1:length(stat_diff_M2)) {
  if (length(unique(stat_diff_M2[,i]))==1) {
    nulcol = c(nulcol, i)
  }
}

nulcol = unique(nulcol)

stat_final = stat_final[,-nulcol]
stat_final[stat_final=="NaN"] <- -1

real_data=cbind(stat_G_Santiago, stat_G_Maio, stat_G_pairwise, stat_L_Santiago, stat_L_Maio, stat_L_pairwise)
real_data=real_data[,-nulcol]

ind1 = matrix(1,dim(stat_same_M1)[1],1)
ind2 = matrix(2,dim(stat_same_M2)[1],1)
ind3 = matrix(3,dim(stat_diff_M1)[1],1)
ind4 = matrix(4,dim(stat_diff_M2)[1],1)

ind = as.factor(c(ind1, ind2, ind3, ind4))

ind1 = matrix("blue",dim(stat_same_M1)[1],1)
ind2 = matrix("green",dim(stat_same_M2)[1],1)
ind3 = matrix("red",dim(stat_diff_M1)[1],1)
ind4 = matrix("orange",dim(stat_diff_M2)[1],1)

ind = as.factor(c(ind1, ind2, ind3, ind4))


modindex <- ind
sumsta <- stat_final
data1 <- data.frame(modindex, sumsta)
model.rf1 <- abcrf(modindex~., data = data1, ntree=500, paral = TRUE, ncore = 40, lda=FALSE)
model.rf1 # confusion matrix

pred = predict(model.rf1, real_data, data1, ntree = 500, paral = TRUE) # model choice

stat_final.pca <- prcomp(~ ., data=stat_final, center = TRUE, scale = TRUE)  # là ou on fait les data PCA
summary(stat_final.pca)
plot(stat_final.pca$x[,1], stat_final.pca$x[,2], col=ind, pch=3, cex=0.7, asp=1, xlab="PC1", ylab="PC2")    
#plot(stat_final.pca$x[,1], stat_final.pca$x[,2], col=ind, pch=3, cex=0.7, asp=1, xlab="PC1", ylab="PC2", xlim=c(-5, -3), ylim=c(3,5))    
legend(-15, -3, legend=c("stat_same_M1", "stat_same_M2", "stat_diff_M1", "stat_diff_M2"), col=c("blue", "green", "red", "orange"), pch=19, cex=0.9)


temp1 <- predict(stat_final.pca,real_data) # add real data to the PCA on simulated data
points(temp1, col="black", pch=8, cex=1.5, lwd=3)

GoodFit <- gfit(real_data, stat_final, nb.replicate= 1000, tol=0.0001)
summary(GoodFit)
plot(GoodFit) # to get gfit plot
