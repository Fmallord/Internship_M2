#install.packages("readxl")
library("readxl")
library("stringr")
library("FactoMineR")
transcrs <- read.table("/media/fmallordy/DATA1/FMallordy/Figures_Paul/TranscriptSwadesh_List_CV2010-2018_TOTAL_PostProcValentin_FINAL_12022019_unaffiliated_NA.txt", sep=";", header=TRUE, row.names = 1)
indivs <- read_excel("/media/fmallordy/DATA1/FMallordy/Figures_Paul/FamQuest_CapeVerde_MasterFile2010-2018_FullFINAL_study_year.xlsx")
indivs=data.frame(indivs)

transcrs$X <- NULL #retire label des indivs (pas de doublon)
total=cbind(indivs, transcrs) #coller les df ensembles dans cet ordre

a=data.frame(str_split_fixed(total$LivingPlaceLoc, " - ", 2))  #spliter les colonnes en ile et localisation
names(a) <- c("LivingPlaceIsland", "LivingPlaceLocality")
total=cbind(total, a)
total$LivingPlaceLocality <- NULL  #retirer localité => on rajoute que l'info de l'île

a=data.frame(str_split_fixed(total$BirthPlaceLoc, " - ", 2))
names(a) <- c("BirthPlaceIsland", "BirthPlaceLocality")
total=cbind(total, a)
total$BirthPlaceLocality <- NULL


######################### MDS

transcrs.MCA <- MCA(transcrs[,1:34], na.method = "Average")
summary(transcrs.MCA, ncp=3) #résumé des stats résumées, par défaut sur 10 variables/individus, et 3 dimensions de coords pour catég variables

##### to colorate the points
test1=transcrs.MCA$ind$coord
col1=rgb(0.11328125, 0.44140625,	0.7187500) #SantoAntao
col2=rgb(0.0000000,	0.0000000,	0.0000000)  #Santiago
col3=rgb(0.55859375,	0.89453125,	0.93359375) #SaoNicolau
col4=rgb(0.1328125,	0.70703125,	0.44921875) #SaoVicente
col5=rgb(0.7421875,	0.0859375,	0.1328125) #Fogo
col6=rgb(0.94921875,	0.5703125,	0.0000000) #Brava
col7=rgb(0.937254902,	0.450980392,	0.835294118) #Sal
col8=rgb(0.91372549,	0.352941176,	0.047058824) #Maio
col9=rgb(0.61176471,	0.258823529,	0.917647059) #BoaVista
col10=rgb(0.792156863,	0.729411765,	0.623529412) #SaoTomé
col11=rgb(0.6, 1, 0.5) # France
#couleurs_base=c(col9, col6, col5, col8, col7, col2, col1, col3, col4) #pour qd c'est en levels
couleurs_base=c(col2, col4, col1, col6, col5, col9, col3, col8, col7)  # pr quand c'est unique
couleur=vector("character", length=length(total$DNACode))
for (i in 1:length(total$DNACode)) {
  for (j in 1:length(unique(total$BirthPlaceIsland))) {
    if (total[i,70]==unique(total$BirthPlaceIsland)[j]) {
      couleur[i]=couleurs_base[j]
    }
  }
}
test1=cbind(couleur, test1)
##### end of color part

#to plot MCA by indivs
plot(transcrs.MCA$ind$coord[,1], transcrs.MCA$ind$coord[,2], main="MCA representation by individuals", xlab="Dim 1", ylab="Dim 2", pch = 19, col=test1[,1], asp=1)
legend(1.3, 0.2, legend=as.character(unique(total$BirthPlaceIsland)), col=couleurs_base, pch=19, cex=0.75)
