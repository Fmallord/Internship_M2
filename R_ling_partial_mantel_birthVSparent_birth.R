#install.packages("readxl")
library("readxl")
library("stringr")
#install.packages("ape") #to use mantel tests
#install.packages("ncf") #to use partial mantel tests
#install.packages("geosphere")
library("geosphere")
detach("package:ncf", unload=TRUE) # to avoid confrontation btwn the 2 functions mantel.test for simple mantel tests
library("ape")

transcrs <- read.table("/media/fmallordy/DATA1/FMallordy/Figures_Paul/TranscriptSwadesh_List_CV2010-2018_TOTAL_PostProcValentin_FINAL_12022019_unaffiliated.txt", sep=";", header=TRUE, row.names = 1)
indivs <- read_excel("/media/fmallordy/DATA1/FMallordy/Figures_Paul/FamQuest_CapeVerde_MasterFile2010-2018_FullFINAL_study_year_parents.xlsx")
indivs=data.frame(indivs)

transcrs=transcrs[(rownames(transcrs) %in% indivs$DNACode),]  # permet de filtrer les mots par les inds qu'on garde du tableau excel avant de tout merger dans un dataframe

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

a=data.frame(str_split_fixed(total$FatherBirthPlaceLoc, " - ", 2))
names(a) <- c("FatherBirthPlaceIsland", "FatherBirthPlaceLocality")
total=cbind(total, a)
total$FatherBirthPlaceLocality <- NULL

a=data.frame(str_split_fixed(total$MotherBirthPlaceLoc, " - ", 2))
names(a) <- c("MotherBirthPlaceIsland", "MotherBirthPlaceLocality")
total=cbind(total, a)
total$MotherBirthPlaceLocality <- NULL

######### Coordinates for all islands
c_Brava=c(-24.7, 14.8666667)
c_Fogo=c(-24.383055555555558, 14.9330556)
c_Santiago=c(-23.62611111111111, 15.0613889)
c_Maio=c(-23.165555555555553, 15.2258333)
c_Boa_Vista=c(-22.803611111111113, 16.1033333)
c_Sal=c(-22.933333333333334, 16.7166667)
c_Sao_Nicolau=c(-24.27111111111111, 16.615)
c_Sao_Vicente=c(-24.95, 16.8333333)
c_Santo_Antao=c(-25.171111111111113, 17.07)
c_Sao_Tome=c(0, 0) #faux, juste pour rajouter une valeur
Birth=list(c_Boa_Vista, c_Brava, c_Fogo, c_Maio, c_Sal, c_Santiago, c_Santo_Antao, c_Sao_Nicolau, c_Sao_Vicente)

########## pour modifier en interne le dataframe total, en remplissant les trous par les coords des îles pr naissance et vie
for (i in 1:length(total$DNACode)) {
  if (total$BirthPlaceLocX[i]=="?") {
    for (j in 1:length(Birth)) {             #tjrs pb de levels, mais l'indiv n'est plus dans le tableau total
      if (levels(total$BirthPlaceIsland)[j]==total$BirthPlaceIsland[i]) {
        total$BirthPlaceLocX[i]=Birth[[j]][1]
        total$BirthPlaceLocY[i]=Birth[[j]][2]
      }
    }
    
  }
}

for (i in 1:length(total$DNACode)) {
  if (total$FatherBirthPlaceLocX[i]=="?") {
    for (j in 1:length(Birth)) {             #tjrs pb de levels, mais l'indiv n'est plus dans le tableau total
      if (levels(total$BirthPlaceIsland)[j]==total$FatherBirthPlaceIsland[i]) {
        total$FatherBirthPlaceLocX[i]=Birth[[j]][1]
        total$FatherBirthPlaceLocY[i]=Birth[[j]][2]
      }
    }
    
  }
}

for (i in 1:length(total$DNACode)) {
  if (total$MotherBirthPlaceLocX[i]=="?") {
    for (j in 1:length(Birth)) {             #tjrs pb de levels, mais l'indiv n'est plus dans le tableau total
      if (levels(total$BirthPlaceIsland)[j]==total$MotherBirthPlaceIsland[i]) {
        total$MotherBirthPlaceLocX[i]=Birth[[j]][1]
        total$MotherBirthPlaceLocY[i]=Birth[[j]][2]
      }
    }
    
  }
}



total[] <- lapply(total, gsub, pattern='°', replacement='') #enlever le signe '°' des coordonnées pr pvoir les utiliser
##### using lapply forces us to use lapply(df$column, levels)[[1]] to know the levels of the dataframe



## pour tous les mots
############################
## calcul matrice des dissimilarités entre individus (quelque soit la différence étymologique entre variants)
m_dist <- matrix(0, nrow=length(total$DNACode), ncol=length(total$DNACode), dimnames=list(total$DNACode, total$DNACode) )
for (i in 1:length(total$DNACode)) {    # chaque ind...
  for (j in 1:length(total$DNACode)) {  # est comparé à tous les autres...
    for (k in 1:34) {                 # pour tous les mots
      if (total[i, 34+k]!= total[j, 34+k]) {
        if (total[i, 34+k]!=-9 & total[j, 34+k]!=-9) {   # on ajoute les distances simples que qd on a les 2 mots
          m_dist[i,j]=m_dist[i,j]+1
        }
        
      }
    }
  }
}
m_dist <- m_dist/34                   # pour passer à une matrice en proportion de dissimilarité


###########################  
## matrice des distances entre lieux de naissance
# (1st longitude, 2nd latitude) in ° (as in FamQuest, X=longitude and Y=latitude)
m_dist_birth <- matrix(0, nrow=length(total$DNACode), ncol=length(total$DNACode), dimnames=list(total$DNACode, total$DNACode) )
for (i in 1:length(total$DNACode)) {    # chaque ind...
  for (j in 1:length(total$DNACode)) {  # est comparé à tous les autres...
    m_dist_birth[i,j]=distCosine(c(as.numeric(total$BirthPlaceLocX[i]), as.numeric(total$BirthPlaceLocY[i])), c(as.numeric(total$BirthPlaceLocX[j]), as.numeric(total$BirthPlaceLocY[j])))
  }
}


###########################  
## matrice des distances entre lieux de naissance des mères
# (1st longitude, 2nd latitude) in ° (as in FamQuest, X=longitude and Y=latitude)
m_dist_mother_birth <- matrix(0, nrow=length(total$DNACode), ncol=length(total$DNACode), dimnames=list(total$DNACode, total$DNACode) )
for (i in 1:length(total$DNACode)) {    # chaque ind...
  for (j in 1:length(total$DNACode)) {  # est comparé à tous les autres...
    m_dist_mother_birth[i,j]=distCosine(c(as.numeric(total$MotherBirthPlaceLocX[i]), as.numeric(total$MotherBirthPlaceLocY[i])), c(as.numeric(total$MotherBirthPlaceLocX[j]), as.numeric(total$MotherBirthPlaceLocY[j])))
  }
}

###########################  
## matrice des distances entre lieux de naissance des pères
# (1st longitude, 2nd latitude) in ° (as in FamQuest, X=longitude and Y=latitude)
m_dist_father_birth <- matrix(0, nrow=length(total$DNACode), ncol=length(total$DNACode), dimnames=list(total$DNACode, total$DNACode) )
for (i in 1:length(total$DNACode)) {    # chaque ind...
  for (j in 1:length(total$DNACode)) {  # est comparé à tous les autres...
    m_dist_father_birth[i,j]=distCosine(c(as.numeric(total$FatherBirthPlaceLocX[i]), as.numeric(total$FatherBirthPlaceLocY[i])), c(as.numeric(total$FatherBirthPlaceLocX[j]), as.numeric(total$FatherBirthPlaceLocY[j])))
  }
}


mantel.test(m_dist, m_dist_father_birth, nperm = 100000)
mantel.test(m_dist, m_dist_mother_birth, nperm = 100000)
plot(m_dist_father_birth, m_dist, pch=19, xlab="Distance between father birth places", ylab="Linguistic distance")
plot(m_dist_mother_birth, m_dist, pch=19, xlab="Distance between mother birth places", ylab="Linguistic distance")

library("ncf")
partial.mantel.test(m_dist, m_dist_father_birth, m_dist_birth, resamp = 10000, method="spearman")   # mantel pour comparer la distance entre mots et la distance en âge (est-ce que gens plus éloignés en âge parlent plus différemment?)
partial.mantel.test(m_dist, m_dist_mother_birth, m_dist_birth, resamp = 10000, method="spearman")   # mantel pour comparer la distance entre mots et la distance en âge (est-ce que gens plus éloignés en âge parlent plus différemment?)

