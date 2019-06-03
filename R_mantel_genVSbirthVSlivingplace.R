#install.packages("readxl")
library("readxl")
#install.packages("geosphere")
library("geosphere")
library("stringr")
detach("package:ncf", unload=TRUE)
#install.packages("ape") #to use mantel tests
library("ape")
transcrs <- read.table("/media/fmallordy/DATA1/FMallordy/Figures_Paul/TranscriptSwadesh_List_CV2010-2018_TOTAL_PostProcValentin_FINAL_12022019_unaffiliated.txt", sep=";", header=TRUE, row.names = 1)
indivs <- read_excel("/media/fmallordy/DATA1/FMallordy/Figures_Paul/FamQuest_CapeVerde_MasterFile2010-2018_FullFINAL_study_year.xlsx")
indivs=data.frame(indivs)

transcrs=transcrs[(rownames(transcrs) %in% indivs$DNACode),]  # permet de filtrer les mots par les inds qu'on garde du tableau excel avant de tout merger dans un dataframe
dist_gen <- read.table("/media/fmallordy/DATA1/FMallordy/Figures_Paul/true_gen.txt", header = TRUE)

# dist_gen <- read.table("/media/fmallordy/DATA1/FMallordy/Figures_Paul/CapeVerde2010-2018_Omni25_PostQCstage3_no_monomorph_FINAL_Pruned50-10-0025.asd.dist", header = TRUE)

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

# dist_gen <- dist_gen[,(colnames(dist_gen) %in% total$DNACode)]
# dist_gen <- dist_gen[(rownames(dist_gen) %in% total$DNACode),]

total=total[match(rownames(dist_gen), total$DNACode),]

###############################
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
Place=list(c_Boa_Vista, c_Brava, c_Fogo, c_Santiago, c_Santo_Antao, c_Sao_Nicolau, c_Sao_Vicente)

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
  if (total$LivingPlaceLocX[i]=="?") {
    for (j in 1:length(Place)) {
      if (levels(total$LivingPlaceIsland)[j]==total$LivingPlaceIsland[i]) {
        total$LivingPlaceLocX[i]=Place[[j]][1]
        total$LivingPlaceLocY[i]=Place[[j]][2]
      }
    }
    
  }
}

total[] <- lapply(total, gsub, pattern='°', replacement='') #enlever le signe '°' des coordonnées pr pvoir les utiliser
##### using lapply forces us to use lapply(df$column, levels)[[1]] to know the levels of the dataframe


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
## matrice des distances entre lieux de résidence
# (1st longitude, 2nd latitude) in ° (as in FamQuest, X=longitude and Y=latitude)
m_dist_place <- matrix(0, nrow=length(total$DNACode), ncol=length(total$DNACode), dimnames=list(total$DNACode, total$DNACode) )
for (i in 1:length(total$DNACode)) {    # chaque ind...
  for (j in 1:length(total$DNACode)) {  # est comparé à tous les autres...
    m_dist_place[i,j]=distCosine(c(as.numeric(total$LivingPlaceLocX[i]), as.numeric(total$LivingPlaceLocY[i])), c(as.numeric(total$LivingPlaceLocX[j]), as.numeric(total$LivingPlaceLocY[j])))
  }
}

dist_genn=as.matrix(dist_gen)

library("ncf")
partial.mantel.test(dist_genn, m_dist_birth, m_dist_place, resamp = 10000, method="spearman")   # mantel pour comparer la distance entre mots et la distance en âge (est-ce que gens plus éloignés en âge parlent plus différemment?)

