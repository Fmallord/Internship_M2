#install.packages("readxl")
library("readxl")
library("stringr")
#install.packages("ape") #to use mantel tests
library("ape")
transcrs <- read.table("/media/fmallordy/DATA1/FMallordy/Figures_Paul/TranscriptSwadesh_List_CV2010-2018_TOTAL_PostProcValentin_FINAL_12022019_unaffiliated.txt", sep=";", header=TRUE)
indivs <- read_excel("/media/fmallordy/DATA1/FMallordy/Figures_Paul/FamQuest_CapeVerde_MasterFile2010-2018_FullFINAL_study_year.xlsx")
indivs=data.frame(indivs)

dist_gen <- read.table("/media/fmallordy/DATA1/FMallordy/Figures_Paul/true_gen.txt", header = TRUE)
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

dist_gen <- dist_gen[,(colnames(dist_gen) %in% total$DNACode)]
dist_gen <- dist_gen[(rownames(dist_gen) %in% total$DNACode),]

total=total[match(rownames(dist_gen), total$DNACode),]

dist_genn=as.matrix(dist_gen)


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
## matrice des différences d'âge, en années
m_age <- matrix(0, nrow=length(total$DNACode), ncol=length(total$DNACode), dimnames=list(total$DNACode, total$DNACode) )
for (i in 1:length(total$DNACode)) {    # chaque ind...
  for (j in 1:length(total$DNACode)) {  # est comparé à tous les autres...
    m_age[i,j]=abs(total$Age[i]-total$Age[j])
        }
}

library("ncf")
mantel.test(m_dist, m_age, resamp = 100000)   # mantel pour comparer la distance lingusitique et la distance en âge (est-ce que gens plus éloignés en âge parlent plus différemment?)
mantel.test(dist_genn, m_age, resamp = 10000)   # mantel pour comparer la distance génétique et la distance en âge (est-ce que gens plus éloignés en âge sont plus différents génétiquement?)
