#install.packages("geosphere")
library("geosphere")
#install.packages("readxl")
library("readxl")
library("stringr")
detach("package:ncf", unload=TRUE)
#install.packages("ape") #to use mantel tests
library("ape")
transcrs <- read.table("/media/fmallordy/DATA1/FMallordy/Figures_Paul/TranscriptSwadesh_List_CV2010-2018_TOTAL_PostProcValentin_FINAL_12022019_unaffiliated.txt", sep=";", header=TRUE, row.names = 1)
indivs <- read_excel("/media/fmallordy/DATA1/FMallordy/Figures_Paul/FamQuest_CapeVerde_MasterFile2010-2018_FullFINAL_study_year.xlsx")
indivs=data.frame(indivs)

dist_gen <- read.table("/media/fmallordy/DATA1/FMallordy/Figures_Paul/CapeVerde2010-2018_Omni25_PostQCstage3_no_monomorph_FINAL_Pruned50-10-0025.asd.dist", header = TRUE)

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

test1=cmdscale(dist_genn, k=5)

a=dist(test1)
b=as.matrix(a)
write.table(b, "/media/fmallordy/DATA1/FMallordy/Figures_Paul/true_gen.txt") 
