library("readxl")
#install.packages("xlsx")
library("xlsx")


transcrs <- read.table("/media/fmallordy/DATA1/FMallordy/Figures_Paul/TranscriptSwadesh_List_CV2010-2018_TOTAL_PostProcValentin_FINAL_12022019_unaffiliated.txt", sep=";", header=TRUE)
indivs <- read_excel("/media/fmallordy/DATA1/FMallordy/Figures_Paul/FamQuest_CapeVerde_MasterFile2010-2018_FullFINAL_study_year.xlsx")
indivs=data.frame(indivs)
transcrs$X <- NULL #retire label des indivs (pas de doublon)
total=cbind(indivs, transcrs) #coller les df ensembles dans cet ordre (!! doivent avoir même ordre!!)
a=data.frame(str_split_fixed(total$LivingPlaceLoc, " - ", 2))  #spliter les colonnes en ile et localisation
names(a) <- c("LivingPlaceIsland", "LivingPlaceLocality")
total=cbind(total, a)
total$LivingPlaceLocality <- NULL  #retirer localité => on rajoute que l'info de l'île
a=data.frame(str_split_fixed(total$BirthPlaceLoc, " - ", 2))
names(a) <- c("BirthPlaceIsland", "BirthPlaceLocality")
total=cbind(total, a)
total$BirthPlaceLocality <- NULL

tronq=cbind(total[35:68])
supertot <- data.frame(matrix(ncol = 34, nrow = 147))
x <- names(total[35:68])
colnames(supertot) <- x
y <- total$DNACode
rownames(supertot) <- y

for (k in 1:length(supertot$green)) {  #pr chaque ind
  for (i in 1:length(tronq)) {         # pour chaque mot
    for (j in 1:length(levels(tronq[,i]))) {  
      if ((levels(tronq[,i])[j]==tronq[k,i]) & (tronq[k,i]==-9)) {
                supertot[k,i] = '?'
      } else if ((levels(tronq[,i])[j]==tronq[k,i]) & (tronq[k,i]!=-9)) {
        supertot[k,i] = j   # on code les différences par des nombres
    }
}
}
}

write.xlsx(supertot, "/media/fmallordy/DATA1/FMallordy/Figures_Paul/ling_AMOVA/supertot_unaffiliated.xlsx") 

################# filter by BirthIsland (to do for every island)
ab=subset(total, total$BirthPlaceIsland=="Santo Antao")
abc=supertot[rownames(supertot) %in% ab$DNACode, ]
write.xlsx(abc, "/media/fmallordy/DATA1/FMallordy/Figures_Paul/ling_AMOVA/Birth_island/supertot_Santo_Antao.xlsx") 

ab=subset(total, total$BirthPlaceIsland=="Sao Vicente")
abc=supertot[rownames(supertot) %in% ab$DNACode, ]
write.xlsx(abc, "/media/fmallordy/DATA1/FMallordy/Figures_Paul/ling_AMOVA/Birth_island/supertot_Sao_Vicente.xlsx") 

ab=subset(total, total$BirthPlaceIsland=="Sao Nicolau")
abc=supertot[rownames(supertot) %in% ab$DNACode, ]
write.xlsx(abc, "/media/fmallordy/DATA1/FMallordy/Figures_Paul/ling_AMOVA/Birth_island/supertot_Sao_Nicolau.xlsx") 

ab=subset(total, total$BirthPlaceIsland=="Sal")
abc=supertot[rownames(supertot) %in% ab$DNACode, ]
write.xlsx(abc, "/media/fmallordy/DATA1/FMallordy/Figures_Paul/ling_AMOVA/Birth_island/supertot_Sal.xlsx") 

ab=subset(total, total$BirthPlaceIsland=="Boa Vista")
abc=supertot[rownames(supertot) %in% ab$DNACode, ]
write.xlsx(abc, "/media/fmallordy/DATA1/FMallordy/Figures_Paul/ling_AMOVA/Birth_island/supertot_Boa_Vista.xlsx") 

ab=subset(total, total$BirthPlaceIsland=="Maio")
abc=supertot[rownames(supertot) %in% ab$DNACode, ]
write.xlsx(abc, "/media/fmallordy/DATA1/FMallordy/Figures_Paul/ling_AMOVA/Birth_island/supertot_Maio.xlsx") 

ab=subset(total, total$BirthPlaceIsland=="Santiago")
abc=supertot[rownames(supertot) %in% ab$DNACode, ]
write.xlsx(abc, "/media/fmallordy/DATA1/FMallordy/Figures_Paul/ling_AMOVA/Birth_island/supertot_Santiago.xlsx") 

ab=subset(total, total$BirthPlaceIsland=="Fogo")
abc=supertot[rownames(supertot) %in% ab$DNACode, ]
write.xlsx(abc, "/media/fmallordy/DATA1/FMallordy/Figures_Paul/ling_AMOVA/Birth_island/supertot_Fogo.xlsx") 

ab=subset(total, total$BirthPlaceIsland=="Brava")
abc=supertot[rownames(supertot) %in% ab$DNACode, ]
write.xlsx(abc, "/media/fmallordy/DATA1/FMallordy/Figures_Paul/ling_AMOVA/Birth_island/supertot_Brava.xlsx") 

################# filter by LivngPlaceIsland (to do for every island)
ab=subset(total, total$LivingPlaceIsland=="Santo Antao")
abc=supertot[rownames(supertot) %in% ab$DNACode, ]
write.xlsx(abc, "/media/fmallordy/DATA1/FMallordy/Figures_Paul/ling_AMOVA/Place_island/supertot_Santo_Antao.xlsx") 

ab=subset(total, total$LivingPlaceIsland=="Sao Vicente")
abc=supertot[rownames(supertot) %in% ab$DNACode, ]
write.xlsx(abc, "/media/fmallordy/DATA1/FMallordy/Figures_Paul/ling_AMOVA/Place_island/supertot_Sao_Vicente.xlsx") 

ab=subset(total, total$LivingPlaceIsland=="Sao Nicolau")
abc=supertot[rownames(supertot) %in% ab$DNACode, ]
write.xlsx(abc, "/media/fmallordy/DATA1/FMallordy/Figures_Paul/ling_AMOVA/Place_island/supertot_Sao_Nicolau.xlsx") 

ab=subset(total, total$LivingPlaceIsland=="Boa Vista")
abc=supertot[rownames(supertot) %in% ab$DNACode, ]
write.xlsx(abc, "/media/fmallordy/DATA1/FMallordy/Figures_Paul/ling_AMOVA/Place_island/supertot_Boa_Vista.xlsx") 

ab=subset(total, total$LivingPlaceIsland=="Santiago")
abc=supertot[rownames(supertot) %in% ab$DNACode, ]
write.xlsx(abc, "/media/fmallordy/DATA1/FMallordy/Figures_Paul/ling_AMOVA/Place_island/supertot_Santiago.xlsx") 

ab=subset(total, total$LivingPlaceIsland=="Brava")
abc=supertot[rownames(supertot) %in% ab$DNACode, ]
write.xlsx(abc, "/media/fmallordy/DATA1/FMallordy/Figures_Paul/ling_AMOVA/Place_island/supertot_Brava.xlsx") 

ab=subset(total, total$LivingPlaceIsland=="Fogo")
abc=supertot[rownames(supertot) %in% ab$DNACode, ]
write.xlsx(abc, "/media/fmallordy/DATA1/FMallordy/Figures_Paul/ling_AMOVA/Place_island/supertot_Fogo.xlsx") 

################# filter by Sex (a tout écrire)
ab=subset(total, total$Sex=="M")
abc=supertot[rownames(supertot) %in% ab$DNACode, ]
write.xlsx(abc, "/media/fmallordy/DATA1/FMallordy/Figures_Paul/ling_AMOVA/Sex/supertot_male.xlsx") 

ab=subset(total, total$Sex=="F")
abc=supertot[rownames(supertot) %in% ab$DNACode, ]
write.xlsx(abc, "/media/fmallordy/DATA1/FMallordy/Figures_Paul/ling_AMOVA/Sex/supertot_female.xlsx") 
