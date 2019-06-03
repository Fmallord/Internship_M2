#install.packages("ape") #to use mantel tests
library("ape")


###############################
##### heatmap, neighbour joining tree and MDS on the matrix of Fst coefficient between 2 islands

Names=c("Brava", "Fogo", "Santiago", "Maio", "Sal", "BoaVista", "SaoNicolau", "SaoVicente", "SantoAntao")

az=matrix(0, nrow=9, ncol=9, dimnames=list(Names, Names))
az[1,1]=0

az= t(matrix(c(0.00000, 0, 0, 0, 0, 0, 0, 0, 0,
0.24423,   0.00000, 0, 0, 0, 0, 0, 0, 0,
0.27437,   0.31716,   0.00000, 0, 0 ,0, 0, 0, 0,
0.73745,   0.72597,   0.46507,   0.00000, 0, 0, 0, 0, 0,
0.72892,   0.71877,   0.48969,   0.58557,   0.00000, 0, 0, 0, 0,
0.65476,   0.66291,   0.46531,   0.50127,   0.30693,   0.00000, 0, 0, 0,
0.62804,   0.63299,   0.44645,   0.53047,   0.38534,   0.19250,   0.00000, 0, 0,
0.64064,   0.66312,   0.48032,   0.54071,   0.42154,   0.42431,   0.44587,   0.00000, 0,
0.64017,   0.65889,   0.49911,   0.53710,   0.46218,   0.48372,   0.49928,   0.06633,   0.00000), nrow=9, ncol=9, dimnames=list(Names, Names)))
for (i in 1:9) {
  for (j in 1:9) {
    if (az[i,j]==0) {
      az[i,j]=az[j,i]     # pour avoir la matrice diagonale remplie des 2 côtés
    }
  }
}

library(lattice)
levelplot(az, col.regions=heat.colors(37))
