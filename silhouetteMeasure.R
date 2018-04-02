silhouetteMeasure <- function(AllDistancesO, Uind)
{
  Asil <- c()
  Bsil <- c()
  for (i in 1:NROW(AllDistancesO))
  {
    Asil[i] <- mean(AllDistancesO[setdiff(which(Uind==Uind[i]),i),i])
    Bsil[i] <- Inf
    for (j in setdiff(unique(Uind),Uind[i]))
    {
      Bsil[i] <- min(Bsil[i],
                     mean(AllDistancesO[which(Uind==j),i]))
    }
  }
  Asil[is.na(Asil)] <- 0
  maxABSil <- Asil
  
  indM <- which(Bsil > Asil)
  maxABSil[indM] <- Bsil[indM]
  Sil  <-  (Bsil - Asil)/maxABSil
  silhouetteInd <- mean(Sil)
  return(silhouetteInd)
}