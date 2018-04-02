

FScore_Wardp <- function (Class, U){
  # calculates F- measure
  Ntotal <- NROW(U)
  totalPairs <- Ntotal*(Ntotal-1)/2
  NClass <- unique(Class)
  NU <- unique(U)
  
  pairsU <- c()
  for (i in NU)
  {
    selInd <- which(U==i)
    if (length(selInd)>1)
    {
      pairsU<- cbind(pairsU,combn(which(U==i),2))
    }
  }
  NpairsU <- NCOL(pairsU)
  pairsC <- c()
  for (i in NClass)
  {
    selInd <- which(Class==i)
    if (length(selInd)>1)
    {
      pairsC<- cbind(pairsC,combn(which(Class==i),2))
    }
  }
  NpairsC <- NCOL(pairsC)
  intersectPairs<-sum(duplicated(t(cbind(pairsC,pairsU))))
  
  
  
  fmeasure <- 2*intersectPairs/(NpairsC + NpairsU)
  
  return(fmeasure)
  
}

FScore_Wardp_Total <- function (Class, U)
{
  FScore_group_value <- 0
  for (i in 1:NCOL(U))
  {
    FScore_group_value <- max(FScore_Wardp(Class, U[,i]),FScore_group_value)
  }
  return(FScore_group_value)
}
