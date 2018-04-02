Ward_p_range <- function(Matrix, reqK,range_p,
                         constraint_type="none",Constrains=c()){
  currentP <- Inf
  currentS <- -Inf
  currentU <- 1:NROW(Matrix)
  for (p in range_p)
  {
    list[U,SI]<- Ward_p(Matrix, reqK,p,constraint_type,Constrains)
    print(c("p=",p,", SI=",SI))
    if (SI > currentS)
    {
      currentS <- SI
      currentU <- U
      currentP <- p
    }
  }
  return(list(currentU,currentP,currentS))
}

Ward_p <- function(Matrix, reqK,p,constraint_type="none",Constrains=c()){
  # function build cluster set
  
  #determine parameters of the input data
  #number of samples
  Ntotal <- NROW(Matrix)
  #number of feautures
  Vtotal <- NCOL(Matrix)
  
  #initialize labels matrix/vector to return
  U <-matrix(0,nrow=Ntotal,ncol=(Ntotal-reqK+1)) 
  U[,1] <- c(1:Ntotal)
  
  #check that there are enough data for analysis
  if (Ntotal <= reqK)
  {
    return(U)
  }
  
  #initialize matrix of centroids with original data points
  Z <- Matrix
  
  #in the begining all the clusters have the size 1
  Nk <-matrix(1,nrow=Ntotal,ncol=1)#(Ntotal-reqK+1))
  
  #initialize weights
  W <-matrix(1.0/Vtotal,nrow=Ntotal,ncol=Vtotal)
  #keep original values to use in constraints algorithm,
  #all the constraints are reference to original data
  ZO <- Z
  WO <- W
  NkO <- Nk
  
  #estimate current distance between centroids
  AllDistances <- 
    Fast_Distance_Function(ZO,WO,NkO,
                           Z,W,Nk,
                           p,U,1,
                           constraint_type,Constrains)
  
  AllDistancesO <- AllDistances

  istep <- 1
  
  K = Ntotal
  while ( (K-istep+1) > reqK ) 
  {
    print(c("merge step: ",istep))
    
    indMin <- which.min(AllDistances)
    distmin <- AllDistances[indMin]
    
    indC1 = U[((indMin-1) %% Ntotal) + 1, istep]
    indC2 = U[floor((indMin-1) / Ntotal) + 1, istep]
    
    mergedIndicies <- c(indC1,indC2)
    
    #select main cluster number to keep
    selectedInd <- min(mergedIndicies)
    
    #update sizes of clusters and current labels
    list[U,Nk] <- Merge_Clusters(mergedIndicies, U, Nk,  istep) 
    
    #update centroid of new cluster
    Z <- Update_Centroid(selectedInd,setdiff(mergedIndicies,selectedInd), U, Matrix, Z ,istep, p)
    
    #update weigths for new cluster
    W <- Update_Weights(selectedInd, U, Matrix, Z, W, istep, p)
    
    #compute number of iterations
    istep = istep + 1
    #estimate current distance between centroids
    AllDistances <- 
      Fast_Distance_Function(ZO,WO,NkO,
                             Z,W,Nk,
                             p,U,istep,
                             constraint_type,Constrains)
   
  } #main cycle
  
  Uind <- U[,NCOL(U)]
  #compute selhouette index
  silhouetteInd <- silhouetteMeasure(AllDistancesO, Uind)
  
  #return partitioning for last iteration
  return (list(U,silhouetteInd))
  
} # end of Ward_p


Merge_Clusters <- function(mergeIndicies, Um, Nk, istep){
  
  # merge two clustres with numbers k1min and k2min in new cluster with number k1min
  
  selectedInd <- min(mergeIndicies)
  unSelectedInd <- setdiff(mergeIndicies,selectedInd)
  
  istep1 <-istep+1
  
  Nk[selectedInd] <- sum(Nk[mergeIndicies])
  Nk[unSelectedInd] <- 0
  
  Um[ ,istep1] <- Um[ ,istep]
  Um[which(Um[,istep1] %in% unSelectedInd),istep1] <- selectedInd
  return(list(U = Um, Nk = Nk))
}

Fast_Distance_Function <- 
  function(ZO,WO,NkO,Z,W,Nk,p,U,istep,constraint_type="none",Constrains=c())
{
  dataTotal <- rbind(ZO,Z[Nk>1,])
  weightTotal <- rbind(WO,W[Nk>1,])
  sizeTotal <- c(NkO,Nk[Nk>1])
  
  
  Ntotal <- NROW(dataTotal)
  AllDistances <- matrix(0,Ntotal,Ntotal)
  #select all possible pair-wise combinations
  indCombinations <- combn(Ntotal,2)
  #get distance to selected pair
  combinedResult <- Get_Cluster_Distance(indCombinations[1,],
                                         indCombinations[2,], dataTotal,
                                         weightTotal, sizeTotal, p)
  
  #upper-triangle indices
  ind = (indCombinations[1,]-1)*Ntotal + indCombinations[2,]
  #lower-triangle indices
  indt = (indCombinations[2,]-1)*Ntotal + indCombinations[1,]
  #assgn results for pairs
  AllDistances[ind] <- combinedResult
  AllDistances[indt] <- combinedResult
  if(constraint_type == "ipoptim")
  {
    Constrains <- 
      constrains_vectorial(matrix_ = Constrains, number_individuals = NROW(AllDistances))
    AllDistances <- ipoptim(dissimilarity = AllDistances, 
                 constrains = Constrains, 
                 number_individuals = NROW(AllDistances), 
                 runs = 10,
                 maxiter = 1000000)
  }else if(constraint_type == "ultratran")
  {
    AllDistances <- 
      floyed_warshal_modified(dissimilarity = AllDistances,
                              constrains = Constrains)
  }
  AllDistances[cbind(1:NROW(AllDistances),1:NROW(AllDistances))] <- Inf
  AllDistancesFull <- matrix(Inf,NROW(Z),NROW(Z))
 
  #extract data for centroids only
  iInd <- 1:NROW(Z)
  if(length(Nk>1)>0)
  {
    iInd[Nk>1] <- ((1:length(which(Nk>1))) + NROW(ZO))
  }
  for (i in 1:NROW(Z))
  {
    for (j in 1:NROW(Z))
    {
      AllDistancesFull[i,j] <- AllDistances[iInd[i],iInd[j]]
    }
  }
  AllDistances <- AllDistancesFull
  removalInd <- which(! ((1:NROW(AllDistances)) %in% U[,istep]))
  AllDistances[removalInd,] <- Inf 
  AllDistances[,removalInd] <- Inf 
  return(AllDistances)
}

Update_Centroid <- function(k1min,k2min, U , Matrix,  Z ,istep, p){
  # update centroid vector Z[k1min, ] for new cluster after cluster merging
  istep1 <-istep+1
  index <- which(U[,istep1]==k1min)
  cluster <- Matrix[index,]
  
  lambdaRate <- 0.1
  decreaseLambdaRate <- 0.9
  centroid  <- Find_Optimal_Centroid(cluster, p, lambdaRate, decreaseLambdaRate)
  Z[k1min, ] <- centroid
  Z[k2min, ]  <- Inf
  return(Z)
}


Select_Cluster <- function(k1min,U, Matrix, istep) {
  # select n-th cluster from data
  
  istep1 <-istep+1
  index <- which(U[,istep1]==k1min)
  cluster <- Matrix[index,]
  
  return (cluster)
  
}

Get_Cluster_Distance <- function(i, j, Z, W, Nk, p){
  
  wij <- 0.5*(W[i, ] + W[j, ] )
  ni<- Nk[i]
  nj<- Nk[j]
  
  d <- (wij*(abs(Z[i, ] - Z[j, ])))^p
 
  dist <- (ni*nj/(ni+nj)) * rowSums(d)

  dist[i==j] = Inf
  return (dist)
}


Update_Weights <- function (k1min, U, Matrix, Z, W, istep, p){
  # update weights matrix after cluster merging
  istep1 =istep+1
  clusterk <- Select_Cluster(k1min, U, Matrix, istep)
  centroid <- as.numeric(Z[k1min,])
  W[k1min,]  <- Weights_p(clusterk, centroid, p)
  return(W)
}
