# Ward p weihgts

Weights_p <-function(cluster,centroid,p){
 # Returns weight vector for input cluster, weight[j] -weight for j-th feature
  
  #compute average distances to center for all features
  distToCenter <- D0(cluster,centroid,p)
  
  #if average distance is always 0, keep initial weigths 
  if (length(distToCenter[distToCenter>0])==0)
  {
    weights <- rep(1/length(distToCenter),length(distToCenter))
    return(weights)
  }
  
  #initialization
  weights <- rep(0.0,NCOL(cluster))
  if( p>1){
    p1 <- 1/(p-1)
    #compute weight for each feature separately
    for( v in 1:NCOL(cluster)) {
      s <- sum((distToCenter[v]/distToCenter[distToCenter>0])^p1) + 10e-20
      weights[v] <- 1.0/s
    } # v
  } # if p > 1
   else{
     #separate case p = 1 to avoid division by 0
     #we select single weight for closest distance are 1.0, others assigned to 0 
     index <- which.min(distToCenter)
     weights[index] <- 1.0
      } 
  
   return(weights)
} 

D0 <- function(cluster, centroid, p){
  c0 <- t(replicate(NROW(cluster),centroid))
  s <- colSums((abs(cluster-c0))^p) + 10e-20
  return (s) 

}
