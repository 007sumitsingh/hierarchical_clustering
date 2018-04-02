
Find_Optimal_Centroid <- function(cluster, p, lambdaRate, decreaseLambdaRate){
  # RETURNS OPTIMAL CENTROID VECTOR FOR GIVEN CLUSTER
  cent_vector <- c()

  for( i in 1 : NCOL(cluster)){
   #sort feature values  
   yreg <- sort(cluster[ ,i])
   #determine initial value minimising the distance
   c0Dist <- sapply(X=yreg,FUN=function(x){return(distance(yreg,p,x))})
   c0Ind <- which.min(c0Dist)
   c0 <- yreg[c0Ind]
   featureRange <- (max(yreg) - min(yreg))
   initLambda <- lambdaRate * featureRange
   
   if (initLambda < 10e-10)
   {
     cent_vector[i] <-c0
     next
   }
   
   # acceptable boundaries in the code
   leftBound <-  yreg[Find_Left_Bound(c0Dist,c0Ind)]
   rightBound <- yreg[Find_Right_Bound(c0Dist,c0Ind)]
   if((rightBound - leftBound) < 10e-10)
   {
     cent_vector[i] <-c0
     next
   }
  center <-  Steepest_Descent_V1(yreg, p, c0, leftBound, rightBound, initLambda, decreaseLambdaRate)
    
  cent_vector[i] <-center
  }
  return(cent_vector)

}

Steepest_Descent_V1 <- function(yreg, p, c0, left,right,initialLambda, decreaseLambdaRate){
  
  niter <- 0
  Max_N_Iter  <- 100
  
  
  lambda <- initialLambda
  h <- 0.01*lambda
  
  while(TRUE){
    
    gradient <- derive1(yreg, p, c0,h)
    del <- lambda*gradient 
    c1 <- c0 - del
    
    # need protection from infinity
    decreaseSteps <- 1
    while((c1 < left) || (c1 > right)){
      lambda <- decreaseLambdaRate*lambda
      del <- lambda*gradient 
      c1 <- c0 - del
      decreaseSteps <- decreaseSteps + 1
      if (decreaseSteps > 10)
      {
        return(c0)
      }
    } # end of step3
    
    niter <- niter + 1
    
    if(abs(c0 - c1) < 10e-10) {
      return(c1)
    }
    c0 <- c1
    if(niter >= Max_N_Iter){
      #warning("Number of iterations for steepest descent exceeded!")
      return(c0)
    }
    
    
 } # END OF MAIN CYCLE - STEP2
  
 
}

  
derive1 <- function(yreg, p, center,h){
  # calculate derivative with c
  dv <- rep(0,length(center))
  for (i in 1:length(center))
  {
    centerh<-center
    centerh[i] <- centerh[i] + h
    dv[i] <- (distance(yreg, p, centerh) - distance(yreg, p, center))/h
  }
  return(dv)
}


distance <- function(yreg, p, c){
  # Calculate the distance from c0 to vector yreg
  d <- sum((abs(yreg-c))^p)
  return (d)
}

  
  Find_Left_Bound <- function(c0Dist,c0Ind){
    index <- max(1,c0Ind-1)
    if (index == c0Ind) return(index)
    while((index>1) && (c0Dist[c0Ind]>=c0Dist[index]))
    {
      index = index - 1
    }
    return(index)
 }
  
 Find_Right_Bound <- function(c0Dist,c0Ind){
   n <- length(c0Dist)
   index <- min(length(c0Dist),c0Ind+1)
   if (index == c0Ind) return(index)
   while((index<n) && (c0Dist[c0Ind]>=c0Dist[index]))
   {
     index = index + 1
   }
   return(index)
 }
    
  
 
