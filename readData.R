testReadData <- function(){
  data(iris)
  database <- iris
  class_database <- 5
  labels_class_database <- c(1,2,3)
  data        <- database[,-class_database]
  class_data  <- database[, class_database]
  label_x <- 10
  Constrains <- constrains_triple_wise(data = data, class = class_data, 
                                       percentaje = label_x)
  Class   <- factor(class_data, labels = labels_class_database)
  
  reqK <- length(unique(class_data))
  p_set <- seq(1,5,0.1)  
  p <- 1.5
  list[U,optimalP,SI] <- Ward_p_range(data, reqK, p_set,
                              constraint_type="none",Constrains=c())
  print(c("optimalP=",optimalP, ", optimalSI=",SI))
  #list[U,SI]<- Ward_p(data, reqK,p,
  #                    constraint_type="ipoptim",Constrains=Constrains)
  #list[U,SI]<- Ward_p(data, reqK,p,
  #                    constraint_type="ultratran",Constrains=Constrains)
  list[U,SI]<- Ward_p(data, reqK,p,
                      constraint_type="none",Constrains=c())
  fWardp <- FScore_Wardp(as.factor(class_data),U[,NCOL(U)])
  AdjWardp <- adj.rand.index(as.factor(class_data),U[,NCOL(U)])
  fWardpG <- FScore_Wardp_Total(as.factor(class_data),U)
  
  print(c("Ward-p single: ", fWardp))
  print(c("Ward-p maximum:", fWardpG))
  B<- dissimilar(mydata = data, method = "euclidean")
  clusters  <- hclust(as.dist(B), method = "ward.D")
  cutree_matrix <- cutree(clusters, 1:length(as.factor(class_data)), NULL)
  fWardpGD <- FScore_Wardp_Total(as.factor(class_data),cutree_matrix)
  
  FscoreWardD <- FScore_hierarchical_clustering(class = Class, tree = clusters)
  print(c("score by name ward-d:",  as.numeric(FscoreWardD)))
  print(c("score independent ward-d:",fWardpGD))
  
  wardDSI <- silhouetteMeasure(B, as.matrix(cutree_matrix)[,reqK])
  
  print(c("silhouetteMeasure ward-p:", SI))
  print(c("silhouetteMeasure ward-d:", wardDSI))
  
}