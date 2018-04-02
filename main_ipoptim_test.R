run_main_ipoptim <- function(){


# CONFIGURATION OF VARIABLES USED. Indicate the methods to be used, the number of times the algorithm is to be executed, etc.
name_ejecution_script <- "iris_ipoptim_10"       # Name for script execution. For example, "complete data iris"
name1_ejecution_script <- "iris_ipoptim_10_p"
RUNS                  <- 2                          # Number of runs of the algorithms for each configuration
database              <- iris                       # Database to use
class_database        <- 5                          # Column in the database containing the class
labels_class_database <- c(1, 2, 3)                 # Class tags. Used to pass classes to the optimization method. Must be a integer
label                 <- c(10)  # Percentage of restrictions
method.diss           <- c("manhattan",
                           "euclidean",                     
                           "std euclidean", 
                           "mahanalobis", 
                           "cosine", 
                           "correlation", 
                           "chebyshev", 
                           "Canberra", 
                           "bray-Curtis",
                           "spearman") #c("manhattan", "euclidean", "std euclidean", "mahanalobis", "cosine", "correlation", "spearman", "chebyshev", "Canberra", "bray-Curtis")
method                <- c("single", "complete", "average", "median", "centroid", "ward.D","ward.P") 

reqK <- length(unique(class_data))
p_set <- 1.5#seq(1,5,0.1) 

# Select the data to use depending on the configuration of the script
data        <- database[,-class_database]
class_data  <- database[, class_database]


# Storing the results in an object of R
table.1 <- data.frame(Dissimilarity.Method = character(),
                      Label = character(),
                      Measurement = character(),
                      single = double(),
                      complete = double(),
                      average = double(),
                      median = double(),
                      centroid = double(),
                      ward.D = double(),
                      stringsAsFactors = FALSE)

table.2 <- data.frame(Dissimilarity.Method = character(),
                      Label = character(),
                      Measurement = character(),
                      single = double(),
                      complete = double(),
                      average = double(),
                      median = double(),
                      centroid = double(),
                      ward.D = double(),
                      stringsAsFactors = FALSE)
aggr <- 0


t <- proc.time() # Start cronometer
for (method.diss in method.diss) 
{
  # Matrix of Dissimilarity
  B <- dissimilar(mydata = data, method = method.diss)
  
  for (label_x in label) 
  {
    file.name <- paste(method.diss, "Label ", label_x, "%.csv")
    file.name1 <- paste(method.diss, "Label ", label_x, "%.csv")
    file.name1 <- paste("adjrand",file.name1)
    
    # making an empty data frame for taking iterative results
    v <- data.frame(Run = character(),
                    single = double(),
                    complete = double(),
                    average = double(),
                    median = double(),
                    centroid = double(),
                    ward.D = double(),
                    stringsAsFactors = FALSE)
    
    v1 <- data.frame(Run = character(),
                    single = double(),
                    complete = double(),
                    average = double(),
                    median = double(),
                    centroid = double(),
                    ward.D = double(),
                    stringsAsFactors = FALSE)
    
    for (rep in 1:RUNS) 
    {
      y <- paste("Run: ", rep)
      y1 <- paste("Run: ", rep)
      ## ModifiedFloyed
      ## input : dissimilarity(B), label
      ## output : B
      
      # Generating Restrictions
      Constrains <- constrains_triple_wise(data = data, class = class_data, percentaje = label_x)
      
      # Generating the Constraint (1,-1,0) matrix
      Constrains_vectorial <- constrains_vectorial(matrix_ = Constrains, number_individuals = nrow(B))
      
      # Iterative Projection to minimize least-square error
      D <- ipoptim(dissimilarity = B, 
                   constrains = Constrains_vectorial, 
                  number_individuals = NROW(B), 
                   runs = 10,
                 maxiter = 1000000)
     
      # Measuring FScore
      Class   <- factor(class_data, labels = labels_class_database)
      Fscore  <- c(single = NA, complete = NA, average = NA, median = NA, centroid = NA, ward.D = NA, ward.P = NA)
      AdjRand <- c(single = NA, complete = NA, average = NA, median = NA, centroid = NA, ward.D = NA, ward.P = NA)
      for (i in 1:length(method)) 
      {
        if(method[i] != "ward.P"){
          clusters  <- hclust(as.dist(D), method = method[i])
          Fscore[i] <- FScore_hierarchical_clustering(class = Class, tree = clusters)
          rezult <- Get_Group_By_Number(class = Class, tree = clusters,(length(Class)-length(unique(Class)) + 1))
          AdjRand[i] <- adj.rand.index(Class,rezult)
        }
        else{
          list[U,p,SI] <- Ward_p_range(data, reqK, p_set,
                                       constraint_type="ipoptim",Constrains=Constrains)
          Fscore[i] <- FScore_Wardp (Class, U)
          AdjRand[i] <- adj.rand.index(Class,U)
        }
      }
  
    }
  }
}

}







