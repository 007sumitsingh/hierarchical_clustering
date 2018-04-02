#' Packages needed to run C ++ code from R
#' You need to install the C ++ compiler (gcc) and language (g ++). If it is in linux, 
#' Can be installed from the package manager "synaptic"
#' In R, install.packages("Rcpp") and install.packages("inline")




# ......................... Clausure transitive (C++ Codes with Rcpp package) ............................
# C++ codes
src <- '
Rcpp::NumericMatrix matrix_adjacency_(matrix_adjacency);
int nrows = matrix_adjacency_.nrow();

for(int y = 0;  y < nrows; ++y)
{
  for(int x = 0; x < nrows; ++x)  
  {
    if(matrix_adjacency_(x,y) == 1)  
    {
      for(int j = 0; j < nrows; ++j) 
      {
        if(matrix_adjacency_(y,j) == 1)    
        {
          matrix_adjacency_(x,j) = 1;
        }
      }  
    }
  }
}

return (matrix_adjacency_);
'
# R function
closure_transitive_cplusplus <- cxxfunction(signature(matrix_adjacency = "numeric"), body = src, plugin="Rcpp")


#' Calculates the transitive closure of a graph from its adjacency matrix. Floyd's Algorithm
#' The relationship between individuals is a binary relation
#'
#' @param matrix Adjacency matrix
#'
#' @return matrix Transitive closure matrix
#'
closure_transitive <- function(matrix = NULL)
{
  if(is.null(matrix) | class(matrix) != "matrix")
  {
    stop("'matrix' must be a matrix")
  }
  closure_transitive_cplusplus(matrix)
  
  return (matrix)
}



#' Function to extract a set of implicit constraints of type c_1 = (x_i, x_j, x_k) and c_2 = (x_i, x_k, x_l)
#' 
#' @param constrains constrains matrix
#'
#' @return matrix constrains matrix expanded
#'
expand_constrains_clausure_transitive <- function(constrains = NULL)
{
  if(is.null(constrains) | class(constrains) != "matrix")
  {
    stop("'constrains' must be a matrix")
  }
  #install.packages("Rcpp")
  #Rcpp::evalCpp("2+2")
  constrains_data_frame <- as.data.frame(constrains)
  
  # We create the set A of all ordered pairs x_i, x_j or x_i, x_k that form part of the constraints x_i, x_j, x_k
  # On this set A is defined the relation R of A in A such that 
  # a R b if and only if a and b form a constraint c of the type x_i, x_j, x_k of the set of constraints C (constrains)
  # For all a, b that belongs to the set A and c that belongs to the set C (constrains)
  A <- rbind(as.matrix(constrains_data_frame[,1:2]), as.matrix(constrains_data_frame[,c(1,3)]))
  A <- unique(A) # There can be no repeated elements in the set of ordered pairs
  A <- as.data.frame(A)
  row.names(A) <- 1:nrow(A)
  
  # Matrix of the relation R on the set A. Used to find the implicit relations
  A_relation_matrix <- matrix(0, nrow = nrow(A), ncol = nrow(A))
  A_relation_matrix <- as.data.frame(A_relation_matrix)
  row.names(A_relation_matrix) <- row.names(A)
  names(A_relation_matrix) <- row.names(A_relation_matrix)
  
  # We traverse the relation matrix, placing 1 'in the position where an ordered pair is related to another ordered pair of the set A
  for(i in 1:nrow(constrains_data_frame))
  {
    row <- row.names(A[ A[,1] == constrains_data_frame[i,1] & A[,2] == constrains_data_frame[i,2], ])
    col <- row.names(A[ A[,1] == constrains_data_frame[i,1] & A[,2] == constrains_data_frame[i,3], ])
    A_relation_matrix[row,col] <- 1
  }
  
  constrains_expanded         <- matrix(0,nrow = 0, ncol=3)
  constrains_expanded         <- as.data.frame(constrains_expanded)
  names(constrains_expanded)  <- c("x_i", "x_j", "x_k")
  
  # We find the transitive closure on the relation matrix
  A_relation_matrix <- closure_transitive(matrix = as.matrix(A_relation_matrix))
  
  # Create the extended constraint set with the implicit constraints from the transitive cluster matrix
  for(i in 1:nrow(A_relation_matrix))  
  {
    for(j in 1:nrow(A_relation_matrix))  
    {
      if(A_relation_matrix[i,j] == 1)
      {
        vector <- c(A[i,][[1]], A[i,][[2]], A[j,2])
        names(vector) <- names(constrains_expanded)
        constrains_expanded <- rbind(constrains_expanded, vector)
      }
    }
  }

  return(as.matrix(constrains_expanded))
}




#' Function that removes the conflicts of the set of constraints
#'
#' @param constrains Matrix with the set of constraints
#'
#' @return Matrix of constraints without conflicts
#'
conflicts_removal <- function(constrains = NULL)
{
  if(is.null(constrains) | class(constrains) != "matrix")
  {
    stop("'constrains' must be a matrix")
  }

  constrains_data_frame <- as.data.frame(constrains)
  
  # Remove Duplicates
  constrains_data_frame <- unique(constrains_data_frame)

  # Explicit conflicts
  # c_1 = (x_i, x_j, x_k) and c_2 = (x_j, x_k, x_i) 
  constrains_data_frame <- cbind(constrains_data_frame, id = 1:nrow(constrains_data_frame)) # Used to remove restrictions
  conflictive <- TRUE
  while(conflictive == TRUE)
  {
    constrains_auxiliar <- constrains_data_frame
    conflicts           <- data.frame() # Save all restrictions are conflicting
    conflicts_auxiliar  <- data.frame()
    
    for (i in 1:nrow(constrains_data_frame)) 
    {
      # We extract the subset of the constraints that form an explicit conflict
      conflicts_auxiliar <- constrains_auxiliar[(constrains_data_frame[i,1] == constrains_auxiliar[,2]) &
                                       (constrains_data_frame[i,2] == constrains_auxiliar[,3]) &
                                       (constrains_data_frame[i,3] == constrains_auxiliar[,1]) ,]  
      if(nrow(conflicts_auxiliar) > 0)
      {
        conflicts <- rbind(conflicts, conflicts_auxiliar, constrains_data_frame[i,])
      }
    }
    
    # remove constrains conflictive
    if(nrow(conflicts) > 0)
    {
      index_constrains_conflictive  <- sample(1:nrow(conflicts), 1)      # We randomly select the constraint to remove
      constrains_to_removal         <- conflicts[index_constrains_conflictive,]
      constrains_data_frame <- constrains_data_frame[constrains_data_frame[,4] != constrains_to_removal[1,4], ]    
    }else
    {
      conflictive <- FALSE
    }
  }
  return (as.matrix(constrains_data_frame[,-4]))
}





#' Function to calculate triple-wise constraints
#' 
#' @param data Data.frame with the data on which the restrictions will be generated. 
#' If the data.frame is a subset of another data.frame it is necessary to assign a new number to the names of the data.frame rows before calling this function.
#' @param class vector with classes of individuals. If the vector is not specified, the constraints are generated in a totally random way
#' @param percentaje Percentage of constraints generated. It is based on the number of individuals.
#' @param proportional logical Indicates whether the sampling of the individuals is done in proportion to the quantity that exists in each class.
#' @param number_constraints_fixed Integer that represents the number of constraints to be generated. 
#' If this value is given, the 'percentage' parameter is ignored.
#' @return matrix Matrix with triple-wise constraints
#'
#' @examples
constrains_triple_wise <- function(data = NULL, class = NULL, percentaje = 10, proportional = FALSE, number_constraints_fixed = NULL)
{
  if(is.null(data) | class(data) != "data.frame")  
  {
    stop("'Data' must be a data.frame with the data to be grouped")
  }

  if((class(percentaje) != "numeric" & class(percentaje) != "integer") | percentaje <= 0 )	# percentaje > 100
  {
    stop("'percentaje' Must be numeric or integer type. Must be a number in the range (0,100]")
  }

  if(!is.null(class) & nrow(data) != length(class))
  {
    stop("The length of the vector 'class' must equal the number of rows of 'data'")
  }
  
  if(!is.null(number_constraints_fixed))
  {
    if(class(number_constraints_fixed) != "integer" & class(number_constraints_fixed) != "numeric")
    {
      stop("'number_constraints_fixed' must be a integer")
    }
  }
  
  number_individuals  <- nrow(data)
  if( is.null(number_constraints_fixed) )
  {
    number_constrains <- ceiling(number_individuals*(percentaje/100))
  }else
  {
    number_constrains <- number_constraints_fixed
  }
  data_auxiliar       <- data
  
  if(is.null(class))
  {
    # When individuals do not have a class tag, the constraints are constructed completely random
    index_individuals <- 1:number_individuals
    x_i               <- sample(index_individuals, number_constrains, replace=TRUE)
    x_j               <- sample(index_individuals, number_constrains, replace=TRUE)
    x_k               <- sample(index_individuals, number_constrains, replace=TRUE)
    constrains <- data.frame(x_i_class = rep(NA, number_individuals), 
                             x_j_class = rep(NA, number_individuals), 
                             x_k_class = rep(NA, number_individuals), 
                             x_i = x_i, 
                             x_j = x_j, 
                             x_k = x_k)
  }else
  {
    # When individuals have a class tag, the constraints are constructed in the form: c = (class x, class x, class y)
    # The original data is added by the sample and class columns to perform the sampling of individuals.
    data_auxiliar <- cbind(sample = row.names(data), data_auxiliar, class) #1:number_individuals
    class_names   <- names(table(class))
    
    if(proportional == FALSE)
    {
      # Vectors that indicate the classes that will form the constraints
      x_i_class <- sample(class_names, number_constrains, replace=TRUE)
      x_j_class <- x_i_class
      x_k_class <- c()
      # The third component of the constraint must be different from the first two
      for(i in 1:length(x_i_class))
      {
        x_k_class[i] <- sample(class_names[class_names != x_i_class[i]], 1, replace=TRUE)
      }
      # Constrains
      constrains <- data.frame(x_i_class, x_j_class, x_k_class, 
                               x_i = rep(0,number_constrains), 
                               x_j = rep(0,number_constrains), 
                               x_k = rep(0,number_constrains))
      
      # Selection of the individuals that are part of each of the restrictions
      for(i in 1:length(names(table(class))))
      {
        sample_individuals_vector         <- data_auxiliar[data_auxiliar$class == names(table(class))[i], "sample"]
        sample_individuals_vector_numeric <- as.numeric(levels(sample_individuals_vector)[sample_individuals_vector])
        
        # For x_i
        if(!is.na(table(constrains[,1])[names(table(class))[i]]))
        {
          constrains[constrains$x_i_class == names(table(class))[i], 4] <- sample( sample_individuals_vector_numeric, table(constrains[,1])[names(table(class))[i]], replace=TRUE ) 
        }
        
        # For x_j
        if(!is.na(table(constrains[,2])[names(table(class))[i]]))
        {
          constrains[constrains$x_j_class == names(table(class))[i], 5] <- sample( sample_individuals_vector_numeric, table(constrains[,2])[names(table(class))[i]], replace=TRUE )
        }
        
        # For x_k
        if(!is.na(table(constrains[,3])[names(table(class))[i]]))
        {
          constrains[constrains$x_k_class == names(table(class))[i], 6] <- sample( sample_individuals_vector_numeric, table(constrains[,3])[names(table(class))[i]], replace=TRUE )
        }
      }
    }else
    {
      stop("The case when the distribution of individuals per class is not balanced yet has not been programmed")
    }
  }

  # Expand the set of initial constraints by transitive closure. The Warshall algorithm is used. Implied constraints are found.
  constrains_expanded_matrix <- expand_constrains_clausure_transitive(constrains = as.matrix(constrains[,4:6]))
  rm(constrains, data_auxiliar)
  
  # Conflicts removal
  constrains_conflict_removed <- conflicts_removal(constrains_expanded_matrix)

  dimnames(constrains_conflict_removed) <- list(1:nrow(constrains_conflict_removed), c("x_i", "x_j", "x_k"))
  return (constrains_conflict_removed)
}




