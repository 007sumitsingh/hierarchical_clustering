#' Packages needed to run C ++ code from R
#' You need to install the C ++ compiler (gcc) and language (g ++). If it is in linux, 
#' Can be installed from the package manager "synaptic"
#' In R, install.packages("Rcpp") and install.packages("inline")


# ......................... Modified Floyd Warshal Algorithm (C++ Codes with Rcpp package) ............................
# C++ codes
src <- '
Rcpp::NumericMatrix matrix_transitive(data);
Rcpp::NumericMatrix matrix_constrains(constrains);
int nrows   = matrix_transitive.nrow();
int h       = matrix_constrains.nrow();
int mincon  = INT_MAX;  
Function min("min");
Function max("max");

// The parent of a m_ij only changes when m_ij is different from the previous value.
double previus_matrix_transitive_ij = INT_MAX;

// X_i value of the constraint. This value is used to assign as parent when a restriction is applied
int x_i_constraint = INT_MAX;

for (int k = 0; k < nrows; ++k) 
{
  for (int i = 0; i < nrows; ++i) 
  {
    for (int j = 0; j < nrows; ++j)
    {
      for (int l = 0; l < h; ++l)
      {
          if(matrix_constrains(l,0)==(i+1) && matrix_constrains(l,1)==(j+1))
	        {
            mincon = as<int>( min(wrap(mincon), wrap(matrix_transitive(i, matrix_constrains(l,2)-1)) ) );
            x_i_constraint = matrix_constrains(l,0);
	        }
      }

      previus_matrix_transitive_ij = matrix_transitive(i,j);
      // Validate whether a constraint is applied to a cell in the array
      // There is no restriction on the cell
      if(mincon == INT_MAX)
      {
        // Calculate the new m_ij
        // mincon can be removed since it is infinite. That is, there is no restriction for the cell being worked on
        matrix_transitive(i,j) = as<double>( min(wrap(matrix_transitive(i,j)), max(wrap(matrix_transitive(i,k)), wrap(matrix_transitive(k,j))), wrap(mincon)) );  
        
        // We verify if m_ij has changed in order to modify the parent matrix
        if(previus_matrix_transitive_ij != matrix_transitive(i,j))
        {
          // As the matrix is symmetric, we have that if the lower diagonal changes, the upper diagonal must change
          matrix_transitive(j,i)  = matrix_transitive(i,j);
        }
      }else
      {
        // Calculate the new m_ij
        matrix_transitive(i,j) = as<double>( min(wrap(matrix_transitive(i,j)), max(wrap(matrix_transitive(i,k)), wrap(matrix_transitive(k,j))), wrap(mincon)) );

        // We verify if m_ij has changed in order to modify the parent matrix
        if(previus_matrix_transitive_ij != matrix_transitive(i,j))
        {
          // The matrix of distance we are converting it into a matrix of order of creation of the groups.
          matrix_transitive(i,j)  = x_i_constraint;

          // Since the matrix is symmetric and a constraint applied on ji must also be applied on ij, we have to
          matrix_transitive(j,i)  = matrix_transitive(i,j);
        }
      }
       
      mincon = INT_MAX; 
    } // End for to j      
  } // End for to i
} // End for to k

return (matrix_transitive);
'
# R function
floyed_warshal_modified_cplusplus <- cxxfunction(signature(data = "numeric", constrains = "numeric"), body = src, plugin="Rcpp")



#' Modified Floyd Warshal Algorithm
#'
#' @param dissimilarity Matrix of dissimilarity of individuals
#' @param constrains Constrains on clustering
#'
#' @return matrix Contains the minimal transitive dissimilarity from the matrix of dissimilarity
#' @export
#'
#' @examples
floyed_warshal_modified <- function(dissimilarity = NULL, constrains = NULL)
{
  if(is.null(dissimilarity))
  {
    stop("'dissimilarity is NULL'")
  }
  
  if(is.null(constrains))
  {
    stop("'constrains is NULL'")
  }
  ultratran <- floyed_warshal_modified_cplusplus(dissimilarity, constrains)
  ultratran <- round(ultratran, 10)
  ultratran[1,1] <- as.character(ultratran[1,1])
  
  names_table_ultratran <- names(table(ultratran))
  # We convert all distances in the order of the group in which they form
  for(i in 2:length(names_table_ultratran))
  {
    ultratran[ ultratran == names_table_ultratran[i] ] <- as.character(i - 1)
  }
  ultratran <- apply(ultratran, 1, as.numeric)
  return (ultratran)
}


testing<-matrix(c(0,1,7,9,11,11,1,0,7,9,12,12,7,7,0,6,10,10,9,9,6,0,8,8,11,12,10,8,0,2,11,12,10,8,2,0), nrow = 6, ncol=6, byrow = TRUE)
q<-matrix(c(1,2,3,1,3,5,3,1,4,4,5,1), nrow=4, ncol = 3, byrow = TRUE)
q<-matrix(c(1,2,3,1,3,5,3,1,4,4,5,1), nrow=4, ncol = 3, byrow = TRUE)
z<-floyed_warshal_modified(testing, q)
z2<-as.dist(z)
clusters<-hclust(z2, method = "single")
plot(clusters)
