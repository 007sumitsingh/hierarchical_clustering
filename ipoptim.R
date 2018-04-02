#' Packages needed to run C ++ code from R
#' You need to install the C ++ compiler (gcc) and language (g ++). If it is in linux, 
#' Can be installed from the package manager "synaptic"
#' In R, install.packages("Rcpp") and install.packages("inline")

# ......................... IPopTim (C++ Codes with Rcpp and RcppArmadillo packages) ............................
# C++ codes
src <- '
arma::mat C = as<arma::mat>(C_);
arma::mat E = as<arma::mat>(E_);
arma::vec d = as<arma::vec>(d__);
arma::vec u = as<arma::vec>(u_);
int r = C.n_rows;
arma::vec a = d;
arma::vec s = a;
int t = 0;
double euc = 1;   //Distance for convergence
Function max("max");
Function print("print");
while (t < r || euc >= 0.0000000001) 
{
  arma::mat d_ = a;
  int p = t % r;          //Algorithm line 2

  double h = (u(p) / 2);
  s = a + C.row(p).t() * h;    //Algorithm line 3
  for (int q = 0; q <r; ++q) 
  {
    //Algorithm line 4
    if (q == p) 
    {
      //Algorithm line 5
      arma::mat C_transpose = C.t();
      u(q) = as<double>(max( wrap(0), wrap(2 * ((C.row(q) * s) / ((C.row(q)) * E * C_transpose.col(q)))) ));     //Algorithm line 6
    }
  }

  h = u(p) / 2;
  a = s - C.row(p).t() * h;     //Algorithm line 11
  t = t + 1;  
  euc = 0;
  for (int i = 0; i < r; ++i) 
  {
    euc = euc + pow((a(i) - d_(i)), 2);
  }
  euc = sqrt(euc);          //Calulating euclidean distance between old a and new a
}

return ( wrap(a) );
'
# R function
ipoptim_cplusplus <- cxxfunction(signature(C_ = "numeric", E_ = "numeric", d__ = "numeric", u_ = "numeric"),
                                 body = src, 
                                 plugin = "RcppArmadillo")





#' IPopTim Algorithm. Calculates an approximate dissimilarity matrix considering the constraints
#' It raises the problem of finding the approximate matrix as a restricted nonlinear optimization problem
#' It uses the Karush-Kuhn tucker method for optimization and the iterative projection algorithm
#'
#' @param dissimilarity Matrix dissimilarity matrix of individuals
#' @param constrains Vector representation of the constraint matrix
#' @param number_individuals Number of individuals in the original date to be compared 
#' @param runs Number of executions of the iterative projection for the ultrametricity constraints
#' @param maxiter Maximum amount of iterations for the iterative projection for the ultrametricity constraints
#' 
#' @return Estimated matrix for the matrix of dissimilarities
#' @export
#'
#' @examples
ipoptim <- function(dissimilarity = NULL, constrains = NULL, number_individuals = 0, runs = 1, maxiter = 1000)
{
  if(class(dissimilarity) != "matrix")
  {
    stop("'dissimilarity' must be square and symmetrical matrix with dissimilarity of individuals")
  }
  
  # Projection on the triple-wise restriction set
  d <- vector_from_matrix(dissimilarity)
  C <- constrains
  C <- t(C) # C as contraint matrix. rows = constrains, cols = length(d): matrix_dissimilarity as vector
  E <- diag(number_individuals * (number_individuals - 1) / 2) # Elementry marix
  d <- t(d)
  d <- t(d)
  u <- vector(mode = "double", length = nrow(C)) # Kuhn-Tucker vector
  a <- ipoptim_cplusplus(C_ = C, E_ = E, d__ = d, u_ = u)

  a <- matrix_from_vector(vector_ = a, number_individuals = nrow(dissimilarity))
  # Projection on the ultrametricity restriction set 
  {
    options(warn = -1)
    distance_between_a_aultrametric <- .Machine$integer.max   # Maximum value of an integer
    ultrametric <- matrix(NA, ncol = 0, nrow = 0)
    for(i in 1:runs)
    {
      ultrametric2 <- ls_fit_ultrametric(as.dist(a), 
                                      method = "IP", 
                                      control = list(verbose = FALSE, runs = 1, maxiter = maxiter))
      distance <- sqrt(sum((a - ultrametric2)^2)) # Euclidean distance 
      if(distance < distance_between_a_aultrametric)
      {
        ultrametric <- ultrametric2
        distance_between_a_aultrametric <- distance
      }
    }
  }
  ultrametric <- as.matrix(ultrametric)
  ultrametric <- round(ultrametric, 10)
  ultrametric[1,1] <- as.character(ultrametric[1,1])
  
  names_table_ultrametric <- names(table(ultrametric))
  # We convert all distances in the order of the group in which they form
  for(i in 2:length(names_table_ultrametric))
  {
    ultrametric[ ultrametric == names_table_ultrametric[i] ] <- i - 1
  }
  ultrametric <- apply(ultrametric, 1, as.numeric)
  return(ultrametric)
}



#' Generating the Constraint (1,-1,0) matrix. 
#' Vector representation of triple-wise constraints
#'
#' @param matrix_ Matrix with constraints
#' @param number_individuals Number of individuals in the original date to be compared
#'
#' @return Matrix with triple-wise constraints in a vector representation
#'
#' @examples
constrains_vectorial <- function(matrix_ = NULL, number_individuals = 0)
{
  r <- nrow(matrix_)
  mat <- matrix(0, (number_individuals * (number_individuals - 1) / 2), r) # The newly converted constraint matrix
  for (i in 1:r) 
  {
    if (matrix_[i, 1] < matrix_[i, 2]) 
    {
      min <- matrix_[i, 1]
      max <- matrix_[i, 2]
    }
    else 
    {
      min <- matrix_[i, 2]
      max <- matrix_[i, 1]
    }
    ind_min <- (min - 1) * number_individuals - ((min - 1) * min / 2)
    ind_min = ind_min + (max - min)
    mat[ind_min, i] <- 1
    if (matrix_[i, 1] < matrix_[i, 3]) 
    {
      min <- matrix_[i, 1]
      max <- matrix_[i, 3]
    }else 
    {
      min <- matrix_[i, 3]
      max <- matrix_[i, 1]
    }
    ind_max <- (min - 1) * number_individuals - ((min - 1) * min / 2)
    ind_max <- ind_max + (max - min)
    mat[ind_max, i] <- -1
  }
  return (mat)
}




# ......................... vector_from_matrix function (C++ Codes with Rcpp package) ............................
# C++ codes
src <- '
Rcpp::NumericMatrix matrix_(data);
int nrows = matrix_.nrow();
int ncolumns = matrix_.ncol();
int k = 0;
int length_data = nrows * (nrows - 1) / 2;
Rcpp::NumericVector d_vector(length_data);

int stop_rows = (nrows - 1);
for (int i = 0; i < stop_rows; ++i) 
{
  for (int j = (i + 1); j < ncolumns; ++j) 
  {
  d_vector[k] = matrix_(i, j);
  k = k + 1;
  }
}
return (d_vector);
'
# R function
vector_from_matrix_cpluplus <- cxxfunction(signature(data = "numeric"), body = src, plugin="Rcpp")


#' Function that represents in vector form a matrix. The matrix must be a square and symmetric matrix.
#' The upper diagonal is used to calculate the vector and the path is done by rows.
#'
#' @param matrix square matrix 
#'
#' @return Vector containing all the elements of the upper diagonal of the matrix
#' @export
#'
#' @examples
vector_from_matrix <- function(matrix_ = NULL)
{
  if(nrow(matrix_) != ncol(matrix_))
  {
    stop("The matrix must be square and symmetrical")
  }
  d_vector <- vector_from_matrix_cpluplus(data = matrix_)
  return(d_vector)
}



# ......................... matrix_from_vector function (C++ Codes with Rcpp package) ............................
# C++ codes
src <- '
Rcpp::NumericVector vector_(data);
Rcpp::NumericMatrix matrix_(as<int>(number_individuals));
int k = 0;
int number_individuals_ = as<int>(number_individuals);

for (int i = 0; i < number_individuals_; ++i) 
{
  for (int j = 0; j < number_individuals_; ++j) 
  {
    if (i < j) 
    {
      matrix_(i, j) = vector_(k);
      matrix_(j, i) = vector_(k);
      k = k + 1;
    }
  }
}

return (matrix_);
'
# R function
matrix_from_vector_cpluplus <- cxxfunction(signature(data = "numeric", number_individuals = "integer"), body = src, plugin="Rcpp")



#' A function that creates a square and symmetric matrix from a vector. 
#'
#' @param matrix matriz cuadrada 
#'
#' @return Vector containing all the elements of the upper diagonal of the matrix
#'
matrix_from_vector <- function(vector_ = NULL, number_individuals = 0)
{
  if( (number_individuals * (number_individuals - 1) / 2) != length(vector_))
  {
    stop("The vector must represent the upper diagonal of a square matrix")
  }
  matrix_ <- matrix_from_vector_cpluplus(data = vector_, number_individuals = number_individuals)
  return(matrix_)
}

